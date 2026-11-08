"""Query every provider in parallel, rank the answers, then write the brief.

Cache hits skip the network. With a time budget, providers that have not
answered by the deadline are reported as "pending"; their threads keep running
and save to the cache, so the next lookup has the fuller answer.
"""

from __future__ import annotations

import sys
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from dataclasses import dataclass, field
from pathlib import Path

from snp_brief.brief import assemble_brief, brief_sentences
from snp_brief.cache_store import cache_path, load_cache, save_cache
from snp_brief.models import ClinVarRecord, DbSnpRecord, GnomadVariant, Paper, UniProtHit
from snp_brief.providers import (
    BY_KEY,
    GeneSummary,
    GwasResult,
    VepAnnotation,
    extras_enabled,
    fetch_citations,
    fetch_gene_summary,
    fetch_gwas,
    fetch_vep,
    parse_citations,
    parse_gene_summary,
    parse_gwas,
    parse_vep,
    Eqtl,
    L2G,
    LdProxy,
    fetch_gtex,
    fetch_ld,
    fetch_litvar,
    fetch_open_targets,
    gtex_variant_id,
    ot_variant_id,
    parse_gtex,
    parse_ld,
    parse_litvar,
    parse_open_targets,
)
from snp_brief.sources.clinvar import fetch_clinvar, parse_clinvar
from snp_brief.sources.dbsnp import fetch_dbsnp, parse_dbsnp
from snp_brief.sources.europepmc import fetch_europepmc, parse_europepmc
from snp_brief.sources.gnomad import fetch_gnomad, parse_gnomad
from snp_brief.sources.uniprot import fetch_uniprot, parse_uniprot

# Shared so a provider that misses the deadline keeps running after gather returns.
_POOL = ThreadPoolExecutor(max_workers=12, thread_name_prefix="provider")
_LOCK = threading.Lock()


@dataclass
class SourceCall:
    """How one source answered: origin is "cache", "live", "no cache", "failed", or "pending"."""

    origin: str
    retrieved: str | None = None
    error: str | None = None
    seconds: float | None = None


@dataclass
class Gathered:
    """One brief and everything it was written from."""

    rsid: str
    db: DbSnpRecord
    clinvar_retrieved: bool
    records: list[ClinVarRecord]
    has_significance: bool
    variants: list[GnomadVariant]
    papers: list[Paper]
    fallback_gene: str | None
    uniprot: UniProtHit | None
    calls: dict[str, SourceCall]
    sentences: list[str]
    text: str
    cache_dir: Path | None = None
    vep: VepAnnotation = field(default_factory=lambda: VepAnnotation(False))
    gwas: GwasResult = field(default_factory=lambda: GwasResult(False))
    citations: list[str] = field(default_factory=list)
    gene_summary: GeneSummary | None = None
    eqtls: list[Eqtl] = field(default_factory=list)
    ld_proxies: list[LdProxy] = field(default_factory=list)
    l2g: list[L2G] = field(default_factory=list)
    litvar: list[str] = field(default_factory=list)

    @property
    def pending(self) -> list[str]:
        return [key for key, call in self.calls.items() if call.origin == "pending"]


def build_brief(rsid: str, *, cache_dir: Path, refresh: bool, offline: bool) -> str:
    return gather(rsid, cache_dir=cache_dir, refresh=refresh, offline=offline).text


def gather(
    rsid: str,
    *,
    cache_dir: Path,
    refresh: bool,
    offline: bool,
    extras_cache: Path | None = None,
    budget: float | None = None,
) -> Gathered:
    """All providers at once. budget (seconds) caps the wait; None waits for every provider."""
    start = time.monotonic()
    deadline = None if budget is None else start + budget
    calls: dict[str, SourceCall] = {}

    def run(source: str, fetcher, *, dirs: tuple[Path, ...] = (cache_dir,), write: Path = cache_dir):
        return _cached_multi(rsid, source, fetcher, dirs=dirs, write=write, refresh=refresh, offline=offline)

    extra_dirs = (cache_dir,) + ((extras_cache,) if extras_cache and extras_cache != cache_dir else ())
    extra_write = extras_cache or cache_dir

    def dbsnp_task():
        return run("dbsnp", lambda: fetch_dbsnp(rsid))

    futures = {
        "dbsnp": _POOL.submit(dbsnp_task),
        "clinvar": _POOL.submit(run, "clinvar", lambda: fetch_clinvar(rsid)),
        "europepmc": _POOL.submit(run, "europepmc", lambda: fetch_europepmc(rsid)),
    }
    if extras_enabled():
        futures["gwas"] = _POOL.submit(
            run, "gwas", lambda: fetch_gwas(rsid, BY_KEY["gwas"].timeout), dirs=extra_dirs, write=extra_write
        )
        futures["ensembl_vep"] = _POOL.submit(
            run, "ensembl_vep", lambda: fetch_vep(rsid, BY_KEY["ensembl_vep"].timeout), dirs=extra_dirs, write=extra_write
        )
        futures["dbsnp_citations"] = _POOL.submit(
            run, "dbsnp_citations", lambda: fetch_citations(rsid), dirs=extra_dirs, write=extra_write
        )
        futures["litvar"] = _POOL.submit(run, "litvar", lambda: fetch_litvar(rsid), dirs=extra_dirs, write=extra_write)
        futures["ensembl_ld"] = _POOL.submit(run, "ensembl_ld", lambda: fetch_ld(rsid), dirs=extra_dirs, write=extra_write)

        def coordinate_ids(builder):
            db_payload, _ = db_future.result()
            record = parse_dbsnp(db_payload, rsid)
            if record.assembly != "GRCh38":
                return []
            return [builder(record.chrom, record.pos, ref, alt) for ref, alt in record.alleles]

        def gtex_task():
            return run("gtex", lambda: fetch_gtex(rsid, coordinate_ids(gtex_variant_id)), dirs=extra_dirs, write=extra_write)

        def ot_task():
            return run(
                "open_targets", lambda: fetch_open_targets(coordinate_ids(ot_variant_id)), dirs=extra_dirs, write=extra_write
            )

        futures["gtex"] = _POOL.submit(gtex_task)
        futures["open_targets"] = _POOL.submit(ot_task)
    db_future = futures["dbsnp"]

    def gnomad_task():
        db_payload, _ = db_future.result()
        return run("gnomad", lambda: fetch_gnomad(rsid, parse_dbsnp(db_payload, rsid)))

    futures["gnomad"] = _POOL.submit(gnomad_task)

    results: dict[str, dict | None] = {}
    _collect(futures, results, calls, deadline)

    db = parse_dbsnp(results.get("dbsnp"), rsid)
    clinvar_retrieved, records = parse_clinvar(results.get("clinvar"), rsid)
    _gnomad_retrieved, variants = parse_gnomad(results.get("gnomad"))
    _papers_retrieved, papers = parse_europepmc(results.get("europepmc"))
    vep = parse_vep(results.get("ensembl_vep")) if "ensembl_vep" in results else VepAnnotation(False)
    gwas = parse_gwas(results.get("gwas")) if "gwas" in results else GwasResult(False)
    citations = parse_citations(results.get("dbsnp_citations"))
    eqtls = parse_gtex(results.get("gtex"))
    ld_proxies = parse_ld(results.get("ensembl_ld"), rsid)
    l2g = parse_open_targets(results.get("open_targets"))
    litvar = parse_litvar(results.get("litvar"))

    used = [record for record in records if not record.no_change] or records
    has_significance = any(record.has_significance for record in used)
    uniprot = None
    gene = None
    if clinvar_retrieved and not has_significance:
        for record in records:
            if record.genes:
                gene = record.genes[0]
                break
        if gene is None and db.genes:
            gene = db.genes[0]
        if gene is None and vep.canonical_gene:
            gene = vep.canonical_gene
        if gene:
            target = gene
            uni = {"uniprot": _POOL.submit(run, "uniprot", lambda: fetch_uniprot(target))}
            _collect(uni, results, calls, deadline)
            uniprot_payload = results.get("uniprot")
            cached_gene = ((uniprot_payload or {}).get("request") or {}).get("gene")
            if uniprot_payload and cached_gene in {None, gene}:
                uniprot = parse_uniprot(uniprot_payload, gene)
    gene_summary = None
    if extras_enabled() and not has_significance:
        symbol = gene or (db.genes[0] if db.genes else None) or vep.canonical_gene or None
        if symbol:
            gs = {"ncbi_gene": _POOL.submit(
                run, "ncbi_gene", lambda: fetch_gene_summary(symbol), dirs=extra_dirs, write=extra_write
            )}
            _collect(gs, results, calls, deadline)
            gene_summary = parse_gene_summary(results.get("ncbi_gene"), symbol)
    return Gathered(
        rsid=rsid,
        db=db,
        clinvar_retrieved=clinvar_retrieved,
        records=records,
        has_significance=has_significance,
        variants=variants,
        papers=papers,
        fallback_gene=gene,
        uniprot=uniprot,
        calls=calls,
        sentences=brief_sentences(db, records, variants, papers, uniprot, clinvar_retrieved=clinvar_retrieved),
        text=assemble_brief(db, records, variants, papers, uniprot, clinvar_retrieved=clinvar_retrieved),
        cache_dir=cache_dir,
        vep=vep,
        gwas=gwas,
        citations=citations,
        gene_summary=gene_summary,
        eqtls=eqtls,
        ld_proxies=ld_proxies,
        l2g=l2g,
        litvar=litvar,
    )


def _collect(futures: dict, results: dict, calls: dict, deadline: float | None) -> None:
    """Wait for providers until the deadline; anything unfinished is marked pending."""
    waiting = {future: key for key, future in futures.items()}
    while waiting:
        timeout = None if deadline is None else max(0.0, deadline - time.monotonic())
        done, _ = wait(list(waiting), timeout=timeout, return_when=FIRST_COMPLETED)
        if not done:
            break
        for future in done:
            key = waiting.pop(future)
            try:
                payload, call = future.result()
            except Exception as exc:  # a provider bug never breaks the brief
                payload, call = None, SourceCall("failed", error=str(exc))
            results[key] = payload
            calls[key] = call
    for key in waiting.values():
        results[key] = None
        calls[key] = SourceCall("pending", error="still loading when the time budget ran out; refresh to see it")


def _cached_multi(
    rsid: str,
    source: str,
    fetcher,
    *,
    dirs: tuple[Path, ...],
    write: Path,
    refresh: bool,
    offline: bool,
) -> tuple[dict | None, SourceCall]:
    started = time.monotonic()
    if not refresh:
        for directory in dirs:
            path = cache_path(directory, rsid, source)
            if path.is_file():
                payload = load_cache(path)
                return payload, SourceCall("cache", _stamp(payload), seconds=time.monotonic() - started)
    if offline:
        print(f"{rsid} {source}: no cached response", file=sys.stderr)
        return None, SourceCall("no cache")
    try:
        payload = fetcher()
    except Exception as exc:
        print(f"{rsid} {source} failed: {exc}", file=sys.stderr)
        return None, SourceCall("failed", error=str(exc), seconds=time.monotonic() - started)
    with _LOCK:
        save_cache(cache_path(write, rsid, source), payload)
    return payload, SourceCall("live", _stamp(payload), seconds=time.monotonic() - started)


def _cached(cache_dir: Path, rsid: str, source: str, fetcher, *, refresh: bool, offline: bool, calls: dict) -> dict | None:
    """Single-source helper kept for callers outside gather."""
    payload, call = _cached_multi(rsid, source, fetcher, dirs=(cache_dir,), write=cache_dir, refresh=refresh, offline=offline)
    calls[source] = call
    return payload


def _stamp(payload: dict | None) -> str | None:
    return (payload or {}).get("retrieved")


def default_user_cache() -> Path:
    """Where live responses for SNPs outside the checked-in cache are saved."""
    import os

    configured = os.environ.get("SNP_BRIEF_USER_CACHE")
    return Path(configured) if configured else Path.home() / ".cache" / "snp-brief"


def load_variant(
    rsid: str,
    *,
    repo_cache: Path,
    user_cache: Path,
    offline: bool = False,
    refresh: bool = False,
    budget: float | None = None,
) -> Gathered:
    """Cache first, live only for what is missing.

    - A SNP in the checked-in cache (repo_cache) is read from there.
    - Any other SNP is read from user_cache; sources missing there are fetched
      live and saved, so the next lookup is instant.
    - refresh re-fetches every source into user_cache and leaves the checked-in
      cache unchanged.
    - offline never calls the network; missing sources are reported as missing.
    """
    if refresh and not offline:
        return gather(rsid, cache_dir=user_cache, refresh=True, offline=False, budget=budget)
    # A SNP fetched or refreshed into user_cache is newer than the checked-in copy.
    # (Extras-only files there, from a checked-in SNP, do not count.)
    if (user_cache / rsid / "dbsnp.json").is_file():
        return gather(rsid, cache_dir=user_cache, refresh=False, offline=offline, budget=budget)
    if (repo_cache / rsid).is_dir():
        # Checked-in SNP: new providers missing from cache/ are fetched into the user cache.
        return gather(rsid, cache_dir=repo_cache, refresh=False, offline=offline, extras_cache=user_cache, budget=budget)
    return gather(rsid, cache_dir=user_cache, refresh=False, offline=offline, budget=budget)
