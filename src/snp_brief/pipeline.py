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
from snp_brief.sources import extra
from snp_brief.sources.extra import BorzoiResult, Constraint, PgxResult, PredictionScores, RegulomeResult, SpliceResult, Timeline
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
    orient_risk_alleles,
    source_enabled,
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
# Fetches in flight, keyed by (rsid, source). A second request for the same source waits for the
# first instead of starting another one, so "Check again" picks up a slow answer when it lands.
_INFLIGHT: dict[tuple[str, str], "_Flight"] = {}
# Sources that just failed, keyed by (rsid, source) -> (monotonic time, SourceCall). Picking up
# late answers reads this instead of trying the same broken server again; Refresh data ignores it.
_RECENT_FAIL: dict[tuple[str, str], tuple[float, "SourceCall"]] = {}
FAIL_MEMORY_SECONDS = 600


class _Flight:
    def __init__(self) -> None:
        self.done = threading.Event()
        self.result: tuple[dict | None, "SourceCall"] | None = None


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
    merged_from: str = ""  # the rsID that was asked for, when dbSNP merged it into rsid
    eqtls: list[Eqtl] = field(default_factory=list)
    ld_proxies: list[LdProxy] = field(default_factory=list)
    l2g: list[L2G] = field(default_factory=list)
    litvar: list[str] = field(default_factory=list)
    findings: list = field(default_factory=list)
    scores: PredictionScores = field(default_factory=lambda: PredictionScores(False))
    pgx: PgxResult = field(default_factory=lambda: PgxResult(False))
    splice: SpliceResult = field(default_factory=lambda: SpliceResult(False))
    regulome: RegulomeResult = field(default_factory=lambda: RegulomeResult(False))
    constraint: Constraint = field(default_factory=lambda: Constraint(False))
    timeline: Timeline = field(default_factory=lambda: Timeline(False))
    borzoi: BorzoiResult = field(default_factory=lambda: BorzoiResult(False))

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
    merged_from: str = "",
) -> Gathered:
    """All providers at once. budget (seconds) caps the wait; None waits for every provider.

    If dbSNP reports that the rsID was merged into another, everything is gathered again for the
    current rsID (once), and the result records the original ID in merged_from.
    """
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
    # Bound before any task that waits on it is submitted (GTEx, Open Targets, gnomAD).
    db_future = futures["dbsnp"]
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
        for key, fetcher in (
            ("myvariant", lambda: extra.fetch_myvariant(rsid)),
            ("pharmgkb", lambda: extra.fetch_pharmgkb(rsid)),
            ("timeline", lambda: extra.fetch_timeline(rsid)),
        ):
            futures[key] = _POOL.submit(run, key, fetcher, dirs=extra_dirs, write=extra_write)

        def splice_task():
            return run(
                "spliceai", lambda: extra.fetch_spliceai(coordinate_ids(extra.splice_variant_ids)), dirs=extra_dirs, write=extra_write
            )

        def constraint_task():
            db_payload, _ = db_future.result()
            genes = parse_dbsnp(db_payload, rsid).genes
            if not genes:
                return None, SourceCall("no data", error="no gene at this position")
            return run("constraint", lambda: extra.fetch_constraint(genes[0]), dirs=extra_dirs, write=extra_write)

        def borzoi_task():
            # Precomputed on a GPU by brev/score_borzoi.py; never fetched over the network.
            for directory in extra_dirs:
                path = cache_path(directory, rsid, "borzoi")
                if path.is_file():
                    payload = load_cache(path)
                    return payload, SourceCall("cache", _stamp(payload))
            return None, SourceCall("no data", error="not computed for this variant yet")

        futures["borzoi"] = _POOL.submit(borzoi_task)
        if source_enabled("spliceai"):
            futures["spliceai"] = _POOL.submit(splice_task)
        if source_enabled("regulome"):
            futures["regulome"] = _POOL.submit(run, "regulome", lambda: extra.fetch_regulome(rsid), dirs=extra_dirs, write=extra_write)
        futures["constraint"] = _POOL.submit(constraint_task)
    def gnomad_task():
        db_payload, _ = db_future.result()
        return run("gnomad", lambda: fetch_gnomad(rsid, parse_dbsnp(db_payload, rsid)))

    futures["gnomad"] = _POOL.submit(gnomad_task)

    results: dict[str, dict | None] = {}
    _collect(futures, results, calls, deadline)

    db = parse_dbsnp(results.get("dbsnp"), rsid)
    if db.current_rsid and not merged_from:
        remaining = None if deadline is None else max(1.0, deadline - time.monotonic())
        return gather(
            db.current_rsid,
            cache_dir=cache_dir,
            refresh=refresh,
            offline=offline,
            extras_cache=extras_cache,
            budget=remaining,
            merged_from=rsid,
        )
    clinvar_retrieved, records = parse_clinvar(results.get("clinvar"), rsid)
    _gnomad_retrieved, variants = parse_gnomad(results.get("gnomad"))
    _papers_retrieved, papers = parse_europepmc(results.get("europepmc"))
    vep = parse_vep(results.get("ensembl_vep")) if "ensembl_vep" in results else VepAnnotation(False)
    gwas = parse_gwas(results.get("gwas"), rsid) if "gwas" in results else GwasResult(False)
    orient_risk_alleles(gwas.hits, db.alleles)
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
        findings=_findings(rsid, cache_dir, extras_cache),
        scores=extra.parse_myvariant(results.get("myvariant")) if "myvariant" in results else PredictionScores(False),
        pgx=extra.parse_pharmgkb(results.get("pharmgkb")) if "pharmgkb" in results else PgxResult(False),
        splice=extra.parse_spliceai(results.get("spliceai")) if "spliceai" in results else SpliceResult(False),
        regulome=extra.parse_regulome(results.get("regulome")) if "regulome" in results else RegulomeResult(False),
        constraint=extra.parse_constraint(results.get("constraint")) if "constraint" in results else Constraint(False),
        timeline=extra.parse_timeline(results.get("timeline")) if "timeline" in results else Timeline(False),
        borzoi=extra.parse_borzoi(results.get("borzoi")) if "borzoi" in results else BorzoiResult(False),
        merged_from=merged_from,
    )


def _findings(rsid: str, cache_dir: Path, extras_cache: Path | None) -> list:
    from snp_brief.literature import cached_findings

    return cached_findings([extras_cache, cache_dir], rsid)


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
                if not offline and _stale(directory, payload):
                    break  # older than SNP_BRIEF_MAX_AGE_DAYS in the user cache: fetch again
                return payload, SourceCall("cache", _stamp(payload), seconds=time.monotonic() - started)
    if offline:
        print(f"{rsid} {source}: no cached response", file=sys.stderr)
        return None, SourceCall("no cache")
    key = (rsid, source)
    with _LOCK:
        failed = _RECENT_FAIL.get(key)
        if failed and not refresh and time.monotonic() - failed[0] < FAIL_MEMORY_SECONDS:
            return None, failed[1]
        flight = _INFLIGHT.get(key)
        mine = flight is None
        if mine:
            flight = _INFLIGHT[key] = _Flight()
    if not mine:
        flight.done.wait(180)
        return flight.result or (None, SourceCall("failed", error="no answer"))
    try:
        flight.result = _fetch_and_save(rsid, source, fetcher, write, started)
        with _LOCK:
            if flight.result[1].origin == "failed":
                _RECENT_FAIL[key] = (time.monotonic(), flight.result[1])
            else:
                _RECENT_FAIL.pop(key, None)
        return flight.result
    finally:
        with _LOCK:
            _INFLIGHT.pop(key, None)
        flight.done.set()


def _fetch_and_save(rsid: str, source: str, fetcher, write: Path, started: float) -> tuple[dict | None, "SourceCall"]:
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


def _stale(directory: Path, payload: dict | None) -> bool:
    """Only the user cache ages out; the checked-in cache/ is a fixed demo snapshot."""
    import os
    from datetime import datetime, timezone

    try:
        if Path(directory).resolve() != default_user_cache().resolve():
            return False
        days = float(os.environ.get("SNP_BRIEF_MAX_AGE_DAYS", "30"))
        stamp = _stamp(payload)
        if not stamp or days <= 0:
            return False
        age = datetime.now(timezone.utc) - datetime.fromisoformat(stamp)
        return age.total_seconds() > days * 86400
    except (ValueError, OSError):
        return False


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


def in_flight(rsid: str, sources) -> bool:
    """True while any of these sources for this rsID is still being fetched in the background."""
    with _LOCK:
        return any((rsid, source) in _INFLIGHT for source in sources)
