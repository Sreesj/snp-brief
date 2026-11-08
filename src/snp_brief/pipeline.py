"""Fetch each source in order, then write the brief. Cache hits skip the network."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

from snp_brief.brief import assemble_brief, brief_sentences
from snp_brief.cache_store import cache_path, load_cache, save_cache
from snp_brief.models import ClinVarRecord, DbSnpRecord, GnomadVariant, Paper, UniProtHit
from snp_brief.sources.clinvar import fetch_clinvar, parse_clinvar
from snp_brief.sources.dbsnp import fetch_dbsnp, parse_dbsnp
from snp_brief.sources.europepmc import fetch_europepmc, parse_europepmc
from snp_brief.sources.gnomad import fetch_gnomad, parse_gnomad
from snp_brief.sources.uniprot import fetch_uniprot, parse_uniprot


@dataclass
class SourceCall:
    """How one source answered: origin is "cache", "live", "no cache", or "failed"."""

    origin: str
    retrieved: str | None = None
    error: str | None = None


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


def build_brief(rsid: str, *, cache_dir: Path, refresh: bool, offline: bool) -> str:
    return gather(rsid, cache_dir=cache_dir, refresh=refresh, offline=offline).text


def gather(rsid: str, *, cache_dir: Path, refresh: bool, offline: bool) -> Gathered:
    calls: dict[str, SourceCall] = {}
    dbsnp_payload = _cached(
        cache_dir, rsid, "dbsnp", lambda: fetch_dbsnp(rsid), refresh=refresh, offline=offline, calls=calls
    )
    db = parse_dbsnp(dbsnp_payload, rsid)
    clinvar_payload = _cached(
        cache_dir, rsid, "clinvar", lambda: fetch_clinvar(rsid), refresh=refresh, offline=offline, calls=calls
    )
    clinvar_retrieved, records = parse_clinvar(clinvar_payload, rsid)
    gnomad_payload = _cached(
        cache_dir,
        rsid,
        "gnomad",
        lambda: fetch_gnomad(rsid, db),
        refresh=refresh,
        offline=offline,
        calls=calls,
    )
    _gnomad_retrieved, variants = parse_gnomad(gnomad_payload)
    paper_payload = _cached(
        cache_dir, rsid, "europepmc", lambda: fetch_europepmc(rsid), refresh=refresh, offline=offline, calls=calls
    )
    _papers_retrieved, papers = parse_europepmc(paper_payload)
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
        if gene:
            uniprot_payload = _cached(
                cache_dir,
                rsid,
                "uniprot",
                lambda: fetch_uniprot(gene),
                refresh=refresh,
                offline=offline,
                calls=calls,
            )
            cached_gene = ((uniprot_payload or {}).get("request") or {}).get("gene")
            if cached_gene in {None, gene}:
                uniprot = parse_uniprot(uniprot_payload, gene)
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
    )


def _cached(
    cache_dir: Path,
    rsid: str,
    source: str,
    fetcher,
    *,
    refresh: bool,
    offline: bool,
    calls: dict[str, SourceCall],
) -> dict | None:
    path = cache_path(cache_dir, rsid, source)
    if path.is_file() and not refresh:
        payload = load_cache(path)
        calls[source] = SourceCall("cache", _stamp(payload))
        return payload
    if offline:
        print(f"{rsid} {source}: no cached response", file=sys.stderr)
        calls[source] = SourceCall("no cache")
        return None
    try:
        payload = fetcher()
    except Exception as exc:
        print(f"{rsid} {source} failed: {exc}", file=sys.stderr)
        calls[source] = SourceCall("failed", error=str(exc))
        return None
    save_cache(path, payload)
    calls[source] = SourceCall("live", _stamp(payload))
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
        return gather(rsid, cache_dir=user_cache, refresh=True, offline=False)
    # A SNP fetched or refreshed into user_cache is newer than the checked-in copy.
    if any((user_cache / rsid).glob("*.json")):
        return gather(rsid, cache_dir=user_cache, refresh=False, offline=offline)
    if (repo_cache / rsid).is_dir():
        return gather(rsid, cache_dir=repo_cache, refresh=False, offline=offline)
    return gather(rsid, cache_dir=user_cache, refresh=False, offline=offline)
