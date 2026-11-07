"""Fetch each source in order, then write the brief. Cache hits skip the network."""

from __future__ import annotations

import sys
from pathlib import Path

from snp_brief.brief import assemble_brief
from snp_brief.cache_store import cache_path, load_cache, save_cache
from snp_brief.sources.clinvar import fetch_clinvar, parse_clinvar
from snp_brief.sources.dbsnp import fetch_dbsnp, parse_dbsnp
from snp_brief.sources.europepmc import fetch_europepmc, parse_europepmc
from snp_brief.sources.gnomad import fetch_gnomad, parse_gnomad
from snp_brief.sources.uniprot import fetch_uniprot, parse_uniprot


def build_brief(rsid: str, *, cache_dir: Path, refresh: bool, offline: bool) -> str:
    dbsnp_payload = _cached(cache_dir, rsid, "dbsnp", lambda: fetch_dbsnp(rsid), refresh=refresh, offline=offline)
    db = parse_dbsnp(dbsnp_payload, rsid)
    clinvar_payload = _cached(
        cache_dir, rsid, "clinvar", lambda: fetch_clinvar(rsid), refresh=refresh, offline=offline
    )
    clinvar_retrieved, records = parse_clinvar(clinvar_payload, rsid)
    gnomad_payload = _cached(
        cache_dir,
        rsid,
        "gnomad",
        lambda: fetch_gnomad(rsid, db),
        refresh=refresh,
        offline=offline,
    )
    _gnomad_retrieved, variants = parse_gnomad(gnomad_payload)
    paper_payload = _cached(
        cache_dir, rsid, "europepmc", lambda: fetch_europepmc(rsid), refresh=refresh, offline=offline
    )
    _papers_retrieved, papers = parse_europepmc(paper_payload)
    used = [record for record in records if not record.no_change] or records
    has_significance = any(record.has_significance for record in used)
    uniprot = None
    if clinvar_retrieved and not has_significance:
        gene = None
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
            )
            cached_gene = ((uniprot_payload or {}).get("request") or {}).get("gene")
            if cached_gene in {None, gene}:
                uniprot = parse_uniprot(uniprot_payload, gene)
    return assemble_brief(
        db,
        records,
        variants,
        papers,
        uniprot,
        clinvar_retrieved=clinvar_retrieved,
    )


def _cached(cache_dir: Path, rsid: str, source: str, fetcher, *, refresh: bool, offline: bool) -> dict | None:
    path = cache_path(cache_dir, rsid, source)
    if path.is_file() and not refresh:
        return load_cache(path)
    if offline:
        print(f"{rsid} {source}: no cached response", file=sys.stderr)
        return None
    try:
        payload = fetcher()
    except Exception as exc:
        print(f"{rsid} {source} failed: {exc}", file=sys.stderr)
        return None
    save_cache(path, payload)
    return payload
