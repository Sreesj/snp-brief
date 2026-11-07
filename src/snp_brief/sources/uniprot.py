"""Reviewed UniProt function text. Used only when ClinVar has no clinical significance.

The disease comment is never read. UniProt describes the protein, not this SNP.
"""

from __future__ import annotations

from datetime import datetime, timezone

from snp_brief.http_client import request_json
from snp_brief.models import UniProtHit

_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot(gene: str) -> dict:
    query = f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)"
    body = request_json(
        _SEARCH,
        params={
            "query": query,
            "fields": "accession,gene_names,cc_function,reviewed",
            "format": "json",
            "size": "5",
        },
    )
    return {
        "source": "uniprot",
        "retrieved": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "request": {"gene": gene, "query": query},
        "response": body,
    }


def parse_uniprot(payload: dict | None, gene: str) -> UniProtHit | None:
    if not payload or not gene:
        return None
    for entry in (payload.get("response") or {}).get("results") or []:
        if not _reviewed(entry):
            continue
        if not _names(entry, gene):
            continue
        text = _function_text(entry)
        accession = (entry.get("primaryAccession") or "").strip()
        if text and accession:
            return UniProtHit(accession, _short(text))
    return None


def _reviewed(entry: dict) -> bool:
    entry_type = str(entry.get("entryType") or "")
    if entry_type:
        return "reviewed" in entry_type.lower() and "unreviewed" not in entry_type.lower()
    return True


def _names(entry: dict, gene: str) -> bool:
    wanted = gene.lower()
    found = False
    for gene_block in entry.get("genes") or []:
        primary = ((gene_block.get("geneName") or {}).get("value") or "").lower()
        if primary == wanted:
            return True
        for synonym in gene_block.get("synonyms") or []:
            if (synonym.get("value") or "").lower() == wanted:
                found = True
    # A reviewed gene_exact search can omit the gene block. Trust it only then.
    return found or not entry.get("genes")


def _function_text(entry: dict) -> str:
    for comment in entry.get("comments") or []:
        if comment.get("commentType") != "FUNCTION":
            continue
        for text in comment.get("texts") or []:
            value = " ".join(str(text.get("value") or "").split())
            if value:
                return value
    return ""


def _short(text: str, limit: int = 480) -> str:
    if len(text) <= limit:
        return text
    cut = text[:limit]
    dot = cut.rfind(". ")
    if dot > 80:
        return cut[: dot + 1]
    return cut.rstrip() + "…"
