"""ClinVar: gene, molecular consequence, and clinical significance for an rsID.

Records are kept only when their dbSNP cross-reference is this rsID. A
multi-variant haplotype (an HGVS allele set) is not this SNP and is dropped.
"""

from __future__ import annotations

import re
from datetime import datetime, timezone

from snp_brief.http_client import request_json
from snp_brief.models import Classification, ClinVarRecord

_TITLE_GENE = re.compile(r"\(([A-Za-z][A-Za-z0-9-]*)\):[cng]\.")
_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_SKIP_TRAITS = {"not provided", "not specified", "none", "see cases", ""}


def fetch_clinvar(rsid: str) -> dict:
    esearch = request_json(
        _ESEARCH,
        params={
            "db": "clinvar",
            "term": rsid,
            "retmode": "json",
            "retmax": "50",
            "tool": "snp_brief",
        },
    )
    result = esearch.get("esearchresult") or {}
    ids = list(result.get("idlist") or [])
    try:
        count = int(result.get("count") or 0)
    except ValueError:
        count = len(ids)
    if count > len(ids) and count <= 100:
        esearch = request_json(
            _ESEARCH,
            params={
                "db": "clinvar",
                "term": rsid,
                "retmode": "json",
                "retmax": str(count),
                "tool": "snp_brief",
            },
        )
        ids = list((esearch.get("esearchresult") or {}).get("idlist") or [])
    esummary = None
    if ids:
        esummary = request_json(
            _ESUMMARY,
            params={
                "db": "clinvar",
                "id": ",".join(ids),
                "retmode": "json",
                "tool": "snp_brief",
            },
        )
    return {
        "source": "clinvar",
        "retrieved": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "request": {"term": rsid},
        "response": {"esearch": esearch, "esummary": esummary},
    }


def parse_clinvar(payload: dict | None, rsid: str) -> tuple[bool, list[ClinVarRecord]]:
    if not payload:
        return False, []
    numeric = rsid[2:]
    esummary = (payload.get("response") or {}).get("esummary") or {}
    result = esummary.get("result") or {}
    records: list[ClinVarRecord] = []
    for uid in result.get("uids") or []:
        doc = result.get(str(uid))
        if not isinstance(doc, dict):
            continue
        if numeric not in _dbsnp_ids(doc):
            continue
        cdna = _cdna_change(_cdna(doc))
        title = (doc.get("title") or "").strip()
        if _haplotype(title, cdna):
            continue
        consequences = [str(item).strip() for item in doc.get("molecular_consequence_list") or [] if str(item).strip()]
        genes = _genes(doc, title)
        classifications = _classifications(doc)
        no_change = _no_change(consequences, cdna)
        records.append(
            ClinVarRecord(
                accession=(doc.get("accession") or f"uid{uid}").strip(),
                title=title,
                genes=genes,
                consequences=consequences,
                cdna=cdna,
                classifications=classifications,
                no_change=no_change,
            )
        )
    return True, records


def _dbsnp_ids(doc: dict) -> set[str]:
    found: set[str] = set()
    for variation in doc.get("variation_set") or []:
        for xref in variation.get("variation_xrefs") or []:
            if str(xref.get("db_source") or "").lower() == "dbsnp":
                found.add(str(xref.get("db_id")))
    return found


def _cdna(doc: dict) -> str:
    for variation in doc.get("variation_set") or []:
        cdna = (variation.get("cdna_change") or "").strip()
        if cdna:
            return cdna
    return ""


def _cdna_change(cdna: str) -> str:
    """Keep the c. or n. change when ClinVar prefixes a transcript."""
    for token in ("c.", "n."):
        index = cdna.find(token)
        if index >= 0:
            return cdna[index:].split()[0]
    return cdna


def _genes(doc: dict, title: str) -> list[str]:
    """Prefer ClinVar's gene sort. The symbol in an HGVS title can be the transcript's gene."""
    sort = (doc.get("gene_sort") or "").strip()
    if sort and not sort.startswith("LOC"):
        return [sort]
    match = _TITLE_GENE.search(title)
    if match and not match.group(1).startswith("LOC"):
        return [match.group(1)]
    listed: list[str] = []
    for gene in doc.get("genes") or []:
        symbol = (gene.get("symbol") or "").strip()
        if symbol and symbol not in listed:
            listed.append(symbol)
    non_loc = [symbol for symbol in listed if not symbol.startswith("LOC")]
    if non_loc:
        return non_loc
    if sort:
        return [sort]
    return listed[:1]


def _haplotype(title: str, cdna: str) -> bool:
    blob = f"{title} {cdna}"
    return ";" in blob or "[" in blob


def _no_change(consequences: list[str], cdna: str) -> bool:
    folded = [item.lower() for item in consequences]
    if folded and all(item == "no sequence alteration" for item in folded):
        return True
    if ">" not in cdna and "=" in cdna:
        return True
    return False


def _classifications(doc: dict) -> list[Classification]:
    found: list[Classification] = []
    for key in ("germline_classification", "oncogenicity_classification", "clinical_impact_classification"):
        block = doc.get(key) or {}
        if not isinstance(block, dict):
            continue
        description = (block.get("description") or "").strip()
        if not description:
            continue
        review = (block.get("review_status") or "").strip()
        found.append(Classification(description, review, _traits(block)))
    return found


def _traits(block: dict) -> list[str]:
    names: list[str] = []
    for trait in block.get("trait_set") or []:
        name = (trait.get("trait_name") or "").strip()
        if name.lower() in _SKIP_TRAITS or name in names:
            continue
        names.append(name)
    return names
