"""ClinVar: gene, molecular consequence, and clinical significance for an rsID.

Records are kept only when their dbSNP cross-reference is this rsID. A
multi-variant haplotype (an HGVS allele set) is not this SNP and is dropped.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

from snp_brief.http_client import request_json, request_text
from snp_brief.models import Classification, ClinVarRecord

_TITLE_GENE = re.compile(r"\(([A-Za-z][A-Za-z0-9-]*)\):[cng]\.")
_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
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
    submissions = fetch_submission_map(esummary, rsid) if esummary else {}
    response = {"esearch": esearch, "esummary": esummary}
    if submissions:
        response["submissions"] = submissions
    return {
        "source": "clinvar",
        "retrieved": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "request": {"term": rsid},
        "response": response,
    }


def parse_clinvar(payload: dict | None, rsid: str) -> tuple[bool, list[ClinVarRecord]]:
    if not payload:
        return False, []
    numeric = rsid[2:]
    response = payload.get("response") or {}
    esummary = response.get("esummary") or {}
    submission_map = response.get("submissions") or {}
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
                submissions=_classifications_from_submissions(submission_map.get(str(uid)) or []),
                aliases=_aliases(doc),
            )
        )
    return True, records


def _aliases(doc: dict) -> list[str]:
    found: list[str] = []
    for item in doc.get("variation_set") or []:
        for alias in item.get("aliases") or []:
            text = str(alias).strip()
            if text and text not in found:
                found.append(text)
    return found[:5]


def fetch_submission_map(esummary: dict | None, rsid: str) -> dict[str, list[dict]]:
    """Submitted interpretations for each VCV that is this rsID.

    esummary only carries ClinVar's aggregate call. The VCV XML lists each
    current SCV. Haplotypes are not this SNP and are not fetched.
    """
    if not esummary:
        return {}
    numeric = rsid[2:] if rsid.lower().startswith("rs") else rsid
    result = esummary.get("result") or {}
    found: dict[str, list[dict]] = {}
    for uid in result.get("uids") or []:
        doc = result.get(str(uid))
        if not isinstance(doc, dict):
            continue
        if numeric not in _dbsnp_ids(doc):
            continue
        title = (doc.get("title") or "").strip()
        cdna = _cdna_change(_cdna(doc))
        if _haplotype(title, cdna):
            continue
        accession = (doc.get("accession") or "").strip()
        if not accession:
            continue
        xml_text = request_text(
            _EFETCH,
            params={
                "db": "clinvar",
                "rettype": "vcv",
                "id": accession,
                "tool": "snp_brief",
            },
            timeout=120,
        )
        found[str(uid)] = parse_submission_xml(xml_text)
    return found


def parse_submission_xml(xml_text: str) -> list[dict]:
    """Current SCV germline interpretations. A removed submission is not a call."""
    root = ET.fromstring(xml_text)
    found: list[dict] = []
    for assertion in root.iter("ClinicalAssertion"):
        status = (assertion.findtext("RecordStatus") or "").strip().lower()
        if status and status != "current":
            continue
        acc_el = assertion.find("ClinVarAccession")
        accession = (acc_el.attrib.get("Accession") if acc_el is not None else "") or ""
        accession = accession.strip()
        description = (assertion.findtext("Classification/GermlineClassification") or "").strip()
        review = (assertion.findtext("Classification/ReviewStatus") or "").strip()
        if not accession or not description:
            continue
        found.append(
            {
                "accession": accession,
                "description": description,
                "review_status": review,
                "traits": _assertion_traits(assertion),
            }
        )
    return found


def _assertion_traits(assertion: ET.Element) -> list[str]:
    """Preferred condition on the submission itself, not an observed phenotype."""
    names: list[str] = []
    trait_set = assertion.find("TraitSet")
    if trait_set is None:
        return names
    for trait in trait_set.findall("Trait"):
        values = trait.findall("Name/ElementValue")
        chosen = _trait_name(values, "Preferred") or _trait_name(values, "included")
        if chosen and chosen not in names:
            names.append(chosen)
    return names


def _trait_name(values: list[ET.Element], kind: str) -> str:
    for value in values:
        if value.attrib.get("Type") == kind and (value.text or "").strip():
            return value.text.strip()
    return ""


def _classifications_from_submissions(items: list) -> list[Classification]:
    found: list[Classification] = []
    if not isinstance(items, list):
        return found
    for item in items:
        if not isinstance(item, dict):
            continue
        description = str(item.get("description") or "").strip()
        review = str(item.get("review_status") or "").strip()
        accession = str(item.get("accession") or "").strip()
        if not description or not accession:
            continue
        traits: list[str] = []
        for trait in item.get("traits") or []:
            name = str(trait).strip()
            if name and name not in traits:
                traits.append(name)
        found.append(Classification(description, review, traits, accession))
    return found


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
