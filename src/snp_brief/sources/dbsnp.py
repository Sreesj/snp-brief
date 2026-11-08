"""dbSNP: rsID, GRCh38 coordinate when the RefSeq accession is current, and alleles."""

from __future__ import annotations

from collections import OrderedDict
from datetime import datetime, timezone

from snp_brief.http_client import request_json
from snp_brief.models import DbSnpRecord

# Current chromosome accessions in GRCh38.p14. Used only to name the assembly.
_GRCH38 = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.10": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}

_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def fetch_dbsnp(rsid: str) -> dict:
    numeric = rsid[2:]
    body = request_json(
        _ESUMMARY,
        params={"db": "snp", "id": numeric, "retmode": "json", "tool": "snp_brief"},
    )
    return {
        "source": "dbsnp",
        "retrieved": _now(),
        "request": {"url": _ESUMMARY, "id": numeric},
        "response": body,
    }


def parse_dbsnp(payload: dict | None, rsid: str) -> DbSnpRecord:
    if not payload:
        return DbSnpRecord(retrieved=False, rsid=rsid)
    result = (payload.get("response") or {}).get("result") or {}
    uids = result.get("uids") or []
    doc = result.get(str(uids[0])) if uids else None
    if not isinstance(doc, dict) or "error" in doc and "chrpos" not in doc:
        return DbSnpRecord(retrieved=True, rsid=rsid)
    accession = (doc.get("acc") or "").strip() or None
    chrom, pos = _chrpos(doc.get("chrpos") or "")
    assembly = "GRCh38" if accession in _GRCH38 else None
    alleles = _alleles(doc.get("spdi") or "", accession, pos)
    genes = []
    for gene in doc.get("genes") or []:
        name = (gene.get("name") or "").strip()
        if name and name not in genes:
            genes.append(name)
    return DbSnpRecord(
        retrieved=True,
        rsid=rsid,
        assembly=assembly,
        accession=accession,
        chrom=chrom,
        pos=pos,
        alleles=alleles,
        genes=genes,
        coding_by_alt=coding_by_alt(doc.get("docsum") or "", accession, [alt for _, alt in alleles]),
        fxn_class=[part.strip() for part in str(doc.get("fxn_class") or "").split(",") if part.strip()],
        current_rsid=_current(doc, rsid),
    )


def _current(doc: dict, rsid: str) -> str:
    """dbSNP answers a merged rsID with the surviving record; its snp_id is the current rsID."""
    snp_id = str(doc.get("snp_id") or "").strip()
    current = f"rs{snp_id}" if snp_id.isdigit() else ""
    return current if current and current != rsid else ""


def coding_by_alt(docsum: str, accession: str | None, alts: list[str]) -> dict[str, str]:
    """Pair each genomic alt with the coding change listed in the same order.

    dbSNP prints HGVS as genomic alleles for an accession, then the coding
    changes for each transcript in that same allele order. The map is kept
    only when those two lists have the same length.
    """
    if not docsum or not accession or not alts:
        return {}
    hgvs = ""
    for part in docsum.split("|"):
        if part.startswith("HGVS="):
            hgvs = part[len("HGVS=") :]
            break
    if not hgvs:
        return {}
    tokens = [token for token in hgvs.split(",") if token]
    prefix = f"{accession}:g."
    genomic: list[str] = []
    for token in tokens:
        if not token.startswith(prefix):
            continue
        change = token.split(":g.", 1)[1]
        if ">" not in change:
            continue
        genomic.append(change.split(">")[-1])
    if genomic != alts:
        return {}
    by_transcript: OrderedDict[str, list[str]] = OrderedDict()
    for token in tokens:
        if ":c." not in token and ":n." not in token:
            continue
        transcript, change = token.split(":", 1)
        if not transcript.startswith(("NM_", "NR_")):
            continue
        by_transcript.setdefault(transcript, []).append(change)
    for transcript, changes in by_transcript.items():
        if transcript.startswith("NM_") and len(changes) == len(genomic):
            return dict(zip(genomic, changes, strict=True))
    return {}


def _chrpos(value: str) -> tuple[str | None, int | None]:
    if ":" not in value:
        return None, None
    chrom, _, pos = value.partition(":")
    try:
        return chrom.strip() or None, int(pos)
    except ValueError:
        return None, None


def _alleles(spdi: str, accession: str | None, pos: int | None) -> list[tuple[str, str]]:
    if not spdi or not accession or pos is None:
        return []
    found: list[tuple[str, str]] = []
    for token in spdi.split(","):
        parts = token.split(":")
        if len(parts) != 4 or parts[0] != accession:
            continue
        try:
            one_based = int(parts[1]) + 1
        except ValueError:
            continue
        # SNVs sit exactly at chrpos. Insertions, duplications and deletions in repeats are often
        # placed a few bases away by SPDI normalisation, so allow a small window for those.
        deleted, inserted = parts[2], parts[3]
        is_snv = len(deleted) == 1 and len(inserted) == 1
        if (is_snv and one_based != pos) or (not is_snv and abs(one_based - pos) > 50):
            continue
        pair = (deleted, inserted)
        if deleted and inserted and deleted != inserted and pair not in found:
            found.append(pair)
    return found


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()
