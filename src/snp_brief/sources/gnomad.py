"""gnomAD allele counts. Population figures are the API counts, not estimates.

Comparison rule, applied separately to each alternate allele: two frequencies
disagree when the higher is at least twice the lower, or the absolute gap is
at least 0.1, or one is 0 and the other is at least 0.001. Ancestry groups
with fewer than 1000 alleles are left out of that comparison. Exome and genome
call sets use the same rule on the allele frequency the API reports. The rule
decides whether to keep one figure or both. It does not create a new figure.
"""

from __future__ import annotations

import sys
import time
from datetime import datetime, timezone

from snp_brief.http_client import request_json
from snp_brief.models import CallSet, DbSnpRecord, GnomadVariant, PopCount

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
MIN_AN = 1000
FOLD = 2.0
ABS_GAP = 0.1
ZERO_FLOOR = 0.001

# gnomAD genetic-ancestry group ids and the labels gnomAD publishes for them.
ANCESTRY = {
    "afr": "African/African American",
    "ami": "Amish",
    "amr": "Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "fin": "European (Finnish)",
    "mid": "Middle Eastern",
    "nfe": "European (non-Finnish)",
    "sas": "South Asian",
    "remaining": "Remaining",
}

_DATASETS = ("gnomad_r4", "gnomad_r3")

_TOTALS_RSID = """
query ($rsid: String!, $dataset: DatasetId!) {
  variant(rsid: $rsid, dataset: $dataset) {
    variantId
    reference_genome
    chrom
    pos
    ref
    alt
    exome { ac an af }
    genome { ac an af }
  }
}
"""

_TOTALS_ID = """
query ($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    reference_genome
    chrom
    pos
    ref
    alt
    exome { ac an af }
    genome { ac an af }
  }
}
"""

_POPS = {
    "exome": """
query ($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    exome { populations { id ac an } }
  }
}
""",
    "genome": """
query ($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variantId
    genome { populations { id ac an } }
  }
}
""",
}


def fetch_gnomad(rsid: str, db: DbSnpRecord) -> dict:
    for dataset in _DATASETS:
        found = _lookup_variants(rsid, db, dataset)
        if not found:
            continue
        packaged = []
        for variant in found:
            packaged.append(
                {
                    "totals": {"data": {"variant": variant}},
                    "populations": _populations(variant, dataset),
                }
            )
        return {
            "source": "gnomad",
            "retrieved": _now(),
            "request": {"rsid": rsid, "dataset": dataset},
            "response": {"dataset": dataset, "variants": packaged},
        }
    return {
        "source": "gnomad",
        "retrieved": _now(),
        "request": {"rsid": rsid},
        "response": {"dataset": None, "variants": []},
    }


def parse_gnomad(payload: dict | None) -> tuple[bool, list[GnomadVariant]]:
    if not payload:
        return False, []
    response = payload.get("response") or {}
    label = _dataset_label(response.get("dataset"))
    variants: list[GnomadVariant] = []
    if "variants" in response:
        for item in response.get("variants") or []:
            raw = ((item.get("totals") or {}).get("data") or {}).get("variant")
            if raw:
                variants.append(_variant(raw, item.get("populations"), label))
        return True, variants
    primary = ((response.get("primary") or {}).get("data") or {}).get("variant")
    if primary:
        variants.append(_variant(primary, response.get("populations"), label))
    for body in (response.get("extras") or {}).values():
        if not isinstance(body, dict):
            continue
        extra = ((body.get("totals") or {}).get("data") or {}).get("variant")
        if extra:
            variants.append(_variant(extra, body.get("populations"), label))
    return True, variants


def frequency_phrase(callsets: list[CallSet]) -> tuple[bool, str]:
    """Return whether the frequency field is unresolved, and the sourced text."""
    usable = [call for call in callsets if _freq(call) is not None]
    if not usable:
        return False, "not reported"
    pieces: list[str] = []
    if len(usable) >= 2:
        first, second = usable[0], usable[1]
        left, right = _freq(first), _freq(second)
        if left is not None and right is not None and disagrees(left, right):
            pieces.extend((_overall(first), _overall(second)))
    pop_set = _population_callset(usable)
    extremes = _extremes(pop_set.populations) if pop_set else None
    if extremes and pop_set is not None:
        low, high = extremes
        overall = _overall(pop_set)
        if overall not in pieces:
            pieces.append(overall)
        pieces.append(f"{pop_set.name} {ANCESTRY[low.id]} {low.ac}/{low.an} alleles")
        pieces.append(f"{pop_set.name} {ANCESTRY[high.id]} {high.ac}/{high.an} alleles")
    if pieces:
        unique: list[str] = []
        for piece in pieces:
            if piece not in unique:
                unique.append(piece)
        return True, "; ".join(unique)
    chosen = max(usable, key=lambda call: (call.an or 0, len(_eligible(call.populations))))
    return False, _overall(chosen)


def disagrees(low: float, high: float) -> bool:
    if high < low:
        low, high = high, low
    if high == low:
        return False
    if low == 0:
        return high >= ZERO_FLOOR
    return (high / low) >= FOLD or (high - low) >= ABS_GAP


def _overall(call: CallSet) -> str:
    if call.af is not None and call.ac is not None and call.an:
        return f"{call.name} allele frequency {format(call.af, '.6g')} ({call.ac}/{call.an} alleles)"
    if call.ac is not None and call.an:
        return f"{call.name} {call.ac}/{call.an} alleles"
    return f"{call.name} allele frequency {format(call.af, '.6g')}"


def _freq(call: CallSet) -> float | None:
    if call.af is not None:
        return float(call.af)
    if call.ac is not None and call.an:
        return call.ac / call.an
    return None


def _eligible(pops: list[PopCount]) -> list[PopCount]:
    return [pop for pop in pops if pop.id in ANCESTRY and pop.an >= MIN_AN]


def _extremes(pops: list[PopCount]) -> tuple[PopCount, PopCount] | None:
    eligible = _eligible(pops)
    if len(eligible) < 2:
        return None
    low = min(eligible, key=lambda pop: (pop.ac / pop.an, -pop.an, pop.id))
    high = max(eligible, key=lambda pop: (pop.ac / pop.an, pop.an, pop.id))
    if low.id == high.id or not disagrees(low.ac / low.an, high.ac / high.an):
        return None
    return low, high


def _population_callset(usable: list[CallSet]) -> CallSet | None:
    """Use the call set whose ancestry extremes disagree by the widest gap."""
    ranked: list[tuple[float, int, CallSet]] = []
    for call in usable:
        pair = _extremes(call.populations)
        if pair is None:
            continue
        low, high = pair
        gap = (high.ac / high.an) - (low.ac / low.an)
        ranked.append((gap, call.an or 0, call))
    if not ranked:
        return None
    return max(ranked, key=lambda item: (item[0], item[1]))[2]


def _variant(raw: dict, populations: dict | None, label: str) -> GnomadVariant:
    callsets: list[CallSet] = []
    for side, nice in (("exome", "exomes"), ("genome", "genomes")):
        built = _callset(raw.get(side), f"{label} {nice}")
        if built is not None:
            callsets.append(built)
    _attach_populations(callsets, populations)
    return GnomadVariant(
        variant_id=str(raw.get("variantId") or ""),
        reference_genome=(raw.get("reference_genome") or None),
        chrom=str(raw.get("chrom") or ""),
        pos=int(raw.get("pos") or 0),
        ref=str(raw.get("ref") or ""),
        alt=str(raw.get("alt") or ""),
        callsets=callsets,
    )


def _callset(block: dict | None, name: str) -> CallSet | None:
    if not block:
        return None
    if block.get("ac") is None and block.get("an") is None and block.get("af") is None:
        return None
    return CallSet(
        name=name,
        ac=block.get("ac"),
        an=block.get("an"),
        af=block.get("af"),
        populations=[],
    )


def _lookup_variants(rsid: str, db: DbSnpRecord, dataset: str) -> list[dict]:
    primary = _gql(_TOTALS_RSID, {"rsid": rsid, "dataset": dataset}, quiet=True)
    variant = ((primary or {}).get("data") or {}).get("variant")
    messages = " ".join(error.get("message", "") for error in primary.get("errors") or [])
    if variant and "Multiple variants" not in messages:
        return [variant, *_other_alleles(variant, db, dataset)]
    if db.assembly != "GRCh38" or not db.chrom or not db.pos:
        return []
    found: list[dict] = []
    for ref, alt in db.alleles[:6]:
        extra = _by_id(f"{_norm(db.chrom)}-{db.pos}-{ref}-{alt}", dataset)
        if extra:
            found.append(extra)
    return found


def _other_alleles(variant: dict, db: DbSnpRecord, dataset: str) -> list[dict]:
    if not _same_site(variant, db):
        return []
    found: list[dict] = []
    for ref, alt in db.alleles[:6]:
        if ref != variant.get("ref") or alt == variant.get("alt"):
            continue
        extra = _by_id(f"{variant['chrom']}-{variant['pos']}-{ref}-{alt}", dataset)
        if extra:
            found.append(extra)
    return found


def _by_id(variant_id: str, dataset: str) -> dict | None:
    body = _gql(_TOTALS_ID, {"variantId": variant_id, "dataset": dataset}, quiet=True)
    variant = ((body or {}).get("data") or {}).get("variant")
    return variant if isinstance(variant, dict) else None


def _attach_populations(callsets: list[CallSet], populations: dict | None) -> None:
    if not populations:
        return
    blocks = populations
    if "body" in populations and "side" in populations:
        blocks = {str(populations["side"]): populations["body"]}
    for side, body in blocks.items():
        if side not in {"exome", "genome"} or not isinstance(body, dict):
            continue
        raw_pops = (((body.get("data") or {}).get("variant") or {}).get(side) or {}).get("populations") or []
        pops = [
            PopCount(str(pop["id"]), int(pop["ac"]), int(pop["an"]))
            for pop in raw_pops
            if pop.get("ac") is not None and pop.get("an") is not None
        ]
        for call in callsets:
            if side == "exome" and call.name.endswith("exomes"):
                call.populations = pops
            if side == "genome" and call.name.endswith("genomes"):
                call.populations = pops


def _populations(variant: dict, dataset: str) -> dict:
    found: dict[str, dict] = {}
    for side in ("exome", "genome"):
        if not variant.get(side):
            continue
        try:
            body = _gql(_POPS[side], {"variantId": variant["variantId"], "dataset": dataset}, quiet=True)
        except RuntimeError as exc:
            print(f"gnomAD populations {variant.get('variantId')} {side}: {exc}", file=sys.stderr)
            continue
        popped = (((body.get("data") or {}).get("variant") or {}).get(side) or {}).get("populations")
        if popped:
            found[side] = body
    return found


def _gql(query: str, variables: dict, *, quiet: bool = False) -> dict:
    del quiet
    time.sleep(0.15)
    body = request_json(GNOMAD_URL, payload={"query": query, "variables": variables})
    variant = ((body.get("data") or {}).get("variant")) if isinstance(body, dict) else None
    if body.get("errors") and not variant:
        message = body["errors"][0].get("message", "gnomAD error")
        if message not in {
            "Variant not found",
            "Multiple variants found, query using variant ID to select one.",
        }:
            print(f"gnomAD: {message}", file=sys.stderr)
    return body


def _same_site(variant: dict, db: DbSnpRecord) -> bool:
    if not db.chrom or not db.pos or not variant.get("chrom") or not variant.get("pos"):
        return False
    if _norm(db.chrom) != _norm(str(variant["chrom"])):
        return False
    if db.pos != int(variant["pos"]):
        return False
    if db.assembly and variant.get("reference_genome") and db.assembly != variant.get("reference_genome"):
        return False
    return True


def _norm(chrom: str) -> str:
    text = chrom.removeprefix("chr")
    if text in {"M", "MT"}:
        return "M"
    return text


def _dataset_label(dataset: str | None) -> str:
    if not dataset:
        return "gnomAD"
    if dataset.startswith("gnomad_r4"):
        return "gnomAD v4"
    if dataset.startswith("gnomad_r3"):
        return "gnomAD v3"
    return dataset


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()
