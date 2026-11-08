"""Newer evidence sources: prediction scores, drug response, cancer evidence, splicing,
regulation, gene constraint and the literature timeline.

Each source has fetch_* (network, returns a cache payload) and parse_* (pure, tolerant of
missing or renamed fields). Response shapes follow each service's public documentation;
`snp-brief --probe RSID` prints what actually comes back so parsers can be checked.

- MyVariant.info (BioThings): dbNSFP scores (REVEL, AlphaMissense, SIFT, PolyPhen-2), CADD,
  SnpEff consequence, and CIViC cancer evidence, in one call.
- ClinPGx (formerly PharmGKB): clinical annotations linking the variant to drug response.
- SpliceAI lookup (Broad Institute): predicted gain or loss of splice sites.
- RegulomeDB: how likely the position is a functional regulatory element.
- gnomAD gene constraint: LOEUF and pLI for the gene.
- Europe PMC: papers per year naming the rsID.
"""

from __future__ import annotations

import urllib.parse
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone

from snp_brief.http_client import request_json

MYVARIANT = "https://myvariant.info/v1/query"
MYVARIANT_FIELDS = ",".join(
    (
        "dbnsfp.revel.score",
        "dbnsfp.alphamissense.score",
        "dbnsfp.alphamissense.pred",
        "dbnsfp.sift.score",
        "dbnsfp.sift.pred",
        "dbnsfp.polyphen2.hdiv.score",
        "dbnsfp.polyphen2.hdiv.pred",
        "dbnsfp.genename",
        "dbnsfp.aa",
        "cadd.phred",
        "cadd.consequence",
        "snpeff.ann",
        "civic",
        "vcf",
        "chrom",
    )
)
CLINPGX_HOSTS = ("https://api.clinpgx.org/v1", "https://api.pharmgkb.org/v1")
SPLICEAI = "https://spliceailookup-api.broadinstitute.org/spliceai/"
REGULOMEDB = "https://regulomedb.org/regulome-search/"
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
EUROPEPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _payload(source: str, request: dict, response) -> dict:
    return {"source": source, "retrieved": _now(), "request": request, "response": response}


def _num(value) -> float | None:
    """A number from a value that may be a string, a list (one per transcript) or missing."""
    if isinstance(value, list):
        numbers = [n for n in (_num(v) for v in value) if n is not None]
        return max(numbers) if numbers else None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _min(value) -> float | None:
    if isinstance(value, list):
        numbers = [n for n in (_num(v) for v in value) if n is not None]
        return min(numbers) if numbers else None
    return _num(value)


def _first(value):
    if isinstance(value, list):
        return next((v for v in value if v not in (None, "", ".")), None)
    return value


def _dig(obj, *path):
    for key in path:
        if isinstance(obj, list):
            obj = obj[0] if obj else None
        if not isinstance(obj, dict):
            return None
        obj = obj.get(key)
    return obj


# ---------------------------------------------------------------- MyVariant.info


@dataclass
class PredictionScores:
    retrieved: bool
    hgvs: str = ""
    gene: str = ""
    consequence: str = ""
    impact: str = ""
    protein_change: str = ""
    cadd_phred: float | None = None
    revel: float | None = None
    alphamissense: float | None = None
    alphamissense_class: str = ""
    sift: float | None = None
    sift_pred: str = ""
    polyphen: float | None = None
    polyphen_pred: str = ""
    civic: list["CivicItem"] = field(default_factory=list)

    @property
    def has_scores(self) -> bool:
        return any(v is not None for v in (self.cadd_phred, self.revel, self.alphamissense, self.sift, self.polyphen))


@dataclass
class CivicItem:
    disease: str
    therapies: list[str]
    evidence_type: str
    significance: str
    level: str
    direction: str
    pmid: str
    statement: str


def fetch_myvariant(rsid: str, timeout: float = 20) -> dict:
    body = request_json(
        MYVARIANT,
        params={"q": f"dbsnp.rsid:{rsid}", "fields": MYVARIANT_FIELDS, "assembly": "hg38", "size": "20"},
        timeout=timeout,
        attempts=2,
    )
    return _payload("MyVariant.info", {"rsid": rsid, "assembly": "hg38"}, body)


_SIFT = {"D": "deleterious", "T": "tolerated"}
_POLYPHEN = {"D": "probably damaging", "P": "possibly damaging", "B": "benign"}
_AM = {"P": "likely pathogenic", "B": "likely benign", "A": "ambiguous"}


def parse_myvariant(payload: dict | None) -> PredictionScores:
    if not payload:
        return PredictionScores(False)
    hits = (payload.get("response") or {}).get("hits") or []
    if not hits:
        return PredictionScores(True)
    # Several alternate alleles come back as separate hits; report the most damaging scores.
    best = PredictionScores(True)
    for hit in hits:
        db = hit.get("dbnsfp") or {}
        cadd = hit.get("cadd") or {}
        revel = _num(_dig(db, "revel", "score"))
        am = _num(_dig(db, "alphamissense", "score"))
        sift = _min(_dig(db, "sift", "score"))
        pp = _num(_dig(db, "polyphen2", "hdiv", "score"))
        phred = _num(cadd.get("phred"))
        if best.revel is None or (revel is not None and revel > best.revel):
            best.revel = revel if revel is not None else best.revel
        if am is not None and (best.alphamissense is None or am > best.alphamissense):
            best.alphamissense = am
            best.alphamissense_class = _AM.get(str(_first(_dig(db, "alphamissense", "pred")) or ""), "")
        if sift is not None and (best.sift is None or sift < best.sift):
            best.sift = sift
            best.sift_pred = _SIFT.get(str(_first(_dig(db, "sift", "pred")) or ""), "")
        if pp is not None and (best.polyphen is None or pp > best.polyphen):
            best.polyphen = pp
            best.polyphen_pred = _POLYPHEN.get(str(_first(_dig(db, "polyphen2", "hdiv", "pred")) or ""), "")
        if phred is not None and (best.cadd_phred is None or phred > best.cadd_phred):
            best.cadd_phred = phred
        ann = (hit.get("snpeff") or {}).get("ann")
        ann = ann[0] if isinstance(ann, list) and ann else ann
        if isinstance(ann, dict) and not best.consequence:
            best.consequence = str(ann.get("effect") or "").split("&")[0].replace("_", " ")
            best.impact = str(ann.get("putative_impact") or "")
            best.gene = str(ann.get("genename") or "")
            best.protein_change = str(ann.get("hgvs_p") or "")
        if not best.gene:
            best.gene = str(_first(db.get("genename")) or "")
        if not best.hgvs:
            best.hgvs = str(hit.get("_id") or "")
        best.civic.extend(_civic_items(hit.get("civic")))
    seen = set()
    unique = []
    for item in best.civic:
        key = (item.disease, tuple(item.therapies), item.significance, item.pmid)
        if key not in seen:
            seen.add(key)
            unique.append(item)
    best.civic = unique
    return best


def _civic_items(civic) -> list[CivicItem]:
    """CIViC evidence as MyVariant nests it (molecularProfiles -> evidenceItems), tolerant of older layouts."""
    if not civic:
        return []
    blocks = civic if isinstance(civic, list) else [civic]
    evidence: list[dict] = []
    for block in blocks:
        if not isinstance(block, dict):
            continue
        for profile in block.get("molecularProfiles") or block.get("molecular_profiles") or []:
            evidence.extend(profile.get("evidenceItems") or profile.get("evidence_items") or [])
        evidence.extend(block.get("evidenceItems") or block.get("evidence_items") or [])
    items = []
    for ev in evidence:
        if not isinstance(ev, dict):
            continue
        disease = _dig(ev, "disease", "name") or _dig(ev, "disease", "displayName") or ""
        therapies = [t.get("name", "") for t in (ev.get("therapies") or ev.get("drugs") or []) if isinstance(t, dict)]
        source = ev.get("source") or {}
        pmid = str(source.get("citationId") or source.get("pubmed_id") or source.get("citation_id") or "")
        items.append(
            CivicItem(
                disease=str(disease),
                therapies=[t for t in therapies if t],
                evidence_type=str(ev.get("evidenceType") or ev.get("evidence_type") or "").title(),
                significance=str(ev.get("significance") or ev.get("clinical_significance") or "").replace("_", " ").title(),
                level=str(ev.get("evidenceLevel") or ev.get("evidence_level") or ""),
                direction=str(ev.get("evidenceDirection") or ev.get("evidence_direction") or "").title(),
                pmid=pmid,
                statement=str(ev.get("description") or ev.get("statement") or "")[:400],
            )
        )
    order = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}
    items.sort(key=lambda i: order.get(i.level, 9))
    return items


# ---------------------------------------------------------------- ClinPGx / PharmGKB


@dataclass
class PgxAnnotation:
    level: str
    drugs: list[str]
    phenotypes: list[str]
    genes: list[str]
    categories: list[str]
    summary: str
    url: str


@dataclass
class PgxResult:
    retrieved: bool
    annotations: list[PgxAnnotation] = field(default_factory=list)


def fetch_pharmgkb(rsid: str, timeout: float = 20) -> dict:
    errors = []
    for host in CLINPGX_HOSTS:
        try:
            body = request_json(
                f"{host}/data/clinicalAnnotation",
                params={"location.fingerprint": rsid, "view": "base"},
                timeout=timeout,
                attempts=2,
            )
            return _payload("ClinPGx", {"rsid": rsid, "host": host}, body)
        except RuntimeError as exc:
            if "HTTP 404" in str(exc):
                return _payload("ClinPGx", {"rsid": rsid, "host": host}, {"data": []})
            errors.append(str(exc)[:160])
    raise RuntimeError("; ".join(errors))


_LEVEL_ORDER = {"1A": 0, "1B": 1, "2A": 2, "2B": 3, "3": 4, "4": 5}


def parse_pharmgkb(payload: dict | None) -> PgxResult:
    if not payload:
        return PgxResult(False)
    body = payload.get("response") or {}
    rows = body.get("data") if isinstance(body, dict) else body
    out = []
    for row in rows or []:
        if not isinstance(row, dict):
            continue
        level = row.get("levelOfEvidence") or row.get("evidenceLevel") or ""
        level = level.get("term", "") if isinstance(level, dict) else str(level)

        def names(*keys):
            for key in keys:
                values = row.get(key)
                if values:
                    return [v.get("name") or v.get("symbol") or "" if isinstance(v, dict) else str(v) for v in values]
            return []

        genes = names("relatedGenes") or [g.get("symbol", "") for g in (_dig(row, "location", "genes") or []) if isinstance(g, dict)]
        drugs = names("relatedChemicals", "chemicals", "drugs")
        if not drugs and row.get("drug"):
            drugs = [str(row["drug"])]
        phenos = names("relatedDiseases", "phenotypes")
        if not phenos and row.get("phenotype"):
            phenos = [str(row["phenotype"])]
        cats = [c.get("term") if isinstance(c, dict) else str(c) for c in (row.get("types") or row.get("phenotypeCategories") or [])]
        summary = ""
        alleles = row.get("allelePhenotypes") or []
        if alleles and isinstance(alleles[0], dict):
            summary = "; ".join(
                f"{a.get('allele', '')}: {a.get('phenotype', '')}".strip(": ") for a in alleles[:3] if a.get("phenotype")
            )
        summary = summary or str(row.get("annotation") or row.get("summary") or "")
        ident = row.get("accessionId") or row.get("id") or ""
        out.append(
            PgxAnnotation(
                level=level,
                drugs=[d for d in drugs if d],
                phenotypes=[p for p in phenos if p],
                genes=[g for g in genes if g],
                categories=[c for c in cats if c],
                summary=summary[:500],
                url=f"https://www.clinpgx.org/clinicalAnnotation/{ident}" if ident else "",
            )
        )
    out.sort(key=lambda a: _LEVEL_ORDER.get(a.level, 9))
    return PgxResult(True, out)


# ---------------------------------------------------------------- SpliceAI


@dataclass
class SpliceResult:
    retrieved: bool
    max_delta: float | None = None
    kind: str = ""  # acceptor gain/loss, donor gain/loss
    gene: str = ""
    allele: str = ""


def fetch_spliceai(variants: list[str], timeout: float = 30) -> dict:
    """variants: chrom-pos-ref-alt strings on GRCh38, one per alternate allele."""
    answers = []
    for variant in variants[:3]:
        body = request_json(SPLICEAI, params={"hg": "38", "variant": variant, "distance": "500", "mask": "1"}, timeout=timeout, attempts=2)
        answers.append({"variant": variant, "body": body})
    return _payload("SpliceAI", {"variants": variants[:3]}, {"answers": answers})


_SPLICE_KINDS = {"DS_AG": "acceptor gain", "DS_AL": "acceptor loss", "DS_DG": "donor gain", "DS_DL": "donor loss"}


def parse_spliceai(payload: dict | None) -> SpliceResult:
    if not payload:
        return SpliceResult(False)
    best = SpliceResult(True)
    for answer in (payload.get("response") or {}).get("answers") or []:
        body = answer.get("body") or {}
        scores = body.get("scores") or []
        if isinstance(scores, dict):
            scores = [scores]
        for row in scores:
            if not isinstance(row, dict):
                continue
            for key, label in _SPLICE_KINDS.items():
                value = _num(row.get(key))
                if value is not None and (best.max_delta is None or value > best.max_delta):
                    best.max_delta = value
                    best.kind = label
                    best.gene = str(row.get("g_name") or row.get("gene_name") or row.get("SYMBOL") or "")
                    best.allele = str(answer.get("variant") or "")
    return best


# ---------------------------------------------------------------- RegulomeDB


@dataclass
class RegulomeResult:
    retrieved: bool
    ranking: str = ""
    probability: float | None = None
    feature_count: int = 0

    @property
    def meaning(self) -> str:
        rank = self.ranking[:1]
        return {
            "1": "likely to affect binding and linked to expression of a gene target",
            "2": "likely to affect binding",
            "3": "less likely to affect binding",
            "4": "minimal binding evidence",
            "5": "minimal binding evidence",
            "6": "minimal binding evidence",
            "7": "no regulatory data at this position",
        }.get(rank, "")


def fetch_regulome(rsid: str, timeout: float = 20) -> dict:
    body = request_json(
        REGULOMEDB, params={"regions": rsid, "genome": "GRCh38", "format": "json", "limit": "100"}, timeout=timeout, attempts=2
    )
    return _payload("RegulomeDB", {"rsid": rsid, "genome": "GRCh38"}, body)


def parse_regulome(payload: dict | None) -> RegulomeResult:
    if not payload:
        return RegulomeResult(False)
    body = payload.get("response") or {}
    score = body.get("regulome_score") or {}
    if isinstance(score, list):
        score = score[0] if score else {}
    features = body.get("features") or []
    count = len(features) if isinstance(features, list) else sum(len(v) for v in features.values() if isinstance(v, list))
    return RegulomeResult(
        True,
        ranking=str(score.get("ranking") or "").strip(),
        probability=_num(score.get("probability")),
        feature_count=count,
    )


# ---------------------------------------------------------------- gnomAD gene constraint

CONSTRAINT_QUERY = """
query Constraint($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    symbol
    gnomad_constraint { pli oe_lof oe_lof_lower oe_lof_upper oe_mis mis_z syn_z }
  }
}
"""


@dataclass
class Constraint:
    retrieved: bool
    gene: str = ""
    loeuf: float | None = None
    oe_lof: float | None = None
    pli: float | None = None
    mis_z: float | None = None

    @property
    def verdict(self) -> str:
        if self.loeuf is None:
            return ""
        if self.loeuf < 0.35:
            return "highly intolerant of loss-of-function changes"
        if self.loeuf < 0.6:
            return "moderately intolerant of loss-of-function changes"
        return "tolerant of loss-of-function changes"


def fetch_constraint(symbol: str, timeout: float = 20) -> dict:
    body = request_json(GNOMAD_API, payload={"query": CONSTRAINT_QUERY, "variables": {"symbol": symbol}}, timeout=timeout, attempts=2)
    return _payload("gnomAD constraint", {"gene": symbol}, body)


def parse_constraint(payload: dict | None) -> Constraint:
    if not payload:
        return Constraint(False)
    gene = _dig(payload.get("response") or {}, "data", "gene") or {}
    c = gene.get("gnomad_constraint") or {}
    return Constraint(
        True,
        gene=str(gene.get("symbol") or (payload.get("request") or {}).get("gene") or ""),
        loeuf=_num(c.get("oe_lof_upper")),
        oe_lof=_num(c.get("oe_lof")),
        pli=_num(c.get("pli")),
        mis_z=_num(c.get("mis_z")),
    )


# ---------------------------------------------------------------- literature timeline


@dataclass
class Timeline:
    retrieved: bool
    by_year: dict[int, int] = field(default_factory=dict)
    total: int = 0

    @property
    def trend(self) -> str:
        if not self.by_year:
            return ""
        this_year = datetime.now(timezone.utc).year
        recent = sum(n for y, n in self.by_year.items() if this_year - 3 <= y < this_year)
        earlier = sum(n for y, n in self.by_year.items() if this_year - 6 <= y < this_year - 3)
        if recent >= 6 and recent > 1.5 * max(earlier, 1):
            return "growing"
        if earlier >= 6 and recent < 0.6 * earlier:
            return "declining"
        return "steady"


def fetch_timeline(rsid: str, timeout: float = 20) -> dict:
    body = request_json(
        EUROPEPMC_SEARCH,
        params={"query": f'"{rsid}"', "format": "json", "resultType": "lite", "pageSize": "1000"},
        timeout=timeout,
        attempts=2,
    )
    rows = ((body.get("resultList") or {}).get("result") or []) if isinstance(body, dict) else []
    years = Counter(str(r.get("pubYear")) for r in rows if str(r.get("pubYear") or "").isdigit())
    slim = {"hitCount": body.get("hitCount") if isinstance(body, dict) else None, "years": dict(years)}
    return _payload("Europe PMC", {"rsid": rsid, "query": f'"{rsid}"'}, slim)


def parse_timeline(payload: dict | None) -> Timeline:
    if not payload:
        return Timeline(False)
    body = payload.get("response") or {}
    by_year = {int(y): int(n) for y, n in (body.get("years") or {}).items() if str(y).isdigit()}
    total = int(body.get("hitCount") or sum(by_year.values()))
    return Timeline(True, dict(sorted(by_year.items())), total)


def splice_variant_ids(chrom: str, pos: int, ref: str, alt: str) -> str:
    return f"chr{str(chrom).removeprefix('chr')}-{pos}-{ref}-{alt}"


def quote(text: str) -> str:
    return urllib.parse.quote(text, safe="")
