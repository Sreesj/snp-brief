"""Evidence providers, ranked by dependability, plus the two newest: Ensembl VEP and GWAS Catalog.

Every provider is queried in parallel (pipeline.gather). For each field the
highest-ranked provider that answered wins; lower-ranked answers are kept as
supporting or fallback evidence. A time budget trades dependability for speed:
providers that have not answered by the deadline are reported as still loading,
finish in the background, and are cached for the next lookup.

Rank (1 = most dependable):

1. ClinVar        curated clinical classification
2. gnomAD         population measurement
3. GWAS Catalog   curated statistical trait association (association, not causation)
3. GTEx           expression QTLs: genotype linked to a gene's expression in a tissue
4. Ensembl LD     population linkage: SNPs inherited together with this one
4. dbSNP          reference annotation (position, alleles, location class)
5. Ensembl VEP    computed annotation: consequences, regulatory overlap, prediction scores
5. Open Targets   computed locus-to-gene score
6. UniProt        gene-level protein function (describes the protein, not the variant)
6. NCBI Gene      gene summary, also for non-coding genes UniProt does not cover
7. Europe PMC     papers that name the rsID (title/abstract only)
7. dbSNP citations  PubMed records NCBI links to the rsID
7. LitVar2        publications mentioning the rsID anywhere in the full text
8. AI background  the language model's general knowledge (Hypothesis Forge only, labelled, unverified)
"""

from __future__ import annotations

import os
import sys
import urllib.parse
from dataclasses import dataclass, field
from datetime import datetime, timezone

from snp_brief.http_client import request_json


@dataclass(frozen=True)
class Provider:
    key: str
    name: str
    rank: int
    tier: str
    fields: tuple[str, ...]
    timeout: float


PROVIDERS = (
    Provider("clinvar", "ClinVar", 1, "curated clinical", ("function", "disease"), 30),
    Provider("gnomad", "gnomAD", 2, "population data", ("frequency",), 30),
    Provider("gwas", "GWAS Catalog", 3, "statistical association", ("disease",), 60),
    Provider("gtex", "GTEx", 3, "expression association (eQTL)", ("function",), 20),
    Provider("ensembl_ld", "Ensembl LD", 4, "population linkage", ("context",), 20),
    Provider("dbsnp", "dbSNP", 4, "reference annotation", ("identity", "function"), 20),
    Provider("ensembl_vep", "Ensembl VEP", 5, "computed annotation and predictions", ("function", "prediction"), 30),
    Provider("open_targets", "Open Targets", 5, "computed locus-to-gene score", ("function",), 20),
    Provider("uniprot", "UniProt", 6, "gene-level function", ("function",), 20),
    Provider("ncbi_gene", "NCBI Gene", 6, "gene-level function", ("function",), 20),
    Provider("europepmc", "Europe PMC", 7, "literature mention", ("papers",), 20),
    Provider("dbsnp_citations", "dbSNP citations", 7, "literature mention", ("papers",), 20),
    Provider("litvar", "LitVar2", 7, "literature mention", ("papers",), 20),
    Provider("myvariant", "MyVariant.info", 5, "computed prediction scores (dbNSFP, CADD) and CIViC", ("prediction",), 20),
    Provider("pharmgkb", "ClinPGx (PharmGKB)", 2, "curated drug-response annotation", ("drug",), 20),
    Provider("spliceai", "SpliceAI", 5, "computed splice prediction", ("prediction",), 30),
    Provider("regulome", "RegulomeDB", 5, "regulatory evidence score", ("prediction",), 20),
    Provider("constraint", "gnomAD constraint", 2, "gene-level population constraint", ("context",), 20),
    Provider("timeline", "Literature timeline", 7, "papers per year (Europe PMC)", ("papers",), 20),
    Provider("borzoi", "Borzoi (our model run)", 5, "computed regulatory prediction, run by Variant Atlas", ("prediction",), 5),
)
PROVIDERS = tuple(sorted(PROVIDERS, key=lambda p: p.rank))  # stable: keeps the listed order within a rank
BY_KEY = {provider.key: provider for provider in PROVIDERS}
AI_TIER = "AI background knowledge (unverified)"

SPEED_BUDGETS = {"fast": 8.0, "balanced": 20.0, "thorough": 45.0}


OPTIONAL_OFF = {"regulome", "spliceai"}  # public services that failed on every live check; opt in with SNP_BRIEF_ENABLE


def source_enabled(key: str) -> bool:
    """RegulomeDB and SpliceAI's public API were unreachable or timing out on every live check, so they are
    off unless SNP_BRIEF_ENABLE lists them (e.g. SNP_BRIEF_ENABLE=regulome,spliceai)."""
    if key not in OPTIONAL_OFF:
        return True
    return key in {k.strip() for k in os.environ.get("SNP_BRIEF_ENABLE", "").split(",")}


def extras_enabled() -> bool:
    """SNP_BRIEF_EXTRAS=off turns off the newer providers (GWAS Catalog, Ensembl VEP)."""
    return os.environ.get("SNP_BRIEF_EXTRAS", "on").lower() not in {"off", "0", "false", "no"}


def tier_label(source: str) -> str:
    if source == "Paper abstract":
        return "rank 7, literature finding (LLM-extracted, quote-verified)"
    for provider in PROVIDERS:
        if provider.name == source:
            return f"rank {provider.rank}, {provider.tier}"
    return AI_TIER if source == "AI" else ""


# ---------------------------------------------------------------- Ensembl VEP

VEP_URL = "https://rest.ensembl.org/vep/human/id/{rsid}"
# Parameter names as documented at rest.ensembl.org/documentation/info/vep_id_get.
VEP_PARAMS = {"AlphaMissense": "1", "CADD": "1", "SpliceAI": "1", "canonical": "1", "variant_class": "1", "domains": "1"}
VEP_BASIC = {"canonical": "1", "variant_class": "1", "domains": "1"}


@dataclass
class VepAnnotation:
    retrieved: bool
    most_severe: str = ""
    variant_class: str = ""
    genes: list[str] = field(default_factory=list)
    canonical_gene: str = ""
    canonical_terms: list[str] = field(default_factory=list)
    impact: str = ""
    sift: str = ""
    polyphen: str = ""
    cadd_phred: float | None = None
    alphamissense: str = ""
    alphamissense_score: float | None = None
    spliceai_max: float | None = None
    regulatory: list[str] = field(default_factory=list)
    motifs: list[str] = field(default_factory=list)
    pubmed: list[str] = field(default_factory=list)
    amino_acids: str = ""
    protein_position: str = ""
    domains: list[str] = field(default_factory=list)


def fetch_vep(rsid: str, timeout: float = 25) -> dict:
    """With prediction plugins first; they make Ensembl slow, so on failure retry once without them."""
    url = VEP_URL.format(rsid=urllib.parse.quote(rsid))
    try:
        body = request_json(url, params={**VEP_PARAMS, "content-type": "application/json"}, timeout=timeout, attempts=2)
        used = VEP_PARAMS
    except RuntimeError as first:
        if "HTTP 400" in str(first):  # unknown rsID: no point retrying
            raise
        used = VEP_BASIC
        body = request_json(url, params={**VEP_BASIC, "content-type": "application/json"}, timeout=timeout, attempts=2)
    return {"source": "Ensembl VEP", "retrieved": _now(), "request": {"rsid": rsid, "params": used}, "response": body}


def parse_vep(payload: dict | None) -> VepAnnotation:
    if not payload:
        return VepAnnotation(False)
    body = payload.get("response")
    records = body if isinstance(body, list) else [body] if isinstance(body, dict) else []
    if not records or not isinstance(records[0], dict) or records[0].get("error"):
        return VepAnnotation(True)
    record = records[0]
    out = VepAnnotation(True, str(record.get("most_severe_consequence") or ""), str(record.get("variant_class") or ""))
    transcripts = [t for t in record.get("transcript_consequences") or [] if isinstance(t, dict)]
    for item in transcripts:
        gene = str(item.get("gene_symbol") or "")
        if gene and gene not in out.genes:
            out.genes.append(gene)
    canonical = next((t for t in transcripts if t.get("canonical")), None)
    severe = [t for t in transcripts if out.most_severe in (t.get("consequence_terms") or [])]
    chosen = canonical if canonical is not None and (not severe or canonical in severe) else (severe[0] if severe else canonical)
    if chosen is not None:
        out.canonical_gene = str(chosen.get("gene_symbol") or "")
        out.canonical_terms = [str(t) for t in chosen.get("consequence_terms") or []]
        out.impact = str(chosen.get("impact") or "")
        out.amino_acids = str(chosen.get("amino_acids") or "")
        start = chosen.get("protein_start")
        out.protein_position = str(start) if start is not None else ""
        for domain in chosen.get("domains") or []:
            if isinstance(domain, dict):
                label = f"{domain.get('db', '')} {domain.get('name', '')}".strip()
            else:
                label = str(domain)
            if label and label not in out.domains:
                out.domains.append(label)
    for item in transcripts:
        out.sift = out.sift or str(item.get("sift_prediction") or "")
        out.polyphen = out.polyphen or str(item.get("polyphen_prediction") or "")
        if out.cadd_phred is None and item.get("cadd_phred") is not None:
            out.cadd_phred = _num(item.get("cadd_phred"))
        am = item.get("alphamissense")
        if isinstance(am, dict) and not out.alphamissense:
            out.alphamissense = str(am.get("am_class") or "")
            out.alphamissense_score = _num(am.get("am_pathogenicity"))
        splice = item.get("spliceai")
        if isinstance(splice, dict):
            scores = [_num(v) for k, v in splice.items() if str(k).upper().startswith("DS_")]
            scores = [s for s in scores if s is not None]
            if scores:
                out.spliceai_max = max(scores + ([out.spliceai_max] if out.spliceai_max is not None else []))
    if out.cadd_phred is None and record.get("cadd_phred") is not None:
        out.cadd_phred = _num(record.get("cadd_phred"))
    for item in record.get("regulatory_feature_consequences") or []:
        label = str((item or {}).get("biotype") or "").replace("_", " ")
        if label and label not in out.regulatory:
            out.regulatory.append(label)
    for item in record.get("motif_feature_consequences") or []:
        for factor in (item or {}).get("transcription_factors") or []:
            if factor not in out.motifs:
                out.motifs.append(str(factor))
    for variant in record.get("colocated_variants") or []:
        for pmid in (variant or {}).get("pubmed") or []:
            if str(pmid) not in out.pubmed:
                out.pubmed.append(str(pmid))
    return out


# ---------------------------------------------------------------- GWAS Catalog

GWAS_V2 = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
# The v1 API was retired (HTTP 410). v2 per the migration guide: /v2/associations?rs_id=..., page/size, default 20.
# The v2 rs_id query is slow (often over 30 s), so ask for small pages and allow up to 60 s.
GWAS_PAGE_SIZE = 10
GWAS_MAX_PAGES = 2


@dataclass
class GwasHit:
    trait: str
    p_value: float | None
    p_text: str
    risk_allele: str = ""
    effect: str = ""
    genes: list[str] = field(default_factory=list)
    pmid: str = ""
    reported_allele: str = ""  # set when the study reported the risk allele on the opposite strand


_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def orient_risk_alleles(hits: list[GwasHit], alleles: list[tuple[str, str]]) -> None:
    """Studies sometimes report the risk allele on the reverse strand (A for a C>T variant).
    When only the complement is an allele of this variant, report it on the forward strand.
    Palindromic variants (A/T, C/G) are ambiguous and left as reported."""
    letters = {a for pair in alleles for a in pair if len(a) == 1}
    if not letters or any(_COMPLEMENT.get(a) in letters for a in letters):
        return
    for hit in hits:
        risk = hit.risk_allele
        if len(risk) == 1 and risk not in letters and _COMPLEMENT.get(risk) in letters:
            hit.reported_allele, hit.risk_allele = risk, _COMPLEMENT[risk]


@dataclass
class GwasResult:
    retrieved: bool
    hits: list[GwasHit] = field(default_factory=list)


ENSEMBL_VARIATION = "https://rest.ensembl.org/variation/human/{rsid}"
ENSEMBL_PHENOTYPE_VARIANT = "https://rest.ensembl.org/phenotype/variant/human/{rsid}"


def fetch_gwas(rsid: str, timeout: float = 60) -> dict:
    """GWAS Catalog associations: EBI's v2 API first (10 s), then the same catalogue via Ensembl.

    EBI's v2 rs_id query often fails server-side (HTTP 500 after ~30 s). Ensembl imports the
    NHGRI-EBI GWAS Catalog into its variation phenotypes, so it is a dependable fallback.
    """
    from snp_brief.sources import gwas_local

    root = gwas_local.default_root()
    if root is not None:
        try:
            local = gwas_local.lookup(rsid, root)
            if local is not None:
                return local  # the full catalogue on disk: instant, and an empty answer is a real "none"
        except Exception as exc:  # a broken local file falls back to the APIs
            print(f"{rsid} gwas: local catalogue unreadable ({exc}); using the API", file=sys.stderr)
    errors: list[str] = []
    pages: list = []
    for page in range(GWAS_MAX_PAGES):
        try:
            body = request_json(
                GWAS_V2, params={"rs_id": rsid, "page": str(page), "size": str(GWAS_PAGE_SIZE)},
                timeout=min(timeout, 10), attempts=1,
            )
        except RuntimeError as exc:
            if "HTTP 404" in str(exc):
                return {"source": "GWAS Catalog", "retrieved": _now(), "request": {"rsid": rsid, "api": "v2"}, "response": {"pages": []}}
            errors.append(f"EBI v2: {_short_error(exc)}")
            break
        pages.append(body)
        if len(_rows(body)) < GWAS_PAGE_SIZE:
            break
    if pages:
        return {"source": "GWAS Catalog", "retrieved": _now(), "request": {"rsid": rsid, "api": "v2"}, "response": {"pages": pages}}
    # Ensembl, two ways: the variation endpoint with phenotypes, then the lighter phenotype
    # endpoint (the first sometimes returns HTTP 500 for variants with many associations).
    for api, url, params in (
        ("ensembl-phenotypes", ENSEMBL_VARIATION.format(rsid=urllib.parse.quote(rsid)), {"phenotypes": "1", "content-type": "application/json"}),
        ("ensembl-phenotype-variant", ENSEMBL_PHENOTYPE_VARIANT.format(rsid=urllib.parse.quote(rsid)), {"content-type": "application/json"}),
    ):
        try:
            body = request_json(url, params=params, timeout=min(timeout, 30), attempts=2)
        except RuntimeError as exc:
            errors.append(f"Ensembl ({api}): {_short_error(exc)}")
            continue
        rows = body if isinstance(body, list) else (body.get("phenotypes") or [])
        return {
            "source": "GWAS Catalog",
            "retrieved": _now(),
            "request": {"rsid": rsid, "api": api},
            "response": {"ensembl_phenotypes": [row for row in rows if _is_gwas(row)]},
        }
    raise RuntimeError("; ".join(errors))


def _short_error(exc: Exception) -> str:
    """One line: drop HTML error pages that some servers return with HTTP 500."""
    text = str(exc)
    cut = text.find("<")
    return (text[:cut].rstrip(": ") if cut > 0 else text)[:160]


def _p_order(hit) -> float:
    """log10 of the p-value, read from the text when the float underflows (e.g. 2e-9237 -> -9236.7)."""
    import math
    import re

    match = re.fullmatch(r"\s*(\d+(?:\.\d+)?)[eE]([-+]?\d+)\s*", hit.p_text or "")
    if match and float(match.group(1)) > 0:
        return math.log10(float(match.group(1))) + int(match.group(2))
    if hit.p_value is not None and hit.p_value > 0:
        return math.log10(hit.p_value)
    return 0.0


def _is_gwas(row) -> bool:
    return isinstance(row, dict) and "gwas" in str(row.get("source") or "").lower()


def parse_gwas(payload: dict | None, rsid: str = "") -> GwasResult:
    if not payload:
        return GwasResult(False)
    rsid = rsid or str((payload.get("request") or {}).get("rsid") or "")
    body = payload.get("response") or {}
    hits: list[GwasHit] = []
    if isinstance(body, dict) and "ensembl_phenotypes" in body:
        for row in body.get("ensembl_phenotypes") or []:
            if _haplotype_row(row, rsid):
                continue  # an association for a multi-SNP haplotype is not this SNP alone
            hit = _ensembl_gwas_hit(row)
            if hit is not None:
                hits.append(hit)
        hits.sort(key=_p_order)
        return GwasResult(True, _dedupe(hits))
    bodies = body.get("pages") if isinstance(body, dict) and "pages" in body else [body]
    for page in bodies:
        for row in _rows(page):
            hit = _gwas_hit(row)
            if hit is not None:
                hits.append(hit)
    hits.sort(key=_p_order)
    return GwasResult(True, hits)


_ASSOCIATION_KEYS = {
    "p_value", "pValue", "pvalue", "pvalue_mantissa", "pvalueMantissa", "efo_traits", "efoTraits",
    "reported_trait", "reportedTrait", "association_id", "associationId",
}


def _rows(body) -> list[dict]:
    """The association list, wherever the response wraps it (_embedded, content, data, results, or a bare list)."""
    if isinstance(body, list):
        return [r for r in body if isinstance(r, dict)]
    if not isinstance(body, dict):
        return []
    for key in ("_embedded", "content", "data", "results", "associations"):
        inner = body.get(key)
        if isinstance(inner, dict):
            found = _rows(inner)
            if found:
                return found
        if isinstance(inner, list) and inner and isinstance(inner[0], dict) and _ASSOCIATION_KEYS & set(inner[0]):
            return inner
    for value in body.values():
        if isinstance(value, list) and value and isinstance(value[0], dict) and _ASSOCIATION_KEYS & set(value[0]):
            return value
    return []


def _haplotype_row(row, rsid: str) -> bool:
    variants = [v.strip() for v in str((row or {}).get("variants") or "").split(",") if v.strip()]
    return len(variants) > 1 or (bool(variants) and bool(rsid) and variants[0] != rsid)


def _effect_text(odds_raw, beta_raw) -> str:
    """Effect size as written: numeric values tidied, text such as '0.3284 unit decrease' kept."""
    for label, raw in (("odds ratio", odds_raw), ("beta", beta_raw)):
        if raw in (None, ""):
            continue
        number = _num(raw)
        return f"{label} {number:g}" if number is not None else f"{label} {str(raw).strip()}"
    return ""


def _ensembl_gwas_hit(row: dict) -> GwasHit | None:
    """One Ensembl phenotype row sourced from the GWAS Catalog."""
    if not isinstance(row, dict):
        return None
    trait = str(row.get("trait") or row.get("description") or "").strip()
    if not trait:
        return None
    attrs = row.get("attributes") if isinstance(row.get("attributes"), dict) else {}

    def get(*keys):
        for key in keys:
            for source in (row, attrs):
                value = source.get(key)
                if value not in (None, ""):
                    return value
        return None

    p_raw = get("pvalue", "p_value", "p-value")
    p_value = _num(p_raw)
    effect = _effect_text(get("odds_ratio"), get("beta_coefficient", "beta")) or str(get("effect_text") or "")
    genes_raw = get("genes", "associated_gene", "reported_genes") or ""
    genes = [g.strip() for g in str(genes_raw).replace(";", ",").split(",") if g.strip() and g.strip() != "NR"]
    study = str(get("study", "external_reference", "pubmed_id") or "")
    pmid = study.split("PMID:")[-1].strip() if "PMID:" in study else (study if study.isdigit() else "")
    risk = str(get("risk_allele") or "")
    if "-" in risk and risk.startswith("rs"):
        risk = risk.split("-", 1)[1]
    return GwasHit(trait, p_value, str(p_raw) if p_raw is not None else "", risk, effect, genes[:4], pmid)


def _dedupe(hits: list[GwasHit]) -> list[GwasHit]:
    seen: set[tuple[str, str]] = set()
    out = []
    for hit in hits:
        key = (hit.trait.lower(), hit.pmid)
        if key not in seen:
            seen.add(key)
            out.append(hit)
    return out


def _gwas_hit(row: dict) -> GwasHit | None:
    traits: list[str] = []
    for item in row.get("efo_traits") or row.get("efoTraits") or []:
        name = (item.get("efo_trait") or item.get("efoTrait") or item.get("trait") or item.get("label")) if isinstance(item, dict) else item
        if name and str(name) not in traits:
            traits.append(str(name))
    reported = row.get("reported_trait") or row.get("reportedTrait") or ((row.get("study") or {}).get("diseaseTrait") or {}).get("trait")
    if not traits and reported:
        traits = [str(reported[0] if isinstance(reported, list) else reported)]
    if not traits:
        return None
    mantissa = row.get("pvalue_mantissa", row.get("pvalueMantissa"))
    exponent = row.get("pvalue_exponent", row.get("pvalueExponent"))
    p_value = _num(row.get("p_value", row.get("pValue", row.get("pvalue"))))
    if p_value is None and mantissa is not None and exponent is not None:
        try:
            p_value = float(mantissa) * 10 ** int(exponent)
        except (TypeError, ValueError):
            p_value = None
    p_text = f"{mantissa}e{exponent}" if mantissa is not None and exponent is not None else (f"{p_value:.1e}" if p_value else "")
    risk = row.get("snp_effect_allele") or row.get("snpEffectAllele") or row.get("risk_allele") or row.get("riskAllele")
    if isinstance(risk, list):
        risk = risk[0] if risk else ""
    if not risk:
        for locus in row.get("loci") or []:
            for allele in (locus or {}).get("strongestRiskAlleles") or []:
                risk = risk or (allele or {}).get("riskAlleleName")
    effect = ""
    odds = _num(row.get("or_value", row.get("orValue", row.get("orPerCopyNum"))))
    beta = _num(row.get("beta", row.get("betaNum")))
    if odds:
        effect = f"odds ratio {odds:g}"
    elif beta:
        effect = f"beta {beta:g}"
    genes: list[str] = []
    for gene in row.get("mapped_genes") or row.get("mappedGenes") or row.get("reported_genes") or row.get("reportedGenes") or []:
        name = gene.get("geneName") if isinstance(gene, dict) else gene
        if name and str(name) not in genes:
            genes.append(str(name))
    for locus in row.get("loci") or []:
        for gene in (locus or {}).get("authorReportedGenes") or []:
            name = (gene or {}).get("geneName")
            if name and name not in genes:
                genes.append(name)
    pmid = row.get("pubmed_id") or row.get("pubmedId") or (((row.get("study") or {}).get("publicationInfo") or {}).get("pubmedId")) or ""
    return GwasHit(", ".join(traits), p_value, p_text, str(risk or ""), effect, genes[:4], str(pmid))


def _num(value) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


# ---------------------------------------------------------------- NCBI: dbSNP citations and Gene summary

REFSNP_URL = "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{number}"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_citations(rsid: str, timeout: float = 20) -> dict:
    body = request_json(REFSNP_URL.format(number=rsid.removeprefix("rs")), timeout=timeout, attempts=2)
    pmids = [str(p) for p in (body.get("citations") or [])]
    # Keep only what is used; the full RefSNP record is large.
    return {"source": "dbSNP citations", "retrieved": _now(), "request": {"rsid": rsid}, "response": {"citations": pmids}}


def parse_citations(payload: dict | None) -> list[str]:
    if not payload:
        return []
    return [str(p) for p in ((payload.get("response") or {}).get("citations") or [])]


@dataclass
class GeneSummary:
    symbol: str
    name: str
    summary: str
    gene_id: str


def fetch_gene_summary(symbol: str, timeout: float = 20) -> dict:
    search = request_json(
        f"{EUTILS}/esearch.fcgi",
        params={"db": "gene", "term": f"{symbol}[sym] AND Homo sapiens[orgn]", "retmode": "json", "tool": "snp_brief"},
        timeout=timeout,
    )
    ids = ((search.get("esearchresult") or {}).get("idlist")) or []
    body: dict = {}
    if ids:
        body = request_json(
            f"{EUTILS}/esummary.fcgi",
            params={"db": "gene", "id": ids[0], "retmode": "json", "tool": "snp_brief"},
            timeout=timeout,
            attempts=2,
        )
    return {"source": "NCBI Gene", "retrieved": _now(), "request": {"gene": symbol}, "response": body}


def parse_gene_summary(payload: dict | None, symbol: str) -> GeneSummary | None:
    if not payload:
        return None
    if ((payload.get("request") or {}).get("gene")) not in {None, symbol}:
        return None
    result = (payload.get("response") or {}).get("result") or {}
    for uid in result.get("uids") or []:
        doc = result.get(str(uid)) or {}
        if str(doc.get("name") or "") != symbol:
            continue
        return GeneSummary(symbol, str(doc.get("description") or ""), str(doc.get("summary") or ""), str(uid))
    return None


# ---------------------------------------------------------------- GTEx eQTLs

GTEX = "https://gtexportal.org/api/v2"


@dataclass
class Eqtl:
    gene: str
    tissue: str
    p_value: float | None
    nes: float | None
    variant_id: str = ""  # GTEx id, chr_pos_ref_alt_b38: the effect is for this alt allele


def gtex_variant_id(chrom: str | None, pos: int | None, ref: str, alt: str) -> str:
    return f"chr{chrom}_{pos}_{ref}_{alt}_b38" if chrom and pos and ref and alt else ""


def fetch_gtex(rsid: str, fallback_ids: list[str], timeout: float = 20) -> dict:
    """Resolve the rsID to GTEx variant ids, then fetch significant single-tissue eQTLs."""
    ids: list[str] = []
    try:
        lookup = request_json(f"{GTEX}/dataset/variant", params={"snpId": rsid, "datasetId": "gtex_v8"}, timeout=timeout, attempts=2)
        for row in lookup.get("data") or []:
            vid = (row or {}).get("variantId")
            if vid and vid not in ids:
                ids.append(str(vid))
    except RuntimeError:
        pass
    ids = ids or [vid for vid in fallback_ids if vid]
    rows: list[dict] = []
    for vid in ids[:3]:
        body = request_json(
            f"{GTEX}/association/singleTissueEqtl",
            params={"variantId": vid, "datasetId": "gtex_v8", "itemsPerPage": "250"},
            timeout=timeout,
            attempts=2,
        )
        rows.extend(r for r in body.get("data") or [] if isinstance(r, dict))
    return {"source": "GTEx", "retrieved": _now(), "request": {"rsid": rsid, "variant_ids": ids[:3]}, "response": {"data": rows}}


def parse_gtex(payload: dict | None) -> list[Eqtl]:
    if not payload:
        return []
    out: list[Eqtl] = []
    for row in (payload.get("response") or {}).get("data") or []:
        gene = str(row.get("geneSymbol") or row.get("gencodeId") or "")
        tissue = str(row.get("tissueSiteDetailId") or row.get("tissueSiteDetail") or "").replace("_", " ")
        if gene and tissue:
            out.append(Eqtl(gene, tissue, _num(row.get("pValue")), _num(row.get("nes")), str(row.get("variantId") or "")))
    out.sort(key=lambda e: e.p_value if e.p_value is not None else 1.0)
    return out


# ---------------------------------------------------------------- Ensembl LD

LD_URL = "https://rest.ensembl.org/ld/human/{rsid}/1000GENOMES:phase_3:{pop}"


@dataclass
class LdProxy:
    rsid: str
    r2: float


def fetch_ld(rsid: str, timeout: float = 20, population: str = "EUR") -> dict:
    body = request_json(
        LD_URL.format(rsid=urllib.parse.quote(rsid), pop=population),
        params={"r2": "0.8", "window_size": "250", "content-type": "application/json"},
        timeout=timeout,
        attempts=2,
    )
    return {"source": "Ensembl LD", "retrieved": _now(), "request": {"rsid": rsid, "population": population}, "response": body}


def parse_ld(payload: dict | None, rsid: str) -> list[LdProxy]:
    if not payload:
        return []
    body = payload.get("response")
    rows = body if isinstance(body, list) else []
    out: dict[str, float] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        other = row.get("variation2") if row.get("variation1") == rsid else row.get("variation1")
        r2 = _num(row.get("r2"))
        if other and other != rsid and r2 is not None and r2 >= 0.8:
            out[str(other)] = max(r2, out.get(str(other), 0.0))
    return sorted((LdProxy(k, v) for k, v in out.items()), key=lambda p: -p.r2)


# ---------------------------------------------------------------- Open Targets locus-to-gene

OPEN_TARGETS = "https://api.platform.opentargets.org/api/v4/graphql"
OT_QUERY = """query VariantL2G($id: String!) {
  variant(variantId: $id) {
    id
    credibleSets(page: {size: 10, index: 0}) {
      rows {
        study { traitFromSource }
        l2GPredictions(page: {size: 3, index: 0}) { rows { score target { approvedSymbol } } }
      }
    }
  }
}"""


@dataclass
class L2G:
    gene: str
    score: float
    trait: str


def ot_variant_id(chrom: str | None, pos: int | None, ref: str, alt: str) -> str:
    return f"{chrom}_{pos}_{ref}_{alt}" if chrom and pos and ref and alt else ""


def fetch_open_targets(variant_ids: list[str], timeout: float = 20) -> dict:
    rows = []
    for vid in [v for v in variant_ids if v][:3]:
        body = request_json(OPEN_TARGETS, payload={"query": OT_QUERY, "variables": {"id": vid}}, timeout=timeout, attempts=2)
        rows.append({"variant_id": vid, "body": body})
    return {"source": "Open Targets", "retrieved": _now(), "request": {"variant_ids": variant_ids[:3]}, "response": rows}


def parse_open_targets(payload: dict | None) -> list[L2G]:
    if not payload:
        return []
    best: dict[tuple[str, str], float] = {}
    for item in payload.get("response") or []:
        variant = (((item or {}).get("body") or {}).get("data") or {}).get("variant") or {}
        for cs in ((variant.get("credibleSets") or {}).get("rows")) or []:
            trait = str(((cs or {}).get("study") or {}).get("traitFromSource") or "")
            for pred in (((cs or {}).get("l2GPredictions") or {}).get("rows")) or []:
                gene = str(((pred or {}).get("target") or {}).get("approvedSymbol") or "")
                score = _num((pred or {}).get("score"))
                if gene and score is not None:
                    key = (gene, trait)
                    best[key] = max(score, best.get(key, 0.0))
    out = [L2G(gene, score, trait) for (gene, trait), score in best.items()]
    out.sort(key=lambda r: -r.score)
    return out


# ---------------------------------------------------------------- LitVar2

LITVAR = "https://www.ncbi.nlm.nih.gov/research/litvar2-api/variant/get/{vid}/publications"


def fetch_litvar(rsid: str, timeout: float = 20) -> dict:
    vid = urllib.parse.quote(f"litvar@{rsid}##", safe="")
    body = request_json(LITVAR.format(vid=vid), timeout=timeout, attempts=2)
    return {"source": "LitVar2", "retrieved": _now(), "request": {"rsid": rsid}, "response": body}


def parse_litvar(payload: dict | None) -> list[str]:
    if not payload:
        return []
    body = payload.get("response") or {}
    pmids = body.get("pmids") if isinstance(body, dict) else body if isinstance(body, list) else None
    return [str(p) for p in (pmids or []) if str(p).strip()]
