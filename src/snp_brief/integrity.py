"""Checks that run on every loaded variant, cached or live, before anyone reads the page.

Each check re-derives what the page shows from the raw saved responses, or tests that a value
is possible at all, independently of the parsers that produced it:

- gnomAD counts and frequencies on the page equal the raw response (af = ac / an)
- every percentage in the summary and frequency card is computable from the raw counts
- carrier and two-copy rates equal 1 - (1 - p)^2 and p^2 for a raw p
- dbSNP and gnomAD agree on position; gnomAD alleles are among dbSNP's alleles
- prediction scores are inside their published ranges
- GWAS p-values are between 0 and 1; risk alleles are alleles of this variant
- ClinPGx levels and CIViC levels are valid; PMIDs are numbers; paper years are not in the future
- frequency-based ACMG codes never contradict a pathogenic ClinVar call without a warning
- sources that name a gene agree, or the disagreement is reported

A failed check does not hide data; it is listed so a reviewer knows exactly what to distrust.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class CheckReport:
    passed: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)  # worth a look, not necessarily wrong
    failures: list[str] = field(default_factory=list)  # a value on the page does not hold up

    @property
    def ok(self) -> bool:
        return not self.failures

    def add(self, name: str, problems: list[str], *, warn: bool = False) -> None:
        if problems:
            (self.warnings if warn else self.failures).extend(problems)
        else:
            self.passed.append(name)


def _raw_gnomad(item) -> dict | None:
    for directory in [getattr(item, "cache_dir", None)]:
        if directory is None:
            continue
        path = Path(directory) / item.rsid / "gnomad.json"
        if path.is_file():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                return None
    return None


def _raw_counts(raw: dict | None) -> dict[tuple[str, str, str], tuple[int, int]]:
    counts: dict[tuple[str, str, str], tuple[int, int]] = {}
    if not raw:
        return counts
    for v in (raw.get("response") or {}).get("variants") or []:
        var = (((v.get("totals") or {}).get("data") or {}).get("variant")) or {}
        for cs in ("exome", "genome"):
            d = var.get(cs)
            if d and d.get("an"):
                counts[(var.get("variantId"), cs, "all")] = (d.get("ac") or 0, d["an"])
        for cs, block in (v.get("populations") or {}).items():
            pv = ((block or {}).get("data") or {}).get("variant") or {}
            for pop in ((pv.get(cs) or {}).get("populations") or []):
                if pop.get("an"):
                    counts[(pv.get("variantId"), cs, pop["id"])] = (pop.get("ac") or 0, pop["an"])
    return counts


def _percent_forms(counts) -> set[float]:
    forms: set[float] = set()
    for ac, an in counts.values():
        p = ac / an * 100
        for sig in (1, 2, 3):
            forms.add(float(f"{p:.{sig}g}"))
        forms.update({round(p, 1), round(p, 2), round(p, 3), round(p, 4)})
    return forms


def check(item, summary=None) -> CheckReport:
    report = CheckReport()
    raw = _raw_gnomad(item)
    counts = _raw_counts(raw)

    # gnomAD counts on the page equal the raw response
    problems = []
    if raw:
        for v in (raw.get("response") or {}).get("variants") or []:
            var = (((v.get("totals") or {}).get("data") or {}).get("variant")) or {}
            for cs in ("exome", "genome"):
                d = var.get(cs)
                if d and d.get("an") and d.get("af") is not None and abs(d.get("ac", 0) / d["an"] - d["af"]) > 1e-9:
                    problems.append(f"gnomAD {cs} frequency {d['af']} is not allele count ÷ total ({d['ac']}/{d['an']})")
        for variant in getattr(item, "variants", None) or []:
            for cs in variant.callsets:
                key = "exome" if "exome" in cs.name else "genome"
                r = counts.get((variant.variant_id, key, "all"))
                if r and (cs.ac, cs.an) != r:
                    problems.append(f"{cs.name} shows {cs.ac}/{cs.an} but gnomAD returned {r[0]}/{r[1]}")
                for pop in cs.populations:
                    r = counts.get((variant.variant_id, key, pop.id))
                    if r and (pop.ac, pop.an) != r:
                        problems.append(f"{cs.name} {pop.id} shows {pop.ac}/{pop.an} but gnomAD returned {r[0]}/{r[1]}")
        report.add("gnomAD counts match the source", problems)

    # every percentage shown is computable from raw counts
    if raw and summary is not None:
        from snp_brief.demo import theme
        from snp_brief.summary import _plan

        forms = _percent_forms(counts)
        text = summary.text + " " + theme.facts_grid(_plan(item))
        problems = []
        for m in re.finditer(r"(\d+(?:\.\d+)?)%", text):
            value = float(m.group(1))
            if value not in forms:
                problems.append(f"{value}% in the summary cannot be computed from gnomAD's counts")
        report.add("Percentages trace back to allele counts", problems)

    # carrier rates
    from snp_brief.interpret import acmg_draft, genotype_rates

    problems = []
    raw_freqs = [ac / an for ac, an in counts.values()]
    for r in genotype_rates(item):
        if abs(r.carriers - (1 - (1 - r.af) ** 2)) > 1e-12 or abs(r.homozygous - r.af**2) > 1e-12:
            problems.append(f"carrier rate for {r.group} does not follow from its frequency")
        if raw and not any(abs(r.af - f) < 1e-12 for f in raw_freqs):
            problems.append(f"frequency {r.af:.6g} for {r.group} is not one gnomAD reported")
    if getattr(item, "variants", None):
        report.add("Carrier rates recomputed", problems)

    # position and alleles agree between dbSNP and gnomAD
    db = getattr(item, "db", None)
    if db is not None and db.pos and getattr(item, "variants", None):
        problems = []
        if all(abs(v.pos - db.pos) > 5 for v in item.variants):
            problems.append(f"dbSNP places it at {db.pos} but gnomAD at {', '.join(str(v.pos) for v in item.variants)}")
        alleles = {(ref, alt) for ref, alt in db.alleles}
        for v in item.variants:
            if alleles and (v.ref, v.alt) not in alleles and len(v.ref) == 1 and len(v.alt) == 1:
                problems.append(f"gnomAD allele {v.ref}>{v.alt} is not among dbSNP's alleles")
        report.add("dbSNP and gnomAD agree on position and alleles", problems)

    # scores in range
    problems = []
    sc = getattr(item, "scores", None)
    if sc is not None and sc.retrieved:
        for name, value in (("REVEL", sc.revel), ("AlphaMissense", sc.alphamissense), ("SIFT", sc.sift), ("PolyPhen-2", sc.polyphen)):
            if value is not None and not 0 <= value <= 1:
                problems.append(f"{name} {value} is outside 0 to 1")
        if sc.cadd_phred is not None and not 0 <= sc.cadd_phred <= 99:
            problems.append(f"CADD {sc.cadd_phred} is outside 0 to 99")
    sp = getattr(item, "splice", None)
    if sp is not None and sp.max_delta is not None and not 0 <= sp.max_delta <= 1:
        problems.append(f"SpliceAI {sp.max_delta} is outside 0 to 1")
    rg = getattr(item, "regulome", None)
    if rg is not None and rg.ranking and not re.fullmatch(r"[1-7][a-f]?", rg.ranking):
        problems.append(f"RegulomeDB rank '{rg.ranking}' is not a valid rank")
    if rg is not None and rg.probability is not None and not 0 <= rg.probability <= 1:
        problems.append(f"RegulomeDB probability {rg.probability} is outside 0 to 1")
    c = getattr(item, "constraint", None)
    if c is not None and c.loeuf is not None and c.loeuf < 0:
        problems.append(f"LOEUF {c.loeuf} is negative")
    if any(x is not None and x.retrieved for x in (sc, sp, rg, c)):
        report.add("Scores are within their published ranges", problems)

    # GWAS
    gwas = getattr(item, "gwas", None)
    if gwas is not None and gwas.hits:
        problems = []
        letters = {a for pair in (db.alleles if db else []) for a in pair}
        for hit in gwas.hits:
            underflow = hit.p_value == 0 and re.fullmatch(r"\d+(\.\d+)?[eE]-\d+", (hit.p_text or "").strip())
            if hit.p_value is not None and not underflow and not 0 < hit.p_value <= 1:
                problems.append(f"p-value {hit.p_text} for {hit.trait} is not between 0 and 1")
            if hit.risk_allele and letters and len(hit.risk_allele) == 1 and hit.risk_allele not in letters and hit.risk_allele in "ACGT":
                problems.append(f"risk allele {hit.risk_allele} for {hit.trait} is not an allele of this variant on either strand")
            if hit.pmid and not hit.pmid.isdigit():
                problems.append(f"PMID '{hit.pmid}' for {hit.trait} is not a number")
        report.add("GWAS values are valid", problems, warn=True)

    # drug and cancer levels, papers
    problems = []
    pgx = getattr(item, "pgx", None)
    for a in (pgx.annotations if pgx is not None else []):
        if a.level and a.level not in {"1A", "1B", "2A", "2B", "3", "4"}:
            problems.append(f"ClinPGx level '{a.level}' is not a ClinPGx evidence level")
    for it in (sc.civic if sc is not None else []):
        if it.level and it.level not in {"A", "B", "C", "D", "E"}:
            problems.append(f"CIViC level '{it.level}' is not a CIViC evidence level")
    t = getattr(item, "timeline", None)
    if t is not None and t.by_year:
        year = datetime.now(timezone.utc).year
        # Journals date issues ahead, so next year's papers appear before it starts; two years ahead is an error.
        if max(t.by_year) > year + 1 or min(t.by_year.values()) < 0:
            problems.append("papers-per-year has a year more than one ahead, or a negative count")
    if (pgx is not None and pgx.annotations) or (sc is not None and sc.civic) or (t is not None and t.by_year):
        report.add("Evidence levels and paper counts are valid", problems)

    # ACMG frequency codes vs ClinVar
    criteria = acmg_draft(item)
    codes = {x.code for x in criteria}
    from snp_brief.interpret import _clinvar_headline

    headline = _clinvar_headline(item).lower()
    problems = []
    if codes & {"BA1", "BS1"} and "pathogenic" in headline and "conflict" not in headline and "ClinVar" not in codes:
        problems.append("a frequency code says benign while ClinVar says pathogenic, with no warning")
    report.add("ACMG draft does not contradict ClinVar silently", problems)

    # genes named by different sources
    names: dict[str, set[str]] = {}
    if db is not None and db.genes:
        names["dbSNP"] = set(db.genes)
    clin = {g for r in getattr(item, "records", None) or [] for g in r.genes}
    if clin:
        names["ClinVar"] = clin
    if sc is not None and sc.gene:
        names["MyVariant.info"] = {sc.gene}
    vep = getattr(item, "vep", None)
    if vep is not None and getattr(vep, "canonical_gene", ""):
        names["Ensembl VEP"] = {vep.canonical_gene}
    if len(names) > 1:
        common = set.intersection(*names.values())
        problems = []
        if not common:
            problems.append(
                "sources name different genes, usually because the variant sits between or inside overlapping genes ("
                + "; ".join(f"{k}: {', '.join(sorted(v))}" for k, v in names.items())
                + ")"
            )
        report.add("Sources agree on the gene", problems, warn=True)
    return report
