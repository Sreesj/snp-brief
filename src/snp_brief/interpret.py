"""Computed interpretation aids: carrier and genotype rates, and an ACMG evidence draft.

Carrier rates assume Hardy-Weinberg equilibrium from gnomAD allele frequencies:
  at least one copy = 1 - (1 - p)^2,  two copies = p^2.
They are population estimates, not a person's risk.

The ACMG/AMP draft applies only the criteria that can be computed from public data, with the
thresholds ClinGen's Sequence Variant Interpretation group published:
  BA1   allele frequency above 5% in a gnomAD population (>= 2,000 alleles)
  BS1   above 1% (a generic cut-off; the right one depends on the disorder)
  PM2   absent from gnomAD or below 1 in 10,000 (applied at Supporting strength)
  PP3 / BP4  REVEL (Pejaver et al. 2022): >= 0.932 strong, >= 0.773 moderate, >= 0.644 supporting;
             <= 0.016 strong, <= 0.183 moderate, <= 0.290 supporting. For non-missense
             changes, SpliceAI >= 0.2 supports PP3 and <= 0.1 supports BP4.
  BP7   synonymous change with SpliceAI <= 0.1
  PVS1  flagged only as "to review": null variant in a gene intolerant of loss of function.
Every criterion carries its evidence. This is a starting point for a qualified reviewer, never
a classification.
"""

from __future__ import annotations

from dataclasses import dataclass

from snp_brief.sources.gnomad import ANCESTRY

MIN_POP_AN = 2000
NULL_TERMS = ("stop gained", "stop_gained", "frameshift", "splice acceptor", "splice donor", "splice_acceptor", "splice_donor", "start lost", "start_lost")


@dataclass
class GenotypeRate:
    group: str
    allele: str
    af: float
    carriers: float  # at least one copy
    homozygous: float
    an: int

    @property
    def carriers_text(self) -> str:
        return one_in(self.carriers)

    @property
    def homozygous_text(self) -> str:
        return one_in(self.homozygous)


def one_in(rate: float) -> str:
    """0.043 -> '1 in 23'; 0.00047 -> '1 in 2,100' (two significant figures)."""
    if rate <= 0:
        return "none observed"
    if rate >= 0.5:
        return f"{rate * 100:.0f}%"
    n = 1 / rate
    if n < 10:
        return f"1 in {n:.1f}"
    if n < 100:
        return f"1 in {round(n)}"
    return f"1 in {float(f'{n:.2g}'):,.0f}"


def genotype_rates(item) -> list[GenotypeRate]:
    """Overall and per-ancestry carrier and two-copy rates for each alternate allele."""
    rows: list[GenotypeRate] = []
    for variant in getattr(item, "variants", None) or []:
        allele = f"{variant.ref}>{variant.alt}"
        callset = max(variant.callsets, key=lambda c: c.an or 0, default=None)
        if callset is None or not callset.an:
            continue
        p = (callset.ac or 0) / callset.an
        rows.append(GenotypeRate(f"All ({callset.name})", allele, p, 1 - (1 - p) ** 2, p * p, callset.an))
        for pop in sorted(callset.populations, key=lambda x: -(x.ac / x.an if x.an else 0)):
            if not pop.an or pop.an < MIN_POP_AN or "_" in pop.id or pop.id.upper() in {"XX", "XY"}:
                continue  # ancestry groups only; sex splits and subgroups are left out
            q = pop.ac / pop.an
            rows.append(GenotypeRate(ANCESTRY.get(pop.id, pop.id), allele, q, 1 - (1 - q) ** 2, q * q, pop.an))
    return rows


@dataclass
class Criterion:
    code: str
    strength: str  # Stand-alone, Strong, Moderate, Supporting, To review
    direction: str  # pathogenic, benign, review
    evidence: str
    source: str


def acmg_draft(item) -> list[Criterion]:
    out: list[Criterion] = []
    variants = getattr(item, "variants", None) or []
    max_af, max_group = 0.0, ""
    seen_in_gnomad = False
    for variant in variants:
        for callset in variant.callsets:
            if callset.ac:
                seen_in_gnomad = True
            for pop in callset.populations:
                if pop.an and pop.an >= MIN_POP_AN and "_" not in pop.id and pop.id.upper() not in {"XX", "XY"}:
                    af = pop.ac / pop.an
                    if af > max_af:
                        max_af, max_group = af, f"{ANCESTRY.get(pop.id, pop.id)} ({callset.name})"
            if callset.af is not None and callset.af > max_af and callset.an and callset.an >= MIN_POP_AN:
                max_af, max_group = callset.af, f"all ({callset.name})"
    if variants:
        if max_af > 0.05:
            out.append(Criterion("BA1", "Stand-alone", "benign", f"{max_af * 100:.2g}% in {max_group}, above the 5% cut-off", "gnomAD"))
        elif max_af > 0.01:
            out.append(Criterion("BS1", "Strong", "benign", f"{max_af * 100:.2g}% in {max_group}; generic 1% cut-off, check the disorder's prevalence", "gnomAD"))
        elif not seen_in_gnomad or max_af < 0.0001:
            text = "absent from gnomAD" if not seen_in_gnomad else f"highest frequency {max_af * 100:.2g}% ({max_group}), below 1 in 10,000"
            out.append(Criterion("PM2", "Supporting", "pathogenic", text, "gnomAD"))
    elif getattr(item, "calls", {}).get("gnomad") and item.calls["gnomad"].origin in {"cache", "live"}:
        out.append(Criterion("PM2", "Supporting", "pathogenic", "absent from gnomAD", "gnomAD"))

    scores = getattr(item, "scores", None)
    splice = getattr(item, "splice", None)
    consequence = _consequence(item)
    revel = scores.revel if scores is not None else None
    splice_delta = splice.max_delta if splice is not None and splice.retrieved else None
    if revel is not None and "missense" in consequence:
        for cut, strength in ((0.932, "Strong"), (0.773, "Moderate"), (0.644, "Supporting")):
            if revel >= cut:
                out.append(Criterion("PP3", strength, "pathogenic", f"REVEL {revel:.3f} (>= {cut})", "MyVariant.info (dbNSFP)"))
                break
        else:
            for cut, strength in ((0.016, "Strong"), (0.183, "Moderate"), (0.290, "Supporting")):
                if revel <= cut:
                    out.append(Criterion("BP4", strength, "benign", f"REVEL {revel:.3f} (<= {cut})", "MyVariant.info (dbNSFP)"))
                    break
    elif splice_delta is not None:
        if splice_delta >= 0.2:
            out.append(Criterion("PP3", "Supporting", "pathogenic", f"SpliceAI {splice_delta:.2f} ({splice.kind}), >= 0.2", "SpliceAI"))
        elif splice_delta <= 0.1:
            if "synonymous" in consequence:
                out.append(Criterion("BP7", "Supporting", "benign", f"synonymous; SpliceAI {splice_delta:.2f}, <= 0.1", "SpliceAI"))
            else:
                out.append(Criterion("BP4", "Supporting", "benign", f"SpliceAI {splice_delta:.2f}, <= 0.1", "SpliceAI"))
    headline = _clinvar_headline(item)
    frequency_benign = [c.code for c in out if c.code in {"BA1", "BS1"}]
    if frequency_benign and "pathogenic" in headline.lower() and "conflict" not in headline.lower():
        out.append(
            Criterion(
                "ClinVar",
                "Check first",
                "review",
                f"ClinVar's best-reviewed classification is {headline}, which outweighs {' and '.join(frequency_benign)}. "
                "Common pathogenic alleles (founder effects, heterozygote advantage, reduced penetrance, as with HbS or "
                "Factor V Leiden) are known exceptions to the frequency rules; do not use these codes to downgrade it.",
                "ClinVar",
            )
        )
    elif frequency_benign and _trait_links(item):
        out.append(
            Criterion(
                "Context",
                "Check first",
                "review",
                f"{' and '.join(frequency_benign)} means it is too common to cause a rare inherited disease on its own. "
                f"It is still linked to {_trait_links(item)}, so it is not necessarily harmless.",
                "GWAS Catalog / ClinVar",
            )
        )
    constraint = getattr(item, "constraint", None)
    if any(term in consequence for term in NULL_TERMS):
        if constraint is not None and constraint.loeuf is not None and constraint.loeuf < 0.35:
            out.append(Criterion("PVS1", "To review", "review", f"{consequence} in {constraint.gene} (LOEUF {constraint.loeuf:.2f}); confirm loss of function is the disease mechanism", "gnomAD constraint"))
        else:
            out.append(Criterion("PVS1", "To review", "review", f"{consequence}; loss-of-function mechanism not established from constraint data", "VEP / MyVariant.info"))
    return out


def _trait_links(item) -> str:
    """What else the variant is linked to: GWAS traits, or ClinVar risk-factor/association/drug-response assertions."""
    kinds = []
    for record in getattr(item, "records", None) or []:
        for c in list(record.classifications) + list(record.submissions):
            text = (c.description or "").lower()
            for word in ("risk factor", "association", "drug response", "protective"):
                if word in text and word not in kinds:
                    kinds.append(word)
    gwas = getattr(item, "gwas", None)
    parts = []
    if kinds:
        parts.append("ClinVar " + ", ".join(kinds) + " assertions")
    if gwas is not None and getattr(gwas, "hits", None):
        traits = list(dict.fromkeys(h.trait for h in gwas.hits))[:2]
        parts.append("GWAS traits such as " + " and ".join(traits))
    return "; ".join(parts)


def _clinvar_headline(item) -> str:
    """ClinVar's best-reviewed classification, as the summary reports it."""
    if not getattr(item, "records", None):
        return ""
    try:
        from snp_brief.summary import _plan

        return _plan(item).headline or ""
    except Exception:
        return ""


def _consequence(item) -> str:
    vep = getattr(item, "vep", None)
    scores = getattr(item, "scores", None)
    text = ""
    if vep is not None and getattr(vep, "most_severe", ""):
        text = vep.most_severe
    elif scores is not None and scores.consequence:
        text = scores.consequence
    return text.replace("_", " ").lower()


NOT_COMPUTED = (
    "PS1/PM5 (same amino-acid change or codon already classified), PS3/BS3 (functional studies), "
    "PS4 (case-control), PP1/BS4 (segregation), PS2/PM6 (de novo) and PM3 (in trans) need case "
    "or literature review and are not computed."
)
