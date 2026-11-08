"""Plain-English summary written only from retrieved facts, then checked.

The evidence brief (brief.py) is complete but dense. This module turns the same
retrieved data into three or four readable sentences.

1. `build_facts` lists every fact the summary may use. Each fact has an id
   (F1, F2, ...), its text, and the source that reported it. Percentages are
   computed from the counts the API returned; they are not new measurements.
2. A writer produces sentences that end with the fact ids they rely on, e.g.
   "... is a missense variant in F5 [F2]." The writer is either the built-in
   template (deterministic, offline) or an LLM (llm.py).
3. `check` rejects any sentence that cites no fact, cites an unknown fact,
   contains a number that is not in the facts it cites, names a gene or rsID
   that no fact contains, or reads like a diagnosis.
4. `render` swaps the ids for source names, so the reader sees "(ClinVar)".

If the LLM output fails the check, the template summary is used instead.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from snp_brief.brief import _kept_traits, _views
from snp_brief.models import AlleleView, ClinVarRecord, GnomadVariant
from snp_brief.reconcile import OTHER, reconcile, side
from snp_brief.sources.gnomad import ANCESTRY, _freq, disagrees, percent, population_spread


@dataclass
class Fact:
    id: str
    text: str
    source: str


@dataclass
class Summary:
    rsid: str
    sentences: list[str]  # with [F#] tags
    facts: list[Fact]
    writer: str  # "template" or "llm:<model>"
    rejected: list[tuple[str, str]] = field(default_factory=list)  # (sentence, reason)

    @property
    def text(self) -> str:
        return render(self.sentences, self.facts)


# ---------------------------------------------------------------- facts


@dataclass
class _Plan:
    """Facts plus the ids the template needs, so it never guesses."""

    facts: list[Fact] = field(default_factory=list)
    function: str | None = None
    uniprot: str | None = None
    frequency: str | None = None
    callset_conflict: str | None = None
    spread: str | None = None
    disease: str | None = None
    disease_other: str | None = None
    dissent: str | None = None
    papers: list[str] = field(default_factory=list)
    # Plain values the template writes, each copied from the fact it cites.
    gene: str = ""
    consequence: str = ""
    change: str = ""
    alt: str = ""
    overall_pct: str = ""
    overall_set: str = ""
    high: tuple[str, str] | None = None
    low: tuple[str, str] | None = None
    conflict_pair: tuple[str, str] | None = None
    headline: str = ""
    conditions: list[str] = field(default_factory=list)
    unresolved_groups: list[str] = field(default_factory=list)
    other_labels: list[str] = field(default_factory=list)
    dissent_count: int = 0
    uniprot_text: str = ""
    aliases: str | None = None
    alias_names: str = ""
    location: str = ""
    location_fact: str | None = None
    ref_minor: tuple[str, str] | None = None
    predictions: str | None = None
    prediction_text: str = ""
    regulatory: str | None = None
    regulatory_text: str = ""
    gwas: list[str] = field(default_factory=list)
    gwas_traits: list[str] = field(default_factory=list)
    gwas_best_p: str = ""
    gwas_count: int = 0
    gwas_source: str = "GWAS Catalog"
    function_source: str = ""
    uniprot_gene: str = ""
    merged: str | None = None
    merged_text: str = ""
    gene_summary: str | None = None
    gene_summary_text: str = ""
    eqtl: str | None = None
    eqtl_genes: list[str] = field(default_factory=list)
    eqtl_text: str = ""
    ld: str | None = None
    ld_count: int = 0
    l2g: str | None = None
    l2g_text: str = ""
    litvar: str | None = None
    litvar_count: int = 0
    findings: list[str] = field(default_factory=list)
    finding_text: list[str] = field(default_factory=list)
    no_significance: bool = False
    paper_labels: list[str] = field(default_factory=list)


def build_facts(item) -> list[Fact]:
    return _plan(item).facts


def _plan(item) -> _Plan:
    plan = _Plan()

    def add(text: str, source: str) -> str:
        fact_id = f"F{len(plan.facts) + 1}"
        plan.facts.append(Fact(fact_id, text, source))
        return fact_id

    db = item.db
    merged_from = getattr(item, "merged_from", "")
    if merged_from:
        plan.merged = add(
            f"dbSNP has merged {merged_from} into {item.rsid}; ClinVar, gnomAD and the other sources list it under {item.rsid}.",
            "dbSNP",
        )
        plan.merged_text = f"{merged_from} is now {item.rsid} in dbSNP (the IDs were merged); results are for {item.rsid}"
    if db.retrieved and db.chrom and db.pos:
        alleles = ", ".join(f"{ref}>{alt}" for ref, alt in db.alleles)
        add(f"{item.rsid} is at {db.assembly or ''} chr{db.chrom}:{db.pos}, alleles {alleles}.".replace("  ", " "), "dbSNP")
    view = _main_view(item)
    if view is not None and view.records:
        record = view.records[0]
        plan.gene = ", ".join(record.genes) or ", ".join(db.genes)
        plan.consequence = ", ".join(record.consequences)
        inferred = ""
        if not plan.consequence:
            plan.consequence = _consequence_from_hgvs(record.cdna or record.title)
            inferred = " (consequence inferred from the HGVS name; ClinVar states none)" if plan.consequence else ""
        plan.change = _change(record)
        plan.function = add(
            f"{item.rsid} is {_phrase_consequence(plan.consequence)} in {plan.gene or 'an unnamed gene'}{inferred}"
            f"; {record.title or record.cdna} (ClinVar {record.accession}).",
            "ClinVar",
        )
        if record.aliases:
            plan.alias_names = record.aliases[0]
            plan.aliases = add(
                f"ClinVar lists other names for this same variant: {', '.join(record.aliases)}.", "ClinVar"
            )
        plan.function_source = "ClinVar"
    vep = getattr(item, "vep", None)
    if plan.function is None and vep is not None and vep.retrieved and (vep.canonical_terms or vep.most_severe):
        terms = [t.replace("_", " ") for t in (vep.canonical_terms or [vep.most_severe])]
        plan.consequence = ", ".join(terms)
        plan.gene = vep.canonical_gene or ", ".join(vep.genes) or ", ".join(db.genes)
        where = f" in {plan.gene}" if plan.gene else " with no overlapping transcript"
        extra = f"; most severe across transcripts: {vep.most_severe.replace('_', ' ')}" if vep.most_severe else ""
        plan.function = add(f"{item.rsid} is {_phrase_consequence(plan.consequence)}{where} (Ensembl VEP{extra}).", "Ensembl VEP")
        plan.function_source = "Ensembl VEP"
    if plan.function_source != "ClinVar" and db.fxn_class and db.genes:
        plan.gene = plan.gene or ", ".join(db.genes)
        plan.location = _location_from_fxn(db.fxn_class).replace("the gene", ", ".join(db.genes))
        plan.location_fact = add(
            f"dbSNP classifies {item.rsid} as {', '.join(c.replace('_', ' ') for c in db.fxn_class)} "
            f"relative to {', '.join(db.genes)}: {plan.location}.",
            "dbSNP",
        )
    if item.uniprot and item.fallback_gene:
        plan.gene = plan.gene or item.fallback_gene
        plan.uniprot_text = _first_sentence(_strip_refs(item.uniprot.text))
        plan.uniprot_gene = item.fallback_gene
        plan.uniprot = add(
            f"Reviewed {item.fallback_gene} protein function (describes the protein, not this variant): "
            f"{_strip_refs(item.uniprot.text)} (UniProt {item.uniprot.accession}).",
            "UniProt",
        )
    gene_summary = getattr(item, "gene_summary", None)
    if gene_summary is not None and (gene_summary.summary or gene_summary.name):
        text = gene_summary.summary or gene_summary.name
        plan.gene = plan.gene or gene_summary.symbol
        plan.gene_summary_text = _first_sentence(text)
        plan.gene_summary = add(
            f"NCBI Gene description of {gene_summary.symbol} (describes the gene, not this variant): "
            f"{gene_summary.name}{'. ' + text if gene_summary.summary else ''} (Gene ID {gene_summary.gene_id}).",
            "NCBI Gene",
        )
    if vep is not None and vep.retrieved:
        bits = []
        if vep.sift:
            bits.append(f"SIFT {vep.sift.replace('_', ' ')}")
        if vep.polyphen:
            bits.append(f"PolyPhen {vep.polyphen.replace('_', ' ')}")
        if vep.alphamissense:
            score = f" ({vep.alphamissense_score:g})" if vep.alphamissense_score is not None else ""
            bits.append(f"AlphaMissense {vep.alphamissense.replace('_', ' ')}{score}")
        if vep.cadd_phred is not None:
            bits.append(f"CADD PHRED {vep.cadd_phred:g}")
        if vep.spliceai_max is not None:
            bits.append(f"SpliceAI maximum delta score {vep.spliceai_max:g}")
        if bits:
            plan.prediction_text = "; ".join(bits)
            plan.predictions = add(
                f"Computational predictions from Ensembl VEP (not experimental evidence): {plan.prediction_text}.", "Ensembl VEP"
            )
        if vep.regulatory or vep.motifs:
            parts = []
            if vep.regulatory:
                parts.append("overlaps an annotated " + ", ".join(vep.regulatory[:3]))
            if vep.motifs:
                parts.append("lies in a predicted binding motif for " + ", ".join(vep.motifs[:3]))
            plan.regulatory_text = "; ".join(parts)
            plan.regulatory = add(
                f"Ensembl regulatory annotation: {item.rsid} {plan.regulatory_text}.", "Ensembl VEP"
            )
    gwas = getattr(item, "gwas", None)
    if gwas is not None and gwas.hits:
        plan.gwas_count = len(gwas.hits)
        for hit in gwas.hits:
            if hit.trait not in plan.gwas_traits:
                plan.gwas_traits.append(hit.trait)
        plan.gwas_best_p = gwas.hits[0].p_text
        for hit in gwas.hits[:3]:
            extra = [f"p = {hit.p_text}"] if hit.p_text else []
            if hit.risk_allele:
                extra.append(f"risk allele {hit.risk_allele}")
            if hit.effect:
                extra.append(hit.effect)
            if hit.genes:
                extra.append("reported genes " + ", ".join(hit.genes))
            if hit.pmid:
                extra.append(f"PMID:{hit.pmid}")
            plan.gwas.append(add(f"GWAS Catalog association with {hit.trait} ({'; '.join(extra)}).", "GWAS Catalog"))
        if plan.gwas_count > 3:
            plan.gwas.append(
                add(f"GWAS Catalog lists {plan.gwas_count} associations for {item.rsid} in total.", "GWAS Catalog")
            )
    eqtls = getattr(item, "eqtls", None) or []
    if eqtls:
        for eq in eqtls:
            if eq.gene not in plan.eqtl_genes:
                plan.eqtl_genes.append(eq.gene)
        top = eqtls[:3]
        bits = []
        for eq in top:
            direction = ""
            if eq.nes is not None:
                direction = ", alternate allele raises expression" if eq.nes > 0 else ", alternate allele lowers expression"
            p_text = f"p = {eq.p_value:.1e}" if eq.p_value is not None else ""
            bits.append(f"{eq.gene} in {eq.tissue} ({p_text}{direction})")
        plan.eqtl_text = "; ".join(bits)
        plan.eqtl = add(
            f"GTEx expression QTLs for {item.rsid}: {len(eqtls)} significant gene-tissue pairs across "
            f"{len(plan.eqtl_genes)} genes; strongest: {plan.eqtl_text}.",
            "GTEx",
        )
    l2g = getattr(item, "l2g", None) or []
    if l2g:
        top = l2g[:3]
        plan.l2g_text = "; ".join(f"{r.gene} (score {r.score:.2f}{', ' + r.trait if r.trait else ''})" for r in top)
        plan.l2g = add(f"Open Targets locus-to-gene predictions for loci including {item.rsid}: {plan.l2g_text}.", "Open Targets")
    ld = getattr(item, "ld_proxies", None) or []
    if ld:
        plan.ld_count = len(ld)
        shown = ", ".join(f"{p.rsid} (r² {p.r2:.2f})" for p in ld[:3])
        plan.ld = add(
            f"Ensembl LD (1000 Genomes, European): {len(ld)} nearby variants are inherited with {item.rsid} "
            f"at r² ≥ 0.8, e.g. {shown}.",
            "Ensembl LD",
        )
    litvar = getattr(item, "litvar", None) or []
    if litvar:
        plan.litvar_count = len(litvar)
        plan.litvar = add(
            f"LitVar2 finds {len(litvar)} publications mentioning {item.rsid} in their text, e.g. "
            + ", ".join(f"PMID:{p}" for p in litvar[:3])
            + ".",
            "LitVar2",
        )
    for finding in (getattr(item, "findings", None) or [])[:4]:
        detail = "; ".join(x for x in (finding.population, finding.effect) if x)
        plan.finding_text.append(f"{finding.claim.rstrip('.')} (PMID:{finding.pmid})")
        plan.findings.append(
            add(
                f"Finding from the abstract of PMID:{finding.pmid}: {finding.claim.rstrip('.')}"
                + (f" ({detail})" if detail else "")
                + f'. Quote: "{finding.quote}"',
                "Paper abstract",
            )
        )
    l2g_rows = getattr(item, "l2g", None) or []
    if not plan.gwas and l2g_rows:
        # Fallback when the GWAS Catalog did not answer: Open Targets' credible sets come from GWAS studies.
        traits = [t for t in dict.fromkeys(r.trait for r in l2g_rows) if t]
        if traits:
            plan.gwas_traits = traits
            plan.gwas_count = len(traits)
            plan.gwas_source = "Open Targets"
            plan.gwas.append(
                add(
                    f"Open Targets lists GWAS credible sets containing {item.rsid} for {len(traits)} traits, including "
                    + "; ".join(traits[:3]) + ".",
                    "Open Targets",
                )
            )
    variant = view.variant if view is not None else (item.variants[0] if item.variants else None)
    if variant is not None:
        _frequency_facts(plan, variant, add)
    if view is not None and view.records:
        _disease_facts(plan, view, add)
    elif item.clinvar_retrieved:
        plan.no_significance = True
        plan.disease = add(f"ClinVar reports no clinical significance for {item.rsid}.", "ClinVar")
    for paper in item.papers[:3]:
        first = paper.authors.split(",")[0].strip() or "Unknown"
        label = f"{first}{' et al.' if ',' in paper.authors else ''} {paper.year}".strip()
        plan.paper_labels.append(f"{label} (PMID:{paper.pmid})")
        plan.papers.append(add(f"Paper PMID:{paper.pmid}: {label}, {paper.title.rstrip('.')}.", "Europe PMC"))
    if not plan.papers and plan.litvar:
        plan.paper_labels = [f"PMID:{p}" for p in litvar[:3]]
        plan.papers.append(plan.litvar)
    citations = getattr(item, "citations", None) or []
    if not plan.papers and citations:
        shown = citations[:3]
        plan.paper_labels = [f"PMID:{pmid}" for pmid in shown]
        plan.papers.append(
            add(f"dbSNP links {len(citations)} PubMed records to {item.rsid}, e.g. {', '.join(plan.paper_labels)}.", "dbSNP citations")
        )
    if not plan.papers and vep is not None and vep.pubmed:
        shown = vep.pubmed[:3]
        plan.paper_labels = [f"PMID:{pmid}" for pmid in shown]
        plan.papers.append(
            add(f"Ensembl lists {len(vep.pubmed)} PubMed records citing this variant, e.g. {', '.join(plan.paper_labels)}.", "Ensembl VEP")
        )
    if not plan.no_significance and item.clinvar_retrieved and not item.has_significance:
        plan.no_significance = True
        plan.disease = plan.disease or add(f"ClinVar reports no clinical significance for {item.rsid}.", "ClinVar")
    return plan


def _main_view(item) -> AlleleView | None:
    views = _views(item.db, item.records, item.variants)
    for view in views:
        if view.records and any(record.has_significance for record in view.records):
            return view
    for view in views:
        if view.records:
            return view
    for view in views:
        if view.variant is not None:
            return view
    return None


def _change(record: ClinVarRecord) -> str:
    match = re.search(r"\((p\.[^)]+)\)", record.title or "")
    if match:
        return match.group(1)
    return record.cdna


def _consequence_from_hgvs(text: str) -> str:
    """A coding-DNA name such as c.*97G>A is in the 3' UTR. Only that unambiguous case is inferred."""
    match = re.search(r"c\.\*\d+[ACGT]>[ACGT]", text or "")
    return "3 prime UTR variant" if match else ""


def _location_from_fxn(classes: list[str]) -> str:
    joined = " ".join(classes).lower()
    if "upstream" in joined:
        return "upstream of the gene, not inside it"
    if "downstream" in joined:
        return "downstream of the gene, not inside it"
    if "intron" in joined:
        return "inside an intron of the gene"
    if "utr" in joined:
        return "in an untranslated region of the gene"
    if "missense" in joined or "coding" in joined or "synonymous" in joined:
        return "inside the coding sequence of the gene"
    return "near the gene"


def _strip_refs(text: str) -> str:
    """Drop UniProt's inline evidence lists such as "(PubMed:10500182, PubMed:12887909)"."""
    text = re.sub(r"\s*\((?:PubMed|ECO|By similarity)[^)]*\)", "", text)
    return re.sub(r"\s+", " ", text).strip()


def _first_sentence(text: str) -> str:
    head = re.split(r"(?<=[.;])\s", text.strip(), maxsplit=1)[0]
    return head.rstrip(".;")


def _frequency_facts(plan: _Plan, variant: GnomadVariant, add) -> None:
    plan.alt = f"{variant.ref}>{variant.alt}"
    usable = [call for call in variant.callsets if _freq(call) is not None]
    if not usable:
        return
    spread = population_spread(variant.callsets)
    if len(usable) >= 2:
        first, second = usable[0], usable[1]
        left, right = _freq(first), _freq(second)
        if left is not None and right is not None and disagrees(left, right):
            a = percent(first.ac or 0, first.an or 0) if first.an else f"{first.af}"
            b = percent(second.ac or 0, second.an or 0) if second.an else f"{second.af}"
            plan.conflict_pair = (f"{a} in {first.name}", f"{b} in {second.name}")
            plan.callset_conflict = add(
                f"The alternate allele {variant.alt} ({variant.ref}>{variant.alt}) has different frequencies in the two gnomAD call sets: "
                f"{first.name} {a} ({first.ac}/{first.an} alleles) and {second.name} {b} ({second.ac}/{second.an} alleles). "
                "The call sets sample different populations (exomes are mostly European ancestry), so this is not a disagreement between sources.",
                "gnomAD",
            )
    chosen = spread[0] if spread else max(usable, key=lambda call: call.an or 0)
    if chosen.an:
        plan.overall_pct = percent(chosen.ac or 0, chosen.an)
        plan.overall_set = chosen.name
        text = (
            f"The alternate allele {variant.alt} ({variant.ref}>{variant.alt}) is carried on {plan.overall_pct} of alleles in {chosen.name} "
            f"({chosen.ac}/{chosen.an} alleles; reported af {format(chosen.af, '.6g') if chosen.af is not None else 'n/a'})."
        )
        if (chosen.ac or 0) * 2 > chosen.an:
            plan.ref_minor = (variant.ref, percent(chosen.an - (chosen.ac or 0), chosen.an))
            text += f" So the reference allele {variant.ref} is the less common allele ({plan.ref_minor[1]})."
        plan.frequency = add(text, "gnomAD")
    if spread:
        pop_set, low, high = spread
        plan.high = (ANCESTRY[high.id], percent(high.ac, high.an))
        plan.low = (ANCESTRY[low.id], percent(low.ac, low.an))
        plan.spread = add(
            f"In {pop_set.name} the {variant.alt} allele frequency varies by ancestry: highest in {plan.high[0]} {plan.high[1]} "
            f"({high.ac}/{high.an} alleles), lowest in {plan.low[0]} {plan.low[1]} ({low.ac}/{low.an} alleles).",
            "gnomAD",
        )


def _disease_facts(plan: _Plan, view: AlleleView, add) -> None:
    verdict = reconcile(view.records, _kept_traits)
    if verdict is None:
        # Aggregate only, e.g. "association" with LACTASE PERSISTENCE.
        record = view.records[0]
        meaningful = [c for c in record.classifications if c.description and c.description.lower() != "not provided"]
        if not meaningful:
            plan.no_significance = True
            plan.disease = add("ClinVar reports no clinical significance for this allele.", "ClinVar")
            return
        top = meaningful[0]
        plan.headline = top.description
        plan.conditions = _kept_traits(top.traits)[:3]
        plan.disease = add(
            f"ClinVar classification: {top.description}"
            + (f", conditions {'; '.join(plan.conditions)}" if plan.conditions else "")
            + f" (ClinVar {record.accession}, {top.review_status}).",
            "ClinVar",
        )
        return
    if not verdict.resolved:
        for label, group in verdict.tallies:
            plan.unresolved_groups.append(f"{len(group)} {label}")
        plan.disease = add(
            f"ClinVar submissions at the highest review level present ({verdict.tier_review}) disagree: "
            + "; ".join(plan.unresolved_groups)
            + ".",
            "ClinVar",
        )
    elif verdict.headline:
        plan.headline = verdict.headline
        plan.conditions = verdict.traits[:3]
        plan.disease = add(
            f"Best-reviewed ClinVar classification: {verdict.headline}"
            + (f", conditions {'; '.join(plan.conditions)}" if plan.conditions else "")
            + f" ({'; '.join(verdict.headline_cites)}"
            + (f"; {verdict.support}" if verdict.support else "")
            + ").",
            "ClinVar",
        )
    if verdict.other:
        plan.other_labels = [label for label, _ in verdict.other[:3]]
        bits = []
        for label, group in verdict.other[:3]:
            traits = _kept_traits(group[0].classification.traits)[:1]
            bits.append(label + (f" ({traits[0]})" if traits else ""))
        plan.disease_other = add("ClinVar also records these assertions: " + "; ".join(bits) + ".", "ClinVar")
    if verdict.dissent:
        plan.dissent_count = len(verdict.dissent)
        labels = sorted({call.label for call in verdict.dissent})
        plan.dissent = add(
            f"{plan.dissent_count} lower-review ClinVar submissions disagree ({', '.join(labels)}).",
            "ClinVar",
        )


# ---------------------------------------------------------------- template writer


def template_summary(item) -> Summary:
    plan = _plan(item)
    return Summary(item.rsid, [sentence for _, sentence in _template_parts(item.rsid, plan)], plan.facts, "template")


def _template_parts(rsid: str, plan: _Plan) -> list[tuple[str, str]]:
    """Template sentences, each tagged with the field it covers (function, frequency, disease, papers)."""
    parts: list[tuple[str, str]] = []
    if plan.merged:
        parts.append(("identity", f"{plan.merged_text} [{plan.merged}]."))
    if plan.function:
        lead = f"{rsid} is {_phrase_consequence(plan.consequence)}" + (f" in {plan.gene}" if plan.gene else "")
        if plan.change:
            lead += f" ({plan.change})"
        ids = [plan.function] + ([plan.aliases] if plan.aliases else [])
        alias = f", also known as {plan.alias_names}" if plan.aliases and plan.alias_names else ""
        parts.append(("function", f"{lead}{alias} [{', '.join(ids)}]."))
    if plan.location_fact:
        parts.append(("function", f"dbSNP places {rsid} {plan.location} [{plan.location_fact}]."))
    if plan.uniprot:
        parts.append(
            (
                "function",
                f"UniProt describes the {plan.uniprot_gene or plan.gene or 'encoded'} protein, not this variant: "
                f"{plan.uniprot_text[:1].lower() + plan.uniprot_text[1:]} [{plan.uniprot}].",
            )
        )
    if plan.gene_summary and not plan.uniprot:
        parts.append(
            ("function", f"NCBI Gene describes {plan.gene} (the gene, not this variant): {plan.gene_summary_text} [{plan.gene_summary}].")
        )
    if plan.regulatory:
        parts.append(("function", f"It {plan.regulatory_text} [{plan.regulatory}]."))
    if plan.eqtl:
        genes = ", ".join(plan.eqtl_genes[:3])
        parts.append(("function", f"In GTEx it is linked to the expression of {genes} ({plan.eqtl_text}) [{plan.eqtl}]."))
    if plan.l2g:
        parts.append(("function", f"Open Targets' locus-to-gene model points to {plan.l2g_text} [{plan.l2g}]."))
    if plan.predictions:
        parts.append(("prediction", f"Computational predictions (not experimental): {plan.prediction_text} [{plan.predictions}]."))
    if plan.frequency:
        text = f"The alternate {plan.alt.split('>')[-1]} allele ({plan.alt}) is carried on {plan.overall_pct} of alleles in {plan.overall_set}"
        ids = [plan.frequency]
        if plan.spread and plan.high and plan.low:
            text += f", from {plan.high[1]} in {plan.high[0]} to {plan.low[1]} in {plan.low[0]}"
            ids.append(plan.spread)
        if plan.ref_minor:
            text += f"; the reference {plan.ref_minor[0]} allele is the less common one ({plan.ref_minor[1]})"
        parts.append(("frequency", f"{text} [{', '.join(ids)}]."))
    if plan.callset_conflict and plan.conflict_pair:
        parts.append(
            (
                "frequency",
                f"gnomAD exomes and genomes give different frequencies ({plan.conflict_pair[0]} vs {plan.conflict_pair[1]}) "
                f"because they sample different populations [{plan.callset_conflict}].",
            )
        )
    if plan.unresolved_groups and plan.disease:
        parts.append(
            (
                "disease",
                "Disease relevance is unresolved: ClinVar submissions at the same review level disagree ("
                + ", ".join(plan.unresolved_groups)
                + f") [{plan.disease}].",
            )
        )
    elif plan.headline and plan.disease:
        if side(plan.headline) == OTHER:
            cond = f" with {', '.join(plan.conditions)}" if plan.conditions else ""
            parts.append(("disease", f"ClinVar records it as {_article(plan.headline.lower())} {plan.headline.lower()}{cond}, not as a pathogenicity call [{plan.disease}]."))
        else:
            cond = f" for {', '.join(plan.conditions)}" if plan.conditions else ""
            parts.append(("disease", f"ClinVar's best-supported classification is {plan.headline}{cond} [{plan.disease}]."))
    elif plan.no_significance and plan.disease:
        parts.append(("disease", f"ClinVar reports no clinical significance for this variant [{plan.disease}]."))
    if plan.gwas:
        traits = ", ".join(plan.gwas_traits[:3])
        more = f" and {len(plan.gwas_traits) - 3} more traits" if len(plan.gwas_traits) > 3 else ""
        best = f" (strongest p = {plan.gwas_best_p})" if plan.gwas_best_p else ""
        parts.append(
            (
                "disease",
                f"Genome-wide association studies link it to {traits}{more}{best}; this is a statistical association, "
                f"not a clinical classification [{', '.join(plan.gwas)}].",
            )
        )
    extras = []
    extra_ids = []
    if plan.disease_other:
        extras.append("it also lists " + ", ".join(plan.other_labels) + " assertions")
        extra_ids.append(plan.disease_other)
    if plan.dissent:
        extras.append(f"{plan.dissent_count} lower-review submissions disagree")
        extra_ids.append(plan.dissent)
    if extras:
        parts.append(("disease_note", f"ClinVar note: {'; '.join(extras)} [{', '.join(extra_ids)}]."))
    if plan.findings:
        parts.append(("papers", f"Papers report: {'; '.join(plan.finding_text)} [{', '.join(plan.findings)}]."))
    if plan.papers:
        parts.append(("papers", f"Papers that cite this rsID: {'; '.join(plan.paper_labels)} [{', '.join(plan.papers)}]."))
    return parts


def field_ids(plan: _Plan) -> dict[str, set[str]]:
    """Fact ids that count as covering each required summary field."""
    return {
        "function": {i for i in (plan.function, plan.uniprot, plan.location_fact, plan.gene_summary, plan.eqtl) if i},
        "frequency": {i for i in (plan.frequency, plan.spread, plan.callset_conflict) if i},
        "disease": {i for i in (plan.disease, *plan.gwas) if i},
    }


def _article(word: str) -> str:
    return "an" if word[:1].lower() in "aeiou" else "a"


def _phrase_consequence(consequence: str) -> str:
    if not consequence:
        return "a variant"
    first = consequence.split(",")[0].strip().replace("_", " ")
    if first == "nonsense":
        first = "nonsense (stop-gained)"
    if "variant" not in first:
        first += " variant"
    return f"{_article(first)} {first}"


# ---------------------------------------------------------------- checker

_TAGS = re.compile(r"\[(F\d+(?:\s*,\s*F\d+)*)\]")
_NUMBER = re.compile(r"(?<![A-Za-z])\d+(?:[.,]\d+)*(?:e-?\d+)?%?")
_NAMED = re.compile(r"\b(rs\d+|[A-Z][A-Z0-9]{2,}[0-9A-Z]*)\b")
_ALLOWED_CAPS = {"SNP", "DNA", "RNA", "GRCH38", "VCV", "SCV", "PMID", "UTR", "NM", "LLM", "OMIM", "USA", "HGVS"}
_FORBIDDEN = ("diagnos", "you have", "your risk", "patient should", "we recommend", "treatment")


def check(sentences: list[str], facts: list[Fact]) -> tuple[list[str], list[tuple[str, str]]]:
    """Return (kept sentences, rejected (sentence, reason) pairs)."""
    by_id = {fact.id: fact for fact in facts}
    corpus = " ".join(fact.text for fact in facts)
    kept: list[str] = []
    rejected: list[tuple[str, str]] = []
    for sentence in sentences:
        reason = _problem(sentence, by_id, corpus)
        if reason:
            rejected.append((sentence, reason))
        else:
            kept.append(sentence)
    return kept, rejected


def _problem(sentence: str, by_id: dict[str, Fact], corpus: str) -> str | None:
    tags = _TAGS.findall(sentence)
    if not tags:
        return "cites no fact"
    ids = [part.strip() for tag in tags for part in tag.split(",")]
    unknown = [fact_id for fact_id in ids if fact_id not in by_id]
    if unknown:
        return f"cites unknown fact {', '.join(unknown)}"
    cited = " ".join(by_id[fact_id].text for fact_id in ids)
    body = _TAGS.sub("", sentence)
    cited_numbers = {_norm_number(n) for n in _NUMBER.findall(cited)}
    for number in _NUMBER.findall(body):
        if _norm_number(number) not in cited_numbers:
            return f"number {number} is not in the cited facts"
    for name in _NAMED.findall(body):
        if name.upper() in _ALLOWED_CAPS:
            continue
        if name not in corpus:
            return f"{name} is not in any fact"
    lowered = body.lower()
    cited_lower = cited.lower()
    for word in _FORBIDDEN:
        # A word copied from a source (e.g. a ClinVar condition "Cystic fibrosis diagnostic test") is not advice.
        if word in lowered and word not in cited_lower:
            return f"reads like clinical advice ({word.strip()})"
    return None


def _norm_number(text: str) -> str:
    return text.replace(",", "").rstrip("%").rstrip(".")


def render(sentences: list[str], facts: list[Fact]) -> str:
    """Replace [F#] tags with the names of the sources that reported them."""
    by_id = {fact.id: fact for fact in facts}

    def swap(match: re.Match) -> str:
        names: list[str] = []
        for part in match.group(1).split(","):
            fact = by_id.get(part.strip())
            if fact and fact.source not in names:
                names.append(fact.source)
        return f"({', '.join(names)})" if names else ""

    out = []
    for sentence in sentences:
        text = _TAGS.sub(swap, sentence)
        out.append(re.sub(r"\s+([.,;])", r"\1", text).strip())
    return " ".join(out)


def facts_block(facts: list[Fact]) -> str:
    from snp_brief.providers import tier_label

    return "\n".join(f"[{fact.id}] ({fact.source}; {tier_label(fact.source)}) {fact.text}" for fact in facts)


# ---------------------------------------------------------------- entry point


def summarize(item, mode: str = "template", client=None, cache_dir=None, refresh: bool = False) -> Summary:
    """mode is "template" or "llm". LLM output that fails the check falls back to the template.

    With cache_dir, an LLM summary is saved to cache/<rsid>/summary_llm.json and reused while
    the facts are unchanged. A cached LLM summary is used even without an API key, so a
    pre-generated demo runs offline.
    """
    from snp_brief.module_cache import result_key, save_result

    fallback = template_summary(item)
    kept, rejected = check(fallback.sentences, fallback.facts)
    fallback.sentences, fallback.rejected = kept, rejected
    if mode != "llm":
        return fallback
    from snp_brief.llm import LLMClient, write_summary

    client = client or LLMClient.from_env()
    block = facts_block(fallback.facts)
    if not refresh:
        cached = _cached_llm_summary(cache_dir, item.rsid, client, block)
        if cached is not None:
            kept, rejected = check(cached["sentences"], fallback.facts)
            if kept and len(rejected) <= len(kept):
                return Summary(item.rsid, kept, fallback.facts, f"llm:{cached['model']} (cached)", rejected)
    if client is None:
        return fallback
    try:
        sentences = write_summary(client, item.rsid, fallback.facts)
    except Exception as exc:  # network, timeout, bad response
        fallback.rejected.append(("(LLM call)", f"failed: {exc}"))
        return fallback
    kept, rejected = check(sentences, fallback.facts)
    if not kept or len(rejected) > len(kept):
        fallback.rejected.extend(rejected)
        return fallback
    kept, filled = fill_missing_fields(item, kept)
    writer = f"llm:{client.model}" + (f" + template for {', '.join(filled)}" if filled else "")
    save_result(
        cache_dir,
        item.rsid,
        "summary_llm",
        result_key("summary_llm", client.model, block),
        {"model": client.model, "facts_hash": _hash(block), "sentences": kept, "filled": filled},
    )
    return Summary(item.rsid, kept, fallback.facts, writer, rejected)


def fill_missing_fields(item, sentences: list[str]) -> tuple[list[str], list[str]]:
    """Every summary must cover function, frequency, and disease relevance when facts exist.

    An LLM summary that dropped a field (or had that sentence rejected by the
    checker) gets the template sentence for that field, in reading order.
    """
    plan = _plan(item)
    required = field_ids(plan)
    parts = _template_parts(item.rsid, plan)
    cited = [set(_ids(sentence)) for sentence in sentences]
    covered = {name for name, ids in required.items() if ids and any(c & ids for c in cited)}
    missing = [name for name, ids in required.items() if ids and name not in covered]
    if not missing:
        return sentences, []
    order = {"function": 0, "frequency": 1, "disease": 2}
    rank_of = []
    for ids in cited:
        positions = [order[name] for name, req in required.items() if ids & req]
        rank_of.append(min(positions) if positions else 3)
    out = list(sentences)
    for name in sorted(missing, key=order.get):
        extra = [sentence for field, sentence in parts if field == name or (name == "disease" and field == "disease_note")]
        where = next((i for i, r in enumerate(rank_of) if r > order[name]), len(out))
        out[where:where] = extra
        rank_of[where:where] = [order[name]] * len(extra)
    return out, missing


def _ids(sentence: str) -> list[str]:
    return [part.strip() for tag in _TAGS.findall(sentence) for part in tag.split(",")]


def _hash(text: str) -> str:
    import hashlib

    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def _cached_llm_summary(cache_dir, rsid: str, client, block: str) -> dict | None:
    """A stored LLM summary for these exact facts: same model when a client is set, any model offline."""
    import json
    from pathlib import Path

    if cache_dir is None:
        return None
    path = Path(cache_dir) / rsid / "summary_llm.json"
    if not path.is_file():
        return None
    try:
        stored = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if stored.get("facts_hash") != _hash(block):
        return None
    if client is not None and stored.get("model") != client.model:
        return None
    return stored
