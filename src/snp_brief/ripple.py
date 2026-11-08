"""Molecular Ripple, first version: the variant's molecular context and what to examine.

This version uses only evidence the pipeline already has. It classifies the
variant context from the ClinVar molecular consequence and lists what a
mechanistic analysis should prioritize for that context. Every statement is
labelled:

- observed: reported by a source (ClinVar, dbSNP, UniProt, gnomAD)
- predicted: a computational prediction (none yet; Ensembl VEP scores such as
  AlphaMissense, SpliceAI and CADD are the next step)
- unknown: a question the evidence does not answer

The labels exist so a chain of predictions never looks like an established pathway.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from snp_brief.summary import _main_view, _plan

OBSERVED = "observed"
PREDICTED = "predicted"
UNKNOWN = "unknown"

CONTEXTS = {
    "missense": (
        "Missense coding variant",
        ["Amino-acid change and its chemistry", "Protein domain at this position", "Protein structure", "Functional studies of this allele"],
    ),
    "nonsense": (
        "Stop-gained (nonsense) variant",
        ["Whether the truncated transcript escapes nonsense-mediated decay", "Loss of protein function", "Whether loss of function is tolerated (population frequency)"],
    ),
    "synonymous": (
        "Synonymous variant",
        ["Splicing effects", "Transcript regulation", "Linkage to a nearby causal variant"],
    ),
    "splice": (
        "Splice-site variant",
        ["Transcript consequences", "Splice predictions (SpliceAI)", "RNA evidence of altered splicing"],
    ),
    "intron": (
        "Intronic variant",
        ["Splice-related predictions", "Known regulatory elements (enhancers) at this position", "Expression QTLs and the gene they act on"],
    ),
    "utr": (
        "UTR / near-gene variant",
        ["mRNA stability and translation", "Regulatory elements", "Expression QTLs"],
    ),
    "intergenic": (
        "Intergenic or unannotated variant",
        ["Regulatory elements", "Candidate target genes", "Expression QTLs and fine-mapping", "Relevant tissue context"],
    ),
}


@dataclass
class Link:
    label: str  # observed, predicted, unknown
    text: str
    source: str = ""


@dataclass
class RippleResult:
    rsid: str
    context_key: str
    context: str
    gene: str
    change: str
    priorities: list[str]
    links: list[Link] = field(default_factory=list)


def context_key(consequence: str) -> str:
    text = consequence.lower()
    for key, words in (
        ("nonsense", ("nonsense", "stop gained", "stop_gained")),
        ("splice", ("splice donor", "splice acceptor", "splice site", "splice_")),
        ("missense", ("missense",)),
        ("synonymous", ("synonymous",)),
        ("intron", ("intron",)),
        ("utr", ("utr", "upstream", "downstream", "5 prime", "3 prime", "2kb", "500b")),
    ):
        if any(word in text for word in words):
            return key
    return "intergenic"


def ripple(item) -> RippleResult:
    plan = _plan(item)
    view = _main_view(item)
    record = view.records[0] if view is not None and view.records else None
    db_classes = ", ".join(c.replace("_", " ") for c in item.db.fxn_class)
    vep = getattr(item, "vep", None)
    # Fallback by rank: ClinVar consequence, then Ensembl VEP (already in plan.consequence), then dbSNP class.
    consequence = plan.consequence or db_classes
    if record is not None and record.genes:
        gene_source = "ClinVar"
    elif plan.function_source == "Ensembl VEP":
        gene_source = "Ensembl VEP"
    else:
        gene_source = "dbSNP"
    key = context_key(consequence)
    name, priorities = CONTEXTS[key]
    result = RippleResult(item.rsid, key, name, plan.gene, plan.change, priorities)
    if record is None and db_classes:
        result.links.append(
            Link(OBSERVED, f"{item.rsid} → {plan.gene or 'no gene named'}: {db_classes}", "dbSNP")
        )
    if record is not None:
        result.links.append(
            Link(
                OBSERVED,
                f"{item.rsid} → {plan.gene or 'gene not named'}: {plan.consequence or 'no consequence stated'}"
                + (f", {plan.change}" if plan.change else ""),
                f"ClinVar {record.accession}",
            )
        )
    if item.uniprot is not None and item.fallback_gene:
        result.links.append(
            Link(OBSERVED, f"{item.fallback_gene} protein function: {item.uniprot.text}", f"UniProt {item.uniprot.accession}")
        )
    if plan.headline:
        cond = f" ({', '.join(plan.conditions)})" if plan.conditions else ""
        result.links.append(Link(OBSERVED, f"Clinical classification: {plan.headline}{cond}", "ClinVar"))
    if plan.high and plan.low:
        result.links.append(
            Link(OBSERVED, f"Frequency varies by ancestry: {plan.high[1]} in {plan.high[0]}, {plan.low[1]} in {plan.low[0]}", "gnomAD")
        )
    if vep is not None and vep.retrieved and plan.function_source == "Ensembl VEP" and record is None:
        result.links.append(
            Link(OBSERVED, f"{item.rsid} → {plan.gene or 'no overlapping gene'}: {plan.consequence}", "Ensembl VEP annotation")
        )
    if vep is not None and vep.regulatory:
        result.links.append(
            Link(OBSERVED, "Overlaps an annotated " + ", ".join(vep.regulatory[:3]), "Ensembl Regulatory Build")
        )
    gene_summary = getattr(item, "gene_summary", None)
    if gene_summary is not None and item.uniprot is None and (gene_summary.summary or gene_summary.name):
        result.links.append(
            Link(OBSERVED, f"{gene_summary.symbol}: {gene_summary.name}", f"NCBI Gene {gene_summary.gene_id}")
        )
    for eq in (getattr(item, "eqtls", None) or [])[:3]:
        direction = "" if eq.nes is None else (" (alternate allele raises expression)" if eq.nes > 0 else " (alternate allele lowers expression)")
        result.links.append(Link(OBSERVED, f"{item.rsid} → expression of {eq.gene} in {eq.tissue}{direction}", "GTEx eQTL (statistical)"))
    ld = getattr(item, "ld_proxies", None) or []
    if ld:
        result.links.append(
            Link(OBSERVED, f"Inherited together with {len(ld)} nearby variants (r² ≥ 0.8), e.g. " + ", ".join(p.rsid for p in ld[:3]), "Ensembl LD")
        )
    scores = []
    for row in (getattr(item, "l2g", None) or [])[:2]:
        scores.append(Link(PREDICTED, f"Locus-to-gene model favours {row.gene} (score {row.score:.2f})" + (f" for {row.trait}" if row.trait else ""), "Open Targets"))
    if vep is not None and vep.retrieved:
        if vep.motifs:
            scores.append(Link(PREDICTED, "Inside a predicted binding motif for " + ", ".join(vep.motifs[:3]), "Ensembl VEP"))
        if vep.alphamissense:
            value = f" ({vep.alphamissense_score:g})" if vep.alphamissense_score is not None else ""
            scores.append(Link(PREDICTED, f"AlphaMissense: {vep.alphamissense.replace('_', ' ')}{value}", "Ensembl VEP"))
        if vep.sift or vep.polyphen:
            scores.append(
                Link(PREDICTED, "; ".join(x for x in (f"SIFT {vep.sift}" if vep.sift else "", f"PolyPhen {vep.polyphen}" if vep.polyphen else "") if x).replace("_", " "), "Ensembl VEP")
            )
        if vep.cadd_phred is not None:
            scores.append(Link(PREDICTED, f"CADD PHRED {vep.cadd_phred:g} (higher = more likely deleterious; 20 is the top 1%)", "Ensembl VEP"))
        if vep.spliceai_max is not None:
            scores.append(Link(PREDICTED, f"SpliceAI maximum delta score {vep.spliceai_max:g} (above 0.5 suggests a splicing effect)", "Ensembl VEP"))
    if scores:
        result.links.extend(scores)
    else:
        result.links.append(
            Link(PREDICTED, "No prediction scores were returned (Ensembl VEP not answered, or scores not available for this variant).")
        )
    if key in {"intron", "intergenic", "utr", "synonymous"}:
        result.links.append(
            Link(
                UNKNOWN,
                f"Which gene this variant regulates; the gene named by {gene_source} ({plan.gene or 'none'}) "
                "is the nearest annotation, not proof of the target.",
            )
        )
        result.links.append(Link(UNKNOWN, "Whether this SNP is causal or tags a nearby causal variant."))
    elif key == "missense":
        result.links.append(Link(UNKNOWN, "How much the amino-acid change alters protein activity in this allele's studies."))
    elif key == "nonsense":
        result.links.append(Link(UNKNOWN, "Whether the truncated product is made or degraded."))
    return result
