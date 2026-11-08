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
        ("utr", ("utr", "upstream", "downstream", "5 prime", "3 prime")),
    ):
        if any(word in text for word in words):
            return key
    return "intergenic"


def ripple(item) -> RippleResult:
    plan = _plan(item)
    view = _main_view(item)
    record = view.records[0] if view is not None and view.records else None
    key = context_key(plan.consequence)
    name, priorities = CONTEXTS[key]
    result = RippleResult(item.rsid, key, name, plan.gene, plan.change, priorities)
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
    result.links.append(
        Link(PREDICTED, "Pathogenicity and splicing scores (AlphaMissense, SpliceAI, CADD) — not yet retrieved; next step is Ensembl VEP.")
    )
    if key in {"intron", "intergenic", "utr", "synonymous"}:
        result.links.append(Link(UNKNOWN, f"Which gene this variant regulates; the gene named by ClinVar ({plan.gene or 'none'}) is the nearest annotation, not proof of the target."))
        result.links.append(Link(UNKNOWN, "Whether this SNP is causal or tags a nearby causal variant."))
    elif key == "missense":
        result.links.append(Link(UNKNOWN, "How much the amino-acid change alters protein activity in this allele's studies."))
    elif key == "nonsense":
        result.links.append(Link(UNKNOWN, "Whether the truncated product is made or degraded."))
    return result
