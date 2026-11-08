"""How each number or label that Variant Atlas computes itself is derived.

Shown on the page under "How this is computed" beside the result it explains, so a reviewer
can check the arithmetic and the thresholds. Data copied from a source is not listed here;
only what we compute from it.
"""

from __future__ import annotations

METHODS: dict[str, tuple[str, list[str]]] = {
    "zygosity": (
        "Alleles and zygosity",
        [
            "Each alternate allele at the rsID is listed separately (ref>alt and protein change), with ClinVar's classification for "
            "that allele only. ClinVar records are matched to alleles by their cDNA change.",
            "Each condition's mode of inheritance comes from the Human Phenotype Ontology disease annotations, looked up through the "
            "OMIM and Orphanet identifiers ClinVar attaches to the condition; if none matches, from the condition name when it "
            "states the mode (for example 'autosomal recessive …'). Otherwise it is shown as not stated.",
            "One copy / two copies follow from the mode: recessive → one copy is a carrier, two copies cause disease; dominant → "
            "one copy is enough; X-linked → males with one copy are affected, females usually carriers; semidominant → one copy "
            "milder, two more severe; multifactorial → contributes to risk.",
            "Copies observed and people with two copies come from gnomAD (the call set with the most alleles).",
            "Flags: a pathogenic allele whose inheritance is not stated; people with two copies in gnomAD for a recessive call "
            "(gnomAD excludes severe paediatric disease, so check penetrance); 10 or more homozygotes for a dominant call.",
        ],
    ),
    "summary": (
        "Summary",
        [
            "Every statement comes from a numbered list of facts taken from the sources; each fact keeps its source.",
            "Percentages are computed from the allele counts the source returned (count ÷ total alleles × 100); they are not new measurements.",
            "The writer (a fixed template, or the AI model when enabled) must cite fact numbers for every sentence.",
            "A checker then rejects any sentence that cites no fact, contains a number not found in the facts it cites, "
            "names a gene or rsID no fact contains, or reads like a diagnosis. Rejected AI sentences are replaced by the template.",
            "The source tags you see replace the fact numbers.",
        ],
    ),
    "reconcile": (
        "Reconciling sources",
        [
            "Each source has a reliability rank (1 = most dependable): ClinVar 1; gnomAD, ClinPGx and gene constraint 2; "
            "GWAS Catalog and GTEx 3; dbSNP and Ensembl LD 4; prediction tools 5; UniProt and NCBI Gene 6; literature 7.",
            "For each field the highest-ranked source that answered is used; lower-ranked answers are kept as supporting evidence.",
            "ClinVar classification uses the submissions at the highest review level present; flagged submissions and records that "
            "report no sequence change are excluded. When those top-level submissions disagree, the variant is shown as Conflicting.",
        ],
    ),
    "frequency": (
        "Frequency card",
        [
            "Allele frequency = alternate allele count ÷ total alleles, per gnomAD call set (exomes, genomes) and ancestry group.",
            "Highest and lowest groups use only groups with at least 1,000 alleles, so tiny samples do not dominate.",
            "Two figures are treated as different when one is at least twice the other, the gap is at least 0.1, "
            "or one is 0 and the other at least 0.001; otherwise one figure is shown.",
        ],
    ),
    "carriers": (
        "How common it is in people",
        [
            "p = allele frequency in the group (alternate alleles ÷ total alleles in gnomAD).",
            "Carry at least one copy = 1 − (1 − p)²; carry two copies = p². This assumes Hardy-Weinberg equilibrium "
            "(random mating, no strong selection), which is a standard approximation for population estimates.",
            "Shown as '1 in N' with N = 1 ÷ rate, rounded to two significant figures.",
            "The 'All' row uses the gnomAD call set with the most alleles; ancestry groups with fewer than 2,000 alleles and "
            "sex-specific groups are left out.",
            "These are population estimates, not an individual's risk.",
        ],
    ),
    "acmg": (
        "ACMG/AMP evidence draft",
        [
            "Only criteria computable from public data are applied, with ClinGen Sequence Variant Interpretation thresholds.",
            "BA1: allele frequency above 5% in any gnomAD ancestry group with at least 2,000 alleles.",
            "BS1: above 1% (a generic cut-off; the correct one depends on the disorder's prevalence and penetrance).",
            "PM2 (Supporting): absent from gnomAD, or highest group frequency below 1 in 10,000.",
            "PP3 / BP4 for missense changes use REVEL (Pejaver et al., 2022): ≥ 0.932 Strong, ≥ 0.773 Moderate, ≥ 0.644 Supporting; "
            "≤ 0.016 Strong, ≤ 0.183 Moderate, ≤ 0.290 Supporting benign. Between 0.290 and 0.644 no code applies.",
            "For non-missense changes, SpliceAI ≥ 0.2 supports PP3 and ≤ 0.1 supports BP4; a synonymous change with SpliceAI ≤ 0.1 gets BP7.",
            "PVS1 is only flagged for review on nonsense, frameshift, start-loss or canonical splice changes; LOEUF below 0.35 "
            "suggests loss of function may be the mechanism, which a reviewer must confirm.",
            "If BA1 or BS1 applies but ClinVar's best-reviewed classification is Pathogenic or Likely pathogenic, a 'Check first' "
            "line says the expert classification outweighs the frequency codes (common pathogenic alleles are known exceptions).",
            "If BA1 or BS1 applies and the variant has GWAS traits or ClinVar risk-factor, association, drug-response or protective "
            "assertions, a 'Context' line says that being common rules out a rare inherited disease, not every effect.",
            "Criteria needing case data or literature review (PS1, PM5, PS3, BS3, PS4, PP1, BS4, PS2, PM6, PM3) are not computed.",
        ],
    ),
    "scores": (
        "Score colours",
        [
            "Scores are copied from the tools; only the colour is ours.",
            "REVEL: red ≥ 0.644, green ≤ 0.290 (ClinGen calibration), amber between.",
            "AlphaMissense: colour follows its own class (likely pathogenic, ambiguous, likely benign).",
            "CADD: red at phred ≥ 20 (the top 1% most deleterious possible changes).",
            "SIFT: red below 0.05. PolyPhen-2 (HumDiv): red above 0.85, green below 0.15.",
            "SpliceAI: the largest of the four delta scores (acceptor/donor gain/loss); red ≥ 0.5, amber ≥ 0.2.",
            "RegulomeDB: red for ranks 1 and 2 (likely to affect binding).",
            "Gene constraint: LOEUF is the upper bound of the observed ÷ expected loss-of-function ratio; red below 0.35, amber below 0.6.",
            "When several alleles or transcripts are scored, the most damaging value is shown.",
        ],
    ),
    "gwas": (
        "GWAS associations",
        [
            "Associations come from the GWAS Catalog download (single-SNP associations only; haplotypes and interactions are left out).",
            "When a study reports the risk allele on the opposite DNA strand (for example A for a C>T variant), it is shown on the "
            "forward strand (T) with a note. Palindromic variants (A/T, C/G) are left as reported because the strand is ambiguous.",
            "Associations are ordered by p-value, read from the exact text so values smaller than 1e-308 still sort correctly.",
        ],
    ),
    "borzoi": (
        "Predicted regulatory effect (Borzoi)",
        [
            "Borzoi (Linder et al., Nature Genetics 2025) is a deep-learning model that predicts RNA-seq, CAGE, DNase/ATAC and ChIP "
            "coverage from 524,288 bases of DNA. We run it ourselves on NVIDIA GPUs; nothing here comes from a database.",
            "Input: the hg38 sequence centred on the variant, checked against the reference base, with the reference and the "
            "alternate allele. The four published model replicates are run on both, and on their reverse complements, and averaged.",
            "Model outputs are converted back to coverage by undoing Borzoi's target transforms (scale, soft clip, 3/4 power).",
            "Gene effect: RNA-seq coverage is summed over each gene's canonical exons (GENCODE v44) inside the 196,608-base output "
            "window; effect = log2((alt + 1) ÷ (ref + 1)) per track, averaged over tracks with the same tissue or cell type. "
            "'x% higher' = 2^|log2| − 1.",
            "Local effect: the same ratio for DNase/ATAC, ChIP and CAGE tracks over the 512 bases around the variant.",
            "Prediction vs measurement: for each GTEx eQTL of this variant, the sign of GTEx's normalized effect size is compared "
            "with Borzoi's predicted change for the same gene and the same alternate allele, in the closest-named tissue when one "
            "matches, else across all tissues. Only clear predictions (|log2| ≥ 0.005, about 0.35%) are counted; smaller "
            "ones are reported separately as 'no clear change', since a near-zero prediction cannot agree or disagree.",
            "Predictions are hypotheses about mechanism, not measurements; agreement with GTEx raises confidence, disagreement is a flag.",
        ],
    ),
    "timeline": (
        "Papers per year",
        [
            "Counts of Europe PMC records whose text matches the rsID exactly, grouped by publication year (up to 1,000 records).",
            "Trend compares the last three full years with the three before: growing if recent ≥ 6 papers and more than 1.5 × earlier; "
            "declining if earlier ≥ 6 and recent below 0.6 × earlier; otherwise steady.",
        ],
    ),
    "ripple": (
        "Molecular Ripple",
        [
            "The variant context comes from the consequence reported by ClinVar or VEP (for example missense, splice, intronic).",
            "Each link is labelled observed (reported by a source), predicted (a computational score or model) or unknown.",
            "Edges in the graph are styled by that label so a chain of predictions never looks like an established pathway.",
        ],
    ),
    "forge": (
        "Hypothesis Forge",
        [
            "Triage sorts the variant into a category from its ClinVar classification and review level, GWAS and eQTL evidence, "
            "and how many primary papers name it.",
            "Each supporting statement must pass the same fact check as the summary.",
            "Experiments are ranked by: consistency of the hypothesis with the facts, how many other hypotheses the result would "
            "separate, net support (statements for minus against), feasibility, then how many evidence gaps it addresses.",
        ],
    ),
}


def method_items(key: str) -> tuple[str, list[str]]:
    return METHODS[key]
