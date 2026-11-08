# Variant Atlas

**Know which allele matters, whether one copy or two cause disease, and how sure anyone really is, before you order primers.**

Variant Atlas turns any list of SNP IDs (rsIDs) into a sourced, plain-text answer on biological function, population frequency and disease relevance. It reconciles 18 sources, shows each allele and its zygosity first, and re-derives every number from the raw data each time a variant loads.

Built at Nucleate Pittsburgh BioHack 2026 (QIAGEN BioTech track). It is a research tool and does not diagnose.

**Live demo:** https://variant-atlas.streamlit.app (see [Deploy](#deploy-always-on-demo) if the link is not up yet)

---

## Contents

- [What it does](#what-it-does)
- [Run it](#run-it)
- [Data sources](#data-sources)
- [How we compute](#how-we-compute)
- [Borzoi on NVIDIA GPUs](#borzoi-on-nvidia-gpus)
- [Data checks and validation](#data-checks-and-validation)
- [Engineering](#engineering)
- [Deploy (always-on demo)](#deploy-always-on-demo)
- [Repository layout](#repository-layout)
- [Limitations](#limitations)

## What it does

One search box takes one rsID or many (commas, spaces or new lines). Each variant opens with:

1. **Alleles and zygosity, check this first.** Every alternate allele separately (ref>alt and protein change), its own ClinVar classification and conditions, each condition's mode of inheritance, what one copy and two copies mean, gnomAD copies and homozygotes, and warnings.
2. **Variant Observatory.** Function, frequency and clinical cards; a plain-text summary with a source tag on every sentence; the data-checks result; population frequencies by ancestry; literature (papers per year, findings quoted from abstracts); and the status of every source.
3. **Interpretation.** An ACMG/AMP evidence draft with "Check first" safety lines, carrier and two-copy rates by ancestry, drug response (ClinPGx) and cancer evidence (CIViC).
4. **Molecular Ripple.** Prediction and constraint scores, our Borzoi expression predictions compared with measured GTEx eQTLs, and molecular links labelled observed, predicted or unknown, with an evidence graph.
5. **Hypothesis Forge.** For variants that are not fully understood: competing explanations, evidence for and against, and ranked experiments that would tell them apart. Well-characterized and thin-evidence variants get no invented hypotheses.
6. **Notes and export.** Notes per variant and a full report for one variant or the whole list.

Every result we calculate has a **"How this is computed"** note beside it with the formula and thresholds.

## Run it

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then:

```bash
git clone https://github.com/Sreesj/snp-brief.git
cd snp-brief
uv run streamlit run streamlit_app.py      # the web app
uv run snp-brief                           # plain-text briefs for fixture/snps.txt
```

The 20 demo variants in `fixture/snps.txt` are cached in `cache/`, so the demo runs instantly and offline. Any other rsID is fetched live and saved to `~/.cache/snp-brief`; sources that are new since a variant was cached are fetched on top of the cached ones.

Useful commands:

| Command | What it does |
| --- | --- |
| `uv run snp-brief --input ids.txt --save out.txt` | Briefs for any list |
| `uv run snp-brief --check rs334 rs6025` | Load variants (cache first, live otherwise) and run the data checks |
| `uv run snp-brief --probe rs4244285` | Call every source live and print its status and response shape |
| `uv run snp-brief --list-models` | Models your AI provider offers |
| `uv run snp-brief --build-gwas FILE.tsv` | Slim the GWAS Catalog download into `data/gwas_catalog.tsv.gz` |
| `uv run snp-brief --build-inheritance phenotype.hpoa` | Build `data/inheritance.tsv.gz` from HPO |
| `uv run snp-brief --refresh-gnomad` | Re-fetch gnomAD (with homozygote counts) for the demo list |
| `uv run snp-brief --export-variants brev/variants.tsv` | Input file for the Borzoi GPU run |
| `uv run snp-brief --import-borzoi brev/out` | Copy Borzoi results into `cache/` |
| `PYTHONPATH=src python3 scripts/verify_cached.py` | Independent cross-check of the cached demo variants |

Optional settings (environment variables or Streamlit secrets; see `.streamlit/secrets.toml.example`): `NEBIUS_API_KEY` for AI-written summaries, `SNP_BRIEF_LLM_MODEL`, `NCBI_API_KEY`, `SNP_BRIEF_APP_PASSWORD`, `SNP_BRIEF_ENABLE=regulome,spliceai` to switch on the two optional services. Without an AI key the app uses its deterministic template summaries.

## Data sources

All sources are queried in parallel. Each has a reliability rank (1 = most dependable); for every field the highest-ranked source that answered wins and lower-ranked answers are kept as support.

| Rank | Source | What we use |
| --- | --- | --- |
| 1 | ClinVar | Classification per allele, review status, conditions with OMIM/Orphanet IDs, protein change |
| 2 | gnomAD v4 | Allele counts and frequencies by ancestry, homozygote and hemizygote counts |
| 2 | ClinPGx (formerly PharmGKB) | Drug-response annotations, evidence levels 1A to 4 |
| 2 | gnomAD gene constraint | LOEUF and pLI |
| 3 | GWAS Catalog | Traits, p-values, risk alleles, PMIDs (official download file, indexed locally) |
| 3 | GTEx v8 | eQTLs: gene, tissue, normalized effect size, p-value |
| 4 | dbSNP | GRCh38 position, alleles, merged rsIDs |
| 4 | Ensembl LD | Variants inherited together (r² ≥ 0.8) |
| 5 | MyVariant.info | REVEL, AlphaMissense, CADD, SIFT, PolyPhen-2, consequence, CIViC |
| 5 | Ensembl VEP | Consequence, regulatory overlap, protein domains (backup) |
| 5 | Open Targets | Locus-to-gene scores and GWAS credible sets |
| 5 | Borzoi (our model run) | Predicted change in expression and chromatin |
| 6 | UniProt, NCBI Gene | Gene-level function, only when ClinVar has no significance; never used for disease |
| 7 | Europe PMC, dbSNP citations, LitVar2 | Papers naming the rsID; papers per year |

Supporting data: Human Phenotype Ontology disease annotations (8,738 disease-inheritance rows) and GENCODE v44 for Borzoi. RegulomeDB and the public SpliceAI lookup server failed on every live check, so they are off by default.

## How we compute

### Reconciling ClinVar

1. Flagged submissions are excluded and counted.
2. Each submission is put on a side: pathogenic, benign or uncertain.
3. The headline is what the highest review-star tier present agrees on (practice guideline 4, expert panel 3, multiple submitters 2, single submitter 1, no criteria 0). If that tier disagrees across sides, the field is **Conflicting** with the tally.
4. Lower-tier submissions on another side are kept as named dissent.
5. Risk factor, drug response, protective and association are listed as additional assertions, never as a conflict.
6. ClinVar records are matched to alleles by their cDNA change, so every call belongs to one allele. Merged rsIDs are followed.

### Frequencies and carrier rates

```
p = AC / AN                      (allele count / total alleles, from gnomAD)
at least one copy = 1 - (1 - p)^2
two copies        = p^2          (Hardy-Weinberg; shown as "1 in N", N rounded to 2 significant figures)
```

Highest and lowest ancestry groups use groups with ≥ 1,000 alleles; carrier tables use groups with ≥ 2,000 alleles. Example, Factor V Leiden (rs6025) in gnomAD v4 exomes: 31,939 / 1,460,848 = 2.19%; one copy 4.32% (about 1 in 23); two copies 0.048% (about 1 in 2,100).

### Alleles and zygosity

Each condition's mode of inheritance comes from HPO disease annotations, looked up through the OMIM and Orphanet IDs that ClinVar attaches to the condition; failing that, from the condition name if it states the mode. Otherwise it is shown as "inheritance not stated", never guessed. Zygosity is stated per condition, never merged across conditions.

| Mode | One copy | Two copies |
| --- | --- | --- |
| Autosomal recessive | Carrier, usually unaffected | Causes disease |
| Autosomal dominant | Enough to cause disease | Causes disease, often more severe |
| Semidominant | Milder disease or risk | More severe disease |
| X-linked | Males affected; females usually carriers | Females with two copies affected |

Flags: pathogenic with inheritance not stated; conditions on one record that differ in inheritance; homozygotes in gnomAD for a recessive pathogenic call (gnomAD excludes severe paediatric disease); 10 or more homozygotes for a dominant pathogenic call.

### ACMG/AMP evidence draft

Only criteria computable from public data, with ClinGen thresholds. A draft for a qualified reviewer, never a classification.

| Code | Rule |
| --- | --- |
| BA1 | Allele frequency > 5% in a gnomAD ancestry group with ≥ 2,000 alleles |
| BS1 | > 1% (generic cut-off) |
| PM2 (supporting) | Absent from gnomAD or < 1 in 10,000 |
| PP3 | Missense REVEL ≥ 0.932 strong, ≥ 0.773 moderate, ≥ 0.644 supporting (Pejaver et al. 2022); SpliceAI ≥ 0.2 for non-missense |
| BP4 | REVEL ≤ 0.016 / ≤ 0.183 / ≤ 0.290; SpliceAI ≤ 0.1 |
| BP7 | Synonymous with SpliceAI ≤ 0.1 |
| PVS1 | To review only: null variant; LOEUF < 0.35 suggests loss of function may be the mechanism |

Safety lines: if BA1/BS1 applies but ClinVar's best-reviewed call is pathogenic, ClinVar wins (sickle cell, Factor V Leiden, prothrombin G20210A, CFTR F508del); if a common variant has GWAS traits or risk-factor assertions, it is flagged as not necessarily harmless.

### Summary and fact check

1. Every usable fact becomes a numbered entry with its source.
2. A fixed template, or the AI model (open-weight `openai/gpt-oss-120b` on Nebius Token Factory), writes sentences that must cite the fact numbers.
3. A sentence is rejected if it cites no fact, contains a number not in its facts, names a gene or rsID no fact contains, or reads like medical advice.
4. Rejected sentences fall back to the template; fact numbers become source tags.

The clinical sentence names the allele ("For the T>A (p.Glu7Val) allele …"), other alleles are named with their own calls, and GWAS risk alleles reported on the opposite strand are shown on the forward strand with a note.

### Literature

Findings extracted from abstracts are kept only if the PMID is one we supplied, the quote appears verbatim in the abstract, and every number in the claim is in the quote. Papers per year come from Europe PMC; the trend compares the last three full years with the three before.

## Borzoi on NVIDIA GPUs

We run [Borzoi](https://www.nature.com/articles/s41588-024-02053-6) (Linder et al., *Nature Genetics* 2025), a published sequence-to-function model, ourselves on an NVIDIA GPU via Brev (`brev/`). Per single-base variant:

1. 524,288 bases of hg38 centred on the variant; the reference base is checked.
2. Reference and alternate sequences through all 4 published replicates, forward and reverse complement (16 passes), averaged.
3. Target transforms undone; gene effect = log2((alt + 1) / (ref + 1)) of RNA-seq coverage over each gene's canonical GENCODE v44 exons, per tissue; local effect over the 512 bases around the variant.
4. Each GTEx eQTL is compared in direction with Borzoi's prediction for the same gene and allele. Only clear predictions (|log2| ≥ 0.005) are counted.

Results on the demo set: 42 alleles on 18 variants; prothrombin G20210A (rs1799963) is predicted to raise F2 expression by about 2.2%, matching its known mechanism; of 74 GTEx comparisons, 2 were clear and both agree; the other 72 had no clear predicted change and are reported, not counted. To rerun, see `brev/README.md`.

## Data checks and validation

Every time a variant loads, 11 checks run and the result appears under the summary:

1. gnomAD counts equal the raw response, and af = AC / AN.
2. Every percentage shown is computable from the raw counts.
3. Carrier rates follow from a frequency gnomAD reported.
4. Homozygotes ≤ half the copies observed.
5. dbSNP and gnomAD agree on position and alleles.
6. Scores are in range (REVEL, AlphaMissense, SIFT, PolyPhen-2, SpliceAI 0–1; CADD 0–99; LOEUF ≥ 0).
7. GWAS p-values valid (values below 1e-308 read from text), risk alleles on either strand, PMIDs numeric.
8. ClinPGx levels 1A–4, CIViC levels A–E, paper years not more than one year ahead.
9. Borzoi values plausible and for an allele dbSNP lists.
10. ACMG frequency codes never contradict a pathogenic ClinVar call silently.
11. Sources agree on the gene, or the disagreement is shown.

Validation: an independent cross-check of all 20 demo variants (`scripts/verify_cached.py`) found no wrong numbers, and frequencies match published values (sickle cell allele 5.7% in African ancestry, Factor V Leiden 2.6% and CFTR F508del 1.5% in Europeans). Live checks were run on uncached variants. 150 automated tests run in CI on every push.

## Engineering

- All sources in parallel; the page waits 8, 20 or 45 seconds (Settings → Speed). Late answers keep arriving in the background and the page redraws within 3 seconds; a source already being fetched is not fetched twice.
- Optional sources retry once; a failed source is remembered for 10 minutes and shown as Unavailable, never stuck on Loading.
- NCBI calls spaced 0.35 s (0.11 s with a key); gnomAD 1 s, with wait-and-retry on rate limits.
- Cache first: the demo variants are checked in; other variants are cached per user for 30 days.
- Local data: GWAS Catalog (23 MB, indexed in SQLite on first use) and the HPO inheritance table.
- The GPU work runs offline; the app never needs a GPU.
- Stack: Python 3.12, Streamlit, standard-library HTTP; PyTorch and borzoi-pytorch on the GPU side.

## Deploy (always-on demo)

The app runs on [Streamlit Community Cloud](https://streamlit.io/cloud), free, from this repository.

1. Open https://share.streamlit.io/deploy?repository=Sreesj/snp-brief&branch=main&mainModule=streamlit_app.py and sign in with the GitHub account that owns the repo.
2. Set **App URL** to `variant-atlas` (giving https://variant-atlas.streamlit.app) and Python to 3.12 under Advanced settings.
3. Optional, under **Secrets**: `NEBIUS_API_KEY` for AI summaries. Do not set `SNP_BRIEF_APP_PASSWORD` for a public demo.
4. Deploy. The first start takes a few minutes.

Community Cloud puts apps to sleep after a period without visitors. The `keep-awake` GitHub workflow opens the app every 6 hours and presses the wake button if it is asleep; set the repository variable `APP_URL` if your URL differs.

Other hosts: `docker build -t variant-atlas . && docker run -p 8501:8501 variant-atlas`.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/snp_brief/pipeline.py` | Parallel fetching, caching, time budget, late answers |
| `src/snp_brief/providers.py`, `sources/` | One module per source; `sources/gwas_local.py` for the GWAS file |
| `src/snp_brief/reconcile.py`, `brief.py` | ClinVar reconciliation and allele matching |
| `src/snp_brief/summary.py`, `llm.py` | Fact list, writers, fact check |
| `src/snp_brief/zygosity.py`, `interpret.py` | Zygosity check, carrier rates, ACMG draft |
| `src/snp_brief/integrity.py`, `methods.py` | Data checks and "How this is computed" notes |
| `src/snp_brief/ripple.py`, `forge.py`, `graph.py`, `literature.py` | Molecular Ripple, Hypothesis Forge, evidence graph, findings |
| `src/snp_brief/demo/` | The Streamlit app, design system and DNA animation |
| `brev/` | Borzoi GPU pipeline |
| `cache/`, `data/` | Cached demo variants; GWAS and HPO data files |
| `tests/`, `scripts/` | Automated tests; cross-check script |

## Limitations

- A research brief, not a diagnostic tool.
- The ACMG draft covers only criteria public data can support; segregation, de novo status, functional and case-control evidence need case data and an expert.
- Zygosity depends on HPO and ClinVar disease IDs; unmatched conditions show "inheritance not stated".
- Carrier rates assume Hardy-Weinberg equilibrium and are population estimates.
- Borzoi predictions are hypotheses about mechanism, single-base variants only.
- Public services can be slow or down; the app reports this rather than hiding it.
