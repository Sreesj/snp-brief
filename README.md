# SNP brief

A scientist hands in SNP ids. For each id the tool writes one plain-text paragraph a judge can read in a pitch: biological function, allele frequency, and disease relevance. Every claim names its source. The tool does not diagnose and does not invent a number.

This is one repository. Fixture parsing, ClinVar, gnomAD, Europe PMC, and the brief writer are modules, not separate projects.

## Demo

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then:

```bash
git clone https://github.com/Sreesj/snp-brief.git
cd snp-brief
uv run snp-brief
```

On the first run uv creates `.venv/` and installs the package from `uv.lock`. It uses Python 3.12 from `.python-version` and downloads it if it is not installed.

The command prints, for each rsID in `fixture/snps.txt`, a short plain-English summary and the full evidence brief. Cached responses live in `cache/`. When a cache file is present the command does not call that source, so the pitch still runs if a source is down. `--offline` refuses every network call. `--refresh` replaces the cache from the live sources.

The same command runs any list:

```bash
uv run snp-brief --input fixture/held-out.txt --save briefs/held-out.txt
```

`fixture/held-out.txt` is one rsID that is not in the demo list. It is run once and saved. It is not used to tune the brief.

## Page: Variant Atlas

```bash
uv run snp-brief-demo
```

Serves a local page. Run the list in the sidebar, or search any rsID at the top. Loading is cache first: SNPs in the checked-in `cache/` load instantly with no network call; any other rsID is fetched live from dbSNP, ClinVar, gnomAD, Europe PMC, and UniProt (10-30 s the first time) and saved to `~/.cache/snp-brief` (or `SNP_BRIEF_USER_CACHE`), so the next lookup is instant. **Refresh evidence** on a SNP re-fetches it live and keeps the old record if any source fails. **Offline only** in the sidebar turns off all network calls. The selected SNP opens in three modules:

- **Variant Observatory** (default): the checked plain-English summary, the evidence brief, gnomAD populations, papers, and which sources answered.
- **Molecular Ripple** (runs when you press Run): the variant context (missense, nonsense, intronic, ...), what a mechanistic analysis should examine for that context, and molecular links labelled *observed*, *predicted*, or *unknown*. Ensembl VEP scores are the next step.
- **Hypothesis Forge** (runs when you press Run): stage 1 sorts the variant as *well characterized*, *conflicting evidence*, *uncertain mechanism*, or *insufficient evidence*, with the reason and next steps. For conflicting and uncertain variants, an LLM proposes two or three competing hypotheses, each with fact-checked support, missing evidence, the observation that would tell it apart, and an experiment. Experiments are ranked by feasibility and open gaps. Well-characterized and thin-evidence variants get no invented hypotheses.

Module results stay in the session when you switch tabs or SNPs. The checked-in `cache/` is never changed by the page.

## Deploy (Streamlit Community Cloud)

1. Push to GitHub. The repo has `streamlit_app.py` (entry point) and `requirements.txt`.
2. At share.streamlit.io, sign in with GitHub, choose **New app**, repo `Sreesj/snp-brief`, branch `main`, main file `streamlit_app.py`.
3. Under **Advanced settings → Secrets**, paste the contents of `.streamlit/secrets.toml.example` with your values: the Nebius key, the model, an app password, and the per-session limit on AI hypothesis runs.
4. Deploy, then open the link and search a few rsIDs.

The password keeps strangers from spending the Nebius credits. Live responses are saved in the app's own cache; if the disk is read-only, results are shown but not saved.

## LLM: Nebius Token Factory

```bash
export NEBIUS_API_KEY=...                 # from tokenfactory.nebius.com; never commit it
uv run snp-brief --list-models            # models your account can use
export SNP_BRIEF_LLM_MODEL=<model id>     # optional; default openai/gpt-oss-120b
export SNP_BRIEF_FORGE_MODEL=<model id>   # optional; a stronger model for Hypothesis Forge
uv run snp-brief --summary llm --forge    # LLM summaries and Forge hypotheses
```

LLM summaries are saved to `cache/<rsid>/summary_llm.json` and Forge results to `cache/<rsid>/forge.json`, keyed by the facts they were written from. Generate them once with the key, commit them, and the demo shows them offline with no key. `--regenerate` ignores those files. Any other OpenAI-compatible server works with `SNP_BRIEF_LLM_URL`, `SNP_BRIEF_LLM_MODEL`, and `SNP_BRIEF_LLM_KEY`.

## Sources, in order

1. **dbSNP** (`src/snp_brief/sources/dbsnp.py`) turns the id into an rsID, a genomic coordinate, and the alternate alleles.
2. **ClinVar** (`src/snp_brief/sources/clinvar.py`) supplies the gene, molecular consequence, and clinical significance. A record is kept only when it cross-references this rsID. A multi-variant haplotype is not this SNP, so it is left out. Distinct alternate alleles are reported separately. They are not treated as a disagreement about one allele.
3. **gnomAD** (`src/snp_brief/sources/gnomad.py`) supplies allele counts and the reported allele frequency, for each alternate allele the rsID actually has. v4 is used when the variant is there, otherwise v3.
4. **Europe PMC** (`src/snp_brief/sources/europepmc.py`) supplies up to three primary research articles that name the rsID in the title or abstract.

**UniProt** (`src/snp_brief/sources/uniprot.py`) is only a function fallback. If ClinVar has no clinical significance, the brief uses the gene’s reviewed function text and cites UniProt. That text is protein function. It is never copied into disease relevance. Disease relevance then says that no clinical significance was found, or the field is left blank when ClinVar was not retrieved.

## Conflict rule

After the fields are filled, sources are compared on each field.

- If they agree, the brief keeps the claim and cites them.
- If they conflict, the brief keeps both statements and marks the field unresolved.
- If nobody reports it, the field says it was not reported. A missing frequency is not written as zero.

**Disease relevance is reconciled by ClinVar review status** (`src/snp_brief/reconcile.py`):

1. Flagged submissions are excluded and counted.
2. Pathogenicity calls (pathogenic side, benign side, uncertain) are ranked by ClinVar review stars. The headline is what the highest review tier present agrees on. An aggregate classification reviewed by an expert panel or practice guideline is used directly.
3. Lower-review submissions on another side are listed as named dissent, not as a conflict.
4. Other assertion types (risk factor, drug response, protective, association, affects) are reported as additional ClinVar assertions, never as a conflict.
5. The field is unresolved only when the top review tier itself disagrees (MTHFR rs1801133 in the fixture).

So Factor V Leiden (rs6025) reads Pathogenic with its drug-response and risk-factor assertions alongside, and ACTN3 (rs1815739) reads Benign with two unreviewed "Pathogenic" submissions noted as dissent.

**Frequency.** A spread between ancestry groups is the finding, not a disagreement between sources: the brief reports the overall frequency and the highest and lowest ancestry groups. Ancestry groups with fewer than 1000 alleles are not compared. Two frequencies differ when the higher is at least twice the lower, or the gap is at least 0.1, or one is 0 and the other is at least 0.001. When gnomAD's exome and genome call sets differ by that rule, both are reported with a note that they sample different populations (exomes are mostly European ancestry); this is not treated as a conflict. Counts are the ones the API returned; percentages are computed from those counts.

## Plain-English summary

Above each evidence brief the tool prints a three-to-five sentence summary (`src/snp_brief/summary.py`):

1. Every fact the summary may use is listed with an id and its source (dbSNP, ClinVar, gnomAD, UniProt, Europe PMC).
2. A writer produces sentences that end with the fact ids they use.
3. A checker rejects any sentence that cites no fact, uses a number that is not in the facts it cites, names a gene or rsID no fact contains, or reads like clinical advice.
4. Ids are replaced with source names, so each sentence ends with "(ClinVar)", "(gnomAD)", and so on.

The default writer is a deterministic template, so the demo works offline. To use an LLM, point the tool at any OpenAI-compatible server (vLLM or NVIDIA NIM on a Brev or Nebius GPU, the NVIDIA API catalog, or Ollama):

```bash
export SNP_BRIEF_LLM_URL=http://<gpu-host>:8000/v1
export SNP_BRIEF_LLM_MODEL=<model name>
export SNP_BRIEF_LLM_KEY=<token, if needed>
uv run snp-brief --summary llm
```

LLM sentences go through the same checker. If more sentences fail than pass, or the call fails, the template summary is printed instead. `--summary off` prints only the evidence brief. The demo page shows the summary, the facts it was written from, and any rejected sentences.

## What this is not

Not a pathogenicity classifier, and not a copy of the genomic-variant-classifier repository. This demo writes a sourced brief.

## Team split

| Piece | Where |
| --- | --- |
| Ids and the fixture | `fixture/`, `src/snp_brief/fixture.py`, `src/snp_brief/sources/dbsnp.py` |
| ClinVar and function | `src/snp_brief/sources/clinvar.py`, `src/snp_brief/sources/uniprot.py` |
| gnomAD frequency | `src/snp_brief/sources/gnomad.py` |
| Europe PMC | `src/snp_brief/sources/europepmc.py` |
| Brief writer, conflict rules, and the demo | `src/snp_brief/brief.py`, `src/snp_brief/reconcile.py`, `src/snp_brief/cli.py` |
| Plain-English summary and fact checker | `src/snp_brief/summary.py`, `src/snp_brief/llm.py` |
| Hypothesis Forge | `src/snp_brief/forge.py` |
| Molecular Ripple | `src/snp_brief/ripple.py` |
| Module result cache | `src/snp_brief/module_cache.py` |

The product rules are in `prd.md`.

## Example

`uv run snp-brief --input fixture/snps.txt` prints, for rs1815739 (ACTN3 R577X), from the cached fixture:

```text
rs1815739
Summary (template): rs1815739 is a nonsense (stop-gained) variant in ACTN3 (p.Arg577Ter) (ClinVar). The C>T allele is carried on 45% of alleles in gnomAD v4 exomes, from 65% in Admixed American to 16% in African/African American (gnomAD). ClinVar's best-supported classification is Benign for ACTININ, ALPHA-3 POLYMORPHISM (ClinVar). ClinVar note: it also lists Affects assertions; 2 lower-review submissions disagree (ClinVar). Primary literature: Alfred T et al. 2011 (PMID:21542061); Eynon N et al. 2012 (PMID:22916217); Tharabenjasin P et al. 2019 (PMID:31145768) (Europe PMC).
Evidence: rs1815739 is GRCh38 chr11:66560624, alleles C>A, C>T (dbSNP). Function: nonsense in ACTN3. ...
```

The stop-gained consequence is not treated as a disease claim, and the unreviewed "Pathogenic" submissions are kept as dissent rather than overturning the reviewed Benign call.

## Tests

```bash
uv run python -m unittest discover -s tests -v
```

CI (`.github/workflows/ci.yml`) runs `uv sync --frozen`, these tests, and `uv run snp-brief --offline` on every push and pull request.
