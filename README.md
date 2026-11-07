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

The command prints one paragraph per rsID in `fixture/snps.txt`. Cached responses live in `cache/`. When a cache file is present the command does not call that source, so the pitch still runs if a source is down. `--offline` refuses every network call. `--refresh` replaces the cache from the live sources.

The same command runs any list:

```bash
uv run snp-brief --input fixture/held-out.txt --save briefs/held-out.txt
```

`fixture/held-out.txt` is one rsID that is not in the demo list. It is run once and saved. It is not used to tune the brief.

## Page

```bash
uv run snp-brief-demo
```

Serves a local page. Offline by default: it reads `cache/` and makes no network call. The sidebar starts filled with `fixture/snps.txt`. **Run briefs** builds one row per SNP, with function, frequency, and disease relevance marked sourced, unresolved, or blank. Select a row to open that paragraph, the gnomAD population counts when frequency is unresolved, the Europe PMC papers, and which sources answered. **Live fetching** calls the sources and writes a separate temp cache, leaving `cache/` unchanged.

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

When ClinVar submissions for one variant disagree, disease relevance lists each distinct interpretation with its review status and citation and marks that field unresolved.

gnomAD populations are part of that rule. Ancestry groups with fewer than 1000 alleles are not compared. Two frequencies disagree when the higher is at least twice the lower, or the gap is at least 0.1, or one is 0 and the other is at least 0.001. The brief then keeps both counts and names the populations. Exome and genome call sets use the same comparison on the allele frequency gnomAD reports. The comparison does not create a new number. Counts and frequencies in the paragraph are the ones the API returned.

## What this is not

Not a pathogenicity classifier, and not a copy of the genomic-variant-classifier repository. This demo writes a sourced brief.

## Team split

| Piece | Where |
| --- | --- |
| Ids and the fixture | `fixture/`, `src/snp_brief/fixture.py`, `src/snp_brief/sources/dbsnp.py` |
| ClinVar and function | `src/snp_brief/sources/clinvar.py`, `src/snp_brief/sources/uniprot.py` |
| gnomAD frequency | `src/snp_brief/sources/gnomad.py` |
| Europe PMC | `src/snp_brief/sources/europepmc.py` |
| Brief writer, conflict rules, and the demo | `src/snp_brief/brief.py`, `src/snp_brief/cli.py` |

The product rules are in `prd.md`.

## Example

`uv run snp-brief` prints one paragraph per id. This is the rs1815739 paragraph from the cached fixture. The stop-gained consequence is not treated as a disease claim, and the two ancestry counts are both kept because they disagree:

```text
rs1815739 is GRCh38 chr11:66560624, alleles C>A, C>T (dbSNP). Function: nonsense in ACTN3. NM_001104.4(ACTN3):c.1729C>T (p.Arg577Ter) (ClinVar VCV000018312, criteria provided, single submitter). Frequency for C>T: unresolved. gnomAD v4 exomes allele frequency 0.445748 (651404/1461372 alleles); gnomAD v4 exomes African/African American 5239/33480 alleles; gnomAD v4 exomes Admixed American 29090/44712 alleles. Disease relevance: unresolved. Benign (ClinVar SCV000040272, no assertion criteria provided), conditions ACTININ, ALPHA-3 POLYMORPHISM; Affects (ClinVar SCV000040273, no assertion criteria provided), conditions ACTN3 DEFICIENCY; Pathogenic (ClinVar SCV000040274, no assertion criteria provided), conditions SPRINTING PERFORMANCE; INCREASED COLD TOLERANCE; Benign (ClinVar SCV006106436, criteria provided, single submitter). Papers: Alfred T et al., 2011, ACTN3 genotype, athletic status, and life course physical capability: meta-analysis of the published literature and findings from nine studies, PMID:21542061; Eynon N et al., 2012, The ACTN3 R577X polymorphism across three groups of elite male European athletes, PMID:22916217; Tharabenjasin P et al., 2019, Association of the ACTN3 R577X (rs1815739) polymorphism with elite power sports: A meta-analysis, PMID:31145768 (Europe PMC).
```

## Tests

```bash
uv run python -m unittest discover -s tests -v
```

CI (`.github/workflows/ci.yml`) runs `uv sync --frozen`, these tests, and `uv run snp-brief --offline` on every push and pull request.
