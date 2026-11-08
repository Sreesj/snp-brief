# PRD: SNP brief

Nucleate Pittsburgh BioHack 2026, QIAGEN BioTech track. Team of five. Demo by Sunday. One GitHub repo once the team can publish it. Use uv so the demo can be rerun.

This file is the product brief for the single repository. The modules below are the team split. They are not five projects.

## Problem

Given a list of SNPs, information on each one is scattered. Produce a clear plain-text summary of biological function, frequency, and disease relevance.

## User and output

A scientist hands in SNP ids. For each SNP the tool returns one paragraph a judge can read in a pitch: function, frequency, disease relevance. Every claim names its source. The tool does not diagnose and does not invent a number.

Input is one rsID per line. The pinned demo list is `fixture/snps.txt`. `fixture/held-out.txt` is one extra id, run once at the end and not used while tuning.

## Sources, in order

1. dbSNP (`src/snp_brief/sources/dbsnp.py`): turn the input into an rsID and a genomic coordinate.
2. ClinVar (`src/snp_brief/sources/clinvar.py`): gene, molecular consequence, clinical significance. Keep a record only when it cites this rsID. Do not treat a multi-variant haplotype as this SNP. Report different alternate alleles separately.
3. gnomAD (`src/snp_brief/sources/gnomad.py`): allele frequency. If populations disagree, keep both numbers and name the populations. In this demo, ancestry groups with fewer than 1000 alleles are not compared. A disagreement is a twofold difference, an absolute gap of at least 0.1, or a frequency of 0 against one of at least 0.001. Exome and genome call sets use that same comparison. The rule does not invent a frequency.
4. Europe PMC (`src/snp_brief/sources/europepmc.py`): a few primary papers for that variant, taken from research articles that name the rsID in the title or abstract.

UniProt (`src/snp_brief/sources/uniprot.py`) is only a function fallback. If ClinVar has no clinical significance, use the gene’s reviewed UniProt function text and cite UniProt. Leave disease relevance blank, or say no clinical significance was found. UniProt describes the protein, not whether this SNP causes disease. Never use UniProt for disease.

## Conflict layer

After the three fields are filled, compare sources on each field (`src/snp_brief/brief.py`). If they agree, keep the claim and cite them. If they conflict, keep both statements and mark the field unresolved. If nobody reports it, leave the field blank. Do not turn two alternate alleles into a false conflict.

## What we are not building

Not a pathogenicity classifier, and not the genomic-variant-classifier repo. That system predicts Pathogenic through Benign from millions of rows. This demo writes a sourced brief. Borrow only the public source list. Do not vendor or fork that repository.

## Demo

One pinned demo list of 20 SNPs (`fixture/snps.txt`), cached under `cache/`, plus one held-out SNP run once at the end. Tests use a frozen snapshot in `tests/cache/`. Cache fixture responses under `cache/` so the pitch still runs if a source is down. One command produces the briefs:

```bash
uv run snp-brief
```

## Team split

One person on ids and the fixture (`fixture/`, `fixture.py`, `sources/dbsnp.py`). One on ClinVar and function (`sources/clinvar.py`, `sources/uniprot.py`). One on gnomAD frequency (`sources/gnomad.py`). One on Europe PMC (`sources/europepmc.py`). One on the brief writer, the conflict rules, and the demo (`brief.py`, `cli.py`).
