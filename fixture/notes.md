# Why these rsIDs

Pinned inputs for the SNP brief. `snps.txt` is the demo list (10). `held-out.txt` is one extra id, kept out of that list so it is run once at the end and not used while tuning. One id per line. Reasons below are why the id was picked. They are not fetched ClinVar, gnomAD, or paper claims, and they are not frequencies.

## Pathogenic

- **rs6025** — Factor V Leiden (`F5`). The standard pathogenic thrombophilia SNP.
- **rs334** — Sickle-cell allele (`HBB`). A second pathogenic SNP whose frequency is known to differ sharply by ancestry.
- **rs1800562** — HFE C282Y. Pathogenic hemochromatosis allele, common in some European ancestries and rare in others.
- **rs1799963** — Prothrombin G20210A (`F2`). Pathogenic thrombosis allele, less familiar than Factor V Leiden, so the brief is not only the one textbook SNP.

## Common benign

- **rs4988235** — Lactase persistence (`MCM6`). Common benign trait allele, high in Europe and low in East Asia, so population frequencies should disagree.
- **rs1815739** — ACTN3 R577X. Common stop-gained allele that is benign. Included so a severe-looking molecular consequence is not treated as disease.
- **rs1801133** — MTHFR C677T. Very common allele, generally benign, and a usual case of a variant that sounds clinical and is not.
- **rs1426654** — SLC24A5 pigmentation allele. Common benign, with a large ancestry frequency difference.
- **rs12913832** — HERC2 eye-color allele. Common benign variant in a non-disease trait.

## No clinical significance

- **rs6054257** — The ordinary example SNP from the VCF specification. Not a named pathogenic, benign, or risk allele. This is the slot where clinical significance should be missing.

## Held out

- **rs34637584** — LRRK2 G2019S. Well-known rare pathogenic Parkinson allele. Same shape as the pathogenic fixture ids, different gene and disease, and not in `snps.txt`.
