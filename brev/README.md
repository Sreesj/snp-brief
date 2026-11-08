# Borzoi on NVIDIA Brev

Predicts each variant's effect on gene expression and chromatin with Borzoi, run on our GPUs,
and saves the result where the app reads it (Molecular Ripple → "Predicted regulatory effect").

## On your laptop
```
uv run snp-brief --export-variants brev/variants.tsv            # the 20 demo SNPs
# or any list:  uv run snp-brief --input my_ids.txt --export-variants brev/variants.tsv
git add brev && git commit -m "brev: variants to score" && git push
```

## On the Brev machine (terminal or a Jupyter terminal)
```
git clone https://github.com/Sreesj/snp-brief.git && cd snp-brief/brev
bash 01_check.sh        # GPU and internet; paste the output if anything says FAILED or 000
bash 02_setup.sh        # ~10 min: packages, hg38, Borzoi's 4 replicates, GENCODE
source .venv/bin/activate
python score_borzoi.py variants.tsv out/     # a few minutes for 20 SNPs on an H100
tar czf borzoi_results.tgz out/
```
Download `borzoi_results.tgz` (Brev's file browser, or `scp`), then on your laptop:
```
tar xzf borzoi_results.tgz -C brev/
uv run snp-brief --import-borzoi brev/out
git add cache/*/borzoi.json && git commit -m "data: Borzoi predictions" && git push
```
Stop the Brev instance when done.
