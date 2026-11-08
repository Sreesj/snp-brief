#!/usr/bin/env bash
# Step 2 on Brev: Python packages, the hg38 genome, Borzoi's track list and GENCODE genes. About 10 minutes.
set -euo pipefail
cd "$(dirname "$0")"
python3 -m venv .venv
source .venv/bin/activate
pip install -q --upgrade pip
pip install -q torch borzoi-pytorch pyfaidx pandas numpy
mkdir -p ref
[ -f ref/hg38.fa ] || { curl -L -o ref/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip ref/hg38.fa.gz; }
python3 -c "import pyfaidx; pyfaidx.Fasta('ref/hg38.fa'); print('genome indexed')"
[ -f ref/targets_human.txt ] || curl -L -o ref/targets_human.txt https://raw.githubusercontent.com/calico/borzoi/main/examples/targets_human.txt
[ -f ref/genes.gtf.gz ] || curl -L -o ref/genes.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz
python3 - <<'PY'
from borzoi_pytorch import Borzoi
for r in range(4):
    Borzoi.from_pretrained(f"johahi/borzoi-replicate-{r}")
    print(f"downloaded Borzoi replicate {r}")
PY
echo "setup done"
