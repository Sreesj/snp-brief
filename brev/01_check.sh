#!/usr/bin/env bash
# Step 1 on Brev: is there a GPU, and can this machine download what we need?
echo "== GPU =="; nvidia-smi --query-gpu=name,memory.total --format=csv || echo "NO GPU FOUND"
echo "== Python =="; python3 --version
for url in https://huggingface.co https://hgdownload.soe.ucsc.edu https://raw.githubusercontent.com https://ftp.ebi.ac.uk https://pypi.org; do
  code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 15 "$url")
  echo "$url -> HTTP $code"
done
echo "== Disk =="; df -h . | tail -1
