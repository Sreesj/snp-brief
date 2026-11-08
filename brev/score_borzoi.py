"""Predict each variant's regulatory effect with Borzoi (Linder et al., Nature Genetics 2025).

Run on a GPU machine after 02_setup.sh:

    source .venv/bin/activate
    python score_borzoi.py variants.tsv out/

variants.tsv comes from `snp-brief --export-variants brev/variants.tsv` (rsid, chrom, pos, ref, alt;
GRCh38, 1-based). For every SNV it writes out/<rsid>/borzoi.json, which the app reads.

Method (also shown in the app under "How this is computed"):
1. Take 524,288 bases of hg38 centred on the variant; check the reference base matches the
   genome. Build the reference and the alternate sequence (one-hot).
2. Run Borzoi's four trained replicates on both sequences and on their reverse complements
   (reverse-complement outputs are flipped back and strand-paired), and average: 16 passes.
3. Undo Borzoi's target transforms (scale, soft clip, 3/4 power) to get coverage.
4. Gene expression: for each gene whose canonical exons fall in the 196,608-base output
   window, sum RNA-seq coverage over its exon bins (32 bases each) per track, and take
   log2((alt + 1) / (ref + 1)). Tracks with the same description (tissue or cell type) are averaged.
5. Local chromatin and transcription start: for DNase/ATAC, ChIP and CAGE tracks, the same
   log2 ratio over the 16 bins (512 bases) around the variant.
The output keeps the genes and tracks with the largest predicted changes.
"""

from __future__ import annotations

import gzip
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyfaidx
import torch

HERE = Path(__file__).resolve().parent
SEQ_LEN = 524_288
BIN = 32
OUT_BINS = 6_144
CROP = (SEQ_LEN - OUT_BINS * BIN) // 2  # 163,840 bases cropped on each side
LOCAL_BINS = 16
REPLICATES = 4
TOP_GENES = 8
TOP_TRACKS = 12


def one_hot(seq: str) -> np.ndarray:
    lut = np.zeros((256, 4), dtype=np.float32)
    for i, base in enumerate("ACGT"):
        lut[ord(base)] = lut[ord(base.lower())] = np.eye(4, dtype=np.float32)[i]
    return lut[np.frombuffer(seq.encode(), dtype=np.uint8)].T  # (4, L); N -> zeros


def load_targets() -> pd.DataFrame:
    t = pd.read_csv(HERE / "ref" / "targets_human.txt", sep="\t", index_col=0)
    t["assay"] = t["description"].str.split(":").str[0].str.upper()
    t["label"] = t["description"].str.split(":", n=1).str[1].fillna(t["description"])
    t["transform"] = np.where(t["sum_stat"].astype(str).str.contains("sqrt"), 0.75, 1.0)
    return t


def undo_transforms(y: np.ndarray, t: pd.DataFrame) -> np.ndarray:
    """y: (tracks, bins) model output -> coverage, following Borzoi's own scoring scripts."""
    y = y / t["scale"].to_numpy()[:, None]
    soft = t["clip_soft"].to_numpy()[:, None]
    over = y > soft
    y = np.where(over, (y - soft) ** 2 + soft, y)
    power = t["transform"].to_numpy()[:, None]
    return np.where(power != 1.0, np.clip(y, 0, None) ** (1.0 / power), y)


def load_genes(chrom_filter: set[str]) -> dict[str, list[tuple[str, str, list[tuple[int, int]]]]]:
    """chrom -> [(gene, strand, [(exon start, end), ...])] for canonical transcripts."""
    genes: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(HERE / "ref" / "genes.gtf.gz", "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] != "exon" or f[0] not in chrom_filter or "Ensembl_canonical" not in f[8]:
                continue
            gtype = re.search(r'gene_type "([^"]+)"', f[8])
            if not gtype or gtype.group(1) not in {"protein_coding", "lncRNA"}:
                continue
            name = re.search(r'gene_name "([^"]+)"', f[8]).group(1)
            genes[(f[0], name, f[6])].append((int(f[3]), int(f[4])))
    out: dict[str, list] = defaultdict(list)
    for (chrom, name, strand), exons in genes.items():
        out[chrom].append((name, strand, exons))
    return out


def models(device_ids: list[int]):
    from borzoi_pytorch import Borzoi

    loaded = []
    for r in range(REPLICATES):
        device = torch.device(f"cuda:{device_ids[r % len(device_ids)]}") if device_ids else torch.device("cpu")
        m = Borzoi.from_pretrained(f"johahi/borzoi-replicate-{r}").to(device).eval()
        loaded.append((m, device))
    return loaded


@torch.no_grad()
def predict(loaded, x: np.ndarray, strand_pair: np.ndarray) -> np.ndarray:
    """x: (2, 4, L) ref and alt -> (2, tracks, bins) averaged over replicates and strands."""
    fwd = torch.from_numpy(x)
    rc = torch.from_numpy(np.ascontiguousarray(x[:, ::-1, ::-1]))  # complement (ACGT reversed) and reverse
    total = None
    for m, device in loaded:
        with torch.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=device.type == "cuda"):
            a = m(fwd.to(device)).float().cpu().numpy()
            b = m(rc.to(device)).float().cpu().numpy()
        if a.shape[1:] != (len(strand_pair), OUT_BINS):
            raise SystemExit(f"unexpected Borzoi output shape {a.shape}; expected (2, {len(strand_pair)}, {OUT_BINS})")
        b = b[:, strand_pair, ::-1]  # back to forward orientation, strands swapped
        s = a + b
        total = s if total is None else total + s
    return total / (2 * len(loaded))


def main() -> int:
    variants_path, out_dir = Path(sys.argv[1]), Path(sys.argv[2])
    rows = pd.read_csv(variants_path, sep="\t", dtype=str)
    rows = rows[(rows["ref"].str.len() == 1) & (rows["alt"].str.len() == 1)]
    targets = load_targets()
    strand_pair = targets["strand_pair"].to_numpy().astype(int)
    assay = targets["assay"].to_numpy()
    label = targets["label"].to_numpy()
    rna = np.where(assay == "RNA")[0]
    local = np.where(np.isin(assay, ["DNASE", "ATAC", "CHIP", "CAGE"]))[0]
    fasta = pyfaidx.Fasta(str(HERE / "ref" / "hg38.fa"))
    chroms = {f"chr{c}" for c in rows["chrom"]}
    genes = load_genes(chroms)
    gpus = list(range(torch.cuda.device_count()))
    print(f"{len(rows)} SNVs, {len(targets)} tracks, GPUs: {gpus or 'none (CPU, slow)'}")
    loaded = models(gpus)
    results: dict[str, list] = defaultdict(list)
    for _, row in rows.iterrows():
        chrom, pos, ref, alt = f"chr{row['chrom']}", int(row["pos"]), row["ref"].upper(), row["alt"].upper()
        start = pos - 1 - SEQ_LEN // 2  # 0-based start; the variant sits at index SEQ_LEN // 2
        seq = str(fasta[chrom][max(0, start): start + SEQ_LEN]).upper()
        if start < 0:
            seq = "N" * (-start) + seq
        seq = seq.ljust(SEQ_LEN, "N")
        genome_ref = seq[SEQ_LEN // 2]
        entry = {"alt": alt, "ref": ref, "genome_ref": genome_ref, "ref_matches_genome": genome_ref == ref}
        if genome_ref != ref:
            print(f"{row['rsid']} {ref}>{alt}: genome has {genome_ref} at {chrom}:{pos}; skipped")
            results[row["rsid"]].append(entry)
            continue
        alt_seq = seq[: SEQ_LEN // 2] + alt + seq[SEQ_LEN // 2 + 1:]
        x = np.stack([one_hot(seq), one_hot(alt_seq)])
        y = predict(loaded, x, strand_pair)
        cov_ref, cov_alt = undo_transforms(y[0], targets), undo_transforms(y[1], targets)
        if not (np.isfinite(cov_ref).all() and np.isfinite(cov_alt).all()):
            entry["error"] = "non-finite prediction"
            results[row["rsid"]].append(entry)
            continue
        out_start = start + CROP  # 0-based genome position of output bin 0
        gene_rows = []
        for name, strand, exons in genes.get(chrom, []):
            bins = sorted({(p - 1 - out_start) // BIN for s, e in exons for p in range(s, e + 1, BIN)
                           if 0 <= (p - 1 - out_start) // BIN < OUT_BINS})
            if len(bins) < 3:
                continue
            r = cov_ref[rna][:, bins].sum(axis=1)
            a = cov_alt[rna][:, bins].sum(axis=1)
            lfc = np.log2(a + 1) - np.log2(r + 1)
            by_tissue: dict[str, list[float]] = defaultdict(list)
            for value, lab in zip(lfc, label[rna]):
                by_tissue[lab].append(float(value))
            tissues = {k: round(float(np.mean(v)), 5) for k, v in by_tissue.items()}
            mid = np.mean([(s + e) / 2 for s, e in exons])
            gene_rows.append({
                "gene": name, "strand": strand, "distance": int(abs(mid - pos)),
                "mean_lfc": round(float(lfc.mean()), 5),
                "max_abs_lfc": round(float(np.abs(lfc).max()), 5),
                "tissues": tissues,
            })
        gene_rows.sort(key=lambda g: -abs(g["mean_lfc"]))
        centre = OUT_BINS // 2
        window = slice(centre - LOCAL_BINS // 2, centre + LOCAL_BINS // 2)
        r = cov_ref[local][:, window].sum(axis=1)
        a = cov_alt[local][:, window].sum(axis=1)
        lfc = np.log2(a + 1) - np.log2(r + 1)
        order = np.argsort(-np.abs(lfc))[:TOP_TRACKS]
        entry["chromatin"] = [
            {"assay": assay[local][i], "label": label[local][i], "lfc": round(float(lfc[i]), 5)} for i in order
        ]
        entry["assay_summary"] = {
            k: round(float(np.mean(np.abs(lfc[assay[local] == k]))), 5) for k in ("DNASE", "ATAC", "CHIP", "CAGE")
            if (assay[local] == k).any()
        }
        entry["genes"] = gene_rows[:TOP_GENES]
        results[row["rsid"]].append(entry)
        top = gene_rows[0] if gene_rows else None
        print(f"{row['rsid']} {ref}>{alt}: top gene {top['gene'] if top else '-'} "
              f"(mean log2FC {top['mean_lfc'] if top else 0:+.4f}); top local track {entry['chromatin'][0]['label'] if entry['chromatin'] else '-'}")
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    for rsid, alleles in results.items():
        folder = out_dir / rsid
        folder.mkdir(parents=True, exist_ok=True)
        payload = {
            "source": "Borzoi",
            "retrieved": stamp,
            "request": {"rsid": rsid, "model": "johahi/borzoi-replicate-0..3", "assembly": "GRCh38"},
            "response": {"model": f"Borzoi, {REPLICATES} replicates, reverse-complement averaged", "alleles": alleles},
        }
        (folder / "borzoi.json").write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"wrote {len(results)} files under {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
