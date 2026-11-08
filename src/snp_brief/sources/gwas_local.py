"""The GWAS Catalog from its download file, so trait lookups never depend on a slow API.

Build once (the file comes from https://www.ebi.ac.uk/gwas/docs/file-downloads, "All
associations"; a .tsv or the .zip it ships in):

    uv run snp-brief --build-gwas ~/Downloads/gwas-catalog-associations.zip

That writes data/gwas_catalog.tsv.gz: only the columns the app uses, small enough to commit.
On first use the app turns it into an indexed SQLite file in the user cache, so each lookup
is a single indexed query.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import sqlite3
import threading
import zipfile
from datetime import datetime, timezone
from pathlib import Path

SLIM_NAME = "gwas_catalog.tsv.gz"
COLUMNS = ("rsid", "trait", "p_value", "risk_allele", "effect", "genes", "pmid", "accession", "sample")
_LOCK = threading.Lock()


def default_root() -> Path | None:
    for candidate in (Path(os.environ.get("SNP_BRIEF_ROOT", "")) if os.environ.get("SNP_BRIEF_ROOT") else None,
                      Path.cwd(), Path(__file__).resolve().parents[3]):
        if candidate is not None and slim_path(candidate).is_file():
            return candidate
    return None


def slim_path(root: Path) -> Path:
    return root / "data" / SLIM_NAME


def _open_source(path: Path):
    if path.suffix == ".zip":
        archive = zipfile.ZipFile(path)
        name = next(n for n in archive.namelist() if n.endswith(".tsv"))
        return io.TextIOWrapper(archive.open(name), encoding="utf-8", errors="replace")
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def build(source: Path, root: Path) -> tuple[Path, int]:
    """Slim the full download to the columns we use; single-SNP associations only."""
    out = slim_path(root)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    csv.field_size_limit(10_000_000)
    with _open_source(source) as handle, gzip.open(out, "wt", encoding="utf-8", newline="") as sink:
        reader = csv.DictReader(handle, delimiter="\t")
        writer = csv.writer(sink, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in reader:
            snps = (row.get("SNPS") or "").strip()
            if not snps.startswith("rs") or any(sep in snps for sep in (";", " x ", ",")):
                continue  # haplotypes and interactions are not this SNP alone
            odds, beta = (row.get("OR or BETA") or "").strip(), (row.get("95% CI (TEXT)") or "").strip()
            effect = ""
            if odds:
                effect = f"{'beta' if any(w in beta.lower() for w in ('increase', 'decrease', 'unit')) else 'odds ratio'} {odds}"
                if beta:
                    effect += f" {beta}"
            risk = (row.get("STRONGEST SNP-RISK ALLELE") or "").strip()
            if "-" in risk:
                risk = risk.split("-", 1)[1]
            writer.writerow(
                (
                    snps,
                    (row.get("MAPPED_TRAIT") or row.get("DISEASE/TRAIT") or "").strip(),
                    (row.get("P-VALUE") or "").strip(),
                    "" if risk == "?" else risk,
                    effect,
                    (row.get("REPORTED GENE(S)") or row.get("MAPPED_GENE") or "").strip(),
                    (row.get("PUBMEDID") or "").strip(),
                    (row.get("STUDY ACCESSION") or "").strip(),
                    (row.get("INITIAL SAMPLE SIZE") or "").strip()[:120],
                )
            )
            rows += 1
    return out, rows


def _db_path(slim: Path) -> Path:
    base = Path(os.environ.get("SNP_BRIEF_USER_CACHE") or Path.home() / ".cache" / "snp-brief")
    return base / f"gwas_catalog_{int(slim.stat().st_mtime)}.sqlite"


def _connect(root: Path) -> sqlite3.Connection | None:
    slim = slim_path(root)
    if not slim.is_file():
        return None
    db = _db_path(slim)
    with _LOCK:
        if not db.is_file():
            db.parent.mkdir(parents=True, exist_ok=True)
            tmp = db.with_suffix(".building")
            con = sqlite3.connect(tmp)
            con.execute(f"CREATE TABLE gwas ({', '.join(c + ' TEXT' for c in COLUMNS)})")
            with gzip.open(slim, "rt", encoding="utf-8") as handle:
                reader = csv.reader(handle, delimiter="\t")
                next(reader, None)
                con.executemany(f"INSERT INTO gwas VALUES ({', '.join('?' * len(COLUMNS))})", reader)
            con.execute("CREATE INDEX gwas_rsid ON gwas(rsid)")
            con.commit()
            con.close()
            tmp.replace(db)
    return sqlite3.connect(db, check_same_thread=False)


def available(root: Path) -> bool:
    return slim_path(root).is_file()


def lookup(rsid: str, root: Path) -> dict | None:
    """A cache payload in the same shape the Ensembl fallback uses, or None if no local file."""
    con = _connect(root)
    if con is None:
        return None
    try:
        rows = con.execute(f"SELECT {', '.join(COLUMNS)} FROM gwas WHERE rsid = ?", (rsid,)).fetchall()
    finally:
        con.close()
    phenotypes = []
    for row in rows:
        rec = dict(zip(COLUMNS, row))
        phenotypes.append(
            {
                "trait": rec["trait"],
                "source": "NHGRI-EBI GWAS catalog (download)",
                "pvalue": rec["p_value"],
                "risk_allele": rec["risk_allele"],
                "effect_text": rec["effect"],
                "genes": rec["genes"],
                "study": f"PMID:{rec['pmid']}" if rec["pmid"] else "",
                "accession": rec["accession"],
            }
        )
    stamp = datetime.fromtimestamp(slim_path(root).stat().st_mtime, timezone.utc).isoformat(timespec="seconds")
    return {
        "source": "GWAS Catalog",
        "retrieved": stamp,
        "request": {"rsid": rsid, "api": "catalog-download"},
        "response": {"ensembl_phenotypes": phenotypes},
    }
