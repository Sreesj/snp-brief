"""Zygosity check: for each alternate allele, does disease need one copy or two?

Shown first on the page, before anything else, so nobody designs an assay or primer and only
later learns the variant is pathogenic only in the homozygous state (or not at all for this
allele).

For each alternate allele of the rsID:
- the allele itself (ref>alt, protein change) and its own ClinVar classification and conditions,
- each condition's mode of inheritance from the Human Phenotype Ontology disease annotations
  (OMIM and Orphanet identifiers that ClinVar attaches to the condition), or from the condition
  name when it states it ("autosomal recessive ..."),
- what one copy and two copies mean under that mode,
- how many people in gnomAD carry one copy and two copies (homozygotes), with consistency flags:
  homozygotes observed for a severe recessive call, or many homozygotes for a dominant call,
  are worth checking.

Inheritance data: `snp-brief --build-inheritance phenotype.hpoa` writes data/inheritance.tsv.gz
from the HPO release file (https://hpo.jax.org/data/annotations).
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

INHERITANCE_FILE = "inheritance.tsv.gz"

# HPO inheritance terms (subtree of HP:0000005 "Mode of inheritance").
HPO_MODES = {
    "HP:0000006": "autosomal dominant",
    "HP:0000007": "autosomal recessive",
    "HP:0001417": "X-linked",
    "HP:0001419": "X-linked recessive",
    "HP:0001423": "X-linked dominant",
    "HP:0001427": "mitochondrial",
    "HP:0001450": "Y-linked",
    "HP:0001426": "multifactorial",
    "HP:0010982": "polygenic",
    "HP:0003745": "sporadic",
    "HP:0001442": "somatic mosaicism",
    "HP:0001428": "somatic mutation",
    "HP:0012275": "autosomal dominant somatic cell mutation",
    "HP:0001425": "heterogeneous",
    "HP:0032113": "semidominant",
    "HP:0001444": "autosomal dominant somatic cell mutation",
}

NAME_HINTS = (
    ("autosomal recessive", "autosomal recessive"),
    ("autosomal dominant", "autosomal dominant"),
    ("x-linked recessive", "X-linked recessive"),
    ("x-linked dominant", "X-linked dominant"),
    ("x-linked", "X-linked"),
    ("mitochondrial", "mitochondrial"),
)


def build(hpoa: Path, root: Path) -> tuple[Path, int]:
    """phenotype.hpoa -> data/inheritance.tsv.gz (disease id, disease name, mode)."""
    out = root / "data" / INHERITANCE_FILE
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with hpoa.open(encoding="utf-8") as handle, gzip.open(out, "wt", encoding="utf-8", newline="") as sink:
        writer = csv.writer(sink, delimiter="\t", lineterminator="\n")
        writer.writerow(("disease_id", "disease_name", "mode", "hpo_id"))
        header = None
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            row = dict(zip(header, fields))
            if row.get("aspect") != "I" or row.get("qualifier") == "NOT":
                continue
            mode = HPO_MODES.get(row.get("hpo_id", ""))
            if mode:
                writer.writerow((row["database_id"], row.get("disease_name", ""), mode, row["hpo_id"]))
                rows += 1
    return out, rows


def _data_file() -> Path | None:
    candidates = [Path(os.environ["SNP_BRIEF_ROOT"])] if os.environ.get("SNP_BRIEF_ROOT") else []
    candidates += [Path.cwd(), Path(__file__).resolve().parents[2]]
    for root in candidates:
        path = root / "data" / INHERITANCE_FILE
        if path.is_file():
            return path
    return None


@lru_cache(maxsize=1)
def _modes() -> dict[str, list[str]]:
    path = _data_file()
    table: dict[str, list[str]] = {}
    if path is None:
        return table
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            modes = table.setdefault(row["disease_id"], [])
            if row["mode"] not in modes:
                modes.append(row["mode"])
    return table


def inheritance_available() -> bool:
    return bool(_modes())


@dataclass
class Condition:
    name: str
    ids: list[str]
    modes: list[str]
    mode_source: str  # "HPO (OMIM:603903)", "condition name", or ""

    @property
    def one_copy(self) -> str:
        return _copy_meaning(self.modes, 1)

    @property
    def two_copies(self) -> str:
        return _copy_meaning(self.modes, 2)


@dataclass
class AlleleZygosity:
    allele: str  # "T>A"
    change: str  # "p.Glu7Val" when known
    classification: str  # ClinVar's call for this allele, or "" when no record
    review: str
    accession: str
    conditions: list[Condition] = field(default_factory=list)
    carriers: int | None = None  # allele count = copies seen
    homozygotes: int | None = None
    hemizygotes: int | None = None
    callset: str = ""
    flags: list[str] = field(default_factory=list)


def _copy_meaning(modes: list[str], copies: int) -> str:
    m = set(modes)
    if not m:
        return "not stated by the sources"
    if m & {"autosomal dominant", "X-linked dominant", "semidominant"} and "autosomal recessive" in m:
        return ("can cause disease (dominant form reported)" if copies == 1
                else "can cause disease (recessive and dominant forms both reported)")
    if m & {"autosomal dominant", "X-linked dominant"}:
        return "enough to cause disease (dominant)" if copies == 1 else "causes disease, often more severe"
    if "semidominant" in m:
        return "milder disease or risk (semidominant)" if copies == 1 else "more severe disease (semidominant)"
    if "autosomal recessive" in m:
        return "carrier, usually unaffected (recessive)" if copies == 1 else "causes disease (recessive)"
    if m & {"X-linked recessive", "X-linked"}:
        return ("males (one X): affected; females: usually carriers (X-linked)" if copies == 1
                else "females with two copies: affected (X-linked)")
    if "mitochondrial" in m:
        return "mitochondrial inheritance: copy number does not apply; depends on heteroplasmy"
    if m & {"multifactorial", "polygenic"}:
        return "contributes to risk with other genes and environment" if copies == 1 else "may add more risk"
    return "inheritance: " + ", ".join(sorted(m))


def _raw_clinvar(item) -> dict:
    directory = getattr(item, "cache_dir", None)
    if directory is None:
        return {}
    path = Path(directory) / item.rsid / "clinvar.json"
    try:
        return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    except (OSError, ValueError):
        return {}


def _trait_ids(raw: dict, accession: str) -> dict[str, list[str]]:
    """trait name -> ["OMIM:603903", "ORPHA:232", ...] from ClinVar's esummary."""
    result = ((raw.get("response") or {}).get("esummary") or {}).get("result") or {}
    out: dict[str, list[str]] = {}
    for uid in result.get("uids") or []:
        rec = result.get(uid) or {}
        # Every record's conditions are used: a condition's mode of inheritance does not depend on the allele.
        for block in ("germline_classification",):
            for trait in (rec.get(block) or {}).get("trait_set") or []:
                ids = []
                for x in trait.get("trait_xrefs") or []:
                    db, ident = x.get("db_source"), str(x.get("db_id", ""))
                    if db == "OMIM" and ident.isdigit():
                        ids.append(f"OMIM:{ident}")
                    elif db == "Orphanet" and ident.isdigit():
                        ids.append(f"ORPHA:{ident}")
                out[str(trait.get("trait_name", ""))] = ids
    return out


def _norm(name: str) -> str:
    return " ".join(re.findall(r"[a-z0-9]+", name.lower()))


def _closest(name: str, ids_by_name: dict[str, list[str]]) -> list[str]:
    """Submissions spell conditions differently from the record ('THROMBOPHILIA DUE TO FACTOR V LEIDEN' vs
    'Factor V Leiden thrombophilia'); match on shared words when the spelling differs."""
    words = set(_norm(name).split()) - {"due", "to", "of", "the", "and", "disease", "susceptibility", "type"}
    best, score = [], 0.0
    for key, ids in ids_by_name.items():
        other = set(key.split()) - {"due", "to", "of", "the", "and", "disease", "susceptibility", "type"}
        if not words or not other or not ids:
            continue
        overlap = len(words & other) / len(words | other)
        if overlap > score:
            best, score = ids, overlap
    return best if score >= 0.6 else []


def _modes_for(name: str, ids: list[str]) -> tuple[list[str], str]:
    table = _modes()
    for ident in ids:
        if ident in table:
            return table[ident], f"HPO ({ident})"
    lower = name.lower()
    for hint, mode in NAME_HINTS:
        if hint in lower:
            return [mode], "condition name"
    return [], ""


def check(item) -> list[AlleleZygosity]:
    from snp_brief.brief import _kept_traits, _views
    from snp_brief.reconcile import reconcile
    from snp_brief.summary import _change

    out: list[AlleleZygosity] = []
    raw = _raw_clinvar(item)
    views = _views(item.db, item.records, item.variants)
    for view in views:
        rec = next((r for r in view.records if r.has_significance), view.records[0] if view.records else None)
        classification = review = accession = change = ""
        conditions: list[Condition] = []
        names: list[str] = []
        if rec is not None:
            accession = rec.accession
            change = _change(rec)
            verdict = reconcile(view.records, _kept_traits)
            if verdict is not None and verdict.resolved and verdict.headline:
                classification, review, names = verdict.headline, verdict.tier_review, list(verdict.traits)
            elif verdict is not None and not verdict.resolved:
                classification = "Conflicting (" + "; ".join(f"{len(g)} {label}" for label, g in verdict.tallies) + ")"
                review = verdict.tier_review
                names = sorted({t for _, g in verdict.tallies for call in g for t in _kept_traits(call.classification.traits)})
            else:
                top = next((c for c in rec.classifications if c.description and c.description.lower() != "not provided"), None)
                if top is not None:
                    classification, review, names = top.description, top.review_status, _kept_traits(top.traits)
            ids_by_name = {_norm(k): v for k, v in _trait_ids(raw, rec.accession).items()}
            seen_names: set[str] = set()
            for name in names[:8]:
                key = _norm(name).rstrip("s")
                if key in seen_names:
                    continue  # "HBB-related disorder" and "HBB-related disorders"
                seen_names.add(key)
                ids = ids_by_name.get(_norm(name)) or _closest(name, ids_by_name)
                modes, source = _modes_for(name, ids)
                conditions.append(Condition(name, ids, modes, source))
        z = AlleleZygosity(f"{view.ref}>{view.alt}", change, classification, review, accession, conditions)
        if view.variant is not None and view.variant.callsets:
            cs = max(view.variant.callsets, key=lambda c: c.an or 0)
            z.carriers, z.homozygotes, z.hemizygotes, z.callset = cs.ac, cs.hom, cs.hemi, cs.name
        z.flags = _flags(z)
        out.append(z)
    return out


def _flags(z: AlleleZygosity) -> list[str]:
    flags = []
    pathogenic = "pathogenic" in z.classification.lower() and "conflict" not in z.classification.lower()
    modes = {m for c in z.conditions for m in c.modes}
    if pathogenic and not modes:
        flags.append("Inheritance of the condition is not stated by the sources: confirm whether one or two copies cause disease before designing an assay.")
    if pathogenic and z.homozygotes:
        if modes and modes <= {"autosomal recessive"}:
            flags.append(
                f"gnomAD includes {z.homozygotes:,} people with two copies. gnomAD excludes severe paediatric disease, "
                "so check penetrance or a milder phenotype."
            )
        elif modes & {"autosomal dominant"} and z.homozygotes >= 10:
            flags.append(
                f"{z.homozygotes:,} people in gnomAD have two copies of an allele called pathogenic for a dominant condition: "
                "suggests reduced penetrance, a risk allele, or a classification to re-check."
            )
    if z.homozygotes is None and z.carriers is not None:
        flags.append("Homozygote count not in the saved gnomAD data; Refresh data to fetch it.")
    return flags
