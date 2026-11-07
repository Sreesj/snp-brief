"""Records passed from source parsers to the brief writer."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class DbSnpRecord:
    retrieved: bool
    rsid: str
    assembly: str | None = None
    accession: str | None = None
    chrom: str | None = None
    pos: int | None = None
    alleles: list[tuple[str, str]] = field(default_factory=list)
    genes: list[str] = field(default_factory=list)
    # Genomic alt base -> coding change such as c.20A>T, from the dbSNP HGVS list.
    coding_by_alt: dict[str, str] = field(default_factory=dict)


@dataclass
class Classification:
    description: str
    review_status: str
    traits: list[str]
    # Set for one submitted interpretation (an SCV). Empty on the VCV aggregate.
    accession: str = ""


@dataclass
class ClinVarRecord:
    accession: str
    title: str
    genes: list[str]
    consequences: list[str]
    cdna: str
    classifications: list[Classification]
    no_change: bool
    # Submitted interpretations (SCVs). Empty when the cache only has the aggregate.
    submissions: list[Classification] = field(default_factory=list)

    @property
    def has_significance(self) -> bool:
        return any(not empty_significance(item.description) for item in self.classifications)


def empty_significance(text: str) -> bool:
    folded = " ".join(text.lower().split())
    return folded in {
        "",
        "not provided",
        "no interpretation for the single variant",
        "no classifications from unflagged records",
        "no assertion provided",
        "no classification provided",
        "none",
        "no data",
        "no classification for the individual variant",
    }


def interpretation_label(text: str) -> str | None:
    """Folded clinical significance, or None when the text is not an interpretation."""
    folded = " ".join(text.lower().split())
    if empty_significance(folded):
        return None
    return folded


@dataclass
class PopCount:
    id: str
    ac: int
    an: int


@dataclass
class CallSet:
    """One gnomAD sequencing call set. Counts and af are copied from the API."""

    name: str
    ac: int | None
    an: int | None
    af: float | None
    populations: list[PopCount] = field(default_factory=list)


@dataclass
class GnomadVariant:
    variant_id: str
    reference_genome: str | None
    chrom: str
    pos: int
    ref: str
    alt: str
    callsets: list[CallSet]


@dataclass
class Paper:
    authors: str
    year: str
    title: str
    pmid: str


@dataclass
class UniProtHit:
    accession: str
    text: str


@dataclass
class AlleleView:
    ref: str
    alt: str
    records: list[ClinVarRecord]
    variant: GnomadVariant | None
