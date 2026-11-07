"""What the demo page shows for one SNP. Plain data, so it is tested without a browser.

Field statuses are read from the brief's own statements, so the table cannot
disagree with the paragraph. A field is unresolved when the brief marks it
unresolved, blank when it was not reported or ClinVar found no clinical
significance, and sourced otherwise.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import NamedTuple

from snp_brief.fixture import is_rsid
from snp_brief.models import CallSet
from snp_brief.pipeline import Gathered, SourceCall
from snp_brief.sources.gnomad import ANCESTRY, MIN_AN, frequency_phrase

FIELDS = ("Function", "Frequency", "Disease relevance")
SOURCED = "sourced"
UNRESOLVED = "unresolved"
BLANK = "blank"

_STATEMENT = re.compile(
    r"(?:Allele (?P<allele>[^:]+): )?(?P<field>Function|Frequency|Disease relevance)"
    r"(?: for (?P<label>[^:]+))?: (?P<claim>.*)",
    re.S,
)
_SEPARATORS = re.compile(r"[\s,;]+")


class Source(NamedTuple):
    key: str
    name: str
    role: str


SOURCES = (
    Source("dbsnp", "dbSNP", "rsID, coordinate, alleles"),
    Source("clinvar", "ClinVar", "gene, consequence, clinical significance"),
    Source("gnomad", "gnomAD", "allele frequency"),
    Source("europepmc", "Europe PMC", "papers"),
    Source("uniprot", "UniProt (fallback)", "function only, when ClinVar has no clinical significance; never disease"),
)
NOT_CONSULTED = "not consulted"


@dataclass
class Statement:
    field: str
    allele: str
    status: str
    text: str


@dataclass
class PopulationTable:
    """gnomAD ancestry counts for one alternate allele, as the API returned them."""

    allele: str
    variant_id: str
    totals: list[str]
    callsets: list[str]
    rows: list[dict[str, str]] = field(default_factory=list)
    named: set[tuple[str, str]] = field(default_factory=set)
    small: set[tuple[str, str]] = field(default_factory=set)


@dataclass
class PaperLink:
    title: str
    url: str
    citation: str


def parse_rsid_text(text: str) -> tuple[list[str], list[str]]:
    """Return the rsIDs in order without repeats, and the tokens that are not rsIDs. `#` starts a comment."""
    ids: list[str] = []
    rejected: list[str] = []
    for line in text.splitlines():
        for token in _SEPARATORS.split(line.split("#", 1)[0]):
            if not token:
                continue
            rsid = token.lower()
            if not is_rsid(rsid):
                rejected.append(token)
            elif rsid not in ids:
                ids.append(rsid)
    return ids, rejected


def statements(sentences: list[str]) -> list[Statement]:
    found: list[Statement] = []
    for sentence in sentences:
        match = _STATEMENT.fullmatch(sentence)
        if match is None:
            continue
        name = match["field"]
        allele = match["allele"] or match["label"] or ""
        found.append(Statement(name, allele, _status(name, match["claim"]), sentence))
    return found


def field_status(found: list[Statement], name: str) -> str:
    seen = {item.status for item in found if item.field == name}
    if UNRESOLVED in seen:
        return UNRESOLVED
    if SOURCED in seen:
        return SOURCED
    return BLANK


def gene_label(item: Gathered) -> str:
    genes: list[str] = []
    for record in item.records:
        for gene in record.genes:
            if gene not in genes:
                genes.append(gene)
    return ", ".join(genes or item.db.genes)


def summary_row(item: Gathered) -> dict[str, str]:
    found = statements(item.sentences)
    row = {"rsID": item.rsid, "Gene": gene_label(item)}
    for name in FIELDS:
        row[name] = field_status(found, name)
    return row


def unresolved_populations(item: Gathered) -> list[PopulationTable]:
    """One table per allele whose frequency the brief marks unresolved."""
    tables: list[PopulationTable] = []
    for variant in item.variants:
        unresolved, phrase = frequency_phrase(variant.callsets)
        if not unresolved:
            continue
        callsets = [call for call in variant.callsets if any(pop.id in ANCESTRY for pop in call.populations)]
        table = PopulationTable(
            allele=f"{variant.ref}>{variant.alt}",
            variant_id=variant.variant_id,
            totals=[_total(call) for call in variant.callsets],
            callsets=[call.name for call in callsets],
        )
        for pop_id, label in ANCESTRY.items():
            row = {"Population": label}
            for call in callsets:
                pop = next((pop for pop in call.populations if pop.id == pop_id), None)
                if pop is None:
                    row[call.name] = ""
                    continue
                row[call.name] = f"{pop.ac}/{pop.an}"
                if f"{call.name} {label} {pop.ac}/{pop.an} alleles" in phrase:
                    table.named.add((label, call.name))
                if pop.an < MIN_AN:
                    table.small.add((label, call.name))
            if any(row[call.name] for call in callsets):
                table.rows.append(row)
        tables.append(table)
    return tables


def paper_links(item: Gathered) -> list[PaperLink]:
    links: list[PaperLink] = []
    for paper in item.papers:
        first = paper.authors.split(",")[0].strip() or "Unknown authors"
        who = first + (" et al." if "," in paper.authors else "")
        year = f", {paper.year}" if paper.year else ""
        links.append(
            PaperLink(
                paper.title.rstrip("."),
                f"https://europepmc.org/article/MED/{paper.pmid}",
                f"{who}{year}, PMID:{paper.pmid}",
            )
        )
    return links


def source_rows(item: Gathered) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in SOURCES:
        call = item.calls.get(source.key)
        if call is None:
            answered, detail, gave = NOT_CONSULTED, _not_consulted(item), ""
        elif call.origin in {"cache", "live"}:
            answered, detail, gave = "yes", _response(call), _gave(source.key, item)
        else:
            answered, detail, gave = "no", _response(call), ""
        rows.append(
            {"Source": source.name, "Used for": source.role, "Answered": answered, "Detail": detail, "Gave": gave}
        )
    return rows


def briefs_text(items: list[Gathered]) -> str:
    """All briefs as one text, the same as `uv run snp-brief` prints them."""
    return "\n\n".join(item.text for item in items) + "\n"


def _status(name: str, claim: str) -> str:
    if claim.startswith("unresolved."):
        return UNRESOLVED
    if claim == "not reported.":
        return BLANK
    if name == "Disease relevance" and claim.startswith("no clinical significance was found"):
        return BLANK
    return SOURCED


def _total(call: CallSet) -> str:
    parts = [call.name]
    if call.af is not None:
        parts.append(f"allele frequency {format(call.af, '.6g')}")
    if call.ac is not None and call.an:
        parts.append(f"({call.ac}/{call.an} alleles)")
    return " ".join(parts)


def _response(call: SourceCall) -> str:
    when = _when(call.retrieved)
    if call.origin == "cache":
        return f"cached, retrieved {when}" if when else "cached"
    if call.origin == "live":
        return f"live, retrieved {when}" if when else "live"
    if call.origin == "no cache":
        return "no cached response, and live fetching is off"
    return f"request failed: {call.error}"


def _when(stamp: str | None) -> str:
    if not stamp:
        return ""
    try:
        moment = datetime.fromisoformat(stamp)
    except ValueError:
        return stamp
    if moment.tzinfo is not None:
        moment = moment.astimezone(timezone.utc)
    return moment.strftime("%Y-%m-%d %H:%M UTC")


def _not_consulted(item: Gathered) -> str:
    if not item.clinvar_retrieved:
        return "ClinVar was not retrieved"
    if item.has_significance:
        return "ClinVar reports clinical significance"
    return "no gene to look up"


def _gave(key: str, item: Gathered) -> str:
    db = item.db
    if key == "dbsnp":
        if not db.chrom or not db.pos:
            return "no coordinate for this rsID"
        place = f"{db.assembly or db.accession or ''} chr{db.chrom}:{db.pos}".strip()
        alleles = ", ".join(f"{ref}>{alt}" for ref, alt in db.alleles)
        return f"{place}, alleles {alleles}" if alleles else place
    if key == "clinvar":
        if not item.records:
            return "no record names this rsID"
        plural = "s" if len(item.records) > 1 else ""
        return f"{len(item.records)} record{plural}: " + ", ".join(record.accession for record in item.records)
    if key == "gnomad":
        if not item.variants:
            return "no variant for this rsID"
        return "; ".join(
            f"{variant.ref}>{variant.alt} in " + ", ".join(call.name for call in variant.callsets)
            for variant in item.variants
        )
    if key == "europepmc":
        if not item.papers:
            return "no research article names this rsID"
        plural = "s" if len(item.papers) > 1 else ""
        return f"{len(item.papers)} research article{plural}"
    if item.uniprot is not None:
        return f"{item.uniprot.accession}: reviewed {item.fallback_gene} function text, used for function only"
    return f"no reviewed function text for {item.fallback_gene}"
