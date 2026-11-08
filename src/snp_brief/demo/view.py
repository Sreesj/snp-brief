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
from snp_brief.sources.gnomad import ANCESTRY, MIN_AN, frequency_phrase, population_spread

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
    Source("clinvar", "ClinVar", "rank 1 · gene, consequence, clinical significance"),
    Source("gnomad", "gnomAD", "rank 2 · allele frequency"),
    Source("gwas", "GWAS Catalog", "rank 3 · trait associations (statistical, not clinical)"),
    Source("gtex", "GTEx", "rank 3 · expression QTLs (genotype linked to gene expression)"),
    Source("dbsnp", "dbSNP", "rank 4 · rsID, coordinate, alleles, location class"),
    Source("ensembl_ld", "Ensembl LD", "rank 4 · variants inherited together with this one"),
    Source("ensembl_vep", "Ensembl VEP", "rank 5 · consequences, regulatory overlap, prediction scores"),
    Source("open_targets", "Open Targets", "rank 5 · locus-to-gene score"),
    Source("uniprot", "UniProt (fallback)", "rank 6 · function only, when ClinVar has no clinical significance; never disease"),
    Source("ncbi_gene", "NCBI Gene", "rank 6 · gene description, when ClinVar has no clinical significance"),
    Source("europepmc", "Europe PMC", "rank 7 · papers naming the rsID"),
    Source("dbsnp_citations", "dbSNP citations", "rank 7 · PubMed records linked to the rsID"),
    Source("litvar", "LitVar2", "rank 7 · publications mentioning the rsID in full text"),
    Source("pharmgkb", "ClinPGx (PharmGKB)", "rank 2 · drug response, curated clinical annotations"),
    Source("constraint", "gnomAD constraint", "rank 2 · how much the gene tolerates damaging changes"),
    Source("myvariant", "MyVariant.info", "rank 5 · REVEL, AlphaMissense, CADD, SIFT, PolyPhen-2, CIViC cancer evidence"),
    Source("spliceai", "SpliceAI", "rank 5 · predicted splice-site gain or loss"),
    Source("regulome", "RegulomeDB", "rank 5 · regulatory evidence at this position"),
    Source("timeline", "Literature timeline", "rank 7 · papers per year naming the rsID"),
    Source("borzoi", "Borzoi (our model run)", "rank 5 · predicted effect on gene expression and chromatin, run on GPU by Variant Atlas"),
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
    label = f"{item.merged_from} → {item.rsid}" if getattr(item, "merged_from", "") else item.rsid
    row = {"rsID": label, "Gene": gene_label(item)}
    for name in FIELDS:
        row[name] = field_status(found, name)
    return row


def population_tables(item: Gathered) -> list[PopulationTable]:
    """One table per allele whose frequency varies by ancestry or is unresolved between call sets."""
    tables: list[PopulationTable] = []
    for variant in item.variants:
        unresolved, _phrase = frequency_phrase(variant.callsets)
        spread = population_spread(variant.callsets)
        if not unresolved and spread is None:
            continue
        named_pops: set[tuple[str, str]] = set()
        if spread is not None:
            pop_set, low, high = spread
            named_pops = {(ANCESTRY[low.id], pop_set.name), (ANCESTRY[high.id], pop_set.name)}
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
                if (label, call.name) in named_pops:
                    table.named.add((label, call.name))
                if pop.an < MIN_AN:
                    table.small.add((label, call.name))
            if any(row[call.name] for call in callsets):
                table.rows.append(row)
        tables.append(table)
    return tables


# Older name kept for callers written before ancestry spread stopped meaning "unresolved".
unresolved_populations = population_tables


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
        elif call.origin == "pending":
            answered, detail, gave = "loading", "still loading when the time budget ran out", ""
        elif call.origin == "failed":
            answered, detail, gave = "unavailable", _response(call), "The database did not respond; it is retried on the next refresh"
        else:
            answered, detail, gave = "no", _response(call), ""
        rows.append(
            {"Source": source.name, "Used for": source.role, "Answered": answered, "Detail": detail, "Gave": gave}
        )
    return rows


def briefs_text(items: list[Gathered], summary: str = "template") -> str:
    """All briefs as one text, the same as `uv run snp-brief --summary <mode>` prints them."""
    from snp_brief.cli import format_item

    return "\n\n".join(format_item(item, summary) for item in items) + "\n"


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
    if call.origin == "live" and call.seconds is not None:
        return f"live, {call.seconds:.1f} s, retrieved {_when(call.retrieved)}"
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
    if key == "gwas":
        hits = item.gwas.hits
        if not hits:
            return "no association for this rsID"
        traits = list(dict.fromkeys(hit.trait for hit in hits))
        return f"{len(hits)} association{'s' if len(hits) != 1 else ''}: " + ", ".join(traits[:3])
    if key == "ensembl_vep":
        vep = item.vep
        if not vep.retrieved or not (vep.most_severe or vep.genes):
            return "no annotation returned"
        return f"{vep.most_severe.replace('_', ' ')}" + (f" in {vep.canonical_gene}" if vep.canonical_gene else "")
    if key == "myvariant":
        sc = item.scores
        parts = []
        if sc.revel is not None:
            parts.append(f"REVEL {sc.revel:.2f}")
        if sc.alphamissense is not None:
            parts.append(f"AlphaMissense {sc.alphamissense:.2f}")
        if sc.cadd_phred is not None:
            parts.append(f"CADD {sc.cadd_phred:.1f}")
        if sc.civic:
            parts.append(f"{len(sc.civic)} CIViC evidence items")
        return ", ".join(parts) or (sc.consequence or "no scores for this variant")
    if key == "pharmgkb":
        n = len(item.pgx.annotations)
        if not n:
            return "no drug-response annotation"
        drugs = list(dict.fromkeys(d for a in item.pgx.annotations for d in a.drugs))
        return f"{n} annotation{'s' if n != 1 else ''}: " + ", ".join(drugs[:3])
    if key == "spliceai":
        sp = item.splice
        return f"max delta {sp.max_delta:.2f} ({sp.kind})" if sp.max_delta is not None else "no splice score"
    if key == "regulome":
        rg = item.regulome
        return f"rank {rg.ranking}" + (f", probability {rg.probability:.2f}" if rg.probability is not None else "") if rg.ranking else "no regulatory score"
    if key == "constraint":
        c = item.constraint
        return f"{c.gene} LOEUF {c.loeuf:.2f}" if c.loeuf is not None else "no constraint data"
    if key == "borzoi":
        bz = item.borzoi
        genes = [g for a in bz.alleles for g in a.genes]
        if not genes:
            return "no prediction"
        top = max(genes, key=lambda g: abs(g.mean_lfc))
        return f"strongest predicted change: {top.gene} (log2 {top.mean_lfc:+.3f})"
    if key == "timeline":
        t = item.timeline
        return f"{t.total} papers" + (f", {t.trend}" if t.trend else "") if t.total else "no papers"
    if key == "ncbi_gene":
        gs = item.gene_summary
        return f"{gs.symbol}: {gs.name}" if gs is not None else "no gene description"
    if key == "gtex":
        genes = list(dict.fromkeys(e.gene for e in item.eqtls))
        return f"{len(item.eqtls)} eQTLs: " + ", ".join(genes[:3]) if item.eqtls else "no significant eQTL"
    if key == "ensembl_ld":
        return f"{len(item.ld_proxies)} variants at r² ≥ 0.8" if item.ld_proxies else "no strong LD partners"
    if key == "open_targets":
        return ", ".join(f"{r.gene} {r.score:.2f}" for r in item.l2g[:3]) if item.l2g else "no locus-to-gene prediction"
    if key == "litvar":
        return f"{len(item.litvar)} publications" if item.litvar else "no publications found"
    if key == "dbsnp_citations":
        return f"{len(item.citations)} PubMed records" if item.citations else "no citations linked"
    if item.uniprot is not None:
        return f"{item.uniprot.accession}: reviewed {item.fallback_gene} function text, used for function only"
    return f"no reviewed function text for {item.fallback_gene}"
