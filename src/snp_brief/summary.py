"""Plain-English summary written only from retrieved facts, then checked.

The evidence brief (brief.py) is complete but dense. This module turns the same
retrieved data into three or four readable sentences.

1. `build_facts` lists every fact the summary may use. Each fact has an id
   (F1, F2, ...), its text, and the source that reported it. Percentages are
   computed from the counts the API returned; they are not new measurements.
2. A writer produces sentences that end with the fact ids they rely on, e.g.
   "... is a missense variant in F5 [F2]." The writer is either the built-in
   template (deterministic, offline) or an LLM (llm.py).
3. `check` rejects any sentence that cites no fact, cites an unknown fact,
   contains a number that is not in the facts it cites, names a gene or rsID
   that no fact contains, or reads like a diagnosis.
4. `render` swaps the ids for source names, so the reader sees "(ClinVar)".

If the LLM output fails the check, the template summary is used instead.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from snp_brief.brief import _kept_traits, _views
from snp_brief.models import AlleleView, ClinVarRecord, GnomadVariant
from snp_brief.reconcile import OTHER, reconcile, side
from snp_brief.sources.gnomad import ANCESTRY, _freq, frequency_phrase, percent, population_spread


@dataclass
class Fact:
    id: str
    text: str
    source: str


@dataclass
class Summary:
    rsid: str
    sentences: list[str]  # with [F#] tags
    facts: list[Fact]
    writer: str  # "template" or "llm:<model>"
    rejected: list[tuple[str, str]] = field(default_factory=list)  # (sentence, reason)

    @property
    def text(self) -> str:
        return render(self.sentences, self.facts)


# ---------------------------------------------------------------- facts


@dataclass
class _Plan:
    """Facts plus the ids the template needs, so it never guesses."""

    facts: list[Fact] = field(default_factory=list)
    function: str | None = None
    uniprot: str | None = None
    frequency: str | None = None
    callset_conflict: str | None = None
    spread: str | None = None
    disease: str | None = None
    disease_other: str | None = None
    dissent: str | None = None
    papers: list[str] = field(default_factory=list)
    # Plain values the template writes, each copied from the fact it cites.
    gene: str = ""
    consequence: str = ""
    change: str = ""
    alt: str = ""
    overall_pct: str = ""
    overall_set: str = ""
    high: tuple[str, str] | None = None
    low: tuple[str, str] | None = None
    conflict_pair: tuple[str, str] | None = None
    headline: str = ""
    conditions: list[str] = field(default_factory=list)
    unresolved_groups: list[str] = field(default_factory=list)
    other_labels: list[str] = field(default_factory=list)
    dissent_count: int = 0
    uniprot_text: str = ""
    no_significance: bool = False
    paper_labels: list[str] = field(default_factory=list)


def build_facts(item) -> list[Fact]:
    return _plan(item).facts


def _plan(item) -> _Plan:
    plan = _Plan()

    def add(text: str, source: str) -> str:
        fact_id = f"F{len(plan.facts) + 1}"
        plan.facts.append(Fact(fact_id, text, source))
        return fact_id

    db = item.db
    if db.retrieved and db.chrom and db.pos:
        alleles = ", ".join(f"{ref}>{alt}" for ref, alt in db.alleles)
        add(f"{item.rsid} is at {db.assembly or ''} chr{db.chrom}:{db.pos}, alleles {alleles}.".replace("  ", " "), "dbSNP")
    view = _main_view(item)
    if view is not None and view.records:
        record = view.records[0]
        plan.gene = ", ".join(record.genes) or ", ".join(db.genes)
        plan.consequence = ", ".join(record.consequences)
        plan.change = _change(record)
        plan.function = add(
            f"{item.rsid} is a {plan.consequence or 'variant with no stated consequence'} in {plan.gene or 'an unnamed gene'}"
            f"; {record.title or record.cdna} (ClinVar {record.accession}).",
            "ClinVar",
        )
    if item.uniprot and item.fallback_gene:
        plan.gene = plan.gene or item.fallback_gene
        plan.uniprot_text = _first_sentence(item.uniprot.text)
        plan.uniprot = add(
            f"Reviewed {item.fallback_gene} protein function: {item.uniprot.text} (UniProt {item.uniprot.accession}).",
            "UniProt",
        )
    variant = view.variant if view is not None else (item.variants[0] if item.variants else None)
    if variant is not None:
        _frequency_facts(plan, variant, add)
    if view is not None and view.records:
        _disease_facts(plan, view, add)
    elif item.clinvar_retrieved:
        plan.no_significance = True
        plan.disease = add(f"ClinVar reports no clinical significance for {item.rsid}.", "ClinVar")
    for paper in item.papers[:3]:
        first = paper.authors.split(",")[0].strip() or "Unknown"
        label = f"{first}{' et al.' if ',' in paper.authors else ''} {paper.year}".strip()
        plan.paper_labels.append(f"{label} (PMID:{paper.pmid})")
        plan.papers.append(add(f"Paper PMID:{paper.pmid}: {label}, {paper.title.rstrip('.')}.", "Europe PMC"))
    if not plan.no_significance and item.clinvar_retrieved and not item.has_significance:
        plan.no_significance = True
        plan.disease = plan.disease or add(f"ClinVar reports no clinical significance for {item.rsid}.", "ClinVar")
    return plan


def _main_view(item) -> AlleleView | None:
    views = _views(item.db, item.records, item.variants)
    for view in views:
        if view.records and any(record.has_significance for record in view.records):
            return view
    for view in views:
        if view.records:
            return view
    for view in views:
        if view.variant is not None:
            return view
    return None


def _change(record: ClinVarRecord) -> str:
    match = re.search(r"\((p\.[^)]+)\)", record.title or "")
    if match:
        return match.group(1)
    return record.cdna


def _first_sentence(text: str) -> str:
    head = re.split(r"(?<=[.;])\s", text.strip(), maxsplit=1)[0]
    return head.rstrip(".;")


def _frequency_facts(plan: _Plan, variant: GnomadVariant, add) -> None:
    plan.alt = f"{variant.ref}>{variant.alt}"
    usable = [call for call in variant.callsets if _freq(call) is not None]
    if not usable:
        return
    unresolved, _ = frequency_phrase(variant.callsets)
    spread = population_spread(variant.callsets)
    if unresolved and len(usable) >= 2:
        first, second = usable[0], usable[1]
        a = percent(first.ac or 0, first.an or 0) if first.an else f"{first.af}"
        b = percent(second.ac or 0, second.an or 0) if second.an else f"{second.af}"
        plan.conflict_pair = (f"{a} in {first.name}", f"{b} in {second.name}")
        plan.callset_conflict = add(
            f"{variant.ref}>{variant.alt} allele frequency disagrees between gnomAD call sets: "
            f"{first.name} {a} ({first.ac}/{first.an} alleles) and {second.name} {b} ({second.ac}/{second.an} alleles).",
            "gnomAD",
        )
    chosen = spread[0] if spread else max(usable, key=lambda call: call.an or 0)
    if chosen.an:
        plan.overall_pct = percent(chosen.ac or 0, chosen.an)
        plan.overall_set = chosen.name
        plan.frequency = add(
            f"{variant.ref}>{variant.alt} allele frequency in {chosen.name} is {plan.overall_pct} "
            f"({chosen.ac}/{chosen.an} alleles; reported af {format(chosen.af, '.6g') if chosen.af is not None else 'n/a'}).",
            "gnomAD",
        )
    if spread:
        pop_set, low, high = spread
        plan.high = (ANCESTRY[high.id], percent(high.ac, high.an))
        plan.low = (ANCESTRY[low.id], percent(low.ac, low.an))
        plan.spread = add(
            f"In {pop_set.name} the frequency varies by ancestry: highest in {plan.high[0]} {plan.high[1]} "
            f"({high.ac}/{high.an} alleles), lowest in {plan.low[0]} {plan.low[1]} ({low.ac}/{low.an} alleles).",
            "gnomAD",
        )


def _disease_facts(plan: _Plan, view: AlleleView, add) -> None:
    verdict = reconcile(view.records, _kept_traits)
    if verdict is None:
        # Aggregate only, e.g. "association" with LACTASE PERSISTENCE.
        record = view.records[0]
        meaningful = [c for c in record.classifications if c.description and c.description.lower() != "not provided"]
        if not meaningful:
            plan.no_significance = True
            plan.disease = add("ClinVar reports no clinical significance for this allele.", "ClinVar")
            return
        top = meaningful[0]
        plan.headline = top.description
        plan.conditions = _kept_traits(top.traits)[:3]
        plan.disease = add(
            f"ClinVar classification: {top.description}"
            + (f", conditions {'; '.join(plan.conditions)}" if plan.conditions else "")
            + f" (ClinVar {record.accession}, {top.review_status}).",
            "ClinVar",
        )
        return
    if not verdict.resolved:
        for label, group in verdict.tallies:
            plan.unresolved_groups.append(f"{len(group)} {label}")
        plan.disease = add(
            f"ClinVar submissions at the highest review level present ({verdict.tier_review}) disagree: "
            + "; ".join(plan.unresolved_groups)
            + ".",
            "ClinVar",
        )
    elif verdict.headline:
        plan.headline = verdict.headline
        plan.conditions = verdict.traits[:3]
        plan.disease = add(
            f"Best-reviewed ClinVar classification: {verdict.headline}"
            + (f", conditions {'; '.join(plan.conditions)}" if plan.conditions else "")
            + f" ({'; '.join(verdict.headline_cites)}"
            + (f"; {verdict.support}" if verdict.support else "")
            + ").",
            "ClinVar",
        )
    if verdict.other:
        plan.other_labels = [label for label, _ in verdict.other[:3]]
        bits = []
        for label, group in verdict.other[:3]:
            traits = _kept_traits(group[0].classification.traits)[:1]
            bits.append(label + (f" ({traits[0]})" if traits else ""))
        plan.disease_other = add("ClinVar also records these assertions: " + "; ".join(bits) + ".", "ClinVar")
    if verdict.dissent:
        plan.dissent_count = len(verdict.dissent)
        labels = sorted({call.label for call in verdict.dissent})
        plan.dissent = add(
            f"{plan.dissent_count} lower-review ClinVar submissions disagree ({', '.join(labels)}).",
            "ClinVar",
        )


# ---------------------------------------------------------------- template writer


def template_summary(item) -> Summary:
    plan = _plan(item)
    sentences: list[str] = []
    if plan.function:
        lead = f"{item.rsid} is {_phrase_consequence(plan.consequence)} in {plan.gene}"
        if plan.change:
            lead += f" ({plan.change})"
        sentences.append(f"{lead} [{plan.function}].")
    if plan.uniprot:
        sentences.append(
            f"ClinVar gives no clinical significance; UniProt describes the {plan.gene or 'encoded'} protein: "
            f"{plan.uniprot_text[:1].lower() + plan.uniprot_text[1:]} [{plan.uniprot}]."
        )
    if plan.callset_conflict and plan.conflict_pair:
        sentences.append(
            f"gnomAD's exome and genome call sets disagree on its frequency ({plan.conflict_pair[0]} vs {plan.conflict_pair[1]}) [{plan.callset_conflict}]."
        )
    if plan.frequency:
        text = f"The {plan.alt} allele is carried on {plan.overall_pct} of alleles in {plan.overall_set}"
        ids = [plan.frequency]
        if plan.spread and plan.high and plan.low:
            text += f", from {plan.high[1]} in {plan.high[0]} to {plan.low[1]} in {plan.low[0]}"
            ids.append(plan.spread)
        sentences.append(f"{text} [{', '.join(ids)}].")
    if plan.unresolved_groups and plan.disease:
        sentences.append(
            "Disease relevance is unresolved: ClinVar submissions at the same review level disagree ("
            + ", ".join(plan.unresolved_groups)
            + f") [{plan.disease}]."
        )
    elif plan.headline and plan.disease:
        if side(plan.headline) == OTHER:
            cond = f" with {', '.join(plan.conditions)}" if plan.conditions else ""
            sentences.append(f"ClinVar records it as {_article(plan.headline.lower())} {plan.headline.lower()}{cond}, not as a pathogenicity call [{plan.disease}].")
        else:
            cond = f" for {', '.join(plan.conditions)}" if plan.conditions else ""
            sentences.append(f"ClinVar's best-supported classification is {plan.headline}{cond} [{plan.disease}].")
    elif plan.no_significance and plan.disease:
        sentences.append(f"ClinVar reports no clinical significance for this variant [{plan.disease}].")
    extras = []
    extra_ids = []
    if plan.disease_other:
        extras.append("it also lists " + ", ".join(plan.other_labels) + " assertions")
        extra_ids.append(plan.disease_other)
    if plan.dissent:
        extras.append(f"{plan.dissent_count} lower-review submissions disagree")
        extra_ids.append(plan.dissent)
    if extras:
        sentences.append(f"ClinVar note: {'; '.join(extras)} [{', '.join(extra_ids)}].")
    if plan.papers:
        sentences.append(f"Primary literature: {'; '.join(plan.paper_labels)} [{', '.join(plan.papers)}].")
    return Summary(item.rsid, sentences, plan.facts, "template")


def _article(word: str) -> str:
    return "an" if word[:1].lower() in "aeiou" else "a"


def _phrase_consequence(consequence: str) -> str:
    if not consequence:
        return "a variant"
    first = consequence.split(",")[0].strip()
    if first == "nonsense":
        first = "nonsense (stop-gained)"
    if "variant" not in first:
        first += " variant"
    return f"{_article(first)} {first}"


# ---------------------------------------------------------------- checker

_TAGS = re.compile(r"\[(F\d+(?:\s*,\s*F\d+)*)\]")
_NUMBER = re.compile(r"(?<![A-Za-z])\d+(?:[.,]\d+)*(?:e-?\d+)?%?")
_NAMED = re.compile(r"\b(rs\d+|[A-Z][A-Z0-9]{2,}[0-9A-Z]*)\b")
_ALLOWED_CAPS = {"SNP", "DNA", "RNA", "GRCH38", "VCV", "SCV", "PMID", "UTR", "NM", "LLM", "OMIM", "USA", "HGVS"}
_FORBIDDEN = ("diagnos", "you have", "your risk", "patient should", "we recommend", "treatment")


def check(sentences: list[str], facts: list[Fact]) -> tuple[list[str], list[tuple[str, str]]]:
    """Return (kept sentences, rejected (sentence, reason) pairs)."""
    by_id = {fact.id: fact for fact in facts}
    corpus = " ".join(fact.text for fact in facts)
    kept: list[str] = []
    rejected: list[tuple[str, str]] = []
    for sentence in sentences:
        reason = _problem(sentence, by_id, corpus)
        if reason:
            rejected.append((sentence, reason))
        else:
            kept.append(sentence)
    return kept, rejected


def _problem(sentence: str, by_id: dict[str, Fact], corpus: str) -> str | None:
    tags = _TAGS.findall(sentence)
    if not tags:
        return "cites no fact"
    ids = [part.strip() for tag in tags for part in tag.split(",")]
    unknown = [fact_id for fact_id in ids if fact_id not in by_id]
    if unknown:
        return f"cites unknown fact {', '.join(unknown)}"
    cited = " ".join(by_id[fact_id].text for fact_id in ids)
    body = _TAGS.sub("", sentence)
    cited_numbers = {_norm_number(n) for n in _NUMBER.findall(cited)}
    for number in _NUMBER.findall(body):
        if _norm_number(number) not in cited_numbers:
            return f"number {number} is not in the cited facts"
    for name in _NAMED.findall(body):
        if name.upper() in _ALLOWED_CAPS:
            continue
        if name not in corpus:
            return f"{name} is not in any fact"
    lowered = body.lower()
    for word in _FORBIDDEN:
        if word in lowered:
            return f"reads like clinical advice ({word.strip()})"
    return None


def _norm_number(text: str) -> str:
    return text.replace(",", "").rstrip("%").rstrip(".")


def render(sentences: list[str], facts: list[Fact]) -> str:
    """Replace [F#] tags with the names of the sources that reported them."""
    by_id = {fact.id: fact for fact in facts}

    def swap(match: re.Match) -> str:
        names: list[str] = []
        for part in match.group(1).split(","):
            fact = by_id.get(part.strip())
            if fact and fact.source not in names:
                names.append(fact.source)
        return f"({', '.join(names)})" if names else ""

    out = []
    for sentence in sentences:
        text = _TAGS.sub(swap, sentence)
        out.append(re.sub(r"\s+([.,;])", r"\1", text).strip())
    return " ".join(out)


def facts_block(facts: list[Fact]) -> str:
    return "\n".join(f"[{fact.id}] ({fact.source}) {fact.text}" for fact in facts)


# ---------------------------------------------------------------- entry point


def summarize(item, mode: str = "template", client=None) -> Summary:
    """mode is "template" or "llm". LLM output that fails the check falls back to the template."""
    fallback = template_summary(item)
    kept, rejected = check(fallback.sentences, fallback.facts)
    fallback.sentences, fallback.rejected = kept, rejected
    if mode != "llm":
        return fallback
    from snp_brief.llm import LLMClient, write_summary

    client = client or LLMClient.from_env()
    if client is None:
        return fallback
    try:
        sentences = write_summary(client, item.rsid, fallback.facts)
    except Exception as exc:  # network, timeout, bad response
        fallback.rejected.append(("(LLM call)", f"failed: {exc}"))
        return fallback
    kept, rejected = check(sentences, fallback.facts)
    if not kept or len(rejected) > len(kept):
        fallback.rejected.extend(rejected)
        return fallback
    return Summary(item.rsid, kept, fallback.facts, f"llm:{client.model}", rejected)
