"""Write one plain-text paragraph. Agree and cite, keep both and mark unresolved, or leave it blank.

Different alternate alleles are not a conflict. UniProt text is protein function
and is never written into disease relevance.
"""

from __future__ import annotations

from snp_brief.models import (
    AlleleView,
    Classification,
    ClinVarRecord,
    DbSnpRecord,
    GnomadVariant,
    Paper,
    UniProtHit,
    empty_significance,
    interpretation_label,
)
from snp_brief.sources.gnomad import frequency_phrase


def assemble_brief(
    db: DbSnpRecord,
    records: list[ClinVarRecord],
    variants: list[GnomadVariant],
    papers: list[Paper],
    uniprot: UniProtHit | None,
    *,
    clinvar_retrieved: bool,
) -> str:
    return " ".join(
        brief_sentences(db, records, variants, papers, uniprot, clinvar_retrieved=clinvar_retrieved)
    )


def brief_sentences(
    db: DbSnpRecord,
    records: list[ClinVarRecord],
    variants: list[GnomadVariant],
    papers: list[Paper],
    uniprot: UniProtHit | None,
    *,
    clinvar_retrieved: bool,
) -> list[str]:
    """The paragraph in order. A field statement starts with its field name, after any allele prefix."""
    sentences = [_identity(db, variants)]
    note = _no_change_note(records)
    if note:
        sentences.append(note)
    used = _used_records(records)
    has_significance = any(record.has_significance for record in used)
    views = _views(db, records, variants)
    if not has_significance:
        gene = _fallback_gene(records, db)
        wrote_function = False
        if any(view.records for view in views):
            multi = len([view for view in views if view.records]) > 1
            for view in views:
                if not view.records:
                    continue
                sentences.append(_function(view, db, _prefix(view, multi)))
                wrote_function = True
        if uniprot and gene:
            sentences.append(
                f"Function: reviewed {gene} protein function: {uniprot.text} (UniProt {uniprot.accession})."
            )
            wrote_function = True
        if not wrote_function:
            sentences.append("Function: not reported.")
        if views:
            for view in views:
                sentences.append(_frequency(view, prefix=""))
        else:
            sentences.append("Frequency: not reported.")
        if clinvar_retrieved:
            sentences.append("Disease relevance: no clinical significance was found (ClinVar).")
        else:
            sentences.append("Disease relevance: not reported.")
    else:
        if not views:
            views = [AlleleView("", "", used, variants[0] if len(variants) == 1 else None)]
        multi = len(views) > 1
        for view in views:
            prefix = _prefix(view, multi)
            sentences.append(_function(view, db, prefix))
            sentences.append(_frequency(view, prefix))
            sentences.append(_disease(view, clinvar_retrieved=clinvar_retrieved, prefix=prefix))
    paper_line = _papers(papers)
    if paper_line:
        sentences.append(paper_line)
    return [sentence.strip() for sentence in sentences if sentence.strip()]


def _used_records(records: list[ClinVarRecord]) -> list[ClinVarRecord]:
    changing = [record for record in records if not record.no_change]
    return changing or list(records)


def _views(
    db: DbSnpRecord,
    records: list[ClinVarRecord],
    variants: list[GnomadVariant],
) -> list[AlleleView]:
    pool = _used_records(records)
    pairs = list(db.alleles)
    if not pairs and variants:
        pairs = [(variant.ref, variant.alt) for variant in variants]
    views: list[AlleleView] = []
    consumed: set[int] = set()
    attached: set[int] = set()
    for ref, alt in pairs:
        cdna = db.coding_by_alt.get(alt, "")
        matched: list[ClinVarRecord] = []
        for record in pool:
            if id(record) in consumed:
                continue
            if cdna and record.cdna.lower() == cdna.lower():
                matched.append(record)
                consumed.add(id(record))
        if not matched and len(pool) == 1 and len(pairs) == 1 and id(pool[0]) not in consumed:
            matched = [pool[0]]
            consumed.add(id(pool[0]))
        variant = _find_variant(variants, ref, alt)
        if variant is not None:
            attached.add(id(variant))
        if matched or variant is not None:
            views.append(AlleleView(ref, alt, matched, variant))
    unmatched = [record for record in pool if id(record) not in consumed]
    # One ClinVar record and one gnomAD variant are the same allele when the
    # coding-change map could not link them. Do not file them as two alleles.
    if len(unmatched) == 1 and len(variants) == 1:
        target = variants[0]
        placed = False
        for index, view in enumerate(views):
            if view.variant is target:
                views[index] = AlleleView(view.ref or target.ref, view.alt or target.alt, unmatched, target)
                placed = True
                break
        if not placed:
            views.append(AlleleView(target.ref, target.alt, unmatched, target))
            attached.add(id(target))
        unmatched = []
    for record in unmatched:
        views.append(AlleleView("", "", [record], None))
    for variant in variants:
        if id(variant) not in attached:
            views.append(AlleleView(variant.ref, variant.alt, [], variant))
    return views


def _find_variant(variants: list[GnomadVariant], ref: str, alt: str) -> GnomadVariant | None:
    for variant in variants:
        if variant.ref == ref and variant.alt == alt:
            return variant
    return None


def _prefix(view: AlleleView, multi: bool) -> str:
    if not multi:
        return ""
    if view.ref and view.alt:
        return f"Allele {view.ref}>{view.alt}: "
    if view.records and view.records[0].cdna:
        return f"Allele {view.records[0].cdna}: "
    if view.records:
        return f"Allele {view.records[0].accession}: "
    return ""


def _identity(db: DbSnpRecord, variants: list[GnomadVariant]) -> str:
    if not db.retrieved:
        return f"{db.rsid} (dbSNP record not retrieved)."
    coord = _coord(db)
    alleles = ""
    if db.alleles:
        alleles = ", alleles " + ", ".join(f"{ref}>{alt}" for ref, alt in db.alleles)
    if coord:
        base = f"{db.rsid} is {coord}{alleles} (dbSNP)"
    else:
        base = f"{db.rsid} (dbSNP genomic coordinate not reported)"
    conflicts: list[str] = []
    for variant in variants:
        if db.chrom and db.pos and variant.chrom and variant.pos:
            if _norm(db.chrom) != _norm(variant.chrom) or db.pos != variant.pos:
                conflicts.append(
                    f"gnomAD {variant.reference_genome or ''} chr{variant.chrom}:{variant.pos} {variant.ref}>{variant.alt}".strip()
                )
        if db.assembly and variant.reference_genome and db.assembly != variant.reference_genome:
            conflicts.append(f"gnomAD assembly {variant.reference_genome}")
    unique = list(dict.fromkeys(conflicts))
    if unique:
        return f"{base}. Coordinate: unresolved. {'; '.join(unique)}."
    return base + "."


def _coord(db: DbSnpRecord) -> str | None:
    if not db.chrom or not db.pos:
        return None
    if db.assembly:
        return f"{db.assembly} chr{db.chrom}:{db.pos}"
    if db.accession:
        return f"{db.accession} chr{db.chrom}:{db.pos}"
    return f"chr{db.chrom}:{db.pos}"


def _no_change_note(records: list[ClinVarRecord]) -> str | None:
    changing = [record for record in records if not record.no_change]
    idle = [record for record in records if record.no_change]
    if not changing or not idle:
        return None
    bits = []
    for record in idle:
        label = record.cdna or record.title or record.accession
        described = _descriptions(record) or ["no clinical significance"]
        bits.append(f"{record.accession} ({label}; {'; '.join(described)})")
    return (
        "ClinVar also links "
        + "; ".join(bits)
        + ". Those records report no sequence alteration, so they are not applied to the alternate allele."
    )


def _function(view: AlleleView, db: DbSnpRecord, prefix: str) -> str:
    records = view.records
    if not records:
        return f"{prefix}Function: not reported."
    keys = {_consequence_key(record) for record in records}
    gene_text, gene_conflict = _gene_text(records, db)
    if len(keys) == 1 and not gene_conflict:
        details = "; ".join(_function_detail(record) for record in records)
        if records[0].consequences:
            lead = f"{', '.join(records[0].consequences)} in {gene_text}"
        else:
            lead = f"no molecular consequence stated for {gene_text}"
        return f"{prefix}Function: {lead}. {details}."
    bits = []
    for record in records:
        consequence = ", ".join(record.consequences) if record.consequences else "consequence not stated"
        bits.append(f"{consequence} in {_one_gene(record, db)}, {_function_detail(record)}")
    if gene_conflict:
        bits.append(gene_text)
    return f"{prefix}Function: unresolved. " + "; ".join(bits) + "."


def _function_detail(record: ClinVarRecord) -> str:
    label = record.title or record.cdna or record.accession
    return f"{label} ({_cite(record)})"


def _consequence_key(record: ClinVarRecord) -> str:
    return ", ".join(item.lower() for item in record.consequences) or "not stated"


def _gene_text(records: list[ClinVarRecord], db: DbSnpRecord) -> tuple[str, bool]:
    genes: list[str] = []
    for record in records:
        for gene in record.genes:
            if gene not in genes:
                genes.append(gene)
    if genes and db.genes and not set(genes) & set(db.genes):
        return (
            f"gene {', '.join(genes)} (ClinVar) and {', '.join(db.genes)} (dbSNP)",
            True,
        )
    if genes:
        return ", ".join(genes), False
    if db.genes:
        return ", ".join(db.genes) + " (dbSNP)", False
    return "gene not reported", False


def _one_gene(record: ClinVarRecord, db: DbSnpRecord) -> str:
    if record.genes:
        return ", ".join(record.genes)
    if db.genes:
        return ", ".join(db.genes) + " (dbSNP)"
    return "gene not reported"


def _frequency(view: AlleleView, prefix: str) -> str:
    label = ""
    if not prefix:
        if view.variant is not None and view.variant.ref and view.variant.alt:
            label = f" for {view.variant.ref}>{view.variant.alt}"
        elif view.ref and view.alt:
            label = f" for {view.ref}>{view.alt}"
    if view.variant is None:
        return f"{prefix}Frequency{label}: not reported."
    unresolved, phrase = frequency_phrase(view.variant.callsets)
    if unresolved:
        return f"{prefix}Frequency{label}: unresolved. {phrase}."
    return f"{prefix}Frequency{label}: {phrase}."


def _disease(view: AlleleView, *, clinvar_retrieved: bool, prefix: str) -> str:
    # Disagreeing submissions replace the aggregate. A rollup such as
    # "drug response" or "Conflicting classifications" is not a second call.
    submitted = _submitted_conflict(view.records)
    if submitted:
        return f"{prefix}Disease relevance: unresolved. {submitted}."
    statements = _disease_statements(view.records)
    if not statements:
        if clinvar_retrieved:
            return f"{prefix}Disease relevance: no clinical significance was found (ClinVar)."
        return f"{prefix}Disease relevance: not reported."
    keys = {key for key, _, _ in statements}
    forced = any(flag for _, _, flag in statements)
    if len(keys) > 1 or forced:
        return f"{prefix}Disease relevance: unresolved. " + _unresolved_disease(view.records) + "."
    traits: list[str] = []
    cites: list[str] = []
    description = ""
    for record in view.records:
        for classification in _meaningful(record):
            parts = _split_description(classification.description)
            if len(parts) != 1:
                continue
            description = parts[0]
            for trait in _kept_traits(classification.traits):
                if trait not in traits:
                    traits.append(trait)
            cite = _cite_classification(record, classification)
            if cite not in cites:
                cites.append(cite)
    trait_text = _trait_text(traits)
    return f"{prefix}Disease relevance: {description}{trait_text} ({'; '.join(cites)})."


def _submitted_conflict(records: list[ClinVarRecord]) -> str | None:
    """Distinct submission interpretations when those calls disagree.

    None when no submission is stored or every stored interpretation is the same
    call. Review status can differ without making the call itself a conflict;
    the aggregate citation is kept in that case.
    """
    grouped: dict[tuple[str, str], dict] = {}
    order: list[tuple[str, str]] = []
    labels: set[str] = set()
    for record in records:
        for submission in record.submissions:
            for part in _split_description(submission.description):
                label = interpretation_label(part)
                if label is None:
                    continue
                labels.add(label)
                review_key = " ".join(submission.review_status.lower().split())
                key = (label, review_key)
                if key not in grouped:
                    grouped[key] = {"part": part, "record": record, "classification": submission, "traits": []}
                    order.append(key)
                bucket = grouped[key]
                for trait in _kept_traits(submission.traits):
                    if trait not in bucket["traits"]:
                        bucket["traits"].append(trait)
    if len(labels) <= 1:
        return None
    chunks: list[str] = []
    for key in order:
        bucket = grouped[key]
        piece = f"{bucket['part']} ({_cite_classification(bucket['record'], bucket['classification'])})"
        piece += _trait_text(bucket["traits"])
        if piece not in chunks:
            chunks.append(piece)
    return "; ".join(chunks)


def _unresolved_disease(records: list[ClinVarRecord]) -> str:
    """Keep each conflicting classification. Say the condition list once when it is shared."""
    items: list[tuple[str, str, list[str]]] = []
    for record in records:
        for classification in _meaningful(record):
            traits = _kept_traits(classification.traits)
            cite = _cite_classification(record, classification)
            for part in _split_description(classification.description):
                items.append((part, cite, traits))
    shared = len({tuple(traits) for _, _, traits in items}) <= 1
    chunks: list[str] = []
    for part, cite, traits in items:
        piece = f"{part} ({cite})"
        if not shared:
            piece += _trait_text(traits)
        if piece not in chunks:
            chunks.append(piece)
    text = "; ".join(chunks)
    if shared and items and items[0][2]:
        text += _trait_text(items[0][2])
    return text


def _disease_statements(records: list[ClinVarRecord]) -> list[tuple[str, str, bool]]:
    statements: list[tuple[str, str, bool]] = []
    for record in records:
        meaningful = _meaningful(record)
        if not meaningful:
            continue
        for classification in meaningful:
            parts = _split_description(classification.description)
            # Several classifications inside one ClinVar description disagree with each other.
            split = len(parts) > 1
            for part in parts:
                trait_text = _trait_text(_kept_traits(classification.traits))
                text = f"{part}{trait_text} ({_cite_classification(record, classification)})"
                conflict = split or "conflict" in part.lower()
                statements.append((part.lower(), text, conflict))
    return statements


def _meaningful(record: ClinVarRecord) -> list[Classification]:
    return [item for item in record.classifications if not empty_significance(item.description)]


def _descriptions(record: ClinVarRecord) -> list[str]:
    found = []
    for item in _meaningful(record):
        found.extend(_split_description(item.description))
    return found


def _split_description(description: str) -> list[str]:
    return [part.strip() for part in description.split(";") if part.strip()]


def _kept_traits(traits: list[str]) -> list[str]:
    skip = {"not provided", "not specified", "none", "see cases"}
    return [trait for trait in traits if trait.lower() not in skip]


def _trait_text(traits: list[str]) -> str:
    if not traits:
        return ""
    shown = traits[:3]
    text = ", conditions " + "; ".join(shown)
    extra = len(traits) - len(shown)
    if extra:
        text += f", and {extra} more condition names in the record"
    return text


def _cite(record: ClinVarRecord) -> str:
    reviews = [item.review_status for item in record.classifications if item.review_status]
    cite = f"ClinVar {record.accession}"
    if reviews:
        cite += f", {reviews[0]}"
    return cite


def _cite_classification(record: ClinVarRecord, classification: Classification) -> str:
    cite = f"ClinVar {classification.accession or record.accession}"
    if classification.review_status:
        cite += f", {classification.review_status}"
    return cite


def _fallback_gene(records: list[ClinVarRecord], db: DbSnpRecord) -> str | None:
    for record in records:
        if record.genes:
            return record.genes[0]
    if db.genes:
        return db.genes[0]
    return None


def _papers(papers: list[Paper]) -> str | None:
    if not papers:
        return None
    cites = []
    for paper in papers:
        first = paper.authors.split(",")[0].strip() or "Unknown authors"
        who = first + (" et al." if "," in paper.authors else "")
        title = paper.title.rstrip(".")
        year = f", {paper.year}" if paper.year else ""
        cites.append(f"{who}{year}, {title}, PMID:{paper.pmid}")
    return "Papers: " + "; ".join(cites) + " (Europe PMC)."


def _norm(chrom: str) -> str:
    text = chrom.removeprefix("chr")
    return "M" if text in {"M", "MT"} else text
