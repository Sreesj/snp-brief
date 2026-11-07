"""Reconcile ClinVar submissions by review status.

A one-star or unreviewed submission should not overturn what the best-reviewed
submissions agree on. The rule:

1. Flagged submissions are excluded (ClinVar itself flags them as outdated or
   unsupported). They are counted, not used.
2. Pathogenicity calls (pathogenic side, benign side, uncertain) are ranked by
   ClinVar review stars. The headline is what the highest tier present agrees
   on. If that tier itself disagrees across sides, the field is unresolved.
3. Lower-tier submissions on a different side are kept as named dissent.
4. Other assertion types (risk factor, drug response, protective, association,
   affects, established risk allele, ...) are not pathogenicity calls. They are
   reported as additional ClinVar assertions, never as a conflict.

An aggregate (VCV) classification reviewed by an expert panel or practice
guideline is ClinVar's own reconciliation, so it is used as the headline.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from snp_brief.models import Classification, ClinVarRecord, interpretation_label

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, multiple submitters": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "criteria provided": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
    "flagged submission": -1,
}

PATHOGENIC = "pathogenic"
BENIGN = "benign"
UNCERTAIN = "uncertain"
OTHER = "other"

SIDE_LABEL = {
    PATHOGENIC: "Pathogenic/Likely pathogenic",
    BENIGN: "Benign/Likely benign",
    UNCERTAIN: "Uncertain significance",
}


SIDE_WORD = {PATHOGENIC: "pathogenic", BENIGN: "benign", UNCERTAIN: "uncertain"}


def stars(review_status: str) -> int:
    return STARS.get(" ".join(review_status.lower().split()), 0)


def side(label: str) -> str:
    folded = " ".join(label.lower().split())
    if folded.startswith("conflicting"):
        return OTHER
    if "pathogenic" in folded:
        return PATHOGENIC
    if "benign" in folded:
        return BENIGN
    if folded.startswith("uncertain significance"):
        return UNCERTAIN
    return OTHER


@dataclass
class Call:
    label: str  # as submitted, e.g. "Pathogenic, low penetrance"
    side: str
    stars: int
    record: ClinVarRecord
    classification: Classification

    @property
    def cite(self) -> str:
        accession = self.classification.accession or self.record.accession
        review = self.classification.review_status
        return f"ClinVar {accession}, {review}" if review else f"ClinVar {accession}"

    @property
    def accession(self) -> str:
        return self.classification.accession or self.record.accession


@dataclass
class Verdict:
    resolved: bool
    headline: str = ""  # e.g. "Pathogenic"
    headline_cites: list[str] = field(default_factory=list)
    support: str = ""  # e.g. "24 of 24 criteria-provided submissions agree"
    traits: list[str] = field(default_factory=list)
    tallies: list[tuple[str, list[Call]]] = field(default_factory=list)  # when unresolved
    tier_review: str = ""
    dissent: list[Call] = field(default_factory=list)
    other: list[tuple[str, list[Call]]] = field(default_factory=list)
    flagged: int = 0


def reconcile(records: list[ClinVarRecord], kept_traits) -> Verdict | None:
    """Verdict from submitted interpretations, or None when no submission is stored."""
    calls: list[Call] = []
    flagged = 0
    for record in records:
        for submission in record.submissions:
            for part in _split(submission.description):
                if interpretation_label(part) is None:
                    continue
                rank = stars(submission.review_status)
                if rank < 0:
                    flagged += 1
                    continue
                calls.append(Call(part, side(part), rank, record, submission))
    expert = _expert_aggregate(records)
    if not calls and expert is None:
        return None
    pathogenicity = [call for call in calls if call.side != OTHER]
    other = _group_other([call for call in calls if call.side == OTHER])
    if expert is not None:
        record, classification = expert
        label = _split(classification.description)[0]
        head_side = side(label)
        dissent = [call for call in pathogenicity if call.side != head_side]
        agreeing = [call for call in pathogenicity if call.side == head_side]
        traits = _traits(kept_traits, [classification.traits] + [c.classification.traits for c in agreeing])
        return Verdict(
            resolved=True,
            headline=label,
            headline_cites=[f"ClinVar {record.accession}, {classification.review_status}"],
            support="",
            traits=traits,
            tier_review=classification.review_status,
            dissent=dissent,
            other=other,
            flagged=flagged,
        )
    if not pathogenicity:
        return Verdict(resolved=True, other=other, flagged=flagged)
    top = max(call.stars for call in pathogenicity)
    tier = [call for call in pathogenicity if call.stars == top]
    sides = {call.side for call in tier}
    tier_review = _tier_review(tier)
    if len(sides) == 1:
        head_side = next(iter(sides))
        labels = Counter(call.label[:1].upper() + call.label[1:].lower() for call in tier)
        headline = labels.most_common(1)[0][0]
        agreeing = [call for call in pathogenicity if call.side == head_side]
        dissent = [call for call in pathogenicity if call.side != head_side]
        count = "the only submission" if len(tier) == 1 else f"all {len(tier)} submissions"
        support = f"{count} at the highest review level present ({tier_review}) on the {SIDE_WORD[head_side]} side"
        if len(labels) > 1:
            others = ", ".join(f"{n} {label}" for label, n in labels.most_common()[1:])
            support += f"; {labels[headline]} {headline}, {others}"
        return Verdict(
            resolved=True,
            headline=headline,
            headline_cites=["ClinVar " + ", ".join([call.accession for call in tier[:2]] + list(dict.fromkeys(c.record.accession for c in tier)))],
            support=support,
            traits=_traits(kept_traits, [call.classification.traits for call in agreeing]),
            tier_review=tier_review,
            dissent=dissent,
            other=other,
            flagged=flagged,
        )
    tallies: list[tuple[str, list[Call]]] = []
    for name in (PATHOGENIC, BENIGN, UNCERTAIN):
        group = [call for call in tier if call.side == name]
        if group:
            tallies.append((SIDE_LABEL[name], group))
    tallies.sort(key=lambda item: -len(item[1]))
    return Verdict(
        resolved=False,
        tallies=tallies,
        tier_review=tier_review,
        dissent=[],
        other=other,
        flagged=flagged,
        traits=_traits(kept_traits, [call.classification.traits for call in tier]),
    )


def render(verdict: Verdict, prefix: str, trait_text) -> str:
    """One Disease relevance sentence from a verdict."""
    tail: list[str] = []
    if verdict.other:
        bits = []
        for label, group in verdict.other:
            best = max(group, key=lambda call: call.stars)
            traits = trait_text(_first_traits(group))
            count = f", {len(group)} submissions" if len(group) > 1 else ""
            bits.append(f"{label}{traits} ({best.cite}{count})")
        tail.append("Other ClinVar assertions: " + "; ".join(bits) + ".")
    if verdict.dissent:
        grouped: dict[str, list[Call]] = {}
        for call in verdict.dissent:
            grouped.setdefault(call.label, []).append(call)
        bits = []
        for label, group in grouped.items():
            best = max(group, key=lambda call: call.stars)
            count = f", {len(group)} submissions" if len(group) > 1 else ""
            bits.append(f"{label} ({best.cite}{count})")
        tail.append("Lower-review submissions that disagree: " + "; ".join(bits) + ".")
    if verdict.flagged:
        noun = "submission was" if verdict.flagged == 1 else "submissions were"
        tail.append(f"{verdict.flagged} flagged ClinVar {noun} excluded.")
    extra = (" " + " ".join(tail)) if tail else ""
    if not verdict.resolved:
        groups = []
        for label, group in verdict.tallies:
            examples = ", ".join(call.accession for call in group[:2])
            groups.append(f"{label} ({len(group)} submissions, e.g. ClinVar {examples})")
        body = (
            f"unresolved. ClinVar submissions at the highest review level present ({verdict.tier_review}) "
            f"disagree: " + "; ".join(groups) + trait_text(verdict.traits) + "."
        )
        return f"{prefix}Disease relevance: {body}{extra}"
    if not verdict.headline:
        return f"{prefix}Disease relevance: no pathogenicity classification in ClinVar.{extra}"
    cites = "; ".join(dict.fromkeys(verdict.headline_cites))
    support = f"; {verdict.support}" if verdict.support else ""
    return (
        f"{prefix}Disease relevance: {verdict.headline}{trait_text(verdict.traits)} "
        f"({cites}{support}).{extra}"
    )


def _expert_aggregate(records: list[ClinVarRecord]) -> tuple[ClinVarRecord, Classification] | None:
    for record in records:
        for classification in record.classifications:
            parts = _split(classification.description)
            if len(parts) != 1 or side(parts[0]) == OTHER:
                continue
            if stars(classification.review_status) >= 3:
                return record, classification
    return None


def _group_other(calls: list[Call]) -> list[tuple[str, list[Call]]]:
    grouped: dict[str, list[Call]] = {}
    names: dict[str, str] = {}
    for call in calls:
        key = call.label.lower()
        names.setdefault(key, call.label)
        grouped.setdefault(key, []).append(call)
    ordered = sorted(grouped.items(), key=lambda item: (-max(c.stars for c in item[1]), -len(item[1])))
    return [(names[key], group) for key, group in ordered]


def _tier_review(tier: list[Call]) -> str:
    reviews = Counter(" ".join(call.classification.review_status.split()) for call in tier)
    return reviews.most_common(1)[0][0] if reviews else ""


def _traits(kept_traits, lists: list[list[str]]) -> list[str]:
    seen: dict[str, str] = {}
    for traits in lists:
        for trait in kept_traits(traits):
            seen.setdefault(trait.lower(), trait)
    return list(seen.values())


def _first_traits(group: list[Call]) -> list[str]:
    seen: dict[str, str] = {}
    for call in group:
        for trait in call.classification.traits:
            if trait.lower() in {"not provided", "not specified", "none", "see cases"}:
                continue
            seen.setdefault(trait.lower(), trait)
    return list(seen.values())


def _split(description: str) -> list[str]:
    return [part.strip() for part in description.split(";") if part.strip()]
