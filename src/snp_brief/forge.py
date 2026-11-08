"""Hypothesis Forge: what is not yet understood about a variant, and what to test next.

Stage 1 (rules, no LLM) sorts the variant into one category from evidence the
pipeline already has:

- well characterized: a clear, reviewed ClinVar pathogenicity call and primary papers
- conflicting evidence: ClinVar submissions at the top review level disagree
- uncertain mechanism: ClinVar records an association or risk allele, not a
  pathogenicity call, so the trait link is known but how the variant acts is not
- insufficient evidence: no clinical significance, or a call with no papers

Each category has a reason (with fact ids) and rule-based next steps.

Stages 2 and 3 (LLM, only for conflicting and uncertain) propose two or three
competing hypotheses. Each hypothesis must cite the numbered facts that support
it, name what evidence is missing, the observation that would tell it apart
from the others, and one experiment. Supporting statements are validated with
the same checker as the summary; a hypothesis with no valid support is dropped.
Experiments are then ranked by transparent rules, not by model-assigned
probabilities.

Well-characterized and insufficient-evidence variants do not get invented
hypotheses: the Forge says why.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

from snp_brief.brief import _kept_traits
from snp_brief.module_cache import load_result, result_key, save_result
from snp_brief.reconcile import OTHER, reconcile, side, stars
from snp_brief.summary import Fact, _main_view, _phrase_consequence, _plan, check, facts_block, render

WELL = "well characterized"
CONFLICTING = "conflicting evidence"
UNCERTAIN = "uncertain mechanism"
INSUFFICIENT = "insufficient evidence"

NEXT_STEPS = {
    WELL: [
        "No competing mechanism to test: the pathogenicity call is reviewed and agreed.",
        "Open questions, if any, are context-specific: penetrance, modifiers, or effect size in under-studied ancestry groups.",
        "Re-check when ClinVar or gnomAD release new data.",
    ],
    CONFLICTING: [
        "Compare the assertion criteria and conditions behind the disagreeing ClinVar submissions; they may describe different phenotypes.",
        "Look for functional studies of this allele (enzyme activity, protein stability) that could settle the call.",
        "Check whether the disagreement tracks ancestry, using the gnomAD population counts.",
    ],
    UNCERTAIN: [
        "Identify the target gene: the nearest gene in ClinVar is not necessarily the gene the variant regulates.",
        "Check expression QTL (GTEx) and regulatory annotation (Ensembl VEP, ENCODE) at this position.",
        "Check fine-mapping or linkage: the associated SNP may tag a nearby causal variant.",
    ],
    INSUFFICIENT: [
        "Gather prerequisite evidence before proposing mechanisms: association studies (GWAS Catalog), expression QTLs, and functional annotation.",
        "Confirm the allele is the one studied in any paper that names this rsID.",
        "Do not interpret disease relevance from protein function alone.",
    ],
}

FORGE_PROMPT = """You are Hypothesis Forge, a research assistant for genetic variants (SNPs).
You get numbered facts retrieved from dbSNP, ClinVar, gnomAD, UniProt and Europe PMC, and the open problem.
Propose 2 or 3 COMPETING, biologically plausible hypotheses that could explain the open problem.

Rules:
- Every item in "support" must be a statement that is directly true from the facts, with the fact ids it uses. Copy numbers exactly from the facts you cite.
- Mechanisms may go beyond the facts, but they are hypotheses: say "may" or "could", never state them as established.
- "gaps" lists evidence that is missing for this hypothesis.
- "discriminating_observation" is the observation that would support this hypothesis over the others.
- "experiment" is one concrete experiment or analysis to make that observation.
- "feasibility" is "high" (public data or a standard computational analysis), "medium" (a standard lab assay), or "low" (a new cohort, animal model, or large study).
- Do not diagnose, give medical advice, or address patients.

Reply with only this JSON object:
{"open_questions": ["..."],
 "hypotheses": [
   {"title": "short name",
    "mechanism": "one or two sentences",
    "support": [{"text": "statement from the facts", "facts": ["F1"]}],
    "gaps": ["..."],
    "discriminating_observation": "...",
    "experiment": "...",
    "feasibility": "high|medium|low"}
 ]}"""

FEASIBILITY_SCORE = {"high": 3, "medium": 2, "low": 1}


@dataclass
class Triage:
    category: str
    reason: str
    fact_ids: list[str]
    next_steps: list[str]
    open_problem: str
    wants_hypotheses: bool


@dataclass
class Hypothesis:
    title: str
    mechanism: str
    support: list[dict]  # {"text", "facts", "rendered"}
    gaps: list[str]
    discriminating_observation: str
    experiment: str
    feasibility: str
    rejected_support: list[dict] = field(default_factory=list)


@dataclass
class ForgeResult:
    rsid: str
    triage: Triage
    facts: list[Fact]
    open_questions: list[str] = field(default_factory=list)
    hypotheses: list[Hypothesis] = field(default_factory=list)
    ranking: list[dict] = field(default_factory=list)
    writer: str = "rules"
    note: str = ""
    created: str = ""

    def to_dict(self) -> dict:
        return {
            "triage": asdict(self.triage),
            "open_questions": self.open_questions,
            "hypotheses": [asdict(h) for h in self.hypotheses],
            "ranking": self.ranking,
            "writer": self.writer,
            "note": self.note,
        }


# ---------------------------------------------------------------- stage 1


def triage(item) -> tuple[Triage, list[Fact]]:
    plan = _plan(item)
    facts = plan.facts
    disease_ids = [i for i in (plan.disease, plan.disease_other, plan.dissent) if i]
    function_ids = [i for i in (plan.function, plan.uniprot) if i]
    view = _main_view(item)
    verdict = reconcile(view.records, _kept_traits) if view is not None and view.records else None
    gene = plan.gene or "the gene"

    if plan.unresolved_groups:
        return (
            Triage(
                CONFLICTING,
                "ClinVar submissions at the highest review level present disagree ("
                + ", ".join(plan.unresolved_groups)
                + ").",
                disease_ids,
                NEXT_STEPS[CONFLICTING],
                f"Why do equally reviewed ClinVar submissions disagree about {item.rsid} in {gene}, "
                "and which mechanism, if any, makes this allele disease-relevant?",
                True,
            ),
            facts,
        )
    if plan.no_significance or not plan.headline:
        return (
            Triage(
                INSUFFICIENT,
                "ClinVar reports no clinical significance"
                + (", and no primary paper names this rsID." if not plan.papers else "."),
                disease_ids + function_ids,
                NEXT_STEPS[INSUFFICIENT],
                f"There is too little variant-specific evidence to propose a mechanism for {item.rsid}.",
                False,
            ),
            facts,
        )
    if side(plan.headline) == OTHER:
        trait = ", ".join(plan.conditions) or "a trait"
        consequence = _phrase_consequence(plan.consequence)
        return (
            Triage(
                UNCERTAIN,
                f"ClinVar records {plan.headline.lower()} with {trait}, not a pathogenicity call, "
                f"and the variant is {consequence} in {gene}. The trait link is reported; how the variant acts is not.",
                function_ids + disease_ids,
                NEXT_STEPS[UNCERTAIN],
                f"How does {item.rsid} ({consequence} in {gene}) produce its association with {trait}?",
                True,
            ),
            facts,
        )
    if verdict is not None:
        reviewed = verdict.resolved and stars(verdict.tier_review) >= 1
    else:
        reviewed = any(
            stars(c.review_status) >= 1 for record in view.records for c in record.classifications
        ) if view is not None else False
    if reviewed and plan.papers:
        return (
            Triage(
                WELL,
                f"ClinVar's best-reviewed call is {plan.headline}"
                + (f" ({verdict.support})" if verdict and verdict.support else "")
                + f", and {len(plan.papers)} primary papers name this rsID.",
                disease_ids + plan.papers[:1],
                NEXT_STEPS[WELL],
                f"{item.rsid} is well characterized; no competing mechanism is needed.",
                False,
            ),
            facts,
        )
    why = "no primary paper names this rsID" if not plan.papers else "no submission has assertion criteria"
    return (
        Triage(
            INSUFFICIENT,
            f"ClinVar's call is {plan.headline}, but {why}.",
            disease_ids,
            NEXT_STEPS[INSUFFICIENT],
            f"The {plan.headline} call for {item.rsid} lacks variant-specific literature.",
            False,
        ),
        facts,
    )


# ---------------------------------------------------------------- stages 2 and 3


def forge(item, *, use_llm: bool = True, client=None, cache_dir=None, refresh: bool = False) -> ForgeResult:
    found, facts = triage(item)
    result = ForgeResult(item.rsid, found, facts)
    if not found.wants_hypotheses:
        result.note = "No hypotheses generated: " + (
            "the variant is well characterized." if found.category == WELL else "the evidence is too thin to propose a mechanism."
        )
        return result
    if not use_llm:
        result.note = "Stage 1 only. Connect an LLM (NEBIUS_API_KEY) to generate competing hypotheses."
        return result
    from snp_brief.llm import LLMClient, parse_json

    base = client or LLMClient.from_env()
    forge_client = base.for_forge() if base is not None and hasattr(base, "for_forge") else base
    block = facts_block(facts)
    model = getattr(forge_client, "model", "") if forge_client is not None else ""
    key = result_key("forge", model, block + found.open_problem)
    if not refresh:
        stored = load_result(cache_dir, item.rsid, "forge", key) if forge_client is not None else _any_cached(cache_dir, item.rsid, block + found.open_problem)
        if stored is not None:
            return _from_stored(result, stored, cached=True)
    if forge_client is None:
        result.note = "Stage 1 only. Connect an LLM (NEBIUS_API_KEY) to generate competing hypotheses."
        return result
    user = (
        f"Variant: {item.rsid}\nCategory: {found.category}\nOpen problem: {found.open_problem}\n\n"
        f"Facts:\n{block}\n\nReply with the JSON object."
    )
    try:
        raw = parse_json(forge_client.chat(FORGE_PROMPT, user, json_mode=True, max_tokens=4000))
    except Exception as exc:
        result.note = f"Hypothesis generation failed ({exc}). Showing stage 1 only."
        return result
    validated = validate(raw, facts)
    validated["writer"] = f"llm:{forge_client.model}"
    validated["facts_hash"] = _hash(block + found.open_problem)
    stored = save_result(cache_dir, item.rsid, "forge", key, validated)
    return _from_stored(result, stored, cached=False)


def validate(raw: dict, facts: list[Fact]) -> dict:
    """Keep hypotheses whose support passes the fact checker, then rank their experiments."""
    hypotheses: list[dict] = []
    for item in (raw.get("hypotheses") or [])[:3]:
        if not isinstance(item, dict):
            continue
        kept_support, dropped = [], []
        for claim in item.get("support") or []:
            if not isinstance(claim, dict):
                continue
            ids = [str(x).strip() for x in claim.get("facts") or [] if str(x).strip()]
            sentence = f"{str(claim.get('text', '')).strip().rstrip('.')} [{', '.join(ids)}]." if ids else str(claim.get("text", ""))
            ok, rejected = check([sentence], facts)
            if ok:
                kept_support.append({"text": claim.get("text", ""), "facts": ids, "rendered": render(ok, facts)})
            else:
                dropped.append({"text": claim.get("text", ""), "reason": rejected[0][1]})
        if not kept_support:
            continue
        feasibility = str(item.get("feasibility", "medium")).lower()
        if feasibility not in FEASIBILITY_SCORE:
            feasibility = "medium"
        hypotheses.append(
            {
                "title": str(item.get("title", "Hypothesis")).strip(),
                "mechanism": str(item.get("mechanism", "")).strip(),
                "support": kept_support,
                "gaps": [str(g) for g in item.get("gaps") or []][:4],
                "discriminating_observation": str(item.get("discriminating_observation", "")).strip(),
                "experiment": str(item.get("experiment", "")).strip(),
                "feasibility": feasibility,
                "rejected_support": dropped,
            }
        )
    return {
        "open_questions": [str(q) for q in raw.get("open_questions") or []][:4],
        "hypotheses": hypotheses,
        "ranking": rank(hypotheses),
    }


def rank(hypotheses: list[dict]) -> list[dict]:
    """Order experiments by visible rules: feasibility first, then how much evidence is missing.

    A feasible experiment that targets the hypothesis with the most open gaps
    resolves the most uncertainty for the least effort. Ties keep model order.
    """
    rows = []
    for index, hyp in enumerate(hypotheses):
        if not hyp.get("experiment"):
            continue
        feas = FEASIBILITY_SCORE[hyp["feasibility"]]
        gaps = len(hyp.get("gaps") or [])
        rows.append(
            {
                "experiment": hyp["experiment"],
                "tests": hyp["title"],
                "feasibility": hyp["feasibility"],
                "open_gaps": gaps,
                "score": feas * 10 + min(gaps, 5),
                "_order": index,
            }
        )
    rows.sort(key=lambda row: (-row["score"], row["_order"]))
    for position, row in enumerate(rows, 1):
        row["rank"] = position
        row.pop("_order")
    return rows


def _from_stored(result: ForgeResult, stored: dict, *, cached: bool) -> ForgeResult:
    result.open_questions = stored.get("open_questions", [])
    result.hypotheses = [Hypothesis(**h) for h in stored.get("hypotheses", [])]
    result.ranking = stored.get("ranking", [])
    result.writer = stored.get("writer", "llm") + (" (cached)" if cached else "")
    result.created = stored.get("created", "")
    if not result.hypotheses:
        result.note = "The model's hypotheses did not pass the fact check. Showing stage 1 only."
    return result


def _hash(text: str) -> str:
    import hashlib

    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def _any_cached(cache_dir, rsid: str, basis: str) -> dict | None:
    """Offline: reuse a stored Forge result for these exact facts, whatever model wrote it."""
    import json
    from pathlib import Path

    if cache_dir is None:
        return None
    path = Path(cache_dir) / rsid / "forge.json"
    if not path.is_file():
        return None
    try:
        stored = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return stored if stored.get("facts_hash") == _hash(basis) else None
