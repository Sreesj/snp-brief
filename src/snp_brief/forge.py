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
REPORTED = "mechanism proposed in literature"
TRAIT_ASSOCIATION = "trait association, mechanism unknown"
EXPLORATORY = "exploratory (no variant-specific evidence)"

MECHANISM_WORDS = (
    "enhancer", "promoter", "chromatin", "splic", "transcription factor", "binding site",
    "expression", "regulat", "loop", "methylation", "eqtl",
)

NEXT_STEPS = {
    WELL: [
        "Open questions, if any, are context-specific: penetrance, modifiers, or effect size in under-studied ancestry groups.",
        "Re-check when ClinVar or gnomAD release new data.",
    ],
    TRAIT_ASSOCIATION: [
        "Identify the target gene: check expression QTLs (GTEx) and the genes reported by the GWAS studies.",
        "Check whether the GWAS signal is driven by this SNP or a linked one (fine-mapping, credible sets).",
        "Look for regulatory annotation at this position in the relevant tissue.",
    ],
    EXPLORATORY: [
        "Start with free public checks: GTEx expression effects, GWAS Catalog for linked SNPs, and regulatory annotation.",
        "Treat every hypothesis below as a starting point; none rests on evidence about this specific variant.",
        "If nothing turns up, the most likely explanation for a common variant is no functional effect.",
    ],
    REPORTED: [
        "Read the cited paper to confirm the proposed mechanism and the cell type it was tested in.",
        "Check whether independent data (expression QTLs, reporter or editing assays) replicate it.",
        "Ask what remains open: tissue specificity, effect size, or other variants on the same haplotype.",
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
- Paper facts contain only the paper's title. A support statement that cites a paper may only say what the title says; do not describe the study design or findings beyond it.
- ClinVar may list other names for the same variant. Never propose that the variant is linked to, or tags, a variant that is one of its own other names.
- If a fact (for example a paper title) says this variant itself has a mechanism, do not propose that it is only a tag for another variant unless you mark that hypothesis "contradicted".
- Do not name any rsID or gene that is not in the facts; say "a nearby variant" instead.
- "consistency" is "consistent" (no fact argues against it), "partly" (some facts fit poorly), or "contradicted" (a fact argues against it).
- "background" lists up to 3 statements from your general biology knowledge that the hypothesis relies on (for example what the nearby gene is known to do). They are shown to the user labelled as unverified AI background knowledge. Never put a number, rsID, or study result in background.
- If the category is exploratory, include one hypothesis that the variant has no functional effect, and prefer cheap first checks (public data lookups) as experiments.
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
    "feasibility": "high|medium|low",
    "consistency": "consistent|partly|contradicted",
    "background": ["..."]}
 ]}"""

CONSISTENCY_SCORE = {"consistent": 2, "partly": 1, "contradicted": 0}

FEASIBILITY_SCORE = {"high": 3, "medium": 2, "low": 1}


@dataclass
class Triage:
    category: str
    reason: str
    fact_ids: list[str]
    next_steps: list[str]
    open_problem: str
    wants_hypotheses: bool
    exploratory: bool = False


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
    consistency: str = "consistent"
    background: list[str] = field(default_factory=list)


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
        if plan.gwas or plan.eqtl:
            return _trait_association(item, plan, facts, function_ids)
        if _has_context(plan):
            return _exploratory(item, plan, facts, "ClinVar reports no clinical significance")
        return (
            Triage(
                INSUFFICIENT,
                "ClinVar reports no clinical significance"
                + (", and no primary paper names this rsID." if not plan.papers else ".")
                + " No GWAS association, location annotation, or prediction was retrieved either.",
                disease_ids + function_ids,
                NEXT_STEPS[INSUFFICIENT],
                f"There is too little evidence of any kind to propose a mechanism for {item.rsid}.",
                False,
            ),
            facts,
        )
    if side(plan.headline) == OTHER:
        trait = ", ".join(plan.conditions) or "a trait"
        consequence = _phrase_consequence(plan.consequence)
        mechanism_papers = [
            (fact_id, paper)
            for fact_id, paper in zip(plan.papers, item.papers)
            if any(word in paper.title.lower() for word in MECHANISM_WORDS)
        ]
        if mechanism_papers:
            fact_id, paper = mechanism_papers[0]
            return (
                Triage(
                    REPORTED,
                    f"ClinVar records {plan.headline.lower()} with {trait}, and a paper title proposes how the variant acts "
                    f'("{paper.title.rstrip(".")}", PMID:{paper.pmid}). Databases here do not record the mechanism itself.',
                    function_ids + disease_ids + [fact_id],
                    NEXT_STEPS[REPORTED],
                    f"A mechanism for {item.rsid} is proposed in the literature (PMID:{paper.pmid}). What would confirm it, "
                    "in which tissue, and what does it leave unexplained?",
                    True,
                ),
                facts,
            )
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
                _well_steps(plan),
                f"{item.rsid} is well characterized; no competing mechanism is needed.",
                False,
            ),
            facts,
        )
    why = "no primary paper names this rsID" if not plan.papers else "no submission has assertion criteria"
    if plan.gwas or plan.eqtl:
        return _trait_association(item, plan, facts, function_ids)
    if _has_context(plan):
        return _exploratory(item, plan, facts, f"ClinVar's call is {plan.headline}, but {why}")
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


def _has_context(plan) -> bool:
    return any(
        (plan.function, plan.location_fact, plan.regulatory, plan.predictions, plan.uniprot, plan.gene_summary, plan.l2g, plan.ld)
    )


def _trait_association(item, plan, facts, function_ids) -> tuple[Triage, list[Fact]]:
    where = _phrase_consequence(plan.consequence) + (f" in or near {plan.gene}" if plan.gene else "")
    found = []
    if plan.gwas:
        found.append(
            f"the GWAS Catalog links it to {', '.join(plan.gwas_traits[:3])} "
            f"({plan.gwas_count} association{'s' if plan.gwas_count != 1 else ''})"
        )
    if plan.eqtl:
        found.append(f"GTEx links it to expression of {', '.join(plan.eqtl_genes[:3])}")
    target = ", ".join(plan.gwas_traits[:3]) if plan.gwas else f"the expression of {', '.join(plan.eqtl_genes[:3])}"
    extra = [i for i in (plan.eqtl, plan.l2g, plan.ld, plan.regulatory) if i]
    return (
        Triage(
            TRAIT_ASSOCIATION,
            "No ClinVar classification, but " + " and ".join(found) + ". "
            "Statistical associations say nothing yet about mechanism.",
            function_ids + plan.gwas[:3] + extra,
            NEXT_STEPS[TRAIT_ASSOCIATION],
            f"How could {item.rsid} ({where}) affect {target}, and is it the causal variant at this locus"
            + (f" or one of the {plan.ld_count} variants inherited with it" if plan.ld_count else "")
            + "?",
            True,
        ),
        facts,
    )


def _exploratory(item, plan, facts, lead: str) -> tuple[Triage, list[Fact]]:
    known = []
    if plan.function or plan.location_fact:
        known.append("its location" + (f" relative to {plan.gene}" if plan.gene else ""))
    if plan.regulatory:
        known.append("regulatory annotation")
    if plan.predictions:
        known.append("computational predictions")
    if plan.uniprot or plan.gene_summary:
        known.append("the nearby gene's function")
    ids = [i for i in (plan.function, plan.location_fact, plan.regulatory, plan.predictions, plan.uniprot, plan.gene_summary) if i]
    return (
        Triage(
            EXPLORATORY,
            f"{lead}, and no GWAS association was found. Only {', '.join(known)} is known, "
            "so any hypothesis is exploratory.",
            ids,
            NEXT_STEPS[EXPLORATORY],
            f"What, if anything, does {item.rsid} do? Propose exploratory explanations from its genomic context, "
            "including the possibility that it has no functional effect.",
            True,
            exploratory=True,
        ),
        facts,
    )


def _well_steps(plan) -> list[str]:
    first = f"No competing mechanism to test: the ClinVar classification ({plan.headline}) comes from the highest review tier present"
    if plan.dissent_count:
        first += f"; {plan.dissent_count} lower-review submissions disagree and are listed in the evidence."
    else:
        first += " and those submissions agree."
    return [first, *NEXT_STEPS[WELL]]


# ---------------------------------------------------------------- stages 2 and 3


def forge(item, *, use_llm: bool = True, client=None, cache_dir=None, refresh: bool = False) -> ForgeResult:
    found, facts = triage(item)
    result = ForgeResult(item.rsid, found, facts)
    if not found.wants_hypotheses:
        result.note = "No hypotheses generated: " + (
            "the variant is well characterized." if found.category == WELL else "no evidence of any kind was retrieved."
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
        f"Variant: {item.rsid}\nCategory: {found.category}{' (exploratory)' if found.exploratory else ''}\n"
        f"Open problem: {found.open_problem}\n\n"
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
    import re

    hypotheses: list[dict] = []
    corpus = " ".join(fact.text for fact in facts)
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
        consistency = str(item.get("consistency", "consistent")).lower()
        if consistency not in CONSISTENCY_SCORE:
            consistency = "partly"
        scrub = lambda text: _scrub_rsids(str(text or "").strip(), corpus)  # noqa: E731
        hypotheses.append(
            {
                "title": scrub(item.get("title", "Hypothesis")),
                "mechanism": scrub(item.get("mechanism", "")),
                "support": kept_support,
                "gaps": [scrub(g) for g in item.get("gaps") or []][:4],
                "discriminating_observation": scrub(item.get("discriminating_observation", "")),
                "experiment": scrub(item.get("experiment", "")),
                "feasibility": feasibility,
                "rejected_support": dropped,
                "consistency": consistency,
                "background": [
                    scrub(b) for b in (item.get("background") or [])[:3] if b and not re.search(r"\d", str(b))
                ],
            }
        )
    return {
        "open_questions": [str(q) for q in raw.get("open_questions") or []][:4],
        "hypotheses": hypotheses,
        "ranking": rank(hypotheses),
    }


def rank(hypotheses: list[dict]) -> list[dict]:
    """Order experiments by visible rules, most important first:

    1. consistency of the tested hypothesis with the facts (contradicted ones go last),
    2. how many fact-checked statements support it,
    3. feasibility of the experiment,
    4. how many evidence gaps it would address.

    Plausibility comes before ease, so an easy test of a weak idea does not rank first.
    """
    rows = []
    for index, hyp in enumerate(hypotheses):
        if not hyp.get("experiment"):
            continue
        consistency = hyp.get("consistency", "consistent")
        rows.append(
            {
                "experiment": hyp["experiment"],
                "tests": hyp["title"],
                "consistency": consistency,
                "support": len(hyp.get("support") or []),
                "feasibility": hyp["feasibility"],
                "open_gaps": len(hyp.get("gaps") or []),
                "_key": (
                    -CONSISTENCY_SCORE.get(consistency, 1),
                    -len(hyp.get("support") or []),
                    -FEASIBILITY_SCORE[hyp["feasibility"]],
                    -len(hyp.get("gaps") or []),
                    index,
                ),
            }
        )
    rows.sort(key=lambda row: row["_key"])
    for position, row in enumerate(rows, 1):
        row["rank"] = position
        row.pop("_key")
    return rows


def _scrub_rsids(text: str, corpus: str) -> str:
    """Replace rsIDs the facts do not contain: the model must not introduce variants."""
    import re

    return re.sub(
        r"\brs\d+\b",
        lambda m: m.group(0) if m.group(0) in corpus else "a nearby variant (not named in the retrieved evidence)",
        text,
    )


def _from_stored(result: ForgeResult, stored: dict, *, cached: bool) -> ForgeResult:
    result.open_questions = stored.get("open_questions", [])
    result.hypotheses = [Hypothesis(**{k: v for k, v in h.items() if k in Hypothesis.__dataclass_fields__}) for h in stored.get("hypotheses", [])]
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
