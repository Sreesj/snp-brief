"""Variant Atlas: search or run a list of SNPs, then investigate one in three modules.

Variant Observatory is the default module (evidence and summary). Molecular
Ripple and Hypothesis Forge run only when someone presses Run in their tab;
results are kept per SNP in the session and in cache/<rsid>/.
"""

from __future__ import annotations

import hmac
import html
import os
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from snp_brief.demo import motion, theme
from snp_brief.interpret import NOT_COMPUTED, acmg_draft, genotype_rates
from snp_brief.demo.view import (
    PopulationTable,
    briefs_text,
    gene_label,
    paper_links,
    parse_rsid_text,
    source_rows,
    unresolved_populations,
)
from snp_brief.forge import CONFLICTING, EXPLORATORY, INSUFFICIENT, REPORTED, TRAIT_ASSOCIATION, UNCERTAIN, WELL, forge
from snp_brief.graph import evidence_dot
from snp_brief.literature import Finding, extract_findings
from snp_brief.llm import LLMClient
from snp_brief.pipeline import Gathered, SourceCall, default_user_cache, in_flight, load_variant
from snp_brief.providers import BY_KEY, SPEED_BUDGETS
from snp_brief.report import build_report, load_notes, save_notes
from snp_brief.ripple import ripple_cached
from snp_brief.summary import _plan, summarize

USER_CACHE = default_user_cache()
MODULES = ("Variant Observatory", "Interpretation", "Molecular Ripple", "Hypothesis Forge")
CATEGORY_KIND = {
    WELL: "benign",
    CONFLICTING: "uncertain",
    UNCERTAIN: "association",
    TRAIT_ASSOCIATION: "association",
    REPORTED: "brand",
    EXPLORATORY: "neutral",
    INSUFFICIENT: "neutral",
}
NAMED_CSS = "color: #b9beff; font-weight: 600"
SMALL_CSS = "color: #62666d"


@dataclass
class Run:
    number: int
    live: bool
    results: list[Gathered]
    rejected: list[str]
    failed: list[str]


# ---------------------------------------------------------------- setup


def load_secrets() -> None:
    """On Streamlit Cloud, settings live in st.secrets; copy the known ones into the environment."""
    try:
        secrets = dict(st.secrets)
    except Exception:  # no secrets file locally
        return
    for name in (
        "NEBIUS_API_KEY",
        "SNP_BRIEF_LLM_MODEL",
        "SNP_BRIEF_FORGE_MODEL",
        "SNP_BRIEF_LLM_URL",
        "SNP_BRIEF_LLM_KEY",
        "SNP_BRIEF_USER_CACHE",
        "SNP_BRIEF_APP_PASSWORD",
        "SNP_BRIEF_FORGE_LIMIT",
        "NCBI_API_KEY",
        "SNP_BRIEF_EXTRAS",
        "SNP_BRIEF_MAX_AGE_DAYS",
    ):
        if name in secrets and not os.environ.get(name):
            os.environ[name] = str(secrets[name])


def password_ok() -> bool:
    """If SNP_BRIEF_APP_PASSWORD is set, ask for it once per session."""
    expected = os.environ.get("SNP_BRIEF_APP_PASSWORD")
    if not expected or st.session_state.get("unlocked"):
        return True
    st.html(theme.top_bar())
    entered = st.text_input("Password", type="password", key="password")
    if entered:
        if hmac.compare_digest(entered, expected):
            st.session_state.unlocked = True
            st.rerun()
        st.error("Wrong password.")
    return False


def forge_runs_left() -> int:
    limit = int(os.environ.get("SNP_BRIEF_FORGE_LIMIT", "20"))
    return max(0, limit - st.session_state.get("forge_llm_runs", 0))


def main() -> None:
    st.set_page_config(page_title="Variant Atlas", page_icon=":material/genetics:", layout="wide")
    st.html(theme.FONTS + theme.CSS + theme.EXTRA_CSS)
    load_secrets()
    if not password_ok():
        return
    root = repo_root()
    offline = sidebar()
    current: Run | None = st.session_state.get("last_run")
    if current is None or not current.results:
        # Landing: a centred hero and one search box. Cached SNPs still load instantly when searched.
        st.html(theme.top_bar() + theme.hero())
        _, middle, _ = st.columns([1, 3, 1])
        with middle:
            search_bar(root, offline)
            st.html(theme.features())
        st.html(theme.about())
        install_motion("landing")
        if st.session_state.get("last_run") is not None and st.session_state.last_run.results:
            st.rerun()
        if current is not None:
            show_results(current)
        return
    st.html(theme.top_bar())
    install_motion("app")
    search_bar(root, offline)
    show_results(current)


def install_motion(mode: str) -> None:
    """The moving helix needs a script, which st.html does not run; a 1-pixel iframe does."""
    code = motion.injector(mode)
    if hasattr(st, "iframe"):
        with st.container(key="va-motion"):
            st.iframe(code, height=1)
    else:
        components.html(code, height=0)


def repo_root() -> Path:
    for candidate in (Path.cwd(), Path(__file__).resolve().parents[3]):
        if (candidate / "fixture" / "snps.txt").is_file():
            return candidate
    return Path.cwd()


def sidebar() -> bool:
    with st.sidebar:
        st.html(theme.top_bar() + theme.section("What's inside") + theme.sidebar_guide())
        with st.expander("Settings"):
            st.radio(
                "Speed",
                list(SPEED_BUDGETS),
                index=1,
                key="speed",
                horizontal=True,
                help="How long to wait for slower databases. Late answers are added when you press Check again.",
            )
            offline = st.toggle("Offline mode", key="offline", help="Use saved data only; no network calls.")
            if LLMClient.from_env() is not None:
                st.toggle("AI-written summaries", key="use_llm", value=True)
    return offline


# ---------------------------------------------------------------- loading


def fetch(root: Path, rsid: str, *, offline: bool, refresh: bool = False, budget: float | None = None) -> Gathered:
    if budget is None:
        budget = SPEED_BUDGETS.get(st.session_state.get("speed", "balanced"), 20.0)
    return load_variant(
        rsid, repo_cache=root / "cache", user_cache=USER_CACHE, offline=offline, refresh=refresh, budget=budget
    )


def reload_variant(item: Gathered) -> None:
    """Pick up providers that finished in the background (cache first, nothing re-fetched)."""
    current: Run | None = st.session_state.get("last_run")
    fresh = fetch(repo_root(), item.rsid, offline=bool(st.session_state.get("offline")), budget=4.0)
    if current is not None:
        current.results = [fresh if old.rsid == item.rsid else old for old in current.results]
    _forget_modules(item.rsid)


def _request_search() -> None:
    st.session_state.search_go = True


def search_bar(root: Path, offline: bool) -> None:
    """One box for one or many rsIDs (separated by commas, spaces or new lines)."""
    left, right = st.columns([6, 1], vertical_alignment="bottom")
    query = left.text_input(
        "Search variants",
        placeholder="Search one or more rsIDs, e.g. rs6025, rs429358",
        key="search_text",
        label_visibility="collapsed",
        on_change=_request_search,
    )
    clicked = right.button("Search", key="search_go_btn", type="primary", width="stretch")
    if not (clicked or st.session_state.pop("search_go", False)):
        return
    ids, rejected = parse_rsid_text(query or "")
    if not ids:
        st.error("Enter one or more rsIDs, such as rs6025.")
        return
    search(root, ids, rejected, offline)


def search(root: Path, ids: list[str], rejected: list[str], offline: bool) -> None:
    current: Run | None = st.session_state.get("last_run")
    if current is None:
        st.session_state.run_count = st.session_state.get("run_count", 0) + 1
        current = Run(st.session_state.run_count, not offline, [], [], [])
        st.session_state.last_run = current
    current.rejected = rejected
    current.failed = []
    placeholder = st.empty()
    found: list[Gathered] = []
    for index, rsid in enumerate(ids):
        existing = next((item for item in current.results if rsid in {item.rsid, item.merged_from}), None)
        if existing is None:
            if len(ids) > 1:
                placeholder.html(theme.loader(f"Analyzing {rsid}", f"{index + 1} of {len(ids)} variants", index, len(ids), skeleton=False))
            else:
                placeholder.html(theme.loader(f"Analyzing {rsid}", "Querying 19 sources in parallel"))
            try:
                existing = fetch(root, rsid, offline=offline)
            except Exception as exc:
                current.failed.append(f"{rsid} could not be loaded ({exc})")
                continue
        found.append(existing)
    placeholder.empty()
    others = [item for item in current.results if all(item is not new for new in found)]
    current.results = found + others
    if found:
        st.session_state.active = found[0].rsid


def refresh_variant(item: Gathered) -> list[str]:
    """Re-fetch one SNP live; keep the previous evidence if any source fails."""
    current: Run | None = st.session_state.get("last_run")
    st.session_state.setdefault("pickup_tries", {}).pop(item.rsid, None)
    fresh = fetch(repo_root(), item.rsid, offline=False, refresh=True)
    failed = [name for name, call in fresh.calls.items() if call.origin == "failed"]
    if failed:
        return failed
    if current is not None:
        current.results = [fresh if old.rsid == item.rsid else old for old in current.results]
    _forget_modules(item.rsid)
    return []


def _forget_modules(rsid: str) -> None:
    modules = st.session_state.setdefault("modules", {})
    for key in [key for key in modules if key[0] == rsid]:
        modules.pop(key)


def plan_for(item: Gathered):
    store = st.session_state.setdefault("plans", {})
    key = (item.rsid, item.text, len(getattr(item, "findings", []) or []))
    if key not in store:
        store[key] = _plan(item)
    return store[key]


# ---------------------------------------------------------------- results


def show_results(current: Run) -> None:
    if current.rejected:
        st.warning("Skipped (not rsIDs): " + ", ".join(current.rejected))
    for problem in current.failed:
        st.error(problem)
    if not current.results:
        return
    by_id = {item.rsid: item for item in current.results}
    if st.session_state.get("active") not in by_id:
        st.session_state.active = current.results[0].rsid
    if len(current.results) == 1:
        render_detail(current.results[0])
        return
    list_side, detail_side = st.columns([1, 3.2], gap="large")
    with list_side:
        st.html(f'<div class="va-list-title">Variants · {len(current.results)}</div>')
        for item in current.results:
            plan = plan_for(item)
            value, _, _ = theme.clinical(plan)
            gene = gene_label(item).split(",")[0]
            label = f"{item.rsid}  {gene}".strip()
            active = item.rsid == st.session_state.active
            if st.button(
                label,
                key=f"pick-{item.rsid}",
                help=value,
                type="primary" if active else "secondary",
                width="stretch",
            ):
                st.session_state.active = item.rsid
                st.rerun()
        if st.button("Clear list", key="clear_list", width="stretch"):
            st.session_state.last_run = None
            st.session_state.pop("active", None)
            st.rerun()
        st.download_button(
            "Download all briefs",
            briefs_text(current.results),
            file_name="variant-atlas-briefs.txt",
            mime="text/plain",
            on_click="ignore",
            width="stretch",
        )
    with detail_side:
        render_detail(by_id[st.session_state.active])


def render_detail(item: Gathered) -> None:
    plan = plan_for(item)
    gene = gene_label(item)
    st.html(theme.header_card(item, plan, gene, last_retrieved(item)))
    actions = st.columns([1, 1, 4])
    if actions[0].button("Refresh data", key=f"refresh-{item.rsid}", help="Fetch every source again"):
        placeholder = st.empty()
        placeholder.html(theme.loader(f"Refreshing {item.rsid}", "Fetching every source again"))
        failed = refresh_variant(item)
        placeholder.empty()
        if failed:
            st.error("Some sources did not answer; the previous data is kept.")
        else:
            st.rerun()
    if item.pending and actions[1].button("Check again", key=f"again-{item.rsid}", help="Add sources that were still loading"):
        reload_variant(item)
        st.rerun()
    if item.pending:
        pick_up_late_sources(item.rsid)
    observatory, interpretation_tab, ripple_tab, forge_tab = st.tabs(list(MODULES))
    with observatory:
        render_observatory(item, plan)
    with interpretation_tab:
        render_interpretation(item)
    with ripple_tab:
        render_ripple(item)
    with forge_tab:
        render_forge(item)
    render_notebook(item)


@st.fragment(run_every=3)
def pick_up_late_sources(rsid: str) -> None:
    """While slow databases are still answering in the background, check every few seconds and
    redraw the page as soon as they land, so nobody has to press Check again."""
    current: Run | None = st.session_state.get("last_run")
    item = next((x for x in (current.results if current else []) if x.rsid == rsid), None)
    if item is None or not item.pending:
        return
    tries = st.session_state.setdefault("pickup_tries", {})
    if in_flight(item.rsid, item.pending):
        names = ", ".join((BY_KEY[key].name if key in BY_KEY else key) for key in item.pending)
        st.html(f'<div class="va-waiting"><i></i>Still waiting on {html.escape(names)}. This page updates by itself.</div>')
        return
    tries[rsid] = tries.get(rsid, 0) + 1
    if tries[rsid] > 3:
        # Never leave a row on "Loading": whatever has not answered by now is unavailable.
        for key in item.pending:
            item.calls[key] = SourceCall("failed", error="no answer after several checks")
        st.rerun(scope="app")
    reload_variant(item)
    st.rerun(scope="app")


def run_module(item: Gathered, key: str, work) -> None:
    """Run a module, recording failure so the tab can show it."""
    results = st.session_state.setdefault("modules", {})
    failures = st.session_state.setdefault("module_errors", {})
    try:
        results[(item.rsid, key)] = work()
        failures.pop((item.rsid, key), None)
    except Exception as exc:
        failures[(item.rsid, key)] = str(exc)


def last_retrieved(item: Gathered) -> str:
    stamps = sorted(call.retrieved for call in item.calls.values() if call.retrieved)
    return stamps[-1][:10] if stamps else ""


def module_cache(item: Gathered) -> Path:
    return item.cache_dir or repo_root() / "cache"


# ---------------------------------------------------------------- Variant Observatory


def render_observatory(item: Gathered, plan=None) -> None:
    plan = plan or plan_for(item)
    st.html(theme.facts_grid(plan))
    mode = "llm" if st.session_state.get("use_llm") else "template"
    summary = cached_summary(item, mode)
    st.html(theme.section("Summary") + theme.summary_html(summary.sentences, summary.facts))
    tables = unresolved_populations(item)
    if tables:
        with st.expander("Population frequencies"):
            for table in tables:
                render_populations(table)
    with st.expander("Literature"):
        render_timeline(item)
        render_findings(item)
        links = paper_links(item)
        if links:
            st.markdown(
                "\n".join(
                    f"{number}. [{md_escape(link.title)}]({link.url})  \n   <span style='color:#5f6b7a;font-size:.85rem'>{html.escape(link.citation)}</span>"
                    for number, link in enumerate(links, 1)
                ),
                unsafe_allow_html=True,
            )
        else:
            st.caption("No research article names this rsID in its title or abstract.")
    with st.expander("Evidence details"):
        st.html(theme.section("Evidence brief") + f'<div class="va-card" style="font-size:.9rem;line-height:1.6">{html.escape(item.text)}</div>')
        st.html(theme.section("Data sources"))
        rows = source_rows(item)
        frame = pd.DataFrame(
            [
                {
                    "Source": row["Source"],
                    "Kind": row["Used for"].split("·", 1)[-1].strip(),
                    "Status": {"yes": "Loaded", "no": "No data", "loading": "Loading…", "unavailable": "Unavailable", "not consulted": "Not needed"}.get(
                        row["Answered"], row["Answered"]
                    ),
                    "Result": row["Gave"],
                }
                for row in rows
            ]
        )
        st.table(frame, hide_index=True, border="horizontal")


def cached_summary(item: Gathered, mode: str):
    store = st.session_state.setdefault("summaries", {})
    key = (item.rsid, mode, item.text)
    if key not in store:
        store[key] = summarize(item, mode, cache_dir=module_cache(item))
    return store[key]


def render_findings(item: Gathered) -> None:
    store = st.session_state.setdefault("findings", {})
    client = LLMClient.from_env()
    found = store.get(item.rsid)
    if found is None and item.findings:
        found = {"findings": [vars(f) for f in item.findings]}
    if found is None:
        if client is not None and st.button("Extract findings from abstracts", key=f"findings-{item.rsid}"):
            placeholder = st.empty()
            placeholder.html(theme.loader("Reading abstracts", "Every finding must quote its abstract word for word"))
            try:
                found = extract_findings(
                    item, client.for_forge(), module_cache(item), offline=bool(st.session_state.get("offline"))
                )
            except Exception:
                found = {"findings": []}
            placeholder.empty()
            store[item.rsid] = found
            item.findings = [Finding(**f) for f in found.get("findings", [])]
            st.session_state.get("summaries", {}).clear()
            st.session_state.get("plans", {}).clear()
            st.session_state.get("modules", {}).pop((item.rsid, "forge"), None)
            st.rerun()
        return
    findings = found.get("findings", [])
    if findings:
        st.html(theme.section("Findings from abstracts"))
        for f in findings:
            st.markdown(f'**{f["claim"]}** ({_pmid_link(f["pmid"])})  \n> {f["quote"]}', unsafe_allow_html=True)


def _pmid_link(pmid: str) -> str:
    return f'<a href="https://europepmc.org/article/MED/{html.escape(str(pmid))}" target="_blank">PMID {html.escape(str(pmid))}</a>'


def render_populations(table: PopulationTable) -> None:
    st.caption(f"{table.allele} · {' · '.join(table.totals)}")
    frame = pd.DataFrame(table.rows, columns=["Population", *table.callsets])

    def highlight(data: pd.DataFrame) -> pd.DataFrame:
        styles = pd.DataFrame("", index=data.index, columns=data.columns)
        for index, label in data["Population"].items():
            for name in table.callsets:
                if (label, name) in table.named:
                    styles.at[index, name] = NAMED_CSS
                elif (label, name) in table.small:
                    styles.at[index, name] = SMALL_CSS
        return styles

    st.dataframe(frame.style.apply(highlight, axis=None), hide_index=True, width="stretch")
    st.caption("Allele count / allele number. Highlighted: highest and lowest ancestry groups.")


# ---------------------------------------------------------------- Molecular Ripple


def render_timeline(item: Gathered) -> None:
    t = item.timeline
    if not t.by_year:
        return
    years = range(min(t.by_year), max(t.by_year) + 1)
    frame = pd.DataFrame({"Papers": [t.by_year.get(y, 0) for y in years]}, index=[str(y) for y in years])
    trend = {"growing": "interest is growing", "declining": "interest is declining", "steady": "interest is steady"}.get(t.trend, "")
    st.html(theme.section(f"Papers per year · {t.total} in Europe PMC" + (f" · {trend}" if trend else "")))
    st.bar_chart(frame, color="#7c86ff", height=180, x_label="", y_label="")


def render_interpretation(item: Gathered) -> None:
    """Computed aids for the people who act on a variant: ACMG draft, genotype rates, drugs, cancer."""
    st.html(theme.section("ACMG/AMP evidence draft") + theme.acmg_html(acmg_draft(item), NOT_COMPUTED))
    rates = genotype_rates(item)
    st.html(theme.section("How common is it in people"))
    if rates:
        frame = pd.DataFrame(
            [
                {
                    "Group": r.group,
                    "Allele": r.allele,
                    "Allele frequency": f"{r.af * 100:.3g}%",
                    "Carry at least one copy": r.carriers_text,
                    "Carry two copies": r.homozygous_text,
                }
                for r in rates
            ]
        )
        st.table(frame, hide_index=True, border="horizontal")
        st.caption("Estimated from gnomAD allele frequencies assuming Hardy-Weinberg equilibrium; groups with fewer than 2,000 alleles are left out.")
    else:
        st.caption("gnomAD has no frequency for this variant.")
    pgx = theme.pgx_html(item.pgx)
    st.html(theme.section("Drug response (ClinPGx)") + (pgx or '<div class="va-acmg-empty">No curated drug-response annotation for this variant.</div>'))
    civic = theme.civic_html(item.scores.civic)
    st.html(theme.section("Cancer evidence (CIViC)") + (civic or '<div class="va-acmg-empty">No CIViC evidence for this variant.</div>'))


def render_ripple(item: Gathered) -> None:
    grid = theme.score_grid(item)
    if grid:
        st.html(theme.section("Prediction and constraint scores") + grid)
    results = st.session_state.setdefault("modules", {})
    key = (item.rsid, "ripple")
    if key not in results:
        st.html(theme.callout("Molecular context for this variant: what it touches, what is predicted, and what is still unknown.", "info"))
        error = st.session_state.get("module_errors", {}).get(key)
        if error:
            st.error("Molecular Ripple could not run for this variant.")
        if st.button("Run Molecular Ripple", key=f"run-ripple-{item.rsid}", type="primary"):
            placeholder = st.empty()
            placeholder.html(theme.loader("Mapping molecular context"))
            run_module(item, "ripple", lambda: ripple_cached(item, module_cache(item)))
            placeholder.empty()
            st.rerun()
        return
    found = results[key]
    st.markdown(f"**Variant context:** {found.context}" + (f" · {found.gene}" if found.gene else ""))
    st.html(theme.section("What to examine"))
    st.markdown("\n".join(f"- {line}" for line in found.priorities))
    st.html(theme.section("Molecular links") + theme.link_rows(found.links))
    st.html(theme.section("Evidence graph"))
    try:
        st.graphviz_chart(evidence_dot(item, results.get((item.rsid, "forge"))), width="stretch")
        st.html(
            '<div class="va-chips">'
            + theme.chip("observed", "benign")
            + theme.chip("predicted", "association")
            + theme.chip("hypothesis", "neutral")
            + "</div>"
        )
    except Exception:
        st.caption("The graph could not be drawn here.")


# ---------------------------------------------------------------- Hypothesis Forge


def render_forge(item: Gathered) -> None:
    results = st.session_state.setdefault("modules", {})
    key = (item.rsid, "forge")
    client = LLMClient.from_env()
    if client is not None and forge_runs_left() == 0:
        client = None
    if key not in results:
        st.html(
            theme.callout(
                "Finds what is not yet understood about this variant and, where useful, proposes competing hypotheses "
                "and the experiments that would tell them apart.",
                "info",
            )
        )
        if st.session_state.get("module_errors", {}).get(key):
            st.error("Hypothesis Forge could not run for this variant.")
        if st.button("Run Hypothesis Forge", key=f"run-forge-{item.rsid}", type="primary"):
            placeholder = st.empty()
            placeholder.html(theme.loader("Forging hypotheses" if client else "Assessing the evidence"))
            run_module(item, "forge", lambda: forge(item, client=client, cache_dir=module_cache(item)))
            placeholder.empty()
            found = results.get(key)
            if found is not None and found.writer.startswith("llm:") and "(cached)" not in found.writer:
                st.session_state.forge_llm_runs = st.session_state.get("forge_llm_runs", 0) + 1
            st.rerun()
        return
    found = results[key]
    tri = found.triage
    category = tri.category[:1].upper() + tri.category[1:]
    st.html(f'<div class="va-chips" style="margin:0 0 .6rem 0">{theme.chip(category, CATEGORY_KIND.get(tri.category, "neutral"))}</div>')
    st.markdown(f"**Why:** {tri.reason}")
    st.markdown(f"**Open problem:** {tri.open_problem}")
    if tri.exploratory:
        st.html(
            theme.callout(
                "Exploratory: there is no evidence about this specific variant, so these hypotheses come from its genomic "
                "context and may be wrong. Treat them as starting points."
            )
        )
    st.html(theme.section("Next steps"))
    st.markdown("\n".join(f"- {line}" for line in tri.next_steps))
    if found.open_questions:
        st.html(theme.section("Open questions"))
        st.markdown("\n".join(f"- {q}" for q in found.open_questions))
    if found.hypotheses:
        st.html(theme.section("Competing hypotheses"))
        columns = st.columns(len(found.hypotheses))
        for number, (column, hyp) in enumerate(zip(columns, found.hypotheses), 1):
            with column, st.container(border=True):
                st.markdown(f"**H{number}. {hyp.title}**")
                st.html(
                    '<div class="va-chips" style="margin:.1rem 0 .4rem 0">'
                    + theme.chip(f"{hyp.feasibility} feasibility", "neutral")
                    + theme.chip(hyp.consistency, {"consistent": "benign", "partly": "uncertain", "contradicted": "pathogenic"}.get(hyp.consistency, "neutral"))
                    + "</div>"
                )
                st.markdown(hyp.mechanism)
                st.markdown("**For**  \n" + "\n".join(f"- {claim['rendered']}" for claim in hyp.support))
                if hyp.against:
                    st.markdown("**Against**  \n" + "\n".join(f"- {claim['rendered']}" for claim in hyp.against))
                if hyp.gaps:
                    st.markdown("**Missing evidence**  \n" + "\n".join(f"- {gap}" for gap in hyp.gaps))
                if hyp.background:
                    st.markdown("**Background (AI, unverified)**  \n" + "\n".join(f"- {line}" for line in hyp.background))
                st.markdown(f"**Deciding observation:** {hyp.discriminating_observation}")
                st.markdown(f"**Experiment:** {hyp.experiment}")
    if found.ranking:
        st.html(theme.section("Suggested experiments, in order"))
        frame = pd.DataFrame(found.ranking)[["rank", "experiment", "tests", "feasibility"]]
        frame.columns = ["#", "Experiment", "Tests", "Feasibility"]
        st.dataframe(frame, hide_index=True, width="stretch")
    if found.hypotheses and client is not None:
        with st.expander("Refine the hypotheses"):
            rejected = [
                hyp.title
                for number, hyp in enumerate(found.hypotheses, 1)
                if st.checkbox(f"Reject H{number}: {hyp.title}", key=f"reject-{item.rsid}-{number}")
            ]
            notes = st.text_area(
                "Constraints or resources", placeholder="e.g. computational only, we have intestinal organoids",
                key=f"steer-{item.rsid}",
            )
            if st.button("Regenerate", key=f"reforge-{item.rsid}") and (rejected or notes.strip()):
                placeholder = st.empty()
                placeholder.html(theme.loader("Regenerating hypotheses"))
                run_module(
                    item,
                    "forge",
                    lambda: forge(item, client=client, cache_dir=module_cache(item), steer={"rejected": rejected, "notes": notes}),
                )
                placeholder.empty()
                st.session_state.forge_llm_runs = st.session_state.get("forge_llm_runs", 0) + 1
                st.rerun()


# ---------------------------------------------------------------- notes and export


def render_notebook(item: Gathered) -> None:
    with st.expander("Notes and export"):
        notes_key = f"notes-{item.rsid}"
        if notes_key not in st.session_state:
            st.session_state[notes_key] = load_notes(USER_CACHE, item.rsid)
        st.text_area("Notes", key=notes_key, height=120, placeholder="Your notes on this variant")
        left, right = st.columns(2)
        if left.button("Save notes", key=f"save-notes-{item.rsid}", width="stretch"):
            ok = save_notes(USER_CACHE, item.rsid, st.session_state[notes_key])
            st.toast("Notes saved" if ok else "Notes kept for this session")
        results = st.session_state.get("modules", {})
        mode = "llm" if st.session_state.get("use_llm") else "template"
        findings = st.session_state.get("findings", {}).get(item.rsid, {}).get("findings")
        report = build_report(
            item,
            cached_summary(item, mode),
            results.get((item.rsid, "ripple")),
            results.get((item.rsid, "forge")),
            findings,
            st.session_state.get(notes_key, ""),
        )
        right.download_button(
            "Download report", report, file_name=f"{item.rsid}-variant-atlas.md", mime="text/markdown",
            on_click="ignore", key=f"report-{item.rsid}", width="stretch",
        )


def md_escape(text: str) -> str:
    return re.sub(r"([\\`*_{}\[\]()#+\-.!|<>$~:])", r"\\\1", text)


if __name__ == "__main__":
    main()
