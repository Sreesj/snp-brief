"""Streamlit page: run a list of SNPs, then investigate one in three modules.

Variant Observatory is the default module and shows the sourced evidence and
summary. Molecular Ripple and Hypothesis Forge run only when someone presses
Run in their tab; results are kept per SNP in the session and, for the Forge
LLM stage, in cache/<rsid>/forge.json.
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

from snp_brief.demo.view import (
    BLANK,
    FIELDS,
    SOURCED,
    UNRESOLVED,
    PopulationTable,
    briefs_text,
    field_status,
    gene_label,
    paper_links,
    parse_rsid_text,
    source_rows,
    statements,
    summary_row,
    unresolved_populations,
)
from snp_brief.forge import CONFLICTING, EXPLORATORY, INSUFFICIENT, REPORTED, TRAIT_ASSOCIATION, UNCERTAIN, WELL, forge
from snp_brief.providers import SPEED_BUDGETS, tier_label
from snp_brief.llm import LLMClient
from snp_brief.pipeline import Gathered, default_user_cache, load_variant
from snp_brief.graph import evidence_dot
from snp_brief.literature import Finding, extract_findings
from snp_brief.report import build_report, load_notes, save_notes
from snp_brief.ripple import OBSERVED, PREDICTED, UNKNOWN, ripple_cached
from snp_brief.summary import summarize

USER_CACHE = default_user_cache()
BADGE_COLOR = {SOURCED: "green", UNRESOLVED: "orange", BLANK: "gray"}
STATUS_CSS = {
    SOURCED: "background-color: #e6f4ea; color: #137333",
    UNRESOLVED: "background-color: #fef3e0; color: #a35200; font-weight: 600",
    BLANK: "color: #80868b",
}
NAMED_CSS = STATUS_CSS[UNRESOLVED]
SMALL_CSS = "color: #9aa0a6"
FALLBACK_CSS = "background-color: #f3eefc; color: #5b3a99; font-weight: 600"
MODULES = ("Variant Observatory", "Molecular Ripple", "Hypothesis Forge")
CATEGORY_COLOR = {
    WELL: "green",
    CONFLICTING: "orange",
    UNCERTAIN: "blue",
    TRAIT_ASSOCIATION: "blue",
    REPORTED: "violet",
    EXPLORATORY: "gray",
    INSUFFICIENT: "gray",
}
LINK_COLOR = {OBSERVED: "green", PREDICTED: "violet", UNKNOWN: "gray"}
PAGE_CSS = """
<style>
[data-testid="stMainBlockContainer"] { padding-top: 3.5rem; }
.snp-brief {
  font-size: 1rem;
  line-height: 1.65;
  padding: 0.85rem 1.1rem;
  border: 1px solid rgba(49, 51, 63, 0.18);
  border-radius: 0.5rem;
  background: #f8f9fb;
}
</style>
"""


@dataclass
class Run:
    number: int
    live: bool
    results: list[Gathered]
    rejected: list[str]
    failed: list[str]


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
    st.title("Variant Atlas")
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
    st.set_page_config(page_title="Variant Atlas", layout="wide")
    load_secrets()
    if not password_ok():
        return
    st.html(PAGE_CSS)
    root = repo_root()
    if "rsid_text" not in st.session_state:
        fixture = root / "fixture" / "snps.txt"
        st.session_state.rsid_text = fixture.read_text(encoding="utf-8") if fixture.is_file() else ""
    offline, clicked = sidebar()
    st.title("Variant Atlas")
    st.caption(
        "Sourced SNP briefs: biological function, allele frequency, and disease relevance, then mechanism and "
        "open questions. Every claim names its source. The tool does not diagnose and does not invent a number."
    )
    if clicked:
        run_briefs(root, offline)
    search_bar(root, offline)
    current: Run | None = st.session_state.get("last_run")
    if current is None:
        st.info("The list in the sidebar is prefilled with fixture/snps.txt. Press Run briefs.")
        return
    show_results(current)


def repo_root() -> Path:
    for candidate in (Path.cwd(), Path(__file__).resolve().parents[3]):
        if (candidate / "fixture" / "snps.txt").is_file():
            return candidate
    return Path.cwd()


def sidebar() -> tuple[bool, bool]:
    with st.sidebar:
        st.header("SNP list")
        st.text_area("rsIDs, one per line", key="rsid_text", height=300)
        st.file_uploader("Or upload a .txt", type=["txt"], key="upload", on_change=take_upload)
        st.radio(
            "Speed vs completeness",
            list(SPEED_BUDGETS),
            index=1,
            key="speed",
            horizontal=True,
            help="All providers are queried at once. Fast shows what has answered within about 8 s, balanced 20 s, "
            "thorough 45 s. Late providers keep loading in the background; press Check again to add them.",
        )
        offline = st.toggle("Offline only", key="offline")
        if offline:
            st.caption("Reads saved data only and makes no network call. SNPs that were never fetched come back blank.")
        else:
            st.caption(
                "Saved data first; anything missing is fetched live from dbSNP, ClinVar, gnomAD, Europe PMC, and "
                "UniProt and saved for next time. Use Refresh evidence on a SNP to re-fetch it."
            )
        if LLMClient.from_env() is not None:
            st.toggle("LLM summary", key="use_llm", value=True)
            st.caption("Sentences that fail the fact check are dropped; if too many fail, the template is used.")
        clicked = st.button("Run briefs", key="run_briefs", type="primary", width="stretch")
    return offline, clicked


def take_upload() -> None:
    upload = st.session_state.get("upload")
    if upload is not None:
        st.session_state.rsid_text = upload.getvalue().decode("utf-8", errors="replace")


def fetch(root: Path, rsid: str, *, offline: bool, refresh: bool = False) -> Gathered:
    budget = SPEED_BUDGETS.get(st.session_state.get("speed", "balanced"), 20.0)
    return load_variant(
        rsid, repo_cache=root / "cache", user_cache=USER_CACHE, offline=offline, refresh=refresh, budget=budget
    )


def reload_variant(item: Gathered) -> None:
    """Pick up providers that finished in the background (cache first, nothing re-fetched)."""
    current: Run | None = st.session_state.get("last_run")
    fresh = fetch(repo_root(), item.rsid, offline=bool(st.session_state.get("offline")))
    if current is not None:
        current.results = [fresh if old.rsid == item.rsid else old for old in current.results]
    modules = st.session_state.setdefault("modules", {})
    for key in [key for key in modules if key[0] == item.rsid]:
        modules.pop(key)


def run_briefs(root: Path, offline: bool) -> None:
    ids, rejected = parse_rsid_text(st.session_state.rsid_text)
    if not ids:
        st.session_state.last_run = None
        st.error("No rsIDs to run. Enter ids such as rs1815739, one per line.")
        return
    progress = st.progress(0.0)
    results: list[Gathered] = []
    failed: list[str] = []
    for index, rsid in enumerate(ids):
        progress.progress(index / len(ids), text=f"Loading {rsid} ({index + 1} of {len(ids)}); new SNPs take 10-30 s")
        try:
            results.append(fetch(root, rsid, offline=offline))
        except Exception as exc:
            failed.append(f"{rsid} failed: {exc}")
    progress.empty()
    st.session_state.run_count = st.session_state.get("run_count", 0) + 1
    st.session_state.last_run = Run(st.session_state.run_count, not offline, results, rejected, failed)


def search_bar(root: Path, offline: bool) -> None:
    with st.form("search", border=False):
        left, right = st.columns([5, 1], vertical_alignment="bottom")
        query = left.text_input("Search one rsID", placeholder="rs7903146", key="search_text")
        submitted = right.form_submit_button("Search", width="stretch")
    if not submitted:
        return
    ids, _ = parse_rsid_text(query or "")
    if not ids:
        st.error("Enter one rsID, such as rs7903146.")
        return
    rsid = ids[0]
    current: Run | None = st.session_state.get("last_run")
    if current is None:
        st.session_state.run_count = st.session_state.get("run_count", 0) + 1
        current = Run(st.session_state.run_count, not offline, [], [], [])
        st.session_state.last_run = current
    if all(item.rsid != rsid for item in current.results):
        try:
            with st.spinner(f"Loading {rsid}; a SNP that was never fetched takes 10-30 s…"):
                current.results.append(fetch(root, rsid, offline=offline))
        except Exception as exc:
            st.error(f"{rsid} failed: {exc}")
            return
    st.session_state.active = rsid


def refresh_variant(item: Gathered) -> list[str]:
    """Re-fetch one SNP from the live sources and drop its module results.

    If any source fails, the previous evidence is kept, so a network problem
    never replaces a good record with a blank one. Returns the failed sources.
    """
    current: Run | None = st.session_state.get("last_run")
    fresh = fetch(repo_root(), item.rsid, offline=False, refresh=True)
    failed = [name for name, call in fresh.calls.items() if call.origin == "failed"]
    if failed:
        return failed
    if current is not None:
        current.results = [fresh if old.rsid == item.rsid else old for old in current.results]
    modules = st.session_state.setdefault("modules", {})
    for key in [key for key in modules if key[0] == item.rsid]:
        modules.pop(key)
    return []


def origin_line(current: Run) -> str:
    counts = {"cache": 0, "live": 0, "missing": 0}
    for item in current.results:
        for call in item.calls.values():
            if call.origin == "cache":
                counts["cache"] += 1
            elif call.origin == "live":
                counts["live"] += 1
            else:
                counts["missing"] += 1
    bits = [f"{counts['cache']} source answers from saved data", f"{counts['live']} fetched live"]
    if counts["missing"]:
        bits.append(f"{counts['missing']} missing or failed")
    return ", ".join(bits)


def show_results(current: Run) -> None:
    st.markdown(f"**{len(current.results)} SNPs**: {origin_line(current)}.")
    if current.rejected:
        st.warning("Skipped, not rsIDs: " + ", ".join(current.rejected))
    for problem in current.failed:
        st.error(problem)
    if not current.results:
        return
    table_side, detail_side = st.columns([4, 8], gap="large")
    with table_side:
        frame = pd.DataFrame([summary_row(item) for item in current.results])
        event = st.dataframe(
            frame.style.map(lambda status: STATUS_CSS.get(status, ""), subset=list(FIELDS)),
            hide_index=True,
            column_config={
                name: st.column_config.TextColumn(name, help="sourced, unresolved, or blank") for name in FIELDS
            },
            on_select="rerun",
            selection_mode="single-row",
            key=f"briefs-{current.number}",
        )
        st.download_button(
            "Download all briefs (.txt)",
            briefs_text(current.results),
            file_name="snp-briefs.txt",
            mime="text/plain",
            on_click="ignore",
            width="stretch",
        )
        st.caption(
            "sourced: a source reports it (for disease, the best-reviewed ClinVar call). unresolved: sources at the same review level disagree, so the brief keeps both. "
            "blank: not reported, or ClinVar found no clinical significance."
        )
    with detail_side:
        selected = event.selection.rows
        by_id = {item.rsid: item for item in current.results}
        if selected:
            st.session_state.active = current.results[selected[0]].rsid
        active = by_id.get(st.session_state.get("active", ""))
        if active is not None:
            render_detail(active)
        else:
            st.info("Select a row, or search an rsID, to open it in the three modules.")


def render_detail(item: Gathered) -> None:
    gene = gene_label(item)
    head, button = st.columns([5, 2], vertical_alignment="bottom")
    head.subheader(f"{item.rsid} · {gene}" if gene else item.rsid)
    if button.button("Refresh evidence", key=f"refresh-{item.rsid}", help="Re-fetch every source live for this SNP"):
        with st.spinner(f"Re-fetching {item.rsid}…"):
            failed = refresh_variant(item)
        if failed:
            st.error(f"Refresh failed for {', '.join(failed)}; the previous evidence is kept.")
        else:
            st.rerun()
    st.caption(f"Evidence last retrieved: {last_retrieved(item)}")
    if item.pending:
        names = ", ".join(item.pending)
        note, again = st.columns([5, 2], vertical_alignment="center")
        note.info(f"Still loading: {names}. Shown now: the providers that answered in time.")
        if again.button("Check again", key=f"again-{item.rsid}"):
            reload_variant(item)
            st.rerun()
    results = st.session_state.setdefault("modules", {})
    failures = st.session_state.setdefault("module_errors", {})
    with st.container(horizontal=True):
        st.badge("Observatory: ready", color="green")
        for name, key in (("Ripple", "ripple"), ("Forge", "forge")):
            if (item.rsid, key) in failures:
                st.badge(f"{name}: failed", color="red")
            elif (item.rsid, key) in results:
                st.badge(f"{name}: ready", color="green")
            else:
                st.badge(f"{name}: not run", color="gray")
    observatory, ripple_tab, forge_tab = st.tabs(list(MODULES))
    with observatory:
        render_observatory(item)
    with ripple_tab:
        render_ripple(item)
    with forge_tab:
        render_forge(item)
    render_notebook(item)


def run_module(item: Gathered, key: str, work) -> None:
    """Run a module, recording success or failure for the status badges."""
    results = st.session_state.setdefault("modules", {})
    failures = st.session_state.setdefault("module_errors", {})
    try:
        results[(item.rsid, key)] = work()
        failures.pop((item.rsid, key), None)
    except Exception as exc:  # shown as a red "failed" badge with the message
        failures[(item.rsid, key)] = str(exc)


def render_notebook(item: Gathered) -> None:
    with st.expander("Research notebook and full report"):
        notes_key = f"notes-{item.rsid}"
        if notes_key not in st.session_state:
            st.session_state[notes_key] = load_notes(USER_CACHE, item.rsid)
        st.text_area("Your notes on this SNP (saved on this computer)", key=notes_key, height=150)
        left, right = st.columns(2)
        if left.button("Save notes", key=f"save-notes-{item.rsid}"):
            ok = save_notes(USER_CACHE, item.rsid, st.session_state[notes_key])
            st.toast("Notes saved." if ok else "Could not save notes here; they stay in this session.")
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
            "Download full report (.md)", report, file_name=f"{item.rsid}-variant-atlas.md", mime="text/markdown",
            on_click="ignore", key=f"report-{item.rsid}",
        )


def last_retrieved(item: Gathered) -> str:
    stamps = sorted(call.retrieved for call in item.calls.values() if call.retrieved)
    if not stamps:
        return "not retrieved"
    return stamps[-1].replace("T", " ")[:16] + " UTC"


def module_cache(item: Gathered) -> Path:
    return item.cache_dir or repo_root() / "cache"


def render_observatory(item: Gathered) -> None:
    found = statements(item.sentences)
    with st.container(horizontal=True):
        for name in FIELDS:
            status = field_status(found, name)
            st.badge(f"{name}: {status}", color=BADGE_COLOR[status])
    mode = "llm" if st.session_state.get("use_llm") else "template"
    summary = cached_summary(item, mode)
    st.markdown("**Summary**")
    st.html(f'<div class="snp-brief">{html.escape(summary.text)}</div>')
    checked = f"{len(summary.sentences)} sentences passed the fact check"
    if summary.rejected:
        checked += f"; {len(summary.rejected)} rejected"
    st.caption(f"Written by {summary.writer} from {len(summary.facts)} retrieved facts. {checked}.")
    with st.expander("Facts the summary may use"):
        st.code(
            "\n".join(f"[{fact.id}] ({fact.source}; {tier_label(fact.source)}) {fact.text}" for fact in summary.facts),
            language=None,
        )
        for sentence, reason in summary.rejected:
            st.caption(f"Rejected: {sentence} ({reason})")
    st.markdown("**Evidence brief**")
    st.html(f'<div class="snp-brief">{html.escape(item.text)}</div>')
    tables = unresolved_populations(item)
    if tables:
        st.markdown("**gnomAD populations** (frequency varies by ancestry)")
        for table in tables:
            render_populations(table)
        st.caption(
            "Allele count / allele number, as gnomAD returned them. Highlighted: the groups the brief names. "
            "Grey: fewer than 1000 alleles, so not compared."
        )
    render_findings(item)
    st.markdown("**Papers**")
    links = paper_links(item)
    if links:
        st.markdown(
            "\n".join(
                f"{number}. [{md_escape(link.title)}]({link.url}) · {md_escape(link.citation)}"
                for number, link in enumerate(links, 1)
            )
        )
    elif "europepmc" in item.calls and item.calls["europepmc"].origin in {"cache", "live"}:
        st.caption("Europe PMC lists no primary research article that names this rsID.")
    else:
        st.caption("Europe PMC did not answer.")
    st.markdown("**Sources**")
    sources = pd.DataFrame(source_rows(item))
    st.table(
        sources.style.map(lambda name: FALLBACK_CSS if "fallback" in name else "", subset=["Source"]),
        hide_index=True,
        border="horizontal",
    )


def cached_summary(item: Gathered, mode: str):
    store = st.session_state.setdefault("summaries", {})
    key = (item.rsid, mode, item.text)
    if key not in store:
        store[key] = summarize(item, mode, cache_dir=module_cache(item))
    return store[key]


def render_findings(item: Gathered) -> None:
    store = st.session_state.setdefault("findings", {})
    client = LLMClient.from_env()
    st.markdown("**Findings from paper abstracts**")
    found = store.get(item.rsid)
    if found is None and item.findings:
        found = {"findings": [vars(f) for f in item.findings], "rejected": [], "writer": "cached"}
    if found is None:
        st.caption("The LLM reads up to five abstracts and extracts findings about this SNP; each must quote the abstract verbatim.")
        if st.button("Extract findings", key=f"findings-{item.rsid}", disabled=client is None):
            with st.spinner("Reading abstracts…"):
                try:
                    found = extract_findings(
                        item, client.for_forge(), module_cache(item), offline=bool(st.session_state.get("offline"))
                    )
                except Exception as exc:
                    found = {"findings": [], "rejected": [], "writer": "", "note": f"Failed: {exc}"}
            store[item.rsid] = found
            item.findings = [Finding(**f) for f in found.get("findings", [])]
            st.session_state.get("summaries", {}).clear()
            st.session_state.get("modules", {}).pop((item.rsid, "forge"), None)
            st.rerun()
        if client is None:
            st.caption("Connect an LLM (NEBIUS_API_KEY) to extract findings.")
        return
    if found.get("note"):
        st.caption(found["note"])
    for f in found.get("findings", []):
        st.markdown(f'- {f["claim"]} (PMID:{f["pmid"]})  \n  > "{f["quote"]}"')
    if found.get("findings"):
        st.caption("LLM-extracted, verified: each quote appears verbatim in the abstract and every number in the claim appears in the quote. These now feed the summary and the Forge.")
    for r in found.get("rejected", []):
        st.caption(f"Rejected: {r['claim']} ({r['reason']})")


def render_ripple(item: Gathered) -> None:
    results = st.session_state.setdefault("modules", {})
    key = (item.rsid, "ripple")
    if key not in results:
        st.caption("Molecular context and what a mechanistic analysis should examine for this kind of variant.")
        error = st.session_state.get("module_errors", {}).get(key)
        if error:
            st.error(f"Molecular Ripple failed: {error}")
        if st.button("Run Molecular Ripple", key=f"run-ripple-{item.rsid}", type="primary"):
            run_module(item, "ripple", lambda: ripple_cached(item, module_cache(item)))
            st.rerun()
        return
    found = results[key]
    st.markdown(f"**Variant context:** {found.context}" + (f" · {found.gene}" if found.gene else ""))
    st.markdown("**What to examine for this context**")
    st.markdown("\n".join(f"- {line}" for line in found.priorities))
    st.markdown("**Molecular links**")
    for link in found.links:
        tag, text = st.columns([1, 6])
        with tag:
            st.badge(link.label, color=LINK_COLOR[link.label])
        text.markdown(link.text + (f" · *{link.source}*" if link.source else ""))
    st.caption(
        "observed: reported by a source. predicted: a computational score. unknown: not answered by the evidence. "
        "Predictions are never drawn as established links."
    )
    st.markdown("**Evidence graph**")
    forge_result = results.get((item.rsid, "forge"))
    try:
        st.graphviz_chart(evidence_dot(item, forge_result), width="stretch")
    except Exception as exc:
        st.caption(f"Graph unavailable: {exc}")
    st.caption("Solid green: observed. Dashed violet: predicted. Dotted grey: Hypothesis Forge hypotheses (run the Forge to add them).")


def render_forge(item: Gathered) -> None:
    results = st.session_state.setdefault("modules", {})
    key = (item.rsid, "forge")
    client = LLMClient.from_env()
    if key not in results:
        st.caption(
            "Stage 1 sorts the variant by what is not yet understood. For conflicting or uncertain variants, "
            "stages 2 and 3 propose competing hypotheses and experiments to tell them apart."
        )
        if client is None:
            st.caption("No LLM connected (set NEBIUS_API_KEY). Stage 1 runs; cached hypotheses are shown if present.")
        if client is not None and forge_runs_left() == 0:
            st.caption("This session has used its AI hypothesis runs; stage 1 and cached results still work.")
            client = None
        error = st.session_state.get("module_errors", {}).get(key)
        if error:
            st.error(f"Hypothesis Forge failed: {error}")
        if st.button("Run Hypothesis Forge", key=f"run-forge-{item.rsid}", type="primary"):
            with st.spinner("Forging hypotheses…" if client else "Running stage 1…"):
                run_module(item, "forge", lambda: forge(item, client=client, cache_dir=module_cache(item)))
            found = results.get(key)
            if found is not None and found.writer.startswith("llm:") and "(cached)" not in found.writer:
                st.session_state.forge_llm_runs = st.session_state.get("forge_llm_runs", 0) + 1
            st.rerun()
        return
    found = results[key]
    tri = found.triage
    st.badge(tri.category, color=CATEGORY_COLOR[tri.category])
    st.markdown(f"**Why:** {tri.reason}")
    st.markdown(f"**Open problem:** {tri.open_problem}")
    st.markdown("**Next steps**")
    st.markdown("\n".join(f"- {line}" for line in tri.next_steps))
    if tri.exploratory:
        st.warning(
            "Exploratory: there is no evidence about this specific variant. The hypotheses below are AI-generated from "
            "its genomic context and gene-level information, and may be wrong. Use them as starting points only."
        )
    if found.note:
        st.info(found.note)
    if found.open_questions:
        st.markdown("**Open questions**")
        st.markdown("\n".join(f"- {q}" for q in found.open_questions))
    if found.hypotheses:
        st.markdown("**Competing hypotheses**")
        columns = st.columns(len(found.hypotheses))
        for column, hyp in zip(columns, found.hypotheses):
            with column, st.container(border=True):
                st.markdown(f"**{hyp.title}**")
                st.caption("Hypothesis, not an established mechanism.")
                st.markdown(hyp.mechanism)
                st.markdown("*Evidence for*")
                st.markdown("\n".join(f"- {claim['rendered']}" for claim in hyp.support))
                st.markdown("*Evidence against*")
                st.markdown("\n".join(f"- {claim['rendered']}" for claim in hyp.against) or "- none in the retrieved facts")
                if hyp.gaps:
                    st.markdown("*Missing evidence*")
                    st.markdown("\n".join(f"- {gap}" for gap in hyp.gaps))
                if hyp.background:
                    st.markdown("*AI background knowledge (unverified)*")
                    st.markdown("\n".join(f"- {line}" for line in hyp.background))
                st.markdown(f"*Would be supported by:* {hyp.discriminating_observation}")
                st.markdown(f"*Experiment:* {hyp.experiment} ({hyp.feasibility} feasibility)")
                st.caption(f"Consistency with the facts: {hyp.consistency}")
    if found.ranking:
        st.markdown("**Experiments, ranked**")
        st.dataframe(
            pd.DataFrame(found.ranking)[["rank", "experiment", "tests", "consistency", "separates", "support", "feasibility", "open_gaps"]],
            hide_index=True,
        )
        st.caption(
            "Ranked by consistency with the facts, then how many other hypotheses the result would tell apart, then net "
            "support (evidence for minus against), then feasibility, then evidence gaps addressed. No model-assigned probabilities."
        )
    if found.hypotheses or found.writer != "rules":
        st.caption(f"Hypotheses written by {found.writer}. Supporting statements passed the same fact check as the summary.")
    if found.hypotheses and client is not None:
        with st.expander("Steer the Forge: reject hypotheses or add your constraints"):
            rejected = [
                hyp.title for number, hyp in enumerate(found.hypotheses, 1)
                if st.checkbox(f"Reject H{number}: {hyp.title}", key=f"reject-{item.rsid}-{number}")
            ]
            notes = st.text_area(
                "Constraints or resources (e.g. 'no animal work', 'we have intestinal organoids', 'computational only')",
                key=f"steer-{item.rsid}",
            )
            if st.button("Re-forge with my input", key=f"reforge-{item.rsid}") and (rejected or notes.strip()):
                if forge_runs_left() == 0:
                    st.warning("This session has used its AI hypothesis runs.")
                else:
                    with st.spinner("Re-forging…"):
                        run_module(
                            item, "forge",
                            lambda: forge(item, client=client, cache_dir=module_cache(item), steer={"rejected": rejected, "notes": notes}),
                        )
                    st.session_state.forge_llm_runs = st.session_state.get("forge_llm_runs", 0) + 1
                    st.rerun()


def render_populations(table: PopulationTable) -> None:
    st.markdown(f"{table.allele} · gnomAD {table.variant_id}")
    st.caption(" · ".join(table.totals))
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

    st.dataframe(frame.style.apply(highlight, axis=None), hide_index=True)


def md_escape(text: str) -> str:
    return re.sub(r"([\\`*_{}\[\]()#+\-.!|<>$~:])", r"\\\1", text)


if __name__ == "__main__":
    main()
