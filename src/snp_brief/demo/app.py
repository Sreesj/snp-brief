"""Streamlit page: run a list of SNPs, then investigate one in three modules.

Variant Observatory is the default module and shows the sourced evidence and
summary. Molecular Ripple and Hypothesis Forge run only when someone presses
Run in their tab; results are kept per SNP in the session and, for the Forge
LLM stage, in cache/<rsid>/forge.json.
"""

from __future__ import annotations

import html
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
from snp_brief.forge import CONFLICTING, INSUFFICIENT, REPORTED, UNCERTAIN, WELL, forge
from snp_brief.llm import LLMClient
from snp_brief.pipeline import Gathered, default_user_cache, load_variant
from snp_brief.ripple import OBSERVED, PREDICTED, UNKNOWN, ripple
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
CATEGORY_COLOR = {WELL: "green", CONFLICTING: "orange", UNCERTAIN: "blue", REPORTED: "violet", INSUFFICIENT: "gray"}
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


def main() -> None:
    st.set_page_config(page_title="Variant Atlas", layout="wide")
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
    return load_variant(rsid, repo_cache=root / "cache", user_cache=USER_CACHE, offline=offline, refresh=refresh)


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
    results = st.session_state.setdefault("modules", {})
    with st.container(horizontal=True):
        st.badge("Observatory: ready", color="green")
        for name, key in (("Ripple", "ripple"), ("Forge", "forge")):
            done = (item.rsid, key) in results
            st.badge(f"{name}: {'ready' if done else 'not run'}", color="green" if done else "gray")
    observatory, ripple_tab, forge_tab = st.tabs(list(MODULES))
    with observatory:
        render_observatory(item)
    with ripple_tab:
        render_ripple(item)
    with forge_tab:
        render_forge(item)


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
        st.code("\n".join(f"[{fact.id}] ({fact.source}) {fact.text}" for fact in summary.facts), language=None)
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


def render_ripple(item: Gathered) -> None:
    results = st.session_state.setdefault("modules", {})
    key = (item.rsid, "ripple")
    if key not in results:
        st.caption("Molecular context and what a mechanistic analysis should examine for this kind of variant.")
        if st.button("Run Molecular Ripple", key=f"run-ripple-{item.rsid}", type="primary"):
            results[key] = ripple(item)
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
        if st.button("Run Hypothesis Forge", key=f"run-forge-{item.rsid}", type="primary"):
            with st.spinner("Forging hypotheses…" if client else "Running stage 1…"):
                results[key] = forge(item, client=client, cache_dir=module_cache(item))
            st.rerun()
        return
    found = results[key]
    tri = found.triage
    st.badge(tri.category, color=CATEGORY_COLOR[tri.category])
    st.markdown(f"**Why:** {tri.reason}")
    st.markdown(f"**Open problem:** {tri.open_problem}")
    st.markdown("**Next steps**")
    st.markdown("\n".join(f"- {line}" for line in tri.next_steps))
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
                st.markdown("*Supported by*")
                st.markdown("\n".join(f"- {claim['rendered']}" for claim in hyp.support))
                if hyp.gaps:
                    st.markdown("*Missing evidence*")
                    st.markdown("\n".join(f"- {gap}" for gap in hyp.gaps))
                st.markdown(f"*Would be supported by:* {hyp.discriminating_observation}")
                st.markdown(f"*Experiment:* {hyp.experiment} ({hyp.feasibility} feasibility)")
                st.caption(f"Consistency with the facts: {hyp.consistency}")
    if found.ranking:
        st.markdown("**Experiments, ranked**")
        st.dataframe(
            pd.DataFrame(found.ranking)[["rank", "experiment", "tests", "consistency", "support", "feasibility", "open_gaps"]],
            hide_index=True,
        )
        st.caption(
            "Ranked by consistency with the facts, then number of fact-checked supporting statements, then feasibility, "
            "then evidence gaps addressed. No model-assigned probabilities."
        )
    if found.hypotheses or found.writer != "rules":
        st.caption(f"Hypotheses written by {found.writer}. Supporting statements passed the same fact check as the summary.")


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
