"""Streamlit page: paste rsIDs, run the briefs, and open one SNP at a time."""

from __future__ import annotations

import html
import re
import tempfile
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
from snp_brief.pipeline import Gathered, gather

LIVE_CACHE = Path(tempfile.gettempdir()) / "snp-brief-live-cache"
BADGE_COLOR = {SOURCED: "green", UNRESOLVED: "orange", BLANK: "gray"}
STATUS_CSS = {
    SOURCED: "background-color: #e6f4ea; color: #137333",
    UNRESOLVED: "background-color: #fef3e0; color: #a35200; font-weight: 600",
    BLANK: "color: #80868b",
}
NAMED_CSS = STATUS_CSS[UNRESOLVED]
SMALL_CSS = "color: #9aa0a6"
FALLBACK_CSS = "background-color: #f3eefc; color: #5b3a99; font-weight: 600"
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
    st.set_page_config(page_title="SNP brief", layout="wide")
    st.html(PAGE_CSS)
    root = repo_root()
    if "rsid_text" not in st.session_state:
        fixture = root / "fixture" / "snps.txt"
        st.session_state.rsid_text = fixture.read_text(encoding="utf-8") if fixture.is_file() else ""
    live, clicked = sidebar()
    st.title("SNP brief")
    st.caption(
        "One plain-text paragraph per SNP: biological function, allele frequency, and disease relevance. "
        "Every claim names its source. The tool does not diagnose and does not invent a number."
    )
    if clicked:
        run_briefs(root, live)
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
        live = st.toggle("Live fetching", key="live")
        if live:
            st.caption(
                "Calls dbSNP, ClinVar, gnomAD, Europe PMC, and UniProt now. "
                f"Responses are saved to {LIVE_CACHE}; cache/ is not changed."
            )
        else:
            st.caption("Cached: reads the checked-in cache/ and makes no network call.")
        clicked = st.button("Run briefs", key="run_briefs", type="primary", width="stretch")
    return live, clicked


def take_upload() -> None:
    upload = st.session_state.get("upload")
    if upload is not None:
        st.session_state.rsid_text = upload.getvalue().decode("utf-8", errors="replace")


def run_briefs(root: Path, live: bool) -> None:
    ids, rejected = parse_rsid_text(st.session_state.rsid_text)
    if not ids:
        st.session_state.last_run = None
        st.error("No rsIDs to run. Enter ids such as rs1815739, one per line.")
        return
    cache_dir = LIVE_CACHE if live else root / "cache"
    verb = "Fetching" if live else "Reading"
    progress = st.progress(0.0)
    results: list[Gathered] = []
    failed: list[str] = []
    for index, rsid in enumerate(ids):
        progress.progress(index / len(ids), text=f"{verb} {rsid} ({index + 1} of {len(ids)})")
        try:
            results.append(gather(rsid, cache_dir=cache_dir, refresh=live, offline=not live))
        except Exception as exc:
            failed.append(f"{rsid} failed: {exc}")
    progress.empty()
    st.session_state.run_count = st.session_state.get("run_count", 0) + 1
    st.session_state.last_run = Run(st.session_state.run_count, live, results, rejected, failed)


def show_results(current: Run) -> None:
    origin = "live sources" if current.live else "the checked-in cache, offline"
    st.markdown(f"**{len(current.results)} SNPs** from {origin}.")
    if current.rejected:
        st.warning("Skipped, not rsIDs: " + ", ".join(current.rejected))
    for problem in current.failed:
        st.error(problem)
    if not current.results:
        return
    table_side, detail_side = st.columns([5, 7], gap="large")
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
            "sourced: a source reports it. unresolved: sources disagree, so the brief keeps both. "
            "blank: not reported, or ClinVar found no clinical significance."
        )
    with detail_side:
        selected = event.selection.rows
        if selected:
            render_detail(current.results[selected[0]])
        else:
            st.info("Select a row to read its brief, gnomAD populations, papers, and sources.")


def render_detail(item: Gathered) -> None:
    found = statements(item.sentences)
    gene = gene_label(item)
    st.subheader(f"{item.rsid} · {gene}" if gene else item.rsid)
    with st.container(horizontal=True):
        for name in FIELDS:
            status = field_status(found, name)
            st.badge(f"{name}: {status}", color=BADGE_COLOR[status])
    st.markdown("**Brief**")
    st.html(f'<div class="snp-brief">{html.escape(item.text)}</div>')
    tables = unresolved_populations(item)
    if tables:
        st.markdown("**gnomAD populations** (frequency is unresolved)")
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
