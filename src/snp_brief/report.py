"""Research notebook and the full report export (Markdown).

The notebook is a free-text note per SNP, saved to <user cache>/<rsid>/notebook.md
so it survives a restart. The report gathers every module's output with sources.
"""

from __future__ import annotations

from pathlib import Path

from snp_brief.providers import tier_label


def load_notes(user_cache: Path, rsid: str) -> str:
    path = Path(user_cache) / rsid / "notebook.md"
    try:
        return path.read_text(encoding="utf-8") if path.is_file() else ""
    except OSError:
        return ""


def save_notes(user_cache: Path, rsid: str, text: str) -> bool:
    path = Path(user_cache) / rsid / "notebook.md"
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return True
    except OSError:
        return False


def build_report(item, summary, ripple_result=None, forge_result=None, findings=None, notes: str = "") -> str:
    out = [f"# {item.rsid}: Variant Atlas report", ""]
    out += ["## Summary", "", summary.text, "", f"_Written by {summary.writer}; every sentence passed the fact check._", ""]
    out += ["## Evidence brief", "", item.text, ""]
    out += ["## Facts used (with provider rank)", ""]
    out += [f"- [{f.id}] **{f.source}** ({tier_label(f.source)}): {f.text}" for f in summary.facts]
    out.append("")
    if findings:
        out += ["## Findings from paper abstracts (LLM-extracted, quote-verified)", ""]
        for f in findings:
            out.append(f'- PMID:{f["pmid"]}: {f["claim"]}  \n  > "{f["quote"]}"')
        out.append("")
    if ripple_result is not None:
        out += ["## Molecular Ripple", "", f"Variant context: {ripple_result.context}", ""]
        out += [f"- What to examine: {p}" for p in ripple_result.priorities]
        out += [f"- *{link.label}*: {link.text}" + (f" ({link.source})" if link.source else "") for link in ripple_result.links]
        out.append("")
    if forge_result is not None:
        tri = forge_result.triage
        out += ["## Hypothesis Forge", "", f"**Category:** {tri.category}", "", f"**Why:** {tri.reason}", "", f"**Open problem:** {tri.open_problem}", ""]
        if tri.exploratory:
            out += ["> Exploratory: AI-generated from genomic context, not from evidence about this variant.", ""]
        out += [f"- Next: {step}" for step in tri.next_steps]
        for number, hyp in enumerate(forge_result.hypotheses, 1):
            out += ["", f"### Hypothesis {number}: {hyp.title} ({hyp.consistency}, {hyp.feasibility} feasibility)", "", hyp.mechanism, ""]
            out += [f"- Support: {c['rendered']}" for c in hyp.support]
            out += [f"- Missing evidence: {g}" for g in hyp.gaps]
            out += [f"- AI background (unverified): {b}" for b in getattr(hyp, "background", [])]
            out += [f"- Would be supported by: {hyp.discriminating_observation}", f"- Experiment: {hyp.experiment}"]
        if forge_result.ranking:
            out += ["", "**Experiments, ranked**", ""]
            out += [f"{r['rank']}. {r['experiment']} (tests {r['tests']})" for r in forge_result.ranking]
        out.append("")
    out += ["## Sources", ""]
    for key, call in item.calls.items():
        out.append(f"- {key}: {call.origin}" + (f", retrieved {call.retrieved}" if call.retrieved else ""))
    if notes.strip():
        out += ["", "## Notebook", "", notes.strip()]
    out += ["", "_This report does not diagnose and is not medical advice._", ""]
    return "\n".join(out)
