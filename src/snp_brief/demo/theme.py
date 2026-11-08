"""Look and feel of the Variant Atlas page: CSS, and small HTML builders.

Pure functions (no Streamlit), so the markup is testable without a browser.
"""

from __future__ import annotations

import html
import re

from snp_brief.reconcile import BENIGN, OTHER, PATHOGENIC, UNCERTAIN, side

FONTS = (
    '<link rel="preconnect" href="https://fonts.googleapis.com">'
    '<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>'
    '<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&'
    'family=JetBrains+Mono:wght@500&display=swap" rel="stylesheet">'
)

CSS = """
<style>
:root {
  --va-bg: #f6f7f9;
  --va-surface: #ffffff;
  --va-border: #e6e8ec;
  --va-text: #0e1726;
  --va-muted: #5f6b7a;
  --va-faint: #8a94a3;
  --va-brand: #1f4fd8;
  --va-brand-soft: #eef3ff;
  --va-radius: 14px;
  --va-shadow: 0 1px 2px rgba(16, 24, 40, .04), 0 1px 3px rgba(16, 24, 40, .06);
}
html, body, [class*="css"], .stApp, button, input, textarea {
  font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif !important;
}
.stApp { background: var(--va-bg); }
[data-testid="stMainBlockContainer"] { padding-top: 1.6rem; max-width: 1280px; }
#MainMenu, footer, [data-testid="stDecoration"], [data-testid="stDeployButton"], .stDeployButton { display: none !important; }
[data-testid="stHeader"] { background: transparent; }
h1, h2, h3, h4 { letter-spacing: -0.01em; color: var(--va-text); }
p, li { color: var(--va-text); }

/* top bar */
.va-top { display: flex; align-items: center; gap: .7rem; margin-bottom: .2rem; }
.va-logo { width: 34px; height: 34px; border-radius: 10px;
  background: linear-gradient(135deg, #1f4fd8 0%, #18a3a3 100%);
  display: grid; place-items: center; box-shadow: var(--va-shadow); }
.va-brand { font-size: 1.25rem; font-weight: 700; color: var(--va-text); }
.va-tagline { color: var(--va-muted); font-size: .92rem; margin: 0 0 1rem 0; }

/* cards */
.va-card { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: var(--va-radius);
  box-shadow: var(--va-shadow); padding: 1.15rem 1.3rem; margin-bottom: .9rem; }
.va-head { display: flex; flex-wrap: wrap; align-items: baseline; gap: .6rem 1rem; }
.va-rsid { font-family: "JetBrains Mono", ui-monospace, monospace; font-size: 1.55rem; font-weight: 500; color: var(--va-text); }
.va-gene { font-size: 1.05rem; font-weight: 600; color: var(--va-brand); }
.va-loc { color: var(--va-muted); font-size: .88rem; font-family: "JetBrains Mono", ui-monospace, monospace; }
.va-note { color: var(--va-muted); font-size: .85rem; margin-top: .35rem; }
.va-chips { display: flex; flex-wrap: wrap; gap: .4rem; margin-top: .75rem; }
.va-chip { display: inline-flex; align-items: center; gap: .35rem; border-radius: 999px; padding: .22rem .65rem;
  font-size: .78rem; font-weight: 600; border: 1px solid transparent; white-space: nowrap; }
.va-chip::before { content: ""; width: 6px; height: 6px; border-radius: 50%; background: currentColor; opacity: .85; }
.va-chip.pathogenic { color: #b42318; background: #fef3f2; border-color: #fecdca; }
.va-chip.benign { color: #067647; background: #ecfdf3; border-color: #abefc6; }
.va-chip.uncertain { color: #b54708; background: #fffaeb; border-color: #fedf89; }
.va-chip.association { color: #175cd3; background: #eff8ff; border-color: #b2ddff; }
.va-chip.neutral { color: #475467; background: #f2f4f7; border-color: #e4e7ec; }
.va-chip.brand { color: var(--va-brand); background: var(--va-brand-soft); border-color: #d1dcff; }

/* key facts */
.va-facts { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: .8rem; margin-bottom: .9rem; }
@media (max-width: 900px) { .va-facts { grid-template-columns: 1fr; } }
.va-fact { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: 12px; padding: .9rem 1rem;
  box-shadow: var(--va-shadow); border-top: 3px solid var(--va-border); }
.va-fact.pathogenic { border-top-color: #f04438; } .va-fact.benign { border-top-color: #17b26a; }
.va-fact.uncertain { border-top-color: #f79009; } .va-fact.association { border-top-color: #2e90fa; }
.va-fact.brand { border-top-color: var(--va-brand); }
.va-fact-label { text-transform: uppercase; letter-spacing: .06em; font-size: .7rem; font-weight: 600; color: var(--va-faint); }
.va-fact-value { font-size: 1.12rem; font-weight: 600; color: var(--va-text); margin-top: .25rem; line-height: 1.3; }
.va-fact-sub { font-size: .82rem; color: var(--va-muted); margin-top: .25rem; line-height: 1.4; }

/* summary */
.va-summary { font-size: 1.02rem; line-height: 1.75; color: var(--va-text); }
.va-src { display: inline-block; font-size: .68rem; font-weight: 600; color: var(--va-muted); background: #f2f4f7;
  border: 1px solid #e4e7ec; border-radius: 6px; padding: 0 .35rem; margin-left: .2rem; vertical-align: 1px; }
.va-section { font-size: .78rem; font-weight: 600; text-transform: uppercase; letter-spacing: .06em; color: var(--va-faint);
  margin: .2rem 0 .55rem 0; }
.va-callout { border-radius: 12px; padding: .8rem 1rem; font-size: .9rem; margin: .4rem 0 .9rem 0;
  background: #fffaeb; border: 1px solid #fedf89; color: #93370d; }
.va-callout.info { background: var(--va-brand-soft); border-color: #d1dcff; color: #1d3fa8; }

/* ripple links */
.va-link { display: flex; gap: .7rem; align-items: flex-start; padding: .55rem 0; border-bottom: 1px solid #f0f2f5; }
.va-link:last-child { border-bottom: none; }
.va-tag { font-size: .68rem; font-weight: 700; text-transform: uppercase; letter-spacing: .05em; border-radius: 6px;
  padding: .15rem .45rem; min-width: 74px; text-align: center; }
.va-tag.observed { color: #067647; background: #ecfdf3; } .va-tag.predicted { color: #6941c6; background: #f4f3ff; }
.va-tag.unknown { color: #475467; background: #f2f4f7; }
.va-link-text { font-size: .9rem; color: var(--va-text); line-height: 1.45; }
.va-link-src { color: var(--va-faint); font-size: .78rem; }

/* loader */
.va-loader { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: var(--va-radius);
  padding: 1.1rem 1.3rem; box-shadow: var(--va-shadow); margin: .4rem 0 1rem 0; }
.va-loader-row { display: flex; align-items: center; gap: .8rem; }
.va-helix { display: flex; gap: 4px; align-items: center; height: 22px; }
.va-helix span { width: 5px; height: 100%; border-radius: 3px; background: linear-gradient(#1f4fd8, #18a3a3);
  animation: va-wave 1.1s ease-in-out infinite; }
.va-helix span:nth-child(2) { animation-delay: .1s } .va-helix span:nth-child(3) { animation-delay: .2s }
.va-helix span:nth-child(4) { animation-delay: .3s } .va-helix span:nth-child(5) { animation-delay: .4s }
@keyframes va-wave { 0%, 100% { transform: scaleY(.35); opacity: .5 } 50% { transform: scaleY(1); opacity: 1 } }
.va-loader-title { font-weight: 600; color: var(--va-text); }
.va-loader-sub { color: var(--va-muted); font-size: .85rem; }
.va-bar { height: 6px; background: #eef1f5; border-radius: 999px; overflow: hidden; margin-top: .85rem; }
.va-bar > div { height: 100%; background: linear-gradient(90deg, #1f4fd8, #18a3a3); border-radius: 999px; transition: width .4s ease; }
.va-skel { border-radius: 8px; height: 14px; margin: .55rem 0;
  background: linear-gradient(90deg, #eef1f5 25%, #f7f9fb 37%, #eef1f5 63%); background-size: 400% 100%;
  animation: va-shimmer 1.3s ease infinite; }
@keyframes va-shimmer { 0% { background-position: 100% 50% } 100% { background-position: 0 50% } }

/* streamlit widgets */
[data-testid="stTabs"] button[role="tab"] { font-weight: 600; font-size: .95rem; padding: .55rem .2rem; }
[data-testid="stTabs"] [data-baseweb="tab-highlight"] { background-color: var(--va-brand); }
[data-testid="stExpander"] { background: var(--va-surface); border-radius: 12px; border: 1px solid var(--va-border); box-shadow: var(--va-shadow); }
[data-testid="stExpander"] summary { font-weight: 600; }
[data-testid="stVerticalBlockBorderWrapper"] { border-radius: 12px; }
.stButton > button, .stDownloadButton > button, [data-testid="stFormSubmitButton"] > button {
  border-radius: 10px; font-weight: 600; }
[data-testid="stTextInput"] input { border-radius: 10px; font-size: 1rem; padding: .7rem .9rem; }
[data-testid="stSidebar"] { background: var(--va-surface); border-right: 1px solid var(--va-border); }
.va-list-title { font-size: .78rem; font-weight: 600; text-transform: uppercase; letter-spacing: .06em; color: var(--va-faint);
  margin: .2rem 0 .4rem 0; }
div[class*="st-key-pick-"] button { justify-content: flex-start; text-align: left; font-weight: 500;
  font-family: "JetBrains Mono", ui-monospace, monospace !important; font-size: .85rem; padding: .45rem .7rem; }
div[class*="st-key-try-"] button { border-radius: 999px; padding: .1rem .7rem; font-size: .8rem;
  font-family: "JetBrains Mono", ui-monospace, monospace !important; }
</style>
"""

LOGO = (
    '<div class="va-logo"><svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="white" stroke-width="2.2" '
    'stroke-linecap="round"><path d="M7 3c0 6 10 6 10 12s-10 6-10 6"/><path d="M17 3c0 6-10 6-10 12s10 6 10 6" opacity=".55"/>'
    '<path d="M8.5 7h7M8.5 17h7" opacity=".8"/></svg></div>'
)


def e(text) -> str:
    return html.escape(str(text if text is not None else ""))


def top_bar() -> str:
    return (
        f'<div class="va-top">{LOGO}<div class="va-brand">Variant Atlas</div></div>'
        '<p class="va-tagline">Sourced evidence, molecular context and open questions for any SNP.</p>'
    )


def chip(text: str, kind: str = "neutral") -> str:
    return f'<span class="va-chip {kind}">{e(text)}</span>'


def classification_kind(headline: str, unresolved: bool = False) -> str:
    if unresolved:
        return "uncertain"
    if not headline:
        return "neutral"
    found = side(headline)
    return {PATHOGENIC: "pathogenic", BENIGN: "benign", UNCERTAIN: "uncertain", OTHER: "association"}[found]


def clinical(plan) -> tuple[str, str, str]:
    """(value, detail, kind) for the Clinical card and chip."""
    if plan.unresolved_groups:
        return "Conflicting", "ClinVar submissions at the same review level disagree", "uncertain"
    if plan.headline:
        kind = classification_kind(plan.headline)
        cond = ", ".join(plan.conditions[:2])
        label = plan.headline[:1].upper() + plan.headline[1:]
        return label, (cond if cond else "ClinVar"), kind
    if plan.gwas_traits:
        return "Trait association", ", ".join(plan.gwas_traits[:2]), "association"
    if plan.no_significance:
        return "No clinical significance", "Not classified in ClinVar", "neutral"
    return "Not reported", "", "neutral"


def function_card(plan) -> tuple[str, str]:
    consequence = (plan.consequence or "").split(",")[0].strip().replace("_", " ")
    if consequence:
        value = consequence[:1].upper() + consequence[1:]
    elif plan.location:
        value = plan.location[:1].upper() + plan.location[1:]
    else:
        value = "Not reported"
    sub = " · ".join(x for x in (plan.gene, plan.change) if x)
    return value, sub


def frequency_card(plan) -> tuple[str, str]:
    if not plan.overall_pct:
        return "Not reported", ""
    sub = plan.overall_set.replace("gnomAD ", "gnomAD ")
    if plan.high and plan.low:
        sub = f"{plan.high[1]} {plan.high[0]} → {plan.low[1]} {plan.low[0]}"
    return f"{plan.overall_pct} of alleles", sub


def header_card(item, plan, gene: str, updated: str) -> str:
    db = item.db
    loc = ""
    if db.chrom and db.pos:
        alleles = ", ".join(f"{r}>{a}" for r, a in db.alleles[:3])
        loc = f"{db.assembly or ''} · chr{db.chrom}:{db.pos:,}" + (f" · {alleles}" if alleles else "")
    value, _, kind = clinical(plan)
    chips = []
    function_value, _ = function_card(plan)
    if function_value != "Not reported":
        chips.append(chip(function_value, "brand"))
    chips.append(chip(value, kind))
    if plan.overall_pct:
        chips.append(chip(f"{plan.overall_pct} in gnomAD", "neutral"))
    if plan.gwas_traits:
        chips.append(chip(f"{plan.gwas_count} trait associations", "association"))
    merged = ""
    if getattr(item, "merged_from", ""):
        merged = f'<div class="va-note">{e(item.merged_from)} was merged into {e(item.rsid)} in dbSNP.</div>'
    return (
        '<div class="va-card"><div class="va-head">'
        f'<span class="va-rsid">{e(item.rsid)}</span>'
        + (f'<span class="va-gene">{e(gene)}</span>' if gene else "")
        + (f'<span class="va-loc">{e(loc.strip(" ·"))}</span>' if loc else "")
        + f'</div>{merged}<div class="va-chips">{"".join(chips)}</div>'
        + (f'<div class="va-note">Updated {e(updated)}</div>' if updated else "")
        + "</div>"
    )


def facts_grid(plan) -> str:
    f_value, f_sub = function_card(plan)
    q_value, q_sub = frequency_card(plan)
    c_value, c_sub, c_kind = clinical(plan)
    cards = [
        ("Function", f_value, f_sub, "brand"),
        ("Frequency", q_value, q_sub, "brand" if plan.overall_pct else "neutral"),
        ("Clinical", c_value, c_sub, c_kind),
    ]
    return '<div class="va-facts">' + "".join(
        f'<div class="va-fact {kind}"><div class="va-fact-label">{e(label)}</div>'
        f'<div class="va-fact-value">{e(value)}</div><div class="va-fact-sub">{e(sub)}</div></div>'
        for label, value, sub, kind in cards
    ) + "</div>"


_TAG = re.compile(r"\[(F\d+(?:\s*,\s*F\d+)*)\]")


def summary_html(sentences: list[str], facts) -> str:
    """The summary with each sentence's sources as small grey tags."""
    by_id = {fact.id: fact.source for fact in facts}
    parts = []
    for sentence in sentences:
        sources: list[str] = []
        for tag in _TAG.findall(sentence):
            for fact_id in tag.split(","):
                name = by_id.get(fact_id.strip())
                if name and name not in sources:
                    sources.append(name)
        text = _TAG.sub("", sentence).strip()
        text = re.sub(r"\s+([.,;])", r"\1", text)
        pills = "".join(f'<span class="va-src">{e(name)}</span>' for name in sources)
        parts.append(f"{e(text)}{pills}")
    return f'<div class="va-card va-summary">{" ".join(parts)}</div>'


def loader(title: str, sub: str = "", done: int = 0, total: int = 0, skeleton: bool = True) -> str:
    bar = ""
    if total:
        pct = max(4, int(100 * done / total))
        bar = f'<div class="va-bar"><div style="width:{pct}%"></div></div>'
    skel = ""
    if skeleton:
        skel = "".join(f'<div class="va-skel" style="width:{w}%"></div>' for w in (92, 78, 85, 60))
    return (
        '<div class="va-loader"><div class="va-loader-row"><div class="va-helix">'
        + "<span></span>" * 5
        + f'</div><div><div class="va-loader-title">{e(title)}</div>'
        + (f'<div class="va-loader-sub">{e(sub)}</div>' if sub else "")
        + f"</div></div>{bar}{skel}</div>"
    )


def link_rows(links) -> str:
    rows = []
    for link in links:
        src = f' <span class="va-link-src">· {e(link.source)}</span>' if link.source else ""
        rows.append(
            f'<div class="va-link"><span class="va-tag {e(link.label)}">{e(link.label)}</span>'
            f'<div class="va-link-text">{e(link.text)}{src}</div></div>'
        )
    return '<div class="va-card">' + "".join(rows) + "</div>"


def callout(text: str, kind: str = "") -> str:
    return f'<div class="va-callout {kind}">{e(text)}</div>'


def section(title: str) -> str:
    return f'<div class="va-section">{e(title)}</div>'
