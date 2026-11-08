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
  --va-bg: #08090a;
  --va-surface: #101113;
  --va-surface-2: #16171a;
  --va-border: rgba(255, 255, 255, .08);
  --va-border-strong: rgba(255, 255, 255, .14);
  --va-text: #f7f8f8;
  --va-muted: #b4b9c1;
  --va-faint: #8d929b;
  --va-accent: #7c86ff;
  --va-accent-2: #4cc9b0;
  --va-radius: 12px;
}
html, body, [class*="css"], .stApp, button, input, textarea, select {
  font-family: "Inter", -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif !important;
  -webkit-font-smoothing: antialiased;
}
html, body { background: var(--va-bg) !important; }
.stApp, [data-testid="stAppViewContainer"], [data-testid="stMain"] { background: transparent !important; color: var(--va-text); }
[data-testid="stSidebarContent"] { position: relative; z-index: 1; }
.element-container:has(> iframe[height="0"]), [data-testid="stElementContainer"]:has(iframe[height="0"]),
.st-key-va-motion { position: absolute !important; width: 1px !important; height: 1px !important; overflow: hidden !important; opacity: 0; pointer-events: none; }
[data-testid="stMainBlockContainer"] { padding-top: 1.4rem; max-width: 1200px; position: relative; z-index: 1; }
#MainMenu, footer, [data-testid="stDecoration"], [data-testid="stDeployButton"], .stDeployButton,
[data-testid="stStatusWidget"] { display: none !important; }
[data-testid="stHeader"] { background: transparent; }
h1, h2, h3, h4, p, li, label, span { color: var(--va-text); }
a { color: var(--va-accent) !important; text-decoration: none; }
a:hover { text-decoration: underline; }
[data-testid="stMarkdownContainer"] p, [data-testid="stMarkdownContainer"] li { color: #d0d3d8; line-height: 1.6; }
[data-testid="stMarkdownContainer"] strong { color: var(--va-text); font-weight: 600; }
[data-testid="stMarkdownContainer"] blockquote { border-left: 2px solid var(--va-accent); color: var(--va-muted); }
[data-testid="stCaptionContainer"], .stCaption { color: var(--va-muted) !important; }

/* top bar and hero */
.va-top { display: flex; align-items: center; gap: .6rem; margin: .2rem 0 1.2rem 0; }
.va-logo { width: 26px; height: 26px; border-radius: 7px; display: grid; place-items: center;
  background: #f4f5f6 url("data:image/svg+xml,%3Csvg%20xmlns%3D%22http%3A//www.w3.org/2000/svg%22%20viewBox%3D%220%200%2024%2024%22%20fill%3D%22none%22%20stroke%3D%22%230a0a0b%22%20stroke-width%3D%222.3%22%20stroke-linecap%3D%22round%22%3E%3Cpath%20d%3D%22M7%202.5c0%206.5%2010%206.5%2010%2013S7%2021.5%207%2021.5%22/%3E%3Cpath%20d%3D%22M17%202.5c0%206.5-10%206.5-10%2013s10%206%2010%206%22/%3E%3Cpath%20d%3D%22M9%206.5h6M9.6%2012h4.8M9%2017.5h6%22%20stroke-width%3D%221.8%22/%3E%3C/svg%3E") center / 16px 16px no-repeat; box-shadow: 0 0 0 1px rgba(255,255,255,.08), 0 4px 14px rgba(0,0,0,.4); }
.va-brand { font-size: 1rem; font-weight: 600; color: var(--va-text); letter-spacing: -0.01em; }
.va-hero { text-align: center; padding: 9vh 0 1.6rem 0; }
.va-hero h1 { font-size: clamp(2.2rem, 5vw, 3.6rem); font-weight: 600; letter-spacing: -0.035em; line-height: 1.05; margin: 0;
  background: linear-gradient(180deg, #ffffff 30%, #9ea3ab 100%); -webkit-background-clip: text; background-clip: text; color: transparent; }
.va-hero p { color: var(--va-muted); font-size: 1.08rem; margin: 1rem auto 0 auto; max-width: 560px; line-height: 1.55; }
.va-pill { display: inline-flex; align-items: center; gap: .45rem; border: 1px solid var(--va-border-strong); border-radius: 999px;
  padding: .28rem .8rem; font-size: .8rem; color: var(--va-muted); margin-bottom: 1.4rem; background: rgba(255,255,255,.02); }
.va-pill i { width: 6px; height: 6px; border-radius: 50%; background: var(--va-accent-2); display: inline-block; }
.va-features { display: flex; justify-content: center; gap: .7rem; flex-wrap: wrap; margin-top: 2.2rem; }
.va-feature { color: #d0d3d8; font-size: .86rem; display: flex; align-items: center; gap: .45rem;
  background: rgba(16,17,19,.85); border: 1px solid var(--va-border); border-radius: 999px; padding: .35rem .8rem; }
.va-feature b { color: var(--va-accent); font-weight: 600; }

/* cards */
.va-card { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: var(--va-radius);
  padding: 1.15rem 1.3rem; margin-bottom: .9rem; }
.va-head { display: flex; flex-wrap: wrap; align-items: baseline; gap: .5rem 1rem; }
.va-rsid { font-family: "JetBrains Mono", ui-monospace, monospace; font-size: 1.6rem; font-weight: 500; color: var(--va-text); letter-spacing: -0.02em; }
.va-gene { font-size: 1.02rem; font-weight: 600; color: var(--va-accent); }
.va-loc { color: var(--va-muted); font-size: .84rem; font-family: "JetBrains Mono", ui-monospace, monospace; }
.va-note { color: var(--va-faint); font-size: .8rem; margin-top: .5rem; }
.va-chips { display: flex; flex-wrap: wrap; gap: .4rem; margin-top: .8rem; }
.va-chip { display: inline-flex; align-items: center; gap: .4rem; border-radius: 999px; padding: .2rem .62rem;
  font-size: .76rem; font-weight: 500; border: 1px solid var(--va-border-strong); color: #c9ccd1; background: rgba(255,255,255,.03); white-space: nowrap; }
.va-chip::before { content: ""; width: 6px; height: 6px; border-radius: 50%; background: #6b7078; }
.va-chip.pathogenic::before { background: #f26d6d; } .va-chip.benign::before { background: #4cc38a; }
.va-chip.uncertain::before { background: #f0b73c; } .va-chip.association::before { background: #5aa2ff; }
.va-chip.brand::before { background: var(--va-accent); } .va-chip.neutral::before { background: #6b7078; }

/* key facts */
.va-facts { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: .75rem; margin-bottom: .9rem; }
@media (max-width: 900px) { .va-facts { grid-template-columns: 1fr; } }
.va-fact { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: 12px; padding: .95rem 1.05rem; }
.va-fact-label { display: flex; align-items: center; gap: .45rem; font-size: .74rem; font-weight: 500; color: var(--va-muted); }
.va-fact-label::before { content: ""; width: 7px; height: 7px; border-radius: 2px; background: #6b7078; }
.va-fact.pathogenic .va-fact-label::before { background: #f26d6d; } .va-fact.benign .va-fact-label::before { background: #4cc38a; }
.va-fact.uncertain .va-fact-label::before { background: #f0b73c; } .va-fact.association .va-fact-label::before { background: #5aa2ff; }
.va-fact.brand .va-fact-label::before { background: var(--va-accent); }
.va-fact-value { font-size: 1.08rem; font-weight: 600; color: var(--va-text); margin-top: .4rem; line-height: 1.3; letter-spacing: -0.01em; }
.va-fact-sub { font-size: .8rem; color: var(--va-muted); margin-top: .3rem; line-height: 1.45; }

/* summary */
.va-summary { font-size: 1rem; line-height: 1.75; color: #dcdfe3; }
.va-src { display: inline-block; font-size: .66rem; font-weight: 500; color: var(--va-muted); border: 1px solid var(--va-border-strong);
  border-radius: 5px; padding: 0 .32rem; margin-left: .25rem; vertical-align: 1px; }
.va-section { font-size: .74rem; font-weight: 500; color: var(--va-muted); margin: .6rem 0 .5rem 0; letter-spacing: .01em; }
.va-callout { border-radius: 10px; padding: .75rem .95rem; font-size: .88rem; margin: .3rem 0 .9rem 0;
  background: rgba(240, 183, 60, .06); border: 1px solid rgba(240, 183, 60, .22); color: #e9c879; }
.va-callout.info { background: rgba(124, 134, 255, .06); border-color: rgba(124, 134, 255, .22); color: #b9beff; }

/* ripple links */
.va-link { display: flex; gap: .8rem; align-items: flex-start; padding: .6rem 0; border-bottom: 1px solid var(--va-border); }
.va-link:last-child { border-bottom: none; }
.va-tag { font-size: .66rem; font-weight: 600; text-transform: uppercase; letter-spacing: .06em; border-radius: 5px;
  padding: .14rem .45rem; min-width: 76px; text-align: center; border: 1px solid var(--va-border-strong); }
.va-tag.observed { color: #4cc38a; } .va-tag.predicted { color: #a99cff; } .va-tag.unknown { color: var(--va-muted); }
.va-link-text { font-size: .9rem; color: #d0d3d8; line-height: 1.45; }
.va-link-src { color: var(--va-faint); font-size: .78rem; }

/* loader */
.va-loader { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: var(--va-radius);
  padding: 1.1rem 1.3rem; margin: .4rem 0 1rem 0; }
.va-loader-row { display: flex; align-items: center; gap: .85rem; }
.va-helix { display: flex; gap: 3px; align-items: center; height: 18px; }
.va-helix span { width: 3px; height: 100%; border-radius: 2px; background: var(--va-accent); animation: va-wave 1.1s ease-in-out infinite; }
.va-helix span:nth-child(2) { animation-delay: .1s; background: #8f8bff } .va-helix span:nth-child(3) { animation-delay: .2s; background: #6fa7e8 }
.va-helix span:nth-child(4) { animation-delay: .3s; background: #5bbfcf } .va-helix span:nth-child(5) { animation-delay: .4s; background: var(--va-accent-2) }
@keyframes va-wave { 0%, 100% { transform: scaleY(.3); opacity: .45 } 50% { transform: scaleY(1); opacity: 1 } }
.va-loader-title { font-weight: 500; color: var(--va-text); font-size: .95rem; }
.va-loader-sub { color: var(--va-muted); font-size: .82rem; }
.va-bar { height: 3px; background: rgba(255,255,255,.06); border-radius: 999px; overflow: hidden; margin-top: .9rem; }
.va-bar > div { height: 100%; background: linear-gradient(90deg, var(--va-accent), var(--va-accent-2)); transition: width .4s ease; }
.va-skel { border-radius: 6px; height: 12px; margin: .6rem 0;
  background: linear-gradient(90deg, rgba(255,255,255,.04) 25%, rgba(255,255,255,.09) 37%, rgba(255,255,255,.04) 63%);
  background-size: 400% 100%; animation: va-shimmer 1.4s ease infinite; }
@keyframes va-shimmer { 0% { background-position: 100% 50% } 100% { background-position: 0 50% } }

/* streamlit widgets on dark */
[data-testid="stTabs"] [role="tablist"] { gap: 1.4rem; border-bottom: 1px solid var(--va-border); }
[data-testid="stTabs"] button[role="tab"] { font-weight: 500; font-size: .92rem; color: var(--va-muted); padding: .6rem 0; background: transparent; }
[data-testid="stTabs"] button[role="tab"][aria-selected="true"] { color: var(--va-text); }
[data-testid="stTabs"] [data-baseweb="tab-highlight"] { background-color: var(--va-accent); height: 2px; }
[data-testid="stTabs"] [data-baseweb="tab-border"] { display: none; }
[data-testid="stExpander"] { background: var(--va-surface); border: 1px solid var(--va-border) !important; border-radius: 12px; }
[data-testid="stExpander"] summary { font-weight: 500; color: var(--va-text); }
[data-testid="stExpander"] summary:hover { color: var(--va-accent); }
[data-testid="stVerticalBlockBorderWrapper"] { border-color: var(--va-border) !important; border-radius: 12px; background: var(--va-surface); }
.stButton > button, .stDownloadButton > button, [data-testid="stFormSubmitButton"] > button {
  border-radius: 8px; font-weight: 500; border: 1px solid var(--va-border-strong); background: var(--va-surface-2); color: var(--va-text); }
.stButton > button:hover, .stDownloadButton > button:hover { border-color: var(--va-accent); color: var(--va-text); }
.stButton > button[kind="primary"], [data-testid="stFormSubmitButton"] > button[kind="primaryFormSubmit"],
[data-testid="stFormSubmitButton"] > button { background: var(--va-text); color: #08090a; border-color: var(--va-text); }
.stButton > button[kind="primary"]:hover, [data-testid="stFormSubmitButton"] > button:hover { background: #dfe1e4; color: #08090a; }
[data-testid="stTextInput"] input, [data-testid="stTextArea"] textarea {
  background: var(--va-surface) !important; color: var(--va-text) !important; border-radius: 10px; }
[data-testid="stTextInput"] input { font-size: 1rem; padding: .78rem 1rem; }
[data-testid="stTextInput"] > div > div, [data-testid="stTextArea"] > div > div { border-color: var(--va-border-strong) !important; background: var(--va-surface) !important; }
[data-testid="stTextInput"] input::placeholder { color: var(--va-faint); }
[data-testid="stForm"] { border: none; padding: 0; }
[data-testid="stSidebar"] { background: #0b0c0d; border-right: 1px solid var(--va-border); }
[data-testid="stTable"] table, [data-testid="stTable"] th, [data-testid="stTable"] td { color: #eef0f2 !important; border-color: var(--va-border) !important; background: transparent !important; font-size: .84rem; }
[data-testid="stTable"] th { color: var(--va-muted) !important; font-weight: 500; }
/* sidebar guide */
.va-guide { margin: .2rem 0 1.2rem 0; }
.va-guide-item { display: flex; gap: .7rem; padding: .55rem 0; border-bottom: 1px solid var(--va-border); }
.va-guide-item:last-child { border-bottom: none; }
.va-guide-icon { width: 22px; height: 22px; border-radius: 6px; flex: none; display: grid; place-items: center;
  border: 1px solid var(--va-border-strong); background: rgba(255,255,255,.02); }
.va-guide-icon i { width: 7px; height: 7px; border-radius: 2px; display: block; }
.va-guide-name { font-size: .86rem; font-weight: 600; color: var(--va-text); }
.va-guide-text { font-size: .78rem; color: var(--va-muted); line-height: 1.45; margin-top: .1rem; }
.va-list-title { font-size: .74rem; font-weight: 500; color: var(--va-muted); margin: .2rem 0 .5rem 0; }
div[class*="st-key-pick-"] button { justify-content: flex-start; text-align: left; font-weight: 500; background: transparent;
  border-color: transparent; color: var(--va-muted); font-family: "JetBrains Mono", ui-monospace, monospace !important; font-size: .84rem; padding: .42rem .65rem; }
div[class*="st-key-pick-"] button:hover { background: var(--va-surface-2); color: var(--va-text); border-color: var(--va-border); }
div[class*="st-key-pick-"] button[kind="primary"] { background: var(--va-surface-2); color: var(--va-text); border-color: var(--va-border-strong); }
/* labels on light (primary) buttons must be dark; the global p/span colour would make them white-on-white */
.stButton > button[kind="primary"] *, [data-testid="stFormSubmitButton"] > button * , [data-testid="stDownloadButton"] > button[kind="primary"] * { color: #08090a !important; }
div[class*="st-key-pick-"] button * { color: inherit !important; }
.va-waiting { display: flex; align-items: center; gap: .55rem; color: var(--va-muted); font-size: .82rem; margin: .2rem 0 .6rem 0; }
.va-waiting i { width: 7px; height: 7px; border-radius: 50%; background: var(--va-accent-2); animation: va-pulse 1.4s ease-in-out infinite; }
@keyframes va-pulse { 0%,100% { opacity: .35; transform: scale(.8); } 50% { opacity: 1; transform: scale(1.1); } }

/* module icons (match the deck): database stack, concentric rings, flask */
.va-ico { display: inline-block; flex: none; background: center / contain no-repeat; }
.va-ico-observatory { background-image: url("data:image/svg+xml,%3Csvg%20xmlns%3D%22http%3A//www.w3.org/2000/svg%22%20viewBox%3D%220%200%2024%2024%22%20fill%3D%22none%22%20stroke%3D%22%238b7cff%22%20stroke-width%3D%221.8%22%20stroke-linecap%3D%22round%22%20stroke-linejoin%3D%22round%22%3E%3Cellipse%20cx%3D%2212%22%20cy%3D%225%22%20rx%3D%227.5%22%20ry%3D%222.6%22/%3E%3Cpath%20d%3D%22M4.5%205v14c0%201.4%203.4%202.6%207.5%202.6s7.5-1.2%207.5-2.6V5%22/%3E%3Cpath%20d%3D%22M4.5%209.7c0%201.4%203.4%202.6%207.5%202.6s7.5-1.2%207.5-2.6%22/%3E%3Cpath%20d%3D%22M4.5%2014.4c0%201.4%203.4%202.6%207.5%202.6s7.5-1.2%207.5-2.6%22/%3E%3C/svg%3E"); }
.va-ico-ripple { background-image: url("data:image/svg+xml,%3Csvg%20xmlns%3D%22http%3A//www.w3.org/2000/svg%22%20viewBox%3D%220%200%2024%2024%22%20fill%3D%22none%22%20stroke%3D%22%236aa8ff%22%20stroke-width%3D%221.8%22%20stroke-linecap%3D%22round%22%20stroke-linejoin%3D%22round%22%3E%3Ccircle%20cx%3D%2212%22%20cy%3D%2212%22%20r%3D%229.5%22/%3E%3Ccircle%20cx%3D%2212%22%20cy%3D%2212%22%20r%3D%225%22/%3E%3C/svg%3E"); }
.va-ico-forge { background-image: url("data:image/svg+xml,%3Csvg%20xmlns%3D%22http%3A//www.w3.org/2000/svg%22%20viewBox%3D%220%200%2024%2024%22%20fill%3D%22none%22%20stroke%3D%22%239be39b%22%20stroke-width%3D%221.8%22%20stroke-linecap%3D%22round%22%20stroke-linejoin%3D%22round%22%3E%3Cpath%20d%3D%22M9%202.5h6%22/%3E%3Cpath%20d%3D%22M10%202.5v6.8L4.3%2019.2a1.6%201.6%200%200%200%201.4%202.3h12.6a1.6%201.6%200%200%200%201.4-2.3L14%209.3V2.5%22/%3E%3Cpath%20d%3D%22M6.6%2015h10.8%22/%3E%3C/svg%3E"); }
.va-guide-icon .va-ico { width: 16px; height: 16px; }
.va-tile b .va-ico { width: 20px; height: 20px; }
.va-feature .va-ico { width: 16px; height: 16px; }
.va-feature b.observatory { color: #8b7cff; } .va-feature b.ripple { color: #6aa8ff; } .va-feature b.forge { color: #9be39b; }
/* about section and scroll reveal */
.va-reveal { opacity: 0; transform: translateY(26px); transition: opacity .9s cubic-bezier(.2,.7,.2,1), transform .9s cubic-bezier(.2,.7,.2,1); }
.va-reveal.va-in { opacity: 1; transform: none; }
.va-about { margin: 14vh auto 6vh auto; max-width: 980px; }
.va-about-block { margin: 0 0 13vh 0; }
.va-eyebrow { font-family: "JetBrains Mono", ui-monospace, monospace; font-size: .74rem; letter-spacing: .14em; text-transform: uppercase;
  color: var(--va-accent-2); margin-bottom: .8rem; }
.va-about h2 { font-size: clamp(1.7rem, 3.4vw, 2.6rem); font-weight: 600; letter-spacing: -0.03em; line-height: 1.1; margin: 0 0 1rem 0; color: var(--va-text); }
.va-about .va-lede { font-size: 1.05rem; line-height: 1.65; color: #b4b8bf; max-width: 680px; margin: 0 0 2rem 0; }
.va-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: .9rem; }
.va-grid.va-two { grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); }
.va-hero h1 { filter: drop-shadow(0 2px 18px rgba(8,9,10,.95)); }
.va-hero p, .va-about h2, .va-about .va-lede { text-shadow: 0 2px 24px rgba(8,9,10,.9), 0 0 2px rgba(8,9,10,.8); }
.va-tile { background: rgba(16,17,19,.94); backdrop-filter: blur(10px); -webkit-backdrop-filter: blur(10px);
  border: 1px solid var(--va-border); border-radius: 12px; padding: 1.1rem 1.15rem; }
.va-tile b { display: flex; align-items: center; gap: .5rem; color: var(--va-text); font-weight: 600; font-size: .95rem; margin-bottom: .4rem; }
.va-tile b i { width: 7px; height: 7px; border-radius: 2px; display: inline-block; }
.va-tile span { color: #c9cdd3; font-size: .88rem; line-height: 1.55; display: block; }
.va-steps { counter-reset: step; }
.va-steps .va-tile b::before { counter-increment: step; content: "0" counter(step); font-family: "JetBrains Mono", monospace;
  font-size: .72rem; color: var(--va-accent); font-weight: 500; }
.va-note { border-left: 2px solid var(--va-border-strong); padding: .2rem 0 .2rem 1rem; color: var(--va-muted); font-size: .88rem; line-height: 1.6; max-width: 720px; }
.va-scroll-cue { text-align: center; color: var(--va-muted); font-size: .76rem; letter-spacing: .12em; text-transform: uppercase; margin-top: 7vh; }
.va-scroll-cue::after { content: ""; display: block; width: 1px; height: 42px; margin: .7rem auto 0; background: linear-gradient(var(--va-muted), transparent);
  animation: va-cue 2.2s ease-in-out infinite; }
@keyframes va-cue { 0%,100% { opacity: .3; transform: scaleY(.6); transform-origin: top; } 50% { opacity: 1; transform: scaleY(1); } }
</style>
"""

# The DNA mark is a CSS background image: st.html strips inline <svg> elements.
LOGO = '<div class="va-logo" role="img" aria-label="Variant Atlas"></div>'


def e(text) -> str:
    return html.escape(str(text if text is not None else ""))


def top_bar() -> str:
    return f'<div class="va-top">{LOGO}<div class="va-brand">Variant Atlas</div></div>'


def hero() -> str:
    return (
        '<div class="va-hero">'
        '<div class="va-pill"><i></i>18 genomic sources, reconciled</div>'
        "<h1>Understand any variant.</h1>"
        "<p>Clinical significance, population frequency, molecular context and the open questions worth testing, "
        "for any SNP, with every claim traced to its source.</p>"
        "</div>"
    )


MODULE_ICONS = {"Variant Observatory": "observatory", "Molecular Ripple": "ripple", "Hypothesis Forge": "forge",
                "Observatory": "observatory", "Ripple": "ripple", "Forge": "forge"}


def module_icon(name: str) -> str:
    key = MODULE_ICONS.get(name)
    return f'<span class="va-ico va-ico-{key}" role="img" aria-label="{e(name)} icon"></span>' if key else ""


GUIDE = (
    ("#8b7cff", "Variant Observatory", "What is known: clinical classification, population frequency, function and papers, "
     "reconciled across 18 sources with every claim sourced."),
    ("#f2555a", "Interpretation", "A draft of the ACMG evidence codes, how many people carry the variant, "
     "drug response from ClinPGx and cancer evidence from CIViC."),
    ("#6aa8ff", "Molecular Ripple", "What the variant touches: gene and protein effects, expression links, regulatory "
     "overlap and prediction scores, labelled observed, predicted or unknown."),
    ("#9be39b", "Hypothesis Forge", "What is not yet understood: competing explanations, evidence for and against, and "
     "the experiments that would tell them apart."),
    ("#5aa2ff", "Findings from papers", "Key results pulled from abstracts, each quoted word for word."),
    ("#4cc38a", "Several at once", "Search one rsID or many: separate them with commas or spaces. Keep notes and "
     "download a full report for any variant."),
)


def sidebar_guide() -> str:
    return '<div class="va-guide">' + "".join(
        f'<div class="va-guide-item"><div class="va-guide-icon">{module_icon(name) or f"<i style=background:{color}></i>"}</div>'
        f'<div><div class="va-guide-name">{e(name)}</div><div class="va-guide-text">{e(text)}</div></div></div>'
        for color, name, text in GUIDE
    ) + "</div>"


def features() -> str:
    items = (("Observatory", "evidence and summary"), ("Interpretation", "ACMG draft, carriers, drugs"), ("Ripple", "molecular context"), ("Forge", "competing hypotheses"))
    return '<div class="va-features">' + "".join(
        f'<div class="va-feature">{module_icon(name)}<b class="{MODULE_ICONS.get(name, "")}">{e(name)}</b>{e(text)}</div>'
        for name, text in items
    ) + "</div>"


ABOUT_WHAT = (
    ("#8b7cff", "Variant Observatory", "The evidence, reconciled. Clinical classification from ClinVar, allele frequencies "
     "across gnomAD populations, GWAS associations, gene function and the literature, merged into one plain-language "
     "summary where every sentence names its source."),
    ("#6aa8ff", "Molecular Ripple", "What the variant touches. Protein changes, expression effects in GTEx tissues, "
     "regulatory overlap, linked variants and gene prioritisation, each labelled observed, predicted or unknown."),
    ("#9be39b", "Hypothesis Forge", "Where the knowledge runs out. Competing explanations for what the variant does, "
     "the evidence for and against each, and the experiments that would tell them apart."),
)

ABOUT_WHO = (
    ("Clinical and lab geneticists", "Check a variant's classification, review status and conflicts before interpretation."),
    ("Genetic counsellors", "Get a clear, sourced summary of what is and is not known to prepare for a conversation."),
    ("Researchers", "Triage GWAS hits and candidate variants, see mechanisms and find the open question worth testing."),
    ("Students and bioinformaticians", "Learn how the major genomic databases describe the same variant, side by side."),
)

ABOUT_HOW = (
    ("Ask", "Enter one rsID or a whole list, separated by commas or spaces."),
    ("Gather", "Eighteen sources are queried in parallel, from ClinVar, gnomAD and ClinPGx to GTEx, Europe PMC and our own Borzoi predictions."),
    ("Reconcile", "Sources are ranked by reliability, so when they disagree the conflict is shown rather than hidden."),
    ("Explain", "A plain-text summary is written from the collected facts and checked so it says nothing the data does not."),
)


def about() -> str:
    def tiles(rows, dots: bool) -> str:
        out = []
        for row in rows:
            if dots:
                color, name, text = row
                head = (module_icon(name) or f'<i style="background:{color}"></i>') + e(name)
            else:
                name, text = row
                head = e(name)
            out.append(f'<div class="va-tile va-reveal"><b>{head}</b><span>{e(text)}</span></div>')
        return "".join(out)

    return (
        '<div class="va-scroll-cue">About Variant Atlas</div>'
        '<div class="va-about">'
        '<div class="va-about-block">'
        '<div class="va-eyebrow va-reveal">What it does</div>'
        '<h2 class="va-reveal">From an rsID to an answer you can trust.</h2>'
        '<p class="va-lede va-reveal">Information about a single variant is scattered across a dozen databases and '
        "thousands of papers, often in disagreement. Variant Atlas collects it, reconciles it and explains it in plain "
        "text: what the variant does, how common it is, which diseases and traits it is linked to, and how sure anyone "
        "really is.</p>"
        f'<div class="va-grid">{tiles(ABOUT_WHAT, True)}</div>'
        "</div>"
        '<div class="va-about-block">'
        '<div class="va-eyebrow va-reveal">Who it is for</div>'
        '<h2 class="va-reveal">Built for people who read variants for a living.</h2>'
        '<p class="va-lede va-reveal">Whether you are interpreting a clinical report, planning an experiment or '
        "learning genomics, the same question comes up: what do we actually know about this variant?</p>"
        f'<div class="va-grid va-two">{tiles(ABOUT_WHO, False)}</div>'
        "</div>"
        '<div class="va-about-block">'
        '<div class="va-eyebrow va-reveal">How it works</div>'
        '<h2 class="va-reveal">Every claim traced to its source.</h2>'
        f'<div class="va-grid va-two va-steps">{tiles(ABOUT_HOW, False)}</div>'
        "</div>"
        '<div class="va-note va-reveal">Variant Atlas summarises public research data. It is not a diagnostic tool and '
        "does not give medical advice; clinical decisions belong with qualified professionals.</div>"
        "</div>"
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


# ---------------------------------------------------------------- interpretation and scores

def _score_tile(label: str, value: str, note: str, kind: str = "neutral") -> str:
    css = {"pathogenic": "pathogenic", "benign": "benign", "uncertain": "uncertain"}.get(kind, "brand")
    return (
        f'<div class="va-fact {css}"><div class="va-fact-label">{e(label)}</div>'
        f'<div class="va-fact-value">{e(value)}</div><div class="va-fact-sub">{e(note)}</div></div>'
    )


def score_grid(item) -> str:
    """Prediction and context scores; each tile says what the number means."""
    tiles = []
    sc = getattr(item, "scores", None)
    if sc is not None and sc.revel is not None:
        kind = "pathogenic" if sc.revel >= 0.644 else "benign" if sc.revel <= 0.29 else "uncertain"
        tiles.append(_score_tile("REVEL", f"{sc.revel:.3f}", "missense pathogenicity, 0 to 1 (ClinGen-calibrated)", kind))
    if sc is not None and sc.alphamissense is not None:
        kind = {"likely pathogenic": "pathogenic", "likely benign": "benign"}.get(sc.alphamissense_class, "uncertain")
        tiles.append(_score_tile("AlphaMissense", f"{sc.alphamissense:.3f}", sc.alphamissense_class or "missense effect, 0 to 1", kind))
    if sc is not None and sc.cadd_phred is not None:
        kind = "pathogenic" if sc.cadd_phred >= 20 else "neutral"
        tiles.append(_score_tile("CADD", f"{sc.cadd_phred:.1f}", "deleteriousness, phred scale; 20 = top 1% of changes", kind))
    if sc is not None and sc.sift is not None:
        tiles.append(_score_tile("SIFT", f"{sc.sift:.2f}", sc.sift_pred or "below 0.05 is deleterious", "pathogenic" if sc.sift < 0.05 else "benign"))
    if sc is not None and sc.polyphen is not None:
        tiles.append(_score_tile("PolyPhen-2", f"{sc.polyphen:.2f}", sc.polyphen_pred or "HumDiv score", "pathogenic" if sc.polyphen > 0.85 else "benign" if sc.polyphen < 0.15 else "uncertain"))
    sp = getattr(item, "splice", None)
    if sp is not None and sp.max_delta is not None:
        kind = "pathogenic" if sp.max_delta >= 0.5 else "uncertain" if sp.max_delta >= 0.2 else "benign"
        tiles.append(_score_tile("SpliceAI", f"{sp.max_delta:.2f}", f"largest change: {sp.kind}" if sp.kind else "splice change, 0 to 1", kind))
    rg = getattr(item, "regulome", None)
    if rg is not None and rg.ranking:
        kind = "pathogenic" if rg.ranking[:1] in {"1", "2"} else "neutral"
        prob = f", probability {rg.probability:.2f}" if rg.probability is not None else ""
        tiles.append(_score_tile("RegulomeDB", rg.ranking, f"{rg.meaning}{prob}", kind))
    c = getattr(item, "constraint", None)
    if c is not None and c.loeuf is not None:
        kind = "pathogenic" if c.loeuf < 0.35 else "uncertain" if c.loeuf < 0.6 else "benign"
        pli = f"; pLI {c.pli:.2f}" if c.pli is not None else ""
        tiles.append(_score_tile(f"{c.gene} constraint", f"LOEUF {c.loeuf:.2f}", f"{c.verdict}{pli}", kind))
    if not tiles:
        return ""
    return f'<div class="va-facts">{"".join(tiles)}</div>'


_STRENGTH_ORDER = {"Check first": -1, "Stand-alone": 0, "Strong": 1, "Moderate": 2, "Supporting": 3, "To review": 4}


def acmg_html(criteria, not_computed: str) -> str:
    if not criteria:
        rows = '<div class="va-acmg-empty">No criterion can be computed from the public data for this variant.</div>'
    else:
        rows = "".join(
            f'<div class="va-acmg-row"><span class="va-acmg-code {e(c.direction)}">{e(c.code)}</span>'
            f'<span class="va-acmg-strength">{e(c.strength)}</span><span class="va-acmg-ev">{e(c.evidence)}</span>'
            f'<span class="va-src">{e(c.source)}</span></div>'
            for c in sorted(criteria, key=lambda c: _STRENGTH_ORDER.get(c.strength, 9))
        )
    return (
        f'<div class="va-card va-acmg">{rows}'
        f'<div class="va-acmg-note">Draft for a qualified reviewer, not a classification. {e(not_computed)}</div></div>'
    )


def pgx_html(pgx) -> str:
    if not pgx.annotations:
        return ""
    rows = []
    for a in pgx.annotations[:8]:
        drugs = ", ".join(a.drugs) or "drug not named"
        cats = " · ".join(a.categories)
        link = f' <a href="{e(a.url)}" target="_blank">ClinPGx</a>' if a.url else ""
        rows.append(
            f'<div class="va-pgx-row"><span class="va-level l{e(a.level[:1])}">Level {e(a.level or "?")}</span>'
            f'<div><b>{e(drugs)}</b>{" · " + e(", ".join(a.genes)) if a.genes else ""}'
            f'<div class="va-fact-sub">{e(cats)}{" · " if cats and a.phenotypes else ""}{e(", ".join(a.phenotypes[:3]))}</div>'
            f'<div class="va-pgx-sum">{e(a.summary)}{link}</div></div></div>'
        )
    return f'<div class="va-card">{"".join(rows)}<div class="va-acmg-note">Evidence levels 1A (strongest, guideline-backed) to 4 (case reports).</div></div>'


def civic_html(items) -> str:
    if not items:
        return ""
    rows = []
    for it in items[:8]:
        therapy = ", ".join(it.therapies)
        pmid = f' <a href="https://pubmed.ncbi.nlm.nih.gov/{e(it.pmid)}/" target="_blank">PMID:{e(it.pmid)}</a>' if it.pmid.isdigit() else ""
        rows.append(
            f'<div class="va-pgx-row"><span class="va-level civic">Level {e(it.level or "?")}</span>'
            f'<div><b>{e(it.disease or "Cancer")}</b>{" · " + e(therapy) if therapy else ""}'
            f'<div class="va-fact-sub">{e(it.evidence_type)} · {e(it.significance)}{" · " + e(it.direction) if it.direction else ""}</div>'
            f'<div class="va-pgx-sum">{e(it.statement)}{pmid}</div></div></div>'
        )
    return f'<div class="va-card">{"".join(rows)}<div class="va-acmg-note">CIViC levels A (validated) to E (inferential).</div></div>'


EXTRA_CSS = """
<style>
.va-acmg-row { display: grid; grid-template-columns: 64px 96px 1fr auto; gap: .8rem; align-items: center; padding: .6rem 0; border-bottom: 1px solid var(--va-border); }
.va-acmg-row:last-of-type { border-bottom: none; }
.va-acmg-code { font-family: "JetBrains Mono", monospace; font-weight: 600; font-size: .85rem; text-align: center; border-radius: 6px; padding: .2rem 0; }
.va-acmg-code.pathogenic { color: #ff8a8e; background: rgba(242,85,90,.12); }
.va-acmg-code.benign { color: #6fdca5; background: rgba(76,195,138,.12); }
.va-acmg-code.review { color: #f5c96a; background: rgba(240,183,60,.12); }
.va-acmg-strength { color: var(--va-muted); font-size: .82rem; }
.va-acmg-ev { color: #e6e8eb; font-size: .88rem; }
.va-acmg-note, .va-acmg-empty { color: var(--va-muted); font-size: .8rem; line-height: 1.5; margin-top: .7rem; }
.va-pgx-row { display: grid; grid-template-columns: 82px 1fr; gap: .9rem; padding: .7rem 0; border-bottom: 1px solid var(--va-border); color: #e6e8eb; font-size: .9rem; }
.va-pgx-row b { color: var(--va-text); font-weight: 600; }
.va-pgx-sum { color: #c9cdd3; font-size: .85rem; line-height: 1.5; margin-top: .3rem; }
.va-level { font-family: "JetBrains Mono", monospace; font-size: .74rem; border-radius: 6px; padding: .25rem .4rem; text-align: center; height: fit-content;
  color: #c7ccff; background: rgba(124,134,255,.14); }
.va-level.l1 { color: #6fdca5; background: rgba(76,195,138,.14); }
.va-level.l2 { color: #9ad8ff; background: rgba(90,162,255,.14); }
.va-level.civic { color: #f5c96a; background: rgba(240,183,60,.12); }
@media (max-width: 700px) { .va-acmg-row { grid-template-columns: 56px 1fr; } .va-acmg-ev, .va-acmg-row .va-src { grid-column: 1 / -1; } }
</style>
"""


def method(key: str) -> str:
    """A collapsed 'How this is computed' note for a result we compute ourselves."""
    from snp_brief.methods import method_items

    title, steps = method_items(key)
    items = "".join(f"<li>{e(step)}</li>" for step in steps)
    return f'<details class="va-method"><summary>How this is computed · {e(title)}</summary><ol>{items}</ol></details>'


METHOD_CSS = """
<style>
.va-method { margin: .5rem 0 1rem 0; border: 1px solid var(--va-border); border-radius: 10px; background: rgba(16,17,19,.9); }
.va-method summary { cursor: pointer; list-style: none; padding: .55rem .85rem; color: var(--va-muted); font-size: .8rem; }
.va-method summary::-webkit-details-marker { display: none; }
.va-method summary::before { content: "ƒ"; font-family: "JetBrains Mono", monospace; color: var(--va-accent-2); margin-right: .55rem; }
.va-method[open] summary { color: var(--va-text); border-bottom: 1px solid var(--va-border); }
.va-method ol { margin: .6rem 1rem .8rem 2.1rem; padding: 0; color: #c9cdd3; font-size: .84rem; line-height: 1.55; }
.va-method li { margin: .25rem 0; }
</style>
"""


def checks_html(report) -> str:
    """Result of the integrity checks run on this variant's data (cached or live)."""
    n = len(report.passed)
    if report.failures:
        state, label = "fail", f"{len(report.failures)} value{'s' if len(report.failures) != 1 else ''} failed a data check"
    elif report.warnings:
        state, label = "warn", f"{n} data checks passed · {len(report.warnings)} to look at"
    else:
        state, label = "pass", f"All {n} data checks passed"
    rows = "".join(f'<li class="f">{e(x)}</li>' for x in report.failures)
    rows += "".join(f'<li class="w">{e(x)}</li>' for x in report.warnings)
    rows += "".join(f'<li class="p">{e(x)}</li>' for x in report.passed)
    return (
        f'<details class="va-checks {state}"><summary><i></i>{e(label)}</summary>'
        f'<ul>{rows}</ul><div class="va-acmg-note">Every value is re-derived from the raw source responses or tested '
        "against the range it can take, each time the variant loads.</div></details>"
    )


CHECKS_CSS = """
<style>
.va-checks { margin: .2rem 0 1rem 0; border: 1px solid var(--va-border); border-radius: 10px; background: rgba(16,17,19,.9); }
.va-checks summary { cursor: pointer; list-style: none; padding: .55rem .85rem; font-size: .82rem; color: #d0d3d8; display: flex; align-items: center; gap: .55rem; }
.va-checks summary::-webkit-details-marker { display: none; }
.va-checks summary i { width: 8px; height: 8px; border-radius: 50%; background: #4cc38a; }
.va-checks.warn summary i { background: #f0b73c; } .va-checks.fail summary i { background: #f2555a; }
.va-checks.fail { border-color: rgba(242,85,90,.45); }
.va-checks ul { margin: .2rem 1rem .4rem 1rem; padding: 0; list-style: none; font-size: .84rem; line-height: 1.5; }
.va-checks li { padding: .2rem 0 .2rem 1.3rem; position: relative; color: #c9cdd3; }
.va-checks li::before { position: absolute; left: 0; font-family: "JetBrains Mono", monospace; }
.va-checks li.p::before { content: "✓"; color: #4cc38a; } .va-checks li.w::before { content: "!"; color: #f0b73c; }
.va-checks li.f::before { content: "×"; color: #f2555a; } .va-checks li.f { color: #ffb3b5; }
.va-checks .va-acmg-note { margin: 0 1rem .7rem 1rem; }
</style>
"""


def borzoi_html(item) -> str:
    """Predicted regulatory effect from our Borzoi run, and how it compares with measured GTEx eQTLs."""
    from snp_brief.sources.extra import concordance

    bz = getattr(item, "borzoi", None)
    if bz is None or not bz.retrieved:
        return ""
    parts = []
    for allele in bz.alleles:
        if not allele.ref_matches_genome:
            parts.append(f'<div class="va-acmg-empty">{e(allele.ref)}&gt;{e(allele.alt)}: the reference genome has a different base here, so it was not scored.</div>')
            continue
        if not allele.genes and not allele.chromatin:
            continue
        rows = []
        for g in allele.genes[:6]:
            direction = "higher" if g.mean_lfc > 0 else "lower"
            pct = (2 ** abs(g.mean_lfc) - 1) * 100
            change_text = (f"{pct:.2g}% {direction} expression (log2 {g.mean_lfc:+.3f})" if abs(g.mean_lfc) >= 0.005
                           else f"no clear change (log2 {g.mean_lfc:+.4f})")
            top_tissue = max(g.tissues.items(), key=lambda kv: abs(kv[1]))[0] if g.tissues else ""
            rows.append(
                f'<div class="va-bz-row"><b>{e(g.gene)}</b><span class="va-bz-bar {"up" if g.mean_lfc > 0 else "down"}" '
                f'style="width:{min(100, abs(g.mean_lfc) * 400):.0f}%"></span><span>{change_text}</span><span class="va-fact-sub">strongest in {e(top_tissue)}; {g.distance / 1000:.0f} kb away</span></div>'
            )
        local = ", ".join(f'{e(c["assay"])} {e(c["label"])} ({c["lfc"]:+.2f})' for c in allele.chromatin[:4])
        parts.append(
            f'<div class="va-card"><div class="va-fact-label">Predicted for the {e(allele.alt)} allele vs {e(allele.ref)}</div>'
            f'{"".join(rows) or "<div class=va-acmg-empty>No gene in range has a predicted change.</div>"}'
            + (f'<div class="va-pgx-sum" style="margin-top:.6rem">Largest local changes at the variant: {local}</div>' if local else "")
            + "</div>"
        )
    conc = concordance(getattr(item, "eqtls", None), bz)
    if conc:
        clear = [c for c in conc if c.clear]
        faint = len(conc) - len(clear)
        agree = sum(c.agrees for c in clear)
        rows = "".join(
            f'<div class="va-bz-conc"><span class="{"ok" if c.agrees else "no"}">{"agrees" if c.agrees else "disagrees"}</span>'
            f'<span><b>{e(c.gene)}</b> in {e(c.tissue)}: GTEx measured {"higher" if c.gtex_nes > 0 else "lower"} expression '
            f'(NES {c.gtex_nes:+.2f}); Borzoi predicts {"higher" if c.predicted > 0 else "lower"} (log2 {c.predicted:+.3f}, {e(c.matched_track)})</span></div>'
            for c in clear[:8]
        )
        head = (f"{agree} of {len(clear)} clear predictions agree with GTEx in direction" if clear
                else "No clear Borzoi prediction for this variant's GTEx genes")
        tail = (
            f'<div class="va-pgx-sum" style="margin-top:.5rem">For {faint} other GTEx eQTL{"s" if faint != 1 else ""}, Borzoi predicts '
            "no clear change (under 0.35%). That usually means this SNP is inherited together with the variant that changes "
            "expression rather than causing it; those cases are not counted as agreeing or disagreeing.</div>"
            if faint else ""
        )
        parts.append(f'<div class="va-card"><div class="va-fact-label">Prediction vs measurement · {e(head)}</div>{rows}{tail}</div>')
    return "".join(parts)


BORZOI_CSS = """
<style>
.va-bz-row { display: grid; grid-template-columns: 90px 120px 1fr; gap: .2rem .8rem; align-items: center; padding: .5rem 0; border-bottom: 1px solid var(--va-border); color: #e6e8eb; font-size: .88rem; }
.va-bz-row .va-fact-sub { grid-column: 3; margin-top: 0; }
.va-bz-bar { height: 6px; border-radius: 3px; min-width: 3px; }
.va-bz-bar.up { background: #4cc9b0; } .va-bz-bar.down { background: #f2555a; }
.va-bz-conc { display: grid; grid-template-columns: 82px 1fr; gap: .8rem; padding: .45rem 0; font-size: .86rem; color: #c9cdd3; border-bottom: 1px solid var(--va-border); }
.va-bz-conc span.ok { color: #6fdca5; } .va-bz-conc span.no { color: #ff8a8e; }
.va-bz-conc b { color: var(--va-text); }
</style>
"""


def zygosity_html(rows, inheritance_loaded: bool) -> str:
    """First thing on the page: each allele, its own classification, and what one vs two copies mean."""
    if not rows:
        return ""
    cards = []
    for z in rows:
        label = f"{e(z.allele)}" + (f" · {e(z.change)}" if z.change else "")
        cls = z.classification or "No ClinVar classification for this allele"
        kind = classification_kind(z.classification) if z.classification else "neutral"
        cond_rows = []
        for c in z.conditions[:4]:
            mode = " / ".join(c.modes) if c.modes else "inheritance not stated"
            src = f'<span class="va-src">{e(c.mode_source)}</span>' if c.mode_source else ""
            copies = (
                f'<div class="va-zy-copies"><span><i>1 copy</i>{e(c.one_copy)}</span><span><i>2 copies</i>{e(c.two_copies)}</span></div>'
                if c.modes else ""
            )
            cond_rows.append(f'<div class="va-zy-cond"><b>{e(c.name)}</b><span class="va-zy-mode">{e(mode)}</span>{src}{copies}</div>')
        counts = ""
        if z.carriers is not None:
            hom = f"{z.homozygotes:,}" if z.homozygotes is not None else "not in saved data"
            counts = (f'<div class="va-zy-counts">gnomAD ({e(z.callset)}): {z.carriers:,} copies observed · '
                      f'people with two copies: {e(hom)}' + (f" · hemizygous: {z.hemizygotes:,}" if z.hemizygotes else "") + "</div>")
        flags = "".join(f'<div class="va-zy-flag">{e(f)}</div>' for f in z.flags)
        cards.append(
            f'<div class="va-zy-card"><div class="va-zy-head"><span class="va-zy-allele">{label}</span>'
            f'<span class="va-chip {kind}">{e(cls[:80])}</span></div>{"".join(cond_rows)}{counts}{flags}</div>'
        )
    note = "" if inheritance_loaded else (
        '<div class="va-acmg-note">Inheritance data (HPO) is not loaded, so modes come only from condition names.</div>'
    )
    return f'<div class="va-zy">{"".join(cards)}{note}</div>'


ZYGOSITY_CSS = """
<style>
.va-zy { display: grid; gap: .6rem; margin: 0 0 1rem 0; }
.va-zy-card { background: var(--va-surface); border: 1px solid var(--va-border); border-radius: 12px; padding: .8rem 1rem; }
.va-zy-head { display: flex; align-items: center; gap: .7rem; flex-wrap: wrap; margin-bottom: .3rem; }
.va-zy-allele { font-family: "JetBrains Mono", monospace; font-weight: 600; color: var(--va-text); font-size: .95rem; }
.va-zy-cond { padding: .5rem 0; border-top: 1px solid var(--va-border); color: #e6e8eb; font-size: .88rem; }
.va-zy-cond b { font-weight: 600; margin-right: .6rem; }
.va-zy-mode { color: var(--va-accent-2); font-size: .82rem; margin-right: .5rem; }
.va-zy-copies { display: grid; grid-template-columns: 1fr 1fr; gap: .6rem; margin-top: .35rem; }
.va-zy-copies span { background: rgba(255,255,255,.03); border: 1px solid var(--va-border); border-radius: 8px; padding: .35rem .55rem; color: #c9cdd3; font-size: .83rem; }
.va-zy-copies i { display: block; font-style: normal; font-size: .7rem; text-transform: uppercase; letter-spacing: .06em; color: var(--va-muted); }
.va-zy-counts { color: var(--va-muted); font-size: .82rem; margin-top: .45rem; }
.va-zy-flag { margin-top: .45rem; color: #f5c96a; font-size: .84rem; padding-left: 1.1rem; position: relative; }
.va-zy-flag::before { content: "!"; position: absolute; left: .2rem; font-family: "JetBrains Mono", monospace; }
@media (max-width: 700px) { .va-zy-copies { grid-template-columns: 1fr; } }
</style>
"""
