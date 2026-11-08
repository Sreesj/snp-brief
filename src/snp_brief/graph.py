"""Evidence graph as Graphviz DOT: variant → genes → traits, styled by evidence type.

Edge styles keep a chain of predictions from looking like an established pathway:
- solid green: observed (reported by a source)
- dashed violet: predicted (computational score or model)
- dotted grey: hypothesis (Hypothesis Forge)
"""

from __future__ import annotations

from snp_brief.summary import _plan

OBSERVED = 'color="#4cc38a", fontcolor="#8a8f98", style=solid'
PREDICTED = 'color="#a99cff", fontcolor="#8a8f98", style=dashed'
HYPOTHESIS = 'color="#62666d", fontcolor="#62666d", style=dotted'


def evidence_dot(item, forge_result=None) -> str:
    plan = _plan(item)
    nodes: dict[str, str] = {}
    edges: list[tuple[str, str, str, str]] = []

    def node(name: str, kind: str) -> str:
        key = f"{kind}:{name}"
        if key not in nodes:
            shape = {"variant": "doubleoctagon", "gene": "box", "trait": "ellipse", "feature": "note", "hyp": "plaintext"}[kind]
            nodes[key] = f'"{_esc(key)}" [label="{_esc(name)}", shape={shape}];'
        return key

    variant = node(item.rsid, "variant")
    for gene in [g.strip() for g in (plan.gene or "").split(",") if g.strip()][:2]:
        edges.append((variant, node(gene, "gene"), _esc(plan.consequence or "located in/near"), OBSERVED))
    for eq in (getattr(item, "eqtls", None) or [])[:3]:
        edges.append((variant, node(eq.gene, "gene"), _esc(f"eQTL, {eq.tissue}"), OBSERVED))
    for row in (getattr(item, "l2g", None) or [])[:2]:
        edges.append((variant, node(row.gene, "gene"), f"L2G {row.score:.2f}", PREDICTED))
    for trait in plan.conditions[:3]:
        edges.append((variant, node(trait, "trait"), _esc(plan.headline or "ClinVar"), OBSERVED))
    for trait in plan.gwas_traits[:3]:
        edges.append((variant, node(trait, "trait"), "GWAS association", OBSERVED))
    vep = getattr(item, "vep", None)
    if vep is not None and vep.regulatory:
        edges.append((variant, node(vep.regulatory[0], "feature"), "overlaps", OBSERVED))
    if vep is not None and vep.motifs:
        edges.append((variant, node(f"{vep.motifs[0]} motif", "feature"), "predicted motif", PREDICTED))
    ld = getattr(item, "ld_proxies", None) or []
    if ld:
        edges.append((variant, node(f"{len(ld)} linked variants", "feature"), "LD r² ≥ 0.8", OBSERVED))
    if forge_result is not None:
        for number, hyp in enumerate(getattr(forge_result, "hypotheses", [])[:3], 1):
            edges.append((variant, node(f"H{number}: {hyp.title[:40]}", "hyp"), "hypothesis", HYPOTHESIS))
    lines = [
        "digraph evidence {",
        '  rankdir=LR; bgcolor="transparent"; pad=0.2;',
        '  node [fontsize=11, fontname="Inter", color="#2a2c31", fontcolor="#e6e8eb", style=filled, fillcolor="#16171a", penwidth=1];',
        '  edge [fontsize=9, fontname="Inter", arrowsize=0.6];',
    ]
    lines += [f"  {line}" for line in nodes.values()]
    seen = set()
    for a, b, label, style in edges:
        if (a, b, label) in seen:
            continue
        seen.add((a, b, label))
        lines.append(f'  "{_esc(a)}" -> "{_esc(b)}" [label="{label}", {style}];')
    lines.append("}")
    return "\n".join(lines)


def _esc(text: str) -> str:
    return str(text).replace("\\", "\\\\").replace('"', "'")
