"""One command: print a plain-text brief for each rsID in a list."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from snp_brief.fixture import read_rsids
from snp_brief.pipeline import gather
from snp_brief.summary import summarize


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="snp-brief",
        description="Print a sourced plain-text brief for each SNP rsID.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("fixture/snps.txt"),
        help="rsID list, one id per line (default: fixture/snps.txt)",
    )
    parser.add_argument("--save", type=Path, help="also write the briefs to this file")
    parser.add_argument("--cache", type=Path, default=Path("cache"), help="directory of cached source responses")
    parser.add_argument("--refresh", action="store_true", help="call the live sources and replace the cache")
    parser.add_argument("--offline", action="store_true", help="use the cache only")
    parser.add_argument(
        "--summary",
        choices=("template", "llm", "off"),
        default="template",
        help="plain-English summary above each evidence brief: template (default, offline), "
        "llm (Nebius via NEBIUS_API_KEY, or SNP_BRIEF_LLM_URL / SNP_BRIEF_LLM_MODEL; falls back to the "
        "template if it fails the fact check), or off",
    )
    parser.add_argument(
        "--forge",
        action="store_true",
        help="also print Hypothesis Forge: stage 1 category, and LLM hypotheses when an LLM is connected",
    )
    parser.add_argument(
        "--speed",
        choices=("fast", "balanced", "thorough", "complete"),
        default="complete",
        help="time budget for live providers: fast ~8 s, balanced ~20 s, thorough ~45 s, complete waits for all (default)",
    )
    parser.add_argument("--findings", action="store_true", help="extract quote-verified findings from paper abstracts (LLM)")
    parser.add_argument("--regenerate", action="store_true", help="ignore cached LLM summaries and Forge results")
    parser.add_argument("--probe", metavar="RSID", help="call every provider live for one rsID and print status, time and response shape")
    parser.add_argument("--list-models", action="store_true", help="list the models your LLM provider offers, then exit")
    parser.add_argument("--build-gwas", type=Path, metavar="FILE", help="slim the GWAS Catalog download (.tsv or .zip) into data/gwas_catalog.tsv.gz")
    args = parser.parse_args(argv)
    if args.list_models:
        return list_models()
    if args.probe:
        return probe(args.probe)
    if args.build_gwas:
        from snp_brief.sources import gwas_local

        out, rows = gwas_local.build(args.build_gwas, Path.cwd())
        print(f"wrote {out} ({rows:,} single-SNP associations, {out.stat().st_size / 1e6:.1f} MB)")
        return 0
    if args.refresh and args.offline:
        print("Choose either --refresh or --offline.", file=sys.stderr)
        return 2
    try:
        rsids = read_rsids(args.input)
    except (FileNotFoundError, ValueError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    paragraphs = []
    for rsid in rsids:
        if args.refresh or args.summary == "llm" or args.forge:
            print(rsid, file=sys.stderr)
        from snp_brief.providers import SPEED_BUDGETS

        item = gather(
            rsid, cache_dir=args.cache, refresh=args.refresh, offline=args.offline, budget=SPEED_BUDGETS.get(args.speed)
        )
        if item.pending:
            print(f"{rsid}: still loading after the time budget: {', '.join(item.pending)}", file=sys.stderr)
        if args.findings:
            from snp_brief.literature import Finding, extract_findings
            from snp_brief.llm import LLMClient

            client = LLMClient.from_env()
            found = extract_findings(item, client.for_forge() if client else None, args.cache, refresh=args.regenerate, offline=args.offline)
            item.findings = [Finding(**f) for f in found.get("findings", [])]
            for r in found.get("rejected", []):
                print(f"{rsid}: finding rejected ({r['reason']}): {r['claim']}", file=sys.stderr)
        block = format_item(item, args.summary, cache_dir=args.cache, regenerate=args.regenerate)
        if args.forge:
            block += "\n" + format_forge(item, cache_dir=args.cache, regenerate=args.regenerate)
        paragraphs.append(block)
    text = "\n\n".join(paragraphs) + "\n"
    sys.stdout.write(text)
    if args.save:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        args.save.write_text(text, encoding="utf-8")
    return 0


def format_item(item, mode: str, cache_dir: Path | None = None, regenerate: bool = False) -> str:
    """The evidence brief, with the checked plain-English summary above it unless mode is off."""
    if mode == "off":
        return item.text
    summary = summarize(item, mode, cache_dir=cache_dir, refresh=regenerate)
    lines = [f"{item.rsid}", f"Summary ({summary.writer}): {summary.text}"]
    lines += [f"  Rejected by checker: {sentence} ({reason})" for sentence, reason in summary.rejected]
    lines.append(f"Evidence: {item.text}")
    return "\n".join(lines)


def format_forge(item, cache_dir: Path | None = None, regenerate: bool = False) -> str:
    from snp_brief.forge import forge

    result = forge(item, cache_dir=cache_dir, refresh=regenerate)
    tri = result.triage
    lines = [f"Forge: {tri.category}. {tri.reason}", f"Open problem: {tri.open_problem}"]
    if tri.exploratory:
        lines.append("  Exploratory: AI-generated from genomic context, not from evidence about this variant.")
    lines += [f"  Next: {step}" for step in tri.next_steps]
    for number, hyp in enumerate(result.hypotheses, 1):
        lines.append(f"  Hypothesis {number} ({hyp.feasibility} feasibility, {hyp.consistency}): {hyp.title}. {hyp.mechanism}")
        lines += [f"    Support: {claim['rendered']}" for claim in hyp.support]
        lines += [f"    Gap: {gap}" for gap in hyp.gaps]
        lines += [f"    AI background (unverified): {line}" for line in hyp.background]
        lines.append(f"    Distinguish by: {hyp.discriminating_observation}")
        lines.append(f"    Experiment: {hyp.experiment}")
    for row in result.ranking:
        lines.append(
            f"  Rank {row['rank']}: {row['experiment']} (tests {row['tests']}; {row['consistency']}, "
            f"{row['support']} supporting facts, {row['feasibility']} feasibility)"
        )
    if result.note:
        lines.append(f"  Note: {result.note}")
    if result.hypotheses:
        lines.append(f"  Written by {result.writer}.")
    return "\n".join(lines)


def list_models() -> int:
    from snp_brief.llm import LLMClient

    client = LLMClient.from_env()
    if client is None:
        print("No LLM configured. Set NEBIUS_API_KEY (Nebius Token Factory) or SNP_BRIEF_LLM_URL.", file=sys.stderr)
        return 2
    try:
        models = client.list_models()
    except Exception as exc:
        print(f"Could not list models: {exc}", file=sys.stderr)
        return 1
    print(f"{client.provider}: {len(models)} models (current: {client.model})")
    for model in models:
        print(model)
    return 0


def probe(rsid: str) -> int:
    """Live check of every provider: status, seconds, and the shape of what came back (for fixing parsers)."""
    import json
    import time

    from snp_brief import providers as p
    from snp_brief.sources import extra as x
    from snp_brief.sources.dbsnp import fetch_dbsnp, parse_dbsnp

    try:
        db = parse_dbsnp(fetch_dbsnp(rsid), rsid)
        print(f"dbsnp: OK; {db.assembly} chr{db.chrom}:{db.pos}, genes {db.genes}")
    except Exception as exc:
        print(f"dbsnp: FAILED: {str(exc)[:300]}")
        db = parse_dbsnp(None, rsid)
    coords = [(db.chrom, db.pos, ref, alt) for ref, alt in db.alleles] if db.assembly == "GRCh38" else []
    checks = {
        "gwas": lambda: p.fetch_gwas(rsid),
        "ensembl_vep": lambda: p.fetch_vep(rsid),
        "gtex": lambda: p.fetch_gtex(rsid, [p.gtex_variant_id(*c) for c in coords]),
        "ensembl_ld": lambda: p.fetch_ld(rsid),
        "open_targets": lambda: p.fetch_open_targets([p.ot_variant_id(*c) for c in coords]),
        "litvar": lambda: p.fetch_litvar(rsid),
        "dbsnp_citations": lambda: p.fetch_citations(rsid),
        "ncbi_gene": lambda: p.fetch_gene_summary(db.genes[0]) if db.genes else {"response": "no gene in dbSNP"},
        "myvariant": lambda: x.fetch_myvariant(rsid),
        "pharmgkb": lambda: x.fetch_pharmgkb(rsid),
        "spliceai": lambda: x.fetch_spliceai([x.splice_variant_ids(*c) for c in coords]),
        "regulome": lambda: x.fetch_regulome(rsid),
        "constraint": lambda: x.fetch_constraint(db.genes[0]) if db.genes else {"response": "no gene in dbSNP"},
        "timeline": lambda: x.fetch_timeline(rsid),
    }
    for name, call in checks.items():
        started = time.monotonic()
        try:
            payload = call()
            shape = _shape(payload.get("response"))
            print(f"{name}: OK in {time.monotonic() - started:.1f} s; shape {shape}")
            print("   sample: " + json.dumps(payload.get("response"), default=str)[:600])
        except Exception as exc:
            print(f"{name}: FAILED in {time.monotonic() - started:.1f} s: {str(exc)[:300]}")
    return 0


def _shape(value, depth: int = 0) -> str:
    if depth > 3:
        return "..."
    if isinstance(value, dict):
        return "{" + ", ".join(f"{k}: {_shape(v, depth + 1)}" for k, v in list(value.items())[:12]) + "}"
    if isinstance(value, list):
        return f"[{len(value)} x {_shape(value[0], depth + 1)}]" if value else "[]"
    return type(value).__name__
