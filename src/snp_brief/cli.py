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
    parser.add_argument("--regenerate", action="store_true", help="ignore cached LLM summaries and Forge results")
    parser.add_argument("--list-models", action="store_true", help="list the models your LLM provider offers, then exit")
    args = parser.parse_args(argv)
    if args.list_models:
        return list_models()
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
        item = gather(rsid, cache_dir=args.cache, refresh=args.refresh, offline=args.offline)
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
    lines += [f"  Next: {step}" for step in tri.next_steps]
    for number, hyp in enumerate(result.hypotheses, 1):
        lines.append(f"  Hypothesis {number} ({hyp.feasibility} feasibility, {hyp.consistency}): {hyp.title}. {hyp.mechanism}")
        lines += [f"    Support: {claim['rendered']}" for claim in hyp.support]
        lines += [f"    Gap: {gap}" for gap in hyp.gaps]
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
