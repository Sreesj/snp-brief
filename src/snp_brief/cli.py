"""One command: print a plain-text brief for each rsID in a list."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from snp_brief.fixture import read_rsids
from snp_brief.pipeline import build_brief


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
    args = parser.parse_args(argv)
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
        if args.refresh:
            print(rsid, file=sys.stderr)
        paragraphs.append(
            build_brief(rsid, cache_dir=args.cache, refresh=args.refresh, offline=args.offline)
        )
    text = "\n\n".join(paragraphs) + "\n"
    sys.stdout.write(text)
    if args.save:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        args.save.write_text(text, encoding="utf-8")
    return 0
