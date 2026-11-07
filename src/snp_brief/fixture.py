"""Load the pinned rsID list. One id per line."""

from __future__ import annotations

import re
from pathlib import Path

_RSID = re.compile(r"^rs[1-9][0-9]*$")


def is_rsid(text: str) -> bool:
    return bool(_RSID.fullmatch(text))


def read_rsids(path: Path) -> list[str]:
    if not path.is_file():
        raise FileNotFoundError(f"missing rsID list: {path}")
    ids: list[str] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if not _RSID.fullmatch(line):
            raise ValueError(f"not an rsID: {line}")
        ids.append(line)
    if not ids:
        raise ValueError(f"no rsIDs in {path}")
    return ids
