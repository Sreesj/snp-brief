"""Read and write cached source responses."""

from __future__ import annotations

import json
from pathlib import Path


def cache_path(cache_dir: Path, rsid: str, source: str) -> Path:
    return cache_dir / rsid / f"{source}.json"


def load_cache(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def save_cache(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
