"""Per-module result cache, next to the shared source cache.

    cache/<rsid>/dbsnp.json, clinvar.json, ...   shared evidence (sources)
    cache/<rsid>/summary_llm.json                LLM summary
    cache/<rsid>/forge.json                      Hypothesis Forge result

Each module result is stored with a key built from what it depends on: the
pipeline version, the model, and a hash of the facts it was written from. A
result is reused only when that key matches, so a new model or changed
evidence never silently reuses an old answer.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

PIPELINE_VERSION = "0.2"


def result_key(module: str, model: str, facts_text: str) -> str:
    digest = hashlib.sha256(facts_text.encode("utf-8")).hexdigest()[:16]
    return f"{module}|{PIPELINE_VERSION}|{model}|{digest}"


def load_result(cache_dir: Path | None, rsid: str, module: str, key: str) -> dict | None:
    if cache_dir is None:
        return None
    path = Path(cache_dir) / rsid / f"{module}.json"
    if not path.is_file():
        return None
    try:
        stored = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if stored.get("key") != key:
        return None
    return stored


def save_result(cache_dir: Path | None, rsid: str, module: str, key: str, result: dict) -> dict:
    stored = {"key": key, "created": datetime.now(timezone.utc).replace(microsecond=0).isoformat(), **result}
    if cache_dir is None:
        return stored
    path = Path(cache_dir) / rsid / f"{module}.json"
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(stored, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    except OSError:
        pass  # a read-only deployment still shows the result; it is just not saved
    return stored
