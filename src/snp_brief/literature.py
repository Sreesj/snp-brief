"""Findings from paper abstracts, extracted by the LLM and verified against the abstract text.

1. Abstracts for up to five PMIDs (Europe PMC papers first, then LitVar2, then
   dbSNP citations) are fetched from Europe PMC (resultType=core) and cached as
   cache/<rsid>/abstracts.json.
2. The LLM extracts at most two findings per abstract that concern this rsID:
   the claim, an exact supporting quote, the population, and an effect size.
3. A finding is kept only if its quote appears verbatim in the abstract and
   every number in the claim appears in the quote. Kept findings are cached as
   cache/<rsid>/findings.json and become facts (rank 7, literature) for the
   summary and Hypothesis Forge.
"""

from __future__ import annotations

import html
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from snp_brief.http_client import request_json
from snp_brief.module_cache import load_result, result_key, save_result

SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
_TAG = re.compile(r"<[^>]+>")
_NUMBER = re.compile(r"(?<![A-Za-z])\d+(?:[.,]\d+)*(?:e-?\d+)?%?")

PROMPT = """You extract findings about one genetic variant from paper abstracts.
For each abstract, return at most two findings that are specifically about the given rsID (or an alias listed).
Each finding needs:
- "pmid": the abstract's PMID
- "claim": one sentence stating the finding in plain English
- "quote": an exact, verbatim span copied from the abstract that supports the claim (no paraphrase, no ellipsis)
- "population": who was studied, if stated (else "")
- "effect": the effect size with its unit as written (odds ratio, beta, p-value), if stated (else "")
Every number in "claim" must also appear in "quote". If an abstract says nothing specific about the variant, return nothing for it.
Do not infer study design beyond what the abstract states. Do not give medical advice.
Reply with only: {"findings": [{"pmid": "...", "claim": "...", "quote": "...", "population": "...", "effect": "..."}]}"""


@dataclass
class Abstract:
    pmid: str
    title: str
    year: str
    text: str


@dataclass
class Finding:
    pmid: str
    claim: str
    quote: str
    population: str
    effect: str


def candidate_pmids(item, limit: int = 5) -> list[str]:
    pmids: list[str] = []
    for source in (
        [p.pmid for p in item.papers],
        getattr(item, "litvar", []) or [],
        getattr(item, "citations", []) or [],
        getattr(getattr(item, "vep", None), "pubmed", []) or [],
    ):
        for pmid in source:
            if pmid and str(pmid) not in pmids:
                pmids.append(str(pmid))
    return pmids[:limit]


def fetch_abstracts(pmids: list[str], timeout: float = 20) -> dict:
    query = " OR ".join(f"EXT_ID:{pmid}" for pmid in pmids) + " AND SRC:MED"
    body = request_json(
        SEARCH, params={"query": query, "format": "json", "resultType": "core", "pageSize": str(len(pmids))}, timeout=timeout
    )
    return {"source": "Europe PMC abstracts", "request": {"pmids": pmids}, "response": body}


def parse_abstracts(payload: dict | None) -> list[Abstract]:
    if not payload:
        return []
    rows = (((payload.get("response") or {}).get("resultList") or {}).get("result")) or []
    out = []
    for row in rows:
        text = html.unescape(_TAG.sub(" ", str(row.get("abstractText") or ""))).strip()
        text = re.sub(r"\s+", " ", text)
        if row.get("pmid") and text:
            out.append(
                Abstract(
                    str(row.get("pmid")),
                    html.unescape(_TAG.sub("", str(row.get("title") or ""))).strip(),
                    str(row.get("pubYear") or ""),
                    text,
                )
            )
    return out


def load_abstracts(item, cache_dir: Path | None, *, offline: bool = False) -> list[Abstract]:
    pmids = candidate_pmids(item)
    if not pmids:
        return []
    path = Path(cache_dir) / item.rsid / "abstracts.json" if cache_dir else None
    if path is not None and path.is_file():
        try:
            stored = json.loads(path.read_text(encoding="utf-8"))
            if stored.get("request", {}).get("pmids") == pmids:
                return parse_abstracts(stored)
        except (OSError, ValueError):
            pass
    if offline:
        return []
    payload = fetch_abstracts(pmids)
    if path is not None:
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        except OSError:
            pass
    return parse_abstracts(payload)


def verify(raw: dict, abstracts: list[Abstract], aliases: list[str] | None = None) -> tuple[list[Finding], list[dict]]:
    """Keep findings whose quote is verbatim in the abstract and whose numbers are in the quote."""
    by_pmid = {a.pmid: a for a in abstracts}
    kept: list[Finding] = []
    rejected: list[dict] = []
    for row in (raw.get("findings") or [])[:10]:
        if not isinstance(row, dict):
            continue
        pmid = str(row.get("pmid") or "").replace("PMID:", "").strip()
        claim = str(row.get("claim") or "").strip()
        quote = str(row.get("quote") or "").strip().strip('"')
        abstract = by_pmid.get(pmid)
        reason = ""
        if abstract is None:
            reason = "PMID not among the abstracts"
        elif not quote or _norm(quote) not in _norm(abstract.text):
            reason = "quote not found verbatim in the abstract"
        else:
            quote_numbers = {_clean(n) for n in _NUMBER.findall(quote)}
            missing = [n for n in _NUMBER.findall(claim) if _clean(n) not in quote_numbers]
            if missing:
                reason = f"number {missing[0]} is not in the quote"
        if reason:
            rejected.append({"pmid": pmid, "claim": claim, "reason": reason})
            continue
        kept.append(Finding(pmid, claim, quote, str(row.get("population") or ""), str(row.get("effect") or "")))
    return kept, rejected


def extract_findings(item, client, cache_dir: Path | None, *, refresh: bool = False, offline: bool = False) -> dict:
    """Returns {"findings": [...], "rejected": [...], "writer": "...", "abstracts": n}; cached per abstract set."""
    abstracts = load_abstracts(item, cache_dir, offline=offline)
    if not abstracts:
        return {"findings": [], "rejected": [], "writer": "", "abstracts": 0, "note": "No abstracts available."}
    basis = "\n".join(f"{a.pmid}:{a.text}" for a in abstracts)
    model = getattr(client, "model", "") if client is not None else ""
    key = result_key("findings", model, basis)
    if not refresh:
        stored = load_result(cache_dir, item.rsid, "findings", key) if client is not None else _any(cache_dir, item.rsid)
        if stored is not None:
            return {**stored, "writer": stored.get("writer", "") + " (cached)"}
    if client is None:
        return {"findings": [], "rejected": [], "writer": "", "abstracts": len(abstracts), "note": "Connect an LLM to extract findings."}
    from snp_brief.llm import parse_json

    aliases = [a for record in item.records for a in getattr(record, "aliases", [])]
    user = (
        f"Variant: {item.rsid}" + (f" (also known as {', '.join(aliases[:3])})" if aliases else "") + "\n\n"
        + "\n\n".join(f"PMID {a.pmid} ({a.year}) {a.title}\n{a.text}" for a in abstracts)
    )
    try:
        raw = parse_json(client.chat(PROMPT, user, json_mode=True, max_tokens=3000))
    except Exception as exc:
        return {"findings": [], "rejected": [], "writer": "", "abstracts": len(abstracts), "note": f"Extraction failed ({exc})."}
    kept, rejected = verify(raw, abstracts, aliases)
    result = {
        "findings": [asdict(f) for f in kept],
        "rejected": rejected,
        "writer": f"llm:{client.model}",
        "abstracts": len(abstracts),
    }
    return save_result(cache_dir, item.rsid, "findings", key, result)


def cached_findings(cache_dirs: list[Path | None], rsid: str) -> list[Finding]:
    """Verified findings already on disk, for the summary and Forge (no network, no LLM)."""
    for directory in cache_dirs:
        stored = _any(directory, rsid)
        if stored:
            return [Finding(**f) for f in stored.get("findings", [])]
    return []


def _any(cache_dir: Path | None, rsid: str) -> dict | None:
    if cache_dir is None:
        return None
    path = Path(cache_dir) / rsid / "findings.json"
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _norm(text: str) -> str:
    return re.sub(r"\s+", " ", text.replace(" ", " ").replace(" ", " ")).strip().lower()


def _clean(number: str) -> str:
    return number.replace(",", "").rstrip("%").rstrip(".")
