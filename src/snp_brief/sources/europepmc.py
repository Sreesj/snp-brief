"""Europe PMC: a few primary research papers that name the rsID in the title or abstract."""

from __future__ import annotations

import html
import re
from datetime import datetime, timezone

from snp_brief.http_client import request_json
from snp_brief.models import Paper

_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
_TAG = re.compile(r"<[^>]+>")
_REJECT = {"review", "review article", "systematic review", "editorial", "letter", "comment"}


def fetch_europepmc(rsid: str) -> dict:
    query = (
        f"(TITLE:{rsid} OR ABSTRACT:{rsid}) AND PUB_TYPE:\"research-article\" NOT PUB_TYPE:review"
    )
    body = request_json(
        _SEARCH,
        params={
            "query": query,
            "format": "json",
            "pageSize": "15",
            "resultType": "lite",
            "sort": "CITED desc",
        },
    )
    return {
        "source": "europepmc",
        "retrieved": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "request": {"query": query},
        "response": body,
    }


def parse_europepmc(payload: dict | None, limit: int = 3) -> tuple[bool, list[Paper]]:
    if not payload:
        return False, []
    hits = ((payload.get("response") or {}).get("resultList") or {}).get("result") or []
    papers: list[Paper] = []
    ranked = sorted(hits, key=lambda hit: int(hit.get("citedByCount") or 0), reverse=True)
    seen: set[str] = set()
    for hit in ranked:
        if not _primary(hit.get("pubType") or ""):
            continue
        pmid = str(hit.get("pmid") or "").strip()
        title = _plain(hit.get("title") or "")
        if not pmid or not title or pmid in seen:
            continue
        seen.add(pmid)
        papers.append(
            Paper(
                authors=_plain(hit.get("authorString") or ""),
                year=str(hit.get("pubYear") or "").strip(),
                title=title,
                pmid=pmid,
            )
        )
        if len(papers) >= limit:
            break
    return True, papers


def _primary(pub_type: str) -> bool:
    tokens = [token.strip().lower() for token in pub_type.split(";") if token.strip()]
    if any(token in _REJECT for token in tokens):
        return False
    return any(token in {"research-article", "journal article"} for token in tokens)


def _plain(value: str) -> str:
    return " ".join(html.unescape(_TAG.sub("", value)).split())
