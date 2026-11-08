"""JSON HTTP with a short retry. NCBI calls are spaced for the public rate limit."""

from __future__ import annotations

import json
import os
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

USER_AGENT = "snp-brief/0.1 (research demo)"
_last_ncbi = 0.0
_ncbi_lock = threading.Lock()


def _pace_ncbi(url: str) -> None:
    """Space NCBI calls 0.35 s apart, across threads (public limit is 3 per second)."""
    global _last_ncbi
    host = urllib.parse.urlparse(url).netloc
    if not host.endswith("ncbi.nlm.nih.gov"):
        return
    with _ncbi_lock:
        spacing = 0.11 if os.environ.get("NCBI_API_KEY") else 0.35  # 10/s with a key, 3/s without
        wait = spacing - (time.monotonic() - _last_ncbi)
        if wait > 0:
            time.sleep(wait)
        _last_ncbi = time.monotonic()


def request_text(
    url: str,
    *,
    params: dict[str, str] | None = None,
    payload: dict | None = None,
    timeout: float = 60,
) -> str:
    return _request(url, params=params, payload=payload, timeout=timeout, decode_json=False)


def request_json(
    url: str,
    *,
    params: dict[str, str] | None = None,
    payload: dict | None = None,
    timeout: float = 60,
) -> dict:
    return _request(url, params=params, payload=payload, timeout=timeout, decode_json=True)


def _request(
    url: str,
    *,
    params: dict[str, str] | None,
    payload: dict | None,
    timeout: float,
    decode_json: bool,
):
    host = urllib.parse.urlparse(url).netloc
    if host == "eutils.ncbi.nlm.nih.gov" and os.environ.get("NCBI_API_KEY"):
        params = {**(params or {}), "api_key": os.environ["NCBI_API_KEY"]}
    if params:
        joiner = "&" if "?" in url else "?"
        url = f"{url}{joiner}{urllib.parse.urlencode(params)}"
    data = None
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    last_error: Exception | None = None
    for attempt in range(3):
        _pace_ncbi(url)
        req = urllib.request.Request(url, data=data, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                text = resp.read().decode()
                if not decode_json:
                    return text
                return json.loads(text)
        except urllib.error.HTTPError as exc:
            body = exc.read().decode(errors="replace")
            last_error = exc
            if exc.code not in {429, 500, 502, 503, 504} or attempt == 2:
                raise RuntimeError(f"HTTP {exc.code} for {url}: {body[:400]}") from exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt == 2:
                raise RuntimeError(f"request failed for {url}: {exc}") from exc
        time.sleep(2**attempt)
    raise RuntimeError(f"request failed for {url}: {last_error}")
