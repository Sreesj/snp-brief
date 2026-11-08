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


_gnomad_lock = threading.Lock()
_last_gnomad = 0.0


def _pace_gnomad(url: str) -> None:
    """gnomAD's API rate-limits bursts; keep calls at least 1 s apart across threads."""
    global _last_gnomad
    if "gnomad.broadinstitute.org" not in url:
        return
    with _gnomad_lock:
        wait = 1.0 - (time.monotonic() - _last_gnomad)
        if wait > 0:
            time.sleep(wait)
        _last_gnomad = time.monotonic()


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
    attempts: int = 3,
) -> str:
    return _request(url, params=params, payload=payload, timeout=timeout, decode_json=False, attempts=attempts)


def request_json(
    url: str,
    *,
    params: dict[str, str] | None = None,
    payload: dict | None = None,
    timeout: float = 60,
    attempts: int = 3,
) -> dict:
    """attempts=1 for optional providers: a timeout should not triple the wait."""
    return _request(url, params=params, payload=payload, timeout=timeout, decode_json=True, attempts=attempts)


def _request(
    url: str,
    *,
    params: dict[str, str] | None,
    payload: dict | None,
    timeout: float,
    decode_json: bool,
    attempts: int = 3,
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
    for attempt in range(attempts):
        _pace_ncbi(url)
        _pace_gnomad(url)
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
            rate_limited = exc.code == 429 or "rate limit" in body.lower()
            if (exc.code not in {429, 500, 502, 503, 504} and not rate_limited) or attempt == attempts - 1:
                raise RuntimeError(f"HTTP {exc.code} for {url}: {body[:400]}") from exc
            if rate_limited:
                time.sleep(20 * (attempt + 1))  # the server asked us to slow down
                continue
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = exc
            if attempt == attempts - 1:
                raise RuntimeError(f"request failed for {url}: {exc}") from exc
        time.sleep(2**attempt)
    raise RuntimeError(f"request failed for {url}: {last_error}")
