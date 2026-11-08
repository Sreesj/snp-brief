"""OpenAI-compatible chat client. Standard library only.

Nebius Token Factory is the default when NEBIUS_API_KEY is set:

    export NEBIUS_API_KEY=...                       # from tokenfactory.nebius.com
    uv run snp-brief --list-models                  # see which models your account can use
    export SNP_BRIEF_LLM_MODEL=openai/gpt-oss-120b  # optional; this is the default
    export SNP_BRIEF_FORGE_MODEL=...                # optional, a stronger model for Hypothesis Forge

Any other OpenAI-compatible server works too (vLLM or NVIDIA NIM on a GPU box,
Ollama on a laptop) with SNP_BRIEF_LLM_URL, SNP_BRIEF_LLM_MODEL and
SNP_BRIEF_LLM_KEY.

The model only sees numbered facts. Everything it writes is checked before it
is shown (summary.check, forge.validate).
"""

from __future__ import annotations

import json
import os
import re
import urllib.error
import urllib.request
from dataclasses import dataclass, replace

from snp_brief.summary import Fact, facts_block

NEBIUS_URL = "https://api.tokenfactory.nebius.com/v1"
NEBIUS_DEFAULT_MODEL = "openai/gpt-oss-120b"

SYSTEM_PROMPT = """You write short plain-English summaries of a single genetic variant (SNP) for scientists.
Rules:
- Use ONLY the numbered facts provided. Do not add knowledge from memory.
- Write 3 or 4 sentences: what the variant does (function), how common it is (frequency, including ancestry differences if given), and its disease relevance per ClinVar. Mention a paper only if facts list one.
- End every sentence with the ids of the facts it uses, in square brackets, e.g. [F2] or [F3, F4].
- Copy numbers exactly as they appear in the facts you cite, including the % sign. Do not compute, round, or convert numbers.
- Copy gene symbols and rsIDs exactly as written in the facts.
- Say which allele a frequency belongs to exactly as the fact does (the alternate allele unless the fact says otherwise).
- Do not say a variant is "in" a gene unless a fact says so; if a fact says it is upstream or downstream, say that.
- UniProt text describes the protein, never the variant. Do not imply the variant changes that function.
- Paper facts contain only the title. Do not describe a paper's design (meta-analysis, cohort, GWAS) or findings unless the title states it.
- Name the gnomAD call set (exomes or genomes) for every frequency you give.
- If a fact says sources disagree or a field is unresolved, say so plainly. Do not pick a side.
- Do not diagnose, give medical advice, or address a patient.
- Output only the sentences, no heading, no list."""


@dataclass
class LLMClient:
    base_url: str
    model: str
    api_key: str | None = None
    timeout: float = 120

    @classmethod
    def from_env(cls) -> "LLMClient | None":
        nebius_key = os.environ.get("NEBIUS_API_KEY")
        url = os.environ.get("SNP_BRIEF_LLM_URL") or (NEBIUS_URL if nebius_key else None)
        if not url:
            return None
        model = os.environ.get("SNP_BRIEF_LLM_MODEL") or (NEBIUS_DEFAULT_MODEL if url.rstrip("/") == NEBIUS_URL else None)
        if not model:
            return None
        key = os.environ.get("SNP_BRIEF_LLM_KEY") or nebius_key
        return cls(url.rstrip("/"), model, key)

    @property
    def provider(self) -> str:
        return "Nebius Token Factory" if self.base_url == NEBIUS_URL else self.base_url

    def for_forge(self) -> "LLMClient":
        """The client Hypothesis Forge uses: SNP_BRIEF_FORGE_MODEL if set, else the same model."""
        forge_model = os.environ.get("SNP_BRIEF_FORGE_MODEL")
        return replace(self, model=forge_model) if forge_model else self

    def chat(self, system: str, user: str, *, json_mode: bool = False, max_tokens: int = 2000) -> str:
        payload: dict = {
            "model": self.model,
            "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}],
            "temperature": 0,
            "max_tokens": max_tokens,
        }
        if json_mode:
            payload["response_format"] = {"type": "json_object"}
        body = self._request("/chat/completions", payload)
        message = body["choices"][0]["message"]
        if message.get("refusal"):
            raise RuntimeError(f"model refused: {message['refusal']}")
        content = message.get("content") or ""
        if not content.strip():
            raise RuntimeError("model returned no content (try a larger max_tokens or another model)")
        return content

    def list_models(self) -> list[str]:
        body = self._request("/models", None)
        return sorted(str(item.get("id")) for item in body.get("data", []) if item.get("id"))

    def _request(self, path: str, payload: dict | None) -> dict:
        headers = {"Accept": "application/json", "User-Agent": "snp-brief/0.2"}
        data = None
        if payload is not None:
            headers["Content-Type"] = "application/json"
            data = json.dumps(payload).encode()
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        request = urllib.request.Request(f"{self.base_url}{path}", data=data, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=self.timeout) as response:
                return json.loads(response.read().decode())
        except urllib.error.HTTPError as exc:
            detail = exc.read().decode(errors="replace")[:300]
            raise RuntimeError(f"HTTP {exc.code} from {self.provider}: {detail}") from exc


def write_summary(client, rsid: str, facts: list[Fact]) -> list[str]:
    user = f"Variant: {rsid}\n\nFacts:\n{facts_block(facts)}\n\nWrite the summary."
    text = client.chat(SYSTEM_PROMPT, user, max_tokens=2000)
    return split_sentences(text)


def split_sentences(text: str) -> list[str]:
    """Split after a closing citation tag (and its period), which ends each sentence."""
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.S).strip()
    parts = re.split(r"(?<=\])\.?\s+|(?<=\]\.)\s*", text)
    out: list[str] = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if not part.endswith("."):
            part += "."
        if out and part[:1].islower() and not part.startswith(("rs", "c.", "p.")):
            # "... [F1]. and is recorded ..." continues the previous sentence.
            out[-1] = out[-1].rstrip(".") + " " + part
            continue
        out.append(part)
    return out


def parse_json(text: str) -> dict:
    """The first JSON object in a reply, tolerating code fences or a thinking block."""
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.S).strip()
    fenced = re.search(r"```(?:json)?\s*(\{.*\})\s*```", text, flags=re.S)
    if fenced:
        text = fenced.group(1)
    start, end = text.find("{"), text.rfind("}")
    if start < 0 or end < start:
        raise ValueError("no JSON object in the reply")
    return json.loads(text[start : end + 1])
