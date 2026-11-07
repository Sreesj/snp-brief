"""OpenAI-compatible chat client for the summary writer. Standard library only.

Works with any server that speaks the OpenAI chat-completions API:
- vLLM or NVIDIA NIM on a Brev / Nebius GPU box: SNP_BRIEF_LLM_URL=http://<host>:8000/v1
- NVIDIA API catalog: SNP_BRIEF_LLM_URL=https://integrate.api.nvidia.com/v1
- Ollama on a laptop: SNP_BRIEF_LLM_URL=http://localhost:11434/v1

Environment:
  SNP_BRIEF_LLM_URL    base URL ending in /v1 (required to enable the LLM)
  SNP_BRIEF_LLM_MODEL  model name the server expects (required)
  SNP_BRIEF_LLM_KEY    bearer token, if the server needs one

The model only sees the numbered facts. Its sentences are checked by
summary.check before anything is shown.
"""

from __future__ import annotations

import json
import os
import re
import urllib.request
from dataclasses import dataclass

from snp_brief.summary import Fact, facts_block

SYSTEM_PROMPT = """You write short plain-English summaries of a single genetic variant (SNP) for scientists.
Rules:
- Use ONLY the numbered facts provided. Do not add knowledge from memory.
- Write 3 or 4 sentences: what the variant does (function), how common it is (frequency, including ancestry differences if given), and its disease relevance per ClinVar. Mention a paper only if facts list one.
- End every sentence with the ids of the facts it uses, in square brackets, e.g. [F2] or [F3, F4].
- Copy numbers exactly as they appear in the facts you cite. Do not compute, round, or convert numbers.
- If a fact says sources disagree or a field is unresolved, say so plainly. Do not pick a side.
- Do not diagnose, give medical advice, or address a patient.
- Output only the sentences, no heading, no list."""


@dataclass
class LLMClient:
    base_url: str
    model: str
    api_key: str | None = None
    timeout: float = 60

    @classmethod
    def from_env(cls) -> "LLMClient | None":
        url = os.environ.get("SNP_BRIEF_LLM_URL")
        model = os.environ.get("SNP_BRIEF_LLM_MODEL")
        if not url or not model:
            return None
        return cls(url.rstrip("/"), model, os.environ.get("SNP_BRIEF_LLM_KEY"))

    def chat(self, system: str, user: str) -> str:
        payload = {
            "model": self.model,
            "messages": [{"role": "system", "content": system}, {"role": "user", "content": user}],
            "temperature": 0,
            "max_tokens": 400,
        }
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        if self.api_key:
            headers["Authorization"] = f"Bearer {self.api_key}"
        request = urllib.request.Request(
            f"{self.base_url}/chat/completions", data=json.dumps(payload).encode(), headers=headers
        )
        with urllib.request.urlopen(request, timeout=self.timeout) as response:
            body = json.loads(response.read().decode())
        return body["choices"][0]["message"]["content"]


def write_summary(client, rsid: str, facts: list[Fact]) -> list[str]:
    user = f"Variant: {rsid}\n\nFacts:\n{facts_block(facts)}\n\nWrite the summary."
    text = client.chat(SYSTEM_PROMPT, user)
    return split_sentences(text)


def split_sentences(text: str) -> list[str]:
    """Split after a closing citation tag (and its period), which ends each sentence."""
    text = re.sub(r"<think>.*?</think>", "", text, flags=re.S).strip()
    parts = re.split(r"(?<=\])\.?\s+|(?<=\]\.)\s*", text)
    out = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if not part.endswith("."):
            part += "."
        out.append(part)
    return out
