"""Clickable demo page. `uv run snp-brief-demo` serves it with Streamlit."""

from __future__ import annotations

import sys
from pathlib import Path


def main() -> int:
    from streamlit.web import cli

    page = Path(__file__).with_name("app.py")
    sys.argv = [
        "streamlit",
        "run",
        str(page),
        "--server.showEmailPrompt=false",
        "--browser.gatherUsageStats=false",
        "--client.toolbarMode=minimal",
        *sys.argv[1:],
    ]
    return cli.main()
