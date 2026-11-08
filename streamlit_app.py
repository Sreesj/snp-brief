"""Entry point for Streamlit Community Cloud: `streamlit run streamlit_app.py`."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from snp_brief.demo.app import main  # noqa: E402

main()
