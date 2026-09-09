from pathlib import Path
import ast


def test_streamlit_ui_source_contract():
    app = Path(__file__).resolve().parents[1] / "app.py"
    text = app.read_text(encoding="utf-8")
    ast.parse(text)
    required = [
        "Sequence provenance, reproducibility and warnings",
        "Mismatch-aware search strategy and computational workload",
        "Validate a custom 20-nt guide against these genes",
        "Recommended interactive range",
        "sequence SHA-256",
    ]
    for marker in required:
        assert marker in text
