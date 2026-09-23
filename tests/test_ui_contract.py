# ============================================================================
# TEST UI CONTRACT
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Contains automated tests that protect the Plant MultiGene gRNA
# Designer from accidental behaviour changes.
#
# HOW TO READ THIS FILE:
# 1. A test prepares sample input.
# 2. It calls the function being checked.
# 3. Assertions compare the result with the expected behaviour.
# 4. Test expectations are kept unchanged; only explanatory structure is added.
#
# MAIN TOP-LEVEL PARTS:
# - function: test_streamlit_ui_source_contract
# ============================================================================

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
        "Saved design settings",
        "Provenance field",
    ]
    for marker in required:
        assert marker in text
    assert "st.json(" not in text
