# ============================================================================
# TEST UI CONTRACT
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: This file contains automated checks for the OpenCRISPR-1 gRNA
# Designer. The tests help make sure scientific and interface behaviour does
# not change accidentally when the program is edited.
#
# HOW TO READ A TEST:
# 1. Prepare sample input or a reusable fixture.
# 2. Call the function or workflow being tested.
# 3. Use assertions to compare the actual result with the expected result.
# 4. Existing test logic and expected scientific outcomes are unchanged.
#
# MAIN TOP-LEVEL TESTS / HELPERS:
# - function: test_streamlit_ui_source_contract
# ============================================================================

from pathlib import Path
import ast


def test_streamlit_ui_source_contract():
    app = Path(__file__).resolve().parents[1] / "app.py"
    text = app.read_text(encoding="utf-8")
    ast.parse(text)
    required = [
        "Evidence status (2026)",
        "conflicting",
        "Sequence provenance and reproducibility",
        "MIT/Hsu legacy specificity",
        "sequence SHA-256",
        "Quick validation of an existing OpenCRISPR guide",
        "Provenance field",
    ]
    for marker in required:
        assert marker in text
    assert "st.json(" not in text
