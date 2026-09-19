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
