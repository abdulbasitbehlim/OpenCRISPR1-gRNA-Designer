# GitHub Deployment Guide — OpenCRISPR-1 gRNA Designer

## Recommended repository

**Repository name:** `OpenCRISPR-gRNA-Designer`

**Description:** Research-oriented Python/Streamlit workbench for conservative OpenCRISPR-1-compatible 20-nt + NGG guide discovery, validation, provenance, guide-format review and local specificity screening.

## Before making the repository public

1. Run `python -m pip install -r requirements-dev.txt` in a clean Python 3.11–3.13 environment.
2. Run `pytest -q --cov=. --cov-config=.coveragerc --cov-report=term-missing --cov-fail-under=80`.
3. Run `streamlit run app.py` and complete `UI_TEST_CHECKLIST.md`.
4. Test a manual sequence, a live NCBI/Ensembl lookup, GX19/gX19/gX20 display, and the custom-guide validator.
5. Test the supplied-FASTA specificity workflow and verify exports.
6. Keep `OPENCRISPR_TERMS_NOTICE.md` and the README licensing boundary intact.

## Upload to GitHub

Create an empty repository and upload **the contents of this folder at repository root**. Do not create an extra nested project folder inside the repository.

The included `.github/workflows/ci.yml` runs the test/coverage matrix on pushes and pull requests. No GitHub secret is required for the normal test suite.

## Streamlit Community Cloud

After the repository is public, create an app from the repository and set the entrypoint to:

`app.py`

If you want to identify your NCBI requests with your own contact email, add `NCBI_EMAIL` as an app secret/environment value where supported. Never commit private credentials.

## Docker

```bash
docker compose up --build
```

Then open `http://localhost:8501`.

## First release

For the first public research-preview release, a conservative tag is `v0.1.0`. Do not claim a Zenodo DOI until the exact public tagged release has actually been archived.

## OpenCRISPR-1 boundary

The repository's MIT license covers this project's original software code only. It does not grant rights to OpenCRISPR-1 protein/scaffold/patents/trademarks or other Profluent materials. Keep `OPENCRISPR_TERMS_NOTICE.md` visible and verify current upstream terms before experimental or commercial use.

## Scientific boundary

The app reports conservative sequence compatibility and transparent validation. It does not provide a validated OpenCRISPR-1 efficacy probability, and MIT/Hsu local-panel specificity is a legacy baseline rather than a state-of-the-art whole-genome guarantee.
