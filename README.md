# OpenCRISPR-1 gRNA Designer v1.4.0

A Python/Streamlit workbench for discovering conservative 20-nt + NGG targeting spacers, reviewing sequence quality, and checking near matches in a supplied FASTA panel. Scores are **not OpenCRISPR-specific editing probabilities**.

## Run locally

Python 3.11 or newer:

```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS / Linux instead: source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open `http://localhost:8501`. Click **Load offline example**, then **Design OpenCRISPR guides**, to try both strands without database access. The bundled sequences are synthetic software examples, not experimentally validated guides.

Docker is also supported with `docker compose up --build`.

The existing hosted endpoint is https://opencrispr1-grna-designer.onrender.com/ . A source update on a development branch does not establish that this hosted instance is running v1.4.0.

## What changed

- An exact match is excluded only at a user-specified, sequence-verified contig/start/strand. No first match is automatically assumed to be the intended locus.
- Scores and total counts include **all matching sites within the mismatch radius**, even if the hit display is capped.
- Duplicate or empty FASTA records are rejected instead of silently losing data.
- NCBI RNA records without exon boundaries, records with multiple CDS annotations, and short-exon fallback joins are blocked. Compound CDS parts are handled independently.
- Requested accession versions are checked; Ensembl version suffixes are not silently ignored.
- Changing target, panel, mismatch radius or exclusion settings invalidates incompatible stored results. New failed designs clear old results.
- JSON exports include reference hashes, search parameters, total mismatch counts, hit truncation, verified intended coordinates, expression formats and validation checks.
- Bounded panel scanning, explicit input limits, a native X20 format option, source links and executable examples improve practical use.

## Scope and limits

The app finds 20-base spacers next to NGG on either strand. Gene lookup selects one transcript; use a specific accession to control the isoform. Manual FASTA records are independent sequence segments. Coordinates are 1-based inclusive and relative to each input segment, not chromosome positions.

The optional local screen supports up to 2,000,000 supplied bases, 0–4 substitutions, and NGG PAMs. It does not cover bulges, alternate PAMs, chromatin, sample variants or an entire human genome. The legacy MIT/Hsu score is a panel summary, not a calibrated OpenCRISPR-1 risk estimate. **PASS means implemented rules passed, not experimentally established editing success.**

No full sgRNA scaffold generator, base-editing designer, prime-editing designer, learned efficacy model or whole-genome search is included. Retain the original reference FASTA with exported JSON.

## Verification

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q --cov=. --cov-config=.coveragerc --cov-fail-under=80
python benchmarks/run_audit_examples.py
```

Audit on Python 3.12: **79 tests passed**, **90.98% branch-aware core coverage**, and **14/14 deterministic examples passed**. Coverage excludes `app.py`; three Streamlit AppTest tests exercise its actual reruns separately. Database retrieval unit tests use mocked responses. Live Ensembl attempts for HBB, TP53 and VEGFA failed in the audit environment (timeouts/server error); successful live retrieval is not claimed.

## Documentation

- [Updated software guide and audit](docs/UPDATED_SOFTWARE_GUIDE.md)
- [Download the updated Word guide](docs/OpenCRISPR1_Updated_Software_Guide.docx)
- [User guide](USER_GUIDE.md)
- [Architecture and API migration](ARCHITECTURE.md)
- [Change log](CHANGELOG.md)
- [Executable example results](benchmarks/results/audit_examples.json)

## Scientific sources

- [Ruffolo et al., Nature, 2025](https://doi.org/10.1038/s41586-025-09298-z)
- [Hwang et al., Genome Medicine, 2026](https://doi.org/10.1186/s13073-026-01682-2)
- [Tian et al., Science Advances, 2025](https://doi.org/10.1126/sciadv.adu7334)
- [Hsu et al., Nature Biotechnology, 2013](https://doi.org/10.1038/nbt.2647)

Experimental findings are context-dependent. This app is a candidate-design aid and does not reproduce the protein/scaffold generative models in the founding paper.

## Licensing and feedback

Original repository software is MIT licensed. OpenCRISPR system/model terms remain separate; see [the terms notice](OPENCRISPR_TERMS_NOTICE.md).

[Report an issue](https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer/issues/new?template=feedback.yml).
