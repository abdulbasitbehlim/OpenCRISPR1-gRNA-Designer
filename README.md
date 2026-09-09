# OpenCRISPR-1 gRNA Designer v1.3.1

Research software for designing and validating conservative NGG-compatible guide RNAs for OpenCRISPR-1, with tests, provenance, validation and reproducible scientific documentation.

## Live app and related projects

**Run the OpenCRISPR-1 app:** https://opencrispr1-grna-designer.onrender.com/

**Give feedback:** https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer/issues/new?template=feedback.yml

**GitHub repository:** https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer

**Also available — Plant MultiGene gRNA Designer:**
- Live app: https://plant-multigene-grna-designer.onrender.com/
- GitHub: https://github.com/abdulbasitbehlim/Plant-MultiGene-gRNA-Designer

> The Render free tier may sleep after inactivity, so the first load can take longer while the service wakes up.

## Scope

This tool is a conservative OpenCRISPR-1 design workbench. It scans independent sequence segments on both strands for 20-nt spacers adjacent to an NGG PAM, ranks transparent sequence-quality features, validates each candidate, and reports guide-expression formats used in OpenCRISPR-1 studies.

The software deliberately separates:

- **OpenCRISPR compatibility** — 20-nt spacer + NGG design rule;
- **sequence quality** — GC, poly-T, homopolymers and transparent heuristic ranking;
- **specificity** — optional supplied-FASTA near-match screening using a legacy MIT/Hsu baseline.

It does not claim that the heuristic score is an OpenCRISPR efficacy probability, and it does not replace genome-wide off-target assessment or experimental validation.

## Key features

- gene lookup or reviewed manual FASTA input;
- independent exon/segment scanning to avoid synthetic junction targets;
- both-strand 20-nt + NGG discovery;
- GX19, gX19 and gX20 guide-expression format reporting;
- PASS / REVIEW / FAIL validation;
- optional local FASTA near-match specificity screen;
- custom existing-guide validation;
- accession/version, assembly/release, retrieval time and SHA-256 provenance;
- CSV, FASTA and JSON exports;
- automated tests and GitHub Actions CI.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
streamlit run app.py
```

Or with Docker:

```bash
docker compose up --build
```

Then open `http://localhost:8501`.

## Validation and testing

The validated v1.3.1 package passed **34/34 tests** with approximately **88% branch-aware scientific-core coverage**. The test suite covers guide discovery, both strands, ambiguity handling, exon/segment boundaries, provenance, validation, guide-format variants and local specificity behavior.

A PASS means the candidate satisfies the rules implemented here. It does **not** mean experimentally proven editing success.

## Evidence boundary

OpenCRISPR-1 has strong positive founding and follow-up reports, but later independent evaluations have not been uniformly consistent across targets and experimental settings. This project therefore treats OpenCRISPR-1 compatibility as a design rule rather than a guarantee of universal performance.

The optional MIT/Hsu specificity calculation is retained as a transparent legacy baseline only. It is not presented as state-of-the-art genome-wide off-target prediction.

## Licensing note

The original software in this repository is MIT licensed. OpenCRISPR/OpenCRISPR-1 system/model terms are separate and are not relicensed by this repository. See `OPENCRISPR_TERMS_NOTICE.md`.
