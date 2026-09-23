# OpenCRISPR-1 gRNA Designer v1.4.0

OpenCRISPR-1 gRNA Designer is a Python and Streamlit application for finding **20-nucleotide guide spacers next to an NGG PAM**, reviewing their sequence quality, and checking near-matches in a user-supplied FASTA panel.

The program is designed as a **candidate-design and review tool**.

Its scores are **not OpenCRISPR-specific editing probabilities**, and a PASS result does not mean that a guide has been experimentally proven to work.

**Hosted app:** https://opencrispr1-grna-designer.onrender.com/

---

## What this tool does

The basic workflow is:

1. enter or retrieve a DNA sequence;
2. search both strands for 20 nt + NGG targets;
3. collect candidate spacers;
4. check simple sequence-quality features;
5. optionally compare guides against a supplied reference panel;
6. validate the result;
7. export the design with provenance and settings.

The goal is to make the guide-selection process easier to inspect and reproduce.

---

## Important note about OpenCRISPR-1

This repository supports guide design for use in an OpenCRISPR-1-related workflow, but it does **not** reproduce the protein-design or scaffold-generative models described in the founding OpenCRISPR research.

It also does not provide a learned OpenCRISPR-specific editing-efficiency predictor.

The sequence-ranking and local-specificity checks in this application should therefore be treated as separate design aids.

---

## Guide format

The program searches for:

\`\`\`text
20 nt spacer + NGG PAM
\`\`\`

on both DNA strands.

The app can also provide guide-expression formatting options such as an X20-style sequence where appropriate.

---

## Input options

The application supports several ways to provide a target.

### Gene lookup

A supported gene can be retrieved from a biological database.

When several transcript choices exist, one transcript is selected according to the program's documented rules.

Use a specific accession when you need tighter control over the exact isoform.

### Accession

Supported NCBI or Ensembl identifiers can be used to retrieve a specific biological record.

Requested accession versions are checked rather than silently ignored.

### Manual FASTA

You can provide your own FASTA sequence.

Each record is treated as a real sequence segment rather than being silently merged into an artificial sequence.

---

## Why the program is strict about sequence structure

CRISPR design can become misleading if sequence pieces are joined together in ways that do not exist in the genome.

For that reason, the software blocks or carefully handles situations such as:

- missing exon-boundary information;
- multiple ambiguous CDS annotations;
- short-exon fallback joining;
- duplicate FASTA record names;
- empty FASTA records.

Compound CDS parts are handled as separate biological pieces where required.

---

## Local specificity screen

The optional local screen compares guide candidates with a supplied reference panel.

Current scope:

- up to **2,000,000 supplied bases**;
- **0–4 substitution mismatches**;
- **NGG PAMs**;
- both strands;
- explicit intended-site exclusion.

The screen is intentionally limited.

It does **not** model:

- DNA/RNA bulges;
- alternate PAMs;
- chromatin accessibility;
- sample-specific variants;
- an entire human genome unless that reference is actually supplied;
- all possible biological off-target mechanisms.

The MIT/Hsu-style specificity score is therefore a **panel summary**, not an OpenCRISPR-1 risk probability.

---

## Intended target exclusion

The program does not simply remove the first exact match and assume it is the intended target.

An exact site is excluded from off-target counting only when the user provides a specific:

- reference record;
- position;
- strand;

and the sequence at that location can be verified.

This prevents a repeated sequence from being incorrectly treated as the intended locus.

---

## PASS, REVIEW and validation

A PASS result means that the candidate passed the checks implemented in this software.

It does **not** mean:

- guaranteed editing;
- guaranteed specificity;
- experimental validation;
- clinical suitability.

The program also clears or invalidates stored results when important settings change, such as:

- target sequence;
- reference panel;
- mismatch radius;
- intended-site exclusion settings.

This helps prevent old results from being shown under new settings.

---

## Sequence and result provenance

JSON exports can retain information such as:

- exact target/reference sequence fingerprints;
- accession or source information;
- selected search settings;
- mismatch counts;
- intended-site coordinates;
- truncation information;
- validation checks;
- output/guide format.

Keeping these details makes a later rerun easier to compare with the original analysis.

---

## Offline example

The application includes a built-in synthetic example.

To try it without database access:

1. open the application;
2. click **Load offline example**;
3. click **Design OpenCRISPR guides**;
4. inspect candidates from both strands.

The bundled sequence is a software demonstration.

It is **not an experimentally validated guide set**.

---

## Installation

Python 3.11 or newer is recommended.

Clone the repository:

\`\`\`bash
git clone https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer.git
cd OpenCRISPR1-gRNA-Designer
\`\`\`

Create a virtual environment:

\`\`\`bash
python -m venv .venv
\`\`\`

Activate it.

### Windows PowerShell

\`\`\`powershell
.venv\Scripts\Activate.ps1
\`\`\`

### macOS/Linux

\`\`\`bash
source .venv/bin/activate
\`\`\`

Install dependencies:

\`\`\`bash
python -m pip install -r requirements.txt
\`\`\`

Run the application:

\`\`\`bash
python -m streamlit run app.py
\`\`\`

Then open:

\`\`\`text
http://localhost:8501
\`\`\`

---

## Docker

Docker is also supported:

\`\`\`bash
docker compose up --build
\`\`\`

---

## Main project files

- \`app.py\` — Streamlit user interface.
- \`opencrispr_designer.py\` — main guide discovery and ranking logic.
- \`local_screening.py\` — supplied-panel near-match screening.
- \`sequence_sources.py\` — target sequence retrieval.
- \`accession_sources.py\` — accession handling and source checks.
- \`validation.py\` — input and guide validation.
- \`workflow.py\` — workflow helpers.
- \`benchmarks/\` — deterministic examples and audit scripts.
- \`tests/\` — automated tests.

---

## Running the tests

Install development requirements:

\`\`\`bash
python -m pip install -r requirements-dev.txt
\`\`\`

Run the test suite:

\`\`\`bash
python -m pytest -q --cov=. --cov-config=.coveragerc --cov-fail-under=80
\`\`\`

Run the audit examples:

\`\`\`bash
python benchmarks/run_audit_examples.py
\`\`\`

The repository also includes automated CI checks.

---

## Current audit evidence

The current repository documents:

- **79 passing tests**;
- **90.98% branch-aware core coverage** in the reported audit environment;
- **14/14 deterministic example checks passing**.

The numerical coverage excludes \`app.py\`, while separate Streamlit tests exercise real application reruns.

Database retrieval unit tests use mocked responses.

Live database access can also fail because of network timeouts or remote-server problems, so a failed live lookup is not automatically a software-design failure.

---

## Scientific limitations

This application currently does not provide:

- a complete sgRNA scaffold generator;
- a base-editing designer;
- a prime-editing designer;
- a learned editing-efficiency model;
- a complete genome-wide search engine;
- experimental activity validation.

For final experimental work, guide candidates should be evaluated with appropriate genome-aware specificity tools and biological validation.

---

## Documentation

More detailed material is available here:

- [Updated software guide and audit](docs/UPDATED_SOFTWARE_GUIDE.md)
- [Word software guide](docs/OpenCRISPR1_Updated_Software_Guide.docx)
- [User guide](USER_GUIDE.md)
- [Architecture](ARCHITECTURE.md)
- [Change log](CHANGELOG.md)
- [Executable audit results](benchmarks/results/audit_examples.json)

---

## Scientific references

- Ruffolo et al. (2025), *Nature*. DOI: 10.1038/s41586-025-09298-z
- Hwang et al. (2026), *Genome Medicine*. DOI: 10.1186/s13073-026-01682-2
- Tian et al. (2025), *Science Advances*. DOI: 10.1126/sciadv.adu7334
- Hsu et al. (2013), *Nature Biotechnology*. DOI: 10.1038/nbt.2647

Experimental findings in these papers are context-dependent.

---

## License and feedback

The original software in this repository is available under the **MIT License**.

OpenCRISPR system/model terms are separate from this repository's software license. See:

- [OPENCRISPR_TERMS_NOTICE.md](OPENCRISPR_TERMS_NOTICE.md)

To report a problem or suggestion:

https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer/issues/new?template=feedback.yml
