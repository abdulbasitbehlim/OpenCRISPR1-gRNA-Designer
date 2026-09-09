# OpenCRISPR-1 gRNA Designer — user guide

## Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Recommended workflow

1. Choose **Gene lookup** for an online gene record or **Manual sequence / FASTA** for a reviewed locus.
2. For gene lookup, enter the gene and organism and choose NCBI RefSeq or Ensembl REST.
3. Set the maximum guide count and minimum sequence-quality score in the sidebar.
4. Run the design and review provenance warnings before selecting a candidate.
5. Inspect the candidate table, target strand, segment, coordinates, GC and sequence-quality score.
6. For the selected spacer, review GX19/gX19/gX20 expression-format alternatives. Do not treat the displayed order as an efficacy ranking.
7. If you have an intended locus plus suspected near-match sequences, use the optional local FASTA specificity screen.
8. Export CSV, FASTA or JSON for downstream review.

## Interpreting the results

**OpenCRISPR-1 compatible** means the candidate follows the conservative 20-nt + NGG target architecture supported by the published OpenCRISPR-1 literature. It does not mean guaranteed editing.

**Sequence quality** is a transparent heuristic for prioritization. It is not an OpenCRISPR-specific probability or experimentally calibrated efficiency.

**MIT/Hsu specificity** is computed only against the FASTA sequences you supply. A high supplied-panel score does not prove whole-genome specificity.

## Guide formats

- GX19: natural first G in a 20-nt spacer.
- gX19: first base is changed to G; the position-1 targeting mismatch must be considered.
- gX20: a non-templated G is prepended to the intact 20-nt spacer.

OpenCRISPR-1 showed robustness to these formats in a 2026 study, but target- and system-specific behavior can still vary.

## Before experimental use

Confirm genome assembly, gene/transcript and coding context, allele/variant sequence, and genome-wide off-targets with an appropriate genome-aware workflow. Computational ranking does not replace experimental validation.

## Validation

Every generated guide receives an overall PASS, REVIEW or FAIL result. Select a guide to inspect each check and its explanation. The local specificity status remains **NOT SCREENED** until you provide a FASTA reference panel and run the MIT/Hsu screen.

Use **Quick validation of an existing OpenCRISPR guide** near the top of the app to check a 20-nt spacer plus its actual three-base PAM. This quick mode validates sequence/PAM compatibility and expression-format flags only; it cannot prove that the target exists at a genomic locus or that the guide is genome-wide specific.

After a local-reference screen, the validator updates the specificity status for that spacer and the updated result is included in exports.

## 2026 evidence-status note

Do not interpret an OpenCRISPR-1 compatibility PASS as proof that the nuclease will outperform SpCas9 at your locus. The literature now contains conflicting evaluations. Ruffolo et al. (2025) reported favorable performance; Tian et al. (2025) reported less favorable on/off-target behavior in their comparison; Hwang et al. (2026) explicitly discussed the discrepancy and then reported favorable head-to-head results under their own conditions. Benchmark your loci and delivery context experimentally.

## Specificity-screen interpretation

The built-in supplied-FASTA screen uses **MIT/Hsu (2013) only as a legacy baseline**. It is retained because it is transparent and dependency-free. It is not state-of-the-art: CFD and later empirical/ML approaches can perform better on validated off-target datasets.

Use the local screen to inspect a small supplied panel, not to claim genome-wide specificity. Before experimental use, use an appropriate genome-aware search and validated scoring/experimental workflow.

## Reproducibility and provenance

After sequence retrieval, expand **Sequence provenance and reproducibility** and retain:

- accession/source version;
- assembly or genomic accession;
- annotation release/date;
- UTC retrieval time;
- sequence SHA-256.

These values are also included in JSON export. The sequence hash allows you to detect a future public-record change even when the gene name is unchanged.

## Ambiguous IUPAC sequence

Ambiguity codes are normalized to `N`. Spacer windows containing ambiguity are skipped instead of being treated as ordinary mismatches. An unresolved nucleotide is allowed only at the degenerate N position of an otherwise resolved NGG/CCN PAM. Ambiguity count and codes are reported in provenance.

## UI testing note

The package includes a static UI contract test, but Streamlit could not be installed/launched in the build environment. Run `UI_TEST_CHECKLIST.md` locally before deployment or reviewer handoff.

## Reviewer/release workflow (v1.3.1)

Freeze the exact target sequence/assembly and its SHA-256 before comparing designers. Run CHOPCHOP (or another comparator) with the same 20-nt + NGG and region assumptions, export a `guide` column, and calculate concordance with `benchmarks/concordance.py`. Do not treat generic CHOPCHOP ranking as an OpenCRISPR-1 efficacy ground truth.

Automated core QA at artifact build: **34 tests, 88% branch-aware coverage**. Numerical coverage excludes `app.py`; complete `UI_TEST_CHECKLIST.md` interactively before release.
