# OpenCRISPR-1 gRNA Designer v1.3.1 — architecture and scientific design

## 1. Scope and current evidence status

This application designs conservative OpenCRISPR-1-compatible guide candidates using a 20-nt spacer and NGG PAM, while deliberately avoiding an unsupported universal “OpenCRISPR efficacy score.”

The literature is **not treated as settled**. Ruffolo et al. (Nature, 2025) reported strong activity/specificity for OpenCRISPR-1. Tian et al. (Science Advances, 2025) reported lower on-target performance and more off-target activity than comparator systems in their evaluation. Hwang et al. (Genome Medicine, 2026) explicitly described these evaluations as conflicting, then reported Cas9-level on-target activity and substantially reduced off-target activity in their own head-to-head experiments. 2026 rice studies further support activity in plant systems. These differences mean OpenCRISPR-1 performance should be considered **context-dependent and still under independent validation**, not universally generalizable.

Accordingly, a software compatibility PASS means “meets the implemented sequence/PAM rules,” not “OpenCRISPR-1 will edit this locus efficiently or specifically.”

## 2. Modules

```text
app.py
  ├─ opencrispr_designer.py # NGG discovery, ranking, formats, legacy MIT/Hsu panel screen
  ├─ sequence_sources.py    # NCBI / Ensembl / manual FASTA + reproducibility metadata
  └─ validation.py          # PASS / REVIEW / FAIL checks

tests/
  ├─ scientific-core tests
  ├─ source/provenance tests
  ├─ validation tests
  └─ static Streamlit UI contract test
```

## 3. Target compatibility

The conservative default is:

- 20-nt targeting spacer;
- NGG PAM;
- both DNA strands scanned;
- independent sequence segments scanned separately.

This is a compatibility rule derived from published OpenCRISPR-1 targeting behavior and conventional SpCas9-guide compatibility. It is not a claim that all NGG guides have comparable activity.

## 4. Guide-expression formats

`guide_format_variants()` reports GX19, gX19, and gX20 configurations evaluated in the 2026 Genome Medicine study. The app does not assign a universal “best” format.

For a non-G-starting targeting spacer, gX19 changes position 1 and is explicitly flagged; gX20 retains the full 20-nt targeting sequence with an added 5′ G.

## 5. Sequence acquisition and reproducibility

`GeneSequenceRecord` stores:

- source and accession;
- source record/version;
- assembly or genomic accession when available;
- annotation release/date when available;
- UTC retrieval timestamp;
- SHA-256 fingerprint of the exact sequence segments used;
- ambiguity count/codes;
- warnings.

NCBI mode attempts to retain the versioned RefSeq accession plus linked chromosome/genomic accession version and record date. Ensembl mode records a versioned transcript ID when available, the returned assembly name, and current Ensembl REST release. Manual FASTA remains explicitly user-supplied/assembly-unspecified unless the user provides that context.

The SHA-256 fingerprint is the reproducibility anchor: if a public record changes later, the exported hash can reveal that the exact sequence input is no longer identical.

## 6. Ambiguous IUPAC bases

Ambiguous input is not silently treated as mismatch evidence:

1. non-ACGT IUPAC ambiguity codes are normalized to `N`;
2. candidate spacers containing `N` are skipped;
3. ambiguity may be accepted only at the degenerate `N` position of an otherwise resolved NGG/CCN PAM;
4. ambiguity in required G/G or C/C positions prevents PAM recognition;
5. provenance records ambiguity count and codes.

## 7. Candidate ranking

`sequence_quality_score()` is an explainable shortlist heuristic using GC, poly-T/homopolymer review terms and small sequence-context terms. It is **not** trained on OpenCRISPR-1 experimental outcomes and is never labelled as editing probability.

## 8. Local specificity: MIT/Hsu is a legacy baseline

`screen_local_reference()` implements the 2013 MIT/Hsu positional mismatch score because it is transparent and dependency-free. v1.3.1 now labels it everywhere as **“MIT/Hsu 2013 legacy baseline (supplied FASTA panel)”**.

This is intentionally not described as state-of-the-art. Independent benchmarking in CRISPOR found CFD more discriminative than MIT on validated off-target datasets (Haeussler et al., 2016), and later empirical/ML approaches can further improve prediction depending on data and context. The correct interpretation is therefore:

- useful transparent baseline for a small user-supplied panel;
- not a replacement for CFD or validated modern predictors;
- not genome-wide specificity;
- not OpenCRISPR-1-specific calibration.

The software does **not** fabricate a CFD value because a correct CFD implementation requires the published mismatch/PAM parameterization and careful backend validation. A future plugin can add CFD or a validated genome-aware backend without changing the core architecture.

## 9. Validation layer

Hard checks: 20-nt resolved spacer, NGG PAM, OpenCRISPR compatibility flag. Review checks: preferred GC, transparent quality threshold, poly-T/homopolymers, 5′-G format implications, and local specificity findings.

A local additional exact hit after excluding one presumed intended target is a hard supplied-panel failure; close hits or low legacy MIT/Hsu score are review flags.

## 10. Streamlit UI

The UI explicitly shows:

- a 2026 evidence-status warning about conflicting independent evaluations;
- sequence provenance and SHA-256 fingerprint;
- ambiguity warnings;
- guide validation;
- GX19/gX19/gX20 formats;
- legacy MIT/Hsu specificity language;
- CSV/FASTA/JSON exports.

JSON exports include provenance, evidence-status wording, and the specificity-method limitation.

## 11. Testing status

v1.3.1 passes **34 automated tests** with **88% branch-aware core coverage** covering both-strand NGG discovery, segment isolation, ambiguity rules, guide formats, MIT/Hsu calculations, local panel behavior, provenance fingerprints, validation, and export labels.

`tests/test_ui_contract.py` statically parses `app.py` and verifies that the evidence warning, provenance, legacy-specificity label, and validator are present. **This is not an interactive Streamlit browser test.** Streamlit is not installed in the artifact environment and PyPI access is unavailable, so the interface could not be launched here. Run the local `UI_TEST_CHECKLIST.md` before calling the UI release-tested.

## 12. Scientific boundaries

- OpenCRISPR-1 generalizability remains an empirical question across loci, cell types, delivery systems, guide formats, and organisms.
- NGG compatibility is not efficacy prediction.
- MIT/Hsu is a legacy panel baseline, not modern genome-wide specificity.
- Exact assembly/genotype and variant context must be verified before experimental use.
- Computational validation does not replace experimental off-target measurement.

## 13. Key references

- Ruffolo JA et al. *Nature*. 2025. doi:10.1038/s41586-025-09298-z.
- Tian R et al. *Science Advances*. 2025;11:eadu7334. doi:10.1126/sciadv.adu7334.
- Hwang H-Y et al. *Genome Medicine*. 2026;18:109. doi:10.1186/s13073-026-01682-2.
- Doench JG et al. *Nature Biotechnology*. 2016;34:184-191. doi:10.1038/nbt.3437.
- Haeussler M et al. *Genome Biology*. 2016. doi:10.1186/s13059-016-1012-2.
- Hsu PD et al. *Nature Biotechnology*. 2013;31:827-832. doi:10.1038/nbt.2647.

## 14. External comparison protocol

A fresh CHOPCHOP concordance run is intentionally not claimed from this artifact environment. `benchmarks/EXTERNAL_COMPARISON_PROTOCOL.md` requires the comparator to use the same frozen input sequence/assembly, 20-nt + NGG model, target segments and filtering scope. `benchmarks/concordance.py` then reports guides found by both, this project only, and the external tool only. Differences must be attributed to PAM handling, strand/segment boundaries, ranking/filtering, or database/version differences rather than described as unexplained superiority.

## 15. CI, packaging and release reproducibility

- Runtime/developer dependencies are exactly pinned.
- GitHub Actions tests Python 3.11, 3.12 and 3.13 on pushes/pull requests and enforces at least 80% core coverage.
- `coverage.txt` and `coverage.xml` preserve the measured artifact result.
- `pyproject.toml` installs the scientific core; Docker/Compose provides a one-command app environment.
- `CITATION.cff` and `.zenodo.json` are prepared for a later tagged archive; no DOI is claimed before it is minted.
- `OPENCRISPR_TERMS_NOTICE.md` separates this repository's MIT-licensed code from Profluent's third-party OpenCRISPR-1 technology/terms.
