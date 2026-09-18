# Change log

## 1.4.0

Fix unverified exact-hit exclusion, capped-score bias, duplicate FASTA data loss, spliced-transcript fallback and stale UI screening results. Require explicit intended coordinates, preserve full counts, validate accession versions, separate compound annotation intervals, and export reference/settings provenance. Add a bounded local scanner, offline examples and Streamlit workflow tests. Include accession, screening and workflow modules in package metadata.

Behavior changes: missing transcript exon annotations and multiple CDS selections now stop for explicit input correction; `exclude_one_exact=True` is no longer accepted. This version does not add a trained activity predictor or genome-wide search.

Verification: 79 tests pass, 90.98% branch-aware core coverage; 14 deterministic examples pass. Live Ensembl retrieval was unsuccessful during this audit. See the guide for exact scope.
