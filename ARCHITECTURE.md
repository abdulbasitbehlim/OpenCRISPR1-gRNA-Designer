# OpenCRISPR 1 Designer architecture

Version 1.4.0 retains the Python/Streamlit application and separates target discovery from local specificity screening.

| File | Responsibility |
| --- | --- |
| `app.py` | Input forms, offline demo, source warnings, guide inspection, coordinate confirmation and exports |
| `sequence_sources.py` | FASTA validation, NCBI/Ensembl gene lookup, independent annotation intervals and sequence provenance |
| `accession_sources.py` | Single accession lookup and requested-version verification |
| `opencrispr_designer.py` | Both-strand NGG discovery, transparent ranking, format display and legacy tuple wrapper |
| `local_screening.py` | Streaming NGG reference search, coordinate exclusion, complete counts, bounded display and panel fingerprint |
| `validation.py` | Sequence checks and metadata-aware PASS/REVIEW/FAIL rules |
| `workflow.py` | Design identity, stale-state invalidation and consistent export construction |

## Coordinate contract

Start/end are 1-based inclusive coordinates of the spacer within the supplied segment. Minus-strand spacers are reverse-complemented into synthesis orientation, but start is still the leftmost reference spacer base. An NGG motif is tested in the corresponding strand orientation. Chromosome mapping is not inferred from transcript coordinates.

## Screening contract

`screen_reference(guide, panel, max_mismatches=3, intended_locus=None, max_hits=250)` returns `PanelScreen`. The intended locus is a `TargetLocus(contig, start, strand)` referring to the panel FASTA. It is excluded only after an exact spacer and supported PAM are found there. An invalid coordinate raises an error; it never triggers first-match exclusion.

Every resolved NGG site on both strands is considered. The search supports 0–4 substitutions, with no indels or bulges. Counts and risk sums are accumulated across all accepted hits. A bounded heap retains only the display subset. The panel length is capped at 2 million bases and displayed hits at 5000. No sequence-quality ranking is computed for reference sites.

`PanelScreen` includes normalized reference SHA-256, base/record counts, ambiguous-base count, searchable-site count, intended coordinate/status, mismatch histogram, full hit count, retained hits, settings, timestamp and method. No searchable sites yields a null score. Ambiguity and unknown intended context cannot yield a panel validation PASS. The score remains an uncalibrated MIT/Hsu legacy summary.

## State and exports

The design key hashes the input fingerprint/version and candidate segment, coordinates, strand, spacer and PAM. New design submissions clear old results before retrieval. Panel settings have a separate identity; changing them clears every stored panel result. The same export helper supplies the visible table and downloadable data. Counts in exported JSON come from the full screen, never from the display list length.

## Retrieval boundaries

Duplicate and empty FASTA records are rejected. NCBI extraction treats each annotated exon and CDS part independently; it does not use a compound-feature bounding span. RNA records lacking exon boundaries and multi-CDS records are blocked. Ensembl returns one selected transcript's exons, which can include untranslated regions. All source intervals remain local coordinates. These checks do not replace genome mapping or exon coverage validation.

## Migration from v1.3.1

Automatic `exclude_one_exact=True` is removed. Use the structured API instead:

```python
from local_screening import TargetLocus, screen_reference
result = screen_reference(
    'ACGTACGTACGTACGTACGT',
    {'intended': 'ACGTACGTACGTACGTACGTAGG'},
    intended_locus=TargetLocus('intended', 1, '+'),
)
print(result.total_hits, result.specificity_score)
```

The legacy tuple wrapper is retained for callers that need only score and display hits. `len(hits)` is not the full hit count when display truncation occurs. Validation and the UI use the structured result.

## Verification boundary

Core numerical coverage excludes the Streamlit file; three AppTest tests exercise actual form submission, screening, invalidation and offline demo reruns. Retrieval tests use mocked public-service responses. Synthetic examples check correctness, not nuclease efficacy. Primary references and limitations are documented in `docs/UPDATED_SOFTWARE_GUIDE.md`.
