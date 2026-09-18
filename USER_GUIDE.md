# OpenCRISPR 1 gRNA Designer user guide

Version 1.4.0. See [the detailed guide](docs/UPDATED_SOFTWARE_GUIDE.md) for the evidence review, audit examples and remaining limitations.

## Try it without a database

Run `python -m streamlit run app.py`, click **Load offline example**, then **Design OpenCRISPR guides**. The demo includes one forward and one reverse target. These synthetic sequences check software behavior only.

## Design your candidates

1. Choose **Gene lookup**, **Accession ID**, or **Manual sequence / FASTA**.
2. Select the correct organism and database. Ensembl gene lookup retrieves independent exons. Gene lookup selects one transcript, so choose a stable transcript accession if isoform identity matters.
3. In manual mode, paste independently reviewed contiguous genomic regions or individual exon records. Do not paste spliced cDNA as though it were contiguous genomic DNA.
4. Set maximum guides and minimum sequence-quality score, then submit.
5. Inspect the returned source version, assembly information, warnings and sequence SHA-256.
6. Review spacer, PAM, strand, segment, start and end. Coordinates are local to the segment, 1-based and inclusive; they are not chromosomal coordinates.
7. Review expression formats. GX19 has a native first G, gX19 replaces base 1, gX20 adds a G before the intact spacer, and X20 keeps a non-G-starting spacer unchanged. These are spacer formats, not complete sgRNAs or ranked efficacy recommendations.

NCBI RNA records without exon annotations now stop with a message directing you to annotated exons or genomic FASTA. A record with several CDS annotations is rejected; select a single gene explicitly. Requested versions that disagree with returned records are rejected.

## Screen a small reference panel

Paste the contents of `examples/panel_demo.fasta` into **Reference FASTA** and select the ACGT-repeat demo spacer. Leave mismatch radius at 3. Without specifying an intended locus, all exact matches remain and the result requires review.

To identify the demo intended site, enable **Exclude a verified intended locus**. Enter FASTA ID `intended`, start `1`, strand `+`. Run the screen. It verifies the spacer and an NGG PAM at those coordinates, excludes only that site, and retains the duplicate exact match and the one-substitution match.

For a reverse-strand hit, start means the leftmost spacer base on the supplied reference, not the leftmost PAM base and not the spacer's 5-prime base. Check the strand carefully. The screen searches the original 20-base targeting spacer; it does not quantify how 5-prime G expression changes affect activity.

The panel limit is 2 million bases; display limits range from 1 to 5000 hits. Total counts and the score use every match within the radius, even when fewer rows are displayed. All-N references produce no score. Ambiguous reference bases, a missing intended coordinate or a radius below 3 require review. A small clean panel can still miss real genomic off-target sites.

## Understand and retain results

PASS means the implemented sequence/panel rules passed. REVIEW means an issue or unresolved context needs attention. FAIL indicates an unsupported sequence/PAM or an additional exact panel match after the intended site is verified. Specificity remains NOT SCREENED without a completed panel search.

Download CSV for tables, FASTA for targeting spacers and JSON for the full audit details. Keep the original target and reference FASTA with the JSON. A sequence hash identifies the input but does not reconstruct the sequence or prove biological correctness.

Changing panel text, mismatch radius, intended coordinate or display limit clears stored panel results. A new design clears previous screens. Changing sidebar design settings requires resubmitting before results and exports are shown.

## API migration

`screen_local_reference` no longer accepts `exclude_one_exact=True`; this raises an actionable error. Use `local_screening.screen_reference` with an explicit `TargetLocus` to obtain complete result metadata. The old tuple wrapper returns only the score and retained display hits; its list length is not the total number of hits when capped.

## Boundaries

The sequence score is a heuristic, not an editing probability. The panel search is NGG-only, substitution-only and local. Genomic mapping, variant checking, genome-wide search, scaffold choice, editing-outcome prediction, delivery and experimental confirmation remain separate work. The app does not run the founding paper's AI protein or guide-scaffold generator.
