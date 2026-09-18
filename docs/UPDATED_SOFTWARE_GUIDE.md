# OpenCRISPR 1 Guide RNA Designer Updated Software Guide

Software version 1.4.0 | Review date 18 September 2026

## 1 Purpose and main conclusion

The updated software is a more reliable workbench for finding candidate OpenCRISPR-1 targeting spacers and reviewing near matches in a small reference panel. The update corrects errors that could hide an exact sequence match, inflate a specificity score, lose FASTA records, or carry an old screening result into a new analysis.

The practical conclusion is that the software is suitable for sequence-level candidate discovery and preliminary panel review. It cannot establish that a guide will edit efficiently, that it is specific across a whole genome, or that a complete guide RNA construct is ready for use. Those questions need additional data and validation.

This guide explains the paper evidence, the limitations found in the original code, the implemented changes, worked examples, installation, daily use, output interpretation and remaining development needs. It is written for a beginner who needs to explain the software to a supervisor or team.

The audit began from repository commit 2f9e45d2243de02d65e4ff200c21d1edf01f3c51, version 1.3.1. The updated source is version 1.4.0. The repository is https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer . The development branch is improve/validated-local-screening. A branch update does not by itself update an existing hosted deployment.

## 2 Understand the software in simple terms

Think of a targeting spacer as a short address. The software searches your supplied DNA for possible addresses that fit its supported rule. It then checks simple sequence properties and, if you provide another FASTA file, looks for similar addresses there. It does not observe an editing experiment.

| Term | Meaning in this software |
| --- | --- |
| Target sequence | The DNA region in which you want to find candidates |
| Spacer | The 20 bases that specify a target sequence |
| PAM | A nearby three-base DNA motif; this app supports NGG |
| N in NGG | Any base in the first position; the two G positions must be resolved |
| Scaffold | The RNA structure that associates with the nuclease; not generated here |
| Reference panel | The sequences you supply for a local near-match search |
| Off-target candidate | Another matching or similar sequence; not proof of actual cleavage |
| SHA-256 | A fingerprint for detecting a change in the saved input |

The app first retrieves or reads sequence records, keeps segments separate, searches both DNA strands, calculates a transparent sequence-quality score and displays candidates. You may then screen one candidate against a reference panel and export the results.

For example, the synthetic record ACGTACGTACGTACGTACGTAGG contains the 20-base spacer ACGTACGTACGTACGTACGT followed by AGG. The spacer occupies positions 1 through 20 in that record. The PAM is separate and should not be appended to the targeting spacer when interpreting the FASTA export.

The output uses the DNA alphabet, including T, to represent spacer sequences. A displayed spacer or 5-prime format is not a complete sgRNA, an RNA order specification, or a cloning instruction.

## 3 What the papers establish

### The founding study

Ruffolo and colleagues generated CRISPR proteins and compatible guide-scaffold components using language models. Their workflow differs from this app, which searches a user-supplied target for spacers. In their characterization, OpenCRISPR-1 showed comparable activity to SpCas9 at NGG sites and reduced activity at non-NGG sites. The expanded PAM comparison retained 49 NGG and 43 non-NGG targets after quality control. Most experiments used the SpCas9 scaffold; OpenCRISPR-1 also functioned with a designed scaffold. These results support a conservative NGG candidate search, but do not validate this app's hand-written ranking score [1].

### Independent guide format evaluation

Hwang and colleagues compared OpenCRISPR-1 with Cas9 across 28 endogenous loci in HEK293T cells, with further work in iPSCs and fibroblasts. They reported favorable specificity and tested GX19, gX19 and gX20 formats. They also discussed conflicting earlier evaluations. The software therefore displays these formats without claiming a universal winner or transferring their reported editing rates to new sequences [2].

### A different comparative result

Tian and colleagues compared FrCas9, SpCas9 and OpenCRISPR-1 using several experimental assays. Their study favored FrCas9 in the evaluated settings. Together with the other reports, it shows why a single superiority claim is inadequate for software that accepts arbitrary user targets [3].

### Plant evidence

Das and colleagues reported OpenCRISPR-1-based knockout, base editing and prime editing in rice. This supports investigation in plant systems. It does not make this spacer-design app a base-editor or prime-editor design program [4].

The implementation decision is to keep three outputs separate: a supported sequence rule, an uncalibrated quality ranking and a local panel summary. A result from any one of these cannot substitute for the other two or for an experiment. The literature summaries above are published findings; the examples later in this guide are software tests performed for this update.

## 4 Limitations found and changes made

| Original limitation | Why it mattered | Updated behavior |
| --- | --- | --- |
| First exact hit automatically removed | A real match could disappear without proof it was the intended site | Exclusion requires a verified contig, start and strand |
| Hit display capped before scoring | A smaller display could produce a higher score | All accepted hits contribute to counts and score |
| Duplicate FASTA names overwritten | Earlier records disappeared silently | Duplicate identifiers cause an input error |
| Empty FASTA records accepted | Missing sequence could go unnoticed | Empty records are rejected |
| Spliced transcript fallback | A candidate could span a junction absent from genomic DNA | RNA records without exon boundaries are blocked |
| Short exon fallback and compound spans | Separate regions could be treated as one interval | Exons and CDS parts remain independent |
| First CDS chosen automatically | A multi-gene record could select the wrong gene | Multiple CDS annotations require explicit input selection |
| Ensembl version suffix discarded | A different record version could be used | Requested and returned versions must agree |
| Screens stored only by spacer | Results could survive a target or reference change | Design identity and panel settings control reuse |
| Limited screening provenance | A score could not be traced to its reference and settings | JSON includes hashes, coordinates, settings and counts |
| Reference scanning ranked every site | Unnecessary work and memory use | Streaming search with a bounded display heap |
| Accession module absent from package list | Installed core package could miss functionality | Required core modules are included |

These changes fix concrete software behavior. They do not remove every scientific limitation. In particular, the app still does not contain an experimentally trained OpenCRISPR-1 activity predictor or a whole-genome off-target engine.

## 5 Worked examples of corrected screening

All examples in this section use deliberately constructed sequences. They check whether the program follows its stated rules. They do not measure editing activity in cells.

### Example 1 An exact match without a confirmed intended target

Supply the synthetic spacer and a panel containing one identical NGG-adjacent sequence named other. The old default removed this match and returned a score of 100 with no hits. The updated program retains the hit, returns a legacy score of 50 and marks the context for review. The software no longer assumes that the first sequence match is the correct biological locus.

### Example 2 A verified intended target

Use a panel record named intended containing ACGTACGTACGTACGTACGTAGG. Specify intended, start 1, plus strand. The program verifies the 20-base spacer and supported PAM at that location, then excludes only that site. The remaining hit count is zero and the panel score is 100. This means no retained matches in that tiny panel; it says nothing about the rest of a genome.

### Example 3 An additional identical locus

Add a second record with the same sequence. After excluding the verified intended site, one exact match remains. Its count is one and the legacy score is 50. The supplied-panel validation becomes FAIL. A repeated target cannot be made unique by choosing the first matching record.

### Example 4 The display limit must not change the answer

Construct 301 independent records containing the same target. After identifying one intended site, 300 exact matches remain.

| Maximum displayed hits | Old score | Updated score | Updated total hits |
| --- | --- | --- | --- |
| 1 | 50.00 | 0.33 | 300 |
| 250 | 0.40 | 0.33 | 300 |
| 1000 | 0.33 | 0.33 | 300 |

The corrected score uses every retained match within the selected mismatch radius. The visible table may be shorter, but its length no longer controls the calculation. Exported metadata explicitly states whether the table is truncated.

### Example 5 A reverse strand intended site

Take the reverse complement of the 23-base synthetic target and place three A bases before it. The spacer's leftmost reference coordinate is 7, on the minus strand. The program verifies and excludes that coordinate correctly. This test also checks that the PAM position is not confused with the spacer start.

### Example 6 One substitution

Change position 1 of the synthetic spacer from A to T and retain the AGG PAM. With a mismatch radius of one, the program reports one mismatch at position 1. Mismatch positions are numbered along the displayed 20-base spacer. The positions nearest the PAM are at its high-numbered end; the legacy seed summary uses positions 13 through 20 as a reporting convention.

## 6 Worked examples of input and state handling

### Example 7 Repeated FASTA identifiers

Two records whose headers both begin with same previously collapsed into one dictionary entry. The new parser rejects them. Even headers such as same first and same second have the same FASTA identifier because the identifier ends at the first space. Give each record a unique first token.

### Example 8 An ambiguous reference

A panel consisting of 40 N bases produces zero searchable NGG spacer windows. The program returns no specificity score and requires review. It does not turn an unsearchable reference into a score of 100. A partially ambiguous panel also carries a warning and cannot pass the panel review checks.

### Example 9 Missing transcript exon boundaries

An RNA accession without exon annotation is rejected with instructions to use annotated exons or a reviewed genomic region. The program does not infer where the introns were removed. This stricter behavior can reject an accession that the earlier version accepted; the rejection prevents unverified genomic adjacency.

### Example 10 Two short exons

Split a 23-base sequence into exons of 10 and 13 bases. Neither exon independently contains 20 spacer bases plus a PAM. The program does not join them to create a candidate. Exon-boundary-spanning genomic targets may still be valid when supplied as correctly mapped contiguous genomic DNA, but that genomic context must be provided explicitly.

### Example 11 A compound CDS

Use two 30-base CDS parts separated by a 23-base intron containing a target. The updated extraction returns the two 30-base parts independently. The intervening intron is not silently included as part of the coding sequence. This fixes bounding-span extraction; it does not create a general chromosome-coordinate mapper.

### Example 12 A non G starting spacer

For the synthetic A-starting spacer, the format table contains gX19, gX20 and an unmodified X20 option. The table explains which form changes the first base and which adds a base. It does not imply that an extra G has no biological effect, and it does not build the scaffold.

### Examples 13 and 14 State and coordinate checks

Changing the panel text clears the old result. A new successful design also clears old screens, and a failed new design cannot continue displaying the previous target as its result. A forward synthetic target is reported at positions 1 through 20 on the plus strand. These checks complement the reverse-strand example.

Run python benchmarks/run_audit_examples.py to reproduce all 14 examples. The machine-readable outcomes are saved in benchmarks/results/audit_examples.json. Accession-version mismatch and additional invalid-input cases are also tested in the automated suite.

## 7 Installation and first use

Use Python 3.11 or newer. The update was tested on Python 3.12. Download or check out the updated development branch, then open a terminal in the repository directory.

For Windows PowerShell:

```text
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

For macOS or Linux, use source .venv/bin/activate in place of the PowerShell activation command. If activation is restricted on Windows, run .venv\Scripts\python.exe -m pip install -r requirements.txt and .venv\Scripts\python.exe -m streamlit run app.py directly.

Open http://localhost:8501 . Click Load offline example and then Design OpenCRISPR guides. This route does not need NCBI or Ensembl. It should return candidates from both strands. The same input is available in examples/target_demo.fasta.

To try the panel screen, select the ACGT-repeat spacer, open Optional local-reference specificity screen and paste examples/panel_demo.fasta. Keep the mismatch radius at 3. Enable Exclude a verified intended locus and enter intended, start 1, plus strand. Run the screen. The duplicate exact site and the one-substitution site remain visible.

## 8 Use the app for your own analysis

### Choose the right sequence input

Gene lookup accepts a gene and organism. Ensembl mode retrieves individual exons from one selected transcript. Gene lookup uses a transcript-selection heuristic, so it is not an analysis of every isoform. Use Accession ID when you need to specify a particular transcript.

Accession ID accepts an NCBI nucleotide record or an Ensembl stable gene or transcript ID. Choose the matching source. If a requested version differs from the returned version, the app stops. For an older version, obtain the correct archived sequence and annotation rather than accepting the latest record silently.

Manual sequence or FASTA accepts reviewed sequence segments. Paste each independent region as a separate FASTA record. Verify that the sequence belongs to your intended assembly, genotype and region. The app cannot tell from plain sequence alone whether you have supplied genomic DNA, spliced cDNA or the wrong organism.

### Review the candidates

Set the maximum guide count and minimum sequence-quality score, then submit. Read the provenance and warnings before choosing a guide. A large number of candidates can be reduced by the displayed guide-count limit and quality threshold; the displayed table is a ranked subset.

Use the spacer, PAM, strand, segment and local coordinates together. A reverse-strand spacer is shown in its targeting orientation. Its start and end still refer to the leftmost and rightmost spacer bases in the supplied sequence. Do not enter these values as chromosome coordinates without an external mapping step.

### Review guide formats

GX19 is a 20-base spacer whose first G naturally matches the target. gX19 substitutes G for the first base of a non-G spacer. gX20 adds G before the complete 20-base spacer. X20 preserves a non-G-starting spacer without modifying it. Promoter choice and delivery method determine which formats are appropriate.

The local search evaluates the original 20-base targeting spacer. It does not model the biological consequences of selecting gX19 or gX20, and it does not score a full scaffold-containing RNA. Keep the expression choice separate from the search result.

### Run the local screen

Supply the reference FASTA and choose a radius from zero to four substitutions. The default is three. A radius below three is flagged for review because the search excludes additional near-match classes. The value three is a software review convention, not an experimentally established universal threshold.

Leave intended-site exclusion off until you know the exact reference identifier, spacer start and strand. With exclusion off, exact matches are retained and the intended locus remains unconfirmed. With exclusion on, wrong coordinates or a nonmatching spacer cause an error.

The local panel limit is 2 million bases. The target-design limit is 500000 input bases. These are application limits for manageable local work, not guarantees about runtime on every machine. Larger references need an indexed genome-search workflow.

## 9 Interpret the scores and exports

Sequence quality is the existing hand-written heuristic. It uses GC content, poly-T and long homopolymers, with small position-specific bonuses. It helps prioritize sequences according to those rules; a value of 80 does not mean an 80 percent editing rate. GC of 40 to 60 percent is a preferred review band, not a universal pass condition for experiments.

The MIT/Hsu calculation originates from SpCas9 specificity work [5]. This app retains it as a legacy summary of resolved NGG near matches in the supplied panel. It has not been calibrated for OpenCRISPR-1. The aggregate score is 100 divided by one plus the sum of pair-risk values over all matching sites within the selected radius, after excluding a verified intended site if provided. Fewer searched regions or a narrower radius can still make a panel appear cleaner.

| Result | How to read it |
| --- | --- |
| PASS | Implemented checks passed within their stated scope |
| REVIEW | A flag or unresolved context needs attention |
| FAIL | A hard sequence rule failed or a verified panel screen retained an additional exact site |
| NOT SCREENED | No completed panel screen is attached to this candidate |
| No score | The panel contained no searchable resolved NGG spacer windows |

An unscreened guide can pass its sequence checks while specificity remains NOT SCREENED. The app displays that distinction. None of these labels is experimental validation.

CSV provides the candidate table, validation summary and key provenance. FASTA provides targeting spacers, with PAM, segment, coordinates and strand in each header. JSON contains settings, complete validation checks, expression alternatives, reference fingerprint, mismatch counts, intended-locus status and retained hit details.

Keep the original target and reference files with the exports. The normalized sequence fingerprint identifies the sequence representation used by the app. It does not preserve the original raw file byte for byte, prove the source is correct or recover missing sequence data. A capped JSON hit list remains a display subset; total counts are separate fields.

## 10 What was tested

The original checkout passed 37 existing tests, but the reproduced edge cases still exposed defects. Passing tests therefore did not establish that the old implementation was free of those errors.

The updated version passed 79 tests on Python 3.12. Core branch-aware coverage was 90.98 percent, above the repository's 80 percent gate. The new panel-screening and workflow modules each reached 100 percent coverage in this suite. Coverage measures execution of code paths, not biological accuracy.

The suite includes an independent search oracle that reverse-complements the full reference and scans it forward, providing a different route for checking strand coordinates. Three Streamlit AppTest tests exercise real reruns for design, screening, reference changes, setting changes, failed submissions and the offline demo. AppTest is not a browser screenshot review or a production-load test.

All 14 deterministic audit examples passed. A Python wheel for version 1.4.0 built successfully. The tests use mocked responses to check database-retrieval logic. Live Ensembl attempts for HBB and VEGFA timed out; TP53 returned a server error. Successful live retrieval for those genes is not claimed. No guide efficacy, experimental off-target activity, whole-genome sensitivity or wet-lab outcome was measured.

To reproduce the automated checks:

```text
python -m pip install -r requirements-dev.txt
python -m pytest -q --cov=. --cov-config=.coveragerc --cov-fail-under=80
python benchmarks/run_audit_examples.py
```

## 11 Remaining limitations and practical next steps

The most important remaining limitation is biological prediction. The app lacks an OpenCRISPR-specific efficacy model trained and evaluated on appropriate experimental data. Adding a generic SpCas9 score and renaming it would not solve this problem. A future model needs an identified dataset, held-out evaluation, comparison with simple baselines and a clearly stated experimental scope.

The local search does not cover a whole genome, non-NGG PAMs, RNA or DNA bulges, large insertions or deletions, sample-specific variants, chromatin accessibility or editing repair outcomes. It also cannot certify absence of off-target effects. For a research candidate, perform an independent genome-aware analysis using the correct assembly and record the search settings and reference version.

Transcript support still selects one isoform. Ensembl exon input can include untranslated regions. Candidates at exon edges can be missed when the full genomic spacer-plus-PAM context extends outside the isolated exon. Use reviewed contiguous genomic sequence when those boundary targets matter. The app does not automatically map its local coordinates to chromosome positions.

No complete sgRNA scaffold, cloning oligonucleotide, base-editing window, pegRNA, nicking guide or delivery construct is designed. The paper's protein and RNA generative models are not incorporated. The software's feasible role is candidate identification and transparent review, with the downstream design decisions handled separately.

Online retrieval depends on external services and network availability. The offline route is useful for frozen, reviewed inputs, but it cannot verify annotations that you do not supply. The application is intended for local or modest shared use; multi-user workload management and production browser compatibility were not benchmarked in this audit.

For the next development phase, prioritize genomic-coordinate mapping and controlled transcript selection, then an independently benchmarked whole-genome search integration. A learned activity model should follow only when suitable data and validation are available. Complete-scaffold support should use versioned, verified source sequences and explicit format selection.

## 12 Explain the update to your supervisor

The previous version could produce misleading local-screening results because it guessed which exact match to remove and calculated scores from a capped hit list. I updated the software so that intended-site exclusion needs verified coordinates and all matching sites contribute to the score. I also corrected FASTA handling, removed unsafe transcript fallbacks and prevented stale results from being reused.

The updated software can find NGG-compatible targeting spacers, review their sequence properties, check a small supplied reference panel and export reproducible records. The 79 passing tests and 14 worked examples support software correctness within that scope. They do not establish biological efficacy. The next scientific step is genome-aware validation and evaluation against experimentally measured OpenCRISPR-1 data.

## 13 References

[1] Ruffolo, J. A., et al. (2025). Design of highly functional genome editors by modelling CRISPR–Cas sequences. Nature, 645, 518–525. https://doi.org/10.1038/s41586-025-09298-z

[2] Hwang, H. Y., Yi, H., Gwon, Y., Jeon, E., and Kim, D. (2026). High-fidelity genome and prime editing enabled by the AI-designed openCRISPR-1. Genome Medicine, 18, 109. https://doi.org/10.1186/s13073-026-01682-2

[3] Tian, R., et al. (2025). Systematic high-throughput evaluation reveals FrCas9's superior specificity and efficiency for therapeutic genome editing. Science Advances, 11, eadu7334. https://doi.org/10.1126/sciadv.adu7334

[4] Das, P., et al. (2026). AI-designed OpenCRISPR-1 performs robust knockout, base editing, and prime editing in rice. New Phytologist, 251, 2300–2306. https://doi.org/10.1111/nph.71272

[5] Hsu, P. D., et al. (2013). DNA targeting specificity of RNA-guided Cas9 nucleases. Nature Biotechnology, 31, 827–832. https://doi.org/10.1038/nbt.2647

[6] Profluent AI. OpenCRISPR official repository and system resources. https://github.com/Profluent-AI/OpenCRISPR

[7] Behlim, A. B. OpenCRISPR1 gRNA Designer source repository. https://github.com/abdulbasitbehlim/OpenCRISPR1-gRNA-Designer
