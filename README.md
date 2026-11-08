# k-mer evidence and sequence contiguity in sorghum

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24056550.svg)](https://doi.org/10.5281/zenodo.24056550)

This repository reproduces comparisons between pooled k-mer support, support within a read pair, support within one read, and an exact contiguous target match. It includes the derived counters for 663 selected target states in eight same-study holdout accessions, a selected-case sequence-signature census, and a fixed synthetic benchmark.

The catalogue comparisons measure agreement between evidence criteria. They do not identify true genotypes or estimate a genomic false-positive rate. The exact criterion is implemented using literal string matching; this project does not introduce a new pattern-matching algorithm.

## Reproduce the catalogue tables

Python 3.9 or later is required. The table analysis and synthetic benchmark use only the Python standard library. No installation, account, network transfer or cluster path is needed.

From this directory:

```bash
python run_analysis.py --output outputs/catalogue_audit
```

The command checks the supplied input files, reconstructs the pooled minima from the complete word counters, joins the read/pair/exact counters, reconciles the selected signature census, and writes 15 compressed CSV tables. It compares the decompressed content of every table with the study result. Compressed bytes are not compared because zlib and zlib-ng can compress identical text differently. Successful completion prints:

```text
ANALYSIS_STATUS=ANALYSIS_COMPLETE
```

Results are in `outputs/catalogue_audit/tables/`; the denominators and reconciliation counts are in `summary.json`. `checks/reproduction.json` records the content checksums, compressed-file checksums and row counts. The default output directory is relative to the repository, even when the command is invoked from another directory. To repeat the analysis, choose a new output directory. Existing nonempty outputs are never overwritten.

To check the inputs without running calculations:

```bash
python run_analysis.py --verify-only
```

The supplied derived inputs occupy about 24 MB. The catalogue replay took about 70 seconds in a local Linux run with Python 3.12.14. Runtime depends on the host; allow 2 GB of available memory for the table workflow. This estimate does not apply to scanning the original FASTQs.

## What is reproduced

| Route | Starting material | Recomputed result | Scope |
|---|---|---|---|
| `run_analysis.py` | Supplied derived word/group counters and selected-case tables | 15 catalogue, signature, comparison, fraction-rule and repeat-profile tables | No FASTQ rescan or new biological validation |
| `run_benchmark.py` | Target sequences, fixed protocol and seeds | Complete synthetic reads, evidence counts and performance summaries | Constructed sequence models, not genomic prevalence |
| `scan_reads.py` | A user-supplied pair of existing gzip FASTQs and a catalogue | Native word/group counts, prefix checks and bounded examples | One paired library; no automatic download or full study orchestration |
| `data/derived/representative_cases/` | Previously extracted and checked sequence examples | Supplied illustrative sequence records and event annotations | Purposive examples; the case extraction/mapping workflow is outside the table replay |

All required input hashes and their origins are in `data/derived/input_manifest.json`. Study source information and data boundaries are described in `docs/data_provenance.md`. Stable target, query and group identifiers are retained so the derived tables can be traced back to the recorded analyses. Stage labels in some column names and source fields (for example `CODE19` or `was_code15_primary`) identify the analysis step that produced a value. They are retained for traceability and are not commands.

## Definitions and denominators

- **Pooled:** the minimum supporting-pair count among the target's distinct canonical k-mers. Different words may be supplied by different read pairs.
- **Same pair:** every required word appears in the union of the two mates of one observation. Word order and multiplicity are not required.
- **Same read:** every required word appears within at least one mate of one observation. This remains a set criterion, without an order or multiplicity requirement.
- **Exact:** the complete contiguous target or its reverse complement occurs within one mate.
- **Q20:** each base of the particular matching word or target satisfies Phred ≥20. Raw comparisons omit this filter.
- **Sequence signature:** the ordered pair of R1/R2 strings, ignoring base qualities. Identical signatures need not represent independent molecules.

The primary catalogue comparison requires at least two read-pair observations at Q20. The fixed catalogue contains 5,304 accession-state cells. Eligibility requires `k < target length`, so k81 is available for 657 states and 5,256 cells. Six states are ineligible for k81; they are not zero-support observations. Both baseline-specific and common-denominator results are retained. `near_full` means target length minus one and can be an alias of a fixed-k group.

The selected signature census covers 60 accession-state cells spanning 37 states. Signature counts for the remaining 5,244 cells are unknown and remain blank. They are not inferred from a different accession with the same target state. The 663 states share 595 underlying anchor pairs and are not independent random genomic loci.

At the primary setting, pooled/same-pair/same-read k31 exceed exact in 60/8/1 cells. Same-read k61 and k71 counts equal exact across the complete catalogue under both raw and Q20 policies. In this catalogue, 659 of 663 targets have no internal repeat of 20 bases or longer on either strand (`target_repeat_profile`), so the observed equality describes a target set with few long internal repeats; the repeat profile does not guarantee that equality. Same-pair k61 has equal primary calls but not equal counts; the count and threshold comparisons remain separate.

The `fraction_rule` table tests a fraction threshold within these evidence-counting rules: a target passes when at least a stated fraction of its distinct k-mers each reach the support threshold. A fraction of 1.0 is the pooled rule above. At Q20 and support 2, k31 targets pass the fraction rule but remain below the exact-support threshold in 60 of 5,304 cells at fraction 1.0 and in 1,771 of 5,304 cells at fraction 0.4. These groups include 14 and 28 cells, respectively, with one exact-supporting read pair; the remaining cells have zero. The corresponding 0.4 count in the common 5,256-cell set is 1,740. This calculation does not reproduce the indexing or query procedures of REINDEER or BIGSI.

## Regenerate the synthetic benchmark

```bash
python run_benchmark.py --output outputs/synthetic_benchmark
```

This executes the complete fixed grid: 37 targets, six generating models, three noise models, two read lengths, two read-pair counts and 20 replicates, giving 53,280 synthetic samples and 1,065,600 paired observations. Every detection criterion sees the same observed reads within a sample. Noise variants share fragment geometry. The primary synthetic metric uses Q20 and at least two distinct ordered sequence signatures.

The constructed complementary-variant mixtures are mechanism stress tests. They are not estimates of how often such mixtures occur in genomes. The model includes substitutions and a separate quality-mask condition, but no indels, calibrated sequencing-error model or whole-genome background. The benchmark does not label the 60 real-data discrepancies as biological false positives. It creates no new reference or FASTQ downloads.

`config/benchmark_protocol.json` documents the seeds, geometry, truth definition and complete condition grid. `config/target_panel.json` contains the selected sequences. The output preserves per-sample counts, generating-template records, generation audits, per-stratum performance and paired-method comparisons. Eight output files are checked by decompressed content. Two of them summarize the completed grid: `dilution_sampling_ceiling` counts diluted trials that sampled zero, one or at least two target-origin read pairs, and `support_unit_agreement` counts calls that differ between read-pair and sequence-signature units. The complete replay took about 212 seconds in a local Linux run with Python 3.12.14. Runtime depends on the host. Execution records and reproduction checks are listed in `validation.json`.

## Optional native FASTQ scanner

This route is for users who already have the required read files. It can take hours on a large dataset and is not part of the quick-start analysis. It requires Linux or another suitable C++17 environment, a compiler, and zlib development headers. A platform-specific binary is not distributed.

```bash
python build_scanner.py
python scan_reads.py --r1 /path/to/reads_1.fastq.gz --r2 /path/to/reads_2.fastq.gz --output outputs/read_scan
```

The scanner streams both complete gzip inputs, checks mate identifiers and record counts, checks gzip completion/CRC, and compares an initial prefix with a separate direct-search implementation. It writes occurrence, read-support and paired-support counters under `counts/`. It does not calculate the complete sequence-signature census. Memory use depends on the target catalogue; the small table-replay memory estimate must not be used to size a full FASTQ scan.

Inputs must use four lines per FASTQ record, Phred+33 qualities, read lengths of 1–512 bases, lines no longer than 4,096 bytes, and a final newline. Wrapped FASTQ records and long reads are not supported. Ambiguous DNA bases interrupt an exact A/C/G/T match. The two inputs must contain corresponding mates in the same order.

The source is the implementation used for the full-catalogue scan. The portable build and command-line wrapper were checked on small paired-read fixtures. The original approximately 130 GB FASTQ collection was not rescanned when assembling this release. See `docs/reproducibility.md` for the exact validation boundary.

## Tests

```bash
python -m unittest discover -s tests -v
```

Tests include content hashing across compression settings, the fraction rule, repeat profiles, independent literal-search comparisons, reverse complements, repeated/low-complexity targets, quality boundaries, paired evidence, duplicate signatures, missing zero rows, baseline eligibility, selected-case missingness and fixed simulation geometry. Optional native tests compile the scanner and compare its counts with the Python detector on small fixtures; they are skipped if `g++` is unavailable. A compiler without zlib headers cannot run those native tests.

## Repository contents

| Path | Contents |
|---|---|
| `src/kmer_evidence/` | Catalogue calculations, signature reconciliation, matcher and simulation |
| `config/` | Fixed protocols, target panel and expected output hashes |
| `data/derived/` | Minimal derived counters, selected-case evidence and input provenance |
| `native/` | C++17 streaming scanner source |
| `tests/` | Calculation, simulation, interface and scanner checks |
| `docs/` | Definitions, source information and reproduction boundaries |

The repository contains analysis code and machine-readable inputs/results. Manuscripts, finished figures and figure-production scripts are maintained separately in the study's author archive.

## Release and citation

Repository: https://github.com/EJSAHN/sorghum-kmer-evidence. The code is released under the MIT License (`LICENSE`). Release v1.0.0 is archived at Zenodo (https://doi.org/10.5281/zenodo.24056550); citation metadata are in `CITATION.cff`. Third-party source data retain their original terms. Scientific source lineage is recorded in `docs/code_provenance.json`; exact replay checks and measured runtimes are in `validation.json`.
