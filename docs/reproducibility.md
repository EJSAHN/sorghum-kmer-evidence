# Reproduction boundaries

The default command reconstructs evidence summaries from the complete supplied native word and group counters. It is not a replay of the original raw FASTQ collection. The original sequence discovery, sample selection, download verification, reference mapping and selected-case extraction are upstream steps.

## Catalogue replay

`run_analysis.py` reconstructs pooled support, reconciles the native counters, enumerates all eligible accession-state-baseline-quality-threshold combinations, and writes the catalogue and selected-case tables. It verifies 222,480 comparison rows against the supplied earlier comparison table, including all-zero rows. The raw/Q20 settings and support thresholds 1, 2 and 3 are fixed in the protocol. The primary same-pair/read counts are not substituted for the sequence-signature counts.

All 15 generated tables are checked by the SHA-256 of their decompressed content. Thirteen match the recorded Ceres output; the fraction-rule and repeat-profile tables are computed from the same supplied counters and target sequences. Compressed bytes are not compared, because zlib and zlib-ng produce different compressed bytes for identical text. This establishes reproduction of the derived-counter calculation; it is not an independent scanner implementation or external biological replication.

The selected-case overlay reconciles 10,080 raw counters against the catalogue before attaching sequence-signature counts. Case context is joined by accession and state. Missing signature counts outside the 60 selected cells remain unmeasured. Composite rules are conjunctions of cell-level marginal calls, not claims that the same read pairs satisfy both criteria.

## Synthetic replay

`run_benchmark.py` starts with the fixed target sequences and generates all templates and reads. The random seeds are derived with SHA256 and Python `random.Random`; Python's process-randomized hash is not used. The runner checks eight output files by decompressed content: the six original outputs, including per-cell counts and generation audits, and two summary tables computed from them (`dilution_sampling_ceiling` and `support_unit_agreement`).

The same-read detector, direct sequence tests and full-length equivalence checks are covered by unit tests. The generator fixes target truth before read sampling and sequencing errors. Samples with no observed target-origin read are retained in sensitivity denominators. Repeated trials share targets and noise geometry; there is no claim that the combined counts are independent genomic trials.

## Native scanner

The distributed C++ source was used for the catalogue counters. The release validation builds it with a local compiler and checks small fixtures against the Python detector, including malformed gzip streams, mismatched mate identifiers and unequal mate counts that must fail. The new portable wrapper does not reproduce all historical download, scheduling or cache-management steps. It scans one paired library at a time and does not automatically aggregate many accession/run files or deduplicate signatures across them.

## Publication inputs

Representative-case sequences and annotations are supplied as derived illustrative records. Their source identifiers, original-file hashes and transformation descriptions are preserved. Supplying these records does not imply that all unresolved biological mechanisms have been determined. Finished figures and the manuscript are outside this analysis repository.

## Execution records

Successful table replay writes `summary.json`, `status.json` and `checks/reproduction.json`. Synthetic replay writes `status.json`, `reproduction.json` and the detailed files under `analysis/`. Failed runs preserve a `FAILED` status and error text. Each rerun requires a separate empty output directory.

Release-validation measurements and substantive file hashes are recorded in `validation.json`. The packaged tests were executed locally on Linux. The original numerical results were computed on Ceres with Python 3.9. The tests have not been run on every operating system.

Seeds, configuration, inputs, simulation code paths and C++ scanner source are unchanged from the code that produced the published numbers. `validation_v1.json` and `validation_v3.json` preserve earlier execution records. `validation.json` records the most recent checks.
