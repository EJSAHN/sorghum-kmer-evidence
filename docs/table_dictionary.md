# Table guide

All generated CSV files are gzip-compressed and use named columns. Integer support counts represent paired FASTQ observations unless the column explicitly states `signatures`. `pi_accession`, `state_id`, `pair_id`, `group_id` and query identifiers preserve the original study identities.

| Output table | Unit and use |
|---|---|
| `native_evidence` | Accession × native group × quality; reconstructed pooled minima and observed pair/read/exact counts |
| `baseline_coverage` | State × baseline; eligible k values, native-group aliases and exclusion reasons |
| `full_catalog_calls` | Accession × state × eligible baseline × quality × support threshold |
| `baseline_summary` | Counts and evidence classes for full, common and prespecified descriptive subsets; subsets can overlap |
| `count_equivalence` | Numeric agreement with exact counts and disagreements at support thresholds 1, 2 and 3 |
| `residual_cases` | Conditions with same-read counts exceeding exact; repeated rows across settings are not distinct biological cases |
| `simple_rule_comparison` | Same-denominator contingency counts for marginal rules and their cell-level conjunctions |
| `signature_coverage` | All 5,304 cells, identifying the 60 measured signature cells and retaining blanks elsewhere |
| `selected_signature_sensitivity` | Alternative quality, support and counting-unit summaries within the same 60 selected cells |
| `selected_case_baselines` | Detailed raw-pair and signature counts for the selected cells, including the full-span equivalence check |
| `incremental_case_context` | Descriptive context for numeric same-read/exact differences; multiple rows can refer to one cell |
| `simulation_context` | Primary synthetic stratum summaries kept separate from real-data denominators |
| `fraction_rule` | Scope × baseline × quality × support threshold × required k-mer fraction; cells passing the fraction rule, exact-positive cells and cells passing while exact support remains below the same threshold |
| `target_repeat_profile` | One row per target state; longest substring occurring at least twice on either strand |
| `primary_comparison` | Primary raw-pair support ≥2 at Q20, with both all-eligible and common fixed-k denominators |

`CONCORDANT_POSITIVE`, `DISCORDANT` and `BELOW_THRESHOLD` are fixed from primary k31 evidence. They are descriptive evidence classes, not known biological truth labels. The primary fixed strata partition the catalogue; other exported scope views may overlap and should not be added together.

`NOT_MEASURED` means that no signature census is available for that cell. An empty signature field is not a zero. `near_full` uses target length minus one. Repeated labels can refer to one native group and do not constitute independent replications.

The four selected-case evidence labels describe quality-threshold loss, observed sequence complementation, observed extra repeat sequence, or an unresolved mechanism. These labels do not establish genome-wide absence or identify a causal locus. Stage labels in column names and values (for example `CODE20` or `was_code15_primary`) identify the analysis step that produced a value and are retained for traceability.

The synthetic benchmark also writes `dilution_sampling_ceiling` (diluted trials by the number of read pairs drawn from the target-bearing template) and `support_unit_agreement` (calls that differ between read-pair and sequence-signature units, by read length, quality, criterion and threshold).

Columns containing `without_exact` retain their original names and describe calls below the exact-support threshold. For example, `fraction_rule_positive_without_exact_cells` retains its original column name. It counts fraction-rule positives with exact support **below `minimum_support`**, not necessarily zero exact observations. At the primary threshold of two, exact counts of both zero and one fall below the threshold.
