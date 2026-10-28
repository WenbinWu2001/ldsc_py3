# Artifact Metadata Field Inventory

Last updated on: 2026-10-28

Munged data filenames use the filesystem-safe trait label when supplied: `<trait>.parquet` and optional `<trait>.sumstats.gz`. The `sumstats.parquet` and `sumstats.gz` names below describe runs without a trait label. See [munging output artifacts](../munge-sumstats.md#output-artifacts) for naming and overwrite rules.

Downstream identity metadata lives in the `sumstats.parquet` footer (for munged summary statistics) and in `ldscore/metadata.json` (for LD scores).
Any metadata emitted by `annotate`, `build-r2-matrix`, `h2`, `partitioned-h2`, `rg`, query-level, or pair-level outputs is diagnostic provenance for regression.
The post-fit `quantile-h2` command is the explicit exception: it consumes one selected `partitioned-h2` model's metadata and coefficient-delete artifact because those files identify the fitted model.

This document inventories metadata, schema metadata, and audit sidecars written by the current LDSC3 workflows.
The goal is to keep downstream contracts obvious: only reusable input artifacts carry required identity metadata.
Diagnostic metadata lives under `diagnostics/`; a post-fit command may consume it only when the selected result directory itself is its declared input contract.

## Boundary Rule

Downstream-required identity metadata has exactly two carriers:

- The **`sumstats.parquet` footer** (discrete `ldsc:*` keys) for munged summary statistics.
- **`ldscore/metadata.json`** for LD-score directories.

Everything else is **diagnostic provenance:** humans and external tools may
inspect it for reproducibility/debugging, but LDSC3 downstream workflows must not
load it, require it, or change behavior when it is present.

There is no "optional but used if present" metadata. No JSON metadata file
contains a top-level `format`; `artifact_type` is the schema discriminator.

## Output Ownership Rule

Workflow-owned artifacts are the files and directories in the current public
layout for that workflow. Without `--overwrite`, any existing owned artifact
blocks the run before outputs are opened. With `--overwrite`, artifacts written
by the current run are replaced, and stale owned siblings left by incompatible
flag combinations are removed after the current write succeeds. Sharded
workflows may scope this owned package to the shard selected by the current
invocation; for `build-r2-matrix`, concrete chromosome runs own only that
chromosome's files, while full `@` suite runs own the all-chromosome package.

Legacy/root diagnostic names that are no longer in the public layout are not
owned by current workflows. They are ignored by preflight and cleanup rather
than blocked or deleted.

Directory artifacts are owned as whole trees. Optional
`diagnostics/query_annotations/` and `diagnostics/pairs/` directories are staged
under `diagnostics/` and then moved into final position as a unit; reruns that
omit those optional details remove the old tree after a successful overwrite.
Always-written audit sidecars such as workflow logs and dropped-SNP tables are
included in the owned family for ordinary result-directory workflows.
The gene LD-score index builder is the deliberate exception while running: its open log, history, and lock live in hidden sibling `.<index-name>.build-state/`.
Once the handler closes, the successful log moves into the published index diagnostics.

## Public Output Layouts

### `munge-sumstats`

```text
sumstats/
  sumstats.parquet
  sumstats.gz
  diagnostics/
    sumstats.log
    dropped_snps/dropped.tsv.gz
```

`sumstats.gz` is present only when the selected output format writes
the legacy TSV artifact.

The downstream-required identity metadata is embedded in the `sumstats.parquet`
**footer** as discrete `ldsc:*` keys; no `metadata.json` is written. `.sumstats.gz`
carries no embedded metadata. Legacy `.sumstats` and `.sumstats.gz` inputs are
marked for rsID-to-panel projection. Footerless Parquet is rejected because it
is neither a self-describing LDSC3 artifact nor a legacy compatibility format.

| Footer key | Explanation | Downstream usage |
| --- | --- | --- |
| `ldsc:artifact_type` | Must be `sumstats`. | Required and validated. |
| `ldsc:snp_identifier` | SNP identifier mode used when the munged artifact was written. | Recovered for regression compatibility checks. |
| `ldsc:genome_build` | Final output genome build for coordinate-family modes; empty (`null`) for rsID-family modes because rsID identity is build-independent. | Recovered for the regression genome-build check when known. |
| `ldsc:trait_name` | Trait label from the raw summary-statistics config, when supplied; empty otherwise. | Used for output labels, not compatibility. |
| `ldsc:qc_info_min`, `ldsc:qc_maf_min` | Effective INFO and MAF thresholds of the munge run; empty when that filter was not applied. `SumstatsMunger.write_output` keeps them for a table returned by `run` and omits them otherwise. | Provenance only; not read. |

### Summary-statistics `diagnostics/dropped_snps/dropped.tsv.gz`

Always-owned audit sidecar for liftover and identity-cleanup drops. This file is
written even for clean runs so users can distinguish "no drops" from "missing
sidecar".

Detailed source-build inference, requested output build, and liftover method
provenance live in `diagnostics/sumstats.log`, not the parquet footer. The
footer metadata and `SumstatsTable.config_snapshot` intentionally expose
only the final downstream compatibility build.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `CHR` | Chromosome for the dropped row when available. | None. |
| `SNP` | SNP label for the dropped row when available. | None. |
| `source_pos` | Source coordinate before liftover or cleanup. | None. |
| `target_pos` | Target coordinate after liftover when available. | None. |
| `reason` | Drop reason. | None. |
| `base_key` | Base-mode identity key (without alleles) involved in cleanup. | None. |
| `identity_key` | Effective allele-aware identity key when applicable. | None. |
| `allele_set` | Normalized allele set used by allele-aware cleanup. | None. |
| `stage` | Workflow stage that produced the drop. | None. |

### Regression legacy summary-statistics audit

When any regression input is legacy text and an output directory is
supplied, `h2`, `partitioned-h2`, and `rg` write
`diagnostics/dropped_snps/legacy_sumstats.tsv.gz`, including a header-only file
for a clean projection. Its fixed fields are `trait_name`, `source_path`,
`SNP`, `A1`, `A2`, `reason`, and `panel_candidate_count`. It is diagnostic only;
the in-memory projected table, not this sidecar, enters regression.

### `build-r2-matrix`

```text
ref-panel/
  <build>/chr<chrom>_r2.parquet
  <build>/chr<chrom>_meta.tsv.gz
  diagnostics/
    metadata.json
    metadata.chr<chrom>.json
    build-r2-matrix.log
    build-r2-matrix.chr<chrom>.log
    dropped_snps/chr<chrom>_dropped.tsv.gz
```

The concrete diagnostics files are `metadata.json` and `build-r2-matrix.log`
for multi-chromosome runs. Concrete single-chromosome runs write
`metadata.chr<chrom>.json` and `build-r2-matrix.chr<chrom>.log` so parallel
per-chromosome invocations can share an output directory without racing on
shared diagnostics.

`diagnostics/metadata*.json` is provenance only.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `artifact_type` | Must be `r2_matrix_build`. | None. |
| `files` | Relative map of emitted R² Parquet files, sidecars, and diagnostic files. | None. |
| `snp_identifier` | SNP identifier mode used for emitted R² matrix artifacts. | None for root metadata. Per-file schema/header metadata carries runtime validation. |
| `source_genome_build` | Build of source PLINK coordinates. | None for root metadata. |
| `emitted_genome_builds` | Builds emitted by the builder. | None. |
| `chromosomes` | Chromosomes emitted by the run. | None. |

### R² matrix Parquet and sidecar metadata

Runtime compatibility is carried by the Parquet footer and paired SNP sidecar, not the diagnostic JSON. Required columns are `IDX_1`, `IDX_2`, `R2`, and `SIGN_R`; indices reference the exact sidecar row order.

| Field | Location | Downstream usage |
| --- | --- | --- |
| `ldsc:artifact_type` | Parquet footer or sidecar header | Identifies `r2_matrix` or `ref_panel_metadata` |
| `ldsc:snp_identifier` | Both | SNP identifier mode of the build; `ldscore` requires an exact match, and `query-r2-matrix` accepts this mode or its base form |
| `ldsc:genome_build` | Both | Coordinate-build identity |
| `ldsc:sorted_by_build` | Parquet footer | Declares the build defining index/position order |
| `ldsc:n_snps`, `ldsc:sidecar_identity_sha256` | Parquet footer | Bind pair indices to the sidecar's SNP count, identity, and order |
| `ldsc:n_samples` | Parquet footer | Reference sample size (at least 3, equal across an R² matrix's chromosome files) that `query-r2-matrix` uses to recover Pearson r |
| `ldsc:r2_bias` | Parquet footer | R² matrices store `unbiased` R² |
| `ldsc:r2_encoding`, `ldsc:r2_scale` | Parquet footer | Package output uses `int16_symmetric` with divisor 32767 |
| `ldsc:r2_min` | Parquet footer | Required build-time sparsification threshold; nonpositive means no R² filtering; `ldscore` rejects a positive value |
| `ldsc:ld_window_mode`, `ldsc:ld_window_value` | Parquet footer | Build window used for downstream compatibility checks |
| `ldsc:row_group_size` | Parquet footer | Informational intended row-group size |

See the [complete format contract](r2-matrices.md) and `write_r2_parquet` in `src/ldsc/_refpanel/build_kernel.py`.

### `query-r2-matrix`

```text
pair_query/
  query_r2_matrix.tsv
  diagnostics/
    metadata.json
    query-r2-matrix.log
```

The TSV preserves input columns and writes `r2`, nullable `sign_r` (+1/-1), signed Pearson `r`, and `status`. `sign_r` and `r` use query allele orientation; they are not panel-orientation labels. There is no `sign` alias. Missing values are serialized as `NaN`.

Diagnostic metadata records `artifact_type=query_r2_matrix_result`, the relative `files` map, `r2_matrix_dir`, `snp_identifier`, the requested `genome_build`, `n_samples`, `n_pairs`, and `status_counts`. The empty successful status is counted as `resolved`; other counts use `not_in_matrix`, `cross_chromosome`, or `absent`. This JSON is provenance only. Sources: `_query_r2_matrix_metadata` in `src/ldsc/_refpanel/query.py` and `QueryR2MatrixDirectoryWriter` in `src/ldsc/_refpanel/outputs.py`.

### `annotate`

```text
annotate/
  query.<chrom>.annot.gz
  baseline.<chrom>.annot.gz                  # annotation grid only
  diagnostics/
    metadata.json
    annotate.log
    dropped_snps/dropped.tsv.gz
    query_annotation_status.tsv
    gene_list_audit.tsv.gz                   # gene lists only
    gene_list_resolution_summary.tsv         # gene lists only
    gene_catalog_issues.tsv                  # when catalog issues exist
    input_issues.tsv                         # when input issues exist
    chromosome_scope.json                    # when a scope record exists
```

`diagnostics/metadata.json` is provenance only; no downstream workflow reads it.
Fields come from `AnnotationDirectoryWriter.write()` in `src/ldsc/_annotation/outputs.py` and the provenance built by `run_annotate()` in `src/ldsc/_annotation/annotate.py`.

| Field | Explanation |
| --- | --- |
| `artifact_type` | Always `annotation_projection`. |
| `files` | Relative paths of the annotation files and diagnostics. |
| `snp_identifier`, `genome_build` | SNP identifier mode and genome build of the run. |
| `grid_kind` | `annotation` or `bim`. |
| `ref_snp_annot_sources`, `ref_snp_bim_sources` | Grid inputs as supplied; the unused one is empty. |
| `query_source_kind` | `bed` or `gene_list`. |
| `query_sources`, `query_annot_bed_sources`, `query_annot_gene_list_sources` | Query inputs as supplied. |
| `gene_coordinate_file`, `gene_catalog_build` | Gene coordinate catalog and its declared build; empty for BED queries. |
| `projection_genome_build` | Concrete build used to project intervals. |
| `padding_bp` | Padding applied to every interval. |
| `gene_exclude_regions`, `gene_list_resolution_policy` | Gene exclusion policy and resolution policy (`strict` or `resolved-only`). |
| `chromosomes` | Chromosomes written. |
| `n_snps` | SNP rows kept after SNP identity cleanup. |
| `snp_rows` | Input, output and dropped logical SNP rows, in total and under `by_chromosome`. |
| `baseline_columns`, `query_columns` | Exported baseline columns (annotation grid only) and query columns. |

### `ldscore`

```text
ldscore/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet
  ldscore.overlap.parquet
  ref_metadata/chr<chrom>_meta.tsv.gz      # with --export-ref-metadata
  diagnostics/
    ldscore.log
    chromosome-logs/chr<chrom>.log
    dropped_snps/chr<chrom>_dropped.tsv.gz
    query_annotation_status.tsv
    gene_list_audit.tsv.gz
    gene_list_resolution_summary.tsv
    gene_catalog_issues.tsv                 # when catalog issues exist
    chromosome_scope.json
    input_issues.tsv
    ld_window_snp_counts.tsv
    reference_snp_funnel.tsv
```

`ldscore.query.parquet` is present for a single query batch. Multiple batches instead use `ldscore.query.batch00001.parquet` and subsequent ordinals, with `files.query_batch00001` keys. Baseline-only output has no query file. All current directories require the ordered `query_batches` manifest, including an empty list for baseline-only output.
`query_annotation_status.tsv` is present for BED/gene-list query runs;
The row-complete audit and per-source summary are present for gene-list runs.
`ld_window_snp_counts.tsv` is present only for R² matrix runs with `--ld-wind-snps`; it records each chromosome's `panel_snps`, `kept_snps`, `ld_wind_snps`, and `counting_basis` (`panel_rows`).
`reference_snp_funnel.tsv` is present for every native and indexed run; it records each chromosome's reference-SNP funnel with `funnel`, `chrom`, `backend`, `step_index`, `step`, `removed` (`NA` for a step that did not run) and `remaining`, as described in [reference SNPs on the annotation SNP grid](../ldscore-snp-universes.md#reference-snp-funnel).
It is owned by the LD-score writer but not listed in the `files` entry of `metadata.json`.
`ldscore.overlap.parquet` holds the annotation overlap matrix (long form:
`row_annotation`, `col_annotation`, `overlap_all_snps`, `overlap_common_snps`)
that `partitioned-h2` requires. It is written only for runs with two or more
annotation columns; a single-annotation (e.g. base-only) run omits it.
Its storage grows as $O(B(B+Q))$ for $B$ baseline and $Q$ query columns.
Only the baseline-row block $A_B^\top A$ and each query's self-overlap are stored, which is enough to rebuild the overlap matrix of every baseline-plus-one-query model; `overlap_config.stored_block` is `baseline_rows_plus_query_diagonal` with queries and `baseline_by_baseline` without.
`partitioned-h2` requires the file, the `h2` collinearity check reads it when present, and `rg` ignores it.
Query files repeat the SNP identity columns (`CHR`, `SNP`, `POS`, and `A1` and `A2` when available).
Opening an LD-score directory validates every declared schema and allele column, and `LDScoreSource.read_queries` checks that the selected query files share the same SNP rows (`src/ldsc/_ldscore/source.py`).
The name-keyed long format survives zero-variance column drops and reordering, because regression assembles each model's overlap matrix by annotation name.
Embedding the matrix in `metadata.json` was rejected because it would reach tens of MB for cell-type panels; `.npz` or position-aligned records were rejected as fragile.

`metadata.json` is downstream-required. Query workflows also persist `diagnostics/chromosome_scope.json`, including validated baseline/reference chromosome sets, the content-resolved PLINK chromosome-to-prefix mapping when applicable, the effective analysis scope, validation status, and ordinary-glob selection caveat. Structural input failures populate `diagnostics/input_issues.tsv` with `input_role`, `source`, `chrom`, `reason`, `details`, and `repair`; unavailable later assessments remain unknown. These diagnostics can exist without scientific outputs after a failed preflight.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `artifact_type` | Must be `ldscore`. | Required and validated. |
| `files` | Relative data-file map with required `baseline` and optional `query` and `overlap`. | Required to locate parquet data; `overlap` is required by `partitioned-h2`. |
| `snp_identifier` | SNP identifier mode used for LD-score rows. | Required for regression compatibility checks. |
| `genome_build` | Genome build associated with LD-score identity metadata. | Required for regression compatibility checks when known. |
| `chromosomes` | Chromosome labels represented in the output. | Reporting/navigation. |
| `baseline_columns` | Ordered baseline annotation LD-score columns. | Required to assemble regression covariates. |
| `query_columns` | Ordered query annotation LD-score columns. | Required for `partitioned-h2` query selection. |
| `counts` | Per-annotation count records, keyed by name: `group`, `column`, `all_reference_snp_count` and, when common-SNP counts are available, `common_reference_snp_count`. Regression uses the all-SNP count when the common-SNP count is absent. | Required for regression count vectors. |
| `count_config` | Common-SNP count settings `common_reference_snp_maf_min` (default 0.05) and `common_reference_snp_maf_operator`, the actual threshold operator (`>=` for native LDSC3 computation; strict `>` for legacy LDSC conversion). | Required to interpret the common-SNP counts and checked against overlap metadata. |
| `overlap_config` | Overlap-matrix provenance: `total_all_reference_snps`, `total_common_reference_snps`, `common_maf_min`, `common_maf_operator`, `stored_block`. `null` for single-annotation runs that write no overlap matrix. | Provides `M_tot` and the all-SNP and common-SNP definitions for overlap-aware `partitioned-h2`. |
| `snp_universe_policy` | `ld_reference_universe` (`ref_panel_snps_file`, the explicit reference-panel restriction or `null`, and `retained_snp_count`, the reference SNPs on the annotation SNP grid; `diagnostics/reference_snp_funnel.tsv` records how they were selected) and `regression_rows_and_weights` (regression-SNP selection, region exclusions, and counts before and after exclusion). | Scientific provenance; ignored by regression. |
| `annotation_types` | Per-column `binary` or `quantitative` classification. | Interpretation and logging only; never changes fitting. |
| `n_baseline_rows` | Number of rows in the baseline parquet table. | Reporting. |
| `n_query_rows` | Genome-wide SNP row count shared by every query batch file, or zero without queries. | Reporting. |
| `row_group_layout` | Row-group strategy. | Reporting/technical provenance. |
| `baseline_row_groups` | Row-group metadata for `ldscore.baseline.parquet`. | Reporting/technical provenance. |
| `query_row_groups` | Row-group metadata for a single `ldscore.query.parquet`, or `null` for multiple/no query batches. | Reporting/technical provenance. |
| `query_batches` | Ordered entries containing `file`, ordered `query_columns`, and `row_groups` with chromosome/index/offset/row count. | Required canonical query discovery and selected reads; older directories must be regenerated. |
| `gene_list_resolution_policy` | `strict` when `--allow-unresolved-genes` is omitted; `resolved-only` when supplied. The metadata field and Python policy strings are unchanged by the CLI rename. | Scientific provenance; ignored by regression. |
| `gene_list_resolution_counts` | Aggregate submitted, rejected, and unique-resolved counts. | Reproducibility/diagnostics; ignored by regression. |
| `query_diagnostics` | Relative paths to available query status, scope, structural input issues, and, for gene runs, the audit and source summary. | Troubleshooting/navigation; ignored by regression. |
| `chromosome_scope` | Validated input scope and selection semantics, separately from chromosomes with retained output rows. Empty for workflows without query coverage preflight. | Scientific provenance; ignored by regression. |
| `legacy_ldscore_conversion` | Present only for explicit legacy LDSC conversion: profile, selected source directories/prefixes/files, rsID intersection counts, count origins, strict common-frequency rule, coordinate evidence, and diagnostic paths. Source-file hashes are not computed or written. | Provenance and compatibility auditing; regression still consumes the ordinary canonical fields. |

Native LDSC3 computation records an inclusive common operator (`>=`).
Explicit legacy LDSC conversion records the actual strict legacy operator (`>`) at threshold `0.05`; `overlap_config` must agree with `count_config`, and the loader rejects drift.
Converted suites use the same canonical directory layout.
Their workflow log is `diagnostics/convert-legacy-ldscores.log` and their always-written issue audit is `diagnostics/conversion_issues.tsv.gz`; they do not also claim a native `diagnostics/ldscore.log`.

`diagnostics/query_annotation_status.tsv` has fixed columns `query`, `source`,
`input_type`, `status`, `reason`, `n_annotation_snps`, and `details`.
`diagnostics/gene_list_audit.tsv.gz` is row-complete and
`diagnostics/gene_list_resolution_summary.tsv` has one row per focal/control
source. Their fixed schemas and null rules are documented in
[Gene-list diagnostics and repair](../results/gene-list-diagnostics.md).

### `h2`

```text
h2/
  h2.tsv
  diagnostics/
    ld_score_regression_bins.tsv
    metadata.json
    h2.log
```

`diagnostics/metadata.json` is provenance plus the plotting navigation contract.
`--output-dir` is required. The bin table contains the exact final-fit
rank-bin summaries used by `ldsc plot`; it is a required current `h2` artifact.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `artifact_type` | Must be `h2_result`. | Selects `h2` plotting and conversion contracts. |
| `files` | Relative map with `summary: "h2.tsv"` and `ld_score_regression_bins: "diagnostics/ld_score_regression_bins.tsv"`. | Conversion follows `summary`; plotting follows `ld_score_regression_bins`. |
| `trait_name` | Trait label resolved from CLI input or summary-statistics metadata. | None. |
| `sumstats_file` | Source summary-statistics path. | None. |
| `ldscore_dir` | Source LD-score directory. | None. |
| `effective_snp_identifier` | SNP identifier mode actually used after any allowed identity downgrade. | None. |
| `genome_build` | Regression dataset genome build when available. | None. |
| `identity_downgrade_applied` | Whether same-family identity downgrade was used. | None. |
| `count_key_used_for_regression` | Count vector used for regression. | None. |
| `retained_ld_columns` | LD-score columns retained in the fitted `h2` model. | None. |
| `dropped_zero_variance_ld_columns` | LD-score columns removed before fitting. | None. |
| `n_snps` | Number of SNPs retained after the effective chi-square filter. | None. |
| `effective_chisq_max` | Chi-square cap actually applied (`null` when uncapped). | None. |
| `n_blocks_used` | Realized jackknife block count from the fit outcome. | None. |
| `samp_prev` | Sample (case) prevalence applied for liability conversion, or `null` for an observed-scale run. | None. |
| `pop_prev` | Population prevalence applied, or `null`. | None. |
| `scale` | `observed` or `liability` (the latter when both prevalences are supplied). | None. |

### `partitioned-h2`

```text
partitioned-h2/
  partitioned_h2.tsv
  diagnostics/
    metadata.json
    coefficient_delete_values.parquet
    partitioned-h2.log
    query_status.tsv
    query_annotations/manifest.tsv
    query_annotations/<query>/metadata.json
    query_annotations/<query>/partitioned_h2.tsv
    query_annotations/<query>/partitioned_h2_full.tsv
    query_annotations/<query>/coefficient_delete_values.parquet
```

The `diagnostics/query_annotations/` tree is present for every query-annotation
run and absent for baseline-only runs. `--output-dir` is required; the retired
`--write-per-query-results` flag is rejected.

The root `partitioned_h2.tsv` includes `query_result_dir`, a path relative to the result root that matches the folder in `diagnostics/query_annotations/manifest.tsv`. Baseline-only summaries use `NaN`. The navigation column is absent from per-query tables and in-memory Python summaries.

Root `diagnostics/metadata.json` stores run-wide provenance once.
It describes the analysis through `analysis_type` (`functional_category` | `cell_type_specific`), `headline_metric` (`enrichment` | `coefficient`), `enrichment_p_test` (`two_sided_t`), and `coefficient_p_test` (`one_sided_greater`).
Each per-query `metadata.json` stores fit-specific fields and `global_metadata: "../../metadata.json"`, resolved relative to that query folder.
Its `files` paths remain relative to the aggregate result root.
The post-fit `quantile-h2` reader requires the enclosing shared metadata when this reference is present and merges it with query-specific fields taking precedence.
Keep the complete result directory when copying fitted queries.
Sources: `PartitionedH2DirectoryWriter` in [_regression/outputs.py](../../../src/ldsc/_regression/outputs.py) and `load_fitted_partitioned_model` in [_regression/quantile_h2.py](../../../src/ldsc/_regression/quantile_h2.py).

Model metadata records `n_snps`, `effective_chisq_max`, and `n_blocks_used` from the fitted population: at the root for a baseline-only model and in each per-query metadata file for query models.

Completed query scans write `diagnostics/query_status.tsv` in requested order, including failed queries; strict failures publish this ledger without new scientific results.
Root metadata records the following runtime fields, available to successful per-query model readers through the shared metadata reference.
Baseline-only runs omit these query-worker fields and the status table.

| Field | Meaning |
| --- | --- |
| `query_workers_requested` | Original nonzero integer `workers` request, including a negative CPU offset if supplied. |
| `query_workers_effective` | Resolved process count after the shared worker policy and query/batch work cap; 1 means one worker. |
| `query_worker_native_threads` | The fixed target of 1 BLAS thread per worker, including serial execution; runtime logs separately report verification. |

These fields describe execution, not scientific filtering or jackknife construction. Worker completion order does not change status order or the stable summary/manifest order. Source: `RegressionRunner.estimate_partitioned_h2_batch()` in [_regression/runner.py](../../../src/ldsc/_regression/runner.py); see [query workers and memory](../regression-options.md#43-query-workers-and-memory).

The root coefficient-delete file is present for a baseline-only fitted model. In the cell-type regime, each written per-query directory contains the delete values for that complete baseline-plus-query fit. Metadata records the block count, fitted annotation order, and relative file path. These float64 matrices are required by `quantile-h2` and must never be concatenated across query runs.

### `quantile-h2`

```text
quantile-h2/
  quantile_h2.tsv
  standardized_coefficients.tsv
  diagnostics/
    metadata.json
    quantile-h2.log
    snp_alignment_issues.tsv.gz
```

The root tables are scientific outputs. Diagnostic metadata records selected-model provenance, inherited common-MAF definition, target/reference sources, quantile and missing-token policies, and statistic definitions. The alignment table is always created with a stable schema; row-addressable exclusions and fatal identity/MAF problems are reported there.

### `rg`

```text
rg/
  rg.tsv
  rg_full.tsv
  h2_per_trait.tsv
  diagnostics/
    metadata.json
    rg.log
    pairs/manifest.tsv
    pairs/<pair>/metadata.json
    pairs/<pair>/rg_full.tsv
```

The `diagnostics/pairs/` tree is present only when per-pair detail output is
requested. `--output-dir` is required.

`rg.tsv`, `rg_full.tsv`, and per-pair `rg_full.tsv` report only nominal
two-sided p-values from each genetic-correlation fit. They do not include
package-computed corrected p-value columns; downstream users choose any
multiple-testing correction appropriate for their analysis.

Root `diagnostics/metadata.json` records `trait_names` (the labels after any `label@<suffix>` relabeling, in input order; plotting reads them), `sumstats_files` (each trait's source path, aligned with `trait_names`), `ldscore_dir`, `pair_kind`, `n_pairs` and `n_traits`.

Root and per-pair `diagnostics/metadata.json` files are provenance only.
Per-pair metadata records the prevalences applied to each trait (`samp_prev_1`, `pop_prev_1`, `samp_prev_2`, `pop_prev_2`, each `null` when unset) and a `scale` (`observed` | `liability`) field; the genetic correlation itself is scale-invariant.

### Plot result

```text
<source-result>/plots/
  <fixed-plot-name>.png
  diagnostics/
    metadata.json
    plot.log
```

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `artifact_type` | `plot_result`. | Identifies the derived family. |
| `plot_kind` | Stable plot identifier selected from source metadata. | Human and programmatic provenance. |
| `source_artifact_type` | Artifact type of the source result. | Provenance. |
| `source_result_dir` | Source result root used for dispatch. | Provenance; not a relocatable input contract. |
| `source_table` | Source-relative numerical table plotted. | Provenance. |
| `files.plot` | Fixed PNG filename relative to the plot root. | Locates the figure. |
| `uncertainty` | Uncertainty convention used by the selected plot. | Interpretation. |
| `heritability_annotations` | `rg` plots only: source-relative `source_table` (null if undeclared), `source_available`, `scale=observed`, `estimate_column=total_h2_obs`, `se_column=total_h2_obs_se`, `uncertainty=block_jackknife_standard_error`, and `missing_label=failed`. A declared missing file retains its path and has `source_available=false`. | Provenance and interpretation of the heatmap diagonal or anchor annotations; see [the plotting contract](../../dev/design/plotting.md#heritability-annotations-in-rg-plots). |
| `created_at` | UTC creation timestamp. | Audit only. |

Plot dispatch into the source result checks only its plotting-relevant metadata fields and declared files.
It does not check `schema_version`.
The `rg` heritability source is optional for plotting; absent sources or unusable estimate–SE pairs display `failed` without refitting.

### `h2` scale conversion result

```text
<h2-result>/postprocessing/liability-scale/
  h2_scale_conversion.tsv
  h2_prevalence_sensitivity.png       # sensitivity mode only
  diagnostics/
    metadata.json
    convert-h2-scale.log
```

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `artifact_type` | `h2_scale_conversion_result`. | Identifies the derived family. |
| `source_artifact_type` | Always `h2_result`. | Provenance. |
| `source_result_dir`, `source_table` | Source result and source-relative `h2` summary. | Provenance. |
| `mode` | `exact` or `sensitivity`. | Interprets row count and optional plot. |
| `samp_prev` | Sample case fraction `P`. | Reproducible conversion assumption. |
| `pop_prev` | Exact population prevalence `K`, otherwise `null`. | Exact-mode assumption. |
| `pop_prev_range` | Inclusive `[MIN, MAX]`, otherwise `null`. | Sensitivity-mode grid. |
| `num_points` | Number of conversion rows. | Grid definition. |
| `files.table` | `h2_scale_conversion.tsv`. | Locates the derived table. |
| `files.plot` | Sensitivity PNG in range mode; absent in exact mode. | Locates the optional figure. |
| `uncertainty` | `block_jackknife_standard_error`. | Interpretation. |
| `population_prevalence_uncertainty_propagated` | Always `false`. | Prevents treating the ribbon as uncertainty in `K`. |
| `created_at` | UTC creation timestamp. | Audit only. |

## Regression Output Rule

`h2`, `partitioned-h2`, and `rg` require `--output-dir` and write their public
result files plus diagnostics. Their numerical `RegressionRunner` methods remain
in-memory APIs.

## Compatibility-Critical Fields

Fields that LDSC3 downstream workflows treat as runtime compatibility-critical:

- Munged summary-statistics `sumstats.parquet` footer: `ldsc:artifact_type`,
  `ldsc:snp_identifier`, `ldsc:genome_build`, and `ldsc:trait_name`.
- LD-score `metadata.json`: `artifact_type`,
  `snp_identifier`, `genome_build`, `files`, `baseline_columns`,
  `query_columns`, and `counts[].column`.
- R² matrix per-file schema/header metadata when present:
  `ldsc:artifact_type`, `ldsc:snp_identifier`,
  `ldsc:genome_build`, and `ldsc:sorted_by_build`.

Everything else in JSON metadata is provenance, navigation, or reporting data.
