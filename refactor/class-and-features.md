# Class And Features

This document records the current class model and user-visible feature surface of the refactored package under `src/ldsc/`.

## 1. Design Style

The implemented code follows the hybrid model chosen earlier:

| Concern | Current pattern |
| --- | --- |
| validated config | frozen dataclass |
| aligned in-memory workflow state | frozen dataclass |
| orchestration | service/workflow class |
| numerical methods | internal kernel module |
| output writing | output manager plus artifact producers |

The package is CLI-first, but the public Python API is intentionally usable without touching `_kernel`.

## 2. Layer Model

| Layer | Public? | Location |
| --- | --- | --- |
| CLI dispatcher | yes | `ldsc.cli` |
| workflow/config/output modules | yes | `src/ldsc/*.py` |
| compute kernel | no | `src/ldsc/_kernel/*.py` |
| tests | n/a | `tests/` |

## 3. Public Configuration Classes

Defined in `ldsc.config`.

### `CommonConfig`

| Field | Meaning |
| --- | --- |
| `snp_identifier` | global identifier mode: `rsid` or `chr_pos` |
| `genome_build` | build context for `chr_pos` workflows |
| `global_snp_restriction_path` | global SNP-universe restriction applied across annotation, reference, and regression SNP sets |
| `log_level` | logging verbosity |
| `fail_on_missing_metadata` | strict metadata handling toggle |

### `AnnotationBuildConfig`

| Field | Meaning |
| --- | --- |
| `baseline_annotation_paths` | baseline `.annot` inputs |
| `query_annotation_paths` | query `.annot` inputs |
| `query_bed_paths` | BED inputs for projection workflows |
| `out_prefix` | output prefix |
| `batch_mode` | one output directory per BED source |
| `compression` | output compression |

### `RefPanelConfig`

| Field | Meaning |
| --- | --- |
| `backend` | `auto`, `plink`, or `parquet_r2` |
| `plink_prefix` / `plink_prefix_chr` | PLINK reference input |
| `parquet_r2_paths` / `parquet_r2_paths_chr` | parquet R2 inputs |
| `frequency_paths` / `frequency_paths_chr` | metadata / frequency sidecars |
| `r2_bias_mode` | raw or unbiased R2 handling |
| `r2_sample_size` | reference sample size for parquet corrections |

### `LDScoreConfig`

| Field | Meaning |
| --- | --- |
| `ld_wind_snps`, `ld_wind_kb`, `ld_wind_cm` | exactly one LD-window definition |
| `maf_min` | retained-SNP MAF threshold |
| `chunk_size` | PLINK block-processing chunk size |
| `compute_m5_50` | emit common-SNP count vector |
| `whole_chromosome_ok` | override the safety guard for windows that effectively span a whole chromosome |

Important behavior:

- `compute_m5_50` changes the count artifact only.
- LD scores are still computed on the retained SNP universe for the run.

### `MungeConfig`

Current wrapper around the legacy-compatible munging parameter surface.

Representative fields:

- `N`, `N_cas`, `N_con`
- `info_min`, `maf_min`, `n_min`, `nstudy_min`
- `out_prefix`, `chunk_size`, `merge_alleles_path`
- `signed_sumstats_spec`, `ignore_columns`, `no_alleles`, `a1_inc`, `keep_maf`, `daner`, `daner_n`

### `RegressionConfig`

| Field | Meaning |
| --- | --- |
| `n_blocks` | block-jackknife count |
| `use_m_5_50` | keep LDSC default of preferring `.M_5_50` |
| `use_intercept` | estimate or constrain intercept |
| `intercept_h2` | fixed h2 intercept |
| `intercept_gencov` | fixed covariance intercept |
| `two_step_cutoff` | two-step estimator threshold |
| `chisq_max` | chi-square outlier filter |
| `samp_prev`, `pop_prev` | liability-scale conversion inputs |

Implemented default:

- `use_m_5_50=True`

## 4. Annotation Feature

### Public module

- `ldsc.annotation_builder`

### Internal kernel

- `ldsc._kernel.annotation`

### Main classes

#### `AnnotationSourceSpec`

| Field | Meaning |
| --- | --- |
| `baseline_annot_paths` | baseline `.annot` paths |
| `query_annot_paths` | query `.annot` paths |
| `bed_paths` | BED paths |
| `gene_set_paths` | gene-set paths |
| `allow_missing_query` | whether empty query input is allowed |

#### `AnnotationBundle`

| Field | Meaning |
| --- | --- |
| `metadata` | aligned SNP metadata |
| `baseline_annotations` | baseline annotation matrix |
| `query_annotations` | query annotation matrix |
| `baseline_columns` | baseline annotation names |
| `query_columns` | query annotation names |
| `chromosomes` | chromosomes represented in the bundle |
| `source_summary` | provenance summary |

Methods:

- `validate()`
- `reference_snps()`
- `has_full_baseline_cover()`
- `annotation_matrix()`
- `summary()`

#### `AnnotationBuilder`

Main responsibilities:

- parse `.annot` files
- enforce row alignment across files
- apply global SNP restriction
- project BED intervals to SNP-level annotations

Main public helpers:

- `run()`
- `project_bed_annotations()`
- `run_bed_to_annot()`
- `make_annot_files()`

## 5. Reference-Panel Feature

### Public surface

Reference-panel classes are re-exported from `ldsc`.

### Internal kernel

- `ldsc._kernel.ref_panel`

### Main classes

#### `RefPanelSpec`

| Field | Meaning |
| --- | --- |
| `backend` | `plink` or `parquet_r2` |
| `bfile_prefix` | PLINK prefix |
| `r2_table_paths` | parquet R2 paths |
| `sample_size` | reference sample size |
| `maf_metadata_paths` | frequency / metadata sidecars |
| `chromosomes` | explicit chromosome set when needed |
| `genome_build` | build metadata |

#### `RefPanel`

Abstract interface methods:

- `available_chromosomes()`
- `load_metadata(chrom)`
- `build_reader(chrom, **kwargs)`
- `filter_to_snps(chrom, snps)`
- `summary()`

#### `PlinkRefPanel`

Current role:

- load `.bim` metadata lazily by chromosome
- apply global SNP restriction
- build PLINK-backed readers for LD-score calculation

#### `ParquetR2RefPanel`

Current role:

- load chromosome metadata for parquet R2 inputs
- apply global SNP restriction
- build sorted-R2 readers

#### `RefPanelLoader`

Role:

- choose and construct the correct backend object from config/spec

## 6. LD-Score Feature

### Public module

- `ldsc.ldscore_calculator`

### Internal kernel

- `ldsc._kernel.ldscore`

### Main classes

#### `ChromLDScoreResult`

| Field | Meaning |
| --- | --- |
| `chrom` | chromosome label |
| `reference_metadata` | metadata aligned to `ld_scores` |
| `ld_scores` | chromosome-specific LD-score table |
| `regression_metadata` | metadata aligned to `w_ld` |
| `w_ld` | chromosome-specific regression-weight LD-score table |
| `snp_count_totals` | chromosome-specific count vectors keyed by descriptive name |
| `baseline_columns` | baseline annotation names |
| `query_columns` | query annotation names |
| `reference_snps` | retained reference SNP IDs |
| `regression_snps` | retained regression SNP IDs |
| `output_paths` | written artifact paths when available |

Methods:

- `validate()`
- `to_ldscore_table()`
- `to_weight_table()`
- `summary()`

#### `LDScoreResult`

Aggregated cross-chromosome result used by downstream regression.

| Field | Meaning |
| --- | --- |
| `reference_metadata` | aggregated metadata aligned to `ld_scores` |
| `ld_scores` | aggregated LD-score table |
| `regression_metadata` | aggregated metadata aligned to `w_ld` |
| `w_ld` | aggregated regression-weight table |
| `snp_count_totals` | aggregated count vectors |
| `baseline_columns` | baseline annotation names |
| `query_columns` | query annotation names |
| `reference_snps` | union of retained reference SNPs |
| `regression_snps` | union of retained regression SNPs |
| `chromosome_results` | list of `ChromLDScoreResult` objects |
| `output_paths` | written artifact paths |

Methods:

- `validate()`
- `to_ldscore_table()`
- `to_weight_table()`
- `summary()`

#### `LDScoreCalculator`

Main responsibilities:

- iterate over chromosomes
- compute one chromosome at a time
- aggregate chromosome results
- send outputs to `OutputManager`

Main methods:

- `run()`
- `compute_chromosome()`
- `write_outputs()`

Important workflow behavior:

- the pipeline is chromosome-wise for computational efficiency
- a global regression-SNP input is partitioned internally by chromosome
- `M` and `M_5_50` are stored as named count arrays in `snp_count_totals`
- regression later defaults to the common-SNP count key when available

Count keys currently used:

- `all_reference_snp_counts`
- `common_reference_snp_counts_maf_gt_0_05`

## 7. Summary-Statistics Munging Feature

### Public module

- `ldsc.sumstats_munger`

### Internal kernel

- `ldsc._kernel.sumstats_munger`

### Main classes

#### `RawSumstatsSpec`

| Field | Meaning |
| --- | --- |
| `path` | raw summary-stat file |
| `compression` | compression mode or `auto` |
| `trait_name` | optional trait label |
| `column_hints` | explicit column-name hints |

#### `SumstatsTable`

| Field | Meaning |
| --- | --- |
| `data` | munged summary-stat table |
| `has_alleles` | whether `A1`/`A2` are present |
| `source_path` | original input path |
| `trait_name` | optional trait label |
| `provenance` | run/source metadata |

Methods:

- `validate()`
- `snp_identifiers()`
- `subset_to()`
- `align_to_metadata()`
- `summary()`

#### `MungeRunSummary`

| Field | Meaning |
| --- | --- |
| `n_input_rows` | raw input row count |
| `n_retained_rows` | retained row count |
| `drop_counts` | per-filter drop counts |
| `inferred_columns` | resolved column hints |
| `used_n_rule` | how sample size was determined |
| `output_paths` | written file paths |

#### `SumstatsMunger`

Main methods:

- `run()`
- `write_output()`
- `build_run_summary()`

Behavior:

- this intentionally follows the old munging workflow without adding new features
- the refactor goal here is cleaner structure, not altered munging semantics

## 8. Regression Feature

### Public module

- `ldsc.regression_runner`

### Internal kernel

- `ldsc._kernel.regression`

### Main classes

#### `RegressionDataset`

| Field | Meaning |
| --- | --- |
| `merged` | merged sumstats + LD-score + weight table |
| `ref_ld_columns` | all LD-score columns before variance filtering |
| `weight_column` | regression-weight column name |
| `reference_snp_count_totals` | available count vectors |
| `count_key_used_for_regression` | selected count vector key |
| `retained_ld_columns` | LD-score columns kept after zero-variance filtering |
| `dropped_zero_variance_ld_columns` | dropped LD-score columns |
| `trait_names` | trait labels |
| `chromosomes_aggregated` | chromosomes included in the aggregated LD-score result |

Method:

- `validate()`

This replaced the vague older mask-style exposure with explicit retained and dropped column names.

#### `RegressionRunner`

Main methods:

- `build_dataset()`
- `estimate_h2()`
- `estimate_partitioned_h2()`
- `estimate_partitioned_h2_batch()`
- `estimate_rg()`

File-driven entry helpers:

- `add_h2_arguments()`
- `add_partitioned_h2_arguments()`
- `add_rg_arguments()`
- `run_h2_from_args()`
- `run_partitioned_h2_from_args()`
- `run_rg_from_args()`

Implemented default:

- when multiple query annotations are provided, `estimate_partitioned_h2_batch()` loops over them one at a time and concatenates a clean summary table with one row per query annotation

This feature was not present in the older code path and is now explicit in the refactor.

## 9. Output Feature

### Public module

- `ldsc.outputs`

### Main classes

#### `OutputSpec`

Configuration only. It answers:

- where outputs go
- which built-in outputs are written
- whether per-chromosome artifacts are written
- compression and layout choices
- whether extra registered artifact names are enabled

Key fields:

- `out_prefix`
- `output_dir`
- `artifact_layout`
- `write_ldscore`
- `write_w_ld`
- `write_counts`
- `write_annotation_manifest`
- `write_per_chrom`
- `aggregate_across_chromosomes`
- `compression`
- `overwrite`
- `write_summary_json`
- `write_summary_tsv`
- `write_run_metadata`
- `enabled_artifacts`

#### `ArtifactConfig`

Advanced per-artifact options reserved for future plot/report expansion.

#### `RunSummary`

Small derived run summary used by writers and future post-processing.

#### `ArtifactProducer`

Extension interface for adding new summaries, plots, reports, or manifests without hardcoding more branches into `OutputManager`.

#### `OutputManager`

Coordinates output generation from result objects and registered producers.

#### `ResultWriter`

Serializes artifacts to disk.

#### `ResultFormatter`

Builds clean manifest and summary tables.

#### `PostProcessor`

Reserved extension point for richer downstream summaries and reports.

## 10. Kernel Estimator Classes

Defined in `ldsc._kernel.regression`.

| Class | Role |
| --- | --- |
| `Jackknife` | base jackknife utility |
| `LstsqJackknifeFast` | fast block jackknife |
| `LstsqJackknifeSlow` | slow delete-block jackknife |
| `RatioJackknife` | ratio-estimator jackknife |
| `IRWLS` | iterative reweighted least squares wrapper |
| `LD_Score_Regression` | estimator base class |
| `Hsq` | heritability estimator |
| `Gencov` | genetic covariance estimator |
| `RG` | combines two heritability fits plus one covariance fit into a genetic-correlation estimate |

## 11. Mapping From Prior Structure

| Current module | Replaces or consolidates |
| --- | --- |
| `ldsc.annotation_builder` | old root `make_annot.py` and `utils/run_bed_to_annot.py` wrappers plus annotation-loading logic |
| `ldsc.ldscore_calculator` | old `ldscore_workflow.py` public role and root `ldsc_new.py` orchestration |
| `ldsc.sumstats_munger` | old root `munge_sumstats.py` public role |
| `ldsc.regression_runner` | old `regression_workflow.py` public role and the file-driven wrapper part of the old regression path |
| `ldsc._kernel.ldscore` | consolidated reusable parts of the previous PLINK kernel and newer chromosome/parquet engine |
| `ldsc._kernel.regression` | consolidated estimator logic and reusable allele helpers from the old regression stack |

The old root wrapper scripts and the old `refactor/ldscore/` package are no longer part of the active codebase.
