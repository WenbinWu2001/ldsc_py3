# Artifact Metadata Field Inventory

This document inventories metadata, manifest, schema-metadata, and audit-sidecar fields written by the current LDSC workflows for:

- `munge-sumstats`
- `build-ref-panel`
- `annotation-builder`
- `ldscore`
- regression modules

The goal is to distinguish fields that define a downstream compatibility contract from fields that are only provenance, audit, navigation, or reporting metadata.

## Compatibility Legend

| Value | Meaning |
| --- | --- |
| Yes | Downstream code reads this field and may reject the artifact, choose a merge mode, choose a runtime behavior, or require the field to proceed. |
| Structural | The field is required to locate or assemble files/tables, but it is not provenance compatibility. |
| Indirect | The field itself is not validated for value compatibility, but its presence helps identify a package-written artifact or a supported format. |
| No | The field is recorded for audit, reporting, navigation, or user labeling only. |

## Shared Identity Metadata

Several artifacts repeat the same minimal identity metadata contract produced by `identity_artifact_metadata()`:

| Field | Explanation | Usage |
| --- | --- | --- |
| `schema_version` | Current LDSC artifact identity schema version. | Compatibility: validates that an artifact was written under the current identity contract. |
| `artifact_type` | Artifact family marker such as `sumstats`, `ref_panel_r2`, `ref_panel_metadata`, or `ldscore`. | Compatibility: validates that the metadata belongs to the expected artifact type. |
| `snp_identifier` | Public SNP identity mode: `rsid`, `rsid_allele_aware`, `chr_pos`, or `chr_pos_allele_aware`. | Compatibility: reconstructs artifact identity mode and rejects incompatible combinations. |
| `genome_build` | Genome-build assumption associated with coordinate-family identity modes. May be `null` for rsID-family workflows. | Compatibility: reconstructs config provenance and rejects incompatible coordinate-build combinations when both sides are known. |

This quartet is repeated, but most downstream compatibility checks depend on it.

## `munge-sumstats`

### `sumstats.metadata.json`

Written beside curated sumstats artifacts such as `sumstats.parquet` and/or `sumstats.sumstats.gz`.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `schema_version` | Identity metadata schema version. | Yes | Yes. `load_sumstats()` requires this through identity metadata validation. |
| `artifact_type` | Must be `sumstats`. | Yes | Yes. `load_sumstats()` rejects artifacts whose type is not `sumstats`. |
| `snp_identifier` | SNP identity mode used when the munged artifact was written. | Yes | Yes. Reconstructs a `GlobalConfig` snapshot for sumstats identity and later regression compatibility. |
| `genome_build` | Genome build implied by the effective munging coordinate provenance. | Yes | Yes. Required in the sidecar and compared with LD-score provenance when both are known. |
| `trait_name` | Optional trait label from the raw sumstats config. | Yes | No. Used for labels and regression output naming, not compatibility. |

### `dropped_snps/dropped.tsv.gz`

Always-owned audit sidecar for liftover and identity-cleanup drops. This sidecar is written even for clean runs so users can distinguish "no drops" from "missing sidecar".

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome for the dropped row when available. | Yes | No. Audit/debugging only. |
| `SNP` | SNP label for the dropped row when available. | Yes | No. Audit/debugging only. |
| `source_pos` | Source coordinate before liftover or cleanup. | Yes | No. Audit/debugging only. |
| `target_pos` | Target coordinate after liftover when available. | Yes | No. Audit/debugging only. |
| `reason` | Drop reason, for example invalid identity, duplicate identity, unmapped liftover, or duplicate target coordinate. | Yes | No. Audit/debugging only. |
| `base_key` | Allele-blind identity key involved in cleanup. | Yes | No. Audit/debugging only. |
| `identity_key` | Effective allele-aware identity key when applicable. | Yes | No. Audit/debugging only. |
| `allele_set` | Normalized allele set used by allele-aware identity cleanup. | Yes | No. Audit/debugging only. |
| `stage` | Workflow stage that produced the drop, such as liftover or identity cleanup. | Yes | No. Audit/debugging only. |

## `build-ref-panel`

### `chr{chrom}_r2.parquet` Arrow schema metadata

Written into package-built canonical R2 parquet files.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `ldsc:schema_version` | Identity metadata schema version. | Yes | Yes. If package schema is detected, the reader validates this through identity metadata validation. |
| `ldsc:artifact_type` | Must be `ref_panel_r2`. | Yes | Yes. Validated before using package-written R2 identity metadata. |
| `ldsc:snp_identifier` | SNP identity mode used when the R2 artifact was built. | Yes | Yes. Compared with runtime `GlobalConfig.snp_identifier`. |
| `ldsc:genome_build` | Build associated with the R2 artifact identity metadata. | Yes | Yes. Compared with runtime `GlobalConfig.genome_build`. |
| `ldsc:sorted_by_build` | Build whose positions are used by `POS_1` and `POS_2`. | Yes | Yes. Used to infer or validate parquet R2 genome build. Conflicts raise. |
| `ldsc:row_group_size` | Requested parquet row-group size used during R2 writing. | Yes, technical | Indirect. Presence helps identify package-written schema, but the value is not used for compatibility. |
| `ldsc:n_samples` | LD-reference sample size from the PLINK genotype reader. | Yes | No compatibility check. Used to auto-fill sample size when R2 values are treated as raw. |
| `ldsc:r2_bias` | Whether R2 values are `raw` or already `unbiased`. Package-built panels write `unbiased`. | Yes | Runtime behavior, not compatibility. Used to resolve R2 bias correction mode. |

### `chr{chrom}_r2.parquet` table columns

These are data columns rather than sidecar fields, but several behave like format-contract fields for downstream readers.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome for pairwise R2 row. | No | Structural. Required by canonical R2 reader. |
| `POS_1` | Position of first endpoint in `ldsc:sorted_by_build`. | No | Structural. Required for row-group pruning and R2 lookup. |
| `POS_2` | Position of second endpoint in `ldsc:sorted_by_build`. | No | Structural. Required for R2 lookup. |
| `SNP_1` | SNP label for first endpoint. | No | Structural. Used in rsID-family matching. |
| `SNP_2` | SNP label for second endpoint. | No | Structural. Used in rsID-family matching. |
| `A1_1` | Allele 1 for first endpoint. | No | Structural in allele-aware modes. |
| `A2_1` | Allele 2 for first endpoint. | No | Structural in allele-aware modes. |
| `A1_2` | Allele 1 for second endpoint. | No | Structural in allele-aware modes. |
| `A2_2` | Allele 2 for second endpoint. | No | Structural in allele-aware modes. |
| `R2` | Pairwise LD value. | No | Structural. Core runtime data. |

### `chr{chrom}_meta.tsv.gz` header metadata

Written as `# ldsc:<key>=<value>` comments before the runtime metadata table.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `ldsc:schema_version` | Identity metadata schema version. | Yes | Yes when sidecar identity metadata is present or required. |
| `ldsc:artifact_type` | Must be `ref_panel_metadata`. | Yes | Yes. Validated before using metadata identity fields. |
| `ldsc:snp_identifier` | SNP identity mode used when metadata was written. | Yes | Yes. Compared with runtime `GlobalConfig.snp_identifier`. |
| `ldsc:genome_build` | Build for the metadata coordinates. | Yes | Yes. Compared with runtime `GlobalConfig.genome_build`. |

### `chr{chrom}_meta.tsv.gz` table columns

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome for retained SNP row. | No | Structural. Used for matching, filtering, and coordinate identity. |
| `POS` | Position in the emitted build. | No | Structural. Used for chr-pos matching and window construction. |
| `SNP` | SNP label. | No | Structural. Used for rsID matching. |
| `A1` | First allele. | No | Structural in allele-aware modes. |
| `A2` | Second allele. | No | Structural in allele-aware modes. |
| `CM` | Genetic-map coordinate, nullable when no map was provided. | No | Runtime data. Used for cM windows when requested. |
| `MAF` | Minor allele frequency from the reference panel. | No | Runtime data. Used for MAF filters and common-SNP count logic. |

### `dropped_snps/chr{chrom}_dropped.tsv.gz`

Always-owned per-chromosome audit sidecar for source identity cleanup and liftover-related drops.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome for the dropped SNP. | Yes | No. Audit/debugging only. |
| `SNP` | SNP label for the dropped row. | Yes | No. Audit/debugging only. |
| `source_pos` | Position before liftover or collision handling. | Yes | No. Audit/debugging only. |
| `target_pos` | Position after liftover when available. | Yes | No. Audit/debugging only. |
| `reason` | Drop reason. | Yes | No. Audit/debugging only. |
| `base_key` | Allele-blind identity key involved in cleanup. | Yes | No. Audit/debugging only. |
| `identity_key` | Effective identity key when applicable. | Yes | No. Audit/debugging only. |
| `allele_set` | Normalized allele set used by allele-aware cleanup. | Yes | No. Audit/debugging only. |
| `stage` | Stage that produced the drop. | Yes | No. Audit/debugging only. |

## `annotation-builder`

### `query.<chrom>.annot.gz`

The annotation writer materializes projected BED query annotations as ordinary LDSC annotation shards. It does not write a manifest or identity metadata sidecar for these query shards.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome. | No | Structural. Parsed as annotation row metadata. |
| `POS` | Base-pair position. | No | Structural. Used for chr-pos identity and row alignment. |
| `SNP` | SNP label. | No | Structural. Used for rsID identity and row alignment. |
| `CM` | Genetic-map coordinate. | No | Runtime data for cM windows if used downstream. |
| `A1` | Optional allele 1. | No | Structural in allele-aware modes when present. |
| `A2` | Optional allele 2. | No | Structural in allele-aware modes when present. |
| `MAF` | Optional allele frequency metadata if present in source annotation. | No | Runtime data. May support MAF filtering/counts in LD-score workflows. |
| Query annotation columns | One column per projected BED query annotation. | No | Structural/runtime data. Consumed as annotation matrix columns. |

### In-memory `AnnotationBundle` provenance

`AnnotationBundle` is not a file manifest, but it is the metadata object shared with LD-score and partitioned-regression workflows.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `metadata` | Row metadata table with `CHR`, `POS`, `SNP`, `CM`, and optional `A1`, `A2`, `MAF`. | No | Structural. Used for reference SNP universe and row alignment. |
| `baseline_annotations` | Baseline annotation matrix aligned to `metadata`. | No | Structural/runtime data. Used in LD-score calculation. |
| `query_annotations` | Query annotation matrix aligned to `metadata`. | No | Structural/runtime data. Used in LD-score calculation and partitioned h2. |
| `baseline_columns` | Ordered names of baseline annotation columns. | No | Structural. Used to label LD-score columns and counts. |
| `query_columns` | Ordered names of query annotation columns. | No | Structural. Used to label query LD-score columns and partitioned h2 output. |
| `chromosomes` | Chromosomes represented in the bundle. | Light provenance | Structural for sharded output and chromosome iteration. |
| `source_summary` | Input source paths/tokens used to construct the bundle. | Yes | No. Logging/reporting only. |
| `config_snapshot` | `GlobalConfig` active when the bundle was built. | Yes | Yes indirectly. Used by LD-score/regression workflows to preserve identity and genome-build assumptions. |

### `dropped_snps/dropped.tsv.gz`

Aggregate annotation identity-cleanup audit sidecar written by the BED-to-annotation output workflow.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome for the dropped row. | Yes | No. Audit/debugging only. |
| `SNP` | SNP label for the dropped row. | Yes | No. Audit/debugging only. |
| `source_pos` | Source position when available. | Yes | No. Audit/debugging only. |
| `target_pos` | Target position when available. | Yes | No. Audit/debugging only. |
| `reason` | Identity cleanup drop reason. | Yes | No. Audit/debugging only. |
| `base_key` | Allele-blind identity key. | Yes | No. Audit/debugging only. |
| `identity_key` | Effective identity key. | Yes | No. Audit/debugging only. |
| `allele_set` | Normalized allele set. | Yes | No. Audit/debugging only. |
| `stage` | Cleanup stage, typically annotation identity cleanup. | Yes | No. Audit/debugging only. |

## `ldscore`

### `manifest.json`

Written with `ldscore.baseline.parquet` and optional `ldscore.query.parquet`.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `format` | Directory-format marker, currently `ldsc.ldscore_result.v1`. | Yes, schema provenance | Yes. `load_ldscore_from_dir()` rejects unsupported formats. |
| `schema_version` | Identity metadata schema version. | Yes | Yes. Validated when reconstructing `GlobalConfig` from the manifest. |
| `artifact_type` | Must be `ldscore`. | Yes | Yes. Validated when reconstructing `GlobalConfig`. |
| `snp_identifier` | SNP identity mode used for LD-score rows. | Yes | Yes. Used to reconstruct config; explicit loader overrides must match. |
| `genome_build` | Genome build associated with LD-score identity metadata. | Yes | Yes. Compared with sumstats provenance during regression when both sides are known. |
| `files` | Relative data-file map, with required `baseline` and optional `query`. | No | Structural. Loader requires `files.baseline` and reads `files.query` if present. |
| `chromosomes` | Chromosome labels represented in the output. | Light provenance | No current compatibility check. |
| `baseline_columns` | Ordered baseline annotation LD-score columns. | No | Structural. Used to assemble regression LD-score table and count vectors. |
| `query_columns` | Ordered query annotation LD-score columns. | No | Structural. Used for partitioned h2 query selection and count vectors. |
| `counts` | List of per-annotation count records. | Yes, analysis provenance | Yes structurally. Regression requires records for retained LD-score columns. |
| `count_config` | Settings used to define common-SNP counts, such as MAF threshold and operator. | Yes | No hard compatibility check. Used for reporting/context. |
| `n_baseline_rows` | Number of rows in the baseline parquet table. | Yes, summary provenance | No. |
| `n_query_rows` | Number of rows in the query parquet table, or zero. | Yes, summary provenance | No. |
| `row_group_layout` | Row-group strategy, currently one row group per chromosome. | Yes, technical | No current compatibility check. |
| `baseline_row_groups` | Row-group metadata for `ldscore.baseline.parquet`. | Yes, technical | No current compatibility check. |
| `query_row_groups` | Row-group metadata for `ldscore.query.parquet`, or `null` when absent. | Yes, technical | No current compatibility check. |

### `counts[]` records in `manifest.json`

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `group` | Whether the annotation column belongs to `baseline` or `query`. | Yes | No hard compatibility check. Helps explain count record origin. |
| `column` | Annotation column name that the count applies to. | No | Yes structurally. Regression aligns count records to retained LD-score columns by this field. |
| `all_reference_snp_count` | Number of reference SNPs in the annotation. | Yes, analysis provenance | Yes. Required for count vectors. |
| `common_reference_snp_count` | Optional number of common reference SNPs in the annotation. | Yes, analysis provenance | Yes when present for all retained columns. Regression prefers common counts for `count_kind=common`. |

### `baseline_row_groups[]` and `query_row_groups[]` records

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `chrom` | Chromosome represented by the row group. | Yes, technical | No current compatibility check. |
| `row_group_index` | Row-group index in the parquet file. | Yes, technical | No current compatibility check. |
| `row_offset` | Starting row offset in the logical table. | Yes, technical | No current compatibility check. |
| `n_rows` | Number of rows in the row group. | Yes, technical | No current compatibility check. |

### `ldscore.baseline.parquet` and `ldscore.query.parquet`

These are primary data tables. Their columns are not sidecar fields, but they enforce part of the output contract.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `CHR` | Chromosome. | No | Structural. Used in row alignment and SNP identity. |
| `POS` | Position. | No | Structural. Used in row alignment and chr-pos identity. |
| `SNP` | SNP label. | No | Structural. Used in row alignment and rsID identity. |
| `A1` | Optional allele 1, required for allele-aware LD-score artifacts. | No | Structural. Loader rejects allele-aware artifacts missing usable alleles. |
| `A2` | Optional allele 2, required for allele-aware LD-score artifacts. | No | Structural. Loader rejects allele-aware artifacts missing usable alleles. |
| `regression_ld_scores` | Weight LD-score column used for regression weights. | No | Runtime data. Required for regression dataset construction. |
| Baseline LD-score columns | One column per baseline annotation. | No | Runtime data. Used as regression covariates. |
| Query LD-score columns | One column per query annotation in `ldscore.query.parquet`. | No | Runtime data. Used for partitioned h2 query analyses. |

## Regression Modules

Regression output metadata is mostly result provenance and navigation. Compatibility checks happen before writing outputs, while building regression datasets from loaded sumstats and LD-score artifacts.

### Unpartitioned h2 `h2.tsv`

The unpartitioned h2 workflow writes a compact summary table, `h2.tsv`, and a
provenance sidecar, `h2.metadata.json`, when `output_dir` is supplied. Input
compatibility is enforced before writing, while loading and merging
`sumstats.metadata.json` and LD-score `manifest.json` provenance.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `trait_name` | Trait label resolved from CLI input or the sumstats metadata sidecar. | Yes | No. Reporting/display only. |
| `n_snps` | Number of SNPs in the final regression dataset after sumstats, LD scores, filters, and identity matching. | Yes | No. Reporting only. |
| `total_h2` | Estimated total SNP heritability from the fitted LDSC h2 model. | No | Result value, not compatibility metadata. |
| `total_h2_se` | Standard error for `total_h2`. | No | Result value, not compatibility metadata. |
| `intercept` | Estimated or fixed LDSC intercept used in the h2 model. | Yes, model provenance | No. Reporting/model interpretation only. |
| `intercept_se` | Standard error for the intercept when estimated; may be unavailable for fixed intercepts. | Yes, model provenance | No. Reporting/model interpretation only. |
| `mean_chisq` | Mean chi-square statistic among regression SNPs. | Yes, model summary provenance | No. Reporting/model interpretation only. |
| `lambda_gc` | Genomic-control lambda computed from regression SNP chi-square statistics. | Yes, model summary provenance | No. Reporting/model interpretation only. |
| `ratio` | LDSC ratio statistic derived from intercept and mean chi-square when defined. | Yes, model summary provenance | No. Reporting/model interpretation only. |
| `ratio_se` | Standard error for `ratio` when defined. | Yes, model summary provenance | No. Reporting/model interpretation only. |

### Unpartitioned h2 `h2.metadata.json`

Written beside `h2.tsv` only when the h2 workflow owns an output directory.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `schema_version` | h2 metadata schema version. | Yes | No current reader compatibility check found. |
| `artifact_type` | Result marker, currently `h2_result`. | Yes, schema provenance | No current reader compatibility check found. |
| `trait_name` | Trait label resolved from CLI input or the sumstats metadata sidecar. | Yes | No. Reporting/display only. |
| `sumstats_file` | Source sumstats file path from CLI args. | Yes | No. Audit only. |
| `ldscore_dir` | Source LD-score directory path from CLI args. | Yes | No. Audit only. |
| `effective_snp_identifier` | SNP identity mode actually used after any allowed identity downgrade. | Yes | No current downstream check. Useful provenance. |
| `genome_build` | Genome build from the regression dataset provenance when available. | Yes | No current downstream check. |
| `identity_downgrade_applied` | Whether same-family allele-aware/base identity downgrade was used. | Yes | No current downstream check. Useful provenance. |
| `count_key_used_for_regression` | Count vector used for regression, for example common or all reference SNP counts. | Yes | No current downstream check. Important run provenance. |
| `retained_ld_columns` | LD-score columns retained in the fitted h2 model. | Yes | No. Audit/reporting only. |
| `dropped_zero_variance_ld_columns` | LD-score columns removed before fitting because they had zero variance. | Yes | No. Audit/reporting only. |
| `n_snps` | Number of SNPs in the fitted regression dataset. | Yes | No. Reporting only. |

### Partitioned h2 `query_annotations/manifest.tsv`

Written only when per-query result output is enabled.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `ordinal` | Stable 1-based query order. | No | Structural/navigation only. |
| `query_annotation` | Original query annotation name. | Yes | No. Used for navigation and display. |
| `slug` | Filesystem-safe slug derived from `query_annotation`. | No | Structural/navigation only. |
| `folder` | Per-query result folder name. | No | Structural/navigation only. |
| `summary_path` | Relative path to per-query compact summary. | No | Structural/navigation only. |
| `partitioned_h2_full_path` | Relative path to per-query full summary. | No | Structural/navigation only. |
| `metadata_path` | Relative path to per-query metadata JSON. | No | Structural/navigation only. |

### Partitioned h2 `query_annotations/*/metadata.json`

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `format` | Result-family marker, currently `ldsc.partitioned_h2_result.v1`. | Yes, schema provenance | No current reader compatibility check found. |
| `trait_name` | Trait label for the partitioned h2 run. | Yes | No. Reporting only. |
| `count_kind` | Requested count kind, for example `common`. | Yes | No. Reporting only. |
| `ldscore_dir` | Source LD-score directory path from CLI args. | Yes | No. Reporting/audit only. |
| `dropped_zero_variance_ld_columns` | LD-score columns removed before fitting because they had zero variance. | Yes | No. Audit/reporting only. |
| `retained_ld_columns` | LD-score columns retained in the fitted model. | Yes | No. Audit/reporting only. |
| `n_snps` | Number of SNPs in the fitted regression dataset. | Yes | No. Reporting only. |
| `effective_snp_identifier` | SNP identity mode actually used after any allowed identity downgrade. | Yes | No current downstream check. Useful provenance. |
| `identity_downgrade_applied` | Whether same-family allele-aware/base identity downgrade was used. | Yes | No current downstream check. Useful provenance. |
| `ordinal` | Stable 1-based query order. | No | Navigation only. |
| `query_annotation` | Original query annotation name. | Yes | Navigation/display only. |
| `slug` | Filesystem-safe slug. | No | Navigation only. |
| `folder` | Per-query folder name. | No | Navigation only. |

### Genetic correlation `pairs/manifest.tsv`

Written only when per-pair detail output is enabled.

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `ordinal` | Stable 1-based pair order. | No | Structural/navigation only. |
| `trait_1` | First trait label. | Yes | No. Display/navigation only. |
| `trait_2` | Second trait label. | Yes | No. Display/navigation only. |
| `slug` | Filesystem-safe slug for the pair. | No | Structural/navigation only. |
| `folder` | Per-pair folder name. | No | Structural/navigation only. |
| `rg_full_path` | Relative path to per-pair full rg row. | No | Structural/navigation only. |
| `metadata_path` | Relative path to per-pair metadata JSON. | No | Structural/navigation only. |

### Genetic correlation `pairs/*/metadata.json`

| Field | Explanation | Provenance field? | Downstream compatibility usage |
| --- | --- | --- | --- |
| `format` | Result-family marker, currently `ldsc.rg_result_family.v1`. | Yes, schema provenance | No current reader compatibility check found. |
| `trait_1` | First trait label. | Yes | No. Reporting/display only. |
| `trait_2` | Second trait label. | Yes | No. Reporting/display only. |
| `source_1` | Source path for first sumstats table. | Yes | No. Audit only. |
| `source_2` | Source path for second sumstats table. | Yes | No. Audit only. |
| `pair_kind` | Whether this pair came from all-pairs or anchor mode. | Yes | No. Reporting only. |
| `status` | Pair status, such as `ok` or `failed`. | Yes | No compatibility check. Used to describe result state. |
| `error` | Error text when pair fitting failed. | Yes | No. Audit/debugging only. |
| `n_snps_used` | Number of SNPs used by the fitted pair. | Yes | No. Reporting only. |
| `n_blocks_used` | Number of jackknife blocks used. | Yes | No. Reporting only. |
| `count_key_used_for_regression` | Count vector used for regression, for example common or all reference SNP counts. | Yes | No current downstream check. Important run provenance. |
| `retained_ld_columns` | LD-score columns retained in the fitted pair. | Yes | No. Audit/reporting only. |
| `dropped_zero_variance_ld_columns` | LD-score columns removed before fitting. | Yes | No. Audit/reporting only. |
| `effective_snp_identifier` | SNP identity mode used after any allowed identity downgrade. | Yes | No current downstream check. Important run provenance. |
| `identity_downgrade_applied` | Whether identity downgrade was applied. | Yes | No current downstream check. Important run provenance. |
| `intercept_h2_policy` | Effective h2 intercept policy. | Yes | No. Reporting/audit only. |
| `intercept_gencov_policy` | Effective genetic-covariance intercept policy. | Yes | No. Reporting/audit only. |
| `ordinal` | Stable 1-based pair order. | No | Navigation only. |
| `slug` | Filesystem-safe pair slug. | No | Navigation only. |
| `folder` | Per-pair folder name. | No | Navigation only. |

Failed rg pairs omit fields that are only available after a regression dataset is built, such as retained LD columns and effective SNP identifier.

## Remarks

The identity quartet is repeated across reloadable artifacts, but it is the core compatibility contract.

Fields that should be treated as compatibility-critical:

- `schema_version`
- `artifact_type`
- `snp_identifier`
- `genome_build`
- LD-score `format`
- LD-score `files`
- LD-score `baseline_columns`, `query_columns`, and `counts[].column`
- Ref-panel R2 `ldsc:sorted_by_build`
- Ref-panel metadata/R2 identity metadata when present
