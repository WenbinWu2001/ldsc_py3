# Artifact Metadata Field Inventory

metadata.json is a downstream contract only for sumstats and ldscore. Downstream workflows validate and consume sumstats/metadata.json and ldscore/metadata.json. Any metadata emitted by annotate, ref-panel, h2, partitioned-h2, rg, query-level, or pair-level outputs is diagnostic provenance only and must not be required to run downstream analysis.

This document inventories metadata, schema metadata, and audit sidecars written
by the current LDSC workflows. The goal is to keep downstream contracts obvious:
only input artifacts required by another LDSC command have required root
`metadata.json` files. Diagnostic metadata lives under `diagnostics/`.

## Boundary Rule

`metadata.json` has exactly two roles:

- **Required contract:** downstream LDSC code validates and consumes it. This
  role is limited to `sumstats/metadata.json` and `ldscore/metadata.json`.
- **Diagnostic provenance:** humans and external tools may inspect it for
  reproducibility/debugging. LDSC downstream workflows must not load it, require
  it, or change behavior when it is present.

There is no "optional but used if present" metadata. No JSON metadata file
contains a top-level `format`; `schema_version` plus `artifact_type` is the
schema discriminator.

## Public Output Layouts

### `munge-sumstats`

```text
sumstats/
  metadata.json
  sumstats.parquet
  sumstats.sumstats.gz
  diagnostics/
    sumstats.log
    dropped_snps/dropped.tsv.gz
```

`sumstats.sumstats.gz` is present only when the selected output format writes
the legacy TSV artifact.

`metadata.json` is downstream-required.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `schema_version` | Current LDSC artifact identity schema version. | Required and validated. |
| `artifact_type` | Must be `sumstats`. | Required and validated. |
| `files` | Relative data-file map for `sumstats.parquet` and/or `sumstats.sumstats.gz`. | Required to locate the data artifact. |
| `snp_identifier` | SNP identity mode used when the munged artifact was written. | Required for regression compatibility checks. |
| `genome_build` | Genome build associated with coordinate-family identity modes. | Required for regression compatibility checks when known. |
| `trait_name` | Trait label from the raw sumstats config, when supplied. | Used for output labels, not compatibility. |

### Sumstats `diagnostics/dropped_snps/dropped.tsv.gz`

Always-owned audit sidecar for liftover and identity-cleanup drops. This file is
written even for clean runs so users can distinguish "no drops" from "missing
sidecar".

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `CHR` | Chromosome for the dropped row when available. | None. |
| `SNP` | SNP label for the dropped row when available. | None. |
| `source_pos` | Source coordinate before liftover or cleanup. | None. |
| `target_pos` | Target coordinate after liftover when available. | None. |
| `reason` | Drop reason. | None. |
| `base_key` | Allele-blind identity key involved in cleanup. | None. |
| `identity_key` | Effective allele-aware identity key when applicable. | None. |
| `allele_set` | Normalized allele set used by allele-aware cleanup. | None. |
| `stage` | Workflow stage that produced the drop. | None. |

### `build-ref-panel`

```text
ref-panel/
  <build>/chr<chrom>_r2.parquet
  <build>/chr<chrom>_meta.tsv.gz
  diagnostics/
    metadata.json
    build-ref-panel.log
    build-ref-panel.chr<chrom>.log
    dropped_snps/chr<chrom>_dropped.tsv.gz
```

The concrete log file is `build-ref-panel.log` for multi-chromosome runs or
`build-ref-panel.chr<chrom>.log` for chromosome-scoped runs.

`diagnostics/metadata.json` is provenance only.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `schema_version` | Diagnostic metadata schema version. | None. |
| `artifact_type` | Must be `ref_panel`. | None. |
| `files` | Relative map of emitted R2, metadata, and diagnostic sidecars. | None. |
| `snp_identifier` | SNP identity mode used for emitted panel artifacts. | None for root metadata. Per-file schema/header metadata carries runtime validation. |
| `source_genome_build` | Build of source PLINK coordinates. | None for root metadata. |
| `emitted_genome_builds` | Builds emitted by the builder. | None. |
| `chromosomes` | Chromosomes emitted by the run. | None. |

### Ref-panel parquet and sidecar metadata

Ref-panel runtime compatibility is carried by per-file metadata, not the root
diagnostic JSON.

| Field | Location | Downstream usage |
| --- | --- | --- |
| `ldsc:schema_version` | `chr<chrom>_r2.parquet` Arrow schema and `chr<chrom>_meta.tsv.gz` header. | Validated when package metadata is present. |
| `ldsc:artifact_type` | Same. Must be `ref_panel_r2` or `ref_panel_metadata`. | Validated when package metadata is present. |
| `ldsc:snp_identifier` | Same. | Compared with runtime identity mode. |
| `ldsc:genome_build` | Same. | Compared with runtime genome build when known. |
| `ldsc:sorted_by_build` | R2 parquet schema metadata. | Used to infer/validate R2 coordinate build. |
| `ldsc:n_samples` | R2 parquet schema metadata. | Used to auto-fill sample size for raw R2 values. |
| `ldsc:r2_bias` | R2 parquet schema metadata. | Used to resolve R2 bias-correction behavior. |

### `annotate`

```text
annotate/
  query.<chrom>.annot.gz
  diagnostics/
    metadata.json
    annotate.log
    dropped_snps/dropped.tsv.gz
```

`diagnostics/metadata.json` is provenance only.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `schema_version` | Diagnostic metadata schema version. | None. |
| `artifact_type` | Must be `annotation_projection`. | None. |
| `files` | Relative map of query annotation shards and diagnostics. | None. |
| `snp_identifier` | SNP identity mode used when writing projected annotations. | None for diagnostic JSON. Query shard columns are runtime data. |
| `genome_build` | Genome-build assumption for projected annotations. | None for diagnostic JSON. |
| `baseline_annot_sources` | Source annotation paths/tokens. | None. |
| `query_annot_bed_sources` | Source BED paths/tokens. | None. |
| `chromosomes` | Chromosomes written. | None. |
| `query_columns` | Projected query annotation columns. | None. |

### `ldscore`

```text
ldscore/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet
  diagnostics/
    ldscore.log
```

`ldscore.query.parquet` is present only when query LD scores are written.

`metadata.json` is downstream-required.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `schema_version` | Current LDSC artifact identity schema version. | Required and validated. |
| `artifact_type` | Must be `ldscore`. | Required and validated. |
| `files` | Relative data-file map with required `baseline` and optional `query`. | Required to locate parquet data. |
| `snp_identifier` | SNP identity mode used for LD-score rows. | Required for regression compatibility checks. |
| `genome_build` | Genome build associated with LD-score identity metadata. | Required for regression compatibility checks when known. |
| `chromosomes` | Chromosome labels represented in the output. | Reporting/navigation. |
| `baseline_columns` | Ordered baseline annotation LD-score columns. | Required to assemble regression covariates. |
| `query_columns` | Ordered query annotation LD-score columns. | Required for partitioned h2 query selection. |
| `counts` | Per-annotation count records. | Required for regression count vectors. |
| `count_config` | Common-SNP count settings. | Reporting/context. |
| `n_baseline_rows` | Number of rows in the baseline parquet table. | Reporting. |
| `n_query_rows` | Number of rows in the query parquet table, or zero. | Reporting. |
| `row_group_layout` | Row-group strategy. | Reporting/technical provenance. |
| `baseline_row_groups` | Row-group metadata for `ldscore.baseline.parquet`. | Reporting/technical provenance. |
| `query_row_groups` | Row-group metadata for `ldscore.query.parquet`, or `null`. | Reporting/technical provenance. |

### `h2`

```text
h2/
  h2.tsv
  diagnostics/
    metadata.json
    h2.log
```

`diagnostics/metadata.json` is provenance only. If `--output-dir` is omitted,
the CLI prints compact TSV output to stdout and writes no diagnostics.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `schema_version` | Diagnostic metadata schema version. | None. |
| `artifact_type` | Must be `h2_result`. | None. |
| `files` | Relative map with `summary: "h2.tsv"`. | None. |
| `trait_name` | Trait label resolved from CLI input or sumstats metadata. | None. |
| `sumstats_file` | Source sumstats path. | None. |
| `ldscore_dir` | Source LD-score directory. | None. |
| `effective_snp_identifier` | SNP identity mode actually used after any allowed identity downgrade. | None. |
| `genome_build` | Regression dataset genome build when available. | None. |
| `identity_downgrade_applied` | Whether same-family identity downgrade was used. | None. |
| `count_key_used_for_regression` | Count vector used for regression. | None. |
| `retained_ld_columns` | LD-score columns retained in the fitted h2 model. | None. |
| `dropped_zero_variance_ld_columns` | LD-score columns removed before fitting. | None. |
| `n_snps` | Number of SNPs in the fitted regression dataset. | None. |

### `partitioned-h2`

```text
partitioned-h2/
  partitioned_h2.tsv
  diagnostics/
    metadata.json
    partitioned-h2.log
    query_annotations/manifest.tsv
    query_annotations/<query>/metadata.json
    query_annotations/<query>/partitioned_h2.tsv
    query_annotations/<query>/partitioned_h2_full.tsv
```

The `diagnostics/query_annotations/` tree is present only when per-query detail
output is requested. If `--output-dir` is omitted, the CLI prints compact TSV
output to stdout and writes no diagnostics.

Root and per-query `diagnostics/metadata.json` files are provenance only.

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
requested. If `--output-dir` is omitted, the CLI prints compact `rg.tsv` output
to stdout and writes no diagnostics.

Root and per-pair `diagnostics/metadata.json` files are provenance only.

## Regression No-Output Rule

Regression commands share one user-facing rule:

- With `--output-dir`, write the public result files plus diagnostics.
- Without `--output-dir`, print the compact public TSV table to stdout and
  write no files or diagnostics.

This applies to `h2`, `partitioned-h2`, and `rg`.

## Compatibility-Critical Fields

Fields that LDSC downstream workflows treat as runtime compatibility-critical:

- Sumstats `metadata.json`: `schema_version`, `artifact_type`,
  `snp_identifier`, `genome_build`, and `files`.
- LD-score `metadata.json`: `schema_version`, `artifact_type`,
  `snp_identifier`, `genome_build`, `files`, `baseline_columns`,
  `query_columns`, and `counts[].column`.
- Ref-panel per-file schema/header metadata when present:
  `ldsc:schema_version`, `ldsc:artifact_type`, `ldsc:snp_identifier`,
  `ldsc:genome_build`, and `ldsc:sorted_by_build`.

Everything else in JSON metadata is provenance, navigation, or reporting data.
