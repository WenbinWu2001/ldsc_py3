# Canonical Artifact metadata.json Scheme

**Date:** 2026-05-16
**Scope:** Public LDSC workflow outputs and their JSON metadata files.

---

## Problem

Current and recent outputs mixed `manifest.json`, `sumstats.metadata.json`,
`h2.metadata.json`, root diagnostic `metadata.json`, nested diagnostic
`metadata.json`, and old `format` fields. That made it hard for users to know
which files are required for downstream analysis and which files are diagnostic
provenance.

## Decision

Only `sumstats/metadata.json` and `ldscore/metadata.json` are downstream
contracts. Downstream LDSC workflows validate and consume those files.

All other metadata is diagnostic provenance only. Diagnostic metadata is written
under `diagnostics/metadata.json` or below a diagnostics detail tree. It must
not be required to run downstream analysis, and downstream workflows must not
change behavior when it is present.

The legacy top-level `format` field is removed from all JSON metadata.
`schema_version` plus `artifact_type` is the discriminator.

There is no backward compatibility for old metadata names or old root
diagnostic metadata paths. Package-written artifacts must be regenerated under
the current scheme.

## Boundary Rule

There is no "optional but used if present" metadata.

- Required metadata is loaded and validated by downstream code.
- Diagnostic metadata is for debugging, audit, and external tooling only.

If a future workflow needs metadata to run, that metadata must be promoted to a
documented required contract rather than opportunistically consumed.

## Required For Downstream

Only these metadata files are required for downstream LDSC analysis:

- `sumstats/metadata.json`
- `ldscore/metadata.json`

Regression workflows use them to recover SNP identity mode, genome-build
assumptions, trait labels, LD-score parquet paths, LD-score columns, and
annotation count records.

## Provenance Only

These metadata files are diagnostic only:

- `annotate/diagnostics/metadata.json`
- `ref-panel/diagnostics/metadata.json`
- `h2/diagnostics/metadata.json`
- `partitioned-h2/diagnostics/metadata.json`
- `partitioned-h2/diagnostics/query_annotations/*/metadata.json`
- `rg/diagnostics/metadata.json`
- `rg/diagnostics/pairs/*/metadata.json`

## Public Layouts

### `annotate`

```text
annotate/
  query.<chrom>.annot.gz
  diagnostics/
    metadata.json
    annotate.log
    dropped_snps/dropped.tsv.gz
```

Diagnostic metadata fields:

- `schema_version`
- `artifact_type: "annotation_projection"`
- `files`
- `snp_identifier`
- `genome_build`
- `baseline_annot_sources`
- `query_annot_bed_sources`
- `chromosomes`
- `query_columns`

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
the legacy TSV artifact. `metadata.json` and `sumstats.parquet` are
downstream-required for parquet output.

Required metadata fields:

- `schema_version`
- `artifact_type: "sumstats"`
- `files`
- `snp_identifier`
- `genome_build`
- `trait_name`

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

The concrete log file may be `build-ref-panel.log` or
`build-ref-panel.chr<chrom>.log`, depending on the builder entry point.

Diagnostic metadata fields:

- `schema_version`
- `artifact_type: "ref_panel"`
- `files`
- `snp_identifier`
- `source_genome_build`
- `emitted_genome_builds`
- `chromosomes`

Per-file validation metadata remains in parquet Arrow schema metadata and
`chr<chrom>_meta.tsv.gz` `# ldsc:*` header comments. The diagnostics root
metadata is not a runtime prerequisite.

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

Required metadata fields:

- `schema_version`
- `artifact_type: "ldscore"`
- `files`
- `snp_identifier`
- `genome_build`
- `chromosomes`
- `baseline_columns`
- `query_columns`
- `counts`
- `count_config`
- `n_baseline_rows`
- `n_query_rows`
- `row_group_layout`
- `baseline_row_groups`
- `query_row_groups`

### `h2`

```text
h2/
  h2.tsv
  diagnostics/
    metadata.json
    h2.log
```

Diagnostic metadata fields:

- `schema_version`
- `artifact_type: "h2_result"`
- `files`
- `trait_name`
- `sumstats_file`
- `ldscore_dir`
- `effective_snp_identifier`
- `genome_build`
- `identity_downgrade_applied`
- `count_key_used_for_regression`
- `retained_ld_columns`
- `dropped_zero_variance_ld_columns`
- `n_snps`

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
output is requested.

Root diagnostic metadata fields:

- `schema_version`
- `artifact_type: "partitioned_h2_result"`
- `files`
- `trait_name`
- `sumstats_file`
- `ldscore_dir`
- `count_kind`
- `query_annotations`

Per-query diagnostic metadata fields:

- `schema_version`
- `artifact_type: "partitioned_h2_query_result"`
- `files`
- `ordinal`
- `query_annotation`
- `slug`
- `folder`
- retained and dropped LD-score columns
- effective SNP identifier
- identity downgrade status
- SNP count

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
requested.

Root diagnostic metadata fields:

- `schema_version`
- `artifact_type: "rg_result"`
- `files`
- `ldscore_dir`
- `sumstats_files`
- `trait_names`
- `pair_kind`

Per-pair diagnostic metadata fields:

- `schema_version`
- `artifact_type: "rg_pair_result"`
- `files`
- `ordinal`
- `trait_1`
- `trait_2`
- `slug`
- `folder`
- status and error
- retained and dropped LD-score columns
- effective SNP identifier
- identity downgrade status
- intercept policies
- block count
- SNP count

## Regression No-Output Rule

Regression commands behave consistently:

- With `--output-dir`, write public result files plus `diagnostics/`.
- Without `--output-dir`, print the compact public TSV table to stdout and
  write no output files or diagnostics.

This applies to `h2`, `partitioned-h2`, and `rg`.

## Validation Rules

Only metadata for artifacts consumed as downstream inputs is validated as a
runtime prerequisite:

- `sumstats/metadata.json` must validate as `artifact_type == "sumstats"`.
- `ldscore/metadata.json` must validate as `artifact_type == "ldscore"`.

For these artifacts, `schema_version` must equal the current package schema
version and `artifact_type` must match the expected artifact type. If validation
fails, callers should report that the artifact must be regenerated with the
current LDSC package.

Diagnostic metadata for result directories, annotation outputs, and ref-panel
directory summaries is not validated by downstream workflows.

No JSON metadata file writes or reads a top-level `format` field.

## Implementation Notes

The migration is direct and intentionally has no fallback:

- `manifest.json` remains removed; LD-score uses `metadata.json`.
- `sumstats.metadata.json` remains removed; sumstats uses `metadata.json`.
- `h2.metadata.json` remains removed; h2 diagnostic metadata uses
  `diagnostics/metadata.json`.
- Diagnostic metadata for annotate, ref-panel, h2, partitioned-h2, and rg moves
  under `diagnostics/`.
- Per-query and per-pair detail trees move under `diagnostics/`.
- Loader tests must construct only canonical `sumstats/metadata.json` and
  `ldscore/metadata.json` contracts.
- Output tests must assert old root metadata names are not written.
- JSON metadata tests must assert that no top-level `format` field is present.
