# Canonical Artifact metadata.json Scheme

**Date:** 2026-05-16
**Scope:** Public LDSC workflow output directories and their JSON metadata entry
points.

---

## Problem

Current public outputs use several naming and schema conventions for the same
concept:

- LD-score directories use `manifest.json`.
- Munged sumstats use `sumstats.metadata.json`.
- Unpartitioned h2 uses `h2.metadata.json`.
- Partitioned h2 and rg detail outputs use nested `metadata.json`.
- Some JSON payloads still write a legacy `format` field such as
  `ldsc.ldscore_result.v1`.

This makes output contracts harder to remember and pushes version/type
information into two competing mechanisms: `format` and
`schema_version`/`artifact_type`.

## Decision

Every public LDSC workflow output directory uses a root `metadata.json` file.
The root metadata file is the only JSON entry point for that output directory.

The legacy `format` field is removed from all JSON metadata. Artifact identity
is represented by:

- `schema_version`
- `artifact_type`

There is no backward compatibility for old metadata names or old `format`
payloads. Package-written artifacts must be regenerated under the current
scheme.

All root metadata files include:

- `schema_version`
- `artifact_type`
- `files`

Identity-bearing input artifacts also include:

- `snp_identifier`
- `genome_build`

Result metadata records provenance only: source paths, trait names, effective
SNP identifier, count key, retained/dropped LD-score columns, row counts, and
SNP counts as applicable.

## Required For Downstream

Only two metadata files are required for downstream analysis:

- `sumstats/metadata.json`
- `ldscore/metadata.json`

Regression workflows require these metadata files to recover the SNP identity
mode, genome-build assumptions, trait labels, LD-score parquet paths, LD-score
columns, and annotation count records needed to assemble regression datasets.

## Provenance Only

The following metadata files are not required for downstream workflow execution:

- `annotate/metadata.json`
- `ref-panel/metadata.json`
- `h2/metadata.json`
- `partitioned-h2/metadata.json`
- `partitioned-h2/query_annotations/*/metadata.json`
- `rg/metadata.json`
- `rg/pairs/*/metadata.json`

These files are written for auditability, navigation, and reproducibility. They
must not become hidden prerequisites for downstream LDSC commands.

## Public Layouts

### `annotate`

```text
annotate/
  metadata.json
  query.<chrom>.annot.gz
  dropped_snps/dropped.tsv.gz
  annotate.log
```

Required root metadata fields:

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
  sumstats.log
  dropped_snps/dropped.tsv.gz
```

`sumstats.sumstats.gz` is present only when the selected output format writes
the legacy TSV artifact.

Required root metadata fields:

- `schema_version`
- `artifact_type: "sumstats"`
- `files`
- `snp_identifier`
- `genome_build`
- `trait_name`

This metadata file is required for downstream regression.

### `build-ref-panel`

```text
ref-panel/
  metadata.json
  <build>/chr<chrom>_r2.parquet
  <build>/chr<chrom>_meta.tsv.gz
  dropped_snps/chr<chrom>_dropped.tsv.gz
  build-ref-panel.log
```

The concrete log file may be `build-ref-panel.log` or
`build-ref-panel.chr<chrom>.log`, depending on the builder entry point.

Required root metadata fields:

- `schema_version`
- `artifact_type: "ref_panel"`
- `files`
- `snp_identifier`
- `source_genome_build`
- `emitted_genome_builds`
- `chromosomes`

Per-file validation metadata remains in the parquet Arrow schema and the
`chr<chrom>_meta.tsv.gz` `# ldsc:*` header comments. The root metadata file is
an index and provenance record, not a runtime prerequisite.

### `ldscore`

```text
ldscore/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet
  ldscore.log
```

`ldscore.query.parquet` is present only when query LD scores are written.

Required root metadata fields:

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

This metadata file is required for downstream regression.

### `h2`

```text
h2/
  metadata.json
  h2.tsv
  h2.log
```

Required root metadata fields:

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
  metadata.json
  partitioned_h2.tsv
  partitioned-h2.log
  query_annotations/manifest.tsv
  query_annotations/<query>/metadata.json
  query_annotations/<query>/partitioned_h2.tsv
  query_annotations/<query>/partitioned_h2_full.tsv
```

The `query_annotations/` tree is present only when per-query detail output is
requested.

Required root metadata fields:

- `schema_version`
- `artifact_type: "partitioned_h2_result"`
- `files`
- `trait_name`
- `sumstats_file`
- `ldscore_dir`
- `count_kind`
- `query_annotations`

Required per-query metadata fields:

- `schema_version`
- `artifact_type: "partitioned_h2_query_result"`
- `files`
- `ordinal`
- `query_annotation`
- `slug`
- `folder`
- result provenance copied from the fitted per-query model, such as retained and
  dropped LD-score columns, effective SNP identifier, identity downgrade status,
  and SNP count.

### `rg`

```text
rg/
  metadata.json
  rg.tsv
  rg_full.tsv
  h2_per_trait.tsv
  rg.log
  pairs/manifest.tsv
  pairs/<pair>/metadata.json
  pairs/<pair>/rg_full.tsv
```

The `pairs/` tree is present only when per-pair detail output is requested.

Required root metadata fields:

- `schema_version`
- `artifact_type: "rg_result"`
- `files`
- `ldscore_dir`
- `sumstats_files`
- `trait_names`
- `pair_kind`

Required per-pair metadata fields:

- `schema_version`
- `artifact_type: "rg_pair_result"`
- `files`
- `ordinal`
- `trait_1`
- `trait_2`
- `slug`
- `folder`
- pair provenance copied from the fitted pair, such as status, error, retained
  and dropped LD-score columns, effective SNP identifier, identity downgrade
  status, intercept policies, block count, and SNP count.

## Validation Rules

Only metadata for artifacts consumed as downstream inputs is validated as a
runtime prerequisite:

- `sumstats/metadata.json` must validate as `artifact_type == "sumstats"`.
- `ldscore/metadata.json` must validate as `artifact_type == "ldscore"`.

For these artifacts, `schema_version` must equal the current package schema
version and `artifact_type` must match the expected artifact type. If validation
fails, callers should report that the artifact must be regenerated with the
current LDSC package.

Metadata for result directories and annotation/ref-panel directory summaries is
not validated by downstream workflows. Those files are still schema-versioned so
humans and external tools can interpret them safely.

No JSON metadata file writes or reads a top-level `format` field.

## Implementation Notes

The migration is a direct rename with no compatibility fallback:

- `manifest.json` becomes `metadata.json`.
- `sumstats.metadata.json` becomes `metadata.json`.
- `h2.metadata.json` becomes `metadata.json`.
- Existing nested `metadata.json` files keep their name and gain
  `schema_version` plus `artifact_type` where missing.
- Loader tests must construct only canonical `metadata.json` files.
- Output tests must assert that old root metadata names are not written.
- JSON metadata tests must assert that no top-level `format` field is present.

