# Artifact Metadata Field Inventory

Downstream identity metadata lives in the `sumstats.parquet` footer (for munged
sumstats) and in `ldscore/metadata.json` (for LD scores). Any metadata emitted by
annotate, ref-panel, h2, partitioned-h2, rg, query-level, or pair-level outputs is
diagnostic provenance only and must not be required to run downstream analysis.

This document inventories metadata, schema metadata, and audit sidecars written
by the current LDSC workflows. The goal is to keep downstream contracts obvious:
only input artifacts required by another LDSC command carry required identity
metadata. Diagnostic metadata lives under `diagnostics/`.

## Boundary Rule

Downstream-required identity metadata has exactly two carriers:

- The **`sumstats.parquet` footer** (discrete `ldsc:*` keys) for munged sumstats.
- **`ldscore/metadata.json`** for LD-score directories.

Everything else is **diagnostic provenance:** humans and external tools may
inspect it for reproducibility/debugging, but LDSC downstream workflows must not
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
invocation; for `build-ref-panel`, concrete chromosome runs own only that
chromosome's files, while full `@` suite runs own the all-chromosome package.

Legacy/root diagnostic names that are no longer in the public layout are not
owned by current workflows. They are ignored by preflight and cleanup rather
than blocked or deleted.

Directory artifacts are owned as whole trees. Optional
`diagnostics/query_annotations/` and `diagnostics/pairs/` directories are staged
under `diagnostics/` and then moved into final position as a unit; reruns that
omit those optional details remove the old tree after a successful overwrite.
Always-written audit sidecars such as workflow logs and dropped-SNP tables are
included in the owned family.

## Public Output Layouts

### `munge-sumstats`

```text
sumstats/
  sumstats.parquet
  sumstats.sumstats.gz
  diagnostics/
    sumstats.log
    dropped_snps/dropped.tsv.gz
```

`sumstats.sumstats.gz` is present only when the selected output format writes
the legacy TSV artifact.

The downstream-required identity metadata is embedded in the `sumstats.parquet`
**footer** as discrete `ldsc:*` keys; no `metadata.json` is written. `.sumstats.gz`
carries no embedded metadata. An input without footer metadata (legacy `.sumstats.gz`
or footer-less parquet) loads with an inferred identifier mode rather than being
rejected.

| Footer key | Explanation | Downstream usage |
| --- | --- | --- |
| `ldsc:artifact_type` | Must be `sumstats`. | Required and validated. |
| `ldsc:snp_identifier` | SNP identity mode used when the munged artifact was written. | Recovered for regression compatibility checks. |
| `ldsc:genome_build` | Final output genome build for coordinate-family identity modes; empty (`null`) for rsID-family modes because rsID identity is build-independent. | Recovered for the regression genome-build check when known. |
| `ldsc:trait_name` | Trait label from the raw sumstats config, when supplied; empty otherwise. | Used for output labels, not compatibility. |

### Sumstats `diagnostics/dropped_snps/dropped.tsv.gz`

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
    metadata.chr<chrom>.json
    build-ref-panel.log
    build-ref-panel.chr<chrom>.log
    dropped_snps/chr<chrom>_dropped.tsv.gz
```

The concrete diagnostics files are `metadata.json` and `build-ref-panel.log`
for multi-chromosome runs. Concrete single-chromosome runs write
`metadata.chr<chrom>.json` and `build-ref-panel.chr<chrom>.log` so parallel
per-chromosome invocations can share an output directory without racing on
shared diagnostics.

`diagnostics/metadata*.json` is provenance only.

| Field | Explanation | Downstream usage |
| --- | --- | --- || `artifact_type` | Must be `ref_panel`. | None. |
| `files` | Relative map of emitted R2, metadata, and diagnostic sidecars. | None. |
| `snp_identifier` | SNP identity mode used for emitted panel artifacts. | None for root metadata. Per-file schema/header metadata carries runtime validation. |
| `source_genome_build` | Build of source PLINK coordinates. | None for root metadata. |
| `emitted_genome_builds` | Builds emitted by the builder. | None. |
| `chromosomes` | Chromosomes emitted by the run. | None. |

### Ref-panel parquet and sidecar metadata

Ref-panel runtime compatibility is carried by per-file metadata, not the root
diagnostic JSON.

| Field | Location | Downstream usage |
| --- | --- | --- || `ldsc:artifact_type` | Same. Must be `ref_panel_r2` or `ref_panel_metadata`. | Validated when package metadata is present. |
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
| --- | --- | --- || `artifact_type` | Must be `annotation_projection`. | None. |
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
  ldscore.overlap.parquet
  diagnostics/
    ldscore.log
```

`ldscore.query.parquet` is present only when query LD scores are written.
`ldscore.overlap.parquet` holds the annotation overlap matrix (long form:
`row_annotation`, `col_annotation`, `overlap_all_snps`, `overlap_common_snps`)
that `partitioned-h2` requires.

`metadata.json` is downstream-required.

| Field | Explanation | Downstream usage |
| --- | --- | --- || `artifact_type` | Must be `ldscore`. | Required and validated. |
| `files` | Relative data-file map with required `baseline` and optional `query` and `overlap`. | Required to locate parquet data; `overlap` is required by `partitioned-h2`. |
| `snp_identifier` | SNP identity mode used for LD-score rows. | Required for regression compatibility checks. |
| `genome_build` | Genome build associated with LD-score identity metadata. | Required for regression compatibility checks when known. |
| `chromosomes` | Chromosome labels represented in the output. | Reporting/navigation. |
| `baseline_columns` | Ordered baseline annotation LD-score columns. | Required to assemble regression covariates. |
| `query_columns` | Ordered query annotation LD-score columns. | Required for partitioned h2 query selection. |
| `counts` | Per-annotation count records. | Required for regression count vectors. |
| `count_config` | Common-SNP count settings (`common_reference_snp_maf_min`, `common_reference_snp_maf_operator: ">"`). | Reporting/context. |
| `overlap_config` | Overlap-matrix provenance: `total_all_reference_snps`, `total_common_reference_snps`, `common_maf_min`, `common_maf_operator`, `stored_block`. | Provides `M_tot` and the universe definition for overlap-aware partitioned-h2. |
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
| --- | --- | --- || `artifact_type` | Must be `h2_result`. | None. |
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

Root and per-query `diagnostics/metadata.json` files are provenance only. The
root file is self-describing about the analysis: `analysis_type`
(`functional_category` | `cell_type_specific`), `headline_metric` (`enrichment` |
`coefficient`), `enrichment_p_test` (`two_sided_t`), and `coefficient_p_test`
(`one_sided_greater`).

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

`rg.tsv`, `rg_full.tsv`, and per-pair `rg_full.tsv` report only nominal
two-sided p-values from each genetic-correlation fit. They do not include
package-computed corrected p-value columns; downstream users choose any
multiple-testing correction appropriate for their analysis.

Root and per-pair `diagnostics/metadata.json` files are provenance only.

## Regression No-Output Rule

Regression commands share one user-facing rule:

- With `--output-dir`, write the public result files plus diagnostics.
- Without `--output-dir`, print the compact public TSV table to stdout and
  write no files or diagnostics.

This applies to `h2`, `partitioned-h2`, and `rg`.

## Compatibility-Critical Fields

Fields that LDSC downstream workflows treat as runtime compatibility-critical:

- Sumstats `sumstats.parquet` footer: `ldsc:artifact_type`,
  `ldsc:snp_identifier`, `ldsc:genome_build`, and `ldsc:trait_name`.
- LD-score `metadata.json`: `artifact_type`,
  `snp_identifier`, `genome_build`, `files`, `baseline_columns`,
  `query_columns`, and `counts[].column`.
- Ref-panel per-file schema/header metadata when present:
  `ldsc:artifact_type`, `ldsc:snp_identifier`,
  `ldsc:genome_build`, and `ldsc:sorted_by_build`.

Everything else in JSON metadata is provenance, navigation, or reporting data.
