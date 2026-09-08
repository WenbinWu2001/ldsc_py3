# Artifact Metadata Field Inventory

Last updated on: 2026-09-07

Downstream identity metadata lives in the `sumstats.parquet` footer (for munged
sumstats) and in `ldscore/metadata.json` (for LD scores). Any metadata emitted by
annotate, ref-panel, h2, partitioned-h2, rg, query-level, or pair-level outputs is
diagnostic provenance for regression. The post-fit `quantile-h2` command is the
explicit exception: it consumes one selected partitioned-h2 model's metadata and
coefficient-delete artifact because those files identify the fitted model.

This document inventories metadata, schema metadata, and audit sidecars written
by the current LDSC workflows. The goal is to keep downstream contracts obvious:
only reusable input artifacts carry required identity metadata. Diagnostic
metadata lives under `diagnostics/`; a post-fit command may consume it only when
the selected result directory itself is its declared input contract.

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
included in the owned family for ordinary result-directory workflows. The
gene-index builder is the deliberate exception while running: its open log,
history, and lock live in hidden sibling `.<index-name>.build-state/`. Once the
handler closes, the successful log moves into the published index diagnostics.

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
carries no embedded metadata. Legacy `.sumstats` and `.sumstats.gz` inputs are
marked for rsID-to-panel projection. Footerless Parquet is rejected because it
is neither a self-describing LDSC3 artifact nor an LDSC2 compatibility format.

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

### Regression legacy-sumstats audit

When any regression input is legacy LDSC2 text and an output directory is
supplied, `h2`, `partitioned-h2`, and `rg` write
`diagnostics/dropped_snps/legacy_sumstats.tsv.gz`, including a header-only file
for a clean projection. Its fixed fields are `trait_name`, `source_path`,
`SNP`, `A1`, `A2`, `reason`, and `panel_candidate_count`. It is diagnostic only;
the in-memory projected table, not this sidecar, enters regression.

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
| --- | --- | --- |
| `artifact_type` | Must be `annotation_projection`. | None. |
| `files` | Relative map of query annotation shards and diagnostics. | None. |
| `snp_identifier` | SNP identity mode used when writing projected annotations. | None for diagnostic JSON. Query shard columns are runtime data. |
| `genome_build` | Genome-build assumption for projected annotations. | None for diagnostic JSON. |
| `baseline_annot_sources` | Source annotation paths/tokens. | None. |
| `query_annot_bed_sources` | Source BED paths/tokens. | None. |
| `chromosomes` | Chromosomes written. | None. |
| `query_columns` | Projected query annotation columns. | None. |
| `baseline_columns` | Baseline annotation columns used as the SNP template. | None. |
| `n_snps` | Number of retained SNP rows after annotation identity cleanup. | None. |
| `padding_bp` | BED or gene interval padding applied before projection. | None. |

### `ldscore`

```text
ldscore/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet
  ldscore.overlap.parquet
  diagnostics/
    ldscore.log
    query_annotation_status.tsv
    gene_list_audit.tsv.gz
    gene_list_resolution_summary.tsv
```

`ldscore.query.parquet` is present only when query LD scores are written.
`query_annotation_status.tsv` is present for BED/gene-list query runs;
The row-complete audit and per-source summary are present for gene-list runs.
`ldscore.overlap.parquet` holds the annotation overlap matrix (long form:
`row_annotation`, `col_annotation`, `overlap_all_snps`, `overlap_common_snps`)
that `partitioned-h2` requires. It is written only for runs with two or more
annotation columns; a single-annotation (e.g. base-only) run omits it.

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
| `count_config` | Common-SNP count settings, including the actual threshold operator (`>=` for native LDSC3 computation; strict `>` for LDSC2 conversion). | Required to interpret count universes and checked against overlap metadata. |
| `overlap_config` | Overlap-matrix provenance: `total_all_reference_snps`, `total_common_reference_snps`, `common_maf_min`, `common_maf_operator`, `stored_block`. `null` for single-annotation runs that write no overlap matrix. | Provides `M_tot` and the universe definition for overlap-aware partitioned-h2. |
| `annotation_types` | Per-column `binary` or `quantitative` classification. | Interpretation and logging only; never changes fitting. |
| `annotation_fingerprints` | SHA256 hashes for the ordered common reference-SNP universe and each fitted annotation's canonical float32 values. | Lets `quantile-h2` verify resupplied sources without storing annotation matrices. |
| `n_baseline_rows` | Number of rows in the baseline parquet table. | Reporting. |
| `n_query_rows` | Number of rows in the query parquet table, or zero. | Reporting. |
| `row_group_layout` | Row-group strategy. | Reporting/technical provenance. |
| `baseline_row_groups` | Row-group metadata for `ldscore.baseline.parquet`. | Reporting/technical provenance. |
| `query_row_groups` | Row-group metadata for `ldscore.query.parquet`, or `null`. | Reporting/technical provenance. |
| `gene_list_resolution_policy` | `strict` or `resolved-only` for gene-list runs. | Scientific provenance; ignored by regression. |
| `gene_list_resolution_counts` | Aggregate submitted, rejected, and unique-resolved counts. | Reproducibility/diagnostics; ignored by regression. |
| `query_diagnostics` | Relative paths to query status and, for gene runs, the audit and source summary. | Troubleshooting/navigation; ignored by regression. |
| `legacy_ldsc2_import` | Present only for explicit LDSC2 conversion: profile, selected source directories/prefixes/files, streaming SHA-256 hashes, rsID intersection counts, count origins, strict common-frequency rule, coordinate evidence, and diagnostic paths. | Provenance and compatibility auditing; regression still consumes the ordinary canonical fields. |

Native LDSC3 computation records an inclusive common operator (`>=`). Explicit
LDSC2 conversion records the actual strict legacy operator (`>`) at threshold
`0.05`; `overlap_config` must agree with `count_config`, and the loader rejects
drift. Converted suites use the same canonical directory layout. Their workflow
log is `diagnostics/convert-ldsc2-ldscores.log` and their always-written issue
audit is `diagnostics/conversion_issues.tsv.gz`; they do not also claim a native
`diagnostics/ldscore.log`.

`diagnostics/query_annotation_status.tsv` has fixed columns `query`, `source`,
`input_type`, `status`, `reason`, `n_annotation_snps`, and `details`.
`diagnostics/gene_list_audit.tsv.gz` is row-complete and
`diagnostics/gene_list_resolution_summary.tsv` has one row per focal/control
source. Their fixed schemas and null rules are documented in
[Gene-list diagnostics and repair](gene-list-diagnostics-and-repair.md).

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
rank-bin summaries used by `ldsc plot`; it is a required current h2 artifact.

| Field | Explanation | Downstream usage |
| --- | --- | --- |
| `artifact_type` | Must be `h2_result`. | Selects h2 plotting and conversion contracts. |
| `files` | Relative map with `summary: "h2.tsv"` and `ld_score_regression_bins: "diagnostics/ld_score_regression_bins.tsv"`. | Conversion follows `summary`; plotting follows `ld_score_regression_bins`. |
| `trait_name` | Trait label resolved from CLI input or sumstats metadata. | None. |
| `sumstats_file` | Source sumstats path. | None. |
| `ldscore_dir` | Source LD-score directory. | None. |
| `effective_snp_identifier` | SNP identity mode actually used after any allowed identity downgrade. | None. |
| `genome_build` | Regression dataset genome build when available. | None. |
| `identity_downgrade_applied` | Whether same-family identity downgrade was used. | None. |
| `count_key_used_for_regression` | Count vector used for regression. | None. |
| `retained_ld_columns` | LD-score columns retained in the fitted h2 model. | None. |
| `dropped_zero_variance_ld_columns` | LD-score columns removed before fitting. | None. |
| `n_snps` | Number of SNPs retained after the effective chi-square filter. | None. |
| `effective_chisq_max` | Chi-square cap actually applied (`null` when uncapped). | None. |
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
    query_annotations/manifest.tsv
    query_annotations/<query>/metadata.json
    query_annotations/<query>/partitioned_h2.tsv
    query_annotations/<query>/partitioned_h2_full.tsv
    query_annotations/<query>/coefficient_delete_values.parquet
```

The `diagnostics/query_annotations/` tree is present for every query-annotation
run and absent for baseline-only runs. `--output-dir` is required; the deprecated
`--write-per-query-results` flag is an accepted no-op.

Root and per-query `diagnostics/metadata.json` files are provenance only. The
root file is self-describing about the analysis: `analysis_type`
(`functional_category` | `cell_type_specific`), `headline_metric` (`enrichment` |
`coefficient`), `enrichment_p_test` (`two_sided_t`), and `coefficient_p_test`
(`one_sided_greater`).

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

The root tables are scientific outputs. Diagnostic metadata records selected-model provenance, inherited common-MAF definition, target/reference sources, quantile and missing-token policies, verification level, and statistic definitions. The alignment table is always created with a stable schema; row-addressable exclusions and fatal identity/MAF problems are reported there.

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

Root and per-pair `diagnostics/metadata.json` files are provenance only. Per-pair
metadata records the prevalences applied to each trait (`samp_prev_1`, `pop_prev_1`,
`samp_prev_2`, `pop_prev_2`, each `null` when unset) and a `scale`
(`observed` | `liability`) field; the rg ratio itself is scale-invariant.

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
| `created_at` | UTC creation timestamp. | Audit only. |

Plot dispatch into the source result checks only its plotting-relevant metadata
fields and declared file. It does not check `schema_version`.

### h2 scale conversion result

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
| `source_result_dir`, `source_table` | Source result and source-relative h2 summary. | Provenance. |
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
