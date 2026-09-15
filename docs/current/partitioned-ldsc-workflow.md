# Partitioned LDSC Workflow: Technical Reference

Last updated on: 2026-09-15

This document describes the refactored workflow for computing LD scores and
running h2, partitioned-h2, and rg regression from one canonical LD-score result
directory.

Continuous fitted annotations use the same LD-score regression. LDSC3 classifies finite annotation values as binary or quantitative for logs and metadata only; classification never changes coefficients. Quantitative rows retain legacy numerical proportion/enrichment summaries with an interpretation warning. After fitting, `ldsc quantile-h2` can project one complete fitted joint model onto target-value quantiles; see [Continuous-Annotation Quantile Heritability](continuous-annotation-quantile-h2.md).

## 1. Overview

The pipeline has two phases:

1. `ldsc ldscore --output-dir <ldscore_dir>` computes annotation LD scores and
   writes one self-contained directory.
2. Regression commands read that directory through `--ldscore-dir <ldscore_dir>`.

The public regression CLI no longer accepts fragmented LD-score artifacts:
`--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`, and
`--query-columns` are removed.

Gene-set analyses may construct the same canonical LD-score directory through
an explicit exact index. By default, each regression model contains the
supplied baseline plus one focal gene set. If LD-score calculation was given a
`--control-gene-list-file`, its `gene_control` column is also part of the fixed
baseline block. See [the indexed gene-list tutorial](../wiki/main-functionalities/ldscore.md).

## 2. LD-Score Result Directory

An LD-score run writes:

```text
<ldscore_dir>/
  metadata.json
  ldscore.baseline.parquet
  ldscore.query.parquet        # one query batch; numbered files for multiple batches; omitted without queries
  ldscore.overlap.parquet      # overlap matrix; omitted for single-annotation (e.g. base-only) runs; consumed by partitioned-h2
  diagnostics/
    ldscore.log
```

`metadata.json` is the downstream metadata contract. It contains:

- `artifact_type: "ldscore"`
- relative file paths for `baseline`, optional `query` or `query_batchNNNNN` entries, and optional `overlap`
- required ordered `query_batches` entries containing `file`, `query_columns`, and chromosome `row_groups`; empty for baseline-only output
- `snp_identifier`, `genome_build`, and processed `chromosomes`
- ordered `baseline_columns` and `query_columns`
- one count record per LD-score annotation column
- `count_config` and `overlap_config` (total SNP-universe sizes and the
  common-MAF threshold/operator used for the overlap matrix)
- `config_snapshot`
- basic row counts

`ldscore.overlap.parquet` stores the annotation overlap matrix `O = AᵀA` in long
form (`row_annotation`, `col_annotation`, `overlap_all_snps`,
`overlap_common_snps`). Only the baseline-rows block `A_Bᵀ·A` plus each query's
self-overlap are stored — enough to reconstruct every baseline-plus-one-query
model overlap matrix without keeping the full annotation matrix. `partitioned-h2`
requires this file; the shared `h2` collinearity guard reads it when present
(two or more LD-score columns) but does not require it; `rg` ignores it. The
common-SNP universe uses
`MAF >= common_maf_min` (inclusive), matching the `.M_5_50`-style common counts.

All scientific artifact paths inside the metadata are relative to `ldscore_dir`. A single query batch uses `ldscore.query.parquet`; multiple batches use `ldscore.query.batch00001.parquet` and subsequent ordinals. Each file covers the same genome-wide regression SNP rows, with chromosome row groups. Baseline values, counts, and overlap statistics are shared. Current readers require the manifest; regenerate older directories. Source: `LDScoreDirectoryWriter.write_batches` in [outputs.py](../../src/ldsc/outputs.py).

`ldscore.baseline.parquet` columns:

```text
CHR, SNP, POS, regression_ld_scores, <baseline LD-score columns...>
```

Columns in each query file:

```text
CHR, SNP, POS, <query LD-score columns...>
```

Query files duplicate SNP key columns intentionally, including available `A1`/`A2`. Opening a source validates all declared schemas and allele metadata. `LDScoreSource.read_queries` checks effective SNP-row alignment across selected files; regression checks selected query rows against the baseline. See [ldscore_source.py](../../src/ldsc/ldscore_source.py) and `load_ldscore_from_dir` in [regression_runner.py](../../src/ldsc/regression_runner.py).

The public `SNP` column in LD-score outputs is a carried label, not necessarily
a dbSNP rsID. For package-built parquet reference panels it comes from the
original PLINK `.bim` `SNP` field through the `chr*_meta.tsv.gz` sidecar. Thus a
panel built from `.bim` IDs such as `22:10684250:C:G` will show those strings in
synthetic unpartitioned LD-score rows. When explicit baseline annotations are
supplied, the LD-score workflow keeps the annotation file's `SNP` labels and
uses the reference panel only to filter/align rows. Public SNP identifier modes
are exactly `rsid`, `rsid_allele_aware`, `chr_pos`, and
`chr_pos_allele_aware`; the default is `chr_pos_allele_aware`. Mode names are
exact. `snp_identifier="rsid"` means "match on the literal `SNP` strings";
`snp_identifier="chr_pos"` uses private `CHR:POS` keys for matching and does
not rewrite the visible `SNP` column. The allele-aware variants add the
unordered `A1/A2` allele set to those base keys.

## 3. Inputs

### Baseline Annotations

Baseline `.annot[.gz]` files are passed to `ldsc ldscore` with
`--baseline-annot-sources`. They define the baseline annotation LD-score columns stored
in `ldscore.baseline.parquet`.

Example:

```text
resources/baseline_v1.2/baseline.@.annot.gz
```

For ordinary unpartitioned LD-score generation, `--baseline-annot-sources` may be
omitted when no query inputs are supplied. The public workflow then creates a
synthetic all-ones baseline annotation named exactly `base` over the retained
reference-panel metadata and writes the same canonical LD-score directory. This
is not a separate h2 workflow.

### Query Annotations

Query annotations are optional and become columns in the saved query batch files.

- `--query-annot-bed-sources`: BED intervals projected onto the baseline SNP universe one execution batch at a time, using private staging under the output directory.
- `--query-annot-sources`: pre-built query `.annot[.gz]` files.

For pre-built chromosome-sharded inputs, use one query annotation file per
chromosome containing all query annotation columns; each column is one query
annotation. This is also the format written by `ldsc annotate`:

```text
query.1.annot.gz   # CHR POS SNP CM query_A query_B ...
query.2.annot.gz   # CHR POS SNP CM query_A query_B ...
...
```

Pass the suite as `--query-annot-sources 'query.@.annot.gz'`. Do not pass one
chromosome-sharded suite per query, because `ldscore` accepts only one query
file for each chromosome in sharded mode. The query headers must be identical
across chromosome files.

These arguments are mutually exclusive and require explicit
`--baseline-annot-sources`. If users intentionally want to test query annotations
against an all-ones universe, they should create an explicit all-ones `base`
baseline annotation over the query annotation universe and run the partitioned
workflow with both baseline and query inputs.

For gene-set or pathway BEDs that need flanking sequence, `--padding-bp`
expands each query BED interval by the requested number of base pairs on both
sides before projection; starts are clipped at zero. The default is `0`, so BED
files are used exactly as supplied. Do not set this option for BED files already
expanded during upstream gene-set preparation, because that would double-count
the flank. The flag is valid only for live BED or gene-list queries; remove it
when using prebuilt query annotations or an exact gene index.

### Reference Panel

LD scores can be computed from:

- PLINK input through `--plink-prefix`
- package-built parquet R2 input through `--r2-dir`

In `--r2-dir` mode, matching `chr*_meta.tsv.gz` sidecars are discovered
automatically. The sidecar is mandatory for the canonical index-format R2
parquet — it defines the index space and provides MAF and cM metadata. A missing
sidecar is a hard error. The same R2 parquet serves all four identifier modes.
External R2 parquet formats are not supported.

### Summary Statistics

Regression commands still read munged summary statistics:

- `ldsc h2 --sumstats-file <file>`
- `ldsc partitioned-h2 --sumstats-file <file>`
- `ldsc rg --sumstats-sources <file-or-glob> <file-or-glob> [...]`

Current `ldsc munge-sumstats` outputs include canonical `CHR` and `POS` columns
beside `SNP`, `Z`, and `N`, write `sumstats.parquet` by default, embed the thin
compatibility payload in the parquet footer (`artifact_type`, `snp_identifier`,
`genome_build`, and optional `trait_name`), and write
`diagnostics/dropped_snps/dropped.tsv.gz` for row-level liftover-drop auditing.
Legacy LDSC2 `.sumstats` and `.sumstats.gz` files may be used directly and are
projected by rsID onto the canonical LD-score panel. Footerless Parquet must be
regenerated. Legacy partitioned LD-score suites are not read directly; only a
complete baseline suite may be explicitly converted, without query annotations.
That converted directory runs the baseline-only functional-category regime:
all imported baseline columns are fitted jointly. It is not a converted
cell-type/query analysis and cannot be extended with legacy query annotations.
In allele-aware modes, current sumstats artifacts require usable `A1/A2`. To
run without allele-aware SNP identity, set `--snp-identifier chr_pos` or
`--snp-identifier rsid` intentionally.

## 4. SNP Universes

The LD-score phase tracks these SNP sets:

| Symbol | Name | Meaning |
|--------|------|---------|
| B | Baseline annotation SNPs | Rows in loaded baseline annotation files, or retained reference-panel metadata rows for synthetic `base` unpartitioned runs |
| A | Raw reference panel SNPs | Rows in PLINK `.bim` or parquet metadata sidecar |
| A' | Prepared reference panel | `A` after an optional explicit reference-panel SNP restriction, retained-panel `maf_min`, and PLINK `keep_indivs_file`; equals A when absent |
| `ld_reference_snps` | LD computation universe | `B ∩ A'` |
| C | Regression SNP mask | Optional identity-only mask from explicit `regr_snps_file` or packaged HM3 regression SNP restriction |
| `ld_regression_snps` | Persisted row set | `B ∩ A' ∩ C`; equals `ld_reference_snps` when C is absent |

LD-score column counts in root metadata are computed over `ld_reference_snps`.
Persisted parquet rows are `ld_regression_snps`.
Explicit reference-panel and regression SNP restriction files are filters only:
duplicate restriction keys collapse to one retained key, and non-identity
columns such as `CM` or `MAF` are ignored. Optional frequency metadata fills
missing `CM`/`MAF`; duplicate frequency metadata identity clusters are dropped
entirely before fill, leaving those values missing unless annotation metadata
already supplied them.

## 5. Count Records

Root metadata stores counts as records keyed by column name, not by array
position:

```json
{
  "group": "baseline",
  "column": "base",
  "all_reference_snp_count": 100000.0,
  "common_reference_snp_count": 85000.0
}
```

Root `metadata.json` also records the threshold used to compute common counts:

```json
{
  "count_config": {
    "common_reference_snp_maf_min": 0.05,
    "common_reference_snp_maf_operator": ">="
  }
}
```

The `common_reference_snp_count` key is omitted when common-SNP counts are
unavailable. Regression falls back to `all_reference_snp_count` in that case.

## 6. Regression Behavior

`partitioned-h2 --threads N` controls concurrent complete query models, with default 1 preserving inline execution. Parsing and resolution use the same helpers as `ldscore` and `build-gene-ldscore-index`: positive N requests N workers, while negative values use CPU affinity with a machine CPU-count fallback. All requests are capped by query count and `--query-batch-size`, which still bounds loading. Choose positive N within your CPU allocation. Parallel processes borrow immutable numeric maps and use one native numerical thread each. Shared preparation and final publication stay in the coordinator. Baseline-only behavior and model-specific filtering, weights, SNP order, and jackknife blocks are unchanged. See [the effective CPU and memory policy](regression-configuration.md#43-query-workers-and-memory) and `RegressionRunner._fit_partitioned_query()`.

Query scans retain strict publication by default. Explicit `--continue-on-query-error` skips per-query model preparation, fitting, and result calculation exceptions and publishes successful models; all attempted queries are recorded in `diagnostics/query_status.tsv`. The log contains each failure's query/stage/traceback and available jackknife block/rank/support/genomic diagnostics. Shared input/output failures and scans without successful fits remain fatal. The same rule applies to direct and indexed LD scores, with unchanged filtering, weighting, and contiguous block construction. See [query failure handling](partitioned-h2-results.md#query-failures-and-continuation) and `RegressionRunner.estimate_partitioned_h2_batch()`.

For the tunable estimator parameters (intercept policy, two-step cutoff,
`chisq_max`, jackknife blocks, counts) and the per-command defaults, see
[`regression-configuration.md`](regression-configuration.md).

`h2` and `rg` use baseline LD scores only, even when query batch files exist.
They also use the embedded `regression_ld_scores` column from
`ldscore.baseline.parquet`; this is the historical `w_ld` LD score over the
regression SNP universe, not the final model-dependent regression weight.
Regression merges on the effective key for the resolved mode: `SNP` in `rsid`,
`SNP:<allele_set>` in `rsid_allele_aware`, `CHR:POS` in `chr_pos`, and
`CHR:POS:<allele_set>` in `chr_pos_allele_aware`. Base modes are fully
allele-blind; allele columns may be preserved for orientation, but they never
affect base-mode identity, duplicate filtering, retention, or drop reasons.
`--allow-identity-downgrade` is regression-only and permits same-family
allele-aware/base mixes to run under the base mode. rsID-family and
coordinate-family modes never mix.

**Outlier handling.** Extreme-chi-square SNPs are controlled differently by
annotation count, matching legacy LDSC. A **single-annotation** fit (plain `h2`
and `rg`) applies no chi-square cap and instead defaults the two-step cutoff to
`30` (when `--two-step-cutoff` is unset and the h2 intercept is free). A
**multi-annotation** fit (`partitioned-h2`, where the two-step estimator does
not apply) instead applies a default outlier cap of `max(0.001 · N.max(), 80)`
when `--chisq-max` is unset, dropping SNPs above it so a few extreme statistics
cannot dominate the regression. An explicit `--chisq-max` overrides the default
in either case. The cap is inclusive (`chi^2 <= cap`, per the `-max`
convention), and any drop is logged as `Removed N SNPs with chi^2 > C (...)`.
The diagnostics metadata records the **post-filter** SNP count (`n_snps`) and
the `effective_chisq_max` actually applied (the default for a partitioned model,
the explicit value, or `null` when uncapped), so the reported SNP set matches
the fitted set.

`partitioned-h2` produces **overlap-aware** category summaries (legacy
`--overlap-annot` math) and auto-detects one of two regimes from the LD-score
directory, with no user flag. It requires `ldscore.overlap.parquet`; a directory
produced by an older `ldsc ldscore` is rejected with a regenerate message.

- **Functional-category regime** — the directory has **no** query columns. A
  single joint fit of all baseline annotations yields one row per baseline
  category; the headline is `enrichment` (+ two-sided `enrichment_p`). This
  reproduces Finucane-2015 / legacy `--overlap-annot`. If only a single
  annotation is retained (e.g. a one-annotation baseline directory), the fit is
  degenerate and collapses to single-annotation `h2` behavior (two-step
  estimator, no default chi-square cap); the run logs a `WARNING` and you should
  supply multiple baseline annotations for a meaningful partitioned analysis.
- **Cell-type-specific regime** — the directory **has** query columns. For each
  query, a `baseline + one query` model is fit, yielding one row per query; the
  headline is `coefficient` (the conditional `tau`, + one-sided `coefficient_p`
  testing `coefficient > 0`). Every query run also writes a staged
  `diagnostics/query_annotations/` tree (`manifest.tsv` plus one sanitized
  folder per query with a one-row `partitioned_h2.tsv`, the full
  baseline-plus-query `partitioned_h2_full.tsv`, and `metadata.json`). The
  retired `--write-per-query-results` flag is rejected.

Both regimes write **one** column schema to `partitioned_h2.tsv`, differing only
in rows and the default sort:

```text
category, prop_snps, prop_h2, prop_h2_se, enrichment, enrichment_se,
enrichment_p, coefficient, coefficient_se, coefficient_z, coefficient_p,
overlap_annot, total_h2_obs, total_h2_obs_se, total_h2_liab, total_h2_liab_se,
category_h2_obs, category_h2_obs_se, category_h2_liab, category_h2_liab_se,
samp_prev, pop_prev
```

All columns are lowercase snake_case, consistent with `h2` and `rg`.

Interpretation. `enrichment = prop_h2 / prop_snps` is the **marginal**
heritability of the SNPs in a category (overlap-aware, baseline-confounded for a
query); `coefficient` is the **conditional** per-SNP contribution beyond the
baseline. `total_h2_obs` is the fitted model's total heritability (constant across
one model's rows; per-query in the cell-type regime) and is the reference total
for the category-level absolute heritabilities. `category_h2_obs = M_c·tau_c` is
the conditional category contribution (observed scale) and can be negative under
overlap; `category_h2_liab` is its liability-scale counterpart
(`= category_h2_obs · c(samp_prev, pop_prev)`, `NaN` without
`--samp-prev`/`--pop-prev`). Note that under overlap (almost always), `prop_h2` is
**not** `category_h2_obs / total_h2_obs` -- see the overlap caveat in
`partitioned-h2-results.md`. `overlap_annot` flags whether the fitted model's
annotations overlap. `--summary-sort-by` defaults to `auto` (→ `coefficient-p`
for cell-type, `category` for functional); the run logs a regime banner and
records `analysis_type` / `headline_metric` / `enrichment_p_test` /
`coefficient_p_test` in `diagnostics/metadata.json`.

`partitioned-h2` treats `partitioned_h2.tsv`, `diagnostics/metadata.json`,
`diagnostics/query_annotations/`, and `diagnostics/partitioned-h2.log` as one
owned output family. `--output-dir` is required. Without `--overwrite`, any
existing owned sibling rejects the run. With `--overwrite`, a successful
baseline-only run removes a stale `diagnostics/query_annotations/` tree from a
previous query configuration.

## 7. CLI Examples

Compute LD scores:

```bash
ldsc ldscore \
  --output-dir results/my_study_ldscore \
  --baseline-annot-sources resources/baseline_v1.2/baseline.@.annot.gz \
  --query-annot-bed-sources my_peaks.bed \
  --padding-bp 0 \
  --plink-prefix resources/1kg/1KG_EUR_Phase3_chr \
  --snp-identifier rsid \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

Compute ordinary unpartitioned LD scores without baseline annotations:

```bash
ldsc ldscore \
  --output-dir results/my_unpartitioned_ldscore \
  --plink-prefix resources/1kg/1KG_EUR_Phase3_chr \
  --snp-identifier rsid \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

Run ordinary h2:

```bash
ldsc h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-file my_gwas.parquet \
  --output-dir results/my_study_h2
```

Run partitioned h2:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-file my_gwas.parquet \
  --output-dir results/my_study_partitioned_h2
```

Query-annotation runs write their per-query result tree automatically. The
retired `--write-per-query-results` flag is rejected; omit it.

Run rg:

```bash
ldsc rg \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-sources trait1.parquet trait2.parquet \
  --output-dir results/trait1_trait2_rg
```

Run all pairwise rg estimates for a trait panel:

```bash
ldsc rg \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-sources results/traits/*.parquet \
  --output-dir results/trait_panel_rg \
  --write-per-pair-detail
```

Run anchor-vs-rest rg estimates:

```bash
ldsc rg \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-sources results/traits/*.parquet \
  --anchor-trait BMI \
  --output-dir results/bmi_anchor_rg
```

## 8. Python API Summary

```python
from ldsc import LDScoreCalculator, LDScoreOutputConfig, load_ldscore_from_dir
from ldsc import load_sumstats
from ldsc import RgDirectoryWriter, RgOutputConfig
from ldsc import RegressionRunner

result = LDScoreCalculator().run(
    annotation_bundle=bundle,
    ref_panel=ref_panel,
    ldscore_config=ldscore_config,
    global_config=global_config,
    output_config=LDScoreOutputConfig(output_dir="results/my_study_ldscore"),
)

loaded = load_ldscore_from_dir("results/my_study_ldscore")
sumstats_table = load_sumstats("results/traits/BMI.parquet")
dataset = RegressionRunner(global_config).build_dataset(sumstats_table, loaded)
sumstats_tables = [
    load_sumstats("results/traits/BMI.parquet"),
    load_sumstats("results/traits/LDL.parquet"),
    load_sumstats("results/traits/CAD.parquet"),
]
rg_result = RegressionRunner(global_config).estimate_rg_pairs(
    sumstats_tables,
    loaded,
    anchor_index=0,
)
RgDirectoryWriter().write(
    rg_result,
    RgOutputConfig(output_dir="results/bmi_anchor_rg", write_per_pair_detail=True),
)
```

Writing LD-score workflows and `load_ldscore_from_dir` return `LDScoreSource`. It retains shared `baseline_table`, counts, overlap statistics, SNP metadata, provenance, and `output_paths`, with no query LD-score tables or chromosome result tables. Explicit `read_queries(names)` calls preserve requested column order and can span saved files. They do not cache values or enforce a read-width limit; the caller owns the RAM for those selections.

`LDScoreResult` remains the materialized object for one execution batch. The small prepared-input Python exception, `LDScoreCalculator.run(..., output_config=None)`, returns one such batch with complete in-memory diagnostics and creates no files. Prepare annotations with `AnnotationBundle.from_frames` and supply resolved restrictions explicitly. Multiple batches require output configuration. See the [zero-write example](../../tutorials/ld-score-calculation.md#small-python-calculations-without-writes).

Generation `query_batch_size` bounds sequential output batches. Direct calculation repeats reference work across them; indexed calculation reuses one chromosome operator per worker through its batches. `threads` is capped at chromosome count. Regression chooses its own batch width independently and retains one complete genome-wide model per query. Sources: `LDScoreCalculator.run`, `run_indexed_ldscore`, and the [memory design](annotation-memory-design.md).

`output_paths` records scientific artifacts such as metadata, shared baseline/overlap files, and query batch files; workflow logs such as `diagnostics/ldscore.log` are audit files and are not included.
`RgOutputConfig(write_per_pair_detail=True)` controls the optional per-pair
detail tree when writing `RgResultFamily` through `RgDirectoryWriter`.
