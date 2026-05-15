# Partitioned LDSC Workflow: Technical Reference

This document describes the refactored workflow for computing LD scores and
running h2, partitioned-h2, and rg regression from one canonical LD-score result
directory.

## 1. Overview

The pipeline has two phases:

1. `ldsc ldscore --output-dir <ldscore_dir>` computes annotation LD scores and
   writes one self-contained directory.
2. Regression commands read that directory through `--ldscore-dir <ldscore_dir>`.

The public regression CLI no longer accepts fragmented LD-score artifacts:
`--ldscore`, `--counts`, `--w-ld`, `--annotation-manifest`, and
`--query-columns` are removed.

## 2. LD-Score Result Directory

An LD-score run writes:

```text
<ldscore_dir>/
  manifest.json
  ldscore.baseline.parquet
  ldscore.query.parquet        # omitted when no query annotations were supplied
```

`manifest.json` is the only metadata entry point. It contains:

- `format: "ldsc.ldscore_result.v1"`
- relative file paths for `baseline` and optional `query`
- `snp_identifier`, `genome_build`, and processed `chromosomes`
- ordered `baseline_columns` and `query_columns`
- one count record per LD-score annotation column
- `config_snapshot`
- basic row counts

All paths inside the manifest are relative to `ldscore_dir`.

`ldscore.baseline.parquet` columns:

```text
CHR, SNP, POS, regression_ld_scores, <baseline LD-score columns...>
```

`ldscore.query.parquet` columns:

```text
CHR, SNP, POS, <query LD-score columns...>
```

`ldscore.query.parquet` duplicates the SNP key columns intentionally. Loaders validate
that query rows match baseline rows exactly on `CHR/SNP/POS`.

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

Query annotations are optional and become `ldscore.query.parquet` columns.

- `--query-annot-bed-sources`: BED intervals projected to the baseline SNP universe in
  memory.
- `--query-annot-sources`: pre-built query `.annot[.gz]` files.

These arguments are mutually exclusive and require explicit
`--baseline-annot-sources`. If users intentionally want to test query annotations
against an all-ones universe, they should create an explicit all-ones `base`
baseline annotation over the query annotation universe and run the partitioned
workflow with both baseline and query inputs.

### Reference Panel

LD scores can be computed from:

- PLINK input through `--plink-prefix`
- package-built parquet R2 input through `--r2-dir`

In `--r2-dir` mode, matching `chr*_meta.tsv.gz` sidecars are discovered
automatically when present. Sidecars are optional but strongly recommended:
without them, the loader can only synthesize SNP presence from parquet endpoint
columns and cannot recover MAF or complete cM metadata.
External raw R2 parquet inputs are supported only in `rsid` and `chr_pos`.
Allele-aware modes require package-built canonical R2 parquet with endpoint
allele columns `A1_1/A2_1/A1_2/A2_2`.

### Summary Statistics

Regression commands still read munged summary statistics:

- `ldsc h2 --sumstats-file <file>`
- `ldsc partitioned-h2 --sumstats-file <file>`
- `ldsc rg --sumstats-sources <file-or-glob> <file-or-glob> [...]`

Current `ldsc munge-sumstats` outputs include canonical `CHR` and `POS` columns
beside `SNP`, `Z`, and `N`, write `sumstats.parquet` by default, write
`sumstats.metadata.json` beside the table, and write
`dropped_snps/dropped.tsv.gz` for row-level liftover-drop auditing. The metadata
sidecar stores the thin compatibility payload: `schema_version`,
`artifact_type`, `snp_identifier`, `genome_build`, and optional `trait_name`.
Legacy package-written `.sumstats.gz`
files without the current sidecar are rejected and must be regenerated.
In allele-aware modes, current sumstats artifacts require usable `A1/A2`. To
run without allele-aware SNP identity, set `--snp-identifier chr_pos` or
`--snp-identifier rsid` intentionally.

## 4. SNP Universes

The LD-score phase tracks these SNP sets:

| Symbol | Name | Meaning |
|--------|------|---------|
| B | Baseline annotation SNPs | Rows in loaded baseline annotation files, or retained reference-panel metadata rows for synthetic `base` unpartitioned runs |
| A | Raw reference panel SNPs | Rows in PLINK `.bim` or parquet metadata sidecar |
| A' | Prepared reference panel | `A` after optional explicit or HM3 reference-panel SNP restriction, retained-panel `maf_min`, and PLINK `keep_indivs_file`; equals A when absent |
| `ld_reference_snps` | LD computation universe | `B ∩ A'` |
| C | Regression SNP mask | Optional mask from explicit `regression_snps_file` or packaged HM3 regression SNP restriction |
| `ld_regression_snps` | Persisted row set | `B ∩ A' ∩ C`; equals `ld_reference_snps` when C is absent |

LD-score column counts in the manifest are computed over `ld_reference_snps`.
Persisted parquet rows are `ld_regression_snps`.

## 5. Count Records

The manifest stores counts as records keyed by column name, not by array
position:

```json
{
  "group": "baseline",
  "column": "base",
  "all_reference_snp_count": 100000.0,
  "common_reference_snp_count": 85000.0
}
```

The top-level manifest also records the threshold used to compute common counts:

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

`h2` and `rg` use baseline LD scores only, even when `ldscore.query.parquet` exists.
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

`partitioned-h2` requires `ldscore.query.parquet` and a non-empty `query_columns` list in
the LD-score manifest. Baseline-only LD-score directories are valid inputs for
`h2` and `rg`, but `partitioned-h2` rejects them instead of treating the
baseline as a query annotation. For each query column, it builds a regression
dataset with:

```text
baseline LD-score columns + one query LD-score column
```

This keeps baseline annotations fixed while estimating one query annotation at a
time. By default the CLI writes the compact aggregate `partitioned_h2.tsv` table
plus `partitioned-h2.log`. With `--write-per-query-results`, it also writes a
staged `query_annotations/` tree containing `manifest.tsv` and one sanitized
folder per query annotation. Each query folder contains a one-row
`partitioned_h2.tsv`, the full fitted-model `partitioned_h2_full.tsv`, and
`metadata.json`. The log is preflighted with the table outputs, but it is not
part of returned `output_paths`.

`partitioned-h2` treats `partitioned_h2.tsv`, `query_annotations/`, and
`partitioned-h2.log` as one owned output family. Without `--overwrite`, any
existing owned sibling rejects the run. With `--overwrite`, a successful
aggregate-only run removes a stale `query_annotations/` tree from a previous
per-query configuration.

## 7. CLI Examples

Compute LD scores:

```bash
ldsc ldscore \
  --output-dir results/my_study_ldscore \
  --baseline-annot-sources resources/baseline_v1.2/baseline.@.annot.gz \
  --query-annot-bed-sources my_peaks.bed \
  --plink-prefix resources/1kg/1KG_EUR_Phase3_chr@ \
  --use-hm3-ref-panel-snps \
  --use-hm3-regression-snps \
  --snp-identifier rsid \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

Compute ordinary unpartitioned LD scores without baseline annotations:

```bash
ldsc ldscore \
  --output-dir results/my_unpartitioned_ldscore \
  --plink-prefix resources/1kg/1KG_EUR_Phase3_chr@ \
  --use-hm3-ref-panel-snps \
  --use-hm3-regression-snps \
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

Run partitioned h2 and save per-query folders:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-file my_gwas.parquet \
  --output-dir results/my_study_partitioned_h2 \
  --write-per-query-results
```

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

The public `LDScoreResult` shape is split:

- `baseline_table`
- `query_table`
- `count_records`
- `baseline_columns`
- `query_columns`
- `ld_regression_snps`
- `chromosome_results`
- `config_snapshot`
- `output_paths`

`output_paths` records scientific artifacts such as `manifest.json`,
`ldscore.baseline.parquet`, and optional `ldscore.query.parquet`; workflow logs such as
`ldscore.log` are audit files and are not included.
`RgOutputConfig(write_per_pair_detail=True)` controls the optional per-pair
detail tree when writing `RgResultFamily` through `RgDirectoryWriter`.
