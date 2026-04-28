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
  baseline.parquet
  query.parquet        # omitted when no query annotations were supplied
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

`baseline.parquet` columns:

```text
CHR, SNP, BP, regr_weight, <baseline LD-score columns...>
```

`query.parquet` columns:

```text
CHR, SNP, BP, <query LD-score columns...>
```

`query.parquet` duplicates the SNP key columns intentionally. Loaders validate
that query rows match baseline rows exactly on `CHR/SNP/BP`.

## 3. Inputs

### Baseline Annotations

Baseline `.annot[.gz]` files are passed to `ldsc ldscore` with
`--baseline-annot-paths`. They define the baseline annotation LD-score columns stored
in `baseline.parquet`.

Example:

```text
resources/baseline_v1.2/baseline.@.annot.gz
```

### Query Annotations

Query annotations are optional and become `query.parquet` columns.

- `--query-annot-bed-paths`: BED intervals projected to the baseline SNP universe in
  memory.
- `--query-annot-paths`: pre-built query `.annot[.gz]` files.

These arguments are mutually exclusive.

### Reference Panel

LD scores can be computed from:

- PLINK input through `--plink-path`
- parquet R2 input through `--r2-paths`

Frequency or metadata sidecars are passed with `--metadata-paths` when needed for MAF
or centiMorgan windows.

### Summary Statistics

Regression commands still read munged summary statistics:

- `ldsc h2 --sumstats-path <file>`
- `ldsc partitioned-h2 --sumstats-path <file>`
- `ldsc rg --sumstats-1-path <file> --sumstats-2-path <file>`

## 4. SNP Universes

The LD-score phase tracks these SNP sets:

| Symbol | Name | Meaning |
|--------|------|---------|
| B | Baseline annotation SNPs | Rows in loaded baseline annotation files |
| A | Raw reference panel SNPs | Rows in PLINK `.bim` or parquet metadata sidecar |
| A' | Restricted reference panel | `A ∩ ref_panel_snps_path`; equals A when absent |
| `ld_reference_snps` | LD computation universe | `B ∩ A'` |
| C | Regression SNP mask | Optional mask from `regression_snps_path` |
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
  "common_reference_snp_count_maf_gt_0_05": 85000.0
}
```

The `common_reference_snp_count_maf_gt_0_05` key is omitted when common-SNP
counts are unavailable. Regression falls back to `all_reference_snp_count` in
that case.

## 6. Regression Behavior

`h2` and `rg` use baseline LD scores only, even when `query.parquet` exists.
They also use the embedded `regr_weight` column from `baseline.parquet`.

`partitioned-h2` reads both tables and loops over query columns. For each query
column, it builds a regression dataset with:

```text
baseline LD-score columns + one query LD-score column
```

This keeps baseline annotations fixed while estimating one query annotation at a
time.

## 7. CLI Examples

Compute LD scores:

```bash
ldsc ldscore \
  --output-dir results/my_study_ldscore \
  --baseline-annot-paths resources/baseline_v1.2/baseline.@.annot.gz \
  --query-annot-bed-paths my_peaks.bed \
  --plink-path resources/1kg/1KG_EUR_Phase3_chr@ \
  --ref-panel-snps-path resources/w_hm3.snplist \
  --regression-snps-path resources/w_hm3.snplist \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

Run ordinary h2:

```bash
ldsc h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-path my_gwas.sumstats.gz \
  --output-dir results/my_study_h2
```

Run partitioned h2:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-path my_gwas.sumstats.gz \
  --output-dir results/my_study_partitioned_h2
```

Run rg:

```bash
ldsc rg \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-1-path trait1.sumstats.gz \
  --sumstats-2-path trait2.sumstats.gz \
  --output-dir results/trait1_trait2_rg
```

## 8. Python API Summary

```python
from ldsc import LDScoreCalculator, LDScoreOutputConfig, load_ldscore_from_dir
from ldsc import RegressionRunner

result = LDScoreCalculator().run(
    annotation_bundle=bundle,
    ref_panel=ref_panel,
    ldscore_config=ldscore_config,
    global_config=global_config,
    output_config=LDScoreOutputConfig(output_dir="results/my_study_ldscore"),
)

loaded = load_ldscore_from_dir("results/my_study_ldscore")
dataset = RegressionRunner(global_config).build_dataset(sumstats_table, loaded)
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
