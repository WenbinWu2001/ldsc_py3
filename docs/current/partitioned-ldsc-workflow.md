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
CHR, SNP, POS, regr_weight, <baseline LD-score columns...>
```

`query.parquet` columns:

```text
CHR, SNP, POS, <query LD-score columns...>
```

`query.parquet` duplicates the SNP key columns intentionally. Loaders validate
that query rows match baseline rows exactly on `CHR/SNP/POS`.

## 3. Inputs

### Baseline Annotations

Baseline `.annot[.gz]` files are passed to `ldsc ldscore` with
`--baseline-annot-sources`. They define the baseline annotation LD-score columns stored
in `baseline.parquet`.

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

Query annotations are optional and become `query.parquet` columns.

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
- package-built parquet R2 input through `--ref-panel-dir`
- advanced/compatibility parquet R2 input through `--r2-sources`

In `--ref-panel-dir` mode, matching `chr*_meta.tsv.gz` sidecars are discovered
automatically when present. Advanced frequency or metadata sidecars can still be
passed with `--metadata-sources` when needed for MAF or centiMorgan windows.

### Summary Statistics

Regression commands still read munged summary statistics:

- `ldsc h2 --sumstats-file <file>`
- `ldsc partitioned-h2 --sumstats-file <file>`
- `ldsc rg --sumstats-1-file <file> --sumstats-2-file <file>`

Current `ldsc munge-sumstats` outputs include canonical `CHR` and `POS` columns
beside `SNP`, `Z`, and `N`, and write `sumstats.metadata.json` beside the table.
The sidecar stores the effective `snp_identifier` and `genome_build`, including
auto-inferred build details for `chr_pos` runs. Older `.sumstats.gz` files
without the sidecar still load, but their config provenance is treated as
unknown.

## 4. SNP Universes

The LD-score phase tracks these SNP sets:

| Symbol | Name | Meaning |
|--------|------|---------|
| B | Baseline annotation SNPs | Rows in loaded baseline annotation files, or retained reference-panel metadata rows for synthetic `base` unpartitioned runs |
| A | Raw reference panel SNPs | Rows in PLINK `.bim` or parquet metadata sidecar |
| A' | Prepared reference panel | `A ∩ ref_panel_snps_file`, after optional retained-panel `maf_min` and PLINK `keep_indivs_file`; equals A when absent |
| `ld_reference_snps` | LD computation universe | `B ∩ A'` |
| C | Regression SNP mask | Optional mask from `regression_snps_file` |
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

`h2` and `rg` use baseline LD scores only, even when `query.parquet` exists.
They also use the embedded `regr_weight` column from `baseline.parquet`.
When the effective identifier mode is `rsid`, regression merges on `SNP`. When
it is `chr_pos`, regression builds a normalized private `CHR:POS` key from both
the sumstats and LD-score tables and merges on that coordinate key.

`partitioned-h2` reads both tables and loops over query columns. For each query
column, it builds a regression dataset with:

```text
baseline LD-score columns + one query LD-score column
```

This keeps baseline annotations fixed while estimating one query annotation at a
time. By default the CLI writes the aggregate `partitioned_h2.tsv` table only.
With `--write-per-query-results`, it also writes a staged
`query_annotations/` tree containing a manifest and one sanitized folder per
query annotation.

## 7. CLI Examples

Compute LD scores:

```bash
ldsc ldscore \
  --output-dir results/my_study_ldscore \
  --baseline-annot-sources resources/baseline_v1.2/baseline.@.annot.gz \
  --query-annot-bed-sources my_peaks.bed \
  --plink-prefix resources/1kg/1KG_EUR_Phase3_chr@ \
  --ref-panel-snps-file resources/w_hm3.snplist \
  --regression-snps-file resources/w_hm3.snplist \
  --snp-identifier rsid \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

Compute ordinary unpartitioned LD scores without baseline annotations:

```bash
ldsc ldscore \
  --output-dir results/my_unpartitioned_ldscore \
  --plink-prefix resources/1kg/1KG_EUR_Phase3_chr@ \
  --ref-panel-snps-file resources/w_hm3.snplist \
  --regression-snps-file resources/w_hm3.snplist \
  --snp-identifier rsid \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

Run ordinary h2:

```bash
ldsc h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-file my_gwas.sumstats.gz \
  --output-dir results/my_study_h2
```

Run partitioned h2:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-file my_gwas.sumstats.gz \
  --output-dir results/my_study_partitioned_h2
```

Run partitioned h2 and save per-query folders:

```bash
ldsc partitioned-h2 \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-file my_gwas.sumstats.gz \
  --output-dir results/my_study_partitioned_h2 \
  --write-per-query-results
```

Run rg:

```bash
ldsc rg \
  --ldscore-dir results/my_study_ldscore \
  --sumstats-1-file trait1.sumstats.gz \
  --sumstats-2-file trait2.sumstats.gz \
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
