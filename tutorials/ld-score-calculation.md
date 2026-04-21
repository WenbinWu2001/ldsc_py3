# LD Score Calculation

Goal: compute LDSC-compatible reference LD scores and regression weights with the refactored package, first from existing annotation files and then from multiple BED annotations.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.

Input-token rules used below:

- exact path: one concrete file
- glob token: one token that expands to multiple files, for example `"beds/*.bed"`
- explicit chromosome suite: `annotations/baseline.@.annot.gz`

Output prefixes stay literal; only input fields are expanded.

Resolution behavior:

- there is no separate `*_chr` argument anymore; the same public argument now accepts either shared files, globs, or explicit `@` suite tokens
- group inputs such as `--baseline-annot`, `--query-annot`, `--r2-table`, and `--frqfile` may resolve to many files
- when a group token resolves to many chromosome-sharded files, the workflow tries to keep only the files whose names match the current chromosome
- if filename-based chromosome filtering is not possible, the workflow reads the matched files and filters rows by the `CHR` column internally
- scalar inputs still must resolve to exactly one file

## Case 1: Existing SNP-Level Annotation Files

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(GlobalConfig(snp_identifier="rsid"))

result = run_ldscore(
    out="tutorial_outputs/r2_ldscores",
    baseline_annot="annotations/baseline.@.annot.gz",
    r2_table="r2/reference.@.parquet",
    frqfile="r2/reference_metadata.@.tsv.gz",
    ld_wind_cm=1.0,
)

print(result.output_paths["ldscore"])
print(result.output_paths["w_ld"])
print(result.baseline_columns)
```

### CLI

```bash
ldsc ldscore \
  --out tutorial_outputs/r2_ldscores \
  --baseline-annot "annotations/baseline.@.annot.gz" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

## Case 2: Build Query Annotations From BED Files

`ldsc annotate` projects each BED basename into a query annotation column and writes chromosome-matched `.annot.gz` files.
The generated query shards are named `query.<chrom>.annot.gz`, so downstream
inputs should use tokens such as `annotations/query_from_beds/query.@.annot.gz`.

### Python API

```python
from ldsc import GlobalConfig, run_bed_to_annot, run_ldscore, set_global_config

set_global_config(GlobalConfig(snp_identifier="rsid"))

run_bed_to_annot(
    bed_files="beds/*.bed",
    baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
)

result = run_ldscore(
    out="tutorial_outputs/r2_ldscores_with_queries",
    baseline_annot="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot="annotations/query_from_beds/query.@.annot.gz",
    r2_table="r2/reference.@.parquet",
    frqfile="r2/reference_metadata.@.tsv.gz",
    ld_wind_cm=1.0,
)

print(result.output_paths["ldscore"])
print(result.query_columns)
```

### CLI

```bash
ldsc annotate \
  --bed-files "beds/*.bed" \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds

ldsc ldscore \
  --out tutorial_outputs/r2_ldscores_with_queries \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot "annotations/query_from_beds/query.@.annot.gz" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

If `--output-dir` does not exist yet, the workflow warns once and creates it automatically.

For Python workflows, `run_bed_to_annot()` and `run_ldscore()` reuse the
registered `GlobalConfig`; shared settings such as `snp_identifier`,
`genome_build`, `log_level`, and restriction paths are no longer passed
directly into `run_ldscore(...)`.
