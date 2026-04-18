# LD Score Calculation

Goal: compute LDSC-compatible reference LD scores and regression weights with the refactored package, first from existing annotation files and then from multiple BED annotations.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.

Input-token rules used below:

- exact path: one concrete file
- glob token: one token that expands to multiple files, for example `"beds/*.bed"`
- explicit chromosome suite: `baseline.@`
- legacy bare prefix: still accepted, but `@` is preferred

Output prefixes stay literal; only input fields are expanded.

## Case 1: Existing SNP-Level Annotation Files

### Python API

```python
from ldsc import run_ldscore

result = run_ldscore(
    out="tutorial_outputs/r2_ldscores",
    baseline_annot="annotations/baseline.@",
    r2_table="r2/reference.@",
    frqfile="r2/reference_metadata.@",
    snp_identifier="rsid",
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
  --baseline-annot "annotations/baseline.@" \
  --r2-table "r2/reference.@" \
  --frqfile "r2/reference_metadata.@" \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

## Case 2: Build Query Annotations From BED Files

`ldsc annotate` projects each BED basename into a query annotation column and writes chromosome-matched `.annot.gz` files.
The generated query shards are named `query.<chrom>.annot.gz`, so downstream
inputs should use tokens such as `annotations/query_from_beds/query.@`.

### Python API

```python
from ldsc import run_bed_to_annot, run_ldscore

run_bed_to_annot(
    bed_files="beds/*.bed",
    baseline_annot_dir="annotations/baseline_chr",
    output_dir="annotations/query_from_beds",
)

result = run_ldscore(
    out="tutorial_outputs/r2_ldscores_with_queries",
    baseline_annot="annotations/baseline_chr/baseline.@",
    query_annot="annotations/query_from_beds/query.@",
    r2_table="r2/reference.@",
    frqfile="r2/reference_metadata.@",
    snp_identifier="rsid",
    ld_wind_cm=1.0,
)

print(result.output_paths["ldscore"])
print(result.query_columns)
```

### CLI

```bash
ldsc annotate \
  --bed-files "beds/*.bed" \
  --baseline-annot-dir annotations/baseline_chr \
  --output-dir annotations/query_from_beds

ldsc ldscore \
  --out tutorial_outputs/r2_ldscores_with_queries \
  --baseline-annot "annotations/baseline_chr/baseline.@" \
  --query-annot "annotations/query_from_beds/query.@" \
  --r2-table "r2/reference.@" \
  --frqfile "r2/reference_metadata.@" \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

The `--*-chr` flags are still accepted for CLI compatibility, but the generic
input flags can now take chromosome-suite tokens directly.
