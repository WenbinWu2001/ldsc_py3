# LD Score Calculation

Goal: compute LDSC-compatible LD scores from either pre-built SNP-level annotation files or raw BED intervals.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
For the parquet backend, the R2 file should use the canonical six-column schema
(`CHR`, `POS_1`, `POS_2`, `R2`, `SNP_1`, `SNP_2`) written with PyArrow row groups.
The paired metadata sidecar is required and defines the raw reference-panel SNP
universe used during LD-score computation.

Input-token rules used below:

- exact path: one concrete file
- glob token: one token that expands to multiple files, for example `"beds/*.bed"`
- explicit chromosome suite: `annotations/baseline.@.annot.gz`

Output prefixes stay literal; only input fields are expanded.

Resolution behavior:

- there is no separate `*_chr` argument anymore; the same public argument now accepts exact paths, globs, or explicit `@` suite tokens
- group inputs such as `--baseline-annot`, `--query-annot`, `--query-annot-bed`, `--r2-table`, and `--frqfile` may resolve to many files
- when a group token resolves to chromosome-sharded files, the workflow tries to keep only the files whose names match the active chromosome
- if filename-based chromosome filtering is not possible, the workflow reads the matched files and filters rows by `CHR` internally
- scalar inputs still must resolve to exactly one file

Important output behavior:

- the in-memory result is always one merged `LDScoreResult` with a single `ldscore_table`
- the default output layer writes one `.l2.ldscore.gz` file per chromosome, exposed as `output_paths["ldscore.chrom_<chrom>"]`
- regression weights now live in the `regr_weight` column of each LD-score table; there is no separate `.w.l2.ldscore.gz` output
- count files such as `.l2.M` and `.l2.M_5_50` are aggregate run-level artifacts

## Case 1: Existing SNP-Level Annotation Files

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="rsid",
        genome_build="hg38",
    )
)

result = run_ldscore(
    out="tutorial_outputs/r2_ldscores",
    baseline_annot="annotations/baseline.@.annot.gz",
    r2_table="r2/reference.@.parquet",
    frqfile="r2/reference_metadata.@.tsv.gz",
    ref_panel_snps_path="filters/reference_universe.txt",
    regression_snps_path="filters/hapmap3.txt",
    ld_wind_cm=1.0,
)

first_ldscore_path = next(
    path
    for name, path in sorted(result.output_paths.items())
    if name.startswith("ldscore.chrom_")
)

print(result.ldscore_table.columns.tolist())
print(result.ldscore_table.head())
print(first_ldscore_path)
print(result.output_paths["counts.common_reference_snp_counts_maf_gt_0_05"])
print(result.config_snapshot)
```

### CLI

```bash
ldsc ldscore \
  --out tutorial_outputs/r2_ldscores \
  --baseline-annot "annotations/baseline.@.annot.gz" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --ref-panel-snps-path filters/reference_universe.txt \
  --regression-snps-path filters/hapmap3.txt \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

## Case 2: Use BED Files Directly During LD-Score Calculation

`ldsc ldscore` can now project BED files in memory through `--query-annot-bed` without first materializing intermediate query `.annot.gz` files.

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="rsid",
        genome_build="hg38",
    )
)

result = run_ldscore(
    out="tutorial_outputs/r2_ldscores_with_queries",
    baseline_annot="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot_bed="beds/*.bed",
    r2_table="r2/reference.@.parquet",
    frqfile="r2/reference_metadata.@.tsv.gz",
    ref_panel_snps_path="filters/reference_universe.txt",
    regression_snps_path="filters/hapmap3.txt",
    ld_wind_cm=1.0,
)

print(result.query_columns)
print(result.ldscore_table.head())
```

### CLI

```bash
ldsc ldscore \
  --out tutorial_outputs/r2_ldscores_with_queries \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed "beds/*.bed" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --ref-panel-snps-path filters/reference_universe.txt \
  --regression-snps-path filters/hapmap3.txt \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

## Optional: Materialize BED Projections for Reuse

If you want reusable query `.annot.gz` shards on disk, call `run_bed_to_annot(...)` or `ldsc annotate` explicitly.

### Python API

```python
from ldsc import GlobalConfig, run_bed_to_annot, set_global_config

set_global_config(GlobalConfig(snp_identifier="rsid", genome_build="hg38"))

bundle = run_bed_to_annot(
    bed_files="beds/*.bed",
    baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
)

print(bundle.query_columns)
print(bundle.chromosomes)
```

### CLI

```bash
ldsc annotate \
  --bed-files "beds/*.bed" \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds
```

The generated query shards are named `query.<chrom>.annot.gz`, so downstream inputs should use a token such as `annotations/query_from_beds/query.@.annot.gz`.

If `--output-dir` does not exist yet, the workflow warns once and creates it automatically.

## Configuration Notes

For Python workflows, `GlobalConfig` now carries only shared runtime settings such as:

- `snp_identifier`
- `genome_build`
- `log_level`

Per-run SNP-universe controls are owned by the workflow-specific configs instead:

- `ref_panel_snps_path` belongs to the LD-score reference-panel input and is passed through `run_ldscore(...)` into `RefPanelSpec`
- the LD-score workflow intersects each chromosome bundle with `ref_panel.load_metadata(chrom)`, so `ref_panel_snps_path` shrinks the sidecar-defined compute-time universe from `B` to `B ∩ A'`
- `regression_snps_path` belongs to the LD-score calculation config and further restricts the normalized `ldscore_table` rows from `B ∩ A'` to `B ∩ A' ∩ C`

`run_bed_to_annot()` no longer applies a reference-panel SNP restriction. It projects BED intervals onto the baseline annotation rows and returns an `AnnotationBundle`; any reference-panel restriction is applied later, during LD-score calculation, when the workflow aligns each chromosome bundle to the prepared reference-panel metadata.

Each returned result carries a frozen `config_snapshot`. If you later change the registered `GlobalConfig`, existing results keep their original snapshot, and downstream merge points raise `ConfigMismatchError` if you try to combine artifacts produced under incompatible identifier or genome-build assumptions.
