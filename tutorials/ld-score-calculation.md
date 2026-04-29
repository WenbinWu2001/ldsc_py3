# LD Score Calculation

Goal: compute LDSC-compatible LD scores from a reference panel alone, from pre-built SNP-level annotation files, or from raw BED intervals plus an explicit baseline.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
For the parquet backend, the R2 file should use the canonical six-column schema
(`CHR`, `POS_1`, `POS_2`, `R2`, `SNP_1`, `SNP_2`) written with PyArrow row groups.
The paired metadata sidecar is required and defines the raw reference-panel SNP
universe used during LD-score computation.

Input-token rules used below:

- exact path: one concrete file
- glob token: one token that expands to multiple files, for example `"beds/*.bed"`
- explicit chromosome suite: `annotations/baseline.@.annot.gz`

Output directories stay literal; only input fields are expanded.

Resolution behavior:

- there is no separate `*_chr` argument anymore; the same public argument now accepts exact paths, globs, or explicit `@` suite tokens
- group inputs such as `--baseline-annot-sources`, `--query-annot-sources`, `--query-annot-bed-sources`, `--r2-sources`, and `--metadata-sources` may resolve to many files
- when a group token resolves to chromosome-sharded files, the workflow tries to keep only the files whose names match the active chromosome
- if filename-based chromosome filtering is not possible, the workflow reads the matched files and filters rows by `CHR` internally
- scalar inputs still must resolve to exactly one file

Genome-build behavior for `chr_pos` inputs:

- Python callers can import `infer_chr_pos_build()` and `resolve_chr_pos_table()` from `ldsc`
- workflow calls can set `GlobalConfig(snp_identifier="chr_pos", genome_build="auto")`
- CLI calls can pass `--snp-identifier chr_pos --genome-build auto`
- auto mode infers hg19 or hg38 and converts 0-based `POS` values to canonical 1-based coordinates when enough reference SNPs are present
- the CLI has no standalone build-inference command; inference is part of existing workflows

Important output behavior:

- the in-memory result is one merged `LDScoreResult` with split `baseline_table` and optional `query_table`
- `--output-dir` writes a canonical LD-score result directory containing `manifest.json`, `baseline.parquet`, and optional `query.parquet`
- regression weights live in the `regr_weight` column of `baseline.parquet`; there is no separate `.w.l2.ldscore.gz` output
- annotation counts are stored as manifest records, not as separate `.M` files
- if both baseline and query inputs are omitted, the workflow synthesizes an all-ones baseline column named exactly `base` over retained reference-panel metadata
- query `.annot` and BED inputs require explicit baseline annotations; create an explicit all-ones `base` baseline annotation yourself if you intentionally want query annotations tested against that universe
- missing output directories are created and existing directories are reused
- existing fixed output files fail before writing starts; reruns that should
  replace them must pass `--overwrite` or `overwrite=True`

## Case 1: Ordinary Unpartitioned LD Scores

For unpartitioned heritability, no baseline annotation input is needed. The
result directory still has the normal `baseline.parquet` and manifest layout,
with `baseline_columns == ["base"]`.

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
    output_dir="tutorial_outputs/unpartitioned_ldscores",
    r2_sources="r2/reference.@.parquet",
    metadata_sources="r2/reference_metadata.@.tsv.gz",
    r2_bias_mode="unbiased",
    ref_panel_snps_file="filters/reference_universe.txt",
    regression_snps_file="filters/hapmap3.txt",
    ld_wind_cm=1.0,
)

print(result.baseline_columns)
print(result.baseline_table.loc[:, ["CHR", "SNP", "BP", "regr_weight", "base"]].head())
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/unpartitioned_ldscores \
  --r2-sources "r2/reference.@.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
  --ref-panel-snps-file filters/reference_universe.txt \
  --regression-snps-file filters/hapmap3.txt \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
```

## Case 2: Existing SNP-Level Annotation Files

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
    output_dir="tutorial_outputs/r2_ldscores",
    baseline_annot_sources="annotations/baseline.@.annot.gz",
    r2_sources="r2/reference.@.parquet",
    metadata_sources="r2/reference_metadata.@.tsv.gz",
    r2_bias_mode="unbiased",
    ref_panel_snps_file="filters/reference_universe.txt",
    regression_snps_file="filters/hapmap3.txt",
    ld_wind_cm=1.0,
    # overwrite=True,  # enable only when intentionally replacing prior outputs
)

print(result.baseline_table.columns.tolist())
print(result.baseline_table.head())
print(result.output_paths["baseline"])
print(result.output_paths["manifest"])
print(result.config_snapshot)
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/r2_ldscores \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --r2-sources "r2/reference.@.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
  --ref-panel-snps-file filters/reference_universe.txt \
  --regression-snps-file filters/hapmap3.txt \
  --snp-identifier rsid \
  --ld-wind-cm 1.0
# Add --overwrite only when intentionally replacing manifest/baseline/query files.
```

## Case 3: Use BED Files Directly During LD-Score Calculation

`ldsc ldscore` can now project BED files in memory through `--query-annot-bed-sources` without first materializing intermediate query `.annot.gz` files. Because BED inputs are query annotations, this mode requires explicit baseline annotations.

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
    output_dir="tutorial_outputs/r2_ldscores_with_queries",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot_bed_sources="beds/*.bed",
    r2_sources="r2/reference.@.parquet",
    metadata_sources="r2/reference_metadata.@.tsv.gz",
    r2_bias_mode="unbiased",
    ref_panel_snps_file="filters/reference_universe.txt",
    regression_snps_file="filters/hapmap3.txt",
    ld_wind_cm=1.0,
)

print(result.query_columns)
print(result.baseline_table.head())
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/r2_ldscores_with_queries \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --r2-sources "r2/reference.@.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
  --ref-panel-snps-file filters/reference_universe.txt \
  --regression-snps-file filters/hapmap3.txt \
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
    query_annot_bed_sources="beds/*.bed",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
    # overwrite=True,  # enable only when intentionally replacing query shards
)

print(bundle.query_columns)
print(bundle.chromosomes)
```

### CLI

```bash
ldsc annotate \
  --query-annot-bed-sources "beds/*.bed" \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds
```

The generated query shards are named `query.<chrom>.annot.gz`, so downstream inputs should use a token such as `annotations/query_from_beds/query.@.annot.gz`.

If `--output-dir` does not exist yet, the workflow warns once and creates it
automatically. If a `query.<chrom>.annot.gz` shard already exists, the command
fails before writing any shard; add `--overwrite` only for an intentional rerun.

## Optional: Inspect `chr_pos` Genome Build Before A Workflow

Use the top-level Python API when you want to preflight a `CHR`/`POS` table in a
notebook or script. `infer_chr_pos_build()` reports the decision only;
`resolve_chr_pos_table()` also returns a normalized table with canonical
chromosome labels and 1-based positions.

```python
import pandas as pd
from ldsc import infer_chr_pos_build, resolve_chr_pos_table

# The table must contain enough HapMap3-overlapping SNPs to support inference.
restriction = pd.read_csv("filters/hapmap3_chr_pos.tsv.gz", sep="\t")

normalized, inference = resolve_chr_pos_table(
    restriction,
    context="tutorial restriction table",
)
decision_only = infer_chr_pos_build(
    normalized.loc[:, ["CHR", "POS"]],
    context="normalized tutorial restriction table",
)

print(inference.genome_build)
print(inference.coordinate_basis)
print(decision_only.summary_message)
print(normalized.head())
```

For command-line runs, use auto mode on the workflow itself:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/auto_build_ldscores \
  --baseline-annot-sources "annotations/baseline.@.annot.gz" \
  --r2-sources "r2/reference.@.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
  --snp-identifier chr_pos \
  --genome-build auto \
  --ld-wind-cm 1.0
```

## Configuration Notes

For Python workflows, `GlobalConfig` now carries only shared runtime settings such as:

- `snp_identifier`
- `genome_build`
- `log_level`

Per-run SNP-universe controls are owned by the workflow-specific configs instead:

- `ref_panel_snps_file` belongs to the LD-score reference-panel input and is passed through `run_ldscore(...)` into `RefPanelConfig`
- the LD-score workflow intersects each chromosome bundle with `ref_panel.load_metadata(chrom)`, so `ref_panel_snps_file` shrinks the sidecar-defined compute-time universe from `B` to `B ∩ A'`; in the no-annotation unpartitioned case, synthetic `B` is the retained reference-panel metadata itself
- `regression_snps_file` belongs to the LD-score calculation config and further restricts the normalized `baseline_table` rows from `B ∩ A'` to `B ∩ A' ∩ C`

`run_bed_to_annot()` no longer applies a reference-panel SNP restriction. It projects BED intervals onto the baseline annotation rows and returns an `AnnotationBundle`; any reference-panel restriction is applied later, during LD-score calculation, when the workflow aligns each chromosome bundle to the prepared reference-panel metadata.

In-process LD-score results carry a frozen `config_snapshot`. If you later change the registered `GlobalConfig`, existing results keep their original snapshot, and downstream merge points raise `ConfigMismatchError` if you try to combine artifacts produced under incompatible identifier or genome-build assumptions. Legacy LD-score directories whose manifest is missing usable config provenance still load with a warning and `config_snapshot=None`.
