# LD Score Calculation

Goal: compute LDSC-compatible LD scores from a reference panel alone, from pre-built SNP-level annotation files, or from raw BED intervals plus an explicit baseline.

The examples below assume chromosome-pattern annotation inputs such as
`annotations/baseline.1.annot.gz` and a package-built R2 directory such as
`r2_ref_panel_1kg30x_1cM_hm3/hg38`.
For the parquet backend, the R2 file should use the canonical six-column schema
(`CHR`, `POS_1`, `POS_2`, `SNP_1`, `SNP_2`, `R2`) written with PyArrow row groups.
The matching `chr*_meta.tsv.gz` sidecars are optional but strongly recommended
because they define the full reference-panel SNP universe and supply `MAF`/`CM`.
Package-built R2 parquet files also store `ldsc:r2_bias` and `ldsc:n_samples`
in Arrow schema metadata, so the examples omit R2 bias and sample-size
arguments. Pass `r2_bias_mode` / `--r2-bias-mode` and
`r2_sample_size` / `--r2-sample-size` manually only for legacy or external raw
R2 parquet files that lack LDSC schema metadata. External raw R2 parquet inputs
are supported only in `rsid` and `chr_pos`; allele-aware modes require
package-built canonical R2 parquet with endpoint allele columns
`A1_1/A2_1/A1_2/A2_2`.

Input-token rules used below:

- exact path: one concrete file
- glob token: one token that expands to multiple files, for example `"beds/*.bed"`
- explicit chromosome suite: `annotations/baseline.@.annot.gz`

Output directories stay literal; only input fields are expanded.

Resolution behavior:

- there is no separate `*_chr` argument anymore; the same public argument now accepts exact paths, globs, or explicit `@` suite tokens
- group inputs such as `--baseline-annot-sources`, `--query-annot-sources`, and `--query-annot-bed-sources` may resolve to many files; package-built parquet panels are supplied as one build directory with `--r2-dir`
- when a group token resolves to chromosome-sharded files, the workflow tries to keep only the files whose names match the active chromosome
- if filename-based chromosome filtering is not possible, the workflow reads the matched files and filters rows by `CHR` internally
- scalar inputs still must resolve to exactly one file

Genome-build behavior for `chr_pos` inputs:

- Python callers can import `infer_chr_pos_build()` and `resolve_chr_pos_table()` from `ldsc`
- workflow calls can set `GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="auto")`
- CLI calls can pass `--snp-identifier chr_pos_allele_aware --genome-build auto`
- auto mode infers hg19 or hg38 and converts 0-based `POS` values to canonical 1-based coordinates when enough reference SNPs are present
- the CLI has no standalone build-inference command; inference is part of existing workflows

Public SNP identifier modes are exactly `rsid`, `rsid_allele_aware`,
`chr_pos`, and `chr_pos_allele_aware`; the default is
`chr_pos_allele_aware`. Mode names are exact. Column aliases apply only to
input headers. Annotation files may omit alleles in allele-aware modes because
they describe genomic membership, but if annotation alleles are present they
participate in allele-aware matching.

Important output behavior:

- the in-memory result is one merged `LDScoreResult` with split `baseline_table` and optional `query_table`
- `--output-dir` writes a canonical LD-score result directory containing `manifest.json`, `ldscore.baseline.parquet`, optional `ldscore.query.parquet`, and `ldscore.log`
- `ldscore.baseline.parquet` and `ldscore.query.parquet` are still single flat files, but each parquet row group contains exactly one chromosome
- `manifest.json` records `row_group_layout`, `baseline_row_groups`, and `query_row_groups` for readers that want to load one chromosome by row-group index
- `LDScoreResult.output_paths` lists scientific data artifacts only; it does not include `ldscore.log`
- regression-universe LD scores live in the `regression_ld_scores` column of `ldscore.baseline.parquet`; there is no separate `.w.l2.ldscore.gz` output
- annotation counts are stored as manifest records, not as separate `.M` files
- if both baseline and query inputs are omitted, the workflow synthesizes an all-ones baseline column named exactly `base` over retained reference-panel metadata
- query `.annot` and BED inputs require explicit baseline annotations; create an explicit all-ones `base` baseline annotation yourself if you intentionally want query annotations tested against that universe
- missing output directories are created and existing directories are reused
- existing owned LD-score artifacts, including unselected siblings such as a
  stale `ldscore.query.parquet`, fail before writing starts; reruns that should
  replace them must pass `--overwrite` or `overwrite=True`
- with overwrite enabled, successful baseline-only runs remove stale
  `ldscore.query.parquet` so the result directory reflects the current
  manifest

Performance behavior for canonical parquet R2 input:

- `snp_batch_size` / `--snp-batch-size` controls the number of SNPs processed per LD-score sliding batch; the default is `128`
- the canonical parquet reader automatically sizes a decoded row-group cache once per chromosome from the actual LD window and `snp_batch_size`
- cache entries are decoded row groups with numeric endpoint indices and R2 values, not pandas DataFrames or SNP strings
- the cache only avoids rereading overlapping parquet row groups; every query still computes its required row-group set and filters to the current SNP window, so cache state does not affect correctness

## Case 1: Ordinary Unpartitioned LD Scores

For unpartitioned heritability, no baseline annotation input is needed. The
result directory still has the normal `ldscore.baseline.parquet` and manifest
layout, with `baseline_columns == ["base"]`.

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="chr_pos_allele_aware",
        genome_build="hg38",
    )
)

result = run_ldscore(
    output_dir="tutorial_outputs/unpartitioned_ldscores",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    use_hm3_ref_panel_snps=True,
    use_hm3_regression_snps=True,
    common_maf_min=0.05,
    ld_wind_cm=1.0,
    # snp_batch_size=128,  # optional; also controls parquet cache sizing
)

print(result.baseline_columns)
print(result.baseline_table.loc[:, ["CHR", "SNP", "POS", "regression_ld_scores", "base"]].head())
```

### CLI

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/unpartitioned_ldscores \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --use-hm3-ref-panel-snps \
  --use-hm3-regression-snps \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
# Add --snp-batch-size 128 only when intentionally tuning the sliding batch size.
```

## Case 2: Existing SNP-Level Annotation Files

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="chr_pos_allele_aware",
        genome_build="hg38",
    )
)

result = run_ldscore(
    output_dir="tutorial_outputs/r2_ldscores",
    baseline_annot_sources="annotations/baseline.@.annot.gz",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    use_hm3_ref_panel_snps=True,
    use_hm3_regression_snps=True,
    common_maf_min=0.05,
    ld_wind_cm=1.0,
    # overwrite=True,  # also removes stale LD-score siblings not produced by this run
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
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --use-hm3-ref-panel-snps \
  --use-hm3-regression-snps \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
# Add --overwrite only when intentionally replacing the LD-score artifact family.
```

## Case 3: Use BED Files Directly During LD-Score Calculation

`ldsc ldscore` can now project BED files in memory through `--query-annot-bed-sources` without first materializing intermediate query `.annot.gz` files. Because BED inputs are query annotations, this mode requires explicit baseline annotations.

### Python API

```python
from ldsc import GlobalConfig, run_ldscore, set_global_config

set_global_config(
    GlobalConfig(
        snp_identifier="chr_pos_allele_aware",
        genome_build="hg38",
    )
)

result = run_ldscore(
    output_dir="tutorial_outputs/r2_ldscores_with_queries",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot_bed_sources="beds/*.bed",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg38",
    use_hm3_ref_panel_snps=True,
    use_hm3_regression_snps=True,
    common_maf_min=0.05,
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
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --use-hm3-ref-panel-snps \
  --use-hm3-regression-snps \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg38 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

## Optional: Read One Chromosome From A Result Directory

Full-file readers such as `pd.read_parquet("ldscore.baseline.parquet")` still work.
For large outputs, use the manifest row-group metadata to read exactly one
chromosome from `ldscore.baseline.parquet` or `ldscore.query.parquet`.

```python
import json
from pathlib import Path

import pyarrow.parquet as pq

ldscore_dir = Path("tutorial_outputs/r2_ldscores")
manifest = json.loads((ldscore_dir / "manifest.json").read_text())

if manifest["row_group_layout"] != "one_per_chromosome":
    raise ValueError(f"Unsupported row-group layout: {manifest['row_group_layout']!r}")

baseline_rg_by_chrom = {
    entry["chrom"]: entry["row_group_index"]
    for entry in manifest["baseline_row_groups"]
}

pf = pq.ParquetFile(ldscore_dir / "ldscore.baseline.parquet")
baseline_chr22 = pf.read_row_group(baseline_rg_by_chrom["22"]).to_pandas()

print(baseline_chr22["CHR"].unique())
print(baseline_chr22.head())
```

If query annotations were supplied, `manifest["query_row_groups"]` has the same
shape for `ldscore.query.parquet`. It is `None` for baseline-only LD-score
results.

## Optional: Materialize BED Projections for Reuse

If you want reusable query `.annot.gz` shards on disk, call `run_bed_to_annot(...)` or `ldsc annotate` explicitly.

### Python API

```python
from ldsc import GlobalConfig, run_bed_to_annot, set_global_config

set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))

bundle = run_bed_to_annot(
    query_annot_bed_sources="beds/*.bed",
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
    # overwrite=True,  # also removes stale query shards outside the current chromosome set
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
automatically. If any root-level `query.*.annot.gz` shard already exists, the
command fails before writing any shard, even if that shard is outside the
current chromosome set. Add `--overwrite` only for an intentional rerun; stale
query shards not produced by the successful run are removed.

For script-style annotation runs, `ldsc.annotation_builder.main(argv)` is the
parser entry point and returns the produced `AnnotationBundle`. If you already
have an argparse namespace from the unified `ldsc` parser, call
`ldsc.annotation_builder.run_annotate_from_args(args)`; it dispatches the
workflow directly without converting the namespace back to argv.

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
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg38" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build auto \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

## Configuration Notes

For Python workflows, `GlobalConfig` now carries only shared runtime settings such as:

- `snp_identifier`
- `genome_build`
- `log_level`

Per-run SNP-universe controls are owned by the workflow-specific configs instead:

- `ref_panel_snps_file` or `use_hm3_ref_panel_snps` belongs to the LD-score reference-panel input and is passed through `run_ldscore(...)` into `RefPanelConfig`
- the LD-score workflow intersects each chromosome bundle with `ref_panel.load_metadata(chrom)`, so reference-panel SNP restriction shrinks the sidecar-defined compute-time universe from `B` to `B ∩ A'`; in the no-annotation unpartitioned case, synthetic `B` is the retained reference-panel metadata itself
- `regression_snps_file` or `use_hm3_regression_snps` belongs to the LD-score calculation config and further restricts the normalized `baseline_table` rows from `B ∩ A'` to `B ∩ A' ∩ C`

`run_bed_to_annot()` no longer applies a reference-panel SNP restriction. It projects BED intervals onto the baseline annotation rows and returns an `AnnotationBundle`; any reference-panel restriction is applied later, during LD-score calculation, when the workflow aligns each chromosome bundle to the prepared reference-panel metadata.

In-process LD-score results carry a frozen `config_snapshot`. If you later change the registered `GlobalConfig`, existing results keep their original snapshot, and downstream merge points raise `ConfigMismatchError` if you try to combine artifacts produced under incompatible identifier or genome-build assumptions. Legacy LD-score directories whose manifest is missing usable config provenance still load with a warning and `config_snapshot=None`.
