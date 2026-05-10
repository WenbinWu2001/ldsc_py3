# Cell-Specific LDSC

Goal: estimate cell-specific enrichment by running partitioned LDSC with one query annotation per cell type.

In this package, cell-specific LDSC is the `partitioned-h2` workflow applied to cell-type annotations. Baseline annotations stay in the model as covariates, and each cell-type query column is tested in a baseline-plus-one-query model through `RegressionRunner.estimate_partitioned_h2_batch()`.

Cell-type query annotations require explicit baseline annotations. The
synthetic all-ones `base` annotation is reserved for ordinary unpartitioned
LD-score generation when no query inputs are present. `partitioned-h2` rejects
baseline-only LD-score directories rather than treating `base` as a cell-type
query.

The examples below assume chromosome-pattern inputs such as
`baseline.@.annot.gz`, `cell_type_beds/*.bed`, and a package-built
build-specific R2 directory. Package-built R2 parquet files store
`ldsc:r2_bias` and `ldsc:n_samples` in schema metadata, so the examples omit
R2 bias and sample-size arguments.

Output directories are literal destinations. Missing directories are created,
existing directories are reused, and existing owned workflow artifacts are
refused before writing unless you pass `--overwrite` or `overwrite=True`.
Successful overwrites remove stale owned siblings not produced by the current
configuration and preserve unrelated files. CLI workflow logs are part of that
preflight policy, but they are audit files rather than returned data artifacts.

For `chr_pos` workflows, `genome_build="auto"` can infer hg19/hg38 and
0-based/1-based coordinates during annotation or LD-score loading. The same
logic is public in Python as `infer_chr_pos_build()` and
`resolve_chr_pos_table()` from `ldsc`; the CLI keeps it under workflow flags such
as `--genome-build auto` rather than a separate command.

## Python API

```python
from ldsc import (
    AnnotationBuildConfig,
    AnnotationBuilder,
    GlobalConfig,
    LDScoreCalculator,
    LDScoreConfig,
    LDScoreOutputConfig,
    RefPanelLoader,
    RefPanelConfig,
    RegressionConfig,
    RegressionRunner,
    load_sumstats,
    set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(
    snp_identifier="chr_pos",
    genome_build="hg19",
)
set_global_config(GLOBAL_CONFIG)

annotation_bundle = AnnotationBuilder(GLOBAL_CONFIG, AnnotationBuildConfig()).run(
    AnnotationBuildConfig(
        baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
        query_annot_bed_sources="annotations/cell_type_beds/*.bed",
    )
)

ref_panel = RefPanelLoader(GLOBAL_CONFIG).load(
    RefPanelConfig(
        backend="parquet_r2",
        r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg19",
        chromosomes=tuple(annotation_bundle.chromosomes),
        ref_panel_snps_file="filters/reference_universe.tsv.gz",
    )
)

ldscore_result = LDScoreCalculator().run(
    annotation_bundle=annotation_bundle,
    ref_panel=ref_panel,
    ldscore_config=LDScoreConfig(
        ld_wind_cm=1.0,
        regression_snps_file="filters/hapmap3.tsv.gz",
    ),
    global_config=GLOBAL_CONFIG,
    output_config=LDScoreOutputConfig(
        output_dir="tutorial_outputs/cell_specific_ldscores",
        # overwrite=True,  # also removes stale LD-score siblings not produced by this run
    ),
)

# Current disk-loaded sumstats recover config_snapshot from
# sumstats.metadata.json. Older artifacts without that sidecar warn and use
# config_snapshot=None.
sumstats = load_sumstats("tutorial_outputs/trait/sumstats.parquet", trait_name="trait")

runner = RegressionRunner(
    global_config=GLOBAL_CONFIG,
    regression_config=RegressionConfig(),
)
cell_specific = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
)

cell_specific.to_csv("tutorial_outputs/cell_specific_ldsc.tsv", sep="\t", index=False)
print(annotation_bundle.query_columns)
print(cell_specific)
```

`annotation_bundle.query_columns` are derived from the query annotation files. For BED inputs, the column names come from the BED basenames after normalization.

The LD-score result and annotation bundle retain known `GlobalConfig`
snapshots. Current disk-loaded sumstats recover the same provenance from
`sumstats.metadata.json`. If an older sumstats artifact lacks that sidecar,
regression treats the sumstats provenance as unknown. With `snp_identifier="chr_pos"`,
the merge uses normalized `CHR:POS` coordinates rather than rsIDs, and the
`SNP` column remains a label. If raw sumstats coordinates need a build
conversion before this regression step, run `munge-sumstats` with
`--target-genome-build` plus either `--liftover-chain-file` or
`--use-hm3-quick-liftover`.

## CLI

First compute baseline-plus-cell-type LD scores. This can project BED files directly without materializing intermediate query `.annot.gz` files:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/cell_specific_ldscores \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "annotations/cell_type_beds/*.bed" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg19" \
  --ref-panel-snps-file filters/reference_universe.tsv.gz \
  --regression-snps-file filters/hapmap3.tsv.gz \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

When reusable query `.annot.gz` shards are useful, use `ldsc annotate` or
`run_bed_to_annot(...)`; both are public `ldsc.annotation_builder` workflow
entry points and return the same `AnnotationBundle` shape used above.

Then run partitioned h2 over the cell-type query columns:

```bash
ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.parquet \
  --trait-name trait \
  --ldscore-dir tutorial_outputs/cell_specific_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/cell_specific_ldsc \
  --write-per-query-results
```

The regression reads query annotation columns from
`tutorial_outputs/cell_specific_ldscores/manifest.json` and
`ldscore.query.parquet`. The LD-score parquet files are flat files with
chromosome-aligned row groups, and the manifest lists those row groups for
targeted chromosome reads. The output file is
`tutorial_outputs/cell_specific_ldsc/partitioned_h2.tsv`, with
`partitioned-h2.log` in the same directory. Its key columns are
`Category`, `Prop._SNPs`, `Prop._h2`, `Enrichment`, `Enrichment_p`,
`Coefficient`, and `Coefficient_p`.
For full column definitions, see
[partitioned-h2-results.md](../docs/current/partitioned-h2-results.md).
With `--write-per-query-results`, the command also writes
`tutorial_outputs/cell_specific_ldsc/query_annotations/manifest.tsv` and one
sanitized folder per cell-type query annotation. Each folder contains the
one-row query summary, the baseline-plus-query `partitioned_h2_full.tsv`, and
`metadata.json` with the original annotation name.
If the partitioned summary already exists, `ldsc partitioned-h2` fails before
writing; the same is true for `partitioned-h2.log` and any stale
`query_annotations/` tree. Add `--overwrite` only when replacing it is
intentional. If an overwrite rerun omits `--write-per-query-results`, the old
per-query tree is removed after the new aggregate summary is written.
