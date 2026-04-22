# Partitioned LDSC

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
The bundled `baseline_v1.2` annotations are hg19-based, so the parquet example uses `genome_build="hg19"` to align raw parquet coordinates to the annotation bundle.
The workflow also accepts `hg37` and `GRCh37` as aliases for `hg19`, and `GRCh38` as an alias for `hg38`; outputs always normalize back to canonical `hg19` or `hg38`.

Path-token rules used in this tutorial:

- use `@` for chromosome suites such as `baseline.@.annot.gz`
- use globs when the filenames do not follow the simple chromosome-suffix convention
- scalar inputs still resolve to exactly one file
- output prefixes remain literal destinations

Resolution behavior:

- there is no separate `*_chr` public argument anymore; one argument now handles both shared inputs and chromosome-sharded inputs
- group inputs may expand to many files through a glob or an `@` suite token
- if those files encode chromosome labels in their filenames, the workflow selects the subset for the active chromosome before loading them
- otherwise the workflow loads the matched files and keeps only rows whose `CHR` value matches the active chromosome
- if multiple files contribute annotation columns for the same chromosome, their SNP rows must align exactly and their annotation column names must be unique

## Python API

The Python API is the cleanest end-to-end path because it keeps the merged `AnnotationBundle` and `LDScoreResult` in memory across the whole workflow.

```python
from ldsc import (
    AnnotationBuildConfig,
    AnnotationBuilder,
    AnnotationSourceSpec,
    GlobalConfig,
    LDScoreCalculator,
    LDScoreConfig,
    MungeConfig,
    RawSumstatsSpec,
    RefPanelLoader,
    RefPanelSpec,
    RegressionConfig,
    RegressionRunner,
    SumstatsMunger,
    load_sumstats,
    run_bed_to_annot,
    set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(
    snp_identifier="chr_pos",
    genome_build="hg19",
)
set_global_config(GLOBAL_CONFIG)

annotation_bundle = AnnotationBuilder(GLOBAL_CONFIG, AnnotationBuildConfig()).run(
    AnnotationSourceSpec(
        baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
        bed_paths="beds/*.bed",
    )
)

# If you want reusable query .annot.gz shards on disk, materialize them explicitly.
# reusable_bundle = run_bed_to_annot(
#     bed_files="beds/*.bed",
#     baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
#     output_dir="annotations/query_from_beds",
# )

sumstats = SumstatsMunger().run(
    RawSumstatsSpec(
        path="data/trait.tsv.gz",
        trait_name="trait",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
    ),
    MungeConfig(
        out_prefix="tutorial_outputs/trait",
        signed_sumstats_spec="BETA,0",
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated .sumstats.gz artifact on disk, keep the same
# GlobalConfig active before calling load_sumstats(). The loader warns because
# the original munge-time config is not recoverable from disk.
# sumstats = load_sumstats("tutorial_outputs/trait.sumstats.gz", trait_name="trait")

ref_panel = RefPanelLoader(GLOBAL_CONFIG).load(
    RefPanelSpec(
        backend="parquet_r2",
        r2_table_paths="r2/reference.@.parquet",
        maf_metadata_paths="r2/reference_metadata.@.tsv.gz",
        chromosomes=tuple(annotation_bundle.chromosomes),
        genome_build="hg19",
        ref_panel_snps_path="filters/reference_universe.tsv.gz",
    )
)

ldscore_result = LDScoreCalculator().run(
    annotation_bundle=annotation_bundle,
    ref_panel=ref_panel,
    ldscore_config=LDScoreConfig(
        ld_wind_cm=1.0,
        regression_snps_path="filters/hapmap3.tsv.gz",
    ),
    global_config=GLOBAL_CONFIG,
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
partitioned = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
)

partitioned.to_csv("tutorial_outputs/partitioned_h2.tsv", sep="\t", index=False)
print(annotation_bundle.query_columns)
print(ldscore_result.ldscore_table.head())
print(partitioned)
```

The Python workflow registers `GlobalConfig` once, then reuses it across the compatible helper functions and workflow classes. The resulting `AnnotationBundle`, `SumstatsTable`, and `LDScoreResult` objects carry frozen `config_snapshot` values, and the regression step raises `ConfigMismatchError` if you accidentally mix artifacts produced under incompatible `snp_identifier` or `genome_build` assumptions.

Within this design:

- `ref_panel_snps_path` belongs to `RefPanelSpec` and restricts the retained reference-panel rows
- `LDScoreCalculator.compute_chromosome()` intersects each chromosome-local annotation bundle with `ref_panel.load_metadata(chrom)`, so the LD-score compute universe is `B ∩ A'`
- `regression_snps_path` belongs to `LDScoreConfig` and further restricts the normalized `ldscore_table` rows to `B ∩ A' ∩ C`
- regression weights are embedded as `regr_weight`; there is no separate `.w.l2.ldscore.gz` artifact in the new default format

## CLI

The CLI supports BED-driven LD-score generation directly:

```bash
ldsc ldscore \
  --out tutorial_outputs/partitioned_ldscores \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed "beds/*.bed" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --ref-panel-snps-path filters/reference_universe.tsv.gz \
  --regression-snps-path filters/hapmap3.tsv.gz \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --ld-wind-cm 1.0
```

You can still materialize reusable query `.annot.gz` files explicitly:

```bash
ldsc annotate \
  --bed-files "beds/*.bed" \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds
```

The regression CLI has two current constraints that matter for partitioned runs:

- `ldsc ldscore` writes per-chromosome `.l2.ldscore.gz` files by default
- `ldsc partitioned-h2` loads one merged `.l2.ldscore.gz` table at a time and needs either `--query-columns` or an explicit annotation manifest

So the CLI regression example below assumes you already have:

- a merged LD-score table such as `tutorial_outputs/partitioned_ldscores.merged.l2.ldscore.gz`
- aggregate counts such as `tutorial_outputs/partitioned_ldscores.l2.M_5_50`
- query column names derived from your BED basenames, here `enhancers,promoters`

```bash
ldsc munge-sumstats \
  --sumstats data/trait.tsv.gz \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --signed-sumstats BETA,0 \
  --out tutorial_outputs/trait

ldsc partitioned-h2 \
  --sumstats tutorial_outputs/trait.sumstats.gz \
  --ldscore tutorial_outputs/partitioned_ldscores.merged.l2.ldscore.gz \
  --counts tutorial_outputs/partitioned_ldscores.l2.M_5_50 \
  --count-kind m_5_50 \
  --query-columns enhancers,promoters \
  --out tutorial_outputs/partitioned_h2
```

If you have a legacy separate weight file, add `--w-ld <path>` to the regression command. Otherwise leave it unset and let the loader use the embedded `regr_weight` column.
