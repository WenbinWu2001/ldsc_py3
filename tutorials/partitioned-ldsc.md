# Partitioned LDSC

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
The parquet R2 files are expected to use the canonical six-column schema
(`CHR`, `POS_1`, `POS_2`, `R2`, `SNP_1`, `SNP_2`) with row-group statistics.
The paired metadata sidecar is required; it defines the raw reference-panel SNP
universe, while the parquet pair rows are queried only for LD values.
The bundled `baseline_v1.2` annotations are hg19-based, so the parquet example uses `genome_build="hg19"` to align parquet coordinates to the annotation bundle.
The workflow also accepts `hg37` and `GRCh37` as aliases for `hg19`, and `GRCh38` as an alias for `hg38`; outputs always normalize back to canonical `hg19` or `hg38`.

Path-token rules used in this tutorial:

- use `@` for chromosome suites such as `baseline.@.annot.gz`
- use globs when the filenames do not follow the simple chromosome-suffix convention
- scalar inputs still resolve to exactly one file
- output directories remain literal destinations

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
    GlobalConfig,
    LDScoreCalculator,
    LDScoreConfig,
    LDScoreOutputConfig,
    MungeConfig,
    RefPanelLoader,
    RefPanelConfig,
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
    AnnotationBuildConfig(
        baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
        query_annot_bed_paths="beds/*.bed",
    )
)

# If you want reusable query .annot.gz shards on disk, materialize them explicitly.
# reusable_bundle = run_bed_to_annot(
#     query_annot_bed_paths="beds/*.bed",
#     baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
#     output_dir="annotations/query_from_beds",
# )

sumstats = SumstatsMunger().run(
    MungeConfig(
        sumstats_path="data/trait.tsv.gz",
        trait_name="trait",
        column_hints={"snp": "SNP", "a1": "A1", "a2": "A2", "p": "P", "N_col": "N"},
        output_dir="tutorial_outputs/trait",
        signed_sumstats_spec="BETA,0",
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated .sumstats.gz artifact on disk, keep the same
# GlobalConfig active before calling load_sumstats(). The loader warns because
# the original munge-time config is not recoverable from disk.
# sumstats = load_sumstats("tutorial_outputs/trait/sumstats.sumstats.gz", trait_name="trait")

ref_panel = RefPanelLoader(GLOBAL_CONFIG).load(
    RefPanelConfig(
        backend="parquet_r2",
        r2_paths="r2/reference.@.parquet",
        metadata_paths="r2/reference_metadata.@.tsv.gz",
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
    output_config=LDScoreOutputConfig(output_dir="tutorial_outputs/partitioned_ldscores"),
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
partitioned = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
)

partitioned.to_csv("tutorial_outputs/partitioned_h2/partitioned_h2.tsv", sep="\t", index=False)
print(annotation_bundle.query_columns)
print(ldscore_result.baseline_table.head())
print(partitioned)
```

The Python workflow registers `GlobalConfig` once, then reuses it across the compatible helper functions and workflow classes. The resulting `AnnotationBundle`, `SumstatsTable`, and `LDScoreResult` objects carry frozen `config_snapshot` values, and the regression step raises `ConfigMismatchError` if you accidentally mix artifacts produced under incompatible `snp_identifier` or `genome_build` assumptions.

Within this design:

- `ref_panel_snps_path` belongs to `RefPanelConfig` and restricts the retained reference-panel rows
- `LDScoreCalculator.compute_chromosome()` intersects each chromosome-local annotation bundle with `ref_panel.load_metadata(chrom)`, so the LD-score compute universe is `B ∩ A'`; parquet pair rows are not scanned to define SNP presence
- `regression_snps_path` belongs to `LDScoreConfig` and further restricts the normalized `baseline_table` rows to `B ∩ A' ∩ C`
- regression weights are embedded as `regr_weight`; there is no separate `.w.l2.ldscore.gz` artifact in the new default format

## CLI

The CLI supports BED-driven LD-score generation directly:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/partitioned_ldscores \
  --baseline-annot-paths "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-paths "beds/*.bed" \
  --r2-paths "r2/reference.@.parquet" \
  --metadata-paths "r2/reference_metadata.@.tsv.gz" \
  --ref-panel-snps-path filters/reference_universe.tsv.gz \
  --regression-snps-path filters/hapmap3.tsv.gz \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --ld-wind-cm 1.0
```

You can still materialize reusable query `.annot.gz` files explicitly:

```bash
ldsc annotate \
  --query-annot-bed-paths "beds/*.bed" \
  --baseline-annot-paths "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds
```

The regression CLI consumes the LD-score result directory directly. It reads
baseline columns from `baseline.parquet`, query columns from `query.parquet`,
and counts from `manifest.json`.

```bash
ldsc munge-sumstats \
  --sumstats-path data/trait.tsv.gz \
  --snp SNP \
  --a1 A1 \
  --a2 A2 \
  --p P \
  --N-col N \
  --signed-sumstats BETA,0 \
  --output-dir tutorial_outputs/trait

ldsc partitioned-h2 \
  --sumstats-path tutorial_outputs/trait/sumstats.sumstats.gz \
  --ldscore-dir tutorial_outputs/partitioned_ldscores \
  --count-kind m_5_50 \
  --output-dir tutorial_outputs/partitioned_h2
```

The command writes `tutorial_outputs/partitioned_h2/partitioned_h2.tsv`.
