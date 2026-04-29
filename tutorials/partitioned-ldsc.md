# Partitioned LDSC

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
The parquet R2 files are expected to use the canonical six-column schema
(`CHR`, `POS_1`, `POS_2`, `R2`, `SNP_1`, `SNP_2`) with row-group statistics.
The paired metadata sidecar is required; it defines the raw reference-panel SNP
universe, while the parquet pair rows are queried only for LD values.
The bundled `baseline_v1.2` annotations are hg19-based, so the parquet example uses `genome_build="hg19"` to align parquet coordinates to the annotation bundle.
The workflow also accepts `hg37` and `GRCh37` as aliases for `hg19`, and `GRCh38` as an alias for `hg38`; outputs always normalize back to canonical `hg19` or `hg38`.
If you are using `chr_pos` inputs and need the package to infer hg19/hg38 and
0-based/1-based coordinates, set `genome_build="auto"` in `GlobalConfig` or pass
`--genome-build auto` to `ldsc annotate` or `ldsc ldscore`. Programmatic
preflight is available with `infer_chr_pos_build()` and
`resolve_chr_pos_table()` from `ldsc`; there is no standalone CLI inference
command.

Path-token rules used in this tutorial:

- use `@` for chromosome suites such as `baseline.@.annot.gz`
- use globs when the filenames do not follow the simple chromosome-suffix convention
- scalar inputs still resolve to exactly one file
- output directories remain literal destinations
- missing output directories are created and existing directories are reused
- existing fixed files are refused before writing starts unless you pass
  `--overwrite` or `overwrite=True`

Resolution behavior:

- there is no separate `*_chr` public argument anymore; one argument now handles both shared inputs and chromosome-sharded inputs
- group inputs may expand to many files through a glob or an `@` suite token
- if those files encode chromosome labels in their filenames, the workflow selects the subset for the active chromosome before loading them
- otherwise the workflow loads the matched files and keeps only rows whose `CHR` value matches the active chromosome
- if multiple files contribute annotation columns for the same chromosome, their SNP rows must align exactly and their annotation column names must be unique

Query annotations require explicit baseline annotations. The LD-score workflow
can synthesize an all-ones `base` column only when both baseline and query
inputs are omitted for ordinary unpartitioned LD scores; it does not use that
synthetic path for partitioned/query LDSC.

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
    PartitionedH2DirectoryWriter,
    PartitionedH2OutputConfig,
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
        baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
        query_annot_bed_sources="beds/*.bed",
    )
)

# If you want reusable query .annot.gz shards on disk, materialize them explicitly.
# reusable_bundle = run_bed_to_annot(
#     query_annot_bed_sources="beds/*.bed",
#     baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
#     output_dir="annotations/query_from_beds",
#     overwrite=True,  # enable only when intentionally replacing query shards
# )

sumstats = SumstatsMunger().run(
    MungeConfig(
        sumstats_file="data/trait.tsv.gz",
        trait_name="trait",
        column_hints={
            "snp": "ID",
            "chr": "#CHROM",
            "pos": "POS",
            "a1": "EA",
            "a2": "NEA",
            "p": "PVAL",
            "N_col": "NEFF",
            "info": "IMPINFO",
        },
        output_dir="tutorial_outputs/trait",
        signed_sumstats_spec="BETA,0",
        # overwrite=True,  # enable only when intentionally replacing sumstats/log files
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated .sumstats.gz artifact on disk, load it directly.
# Current artifacts recover config_snapshot from sumstats.metadata.json. Older
# files without that sidecar warn and use config_snapshot=None.
# sumstats = load_sumstats("tutorial_outputs/trait/sumstats.sumstats.gz", trait_name="trait")

ref_panel = RefPanelLoader(GLOBAL_CONFIG).load(
    RefPanelConfig(
        backend="parquet_r2",
        r2_sources="r2/reference.@.parquet",
        metadata_sources="r2/reference_metadata.@.tsv.gz",
        chromosomes=tuple(annotation_bundle.chromosomes),
        genome_build="hg19",
        r2_bias_mode="unbiased",
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
        output_dir="tutorial_outputs/partitioned_ldscores",
        # overwrite=True,  # enable only when intentionally replacing LD-score outputs
    ),
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
partitioned_result = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
    include_model_categories=True,
)
partitioned = partitioned_result.summary

PartitionedH2DirectoryWriter().write(
    partitioned,
    PartitionedH2OutputConfig(
        output_dir="tutorial_outputs/partitioned_h2",
        # overwrite=True,
        write_per_query_results=True,
    ),
    per_query_category_tables=partitioned_result.per_query_category_tables,
    metadata={"trait_name": "trait", "count_kind": "common", "ldscore_dir": "tutorial_outputs/partitioned_ldscores"},
    per_query_metadata=partitioned_result.per_query_metadata,
)
print(annotation_bundle.query_columns)
print(ldscore_result.baseline_table.head())
print(partitioned)
```

The Python workflow registers `GlobalConfig` once, then reuses it across the compatible helper functions and workflow classes. In-process results such as `AnnotationBundle`, `SumstatsTable` from `SumstatsMunger.run()`, and `LDScoreResult` carry frozen `config_snapshot` values, and the regression step raises `ConfigMismatchError` if you accidentally mix known artifacts produced under incompatible `snp_identifier` or `genome_build` assumptions. A `SumstatsTable` loaded from a current disk artifact recovers this provenance from `sumstats.metadata.json`; older artifacts without the sidecar have unknown provenance (`config_snapshot=None`) and do not trigger sumstats-side compatibility validation.

Munged sumstats written by this workflow include canonical `CHR` and `POS`
columns. The raw munger accepts common coordinate headers such as `#CHROM`,
`CHROM`, `CHR`, `POS`, and `BP`, or explicit `--chr`/`--pos` flags. Leading
raw `##` metadata lines are skipped before the real header is parsed. In
`chr_pos` mode, downstream regression merges by normalized `CHR:POS` rather
than by the literal rsID in `SNP`.

Within this design:

- `ref_panel_snps_file` belongs to `RefPanelConfig` and restricts the retained reference-panel rows
- `LDScoreCalculator.compute_chromosome()` intersects each chromosome-local annotation bundle with `ref_panel.load_metadata(chrom)`, so the LD-score compute universe is `B ∩ A'`; parquet pair rows are not scanned to define SNP presence
- `regression_snps_file` belongs to `LDScoreConfig` and further restricts the normalized `baseline_table` rows to `B ∩ A' ∩ C`
- regression weights are embedded as `regr_weight`; there is no separate `.w.l2.ldscore.gz` artifact in the new default format
- query `.annot` and BED inputs are accepted only when baseline annotations are supplied explicitly

## CLI

The CLI supports BED-driven LD-score generation directly:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/partitioned_ldscores \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --r2-sources "r2/reference.@.parquet" \
  --metadata-sources "r2/reference_metadata.@.tsv.gz" \
  --r2-bias-mode unbiased \
  --ref-panel-snps-file filters/reference_universe.tsv.gz \
  --regression-snps-file filters/hapmap3.tsv.gz \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

You can still materialize reusable query `.annot.gz` files explicitly:

```bash
ldsc annotate \
  --query-annot-bed-sources "beds/*.bed" \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds
```

The regression CLI consumes the LD-score result directory directly. It reads
baseline columns from `baseline.parquet`, query columns from `query.parquet`,
and counts from `manifest.json`. Both parquet files stay flat, but their row
groups are chromosome-aligned and listed in the manifest for targeted reads.

```bash
ldsc munge-sumstats \
  --sumstats-file data/trait.tsv.gz \
  --snp ID \
  --chr '#CHROM' \
  --pos POS \
  --a1 EA \
  --a2 NEA \
  --p PVAL \
  --N-col NEFF \
  --info IMPINFO \
  --signed-sumstats BETA,0 \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --output-dir tutorial_outputs/trait

ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.sumstats.gz \
  --ldscore-dir tutorial_outputs/partitioned_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/partitioned_h2
```

The command writes `tutorial_outputs/partitioned_h2/partitioned_h2.tsv`.
If that TSV already exists, the command fails before writing; add `--overwrite`
only when replacing the previous summary is intentional.

To also materialize one result folder per query annotation, add
`--write-per-query-results`:

```bash
ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.sumstats.gz \
  --ldscore-dir tutorial_outputs/partitioned_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/partitioned_h2 \
  --write-per-query-results
```

This keeps the aggregate `partitioned_h2.tsv` and adds
`query_annotations/manifest.tsv` plus sanitized query folders such as
`query_annotations/0001_enhancer_a/`. Each query folder contains its one-row
`partitioned_h2.tsv`, the fitted baseline-plus-query `model_categories.tsv`,
and `metadata.json` with the original query annotation name.
