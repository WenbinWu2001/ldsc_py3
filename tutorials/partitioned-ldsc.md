# Partitioned LDSC

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

`partitioned-h2` requires explicit query annotations in the LD-score result
directory. Baseline-only LD-score directories, including synthetic all-ones
`base` outputs from ordinary `ldsc ldscore` runs, are valid for `h2` and `rg`
but are rejected by `partitioned-h2`.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
Package-built parquet R2 files use canonical pair columns (`CHR`, `POS_1`,
`POS_2`, `SNP_1`, `SNP_2`, `R2`) plus endpoint allele columns
(`A1_1`, `A2_1`, `A1_2`, `A2_2`) in allele-aware modes, with row-group
statistics. The paired metadata sidecar is required; it defines the raw
reference-panel SNP universe, while the parquet pair rows are queried only for
LD values. Six-column external/raw R2 inputs are supported only in base `rsid`
and `chr_pos` modes.
Package-built panels carry `ldsc:r2_bias` and `ldsc:n_samples` in parquet
schema metadata, so the examples omit R2 bias and sample-size arguments.
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
- existing owned workflow artifacts are refused before writing starts unless
  you pass `--overwrite` or `overwrite=True`; successful overwrites remove
  stale owned siblings not produced by the current configuration and preserve
  unrelated files

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
#     overwrite=True,  # also removes stale query shards outside the current chromosome set
# )

sumstats = SumstatsMunger().run(
    MungeConfig(
        raw_sumstats_file="data/trait.tsv.gz",
        trait_name="trait",
        # Most common columns are inferred automatically. If this file has
        # NEFF and you want to use it as the munger's N, pass this explicitly:
        # column_hints={"N_col": "NEFF"},
        # use_hm3_snps=True,  # packaged HM3 row restriction
        # target_genome_build="hg38",
        # liftover_chain_file="resources/liftover/hg19ToHg38.over.chain",
        output_dir="tutorial_outputs/trait",
        # overwrite=True,  # also removes stale unselected sumstats sibling formats
    ),
    global_config=GLOBAL_CONFIG,
)

# If you already have a curated sumstats artifact on disk, load it directly.
# Current parquet and .sumstats.gz artifacts recover config_snapshot from
# root metadata.json. Old package-written files without current provenance
# must be regenerated with the current LDSC package.
# sumstats = load_sumstats("tutorial_outputs/trait/sumstats.parquet", trait_name="trait")

ref_panel = RefPanelLoader(GLOBAL_CONFIG).load(
    RefPanelConfig(
        backend="parquet_r2",
        r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg19",
        chromosomes=tuple(annotation_bundle.chromosomes),
        use_hm3_ref_panel_snps=True,
    )
)

ldscore_result = LDScoreCalculator().run(
    annotation_bundle=annotation_bundle,
    ref_panel=ref_panel,
    ldscore_config=LDScoreConfig(
        ld_wind_cm=1.0,
        use_hm3_regression_snps=True,
    ),
    global_config=GLOBAL_CONFIG,
    output_config=LDScoreOutputConfig(
        output_dir="tutorial_outputs/partitioned_ldscores",
        # overwrite=True,  # also removes stale LD-score siblings not produced by this run
    ),
)

runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
partitioned_result = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
    include_full_partitioned_h2=True,
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

The Python workflow registers `GlobalConfig` once, then reuses it across the compatible helper functions and workflow classes. In-process results such as `AnnotationBundle`, `SumstatsTable` from `SumstatsMunger.run()`, and `LDScoreResult` carry frozen `config_snapshot` values, and the regression step raises `ConfigMismatchError` if you accidentally mix artifacts produced under incompatible `snp_identifier` or `genome_build` assumptions. A `SumstatsTable` loaded from a current disk artifact recovers this provenance from root `metadata.json`; old package-written sumstats without current provenance must be regenerated with the current LDSC package. `SumstatsMunger.run()` is also the implementation path behind `ldsc munge-sumstats` after CLI parsing, and it owns root `metadata.json`, fixed `sumstats.parquet` output by default, optional `sumstats.sumstats.gz` compatibility output, and diagnostics under `diagnostics/`. Workflow logs are preflighted audit files; returned `output_paths` mappings include data artifacts and the dropped-SNP audit sidecar, but not logs. For `munge-sumstats`, `ldscore`, `partitioned-h2`, and `annotate`, output directories represent coherent artifact families: no-overwrite runs reject any owned sibling, and successful overwrites delete stale owned siblings not produced by the current configuration.

Munged sumstats written by this workflow include canonical `CHR` and `POS`
columns. The raw munger accepts common coordinate headers such as `#CHROM`,
`CHROM`, `CHR`, `POS`, and `BP`, or explicit `--chr`/`--pos` flags. Leading
raw `##` metadata lines are skipped before the real header is parsed. In
`chr_pos`-family modes, downstream regression merges by the active effective
coordinate key (`CHR:POS` or `CHR:POS:<allele_set>`) rather than by the literal
rsID in `SNP`; `SNP` is treated as a label. Optional munger liftover is also
`chr_pos`-family behavior, runs after the source-build keep-list filter, changes
`CHR`/`POS` without rewriting `SNP` or allele sets, drops duplicate
source/target coordinate groups, and requires `--target-genome-build` plus one
method flag.
Drop counts are written to `diagnostics/sumstats.log`, examples appear only at `DEBUG`, and
row-level drops are audited in `diagnostics/dropped_snps/dropped.tsv.gz`; the metadata
sidecar is only the compatibility snapshot.

Within this design:

- `ref_panel_snps_file` or `use_hm3_ref_panel_snps` belongs to `RefPanelConfig` and restricts the retained reference-panel rows
- `LDScoreCalculator.compute_chromosome()` intersects each chromosome-local annotation bundle with `ref_panel.load_metadata(chrom)`, so the LD-score compute universe is `B ∩ A'`; parquet pair rows are not scanned to define SNP presence
- `regression_snps_file` or `use_hm3_regression_snps` belongs to `LDScoreConfig` and further restricts the normalized `baseline_table` rows to `B ∩ A' ∩ C`
- explicit reference-panel and regression SNP restriction files are identity-only filters: duplicate restriction keys collapse to one retained key, and non-identity columns such as `CM` or `MAF` are ignored
- regression-universe LD scores are embedded as `regression_ld_scores`; there is no separate `.w.l2.ldscore.gz` artifact in the new default format
- query `.annot` and BED inputs are accepted only when baseline annotations are supplied explicitly

## CLI

The CLI supports BED-driven LD-score generation directly:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/partitioned_ldscores \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg19" \
  --use-hm3-ref-panel-snps \
  --use-hm3-regression-snps \
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

Both CLI paths preflight their workflow logs with the scientific outputs:
`diagnostics/ldscore.log` under `tutorial_outputs/partitioned_ldscores` and
`diagnostics/annotate.log` under `annotations/query_from_beds`.

The annotate command is implemented in the public `ldsc.annotation_builder`
workflow module. Python code that wants the same parser behavior can call
`ldsc.annotation_builder.main(argv)`, while code that already has parsed
top-level CLI arguments can call `run_annotate_from_args(args)` directly. When
the workflow writes query shards, it also writes `diagnostics/annotate.log`
under the output directory.

The regression CLI consumes the LD-score result directory directly. It reads
baseline columns from `ldscore.baseline.parquet`, query columns from `ldscore.query.parquet`,
and counts from root `metadata.json`. Both parquet files stay flat, but their row
groups are chromosome-aligned and listed in the metadata for targeted reads.
If root `metadata.json` has an empty `query_columns` list, use `ldsc h2`/`ldsc rg`
instead of `ldsc partitioned-h2`.

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --trait-name trait \
  --use-hm3-snps \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --output-dir tutorial_outputs/trait

# Add explicit repair flags only when --infer-only reports that they are needed,
# for example --N-col NEFF if that is appropriate for the analysis.

# Optional if downstream LD scores/reference panels are hg38:
#   --target-genome-build hg38 \
#   --liftover-chain-file resources/liftover/hg19ToHg38.over.chain

ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.parquet \
  --ldscore-dir tutorial_outputs/partitioned_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/partitioned_h2
```

The command writes `tutorial_outputs/partitioned_h2/partitioned_h2.tsv` and
`tutorial_outputs/partitioned_h2/diagnostics/partitioned-h2.log`.
The summary columns are documented in
[partitioned-h2-results.md](../docs/current/partitioned-h2-results.md).
If any partitioned-h2 owned output already exists, including a stale
`diagnostics/query_annotations/` tree from an earlier per-query run, the command fails
before writing; add `--overwrite` only when replacing the previous summary is
intentional.

To also materialize one result folder per query annotation, add
`--write-per-query-results`:

```bash
ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/sumstats.parquet \
  --ldscore-dir tutorial_outputs/partitioned_ldscores \
  --count-kind common \
  --output-dir tutorial_outputs/partitioned_h2 \
  --write-per-query-results
```

This keeps the aggregate `partitioned_h2.tsv` and adds
`diagnostics/query_annotations/manifest.tsv` plus sanitized query folders such as
`diagnostics/query_annotations/0001_enhancer_a/`. Each query folder contains its one-row
`partitioned_h2.tsv`, the fitted baseline-plus-query `partitioned_h2_full.tsv`,
and `metadata.json` with the original query annotation name.
If you later rerun the same output directory without `--write-per-query-results`
and pass `--overwrite`, the old `diagnostics/query_annotations/` tree is removed after the
new aggregate summary is written.
