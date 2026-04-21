# Partitioned LDSC

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
The bundled `baseline_v1.2` annotations are hg19-based, so the parquet example
uses `genome_build="hg19"` to align raw parquet coordinates to the annotation bundle.
The workflow also accepts `hg37` and `GRCh37` as aliases for `hg19`, and
`GRCh38` as an alias for `hg38`; outputs always normalize back to canonical
`hg19` or `hg38`.

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
    get_global_config,
    load_sumstats,
    run_bed_to_annot,
    set_global_config,
)

set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
global_config = get_global_config()

run_bed_to_annot(
    bed_files="beds/*.bed",
    baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
    output_dir="annotations/query_from_beds",
)

annotation_bundle = AnnotationBuilder(global_config, AnnotationBuildConfig()).run(
    AnnotationSourceSpec(
        baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
        query_annot_paths="annotations/query_from_beds/query.@.annot.gz",
    )
)
# `run()` resolves the `@` tokens and automatically bundles the per-chromosome
# annotation shards into one in-memory bundle.
# `run_bed_to_annot()` writes the projected query shards as `query.<chrom>.annot.gz`.

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
)

# If you already have a curated .sumstats.gz artifact on disk, load it directly:
# sumstats = load_sumstats("tutorial_outputs/trait.sumstats.gz", trait_name="trait")

ref_panel = RefPanelLoader(global_config).load(
    RefPanelSpec(
        backend="parquet_r2",
        r2_table_paths="r2/reference.*.parquet",
        chromosomes=tuple(annotation_bundle.chromosomes),
    )
)

ldscore_result = LDScoreCalculator().run(
    annotation_bundle=annotation_bundle,
    ref_panel=ref_panel,
    ldscore_config=LDScoreConfig(ld_wind_cm=1.0),
    global_config=global_config,
)

runner = RegressionRunner(regression_config=RegressionConfig())
partitioned = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
)

partitioned.to_csv("tutorial_outputs/partitioned_h2.tsv", sep="\t", index=False)
print(partitioned)
```

The Python workflow registers `GlobalConfig` once, then reuses it across the
compatible helper functions and workflow classes.

## CLI

```bash
ldsc annotate \
  --bed-files "beds/*.bed" \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --output-dir annotations/query_from_beds

ldsc ldscore \
  --out tutorial_outputs/partitioned_ldscores \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot "annotations/query_from_beds/query.@.annot.gz" \
  --r2-table "r2/reference.*.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --ld-wind-cm 1.0

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
  --ldscore tutorial_outputs/partitioned_ldscores.l2.ldscore.gz \
  --w-ld tutorial_outputs/partitioned_ldscores.w.l2.ldscore.gz \
  --counts tutorial_outputs/partitioned_ldscores.M_5_50 \
  --count-kind m_5_50 \
  --annotation-manifest tutorial_outputs/partitioned_ldscores.annotation_groups.tsv \
  --out tutorial_outputs/partitioned_h2
```

The default CLI summary is intentionally query-focused: one output table with one row per query annotation.
If `--output-dir` does not exist yet during annotation projection, the workflow warns once and creates it automatically.

For the Python workflow-layer API, the annotation bundle already provides SNP metadata such as
`CHR`, `POS`, `SNP`, and `CM`, so the partitioned example only needs the parquet `R2` tables.
If your workflow depends on separate frequency metadata such as `MAF`, load that through the
reference panel spec as well.
