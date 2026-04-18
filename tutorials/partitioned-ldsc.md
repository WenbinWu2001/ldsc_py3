# Partitioned LDSC

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/reference.1.parquet`, and `r2/reference_metadata.1.tsv.gz`.
The bundled `baseline_v1.2` annotations are hg19-based, so the parquet example
uses `genome_build="hg19"` to align raw parquet coordinates to the annotation bundle.
The workflow also accepts `hg37` and `GRCh37` as aliases for `hg19`, and
`GRCh38` as an alias for `hg38`; outputs always normalize back to canonical
`hg19` or `hg38`.

Path-token rules used in this tutorial:

- use `@` for chromosome suites such as `baseline.@`
- use globs when the filenames do not follow the simple chromosome-suffix convention
- scalar inputs still resolve to exactly one file
- output prefixes remain literal destinations

## Python API

```python
from ldsc import (
    AnnotationBuildConfig,
    AnnotationBuilder,
    AnnotationSourceSpec,
    CommonConfig,
    LDScoreCalculator,
    LDScoreConfig,
    MungeConfig,
    RawSumstatsSpec,
    RefPanelLoader,
    RefPanelSpec,
    RegressionConfig,
    RegressionRunner,
    SumstatsMunger,
    run_bed_to_annot,
)

common = CommonConfig(snp_identifier="chr_pos", genome_build="hg19")

run_bed_to_annot(
    bed_files="beds/*.bed",
    baseline_annot_dir="annotations/baseline_chr",
    output_dir="annotations/query_from_beds",
)

annotation_bundle = AnnotationBuilder(common, AnnotationBuildConfig()).run(
    AnnotationSourceSpec(
        baseline_annot_paths="annotations/baseline_chr/baseline.@",
        query_annot_paths="annotations/query_from_beds/query.@",
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
    common,
)

ref_panel = RefPanelLoader(common).load(
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
    common_config=common,
)

runner = RegressionRunner(common, RegressionConfig())
partitioned = runner.estimate_partitioned_h2_batch(
    sumstats,
    ldscore_result,
    annotation_bundle,
)

partitioned.to_csv("tutorial_outputs/partitioned_h2.tsv", sep="\t", index=False)
print(partitioned)
```

## CLI

```bash
ldsc annotate \
  --bed-files "beds/*.bed" \
  --baseline-annot-dir annotations/baseline_chr \
  --output-dir annotations/query_from_beds

ldsc ldscore \
  --out tutorial_outputs/partitioned_ldscores \
  --baseline-annot "annotations/baseline_chr/baseline.@" \
  --query-annot "annotations/query_from_beds/query.@" \
  --r2-table "r2/reference.*.parquet" \
  --frqfile "r2/reference_metadata.@" \
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

For the Python workflow-layer API, the annotation bundle already provides SNP metadata such as
`CHR`, `POS`, `SNP`, and `CM`, so the partitioned example only needs the parquet `R2` tables.
If your workflow depends on separate frequency metadata such as `MAF`, load that through the
reference panel spec as well.
