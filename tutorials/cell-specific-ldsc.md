# Cell-Specific LDSC

Goal: estimate cell-specific enrichment by running partitioned LDSC with one query annotation per cell type.

In this package, cell-specific LDSC is the `partitioned-h2` workflow applied to cell-type annotations. Baseline annotations stay in the model as covariates, and each cell-type query column is tested in a baseline-plus-one-query model through `RegressionRunner.estimate_partitioned_h2_batch()`.

The examples below assume chromosome-pattern inputs such as `baseline.@.annot.gz`, `cell_type_beds/*.bed`, `reference.@.parquet`, and `reference_metadata.@.tsv.gz`.

## Python API

```python
from ldsc import (
    AnnotationBuildConfig,
    AnnotationBuilder,
    AnnotationSourceSpec,
    GlobalConfig,
    LDScoreCalculator,
    LDScoreConfig,
    RefPanelLoader,
    RefPanelSpec,
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
    AnnotationSourceSpec(
        baseline_annot_paths="annotations/baseline_chr/baseline.@.annot.gz",
        bed_paths="annotations/cell_type_beds/*.bed",
    )
)

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

sumstats = load_sumstats("tutorial_outputs/trait.sumstats.gz", trait_name="trait")

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

## CLI

First compute baseline-plus-cell-type LD scores. This can project BED files directly without materializing intermediate query `.annot.gz` files:

```bash
ldsc ldscore \
  --out tutorial_outputs/cell_specific_ldscores \
  --baseline-annot "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed "annotations/cell_type_beds/*.bed" \
  --r2-table "r2/reference.@.parquet" \
  --frqfile "r2/reference_metadata.@.tsv.gz" \
  --ref-panel-snps-path filters/reference_universe.tsv.gz \
  --regression-snps-path filters/hapmap3.tsv.gz \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --ld-wind-cm 1.0
```

Then run partitioned h2 over the cell-type query columns:

```bash
ldsc partitioned-h2 \
  --sumstats tutorial_outputs/trait.sumstats.gz \
  --trait-name trait \
  --ldscore tutorial_outputs/cell_specific_ldscores.merged.l2.ldscore.gz \
  --counts tutorial_outputs/cell_specific_ldscores.l2.M_5_50 \
  --count-kind m_5_50 \
  --query-columns neuron,astrocyte,microglia \
  --out tutorial_outputs/cell_specific_ldsc
```

If the LD-score run wrote an annotation manifest, prefer it over manually listing query names:

```bash
ldsc partitioned-h2 \
  --sumstats tutorial_outputs/trait.sumstats.gz \
  --ldscore tutorial_outputs/cell_specific_ldscores.merged.l2.ldscore.gz \
  --counts tutorial_outputs/cell_specific_ldscores.l2.M_5_50 \
  --count-kind m_5_50 \
  --annotation-manifest tutorial_outputs/cell_specific_ldscores.annotation_manifest.tsv \
  --out tutorial_outputs/cell_specific_ldsc
```

The output file is `tutorial_outputs/cell_specific_ldsc.partitioned_h2.tsv`. Its key columns are `query_annotation`, `coefficient`, `coefficient_p`, `category_h2`, `proportion_h2`, and `enrichment`.
