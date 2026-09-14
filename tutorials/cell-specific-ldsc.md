# Cell-Specific LDSC

Last updated on: 2026-09-14

Goal: estimate cell-specific enrichment by running partitioned LDSC with one query annotation per cell type.

After writing the aggregate query result, use the [plotting results manual](plotting-results.md) to create and interpret the horizontal nominal-p-value summary across the separate baseline-conditional fits.

In this package, cell-specific LDSC is the `partitioned-h2` workflow applied to cell-type annotations. Baseline annotations stay in the model as covariates, and each cell-type query column is tested in a baseline-plus-one-query model through `RegressionRunner.estimate_partitioned_h2_batch()`.

Cell-type query annotations require explicit baseline annotations. The
synthetic all-ones `base` annotation is reserved for ordinary unpartitioned
LD-score generation when no query inputs are present. `partitioned-h2` treats baseline-only directories as one joint functional-category fit; they have no cell-type query tests.

The examples below assume chromosome-pattern inputs such as
`baseline.@.annot.gz`, `cell_type_beds/*.bed`, and a package-built
build-specific R2 directory. Package-built R2 parquet files store
`ldsc:r2_bias` and `ldsc:n_samples` in schema metadata, so the examples omit
R2 bias and sample-size arguments.

Output directories are literal destinations. Missing directories are created,
existing directories are reused, and existing current-contract workflow
artifacts are refused before writing unless you pass `--overwrite` or
`overwrite=True`. Successful overwrites remove stale owned siblings not
produced by the current configuration and preserve unrelated files. CLI workflow
logs are part of that preflight policy, but they are audit files rather than
returned data artifacts. Removed legacy root diagnostic names are ignored by
current preflight and cleanup.

For `chr_pos` workflows, `genome_build="auto"` can infer hg19/hg38 and
0-based/1-based coordinates during annotation or LD-score loading. The same
logic is public in Python as `infer_chr_pos_build()` and
`resolve_chr_pos_table()` from `ldsc`; the CLI keeps it under workflow flags such
as `--genome-build auto` rather than a separate command.

## Python API

The writing workflow returns an `LDScoreSource` with shared baseline values and saved query paths. Completed query tables are released after writing. Regression reads only the columns it needs, including selections spanning several saved query files.

```python
from ldsc import (
    GlobalConfig, RegressionConfig, RegressionRunner,
    load_sumstats, run_ldscore, set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg19")
set_global_config(GLOBAL_CONFIG)
ldscore_dir = "tutorial_outputs/cell_specific_ldscores"

source = run_ldscore(
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot_bed_sources="annotations/cell_type_beds/*.bed",
    r2_dir="r2_ref_panel_1kg30x_1cM_hm3/hg19",
    ld_wind_cm=1.0,
    query_batch_size=1000,
    threads=1,
    output_dir=ldscore_dir,
)
# Later sessions can reopen the handle with load_ldscore_from_dir(ldscore_dir).
sumstats = load_sumstats("tutorial_outputs/trait/trait.parquet", trait_name="trait")
runner = RegressionRunner(global_config=GLOBAL_CONFIG, regression_config=RegressionConfig())
result = runner.estimate_partitioned_h2_batch(
    sumstats, source,
    query_batch_size=1000,
    output_dir="tutorial_outputs/cell_specific_h2",
    metadata={"ldscore_dir": ldscore_dir, "trait_name": "trait"},
)
print(result.summary)
print(result.per_query_artifacts)  # Persistent category/delete-value/metadata paths.
```

Run `munge-sumstats` first if the trait is still in a raw format; see [heritability estimates](heritability-estimates.md). Current curated Parquet inputs carry identity/build provenance in their footer. The effective SNP key controls alignment; reference contributors and regression/output SNPs remain separate universes.

For 1,000 pathways, each query still fits separately against shared baseline categories. The default `query_batch_size=1000` permits all 1,000 queries in one execution batch. Reducing it bounds active query workspace and writes multiple numbered genome-wide query files, with an ordered `query_batches` manifest in root metadata. Direct calculation repeats reference/genotype work between batches; indexed calculation reuses one chromosome operator per worker through its batches. Direct and indexed `threads` default to 1 and are capped at the chromosome count. Freed pages can remain reserved by the allocator, so RSS need not immediately fall. See `LDScoreCalculator.run` and `LDScoreSource.read_queries` in the [memory design](../docs/current/annotation-memory-design.md).

Generation and regression batch widths are independent. To inspect values explicitly, call `source.read_queries(["query_name"])`; reads preserve requested column order and have no cache or width cap. Current readers require the `query_batches` manifest; regenerate older directories. For the small prepared-input Python route that writes nothing, see [LD-score calculation without writes](ld-score-calculation.md#small-python-calculations-without-writes).

For reusable annotations, use `with run_annotate(..., output_dir=...) as bundle:`. For explicit low-level preparation, use `with AnnotationBuilder(config).run(source_config, output_dir=...) as bundle:` and finish all borrowers before closure. Returned standalone bundles reference saved query outputs and may depend on original baseline inputs. See the [developer memory design](../docs/current/annotation-memory-design.md).

## CLI

Intercepts are estimated by default. To request the standard fixed intercepts, add `--intercept-h2 1`. This fixes the h2 intercept to 1. The removed `--no-intercept` shortcut is rejected; see the [legacy flag map](../docs/current/legacy-cli-flag-map.md#regression-intercept-consolidation).

First compute baseline-plus-cell-type LD scores. This can project BED files directly without materializing intermediate query `.annot.gz` files:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/cell_specific_ldscores \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "annotations/cell_type_beds/*.bed" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg19" \
  --ref-panel-snps-file filters/reference_universe.tsv.gz \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --common-maf-min 0.05 \
  --ld-wind-cm 1.0
```

When reusable query `.annot.gz` shards are useful, use `ldsc annotate` or
`run_annotate(...)`; both are public `ldsc.annotation_builder` workflow
entry points and return the same `AnnotationBundle` shape used above.

Then run partitioned h2 over the cell-type query columns:

```bash
ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/trait.parquet \
  --trait-name trait \
  --ldscore-dir tutorial_outputs/cell_specific_ldscores \
  --count-kind common \
  --summary-sort-by enrichment-p \
  --output-dir tutorial_outputs/cell_specific_ldsc
```

The regression reads query annotation columns from
`tutorial_outputs/cell_specific_ldscores/metadata.json` and
the query files in its ordered `query_batches` manifest. One batch uses `ldscore.query.parquet`; multiple batches use numbered query files. Each file has chromosome-aligned row groups. `ldscore --query-batch-size` controls generation width, while `partitioned-h2 --query-batch-size` independently controls read width; both default to 1000. Each pathway remains a separate regression model. The output file is
`tutorial_outputs/cell_specific_ldsc/partitioned_h2.tsv`, with
`diagnostics/partitioned-h2.log` under the same directory. Its key columns are
`category`, `prop_snps`, `prop_h2`, `enrichment`, `enrichment_p`,
`coefficient`, and `coefficient_p`.
For full column definitions, see
[partitioned-h2-results.md](../docs/current/partitioned-h2-results.md).
For query-annotation runs, the command writes by default
`tutorial_outputs/cell_specific_ldsc/diagnostics/query_annotations/manifest.tsv` and one
sanitized folder per cell-type query annotation. Each folder contains the
one-row query summary, the baseline-plus-query `partitioned_h2_full.tsv`,
`coefficient_delete_values.parquet`, and `metadata.json` with the original annotation name.
If the partitioned summary already exists, `ldsc partitioned-h2` fails before
writing; the same is true for `diagnostics/partitioned-h2.log` and any stale
`diagnostics/query_annotations/` tree. Add `--overwrite` only when replacing it is
intentional. The retired `--write-per-query-results` flag is rejected; omit it.
