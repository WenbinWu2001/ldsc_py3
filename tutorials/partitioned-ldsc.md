# Partitioned LDSC

Last updated on: 2026-09-14

Goal: run partitioned LDSC in the refactored package by building query annotations, computing baseline-plus-query LD scores, and fitting one partitioned model per query annotation.

After writing partitioned or quantile results, use the [plotting results manual](plotting-results.md) for the functional-enrichment, query-evidence, and continuous-annotation quantile plots.

`partitioned-h2` accepts both baseline-only and query-annotation LD-score
directories. Baseline-only runs fit the functional-category model and keep the
complete-model artifacts at the result root. Query-annotation runs fit one
baseline-plus-query model per query and also write the per-query diagnostics
tree described below.

The examples below assume chromosome-pattern inputs such as `annotations/baseline.1.annot.gz`, `r2/chr1_r2.parquet`, and `r2/chr1_meta.tsv.gz`.
Package-built parquet R2 files use canonical `IDX_1`, `IDX_2`, `R2`, and `SIGN_R` columns bound to the paired metadata sidecar, with row-group statistics. The paired metadata sidecar is required; it defines the raw
reference-panel SNP universe, while the parquet pair rows are queried only for
LD values. External R2 formats are not supported by this workflow.
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
- existing current-contract workflow artifacts are refused before writing starts unless
  you pass `--overwrite` or `overwrite=True`; successful overwrites remove
  stale owned siblings not produced by the current configuration, preserve
  unrelated files, and ignore removed legacy root diagnostics

Resolution behavior:

- there is no separate `*_chr` public argument anymore; one argument now handles both shared inputs and chromosome-sharded inputs
- group inputs may expand to many files through a glob or an `@` suite token
- chromosome scope comes from validated contents; whole-genome sources are scanned once in bounded chunks and normalized into private chromosome artifacts
- chunk sizing accounts for the baseline/query files aligned together; identity cleanup remains global across chromosomes and reads metadata separately from numeric staging
- if multiple files contribute annotation columns for the same chromosome, their SNP rows must align exactly and their annotation column names must be unique

Query annotations require explicit baseline annotations. The LD-score workflow
can synthesize an all-ones `base` column only when both baseline and query
inputs are omitted for ordinary unpartitioned LD scores; it does not use that
synthetic path for partitioned/query LDSC.

Annotation preparation remains serial before chromosome computation; `--threads` applies to the subsequent LD-score workers. At INFO, `diagnostics/annotate.log` or `diagnostics/ldscore.log` marks input reading, SNP identity checks and chromosome preparation, and preparation completion with retained counts and elapsed time. CM/MAF notices appear once per file read, and intentional gene exclusions use one line per gene set. See [preparation logging](../docs/current/workflow-logging.md#annotation-preparation).

## Python API

The writing workflow returns an `LDScoreSource` with shared baseline values and saved query paths. Completed query tables are released after writing. Regression reads only the columns it needs, including selections spanning several saved query files.

```python
from ldsc import (
    GlobalConfig, RegressionConfig, RegressionRunner,
    load_sumstats, run_ldscore, set_global_config,
)

GLOBAL_CONFIG = GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
set_global_config(GLOBAL_CONFIG)
ldscore_dir = "tutorial_outputs/partitioned_ldscores"

source = run_ldscore(
    baseline_annot_sources="annotations/baseline_chr/baseline.@.annot.gz",
    query_annot_bed_sources="beds/*.bed",
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
    output_dir="tutorial_outputs/partitioned_h2",
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

The CLI supports BED-driven LD-score generation directly:

```bash
ldsc ldscore \
  --output-dir tutorial_outputs/partitioned_ldscores \
  --baseline-annot-sources "annotations/baseline_chr/baseline.@.annot.gz" \
  --query-annot-bed-sources "beds/*.bed" \
  --r2-dir "r2_ref_panel_1kg30x_1cM_hm3/hg19" \
  --ref-panel-snps-file filters/reference_universe.tsv.gz \
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
baseline columns from `ldscore.baseline.parquet`, query columns from the files in `metadata.json.query_batches`,
and counts from root `metadata.json`. Every Parquet file is genome-wide, with chromosome-aligned row groups listed in metadata. An empty `query_columns` list selects one baseline-only functional-category model in `partitioned-h2`; use `h2` for an unpartitioned heritability estimate or `rg` for genetic correlation.

```bash
ldsc munge-sumstats \
  --raw-sumstats-file data/trait.tsv.gz \
  --trait-name trait \
  --snp-identifier chr_pos \
  --source-genome-build hg19 \
  --output-genome-build hg19 \
  --output-dir tutorial_outputs/trait

# Add explicit repair flags only when --infer-only reports that they are needed,
# for example --N-col NEFF if that is appropriate for the analysis.

# Optional if downstream LD scores/reference panels are hg38:
#   --output-genome-build hg38 \
#   --liftover-chain-file resources/liftover/hg19ToHg38.over.chain
# Packaged HM3 mapping is automatic when builds differ; the chain overrides it.

ldsc partitioned-h2 \
  --sumstats-file tutorial_outputs/trait/trait.parquet \
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

For query-annotation runs, the same command keeps the aggregate
`partitioned_h2.tsv` and adds
`diagnostics/query_annotations/manifest.tsv` plus sanitized query folders such as
`diagnostics/query_annotations/0001_enhancer_a/`. Each query folder contains its one-row
`partitioned_h2.tsv`, the fitted baseline-plus-query `partitioned_h2_full.tsv`,
`coefficient_delete_values.parquet`, and `metadata.json` with the original query annotation name.
The retired `--write-per-query-results` flag is rejected; omit it. Baseline-only runs do not create
the per-query tree and keep their complete-model artifacts at the result root.
