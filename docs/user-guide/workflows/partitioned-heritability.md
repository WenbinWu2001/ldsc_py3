# Partitioned heritability

Last updated on: 2026-10-28

Partitioned heritability asks how SNP heritability is distributed across annotations, estimated with stratified LD score regression (S-LDSC; Finucane et al. 2015).
LDSC3 runs it in two steps: `ldsc ldscore` computes LD scores for baseline annotations and any query annotations, and `ldsc partitioned-h2` fits the regression.
The analysis has two forms, chosen automatically from the LD-score directory:

- **Baseline only.** With no query annotations, one joint model fits all baseline annotations and reports one row per baseline annotation; read `enrichment`.
- **Baseline plus queries.** With query annotations, each query is fitted separately in its own baseline-plus-one-query model and reported as one row; read `coefficient`.

Testing one query annotation per cell type is [cell-type enrichment](cell-type-enrichment.md), which uses the same commands.

## Inputs you need

- **Munged summary statistics** for one trait: the Parquet file written by [`munge-sumstats`](munge-sumstats.md), or a legacy `.sumstats(.gz)` file with `SNP`, `A1`, `A2`, `Z` and `N` (see [legacy summary statistics](../../reference/formats/legacy-sumstats.md)).
- **Baseline annotations**: a chromosome suite of `.annot(.gz)` files, for example a baseline-LD release.
  No baseline annotations are packaged with LDSC3; see [reference data](../reference-data.md).
- **A reference panel**: an [R² matrix](r2-matrices.md) built by `build-r2-matrix` (`--r2-matrix-dir`), or PLINK files (`--plink-prefix`).
  An R² matrix needs its per-chromosome sidecars, and other R² formats are not accepted; see [R² matrix format](../../reference/formats/r2-matrices.md).
- **Annotation files that align**: when several files supply annotation columns for one chromosome, their SNP rows must match exactly and their column names must be unique; see [input paths](../../reference/input-paths.md).
- **Query annotations** (optional), from exactly one of three routes: BED files (`--query-annot-bed-sources`), prebuilt annotation files (`--query-annot-sources`), or gene lists (`--query-annot-gene-list-sources`; see [gene-list LD scores](gene-list-ld-scores.md)).

Use the same SNP identifier mode and genome build for the annotations, the reference panel and the munged summary statistics.
`partitioned-h2` has no build override: if its grouped build check fails, rebuild the mismatched upstream input.
See [SNP identity and genome build](../concepts/snp-identity-and-genome-build.md) and the [genome-build policy](../../reference/snp-identity-and-genome-build.md#genome-build-policy).

## Step 1: compute baseline-plus-query LD scores

This example projects BED files directly onto the baseline annotation rows, with no intermediate query `.annot.gz` files:

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-bed-sources "/path/to/beds/*.bed" \
  --r2-matrix-dir /path/to/r2-matrix/hg19 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --output-dir /path/to/results/partitioned_ldscores
```

- `--baseline-annot-sources` uses `@` as the chromosome placeholder; query annotations always require explicit baseline annotations.
- `--query-annot-bed-sources` takes a quoted glob; each BED file becomes one query annotation.
  BED intervals are zero-based half-open; see [BED intervals](../../reference/formats/bed-intervals.md) and the optional `--padding-bp`.
- `--r2-matrix-dir` points to one build-specific R² matrix.
  The LD window must use the unit the R² matrix was built with and be no wider.
- `--snp-identifier` and `--genome-build` must match the munged summary statistics; coordinate-family modes require `--genome-build` on the command line.
- `--ld-wind-cm 1.0` sets a 1 cM LD window; supply exactly one of `--ld-wind-cm`, `--ld-wind-kb` or `--ld-wind-snps`.

Defaults used by this command:

- Regression SNPs are the packaged HM3 map; `--regr-snps-file` replaces it.
- `--regr-snps-exclude-regions mhc-and-centromeres` removes MHC and centromere SNPs from the regression SNPs only; reference SNPs and annotation counts are unchanged.
- `--common-maf-min 0.05` defines the common-SNP counts (MAF at least 0.05).
- All reference SNPs contribute to LD scores; `--ref-panel-snps-file` restricts them.

For a baseline-only analysis, omit the query flag.
For the difference between reference SNPs and regression SNPs, see [reference and regression SNPs](../../reference/ldscore-snp-universes.md); all `ldscore` options are described in [LD scores](ld-scores.md).

### Reuse query annotations across runs

To keep query annotations as files, project them once with `ldsc annotate` on the baseline annotation grid:

```bash
ldsc annotate \
  --ref-snp-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-bed-sources "/path/to/beds/*.bed" \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --output-dir /path/to/annotations/query_from_beds
```

The annotation grid route writes aligned `baseline.<chrom>.annot.gz` and `query.<chrom>.annot.gz` files on the same cleaned rows.
Pass both families together to `ldscore`:

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/annotations/query_from_beds/baseline.@.annot.gz" \
  --query-annot-sources "/path/to/annotations/query_from_beds/query.@.annot.gz" \
  --r2-matrix-dir /path/to/r2-matrix/hg19 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --output-dir /path/to/results/partitioned_ldscores
```

See [annotations](annotations.md) for the BIM-grid alternative and the [`annotate` reference](../../reference/annotate.md) for its outputs.

## Step 2: fit partitioned heritability

```bash
ldsc partitioned-h2 \
  --sumstats-file /path/to/results/trait/trait.parquet \
  --ldscore-dir /path/to/results/partitioned_ldscores \
  --output-dir /path/to/results/partitioned_h2
```

The command reads the baseline, query, count and overlap files recorded in the LD-score directory's `metadata.json`; it never recomputes LD scores.
An LD-score directory converted from legacy LDSC with [`convert-ldsc2-ldscores`](migrating-from-legacy-ldsc.md) is accepted for a baseline-only analysis.

Defaults used by this command:

- `--count-kind common` scales estimates by common-SNP counts and `--count-kind all` by all-SNP counts; the run stops before fitting if the LD-score directory lacks the requested counts for a baseline annotation, or for a query annotation unless `--continue-on-query-error` records that query as failed; the other kind is never substituted.
- The intercept is estimated; add `--intercept-h2 1` to fix it at 1 (the replacement for legacy LDSC `--no-intercept`; see the [flag map](../../reference/legacy-flag-map.md#regression-intercept-consolidation)).
- `--n-blocks 200` sets the block jackknife.
- Results are on the observed scale; for a binary trait, add `--samp-prev` and `--pop-prev` as described in [liability scale](liability-scale.md).
- `--summary-sort-by auto` sorts query rows by ascending `coefficient_p` and baseline-only rows in annotation order.

For unpartitioned SNP heritability, use [`h2`](heritability.md) instead.

## Reading the outputs

A run with query annotations writes:

```text
partitioned_h2/
    partitioned_h2.tsv
    diagnostics/
        metadata.json
        partitioned-h2.log
        query_status.tsv
        query_annotations/
            manifest.tsv
            0001_<query-name>/
                partitioned_h2.tsv
                partitioned_h2_full.tsv
                coefficient_delete_values.parquet
                metadata.json
```

A baseline-only run writes no per-query tree and keeps its complete-model files, including `diagnostics/coefficient_delete_values.parquet`, under the result root.

`partitioned_h2.tsv` has one row per baseline annotation (baseline only) or per successfully fitted query (baseline plus queries); baseline rows are not repeated in the query summary.
Choose the column by question:

- **Is a baseline annotation enriched?** Read `enrichment` (proportion of heritability divided by proportion of SNPs) and its two-sided `enrichment_p`.
- **Does a query annotation add heritability beyond the baseline?** Read the conditional `coefficient` and its one-sided `coefficient_p`.
  A query's `enrichment` is confounded by its overlap with enriched baseline annotations, so prefer `coefficient_p`.

For continuous annotations, the enrichment-type summaries lack their binary interpretation; see [continuous annotations](continuous-annotations.md).
Every column, including the difference between `prop_h2` and `category_h2_obs` under overlap, is defined in the [result reference](../../reference/results/partitioned-h2.md#partitioned_h2tsv).

The root table's `query_result_dir` column gives each query's folder relative to the result directory, so you can open a query's details even after re-sorting the table; `category` keeps the original query name.
Folder numbers follow this run's sorted order and are not stable identifiers across runs.
Each folder's `partitioned_h2_full.tsv` reports the complete fitted baseline-plus-query model.

Keep the whole result directory together when moving or sharing fitted queries.
Run-wide provenance is stored once in `diagnostics/metadata.json`, and each query's `metadata.json` points to it, so a query folder alone is incomplete for [`quantile-h2`](continuous-annotations.md).
To summarize the result as a figure, run `ldsc plot --result-dir /path/to/results/partitioned_h2`; see [plotting](plotting.md).
The full output layout is listed in [artifact metadata](../../reference/formats/artifact-metadata.md#partitioned-h2).

## Options and variations

### Many query annotations

Partitioned heritability can test many query annotations in one run, for example 1,000 gene lists; each query is still fitted separately against the same baseline.
Both commands take `--query-batch-size`, default 100, and the two settings are independent:

- In `ldscore`, it sets how many query annotations are computed together.
  Each completed batch is written as its own genome-wide query file (`ldscore.query.parquet` for one batch, numbered `ldscore.query.batchNNNNN.parquet` files otherwise), recorded in the `query_batches` manifest of `metadata.json`.
  Direct computation repeats the reference-panel pass for every batch, so smaller batches use less memory but take longer.
- In `partitioned-h2`, it sets how many query annotations are loaded at once; one read can span several query files.
  A smaller value loads fewer query values but does not reduce the SNPs or baseline annotations in each fit.

LD-score directories without the `query_batches` manifest come from an older version; regenerate them.

### Fit queries concurrently

With several allocated CPUs, `--workers` fits several complete baseline-plus-query models at once:

```bash
ldsc partitioned-h2 \
  --sumstats-file /path/to/results/trait/trait.parquet \
  --ldscore-dir /path/to/results/partitioned_ldscores \
  --workers 4 \
  --query-batch-size 100 \
  --output-dir /path/to/results/partitioned_h2
```

This loads at most 100 query annotations at a time and runs up to four fits together.
Each worker performs its own filtering, weighting, regression and block jackknife; SNP order, jackknife blocks and the final row order are the same as in a serial run.
Start with the default of one worker and increase it only within your CPU and memory allocation, because every worker needs memory for its own model.
Worker counts, CPU detection and process start-up are described in [parallel execution](../../reference/parallel-execution.md) and [query workers and memory](../../reference/regression-options.md#43-query-workers-and-memory).

### When a query fails

By default, any query error stops the run without publishing new results.
For a large scan in which the successful fits should be kept, add `--continue-on-query-error`:

| Regression command | If a query fails |
| --- | --- |
| Without `--continue-on-query-error` (default) | Collect the per-query errors and fail without publishing new scientific results. |
| With `--continue-on-query-error` | Skip and record the failed query, finish the scan and publish the successful fits. |

Every completed scan writes `diagnostics/query_status.tsv` in the requested order, with columns `query_annotation`, `status`, `stage`, `error_type` and `error_message`.
Status is `success`, `unestimable` (a block jackknife deletion left the model unidentifiable) or `failed`.
A failed query has no row in `partitioned_h2.tsv` and no folder; its absence is not a null association, so keep the full status table when accounting for the number of tests.
`diagnostics/partitioned-h2.log` records each failure's stage, traceback and, for jackknife failures, the affected blocks and SNP support.
Shared input or output failures, interrupts, worker exits and scans with no successful query still stop the run.
The flag does not change filtering, weights, jackknife blocks or the solver, and it is separate from `ldscore --allow-unresolved-genes`, which acts on gene lists before regression.
See [query failures](../../reference/results/partitioned-h2.md#query-failures-and-continuation) and [troubleshooting](../troubleshooting.md#partitioned-h2-a-query-regression-fails).

### Other regression settings

`--count-kind`, `--intercept-h2`, `--chisq-max`, `--n-blocks`, `--summary-sort-by` and `--allow-identity-downgrade` are described in [regression options](../../reference/regression-options.md) and `ldsc partitioned-h2 --help`.

## Caveats

- Query annotations require explicit baseline annotations.
  The synthetic all-ones `base` annotation is produced only for ordinary unpartitioned LD scores, which serve [`h2`](heritability.md) and [`rg`](genetic-correlation.md), not partitioned heritability.
- `partitioned-h2` needs `ldscore.overlap.parquet`, which `ldscore` writes when there are at least two annotation columns; see [missing overlap matrix](../troubleshooting.md#partitioned-h2-missing-overlap-matrix).
- Existing outputs, including a stale `diagnostics/query_annotations/` tree, stop the run before anything is written; add `--overwrite` only when you intend to replace them.
  See [logs and outputs](../../reference/logs-and-outputs.md).
- Query runs always write the per-query tree; the retired `--write-per-query-results` flag is rejected.
- `prop_h2`, `category_h2_obs` and `enrichment` are regression estimates and can be negative or exceed one when estimates are noisy or annotations overlap.

## Reference

- [Partitioned heritability result columns](../../reference/results/partitioned-h2.md)
- [Regression options](../../reference/regression-options.md)
- [Reference and regression SNPs](../../reference/ldscore-snp-universes.md)
- [Artifact metadata](../../reference/formats/artifact-metadata.md)
- [Input paths](../../reference/input-paths.md)
- [Logs and outputs](../../reference/logs-and-outputs.md)
- [Parallel execution](../../reference/parallel-execution.md)
- [Guided tutorial](../guided-tutorial.md)

For the Python API, see [Python usage](../python-usage.md#partitioned-heritability).
