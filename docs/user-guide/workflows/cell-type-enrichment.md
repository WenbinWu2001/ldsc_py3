# Cell-type enrichment

Last updated on: 2026-10-28

Cell-type enrichment asks which cell types' annotations carry more trait heritability than the baseline annotations explain.
Each cell type becomes one query annotation, for example the SNPs near its specifically expressed genes or inside its open-chromatin peaks.
LDSC3 fits every query separately in a model of the form

`[query annotation, baseline annotations]`

with stratified LD score regression (S-LDSC; Finucane et al. 2015), as LDSC-SEG (Finucane et al. 2018) did for specifically expressed genes.
The query's conditional `coefficient` and its one-sided `coefficient_p` measure whether the cell type adds heritability beyond the baseline.
This is the query form of [partitioned heritability](partitioned-heritability.md), run with the same `ldscore` and `partitioned-h2` commands; that page explains each step, the outputs and the batch and worker options in full.

## Inputs you need

- **Munged summary statistics** for the trait: a Parquet file from [`munge-sumstats`](munge-sumstats.md), or a legacy `.sumstats(.gz)` file.
- **Baseline annotations**, kept in every model.
- **One query per cell type**, from one of three routes:
  - one BED file per cell type (Route A below);
  - one gene list per cell type, such as its top specifically expressed genes (Route B below);
  - prebuilt annotation files with one column per cell type (`--query-annot-sources`; see [partitioned heritability](partitioned-heritability.md#reuse-query-annotations-across-runs)).
- **A reference panel** ([R² matrix](r2-matrices.md) or PLINK files), or, for gene lists only, a prebuilt [gene LD-score index](gene-list-ld-scores.md).

Use the same SNP identifier mode and genome build throughout; see [SNP identity and genome build](../concepts/snp-identity-and-genome-build.md).

## Route A: one BED file per cell type

Compute baseline-plus-cell-type LD scores, projecting the BED files directly:

```bash
ldsc ldscore \
  --baseline-annot-sources "/path/to/baseline/baseline.@.annot.gz" \
  --query-annot-bed-sources "/path/to/cell_type_beds/*.bed" \
  --r2-matrix-dir /path/to/r2-matrix/hg19 \
  --snp-identifier chr_pos_allele_aware \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --output-dir /path/to/results/cell_type_ldscores
```

Then fit one model per cell type:

```bash
ldsc partitioned-h2 \
  --sumstats-file /path/to/results/trait/trait.parquet \
  --ldscore-dir /path/to/results/cell_type_ldscores \
  --output-dir /path/to/results/cell_type_h2
```

Each BED file becomes one query annotation named after its file name without the final suffix, so the names must be unique; see [BED intervals](../../reference/formats/bed-intervals.md).
To reuse the cell-type annotations across traits or LD windows, write them once with `ldsc annotate`; see [reuse query annotations](partitioned-heritability.md#reuse-query-annotations-across-runs).

## Route B: one gene list per cell type

This case study tests protein-coding gene lists, one list per cell type.
It omits the optional control gene list, so each query is conditioned on the baseline annotations only.
Finucane et al. (2018) also conditioned on an all-genes annotation; to add a fixed gene control, see [optional control gene list](gene-list-ld-scores.md#optional-control-gene-list).

### Set up inputs

```bash
GENE_LIST_SOURCES="/path/to/query_gene_lists/*.txt"
TRAIT_NAME="trait_name"
SUMSTATS_FILE="/path/to/munged/${TRAIT_NAME}.parquet"
RESULT_ROOT="/path/to/results/cell_type_enrichment"
```

Each gene-list file holds one Ensembl gene ID or gene name per line, without a header, and each matched file is tested as a separate query.
The file format and name rules are in [gene-list LD scores](gene-list-ld-scores.md) and [gene lists](../../reference/formats/gene-lists.md).

### Step 1: compute LD scores

Choose direct mode or indexed mode.
Both write the same kind of LD-score directory for Step 2.

#### Option A: direct mode

Direct mode reads the baseline annotations and the PLINK reference panel and projects the gene intervals chromosome by chromosome:

```bash
BASELINE_ANNOT_SOURCES="/path/to/resources/baseline/baseline.@.annot.gz"
PLINK_PREFIX="/path/to/resources/plink/1000G.EUR.QC."
GENE_COORDINATE_FILE="/path/to/resources/gene-coordinates.hg19.tsv.gz"
LDSCORE_DIR="${RESULT_ROOT}/direct_mode/ldscore"

ldsc ldscore \
  --query-annot-gene-list-sources "${GENE_LIST_SOURCES}" \
  --baseline-annot-sources "${BASELINE_ANNOT_SOURCES}" \
  --plink-prefix "${PLINK_PREFIX}" \
  --gene-coordinate-file "${GENE_COORDINATE_FILE}" \
  --snp-identifier rsid \
  --genome-build hg19 \
  --ld-wind-cm 1.0 \
  --padding-bp 100000 \
  --gene-exclude-regions mhc \
  --output-dir "${LDSCORE_DIR}"
```

- `--plink-prefix` takes the shared stem of the chromosome files, here `1000G.EUR.QC.` for `1000G.EUR.QC.1.bed`, `.bim` and `.fam`; an `@` placeholder also works.
  The PLINK files must correspond to the baseline annotation suite.
- `--gene-coordinate-file` is the only source of gene coordinates: one-based, inclusive, and in the build given by `--genome-build` (hg19 or hg38); nothing is packaged.
- `--genome-build hg19` sets the build for gene projection and the named regions.
- `--padding-bp 100000` adds 100 kb on each side of every gene; direct gene-list runs require an explicit value, `0` for gene bodies.
- `--gene-exclude-regions mhc` drops genes whose unpadded interval overlaps the MHC, before padding; it also applies to a control gene list.
  It removes genes, not reference SNPs; the direct-mode default is `none`.
- Resolution is strict: any unresolved or ambiguous identifier stops the run before LD-score computation.
  Add `--allow-unresolved-genes` only for a deliberate exploratory subset analysis.
- Add `--overwrite` only when you intend to replace earlier results in `--output-dir`.

Defaults used by this command:

- `--regr-snps-exclude-regions mhc-and-centromeres` removes MHC and centromere SNPs from the regression SNPs; it changes neither the reference SNPs nor the annotation counts and overlaps.
- `--common-maf-min 0.05` defines common SNPs (MAF at least 0.05) for the common-SNP counts.

Direct mode writes private query files under the output directory while it computes, but does not publish reusable `.annot.gz` files; allow temporary disk space.
A completed or handled failed run removes them; an interrupted run can leave them behind.
Start with `--workers 1`, and lower `--query-batch-size` from its default of 100 to reduce the active query workspace; direct mode repeats the reference-panel pass for every query batch.

#### Option B: indexed mode

Indexed mode assembles the LD scores, annotation counts and overlaps from a prebuilt gene LD-score index instead of rereading the PLINK genotypes:

```bash
GENE_LDSCORE_INDEX_DIR="/path/to/resources/gene-ldscore-index"
LDSCORE_DIR="${RESULT_ROOT}/indexed_mode/ldscore"

ldsc ldscore \
  --query-annot-gene-list-sources "${GENE_LIST_SOURCES}" \
  --gene-ldscore-index-dir "${GENE_LDSCORE_INDEX_DIR}" \
  --output-dir "${LDSCORE_DIR}"
```

The index fixes the baseline annotations, reference SNPs, SNP identifier mode, genome build (hg19 or hg38), LD window, gene padding, gene exclusion, regression SNPs, region exclusion, MAF rules, annotation counts and regression weights; its metadata is the record of these settings.
Do not pass live settings such as `--baseline-annot-sources`, `--plink-prefix`, `--r2-matrix-dir`, `--snp-identifier`, `--genome-build`, `--ld-wind-*`, `--padding-bp`, `--gene-exclude-regions` or `--regr-snps-file`.
`--workers` and `--query-batch-size` are accepted; a smaller batch reduces temporary space but rereads the index operators for each batch.
An optional `--control-gene-list-file` must be one existing file and inherits the index's gene coordinate catalog, build, padding and gene exclusion.
The finished LD-score directory is self-contained; Step 2 does not need the index.
To build an index, see [gene-list LD scores](gene-list-ld-scores.md) and the [gene LD-score index reference](../../reference/gene-ldscore-index.md).

#### Check the gene-list diagnostics

Both modes write gene-list diagnostics under `diagnostics/` in the LD-score directory.
Start with `gene_list_resolution_summary.tsv`, then filter `gene_list_audit.tsv.gz`; `query_annotation_status.tsv` explains queries skipped after identifier checks.
Without a control gene list, `ldscore.baseline.parquet` contains only the supplied baseline annotations and no `gene_control` column.
See [gene-list diagnostics](../../reference/results/gene-list-diagnostics.md) for the repair workflow.

### Step 2: fit one model per gene list

```bash
PARTITIONED_H2_DIR="${RESULT_ROOT}/partitioned_h2/${TRAIT_NAME}"

ldsc partitioned-h2 \
  --sumstats-file "${SUMSTATS_FILE}" \
  --ldscore-dir "${LDSCORE_DIR}" \
  --output-dir "${PARTITIONED_H2_DIR}"
```

Defaults used by this command:

- `--count-kind common` uses annotation counts from common reference SNPs.
- `--num-blocks 200` sets the block jackknife.
- The intercept is estimated; add `--intercept-h2 1` to fix it at 1.
- Without `--samp-prev` and `--pop-prev`, results are on the observed scale; see [liability scale](liability-scale.md).
- `--summary-sort-by auto` sorts the gene lists by ascending `coefficient_p`.
  The other choices are `category`, `prop-snps`, `prop-h2`, `enrichment`, `enrichment-p`, `coefficient` and `coefficient-p`; p-values sort ascending, other numbers descending, and missing values last.

By default, any failed query stops the run without publishing new results; add `--continue-on-query-error` to keep the successful fits.
For many gene lists, `--workers 4 --query-batch-size 100` fits up to four models at once while loading at most 100 queries.
Both options are explained in [partitioned heritability](partitioned-heritability.md#options-and-variations).

## Reading the results

The root `partitioned_h2.tsv` has one row per successfully fitted query and no baseline rows.
Focus on `coefficient` and its one-sided `coefficient_p`: a positive coefficient with a small p-value indicates that the cell type contributes heritability beyond the baseline annotations.
A query's `enrichment` also reflects its overlap with enriched baseline annotations, so it is not the test of interest here.

Use `query_result_dir` to open each query's folder under `diagnostics/query_annotations/`, relative to the result directory; `category` keeps the original query name, and folder numbers follow this run's sorted order.
Each folder's `partitioned_h2_full.tsv` reports the complete fitted model, `[query annotation, baseline annotations]`.
Queries that failed appear only in `diagnostics/query_status.tsv`, with status `unestimable` or `failed`, and have no folder.
Copy the whole result directory when moving fitted queries, because each query's metadata refers to the shared `diagnostics/metadata.json`.
The output layout and columns are described under [reading the outputs](partitioned-heritability.md#reading-the-outputs) and in the [result reference](../../reference/results/partitioned-h2.md).

For an exploratory summary across all queries, run:

```bash
ldsc plot --result-dir "${PARTITIONED_H2_DIR}"
```

The horizontal scatter plot shows the nominal one-sided $-\log_{10}(P)$ from the separate baseline-conditional fits.
It does not compare raw coefficients, show enrichment or correct for multiple testing; see [plotting](plotting.md).

## Caveats

- Each query is a separate test; correct for the number of cell types or gene lists you tested, counting failed queries from `diagnostics/query_status.tsv`.
  A failed query is not a null association.
- Query annotations require explicit baseline annotations; the synthetic all-ones `base` annotation applies only to ordinary unpartitioned LD scores.
- A baseline-only LD-score directory gives the joint baseline model, not cell-type tests.
- For how this workflow differs from legacy LDSC's cell-type analysis, see [migrating from legacy LDSC](migrating-from-legacy-ldsc.md#partitioned-and-cell-type-analyses).

## Reference

- [Partitioned heritability](partitioned-heritability.md)
- [Gene-list LD scores](gene-list-ld-scores.md)
- [Partitioned heritability result columns](../../reference/results/partitioned-h2.md)
- [Gene lists](../../reference/formats/gene-lists.md) and [gene-list diagnostics](../../reference/results/gene-list-diagnostics.md)
- [Gene LD-score index](../../reference/gene-ldscore-index.md)
- [Regression options](../../reference/regression-options.md)
- [Parallel execution](../../reference/parallel-execution.md)

For the Python API, see [Python usage](../python-usage.md#cell-type-enrichment).
