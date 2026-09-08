# LDSC3 - Guided Analysis Tutorial

Last updated on: 2026-09-07

This tutorial walks through how to use the `ldsc` package for a series of LD score-based heritability analyses.

## Overview

The analysis pipeline involves:

- **Before any analysis:** curate a standardized sumstats file with `ldsc munge-sumstats`.
- **Analysis 1:** estimate the heritability of a trait.
  - Step 1: compute (unpartitioned) LD scores with `ldsc ldscore`.
  - Step 2: estimate trait heritability with `ldsc h2`.
- **Analysis 2:** estimate cross-trait genetic correlation between multiple traits.
  - Step 1: compute (unpartitioned) LD scores with `ldsc ldscore` (same as in Analysis 1).
  - Step 2: estimate cross-trait genetic correlation with `ldsc rg`.
- **Analysis 3:** partition heritability with functional annotations (known as cell-type-specific regression).
  - Step 1: compute partitioned LD scores with `ldsc ldscore`, supplying query annotations.
  - Step 2: estimate the heritability contribution of each query annotation with `ldsc partitioned-h2`.

The files you will need are:

- raw GWAS summary statistics files for the traits of interest (or munged sumstats from the legacy ldsc python2 codebase);
- a reference panel for LD score calculation (either R2 parquets or the PLINK suite);
- a set of baseline annotations;
- raw BED files or one-column gene lists for the pathways / cell types whose h2 contribution you want to test (the "query annotations").


As a motivating example, we study the `mdd2025` trait, using 1000 Genomes Phase 3 as the reference panel and `1000G_EUR_Phase3_baseline` as the baseline annotations, with `Hippocampus_PP1.bed`, `Cerebellum_PC16.bed`, `Cerebellum_PP3.bed`, and `Cerebellum_PP1.bed` as the pathways whose heritability contribution we test after controlling for the baseline annotations.

All of these resources can be found under `ldsc3_test_bundle/resources/`, except for the raw GWAS sumstats files.

Before running the commands below, follow the README to clone the `main` branch of the GitHub repo and install it in a suitable conda environment. Then run `conda activate ldsc3` before the analysis.

We first set up the input and output directories:

```bash
INPUT_ROOT="/path/to/ldsc3_test_bundle"
OUTPUT_ROOT="${INPUT_ROOT}/tutorial_output"
```

### Remarks

1. By default, `ldsc ldscore` writes regression rows from the bundled HapMap3 set, after applying its default MHC-and-centromere exclusion. `h2`, `rg`, and `partitioned-h2` consume those rows and have no HapMap3 flag. The separate `munge-sumstats --use-hm3-snps` flag optionally restricts summary statistics to the bundled HapMap3 map.
2. For flags ending in `-sources`, you can use the glob pattern `*` to match multiple files and `@` as the placeholder for the chromosome number.
3. Output directories are created automatically. Use `--overwrite` to allow overwriting existing output files.
4. `ldsc ldscore` excludes the extended MHC and centromere regions by default.
5. Be careful about the genome build of your input file. Although the package runs a guardrail check on the genome build, it is not comprehensive.
6. To inspect what happened under the hood, check the log file in the `diagnostics/` directory inside each command's output directory.
7. The memory and run-time recommendations below are NOT accurate (codebase not yet benchmarked); treat them as rough upper bounds. If you hit OOM, increase the memory allocation. 24 GB should be sufficient for the initial test runs.

## Config you should follow across all steps

- `snp_identifier=chr_pos`: use `chr` + `pos` as the unique identifier for a SNP.
- ... [TODO]

## Munge-sumstats

This section illustrates the `ldsc munge-sumstats` command.

**Goal:** convert raw summary statistics files into a standardized format for use in downstream ldsc regression.

**Recommended memory allocation:** < 4 GB

**Expected running time:** < 5 min

**Minimal command:**

```bash
TRAIT_NAME="mdd2025"  # used as the output dir name and as a trait label in downstream regression outputs
RAW_SUMSTATS_FILE="/path/to/mdd2025_raw_sumstats.tsv"
SUMSTATS_OUT_DIR="${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAME}"

ldsc munge-sumstats \
  --snp-identifier chr_pos \
  --source-genome-build "auto" \
  --output-genome-build "hg19" \
  --raw-sumstats-file "${RAW_SUMSTATS_FILE}" \
  --use-hm3-snps \
  --trait-name "${TRAIT_NAME}" \
  --output-dir "${SUMSTATS_OUT_DIR}" \
  --overwrite
```

1. The package auto-infers the genome build, the raw sumstats file format, and common column names, and automatically strips any comment lines starting with `#` from the file header. You can still pass legacy flags such as `--snp` to specify a column name if your sumstats file uses non-standard column names. For the table of column aliases the package recognizes, see the detailed `munge-sumstats` wiki.
2. Before running the command for real, we strongly recommend using `--infer-only` to preview the inferred columns and suggested flags without writing any output. This is especially helpful if you are unsure which flags to use or about the genome build of your input file: the program will suggest the correct command to run. See the `munge-sumstats` wiki for details.
   ```bash
   ldsc munge-sumstats \
   ... (your-guessed-flags) \
   --infer-only
   ```
3. The legacy `.sumstats.gz` format is supported. Specify the file via `--raw-sumstats-file` and pass a `--trait-name` so that downstream outputs use stable trait labels.
4. If the source and output builds differ, add either `--use-hm3-quick-liftover` together with `--use-hm3-snps`, or provide `--liftover-chain-file`.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
sumstats_processed/
    mdd2025/
        sumstats.parquet
        diagnostics/
            sumstats.log
            dropped_snps/
                dropped.tsv.gz
```

## Specify the reference panel

Before the regressions in any of the analyses below, we first specify the path to the reference panel:

```bash
# LD reference panel from 1000 Genomes Phase 3. Use either R2 parquet files or PLINK files.
R2_DIR="${INPUT_ROOT}/resources/ldref/ldref_r2_1kgph3_1cM"
PLINK_PREFIX="${INPUT_ROOT}/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC."  # hg19
```

Two reference-panel sources are supported (choose one):

- per-chromosome R2 parquets computed from the PLINK files;
- the PLINK suite itself.

## Analysis 1: estimate the heritability of a trait

**Goal:** compute LD scores from a genomic reference panel, then estimate the trait's heritability using the munged sumstats file.

**Recommended memory allocation:** Step 1: 24 GB (generous, for safety); Step 2: < 4 GB

**Expected running time:** Step 1: < 1 h; Step 2: < 1 min

### Step 1: compute (unpartitioned) LD scores with `ldsc ldscore`

```bash
# Step 1: Calculate unpartitioned LD scores from R2 parquets (via `--r2-dir`).
LDSCORE_OUTPUT_DIR="${OUTPUT_ROOT}/unpartitioned_ldscore"

ldsc ldscore \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --r2-dir "${R2_DIR}/hg19" \
  --output-dir "${LDSCORE_OUTPUT_DIR}" \
  --ld-wind-cm 1.0 \
  --overwrite
```

If you use the PLINK suite as the reference panel, replace `--r2-dir ...` in the command with:

```bash
--plink-prefix "${PLINK_PREFIX}"
```

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
unpartitioned_ldscore/
    metadata.json
    ldscore.baseline.parquet
    diagnostics/
        ldscore.log
```

### Step 2: estimate trait heritability with `ldsc h2`

```bash
# Step 2: Regress sumstats on unpartitioned LD scores.
SUMSTATS_FILE="${SUMSTATS_OUT_DIR}/sumstats.parquet"
LDSCORE_DIR="${LDSCORE_OUTPUT_DIR}"  # the dir storing unpartitioned ldscore.baseline.parquet and metadata.json

H2_OUTPUT_DIR="${OUTPUT_ROOT}/h2/${TRAIT_NAME}"

ldsc h2 \
  --sumstats-file "${SUMSTATS_FILE}" \
  --ldscore-dir "${LDSCORE_DIR}" \
  --output-dir "${H2_OUTPUT_DIR}" \
  --overwrite
```

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
h2/mdd2025/
    h2.tsv
    diagnostics/
        ld_score_regression_bins.tsv
        metadata.json
        h2.log
```

The bin table records the exact fitted-data summary used by the optional binned LD Score regression diagnostic. Install the plotting extra and run `ldsc plot --result-dir "${H2_OUTPUT_DIR}"` to create it.

## Analysis 2: estimate cross-trait genetic correlation between multiple traits

**Goal:** estimate the genetic correlation between traits using the GWAS sumstats and LD scores.

**Recommended memory allocation:** Step 1: reuses Analysis 1 LD scores; Step 2: < 4 GB

**Expected running time:** Step 1: reuses Analysis 1 LD scores; Step 2: < 1 min

### Munge two additional sumstats

```bash
# Munge two other traits.
TRAIT_NAMES=("scz2022" "adhd2019")
RAW_SUMSTATS_FILES=(
  "/path/to/scz2022_raw_sumstats.tsv"
  "/path/to/adhd2019_raw_sumstats.tsv"
)  # change raw sumstats file path here

for i in "${!TRAIT_NAMES[@]}"; do
  TRAIT_NAME="${TRAIT_NAMES[$i]}"
  RAW_SUMSTATS_FILE="${RAW_SUMSTATS_FILES[$i]}"
  SUMSTATS_OUT_DIR="${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAME}"

  ldsc munge-sumstats \
    --snp-identifier chr_pos \
    --source-genome-build "auto" \
    --output-genome-build "hg19" \
    --raw-sumstats-file "${RAW_SUMSTATS_FILE}" \
    --use-hm3-snps \
    --trait-name "${TRAIT_NAME}" \
    --output-dir "${SUMSTATS_OUT_DIR}" \
    --overwrite
done

```

### Step 1: compute (unpartitioned) LD scores with `ldsc ldscore` (same as in Analysis 1)

Reuse the unpartitioned LD scores from the h2 section.

### Step 2: estimate cross-trait genetic correlation with `ldsc rg`

This command computes the cross-trait genetic correlation for *all* pairs of the input sumstats.

```bash
# Step 2: Regress all trait pairs on unpartitioned LD scores.
LDSCORE_DIR="${OUTPUT_ROOT}/unpartitioned_ldscore"
SUMSTATS_SOURCES="${OUTPUT_ROOT}/sumstats_processed/*/sumstats.parquet"  # use glob pattern to match all sumstats

RG_OUTPUT_DIR="${OUTPUT_ROOT}/rg/mdd2025_scz2022_adhd2019"

ldsc rg \
  --sumstats-sources "${SUMSTATS_SOURCES}" \
  --ldscore-dir "${LDSCORE_DIR}" \
  --output-dir "${RG_OUTPUT_DIR}" \
  --overwrite
```

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
rg/mdd2025_scz2022_adhd2019/
    rg.tsv
    rg_full.tsv
    h2_per_trait.tsv
    diagnostics/
        metadata.json
        rg.log
```

## Analysis 3: partition heritability with functional annotations (known as cell-type-specific regression)

**Goal:** partition a trait's heritability into the contributions of a set of pathways / cell types, testing whether each query annotation contributes to heritability after controlling for the baseline annotations.

Many hypotheses in this analysis are gene sets: genes differentially expressed
in a tissue or cell type, genes prioritized by proteomics, or genes sharing a
GO or SynGO term. LDSC-SEG demonstrated that S-LDSC can test whether
heritability is enriched near specifically expressed genes conditional on both
the baseline model and an all-genes annotation
([Finucane et al., 2018](https://doi.org/10.1038/s41588-018-0081-4)). For a
large collection of gene sets, an exact gene LD-score index computes the fixed
PLINK/reference work once and reuses it without changing the downstream model.
The dedicated [LDSC-SEG tutorial for protein-coding gene
lists](LDSC-SEG-PC-genes.md) does not add a control-gene annotation and
fits each query with the baseline annotations only.

**Recommended memory allocation:** Step 1: 24 GB (generous, for safety); Step 2: < 4 GB

**Expected running time:** Step 1: < 2 h; Step 2: < 5 min (depends on the number of query annotations)

### Step 1: compute partitioned LD scores with `ldsc ldscore` using query annotations

We also need to specify the paths to the raw BED files (here, four of them) and to the baseline annotations:

```bash
# annotations
RAW_QUERY_BED_SOURCES="$INPUT_ROOT/resources/example_raw_annot/*.bed"  # put all your query BED files in one directory
BASELINE_ANNOT_SOURCES="${INPUT_ROOT}/resources/1000G_EUR_Phase3_baseline/baseline.@.annot.gz"
```

Place all of your raw query files in the same directory and use a glob pattern (as above) to match them. This lets their LD scores be computed in parallel within `ldsc ldscore`.

```bash
PARTITIONED_LDSCORE_OUTPUT_DIR="${OUTPUT_ROOT}/partitioned_ldscore"

# Step 1: Calculate partitioned LD scores. These can be reused across traits.
ldsc ldscore \
  --snp-identifier chr_pos \
  --genome-build "hg19" \
  --query-annot-bed-sources "${RAW_QUERY_BED_SOURCES}" \
  --baseline-annot-sources "${BASELINE_ANNOT_SOURCES}" \
  --r2-dir "${R2_DIR}/hg19" \
  --output-dir "${PARTITIONED_LDSCORE_OUTPUT_DIR}" \
  --ld-wind-cm 1.0 \
  --overwrite
```

**Remarks:**

- For PLINK input, replace `--r2-dir ...` with `--plink-prefix "${PLINK_PREFIX}"`.
- For partitioned h2 analysis, `ldsc ldscore` accepts BED files directly, so you do not need to build the annotations yourself to run this analysis. If you do want to generate annotations for other purposes, follow the *Make annotations* section below.
- Alternatively, pass one-column plain/gzip lists through
  `--query-annot-gene-list-sources`; direct mode also requires a one-based
  build-aware `--gene-coordinate-file` and an explicit `--padding-bp`. Exact
  IDs and case-sensitive names resolve only against that catalog. The BED and
  gene-list flags are mutually exclusive.
- If any input query gene list or BED file is not found in the scientific
  results, that query hit a failure. Check
  `diagnostics/query_annotation_status.tsv` for the reason. For gene lists,
  start with `diagnostics/gene_list_resolution_summary.tsv` and then filter
  `diagnostics/gene_list_audit.tsv.gz`. See the detailed
  [diagnostics and repair guide](../current/gene-list-diagnostics-and-repair.md).
- For a prebuilt exact gene index, use
  `--gene-ldscore-index-dir <index-dir>` with
  `--query-annot-gene-list-sources` and omit live baseline, PLINK/R²,
  genome-build, window, padding, and region arguments. Remove `--padding-bp`
  entirely rather than passing zero. The indexed run writes
  the same canonical output directory and adds no gene control by default. To
  add one, pass an existing one-column file with `--control-gene-list-file`.
- For the complete indexed workflow, see [Build an exact gene LD-score
  index](https://github.com/WenbinWu2001/ldsc_py3/blob/ldsc3-beta/docs/wiki/utility-functionalities/build-gene-ldscore-index.md) and [Calculate LD
  scores for gene lists with an index](https://github.com/WenbinWu2001/ldsc_py3/blob/ldsc3-beta/docs/wiki/main-functionalities/ldscore-from-gene-list.md).
  Per-chromosome `Finished` lines report durable private staging; the index
  becomes public only after complete reload validation. Stages cannot be
  resumed or used for incremental chromosome updates.

**Caveat:** if you use `--ld-wind-cm`, make sure your PLINK suite has non-missing genetic coordinates (the third column in the `.bim` file). If they are missing (e.g., all zeros), the program will raise an error.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
partitioned_ldscore/
    metadata.json
    ldscore.baseline.parquet
    ldscore.query.parquet
    ldscore.overlap.parquet
    diagnostics/
        ldscore.log
```

### Step 2: estimate each query annotation's heritability contribution with `ldsc partitioned-h2`

```bash
# Step 2: Regress sumstats on partitioned LD scores.
SUMSTATS_FILE="${SUMSTATS_OUT_DIR}/sumstats.parquet"
PARTITIONED_LDSCORE_DIR="${PARTITIONED_LDSCORE_OUTPUT_DIR}"

PARTITIONED_H2_OUTPUT_DIR="${OUTPUT_ROOT}/partitioned-h2/${TRAIT_NAME}"

ldsc partitioned-h2 \
  --sumstats-file "${SUMSTATS_FILE}" \
  --ldscore-dir "${PARTITIONED_LDSCORE_DIR}" \
  --output-dir "${PARTITIONED_H2_OUTPUT_DIR}" \
  --overwrite
```

**Remarks:**

1. `partitioned-h2` reports the total h2 implied by the partitioned model; this is not necessarily identical to the standalone unpartitioned h2 estimate above.
2. Query-annotation runs automatically save each query's complete baseline-plus-query fit under `diagnostics/query_annotations/`; the deprecated `--write-per-query-results` flag is unnecessary.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
partitioned-h2/mdd2025/
    partitioned_h2.tsv
    diagnostics/
        metadata.json
        partitioned-h2.log
        query_annotations/
            manifest.tsv
            0001_<query-annotation-slug>/
                partitioned_h2.tsv
                partitioned_h2_full.tsv
                coefficient_delete_values.parquet
                metadata.json
                ...
```

Here:

- `partitioned_h2.tsv` is a summary over all query annotations, sorted by the p-value of the query annotation's coefficient (`coefficient-p`), from most to least significant.
- `query_annotations/` holds the per-query results, which also include the coefficients and h2 share of the baseline annotations. `manifest.tsv` records each original query annotation name, its results folder name, and the paths to the relevant results files.

For a continuous annotation, continue with [Continuous annotations in partitioned LDSC](continuous-annotation-partitioned-ldsc.md). That post-fit workflow uses one complete fitted model and its saved jackknife coefficient values; it does not combine different per-query regressions.

For concise exploratory figures from h2, rg, functional partitioning, cell-type/query, or continuous-annotation quantile results, follow the [plotting results manual](../../tutorials/plotting-results.md). The cell-type summary plot uses the aggregate partitioned-h2 root because each point is the nominal one-sided p-value from a separate baseline-conditional query fit.

## Make annotations

**Goal:** curate binary annotations from BED files, given a set of baseline annotation files.

**Recommended memory allocation:** < 8 GB for four query BED files

**Expected running time:** < 15 min for four query BED files

This step is not part of the analysis pipeline. We preserve it for convenience in case you need it for other purposes.

```bash
ANNOT_OUT_DIR="${OUTPUT_ROOT}/annot_processed"

ldsc annotate \
    --snp-identifier chr_pos \
    --genome-build hg19 \
    --query-annot-bed-sources "$RAW_QUERY_BED_SOURCES" \
    --baseline-annot-sources "$BASELINE_ANNOT_SOURCES" \
    --output-dir "$ANNOT_OUT_DIR" \
    --overwrite
```

**Remarks:**

- For each chromosome, all input query annotations are combined into one file as separate columns, with each column name being the input's base filename (used as the annotation name).
- The resulting `query.@.annot.gz` suite is the canonical input to downstream `ldsc ldscore`: one query annotation file per chromosome, containing multiple annotation columns; each column is one query annotation. Pass it with `--query-annot-sources "$ANNOT_OUT_DIR/query.@.annot.gz"`.
- Do not create one chromosome-sharded suite per query. In sharded mode, `ldscore` expects one query file per chromosome and reads the separate queries from its columns.
- The CM column is preserved in annotations for backward compatibility but is not read by the package downstream.
- See the `annotate` wiki for more on this functionality.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
annot_processed/
    query.1.annot.gz ... query.22.annot.gz
    diagnostics/
        metadata.json
        annotate.log
        dropped_snps/
            dropped.tsv.gz
```


## Backward compatibility with the legacy ldsc python2 codebase

Munged sumstats, legacy ld ref suite (unpartitioned / baseline ref) -- not tested.

[TODO]

## TODO

- How to reuse previously generated annotations for partitioned LDSC.
- Refine memory and run-time numbers with proper benchmarking rather than guessing from log files. In particular, the SLURM memory figure for `ldscore` is inaccurate (it somehow always reports the allocated memory minus 2 MB).
- complete main functionality wiki. add link in this guided tutorial.
- go with quarto?
- User checklist: snp id, genome build, and etc.
