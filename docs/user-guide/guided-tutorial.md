# Guided tutorial

Last updated on: 2026-11-09

This tutorial runs one complete analysis with LDSC3, from raw summary statistics to SNP heritability, genetic correlation and partitioned heritability.
It follows the main path only; each section links to the workflow page that covers the options, outputs and caveats in full.

## What you will run

| Stage | Command | Result |
| --- | --- | --- |
| Munge summary statistics for three traits | `ldsc munge-sumstats` | `sumstats_processed/<trait>/<trait>.parquet` |
| Analysis 1: SNP heritability | `ldsc ldscore`, then `ldsc h2` | `h2/mdd2025/h2.tsv` |
| Analysis 2: genetic correlation | `ldsc rg`, then `ldsc plot` | `rg/mdd2025_scz2022_adhd2019/rg.tsv` |
| Analysis 3: partitioned heritability (cell-type enrichment) | `ldsc ldscore` with query annotations, then `ldsc partitioned-h2` | `partitioned-h2/mdd2025/partitioned_h2.tsv` |

Analyses 1 and 2 share the same unpartitioned LD scores.
Analysis 3 computes its own LD scores with baseline and query annotations.

## Before you begin

Install the development version as described in [installation](installation.md#development-version), then activate the environment:

```bash
conda activate ldsc3-dev
```

The example studies the trait `mdd2025`, with `scz2022` and `adhd2019` added for genetic correlation.
It uses 1000 Genomes Phase 3 as the reference panel, `1000G_EUR_Phase3_baseline` as the baseline annotations, and four BED files (`Hippocampus_PP1.bed`, `Cerebellum_PC16.bed`, `Cerebellum_PP3.bed` and `Cerebellum_PP1.bed`) as the query annotations.
All of these resources are under `ldsc3_test_bundle/resources/`; the raw summary statistics are not included, so replace each `/path/to/...` below with your own files.

Run every command below in one Bash session, starting with these variables:

```bash
INPUT_ROOT="/path/to/ldsc3_test_bundle"
OUTPUT_ROOT="${INPUT_ROOT}/tutorial_output"

# Reference panel from 1000 Genomes Phase 3: an R² matrix, or the PLINK files it was built from.
R2_MATRIX_DIR="${INPUT_ROOT}/resources/ldref/ldref_r2_1kgph3_1cM"
PLINK_PREFIX="${INPUT_ROOT}/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC."  # hg19
```

The tutorial uses two settings throughout:

- the SNP identifier mode `chr_pos`, which identifies a SNP by chromosome and position without alleles (the package default is `chr_pos_allele_aware`);
- the genome build hg19 for every coordinate input.

Every input you combine must use the same SNP identifier mode and genome build; see [SNP identity and genome build](concepts/snp-identity-and-genome-build.md).

A few conventions apply to every command:

- Quote path tokens such as globs and `@` chromosome suites so the shell passes them to LDSC3 unchanged; see [input paths](../reference/input-paths.md#pattern-support-in-command-help).
- Each command creates its `--output-dir` if needed; `--overwrite` allows it to replace the outputs of an earlier run.
- Each command writes its log under `diagnostics/` in its output directory.
- By default, the regression SNPs are the HM3 map without the MHC and centromere regions, and munging keeps only HM3 SNPs; see [reference and regression SNPs](../reference/ldscore-snp-universes.md).

## Munge summary statistics

**Goal:** convert raw summary statistics into munged summary statistics, the input of every regression command.

Munge all three traits once:

```bash
TRAIT_NAMES=("mdd2025" "scz2022" "adhd2019")
RAW_SUMSTATS_FILES=(
  "/path/to/mdd2025_raw_sumstats.tsv"
  "/path/to/scz2022_raw_sumstats.tsv"
  "/path/to/adhd2019_raw_sumstats.tsv"
)

for i in "${!TRAIT_NAMES[@]}"; do
  ldsc munge-sumstats \
    --snp-identifier chr_pos \
    --source-genome-build "auto" \
    --output-genome-build "hg19" \
    --raw-sumstats-file "${RAW_SUMSTATS_FILES[$i]}" \
    --trait-name "${TRAIT_NAMES[$i]}" \
    --output-dir "${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAMES[$i]}" \
    --overwrite
done
```

LDSC3 infers the file format, common column names and the source genome build.
Before the real run, add `--infer-only` to one command to preview the inferred columns and get a suggested command; it writes nothing.
If the source and output builds differ, munging lifts the positions over to the output build.
See [munge summary statistics](workflows/munge-sumstats.md) for column flags, liftover and the run summary.

Each trait's output directory contains:

```text
sumstats_processed/
    mdd2025/
        mdd2025.parquet
        diagnostics/
            sumstats.log
            dropped_snps/
                dropped.tsv.gz
```

## Analysis 1: SNP heritability

**Goal:** compute LD scores from the reference panel, then estimate the SNP heritability of `mdd2025`.

**Starting here?** Run [Before you begin](#before-you-begin) and munge `mdd2025` in [Munge summary statistics](#munge-summary-statistics) first.

### Step 1: compute unpartitioned LD scores

```bash
LDSCORE_DIR="${OUTPUT_ROOT}/unpartitioned_ldscore"

ldsc ldscore \
  --snp-identifier chr_pos \
  --genome-build hg19 \
  --r2-matrix-dir "${R2_MATRIX_DIR}/hg19" \
  --ld-wind-cm 1.0 \
  --output-dir "${LDSCORE_DIR}" \
  --overwrite
```

To use the PLINK files instead of the R² matrix, replace `--r2-matrix-dir ...` with `--plink-prefix "${PLINK_PREFIX}"`.
Without annotations, `ldscore` computes LD scores for one all-ones annotation named `base`.

```text
unpartitioned_ldscore/
    metadata.json
    ldscore.baseline.parquet
    diagnostics/
        ldscore.log
```

Pass this whole LD-score directory to the regression commands.
See [calculate LD scores](workflows/ld-scores.md) for LD windows, reference-panel options and memory.

### Step 2: estimate SNP heritability

```bash
TRAIT_NAME="mdd2025"
H2_OUTPUT_DIR="${OUTPUT_ROOT}/h2/${TRAIT_NAME}"

ldsc h2 \
  --sumstats-file "${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAME}/${TRAIT_NAME}.parquet" \
  --ldscore-dir "${LDSCORE_DIR}" \
  --output-dir "${H2_OUTPUT_DIR}" \
  --overwrite
```

```text
h2/mdd2025/
    h2.tsv
    diagnostics/
        ld_score_regression_bins.tsv
        metadata.json
        h2.log
```

`h2.tsv` has one row.
Read `total_h2_obs` and `total_h2_obs_se`, the SNP heritability on the observed scale and its block-jackknife standard error.
`intercept` near 1 indicates little confounding.

To draw the binned diagnostic plot, run `ldsc plot --result-dir "${H2_OUTPUT_DIR}"`.
For a binary trait, `ldsc convert-h2-scale` reports the result on the liability scale without refitting.
See [estimate SNP heritability](workflows/heritability.md) for every column, the liability scale and the options.

## Analysis 2: genetic correlation

**Goal:** estimate the genetic correlation between every pair of the three traits.

**Starting here?** Run [Before you begin](#before-you-begin), munge all three traits in [Munge summary statistics](#munge-summary-statistics), and compute the LD scores in [Analysis 1, Step 1](#step-1-compute-unpartitioned-ld-scores).

### Step 1: estimate genetic correlation

```bash
LDSCORE_DIR="${OUTPUT_ROOT}/unpartitioned_ldscore"
RG_OUTPUT_DIR="${OUTPUT_ROOT}/rg/mdd2025_scz2022_adhd2019"

ldsc rg \
  --sumstats-sources \
    "${OUTPUT_ROOT}/sumstats_processed/mdd2025/mdd2025.parquet" \
    "${OUTPUT_ROOT}/sumstats_processed/scz2022/scz2022.parquet" \
    "${OUTPUT_ROOT}/sumstats_processed/adhd2019/adhd2019.parquet" \
  --ldscore-dir "${LDSCORE_DIR}" \
  --output-dir "${RG_OUTPUT_DIR}" \
  --overwrite
```

`rg` fits every pair of the listed traits.

```text
rg/mdd2025_scz2022_adhd2019/
    rg.tsv
    rg_full.tsv
    h2_per_trait.tsv
    diagnostics/
        metadata.json
        rg.log
```

`rg.tsv` has one row per pair.
Read `rg`, the genetic correlation, with its block-jackknife standard error `rg_se` and its nominal two-sided `p`.
`h2_per_trait.tsv` holds each trait's SNP heritability from the same run.

### Step 2: plot the result

```bash
ldsc plot --result-dir "${RG_OUTPUT_DIR}"
```

This writes `plots/rg_heatmap.png` inside the `rg` result directory: genetic correlations below the diagonal, and each trait's observed-scale SNP heritability on the diagonal.
Plotting reads the saved results and does not refit.
See [estimate genetic correlation](workflows/genetic-correlation.md) for the anchor-trait design and failed pairs, and [plotting](workflows/plotting.md) for how to read the figure.

## Analysis 3: partitioned heritability (cell-type enrichment)

**Goal:** test whether each query annotation adds SNP heritability for `mdd2025` beyond the baseline annotations.

The method is stratified LD score regression (S-LDSC; Finucane et al. 2015): each query annotation is fitted in its own model together with the same baseline annotations.

**Starting here?** Run [Before you begin](#before-you-begin) and munge `mdd2025` in [Munge summary statistics](#munge-summary-statistics) first.
This analysis does not use the LD scores from Analysis 1.

### Step 1: compute LD scores with baseline and query annotations

`ldsc ldscore` reads BED files directly, so no separate annotation step is needed.
Put the query BED files in one directory and select them with a quoted glob:

```bash
BASELINE_ANNOT_SOURCES="${INPUT_ROOT}/resources/1000G_EUR_Phase3_baseline/baseline.@.annot.gz"
QUERY_BED_SOURCES="${INPUT_ROOT}/resources/example_raw_annot/*.bed"
PARTITIONED_LDSCORE_DIR="${OUTPUT_ROOT}/partitioned_ldscore"

ldsc ldscore \
  --snp-identifier chr_pos \
  --genome-build "hg19" \
  --baseline-annot-sources "${BASELINE_ANNOT_SOURCES}" \
  --query-annot-bed-sources "${QUERY_BED_SOURCES}" \
  --r2-matrix-dir "${R2_MATRIX_DIR}/hg19" \
  --ld-wind-cm 1.0 \
  --output-dir "${PARTITIONED_LDSCORE_DIR}" \
  --overwrite
```

As in Analysis 1, `--plink-prefix "${PLINK_PREFIX}"` can replace `--r2-matrix-dir ...`.
These LD scores can be reused for any trait.

```text
partitioned_ldscore/
    metadata.json
    ldscore.baseline.parquet
    ldscore.query.parquet
    ldscore.overlap.parquet
    diagnostics/
        ldscore.log
        query_annotation_status.tsv
```

`diagnostics/query_annotation_status.tsv` records whether each query annotation was kept or skipped.
With more than 100 query annotations, the query LD scores are split into numbered files; see [many query annotations and workers](workflows/ld-scores.md#many-query-annotations-and-workers).
To use gene lists instead of BED files, see [cell-type enrichment](workflows/cell-type-enrichment.md).

### Step 2: estimate partitioned heritability

```bash
TRAIT_NAME="mdd2025"
PARTITIONED_H2_OUTPUT_DIR="${OUTPUT_ROOT}/partitioned-h2/${TRAIT_NAME}"

ldsc partitioned-h2 \
  --sumstats-file "${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAME}/${TRAIT_NAME}.parquet" \
  --ldscore-dir "${PARTITIONED_LDSCORE_DIR}" \
  --output-dir "${PARTITIONED_H2_OUTPUT_DIR}" \
  --overwrite
```

```text
partitioned-h2/mdd2025/
    partitioned_h2.tsv
    diagnostics/
        metadata.json
        partitioned-h2.log
        query_status.tsv
        query_annotations/
            manifest.tsv
            0001_<query-annotation-name>/
                partitioned_h2.tsv
                partitioned_h2_full.tsv
                coefficient_delete_values.parquet
                metadata.json
```

`partitioned_h2.tsv` has one row per query annotation, sorted from the smallest to the largest `coefficient_p`.
Read `coefficient`, the query annotation's per-SNP heritability conditional on the baseline, and its one-sided `coefficient_p`, which tests whether the query adds heritability.
Prefer `coefficient_p` over `enrichment` for query annotations, because enrichment is confounded by overlap with enriched baseline annotations.
Each query's complete baseline-plus-query fit is under `diagnostics/query_annotations/`; the `query_result_dir` column points to it.

By default, any failed query stops the run without publishing results; `--continue-on-query-error` keeps the successful queries.
Keep the whole result directory together when you move it.
See [partitioned heritability](workflows/partitioned-heritability.md) for every output, concurrent fitting and query failures.

## Next steps

- Build reusable annotation files from BED files or gene lists: [annotations](workflows/annotations.md).
- Score many gene lists, optionally with a gene LD-score index: [gene-list LD scores](workflows/gene-list-ld-scores.md) and [cell-type enrichment](workflows/cell-type-enrichment.md).
- Analyze a continuous annotation by quantile: [continuous annotations](workflows/continuous-annotations.md).
- Report a binary trait on the liability scale: [liability scale](workflows/liability-scale.md).
- Plot any saved result: [plotting](workflows/plotting.md).
- Use legacy LDSC files: legacy munged `.sumstats.gz` files go directly to `h2`, `rg` and `partitioned-h2`, and legacy LD-score files need `ldsc convert-legacy-ldscores` first; see [migrating from legacy LDSC](workflows/migrating-from-legacy-ldsc.md).
- Find reference panels and baseline annotations: [reference data](reference-data.md).
