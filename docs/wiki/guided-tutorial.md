# LDSC3 - Guided Analysis Tutorial

This tutorial walks through how to use the `ldsc` package for a series of LD score-based heritability analyses.

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
- raw BED files for the pathways / cell types whose h2 contribution you want to test (the "query annotations").


As a motivating example, we study the `mdd2025` trait, using 1000 Genomes Phase 3 as the reference panel and `1000G_EUR_Phase3_baseline` as the baseline annotations, with `Hippocampus_PP1.bed`, `Cerebellum_PC16.bed`, `Cerebellum_PP3.bed`, and `Cerebellum_PP1.bed` as the pathways whose heritability contribution we test after controlling for the baseline annotations.

All of these resources can be found under `ldsc3_test_bundle/resources/`, except for the raw GWAS sumstats files.

Before running the commands below, follow the README to clone the `main` branch of the GitHub repo and install it in a suitable conda environment. Then run `conda activate ldsc3` before the analysis.

We first set up the input and output directories:

```bash
INPUT_ROOT="/users/w/e/wenbinwu/Sullivan/LDSC/ldsc3_test_bundle"
OUTPUT_ROOT="/users/w/e/wenbinwu/Sullivan/LDSC/ldsc3_test_bundle/tutorial_output"
```

### Remarks

1. By default, all regressions use HapMap3 SNPs. These SNPs are bundled with the package and can be enabled with a single flag in each command (see below).
2. For flags ending in `-sources`, you can use the glob pattern `*` to match multiple files and `@` as the placeholder for the chromosome number.
3. Output directories are created automatically. Use `--overwrite` to allow overwriting existing output files.
4. `ldsc ldscore` excludes the extended MHC and centromere regions by default.
5. Be careful about the genome build of your input file. Although the package runs a guardrail check on the genome build, it is not comprehensive.
6. To inspect what happened under the hood, check the log file in the `diagnostics/` directory inside each command's output directory.
7. The memory and run-time recommendations below are NOT accurate (codebase not yet benchmarked); treat them as rough upper bounds. If you hit OOM, increase the memory allocation. 24 GB should be sufficient for the initial test runs.

## Config you should follow across all steps

- `snp_identifier=chr_pos`: use `chr` + `pos` as the unique identifier for a SNP.
- ... [TODO]

### TODO

- Results and plots interpretation (which columns to use in each scenario, and what they mean).
- How to reuse previously generated annotations for partitioned LDSC.
- Refine memory and run-time numbers with proper benchmarking rather than guessing from log files. In particular, the SLURM memory figure for `ldscore` is inaccurate (it somehow always reports the allocated memory minus 2 MB).

## Munge-sumstats

This section illustrates the `ldsc munge-sumstats` command.

**Goal:** convert raw summary statistics files into a standardized format for use in downstream ldsc regression.

**Recommended memory allocation:** < 4 GB

**Expected running time:** < 5 min

**Minimal command:**

```bash
TRAIT_NAME="mdd2025"  # used as the output dir name and as a trait label in downstream regression outputs
RAW_SUMSTATS_FILE="/users/w/e/wenbinwu/Sullivan/LDSC/data/pgc_GWAS/mdd/mdd2025/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv" 
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
  --use-hm3-regression-snps \
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
        metadata.json
        h2.log
```

## Analysis 2: estimate cross-trait genetic correlation between multiple traits

**Goal:** estimate the genetic correlation between traits using the GWAS sumstats and LD scores.

**Recommended memory allocation:** Step 1: reuses Analysis 1 LD scores; Step 2: < 4 GB

**Expected running time:** Step 1: reuses Analysis 1 LD scores; Step 2: < 1 min

### Munge two additional sumstats

```bash
# Munge two other traits.
TRAIT_NAMES=("scz2022" "adhd2019")
RAW_SUMSTATS_FILES=(
  "/users/w/e/wenbinwu/Sullivan/LDSC/data/pgc_GWAS/scz/scz2022/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.cleaned.tsv"
  "/users/w/e/wenbinwu/Sullivan/LDSC/data/pgc_GWAS/adhd/adhd2019/daner_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta"
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
  --use-hm3-regression-snps \
  --output-dir "${PARTITIONED_LDSCORE_OUTPUT_DIR}" \
  --ld-wind-cm 1.0 \
  --overwrite
```

**Remarks:**

- For PLINK input, replace `--r2-dir ...` with `--plink-prefix "${PLINK_PREFIX}"`.
- For partitioned h2 analysis, `ldsc ldscore` accepts BED files directly, so you do not need to build the annotations yourself to run this analysis. If you do want to generate annotations for other purposes, follow the *Make annotations* section below.

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
  --write-per-query-results \
  --overwrite
```

**Remarks:**

1. `partitioned-h2` reports the total h2 implied by the partitioned model; this is not necessarily identical to the standalone unpartitioned h2 estimate above.
2. `--write-per-query-results` saves the results for each query annotation separately, under `diagnostics/query_annotations/`.

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
                metadata.json
                ...
```

Here:

- `partitioned_h2.tsv` is a summary over all query annotations, sorted by the p-value of the query annotation's coefficient (`coefficient-p`), from most to least significant.
- `query_annotations/` holds the per-query results, which also include the coefficients and h2 share of the baseline annotations. `manifest.tsv` records each original query annotation name, its results folder name, and the paths to the relevant results files.

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

- All input query annotations are combined into one file as separate columns, with each column name being the input's base filename (used as the annotation name).
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

Munged sumstats, legacy annotations.

[TODO]
