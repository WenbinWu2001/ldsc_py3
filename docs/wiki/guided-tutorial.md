

# LDSC3 - Guided Analysis Tutorial

In this tutorial, we illustrate how to use our ldsc3 package for a series of LD score-based heritability analysis.

The analysis pipeline involves:

- **Before all analysis:** curate a standard sumstats file with `ldsc munge-sumstats`
- **Analysis 1:** estimate heritabiliy of a trait
  - Step 1: compute (unpartitioned) LD scores with `ldsc ldscore`
  - Step 2: use `ldsc h2` to estimate trait heritability.

- **Analysis 2:** estimate cross-trait genetic correlation between multiple traits
  - Step 1: compute (unpartitioned) LD scores with `ldsc ldscore` (same as in Analysis 1)
  - Step 2: use `ldsc rg` to estimate cross-trait genetic correlation.


- **Analysis 3:** partition heritability with functional annotations (known as cell-type-specific h2 regression in ldsc2)
  - Step 1: compute partitioned LD scores with `ldsc ldscore` with query annotations
  - Step 2: use `ldsc partitioned-h2` to estimate heritability contribtion of each query annotation



The files you would need is:

- raw GWAS summary statistics files for traits of interest (or the munged sumstats from legacy ldsc2 codebase)
- a reference panel for LD score calculation (either R2 parquets or PLINK suite)
- a set of baseline annotations
- raw bed files for the pathways / cell types you want to test for h2 contribution ("query annotations")



Useful rules: for flags ending with `-sources`, You can use glob pattern `*` to match multiple files and `@` as the placeholder for chromosome number.



As a motivating example, we study `mdd2025` trait, use 1000 Genomes Phase 3 as the reference panel, baseline annotations as `1000G_EUR_Phase3_baseline`, with `Hippocampus_PP1.bed`, `Cerebellum_PC16.bed`,  `Cerebellum_PP3.bed`, `Cerebellum_PP1.bed` as the pathways to test heritability contribution after controling baseline pathways. 

All these resources can be found under `ldsc3_test_bundle/resources/`, except for the raw GWAS sumstats files.





Before running the commands below, please follow the readme page to clone the `main` branch of the github repo and install it with a proper conda environment. 

Run `conda activate ldsc3` before the following commands.





We first set up the input and output directory:

```bash
INPUT_ROOT="/users/w/e/wenbinwu/Sullivan/LDSC/ldsc3_test_bundle"
INPUT_ROOT="/users/w/e/wenbinwu/Sullivan/LDSC/ldsc3_test_bundle/tutorial_output"
```




### Remarks:
1. By default, we use Hapmap3 SNPs in all regression. These snps are bundled with the package and can be easily used with a flag in each functionality (see below).
2. Use `*` for glob patterns and `@` as the chromosome placeholder.
3. Output directories are created automatically. Use `--overwrite` to allow overwriting existing output files.
4. `ldsc ldscore` excludes extended MHC and centromere regions by default.
5. Be careful of the genome build of your input file! Although this package does a guardrail check on the genome build, it is not comprehensive.
6. To track the details under the hood, you can check the log file under `diagnostics/` directory under the output directory in each command.
7. The memory usage and run time recommendations below are not accurate as we didn't do any benchmarking yet. Expect them as the "upper bound". If you hit OOM, increase memory allocation. 24GB should be more than enough for our initial test runs.





## Config you should follow across all steps

- `snp_identifier=chr_bp`: Use `chr` + `bp` as the unique identifier for snp.
- ... [TODO]



### TODO

- results and plots interpretation (which columns to use in each scenario, what does it mean)
- how to use the previous annotations for pldsc
- refine memory usage and running time with refined benchmarking rather than guessing from log files. In particular, memory usage from slurm job for `ldscore` is inaccurate (always shows uses allocated memory minus 2MB.)



## Munge-sumstats

This section illustrates the usage of `ldsc munge-sumstats` command. 

**Goal:** Convert raw summary statistics files to a standard format, used in downstream ldsc regression.

**recommended memory allocation:** <4GB

**expected running time:** <5min

**minimal command:**

```bash
TRAIT_NAME="mdd2025"  # used as the output dir name and as a trait label in downstream regression outputs
RAW_SUMSTATS_FILE="/users/w/e/wenbinwu/Sullivan/LDSC/data/pgc_GWAS/mdd/mdd2025/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv" 
SUMSTATS_OUT_DIR="${OUTPUT_ROOT}/sumstats_processed"

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

1. The package auto-infers genome build, raw sumstats file format, and common column names, and automatically strips off any comment lines starting with `#` in the file header. You can still pass the legacy flags like `--snp` to specify the column name, if your sumstats file has a column names. For a table of column aliases that this package auto-infers, see the detailed wiki on `munge-sumstats` module.
2. Before you actually run the command, we strongly recommend you use `--infer-only` to preview inferred columns and suggested flags without writing outputs. This is also very helpful if you are unsure about what flags to use or about the genome build of your input files. The program will give you suggestions on the correct command to use. See more in the `munge-sumstats` wiki.
  ```bash
  ldsc munge-sumstats \
  ... (your-guessed-flags) \
  --infer-only
  ```
3. Legacy `.sumstats.gz` is supported. Specify them via `--raw-sumstats-file` and pass a `--trait-name` so that downstream outputs use stable trait labels.
3. If source and output builds differ, add either `--use-hm3-quick-liftover` with `--use-hm3-snps`, or provide `--liftover-chain-file`.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
sumstats_processed/
#     mdd2025/
#       sumstats.parquet
#       diagnostics/sumstats.log
#       diagnostics/dropped_snps/dropped.tsv.gz
```



## Specify the ref panel

Before the subsequeny regression in any analysis below, we first need to specify path of the reference panel:

```bash
# LD reference panel from 1000 Genomes Phase 3. Use either R2 parquet files or PLINK files.
R2_DIR="${INPUT_ROOT}/resources/ldref/ldref_r2_1kgph3_1cM"
PLINK_PREFIX="${INPUT_ROOT}/resources/1000G_EUR_Phase3_plink/1000G.EUR.QC."  # hg19
```

Here we provide two sources of reference panels: 

- per-chromosome R2 parquets calculated from the PLINK.
- PLINK suite.



## **Analysis 1:** estimate heritabiliy of a trait

**Goal:** Calculate LD scores using genomic reference panel, and estimate heritability of the trait using the munged sumstats file.

**recommended memory allocation:** Step 1: 24GB (enough, for safety only);Step 2: <4GB

**expected running time:** Step 1: <2h; Step 2: <1min

### Step1: compute (unpartitioned) LD scores with `ldsc ldscore`

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

If you are using the PLINK suite as reference panel, please replace `--r2-dir ...` in the command with

```bash
`--plink-prefix "${PLINK_PREFIX}"
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

### Step 2: use `ldsc h2` to estimate trait heritability.

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



## **Analysis 2:** estimate cross-trait genetic correlation between multiple traits

**Goal:** Estimate the genetic correlation between two traits using the gwas sumstats and LD structrue (ld scores).

**recommended memory allocation:** Step 1: 24GB (enough, for safety only);Step 2: <4GB

**expected running time:** Step 1: <2h; Step 2: <1min

### Munge another two sumstats

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

### Step 2: use `ldsc rg` to estimate cross-trait genetic correlation.

This command calculates the cross-trait genetic correlation for *all* pairs of the input sumstats.

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





## **Analysis 3:** partition heritability with functional annotations (known as cell-type-specific h2 regression in ldsc2)

**Goal:** Partitioned the heritability of a trait as sums of contributions of a set of pathway, testing whether the query pathway contributes to heritability when controlling the baseline pathways.

**recommended memory allocation:** Step 1: 24GB (enough, for safety only);Step 2: <4GB

**expected running time:** Step 1: <1h; Step 2: <1min

### Step 1: compute partitioned LD scores with `ldsc ldscore` with query annotations

We need to also specify the path of the raw bed files (here, four bed files) and the baseline annotations:

```bash
# annotations
RAW_QUERY_BED_SOURCES="$INPUT_ROOT/resources/example_raw_annot/*.bed"  # put all your query BED files in one directory
BASELINE_ANNOT_SOURCES="${INPUT_ROOT}/resources/1000G_EUR_Phase3_baseline/baseline.@.annot.gz"
```

Please curate all your query raw files under the same directory and use glob pattern as above to match all of them. This would allow their LD scores to be computed in parallel in `ldsc ldscore`.

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
- For partitioned h2 analysis, `ldsc ldscore` accepts bed files as input, so you do not have to make annotations on your own in order to run this functionality. However, if you want to generate the annotations for other purposes, follow the *make annotations* section below.

**Caveat:** If you are using `--ld-wind-cm`, please make sure your PLINK suite contains non-missing genetic coordinates (third column in `.bim` file). If they are missing (e.g., all zeros), the program will throw an error.

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

### Step 2: use `ldsc partitioned-h2` to estimate heritability contribtion of each query annotation

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

1. `partitioned-h2` reports the total h2 implied by the partitioned model; this is not necessarily identical to the standalone unpartitioned h2 estimate below.
1. `--write-per-query-results` means to save the results for each query annotation separately, under `diagonstics/query_annotations/`.

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

Here,

- `partitioned_h2.tsv` is a summary over all query annotations, sorted by p-value of the coefficient of the query annotation (`coefficient-p`) from most significant to least significant.
- `query_annotations/` contains the per-query results, which also contains the coefficients and h2 portion of the baseline annotations. `manifest.tsv `records each original query annotation name, the results folder name, and the relevant paths of the results files.



## make annotations

**Goal:** Curate binary annotations from bed files, given a set of baseline annotations files.

**recommended memory allocation:** <8GB for four query bed files

**expected running time:** <15min for four query bed files

It is not used in the analysis pipeline. However, we preserve this functionality for convenience in case users need it for other purposes.

``` bash
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

- All input query annotations are summarized together as differerent columns, with columns name being the base filename (as the annotation name).
- in annotations CM is preserved for backward compatibility but not used

- See more in the `annotate` wiki for this functionality.

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



## Backward compatibility to legacy ldsc2 codebase

munged sumstats, legacy annotations

[TODO]


