# LDSC-SEG for Protein-Coding Gene Lists

Last updated on: 2026-08-16

This tutorial tests whether one or more protein-coding gene lists are enriched for trait heritability. For each query gene list, LDSC3 fits the model

`[query annotation, baseline annotations]`.

Use the query category's coefficient to assess whether the gene list contributes additional heritability after accounting for the baseline annotations. The commands below omit the optional control-gene file, so no control annotation is added.

The workflow has two steps:

1. Calculate partitioned LD scores using either *direct mode* or *fast indexed mode*.
2. Run partitioned-heritability regression for the target trait.

## Set up inputs

```bash
PROJECT_ROOT="/users/w/e/wenbinwu/Sullivan/LDSC/2026-08-10_sldsc_pc_genes_index_suite"
GENE_LIST_SOURCES="${PROJECT_ROOT}/ldsc3/query_gene_lists/*.txt"
TRAIT_NAME="trait_name"
SUMSTATS_FILE="/path/to/${TRAIT_NAME}.sumstats.gz"
```

Each query gene-list file must contain one Ensembl gene ID or gene name per line, without a header. Glob patterns are supported, and each matched gene list is tested separately.

The summary statistics may be either an LDSC3 sumstats artifact or a legacy LDSC2 `.sumstats` or `.sumstats.gz` file.

## Step 1: calculate partitioned LD scores

Choose either direct mode or fast mode. Both produce the same canonical LD-score output format for the regression in Step 2.

### Option A: direct mode

Direct mode reads the baseline annotations and PLINK reference panel, constructs the query annotations in memory, and calculates their LD scores and annotation counts. The query annotations are not written to disk.

- **Expected runtime:** approximately 1.5 hours
- **Expected memory usage:** approximately 10 GB

```bash
RESULT_ROOT="${PROJECT_ROOT}/ldsc3/example_output/direct_mode"

BASELINE_ANNOT_SOURCES="/users/w/e/wenbinwu/Sullivan/LDSC/data/resources_from_Jerry/partitioned_LDSC/1000G_EUR_Phase3_baseline/baseline.@.annot.gz"
PLINK_PREFIX="/users/w/e/wenbinwu/Sullivan/LDSC/data/resources_from_Jerry/partitioned_LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.@"
GENE_COORDINATE_FILE="/path/to/gene-coordinates.hg19.tsv.gz"
PARTITIONED_LDSCORE_DIR="${RESULT_ROOT}/pldsc/ldscore"

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
  --gene-list-resolution-policy strict \
  --output-dir "${PARTITIONED_LDSCORE_DIR}" \
  --overwrite
```

Flags used in this command:

- `--query-annot-gene-list-sources` specifies one or more query gene-list files. Glob patterns are supported.
- `--baseline-annot-sources` specifies the baseline annotation suite. Use `@` as the chromosome-number placeholder.
- `--plink-prefix` specifies the PLINK reference-panel prefix. Pass the shared plain stem (for example, `1000G.EUR.QC.` for `1000G.EUR.QC.1.{bed,bim,fam}`); the PLINK suite must correspond to the baseline annotation suite. The older `@` placeholder remains supported.
- `--snp-identifier rsid` matches SNPs by rsID.
- `--genome-build hg19` specifies the build used for gene projection and named region definitions.
- `--ld-wind-cm 1.0` calculates LD within a 1-cM window.
- `--padding-bp 100000` adds 100 kb to either side of each gene interval.
- `--gene-coordinate-file` is the sole one-based hg19 catalog for focal and control lists; there is no packaged fallback.
- `--gene-list-resolution-policy strict` stops before LD-score work if any submitted identifier is unresolved or ambiguous. Use `resolved-only` only for deliberate exploratory subset analysis.
- No control-gene annotation is added by default. To add one, pass a single existing file with `--control-gene-list-file`; it must contain one Ensembl gene ID or gene name per line, without a header.
- `--gene-exclude-regions mhc` removes query genes whose unpadded intervals overlap the MHC before padding and projection. If a control gene list is enabled, the same filter also applies to it. This is a gene-level filter; it does not remove SNPs from the LD reference panel.
- `--output-dir` specifies the LD-score output directory.
- `--overwrite` permits replacement of existing result artifacts. Use it with caution.

Relevant flags omitted because their default values are used:

- `--regr-snps-exclude-regions mhc-and-centromeres` removes MHC and pericentromeric SNPs from the LD-score rows subsequently used in regression. It does not remove these SNPs from the LD reference universe or change annotation count and overlap statistics.
- `--common-maf-min 0.05` defines common SNPs as SNPs with MAF >= 0.05 when calculating common-SNP annotation counts.

### Option B: fast indexed mode

Fast mode uses a precomputed exact gene LD-score index to assemble LD scores, annotation counts, and overlap statistics for the query gene lists. It does not reread the PLINK reference panel.

- **Expected runtime:** approximately 10 minutes
- **Expected memory usage:** approximately 6 GB
- **Supported genome build:** hg19

```bash
RESULT_ROOT="${PROJECT_ROOT}/ldsc3/example_output/fast_mode"

GENE_LDSCORE_INDEX_DIR="${PROJECT_ROOT}/ldsc3/1000G_EUR_Phase3_baseline__hg19__pad100000bp__ld1.0cM"
PARTITIONED_LDSCORE_DIR="${RESULT_ROOT}/pldsc/ldscore"

ldsc ldscore \
  --query-annot-gene-list-sources "${GENE_LIST_SOURCES}" \
  --gene-ldscore-index-dir "${GENE_LDSCORE_INDEX_DIR}" \
  --output-dir "${PARTITIONED_LDSCORE_DIR}" \
  --overwrite
```

Flags used in this command:

- `--query-annot-gene-list-sources` specifies one or more query gene-list files. Glob patterns are supported.
- `--gene-ldscore-index-dir` specifies the precomputed gene LD-score index. The index contains the fixed baseline LD scores, embedded gene catalog, gene-to-SNP projection data, annotation-count information, regression weights, and the operators needed to assemble gene-list LD scores.
- No control-gene annotation is added by default. Fast mode accepts one existing control-gene file through `--control-gene-list-file`; multiple files, glob patterns, and sentinel values are not supported. An enabled control inherits the index's gene catalog, build, padding, and gene-exclusion policy.
- `--output-dir` specifies the LD-score output directory. The completed directory is self-contained; the index is not needed by the regression in Step 2.
- `--overwrite` permits replacement of existing result artifacts. Use it with caution.

Fast mode inherits the baseline annotations, reference SNP universe, SNP identity mode, genome build, LD window, gene padding, gene-region exclusion, regression-SNP selection, SNP-region exclusion, MAF rules, annotation counts, and regression weights from the index. The index metadata is the authoritative record of these settings.

Do not supply live overrides such as `--baseline-annot-sources`, `--plink-prefix`, `--r2-dir`, `--snp-identifier`, `--genome-build`, `--ld-wind-*`, `--padding-bp`, `--gene-exclude-regions`, or `--regr-snps-file` in indexed mode.

### LD-score outputs

Both modes write:

```text
ldscore/
    metadata.json
    ldscore.baseline.parquet
    ldscore.query.parquet
    ldscore.overlap.parquet
    diagnostics/
        ldscore.log
        query_annotation_status.tsv
        gene_list_audit.tsv.gz
        gene_list_resolution_summary.tsv
```

Use `diagnostics/ldscore.log` to monitor progress. Start curation with
`gene_list_resolution_summary.tsv`, then filter `gene_list_audit.tsv.gz`; after
catalog preflight, `query_annotation_status.tsv` explains run-specific focal
skips. See [Gene-list diagnostics and repair](../current/gene-list-diagnostics-and-repair.md)
for the detailed repair workflow.

Because no control-gene file is supplied, `ldscore.baseline.parquet` contains only the supplied baseline categories; no `gene_control` column is added.

## Step 2: run partitioned-heritability regression

Use the LD-score directory produced by either direct or fast mode in Step 1.

```bash
PARTITIONED_H2_DIR="${RESULT_ROOT}/pldsc/regr/${TRAIT_NAME}"

ldsc partitioned-h2 \
  --sumstats-file "${SUMSTATS_FILE}" \
  --ldscore-dir "${PARTITIONED_LDSCORE_DIR}" \
  --output-dir "${PARTITIONED_H2_DIR}" \
  --write-per-query-results \
  --overwrite
```

Flags used in this command:

- `--sumstats-file` specifies the munged summary statistics for one trait. LDSC3 artifacts and legacy LDSC2 `.sumstats` or `.sumstats.gz` files are supported.
- `--ldscore-dir` specifies the output directory from Step 1.
- `--output-dir` specifies the partitioned-heritability result directory.
- `--write-per-query-results` writes a detailed result directory for each query under `diagnostics/query_annotations/`. Without this flag, only the aggregate `partitioned_h2.tsv` is written.
- `--overwrite` permits replacement of existing result artifacts. Use it with caution.

Relevant flags omitted because their default values are used:

- `--count-kind common` uses annotation counts calculated from common reference SNPs.
- `--n-blocks 200` uses 200 blocks for block-jackknife standard-error estimation.
- The heritability intercept is estimated because neither `--intercept-h2` nor `--no-intercept` is specified.
- `--samp-prev` and `--pop-prev` are omitted, so results are reported on the observed scale without liability-scale conversion.
- `--summary-sort-by auto` sorts query gene lists by ascending `coefficient_p`. Other choices are `category`, `prop-snps`, `prop-h2`, `enrichment`, `enrichment-p`, `coefficient`, and `coefficient-p`. P-value sorts are ascending, other numeric sorts are descending, and missing values are placed last.

### Regression outputs

```text
regr/<trait>/
    partitioned_h2.tsv
    diagnostics/
        metadata.json
        partitioned-h2.log
        query_annotations/
            manifest.tsv
            0001_<query-name>/
                partitioned_h2.tsv
                partitioned_h2_full.tsv
                metadata.json
```

The root `partitioned_h2.tsv` contains one summary row per query gene list and omits baseline-category rows. Each query-specific `partitioned_h2_full.tsv` reports the complete fitted model:

`[query annotation, baseline annotations]`.

Focus on `coefficient` and its one-sided `coefficient_p`. A positive coefficient with a small p-value indicates that the query gene list contributes additional heritability after accounting for the baseline annotations.

## Differences from the LDSC2 workflow

- LDSC2 calculates new LD scores only for the query annotation and uses precomputed LD scores, annotations, allele-frequency files, and regression weights for the baseline categories.
- LDSC3 direct mode calculates LD scores for both baseline and query annotations and derives the required count and weight information from the PLINK reference panel in the same run.
- LDSC3 fast mode retrieves the fixed baseline information from the index and assembles query LD scores from precomputed exact components.
- By default, LDSC3 starts with the packaged HapMap3 regression SNP set and removes MHC and pericentromeric SNPs from the rows used in regression. LDSC2 uses the provided `w_hm3.print_snps` and does not automatically apply an additional explicit MHC or pericentromeric-region exclusion step.
