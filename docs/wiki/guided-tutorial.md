# LDSC3 - Guided Analysis Tutorial

Last updated on: 2026-09-15

This tutorial walks through how to use the `ldsc` package for a series of LD score-based heritability analyses.

## Overview

If you already have a reusable legacy LD-score suite, use the [convert-ldsc2-ldscores utility](utility-functionalities/convert-ldsc2-ldscores.md) for the standard filenames and conversion flags. It provides an alternative to computing unpartitioned LD scores for `h2`/`rg`, or imports a baseline-only model for `partitioned-h2`; query/cell-type scores still use the current LD-score workflow.

The analysis pipeline involves:

- **Before any analysis:** curate a standardized sumstats file with `ldsc munge-sumstats`.
- **Analysis 1:** estimate the heritability of a trait.
  - Step 1: compute (unpartitioned) LD scores with `ldsc ldscore`.
  - Step 2: estimate trait heritability with `ldsc h2`.
- **Analysis 2:** estimate cross-trait genetic correlation between multiple traits.
  - Step 1: compute (unpartitioned) LD scores with `ldsc ldscore` (same as in Analysis 1).
  - Step 2: estimate cross-trait genetic correlation with `ldsc rg`.
  - Step 3: plot correlations and single-trait observed heritabilities with `ldsc plot`.
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

These commands target the current `restructure` development branch. Follow the [README development installation](../../README.md#development-version), then run `conda activate ldsc3-dev`. The default installation includes Matplotlib and exposes `ldsc plot`; there is no separate plotting extra. Run the Bash examples below in one Bash session.

We first set up the input and output directories:

```bash
INPUT_ROOT="/path/to/ldsc3_test_bundle"
OUTPUT_ROOT="${INPUT_ROOT}/tutorial_output"
```

### Remarks

1. By default, `ldsc ldscore` writes regression rows from the bundled HapMap3 set, after applying its default MHC-and-centromere exclusion. `h2`, `rg`, and `partitioned-h2` consume those rows and have no HapMap3 flag. `munge-sumstats` also restricts to packaged HapMap3 by default; use `--sumstats-snps-file FILE` to replace its keep-list or `--no-snp-restriction` to disable it while retaining ordinary QC.
2. Quote input patterns so LDSC receives them intact. Baseline and prebuilt annotation sources support `*` and complete-autosomal `@` suites; BED, gene-list, and rg sumstats sources support globs but not `@`. Genetic-map sources for `ldscore` and `build-gene-ldscore-index` accept comma-separated exact paths only. See the [command-specific path rules](../current/path-specification.md#pattern-support-in-command-help).
3. Output directories are created automatically. Use `--overwrite` to allow overwriting existing output files.
4. `ldsc ldscore` excludes the extended MHC and centromere regions by default.
5. Be careful about the genome build of your input file. Although the package runs a guardrail check on the genome build, it is not comprehensive.
6. To inspect what happened under the hood, check the log file in the `diagnostics/` directory inside each command's output directory.
7. Resource use depends on SNP, sample, annotation, and worker counts. Start with `--threads 1`, allow temporary disk space under the output directory, and lower `--query-batch-size` to reduce active query workspace. Direct mode repeats reference work across batches; indexed mode keeps one chromosome operator per worker while writing that chromosome's batches. In the [1,000-query chromosome-22 benchmark](../audits/annotation-memory/sequential-query-batches.md#resource-measurements-and-limits), width 100 used 82.2% less peak private scratch than width 1,000 and took 1.903 times as long, with essentially unchanged peak RAM. The report separates measured results from conditional whole-genome projections. The [earlier local benchmark](../audits/annotation-memory/results.md) predates sequential output batching.

## Config you should follow across all steps

- `snp_identifier=chr_pos`: use `chr` + `pos` as the unique identifier for a SNP.
- Use `hg19` consistently for the coordinate-based examples below. Choose the output build explicitly during munging and provide matching baseline, BED, and reference-panel coordinates.
- These examples explicitly select base `chr_pos`; the general package default is `chr_pos_allele_aware`. Use the same identity mode across materialized inputs. See [global configuration](main-functionalities/global-config.md).

## Munge-sumstats

This section illustrates the `ldsc munge-sumstats` command.

**Goal:** convert raw summary statistics files into a standardized format for use in downstream ldsc regression.

**Minimal command:**

```bash
TRAIT_NAME="mdd2025"  # used for the output directory, data filename, and trait metadata
RAW_SUMSTATS_FILE="/path/to/mdd2025_raw_sumstats.tsv"
SUMSTATS_OUT_DIR="${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAME}"

ldsc munge-sumstats \
  --snp-identifier chr_pos \
  --source-genome-build "auto" \
  --output-genome-build "hg19" \
  --raw-sumstats-file "${RAW_SUMSTATS_FILE}" \
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
4. Always choose the output genome build explicitly. If it differs from the resolved source build, packaged HM3 uses quick liftover automatically from package-bundled reference HM3 metadata. Supply `--liftover-chain-file FILE` to disable quick liftover and use the chain instead. Custom-list or unrestricted cross-build runs require a chain; matching builds need no mapping, and unresolved source builds stop with a request for an explicit source build. See the [method-selection table](main-functionalities/munge-sumstats.md#genome-build-conversion).
5. A successful run prints the selected restriction, mapping method, and mapping/drop counts to stdout and `diagnostics/sumstats.log`, including at `--log-level ERROR`. Mapping counts apply after earlier QC and keep-list filtering; separate whole-run stage counts reconcile all parsed input rows with the final retained rows.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
sumstats_processed/
    mdd2025/
        mdd2025.parquet
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
SUMSTATS_FILE="${SUMSTATS_OUT_DIR}/${TRAIT_NAME}.parquet"
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

For a binary trait, `ldsc convert-h2-scale` converts a saved h2 result using supplied sample and population prevalences and writes below `<h2-result-dir>/postprocessing/liability-scale/`. It does not refit regression. See [heritability scale conversion](main-functionalities/h2.md#convert-a-saved-estimate-to-liability-scale).

The bin table records the exact fitted-data summary used by the binned LD Score regression diagnostic. Matplotlib is included in the default package installation; run `ldsc plot --result-dir "${H2_OUTPUT_DIR}"` to create the figure explicitly.

## Analysis 2: estimate cross-trait genetic correlation between multiple traits

**Goal:** estimate the genetic correlation between traits using the GWAS sumstats and LD scores.

### Munge two additional sumstats

```bash
# Munge two other traits.
RG_TRAIT_NAMES=("scz2022" "adhd2019")
RG_RAW_SUMSTATS_FILES=(
  "/path/to/scz2022_raw_sumstats.tsv"
  "/path/to/adhd2019_raw_sumstats.tsv"
)  # change raw sumstats file path here

for i in "${!RG_TRAIT_NAMES[@]}"; do
  RG_TRAIT_NAME="${RG_TRAIT_NAMES[$i]}"
  RG_RAW_SUMSTATS_FILE="${RG_RAW_SUMSTATS_FILES[$i]}"
  RG_SUMSTATS_OUT_DIR="${OUTPUT_ROOT}/sumstats_processed/${RG_TRAIT_NAME}"

  ldsc munge-sumstats \
    --snp-identifier chr_pos \
    --source-genome-build "auto" \
    --output-genome-build "hg19" \
    --raw-sumstats-file "${RG_RAW_SUMSTATS_FILE}" \
    --trait-name "${RG_TRAIT_NAME}" \
    --output-dir "${RG_SUMSTATS_OUT_DIR}" \
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
SUMSTATS_SOURCES="${OUTPUT_ROOT}/sumstats_processed/*/*.parquet"  # match the trait-named sumstats files

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

### Step 3: plot genetic correlations and trait heritabilities

```bash
ldsc plot --result-dir "${RG_OUTPUT_DIR}"
```

This writes `plots/rg_heatmap.png` below the rg result directory, with plot metadata and a log under `plots/diagnostics/`. Lower-triangle cells show genetic correlation with jackknife SE in parentheses. Light-gray diagonal cells show each trait's observed-scale heritability with its jackknife SE underneath, using two decimal places. Add `--overwrite` when regenerating an existing figure.

The diagonal uses the saved single-trait estimates in `h2_per_trait.tsv`, matched by trait name. It retains finite estimates outside 0–1; missing or unusable estimate–SE pairs display `failed`. Liability columns do not change the displayed h2 scale. Plotting reads the saved results without rerunning regression.

For a result produced with `ldsc rg --anchor-trait ...`, the same plot command instead writes `plots/rg_anchor_forest.png`: correlations have points and one-SE bars, partner h2 estimates appear in a separate text column, and the anchor's h2 appears in the subtitle. See [the rg wiki page](main-functionalities/rg.md) and [the plotting manual](../../tutorials/plotting-results.md) for interpretation and malformed-input behavior.

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

### Step 1: compute partitioned LD scores with `ldsc ldscore` using query annotations

We also need to specify the paths to the raw BED files (here, four of them) and to the baseline annotations:

```bash
# annotations
RAW_QUERY_BED_SOURCES="$INPUT_ROOT/resources/example_raw_annot/*.bed"  # put all your query BED files in one directory
BASELINE_ANNOT_SOURCES="${INPUT_ROOT}/resources/1000G_EUR_Phase3_baseline/baseline.@.annot.gz"
```

Place your raw query files in the same directory and use a quoted glob to select them together. Direct LDSC prepares, computes, writes, and releases one query batch before starting the next, repeating reference work across batches. The saved baseline artifact is shared. Chromosome concurrency is controlled separately by `--threads`, which defaults to 1 and is capped at the chromosome count in both direct and indexed calculation. See [LD-score memory controls](main-functionalities/ldscore.md#memory-for-many-pathways) for the implementation sources and tradeoffs.

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
  index](utility-functionalities/build-gene-ldscore-index.md) and [Calculate LD
  scores for gene lists with an index](main-functionalities/ldscore-from-gene-list.md).
  Per-chromosome `Finished` lines report durable private staging; the index
  becomes public only after complete reload validation. Stages cannot be
  resumed or used for incremental chromosome updates.

**Genetic coordinates:** a cM window with PLINK input needs informative BIM CM values or a matching explicit genetic map through `--genetic-map-hg19-sources` (comma-separated exact paths). With R² input, CM comes from the required matching `chrN_meta.tsv.gz` sidecars. Genetic-map flags do not replace those sidecar values. See [effective CM coordinates](main-functionalities/ldscore.md#effective-cm-coordinates-and-metadata-export).

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

The example's four queries fit in one batch and therefore share `ldscore.query.parquet`. When the query count exceeds `--query-batch-size` (default 1000), the directory instead contains `ldscore.query.batch00001.parquet`, `ldscore.query.batch00002.parquet`, and so on. Each covers the same genome-wide SNP rows. The required `metadata.json.query_batches` manifest records file membership and chromosome row groups. Current readers require this manifest; regenerate older LD-score directories. Regression can select columns across these files with a different batch width.

### Step 2: estimate each query annotation's heritability contribution with `ldsc partitioned-h2`

```bash
# Step 2: Select the focal trait explicitly, including when running this section alone.
TRAIT_NAME="mdd2025"
SUMSTATS_FILE="${OUTPUT_ROOT}/sumstats_processed/${TRAIT_NAME}/${TRAIT_NAME}.parquet"
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
2. Query-annotation runs automatically save each successful query's complete baseline-plus-query fit under `diagnostics/query_annotations/`; `--write-per-query-results` has been removed and must be omitted.
3. Error handling is strict by default: any query error prevents publication of new scientific results. Add `--continue-on-query-error` to skip failed query fits and publish successful ones. Check `diagnostics/query_status.tsv` for every attempted query and `diagnostics/partitioned-h2.log` for tracebacks and available block diagnostics. Shared input/output failures, interrupts, and all-query-failed scans still stop. See [query error handling](main-functionalities/partitioned-h2.md#choose-what-happens-when-a-query-fails), implemented by `RegressionRunner.estimate_partitioned_h2_batch()` in [regression_runner.py](../../src/ldsc/regression_runner.py).

4. The command above uses the default inline fit. With four allocated CPUs, add `--threads 4` to fit the four query models concurrently; each process fits one complete model. For larger scans, `--query-batch-size 100` limits loaded query columns and caps concurrency. Parsing and worker-count resolution match `ldscore` and index construction, with queries/batch width replacing the chromosome work cap. See [concurrent query fitting](main-functionalities/partitioned-h2.md#fit-queries-concurrently) for negative counts, numerical-thread limits, memory costs, and failure behavior.

**Outputs:**

Upon a successful run, you should expect the following files in your output directory:

```
partitioned-h2/mdd2025/
    partitioned_h2.tsv
    diagnostics/
        metadata.json
        partitioned-h2.log
        query_status.tsv
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
- Do not create one chromosome-sharded annotation suite per query. In sharded annotation-input mode, `ldscore` expects one query annotation file per chromosome and reads the separate queries from its columns. This input rule is distinct from the numbered genome-wide query LD-score output files described above.
- The generated annotation header includes `CM` for compatibility, with values written as `NA`. Input annotation CM values are discarded; LD scoring uses reference-panel CM. See [`_annotation_parsing.normalize_annotation_chunk`](../../src/ldsc/_annotation_parsing.py).
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

Regression accepts legacy `.sumstats` and `.sumstats.gz` files containing `SNP`, `A1`, `A2`, `Z`, and `N` directly through `--sumstats-file` or `--sumstats-sources`; remunging is optional. Legacy LD-score suites must first pass through `ldsc convert-ldsc2-ldscores`. That converter accepts complete chromosome 1–22 unpartitioned or baseline-only suites with the required counts and supporting inputs. Regression then reads the converted canonical directory.

See the [sumstats compatibility contract](../current/legacy-sumstats-compatibility.md) and [conversion guide](utility-functionalities/convert-ldsc2-ldscores.md). Dedicated [sumstats tests](../../tests/test_legacy_sumstats_compatibility.py) and [conversion tests](../../tests/test_legacy_ldscore_converter.py) cover these routes; this coverage does not imply equivalence for every possible dataset.

## Further workflow references

- [LD-score calculation](main-functionalities/ldscore.md), [heritability](main-functionalities/h2.md), and [partitioned heritability](main-functionalities/partitioned-h2.md).
- [BED preparation](utility-functionalities/how-to-customize-your-bed-files.md), [reusable annotations](utility-functionalities/annotate.md), and [resource navigation](utility-functionalities/resources.md).
- [Global configuration](main-functionalities/global-config.md) and [interface changes](main-functionalities/changes-from-LDSC2.md).
