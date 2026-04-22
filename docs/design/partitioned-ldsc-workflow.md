# Partitioned LDSC Workflow: Technical Reference

This document describes the end-to-end workflow for computing partitioned LD scores and running
partitioned heritability (h²) regression, from raw input files to the final output. It covers
the SNP universe design, the role of each optional argument, and the file formats at each stage.

---

## 1. Overview

Partitioned LD score regression decomposes trait heritability into contributions from
genomic annotations (e.g., a set of enhancer regions defined by BED files). The pipeline
has two sequential phases:

1. **LD score computation** — for each chromosome, compute annotation-weighted LD scores and
   write `.l2.ldscore.gz` / `.M` / `.M_5_50` files.
2. **Regression** — read the LD score files and summary statistics, run the partitioned LDSC
   regression, and report per-annotation heritability enrichment.

---

## 2. Inputs

### 2.1 Baseline annotations

Pre-built `.annot[.gz]` files with one row per SNP (in the reference panel SNP universe)
and one binary column per annotation category. Passed via `--baseline-annot`.

Example: `resources/baseline_v1.2/baseline.@.annot.gz` (one file per chromosome, `@` is
a token expanded to the chromosome number).

### 2.2 Query annotations (BED or pre-built .annot.gz)

The new annotation of interest, either:

- **BED file(s)** (`--query-annot-bed`): raw genomic intervals. The tool projects these onto
  the reference panel SNPs in memory to produce binary annotation columns. No intermediate
  `.annot.gz` files are written.
- **Pre-built `.annot.gz` files** (`--query-annot`): already-projected annotation files.
  Mutually exclusive with `--query-annot-bed`.

### 2.3 Reference panel

Genotype data used to compute LD between SNPs. Either:

- **PLINK format** (`--bfile`): binary `.bed/.bim/.fam` files.
- **Parquet R² format** (`--r2-table`): pre-computed pairwise R² values.

### 2.4 Summary statistics

GWAS summary statistics in `.sumstats.gz` format, used only in the regression phase.
Passed via `--sumstats` on `ldsc h2`, `ldsc partitioned-h2`, or `ldsc rg`.

---

## 3. SNP Universes

The pipeline tracks four distinct SNP sets. Understanding these is essential to interpreting
output counts and row counts.

### Notation

| Symbol | Name | Meaning |
|--------|------|---------|
| **B** | Baseline annotation SNPs | All SNPs with a row in the loaded baseline `.annot.gz` files |
| **A** | Raw reference panel SNPs | All SNPs in the PLINK `.bim` or parquet R² index |
| **A'** | Restricted reference panel | A ∩ `ref_panel_snps_path`; equals A when `--ref-panel-snps` is absent |
| **`ld_reference_snps`** | LD computation universe | B ∩ A' — SNPs in both annotation and (restricted) reference panel |
| **C** | Regression SNP flags | Binary mask from `regression_snps_path`; flags 1 only for SNPs in `ld_reference_snps` |
| **`ld_regression_snps`** | Regression output SNPs | B ∩ A' ∩ C — rows of the output `.l2.ldscore.gz` file |

### Key invariant

`ld_regression_snps ⊆ ld_reference_snps ⊆ B`

- **LD scores** are computed over `ld_reference_snps` (the full reference universe). This
  is the SNP set used for column sums → **M** and **M_5_50** counts.
- The **output file rows** are only `ld_regression_snps`. When `--regression-snps` is absent,
  `ld_regression_snps = ld_reference_snps` and all reference SNPs appear as rows.
- `ld_reference_snps` is an internal variable. It is **never written as a separate file** and
  is discarded after `.M` / `.M_5_50` are computed.

---

## 4. Optional Arguments: `--ref-panel-snps` and `--regression-snps`

Both flags are valid only on `ldsc ldscore`.

### 4.1 `--ref-panel-snps` → restricts the reference panel (A → A')

```
ldsc ldscore --ref-panel-snps hm3.snplist ...
```

Provides a SNP list file (one rsID per line, or CHR:POS depending on `--snp-identifier` mode).
The `RefPanel` applies this restriction internally at `load_metadata()` time (in
`RefPanel._apply_snp_restriction()`). The restricted set A' is stored on `RefPanelSpec`; no
separate filtering step is visible to the caller.

**Effect:** `ld_reference_snps = B ∩ A'` (smaller than `B ∩ A`). `.M` and `.M_5_50` count
only SNPs in this smaller universe. Improves regression conditioning when the reference panel
contains low-quality or imputed SNPs not present in the annotation.

**When absent:** `A' = A`; all reference panel SNPs are eligible.

### 4.2 `--regression-snps` → restricts the output rows (ld_reference_snps → ld_regression_snps)

```
ldsc ldscore --regression-snps hm3.snplist ...
```

Provides a second SNP list that flags which SNPs from `ld_reference_snps` should appear in the
output `.l2.ldscore.gz` rows. The mask is applied **after** LD computation, so LD scores for
all `ld_reference_snps` are computed but only `ld_regression_snps` rows are written.

**Effect:** rows of `.l2.ldscore.gz` = `ld_regression_snps` ⊊ `ld_reference_snps`. The
`regr_weight` column (see Section 6) covers only these rows.

**When absent:** `ld_regression_snps = ld_reference_snps`; all reference SNPs are output.

### 4.3 Typical usage

```bash
# Common case: restrict both reference panel and regression SNPs to HapMap3
ldsc ldscore \
  --baseline-annot resources/baseline_v1.2/baseline.@ \
  --query-annot-bed my_peaks.bed \
  --bfile resources/1kg/1KG_EUR_Phase3_plink/1KG_EUR_Phase3_chr@ \
  --ref-panel-snps resources/w_hm3.snplist \
  --regression-snps resources/w_hm3.snplist \
  --out my_study
```

With both flags pointing to the same HapMap3 list:
- `A'` = HapMap3 SNPs in the reference panel
- `ld_reference_snps` = baseline annotation SNPs ∩ HapMap3 ref panel SNPs
- `ld_regression_snps` = same (since `--regression-snps` = `--ref-panel-snps` here)

---

## 5. Annotation Bundle

Before LD score computation begins, all annotation sources (baseline + query) are loaded and
intersected into a single `AnnotationBundle` by `AnnotationBuilder.run()`. The bundle holds:

- `baseline_matrix`: DataFrame with rows = B, columns = baseline annotation categories
- `query_matrix`: DataFrame with rows = B, columns = query annotation categories (empty if
  no query annotation provided)
- `metadata`: per-SNP CHR/BP/SNP index aligned to the combined matrix rows

The bundle is the **unified SNP universe** for LD computation. Its row set defines **B**.

When `--query-annot-bed` is used, BED intervals are projected onto B in memory using
`pybedtools` overlap. No `.annot.gz` files are written unless the caller explicitly invokes
`run_bed_to_annot(output_dir=...)` via the Python API.

---

## 6. LD Score Computation (per chromosome)

`LDScoreCalculator.run()` iterates over chromosomes. For each chromosome:

1. `ref_panel.load_metadata(chrom)` returns metadata for `A'` (restriction already applied).
2. `ld_reference_snps = B_chrom ∩ A'_chrom` (frozenset intersection).
3. LD scores are computed for all SNPs in `ld_reference_snps`, weighted by annotation columns.
4. `.M` and `.M_5_50` column sums are accumulated over `ld_reference_snps`.
5. The `regression_snps` set C is applied as a binary mask to select `ld_regression_snps`.
6. Output row = `ld_regression_snps`. A `regr_weight` column is appended (the total LD score
   used as regression weight, previously the `.w.l2.ldscore.gz` value).
7. The single merged file `<out_prefix>.<chrom>.l2.ldscore.gz` is written.

`ld_reference_snps` is an internal variable and is not written to disk.

---

## 7. Output Files

### 7.1 `.l2.ldscore.gz` (per chromosome)

| Column | Description |
|--------|-------------|
| CHR | Chromosome |
| SNP | SNP identifier (rsID or chr:pos depending on `--snp-identifier`) |
| BP | Base position |
| `<annot>_L2` | LD score for each baseline/query annotation column |
| `regr_weight` | Per-SNP regression weight (total LD score); replaces the old `.w.l2.ldscore.gz` |

**Rows:** `ld_regression_snps` = B ∩ A' ∩ C.

### 7.2 `.M` and `.M_5_50` (per annotation, per chromosome)

Tab-separated count files. Each value is the **sum of annotation indicator values over
`ld_reference_snps`** (not `ld_regression_snps`). `M_5_50` restricts the count to SNPs
with MAF between 5% and 50%.

These counts are used as denominators in the heritability partitioning formula and are
tracked over the full reference universe (`ld_reference_snps`), not the regression subset.

### 7.3 No `.w.l2.ldscore.gz` file

The old per-chromosome weight file is no longer written. Its information is now embedded as
the `regr_weight` column in `.l2.ldscore.gz`. Pre-computed weight files from external sources
(e.g., EUR `w_ld_chr/` downloads) remain usable via `--w-ld` on regression subcommands.

---

## 8. Regression Phase

`ldsc partitioned-h2 --ref-ld-chr <prefix>. --sumstats <file> [--w-ld <prefix>]`

### 8.1 Loading LD scores from disk

`load_ldscore_from_files()` rebuilds an `LDScoreResult` from disk artifacts. It detects
the file format automatically:

- **New format (default):** `regr_weight` column present in `.l2.ldscore.gz`. `weight_path`
  argument omitted; weights extracted from the embedded column.
- **Legacy format:** no `regr_weight` column. Pass `weight_path` pointing to a separate
  `.w.l2.ldscore.gz` file (or EUR `w_ld_chr/` directory prefix).

SNP set reconstruction:

```python
ld_reference_snps = frozenset()   # not recoverable from disk; never needed after M/M_5_50 load
ld_regression_snps = frozenset(build_snp_id_series(file_rows_df, snp_identifier))
```

`snp_identifier` is passed from `GlobalConfig` and must match the mode used during LD score
computation (`"rsid"` or `"chr_pos"`).

### 8.2 Regression weights

`_resolve_regression_weights()` applies the following priority:

1. If `--w-ld` is provided, load weights from the specified file (backward compatibility with
   pre-computed EUR weights).
2. If the loaded `.l2.ldscore.gz` has a `regr_weight` column, extract it.
3. Otherwise raise `ValueError` — weights are required.

### 8.3 Regression

The regression runner intersects the loaded LD scores with the GWAS summary statistics
on the shared SNP identifier, then fits the LDSC regression to produce per-annotation
heritability estimates.

---

## 9. Python API Summary

```python
from ldsc import (
    AnnotationBuilder, AnnotationSourceSpec,   # BED → AnnotationBundle
    LDScoreCalculator, LDScoreResult,           # LD score computation
    load_ldscore_from_files,                    # load from disk (resumable)
    RegressionRunner,                           # partitioned h² regression
    run_bed_to_annot,                           # write BED → .annot.gz explicitly
)
from ldsc.config import GlobalConfig, LDScoreConfig
from ldsc._kernel.ref_panel import RefPanelSpec

# 1. Build annotation bundle from BED files
global_config = GlobalConfig(genome_build="hg38", snp_identifier="rsid")
source_spec = AnnotationSourceSpec(
    baseline_annot_paths=("resources/baseline_v1.2/baseline.@.annot.gz",),
    bed_paths=("my_peaks.bed",),
)
bundle = AnnotationBuilder(global_config).run(source_spec)

# 2. Compute LD scores
ref_spec = RefPanelSpec(
    backend="plink",
    bfile_prefix="resources/1kg/1KG_EUR_Phase3_chr@",
    genome_build="hg38",
    ref_panel_snps_path="resources/w_hm3.snplist",   # optional
)
ldscore_config = LDScoreConfig(
    regression_snps_path="resources/w_hm3.snplist",  # optional
)
from ldsc._kernel.ref_panel import PlinkRefPanel
result = LDScoreCalculator().run(
    annotation_bundle=bundle,
    ref_panel=PlinkRefPanel(global_config, ref_spec),
    ldscore_config=ldscore_config,
    global_config=global_config,
    output_spec=...,
)

# 3. Load from disk (after step 2 has already written files)
result = load_ldscore_from_files(
    ldscore_path="my_study.@.l2.ldscore.gz",
    counts_path="my_study.@.l2.M_5_50",
    count_kind="m_5_50",
    snp_identifier="rsid",
    # weight_path omitted: uses embedded regr_weight column
)

# 4. Run regression
runner = RegressionRunner()
h2_result = runner.run_partitioned_h2(
    ldscore_result=result,
    sumstats_path="my_study.sumstats.gz",
    global_config=global_config,
)
```

---

## 10. CLI Reference

```bash
# Step 1: compute LD scores (BED query, in-memory projection)
ldsc ldscore \
  --baseline-annot resources/baseline_v1.2/baseline.@ \
  --query-annot-bed my_peaks.bed \
  --bfile resources/1kg/1KG_EUR_Phase3_chr@ \
  --ref-panel-snps resources/w_hm3.snplist \
  --regression-snps resources/w_hm3.snplist \
  --out my_study

# Step 1 (alternative): pre-built query annotation files
ldsc ldscore \
  --baseline-annot resources/baseline_v1.2/baseline.@ \
  --query-annot my_peaks.@.annot.gz \
  --bfile resources/1kg/1KG_EUR_Phase3_chr@ \
  --out my_study

# Step 1b (optional): explicitly save BED → .annot.gz for reuse
ldsc annotate \
  --bed-file my_peaks.bed \
  --baseline-annot resources/baseline_v1.2/baseline.@ \
  --out my_peaks_annot

# Step 2: partitioned h² regression (new format, no --w-ld needed)
ldsc partitioned-h2 \
  --ref-ld-chr my_study. \
  --sumstats my_gwas.sumstats.gz \
  --out my_study_h2

# Step 2 (backward compat): with separate EUR weight file
ldsc partitioned-h2 \
  --ref-ld-chr my_study. \
  --w-ld resources/eur_w_ld_chr/ \
  --sumstats my_gwas.sumstats.gz \
  --out my_study_h2
```
