# Column Schema: Canonical Names, Data Types, and Ordering

This document is the single source of truth for column conventions across all
Python-written artifacts in this package. It governs `column_inference.py`, all
kernel output routines, and any doc or test that asserts column layout.

---

## 1. Canonical Column Names

All written artifacts use the canonical name. The alias families below are
accepted on **read** only; they are never written.

| Canonical | Input aliases | Context |
|-----------|--------------|---------|
| `CHR` | `#CHROM`, `CHROM`, `CHROMOSOME` | all |
| `POS` | `BP`, `POSITION`, `BASE_PAIR`, `BASEPAIR` | all |
| `SNP` | `RSID`, `RS`, `ID`, `MARKERNAME`, `SNPID`, `MARKER` | all |
| `CM` | `CMBP`, `CENTIMORGAN` | annotation, LD-score |
| `MAF` | `FREQ`, `FREQUENCY` | ref panel, annotation |
| `FRQ` | `EAF`, `MAF`, `FRQ_U` | sumstats only |
| `A1` | `ALLELE1`, `EFFECT_ALLELE`, `INC_ALLELE`, `EA`, `REF` | alleles |
| `A2` | `ALLELE2`, `OTHER_ALLELE`, `DEC_ALLELE`, `NEA`, `ALT` | alleles |
| `N` | `WEIGHT` | sumstats |
| `N_CAS` | `NCASE`, `N_CASE`, `N_CASES`, `NCAS`, `CASES_N` | sumstats |
| `N_CON` | `NCONTROL`, `N_CONTROL`, `N_CONTROLS`, `NCON`, `CONTROLS_N` | sumstats |
| `NSTUDY` | `N_STUDY`, `NSTUDIES`, `N_STUDIES` | sumstats |
| `Z` | `ZSCORE`, `GC_ZSCORE` | sumstats |
| `P` | `PVALUE`, `P_VAL`, `GC_PVALUE` | sumstats |
| `BETA` | `B`, `EFFECT`, `EFFECTS` | sumstats |
| `OR` | *(none)* | sumstats |
| `LOG_ODDS` | *(none)* | sumstats |
| `INFO` | *(none)* | sumstats |
| `R2` | *(none)* | pairwise LD |
| `L2` | *(none)* | LD-score output |

---

## 2. Data Types

### In-memory (canonical) dtypes

These dtypes govern pandas DataFrames during processing and are the source of
truth for any computation or validation.

| Column(s) | dtype | Notes |
|-----------|-------|-------|
| `CHR` | `str` (pandas `object`) | No `chr` prefix: `"1"`, `"22"`, `"X"`, `"MT"` |
| `SNP`, `SNP_1`, `SNP_2` | `str` | |
| `A1`, `A2` | `str` | Always uppercased on read |
| `POS`, `POS_1`, `POS_2` | `np.int64` | Explicit `astype(np.int64)`, not plain `int` |
| `CM` | `float64` | |
| `MAF` | `float64` | Values ∈ [0, 0.5]; always minor-allele frequency |
| `FRQ` | `float64` | Values ∈ [0, 1]; effect-allele frequency from sumstats |
| `N`, `N_CAS`, `N_CON`, `NSTUDY` | `float64` | Non-integer in meta-analysis / effective N |
| `Z`, `P`, `BETA`, `OR`, `LOG_ODDS`, `INFO` | `float64` | |
| `R2` | `float64` | |
| `L2` (and annotation-specific LD columns) | `float64` | |

### On-disk dtypes (parquet only)

Float columns are narrowed to `float32` immediately before any `to_parquet`
write. This halves the storage footprint for float data with no meaningful
precision loss. Text-based formats (`.annot.gz`, `.sumstats.gz`) are unaffected
— their size is governed by the number of decimal digits printed, not numpy
dtype.

| Parquet artifact | Columns cast to `float32` on write | Columns kept as-is |
| --- | --- | --- |
| `baseline.parquet`, `query.parquet` | `CM`, `MAF`, `regr_weight`, all `L2` / annotation columns | `CHR` (str), `POS` (int64), `SNP` (str) |
| Pairwise R² parquet | `R2` | `CHR` (str), `POS_1`, `POS_2` (int64), `SNP_1`, `SNP_2` (str) |

`N`, `N_CAS`, `N_CON`, `NSTUDY` are stored in `.sumstats.gz` (text), never in
parquet. They remain `float64` everywhere and are **not** subject to float32
narrowing.

The casting is applied by a narrow helper (`_cast_parquet_floats`) called at the
write site in `outputs.py`. It selects all `float64` columns and recasts them;
no column names need to be listed explicitly.

### Why `object` dtype for strings?
pandas has no native string dtype in the numpy layer. `str` columns are stored as
Python `str` objects in a numpy pointer array; the column's `.dtype` reports
`object`. This is the standard interoperability convention with numpy/scipy and
requires no special handling.

### Why `float64` for `N`-family columns?
Effective sample sizes from meta-analysis and binary-trait liability-scale
conversions can be non-integer. Using `float64` avoids `NA`-handling issues
with numpy `int64` (which cannot represent `NaN`) and matches legacy LDSC
convention. float32 is unsafe here: it can only represent integers exactly up to
2²⁴ = 16,777,216, and modern mega-GWAS meta-analyses routinely exceed this.

### `MAF` vs `FRQ`
These measure the same quantity but in different contexts:

- **`MAF`** (ref panel, annotation): always ≤ 0.5; represents the minor allele
  in the reference population.
- **`FRQ`** (munged sumstats): effect-allele frequency; may exceed 0.5 depending
  on allele coding. When a raw sumstats file provides a column named `MAF`, the
  munger renames it to `FRQ` on output.

### `CHR` format
`normalize_chromosome` (in `chromosome_inference.py`) strips the `chr` prefix on
read. All Python-written artifacts store bare labels: `"1"`, `"22"`, `"X"`,
`"MT"`. External files (e.g., the hm3 curated map with `chr` prefix) are
normalized on read; their on-disk format is not changed.

---

## 3. Column Ordering in Written Artifacts

**Rule:** SNP metadata columns appear in the order **(CHR, POS, SNP, ...)** in
all written artifacts. CHR and POS are always adjacent and always precede SNP.

This applies to artifacts **written by this package**. External input files
(user-supplied annotation files, the hm3 curated map, PLINK BIM files, etc.) are
normalized on read regardless of their on-disk column order; the input order is
not preserved in any output.

For annotation output files specifically: the kernel always reconstructs the
column layout from its internal representation rather than passing through the
input file unchanged. This means the leading metadata columns (`CHR, POS, SNP,
CM`) follow the rule unconditionally. The annotation-specific columns that follow
`CM` retain their input order (the kernel does not reorder them).

| Artifact | Leading columns | Remaining columns |
|----------|----------------|-------------------|
| Annotation (`.annot.gz`) | `CHR, POS, SNP, CM` | annotation columns (input order preserved) |
| LD-score output (`baseline.parquet`, `query.parquet`) | `CHR, POS, SNP, CM, MAF` | L2 / annotation columns |
| Munged sumstats (`.sumstats.gz`) | `CHR, POS, SNP, A1, A2` | `N, Z, FRQ, ...` |
| Canonical pairwise R² parquet | `CHR, POS_1, POS_2, SNP_1, SNP_2, R2` | |

---

## 4. Implementation Locations

The following locations must be kept consistent with this document.

| Location | What to keep in sync |
|----------|---------------------|
| `src/ldsc/column_inference.py` — `INTERNAL_SUMSTATS_ARTIFACT_SPECS` | Tuple order: `CHR, POS, SNP, ...` |
| `src/ldsc/column_inference.py` — `INTERNAL_LDSCORE_ARTIFACT_SPECS` | Tuple order: `CHR, POS, SNP, CM, MAF` |
| `src/ldsc/column_inference.py` — `INTERNAL_ANNOT_ARTIFACT_SPECS` | Already correct: `CHR, POS, SNP, CM, MAF` |
| `src/ldsc/_kernel/ldscore.py` — `ANNOT_META_COLUMNS` | `("CHR", "POS", "SNP", "CM", "MAF")` |
| `src/ldsc/_kernel/ldscore.py` — explicit reorder list (line ~1904) | `["CHR", "POS", "SNP", "CM", "MAF"]` |
| `src/ldsc/outputs.py` — `required_baseline` / `required_query` validation lists | `["CHR", "POS", "SNP", ...]` |
| Docstrings in `ldscore.py` referencing `CHR, SNP, POS` | Update to `CHR, POS, SNP` |
