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
| `A1` | `ALLELE1`, `ALLELE_1`, `EFFECT_ALLELE`, `REFERENCE_ALLELE`, `INC_ALLELE`, `EA` | alleles |
| `A2` | `ALLELE2`, `ALLELE_2`, `OTHER_ALLELE`, `NON_EFFECT_ALLELE`, `DEC_ALLELE`, `NEA` | alleles |
| `N` | `WEIGHT` | sumstats |
| `N_CAS` | `NCASE`, `N_CASE`, `N_CASES`, `NCAS`, `CASES_N` | sumstats |
| `N_CON` | `NCONTROL`, `N_CONTROL`, `N_CONTROLS`, `NCON`, `CONTROLS_N` | sumstats |
| `NSTUDY` | `N_STUDY`, `NSTUDIES`, `N_STUDIES` | sumstats |
| `Z` | `ZSCORE`, `GC_ZSCORE` | sumstats |
| `P` | `PVALUE`, `P_VALUE`, `PVAL`, `P_VAL`, `GC_PVALUE` | sumstats |
| `BETA` | `B`, `EFFECT`, `EFFECTS` | sumstats |
| `OR` | *(none)* | sumstats |
| `LOG_ODDS` | *(none)* | sumstats |
| `INFO` | `IMPINFO` | sumstats |
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

### On-disk dtypes (parquet)

LD-score and pairwise R² float columns are narrowed to `float32` immediately
before parquet writes. This halves the storage footprint for those dense float
matrices with no meaningful precision loss. Munged sumstats Parquet preserves
the numeric precision produced by the munger because downstream regression
uses those values as primary trait data. Text-based formats (`.annot.gz` and
legacy `.sumstats.gz`) are unaffected; their size is governed by the number of
decimal digits printed, not numpy dtype.

| Parquet artifact | Columns cast to `float32` on write | Columns kept as-is |
| --- | --- | --- |
| `sumstats.parquet` | none | `CHR` (str), `POS` (int64 when complete), `SNP` (str), alleles (str), `Z`/`N`/`FRQ` (numeric precision preserved) |
| `ldscore.baseline.parquet`, `ldscore.query.parquet` | `regression_ld_scores`, all LD-score / annotation columns | `CHR` (str), `POS` (int64), `SNP` (str) |
| Pairwise R² parquet | `R2` | `CHR` (str), `POS_1`, `POS_2` (int64), `SNP_1`, `SNP_2` (str) |

Pairwise R² parquet schema metadata also stores `ldsc:sorted_by_build`,
`ldsc:row_group_size`, `ldsc:n_samples`, `ldsc:r2_bias`, and the minimal
identity provenance keys `ldsc:schema_version`, `ldsc:artifact_type`,
`ldsc:snp_identifier`, and `ldsc:genome_build`. The R2 bias keys let downstream
readers distinguish package-built unbiased R2 from raw sample R2 without
requiring users to repeat sample-size arguments. Old package-written artifacts
without current identity provenance must be regenerated with the current
package.

Reference-panel metadata sidecars (`chr{chrom}_meta.tsv.gz`) remain TSV files,
but current package-written sidecars start with comment lines carrying the same
minimal identity provenance as `ldsc:schema_version`, `ldsc:artifact_type`,
`ldsc:snp_identifier`, and `ldsc:genome_build`.

`N`, `N_CAS`, `N_CON`, and `NSTUDY` remain `float64` in memory and are **not**
subject to float32 narrowing. Curated `sumstats.parquet` keeps `N` as numeric
data; legacy `.sumstats.gz` writes the selected curated columns as text with the
configured float format.

The casting is applied by a narrow helper (`_cast_parquet_floats`) called at the
write site in `outputs.py`. It selects all `float64` columns and recasts them;
no column names need to be listed explicitly. LD-score output parquet files are
written with one row group per chromosome and record that layout in
`manifest.json`; sumstats Parquet is sorted by `CHR`, `POS`, and original row
order, and records its row groups in `sumstats.log`.

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

### `A1` and `A2`

`A1` and `A2` remain the canonical written LDSC column names for compatibility.
`A1` is the allele that the signed statistic is relative to; in most GWAS
formats this is the effect or increasing allele. `A2` is the counterpart
allele. Positive `Z`, positive `BETA`, positive `LOG_ODDS`, and `OR > 1` mean
the effect is increasing with respect to `A1`.

Aliases such as `REFERENCE_ALLELE` are interpreted in this signed-statistic
sense only. They do not claim that `A1` is the genome reference allele. Genome
`REF`/`ALT` columns are therefore not global aliases; they are applied only by
format-aware inference when the workflow can make the semantics explicit.

### SNP identity modes

The public SNP identifier modes are exactly `rsid`, `rsid_allele_aware`,
`chr_pos`, and `chr_pos_allele_aware`. The default is
`chr_pos_allele_aware`. These mode names are exact; aliases such as `RSID`,
`BP`, or `POSITION` are input-column aliases only.

Base modes are fully allele-blind. `rsid` uses only `SNP`; `chr_pos` uses only
`CHR:POS`. If `A1/A2` columns are present in base modes, they may be preserved
for output or signed-statistic orientation, but they do not affect identity,
duplicate filtering, retention, or drop reasons.

Allele-aware modes use `A1/A2` only to build safer merge keys. The unordered,
strand-aware allele set participates in `SNP:<allele_set>` or
`CHR:POS:<allele_set>` identity. Allele-aware artifact cleanup drops rows with
missing alleles, invalid or non-SNP alleles, identical pairs, strand-ambiguous
pairs, package-wide multi-allelic base-key clusters, and duplicate effective
merge-key clusters. The duplicate policy is drop-all after computing the
effective merge key for the active mode.

Restriction files may omit alleles. Allele-free restrictions match by base key
and can keep multiple candidate rows until later artifact cleanup. Packaged HM3
restrictions are allele-free base-key filters. Allele-bearing restrictions in
allele-aware modes match by the effective allele-aware key. Annotation files may
omit alleles even in allele-aware modes because they describe genomic
membership; if annotation files include alleles, those alleles participate in
allele-aware matching.

### `CHR` format
`normalize_chromosome` (in `chromosome_inference.py`) strips the `chr` prefix on
read. All Python-written artifacts store bare labels: `"1"`, `"22"`, `"X"`,
`"MT"`. External files (e.g., the hm3 curated map with `chr` prefix) are
normalized on read; their on-disk format is not changed.

---

## 3. Column Ordering in Written Artifacts

**Rule:** SNP metadata columns appear in a stable workflow-specific order.
Annotation artifacts use **(CHR, POS, SNP, ...)**. Munged sumstats use
**(SNP, CHR, POS, ...)** to keep the legacy leading `SNP` convention while
adding coordinates. Pairwise R2 uses the documented pairwise coordinate
variant. Public LD-score result tables currently use **(CHR, SNP, POS, ...)** to
match `LDScoreResult` and regression-loader validation.

Column order does not define row identity. In `rsid`-family modes, `SNP`
supplies the base identifier. In `chr_pos`-family modes, `SNP` is only a label
and matching uses `CHR` and `POS`; munger liftover updates those coordinate
columns without rewriting `SNP`. Liftover-specific duplicate coordinate
handling is therefore based on `CHR/POS` groups, not rsID labels.

This applies to artifacts **written by this package**. External input files
(user-supplied annotation files, the hm3 curated map, PLINK BIM files, etc.) are
normalized on read regardless of their on-disk column order; the input order is
not preserved in any output.

For annotation output files specifically: the kernel always reconstructs the
column layout from its internal representation rather than passing through the
input file unchanged. This means the leading metadata columns (`CHR, POS, SNP,
CM`) follow the annotation rule unconditionally. The annotation-specific columns that follow
`CM` retain their input order (the kernel does not reorder them).

| Artifact | Leading columns | Remaining columns |
|----------|----------------|-------------------|
| Annotation (`.annot.gz`) | `CHR, POS, SNP, CM` | annotation columns (input order preserved) |
| LD-score output (`ldscore.baseline.parquet`) | `CHR, SNP, POS, regression_ld_scores` | baseline LD-score columns |
| LD-score output (`ldscore.query.parquet`) | `CHR, SNP, POS` | query LD-score columns |
| Munged sumstats (`sumstats.parquet` or `.sumstats.gz`) | `SNP, CHR, POS, A1, A2` | `Z, N, FRQ` |
| Canonical pairwise R² parquet | `CHR, POS_1, POS_2, SNP_1, SNP_2, A1_1, A2_1, A1_2, A2_2, R2` when endpoint alleles are available | endpoint allele columns are required in allele-aware modes and ignored for identity in base modes |
| Dropped-SNP audit sidecar (`dropped_snps/*.tsv.gz`) | `CHR, SNP, source_pos, target_pos, reason, base_key, identity_key, allele_set, stage` | |

Dropped-SNP audit sidecars are always written by liftover-aware public
workflows, including header-only files for clean processed outputs. They use
nullable read dtypes: `CHR`, `SNP`, and `reason` as strings, and `source_pos`
and `target_pos` as nullable integer positions. Identity fields are nullable
when a drop reason has no effective identity key. Sumstats and reference-panel
sidecars may include identity reasons (`missing_allele`, `invalid_allele`,
`strand_ambiguous_allele`, `multi_allelic_base_key`, `duplicate_identity`) and
liftover reasons (`missing_coordinate`, `source_duplicate`,
`unmapped_liftover`, `cross_chromosome_liftover`, and `target_collision`);
reference-panel sidecars do not emit `missing_coordinate`.

---

## 4. Implementation Locations

The following locations must be kept consistent with this document.

| Location | What to keep in sync |
|----------|---------------------|
| `src/ldsc/column_inference.py` — `INTERNAL_SUMSTATS_ARTIFACT_SPECS` | Tuple order: `CHR, POS, SNP, ...` |
| `src/ldsc/column_inference.py` — `INTERNAL_LDSCORE_ARTIFACT_SPECS` | Strict reload specs still include optional metadata fields used by internal/legacy readers: `CHR, POS, SNP, CM, MAF` |
| `src/ldsc/column_inference.py` — `INTERNAL_ANNOT_ARTIFACT_SPECS` | Already correct: `CHR, POS, SNP, CM, MAF` |
| `src/ldsc/_kernel/ldscore.py` — `ANNOT_META_COLUMNS` | `("CHR", "POS", "SNP", "CM", "MAF")` |
| `src/ldsc/ldscore_calculator.py` — `_split_ldscore_table` | Public output order: baseline `CHR, SNP, POS, regression_ld_scores, ...`; query `CHR, SNP, POS, ...` |
| `src/ldsc/outputs.py` — `_write_chromosome_aligned_parquet` | One row group per `CHR` in `ldscore.baseline.parquet` and `ldscore.query.parquet` |
| `src/ldsc/outputs.py` — `required_baseline` / `required_query` validation lists | Baseline requires `CHR`, `POS`, `SNP`, `regression_ld_scores`, and baseline columns; query requires `CHR`, `POS`, `SNP`, and query columns |

`regression_ld_scores` is the historical `w_ld`: an LD score computed over the
regression SNP universe. It is stored with the baseline LD-score table because
both h2 and rg use it as one component when the regression kernel computes
final model-dependent weights.
