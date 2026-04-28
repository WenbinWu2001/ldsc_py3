# Implementation Verification Checklist: Parquet R² Row-Group Pruning

This checklist verifies that the implementation matches the design spec in
`parquet-r2-format-and-read-pipeline.md` and `partitioned-ldsc-workflow.md`.
Items are organized by concern, not by task. Each item is a concrete behavioral
claim that can be confirmed by reading code or running a test.

---

## 1. Writer — Canonical 6-Column Schema

- [ ] Output parquet contains exactly six columns: `CHR`, `POS_1`, `POS_2`, `R2`, `SNP_1`, `SNP_2`. No additional columns written.
- [ ] Legacy columns (`hg19_pos_1`, `hg38_pos_1`, `hg19_pos_2`, `hg38_pos_2`, `hg19_Uniq_ID_1/2`, `hg38_Uniq_ID_1/2`, `Dprime`, `+/-corr`) are absent.
- [ ] `POS_1` and `POS_2` are populated from the column corresponding to the supplied `genome_build` parameter (e.g., `hg19_pos` when `genome_build="hg19"`).
- [ ] `R2` column is written as `float32`.
- [ ] `POS_1` and `POS_2` columns are written as `int64`.
- [ ] `CHR`, `SNP_1`, `SNP_2` columns are written as string.
- [ ] Empty LD outputs preserve the same canonical dtypes rather than falling back to object-typed empty pandas columns.
- [ ] Default `row_group_size` is `50_000` when the caller does not supply one.
- [ ] Canonical LD parquet writing requires PyArrow; no pandas/fastparquet fallback is used for the LD pair table because schema metadata and row-group sizing are part of the format contract.

## 2. Writer — Sort Invariant

- [ ] `POS_1 < POS_2` holds for every row (canonical pair orientation enforced at write time).
- [ ] Rows are written in non-decreasing `POS_1` order.
- [ ] If a pair arrives with `POS_1` less than the previous pair's `POS_1`, `write_ld_parquet` raises `ValueError` immediately.
- [ ] The `ValueError` message includes both the offending `POS_1` value and the previous `POS_1` value.

## 3. Writer — Parquet Schema Metadata

- [ ] `ldsc:sorted_by_build` is written into `schema.metadata` as a UTF-8 byte string matching the `genome_build` argument.
- [ ] `ldsc:row_group_size` is written into `schema.metadata` as a UTF-8 decimal string matching the row group size used.
- [ ] Both keys are readable via `pq.ParquetFile(path).schema_arrow.metadata` without loading any data rows.

## 4. Reader — Schema Detection at Open Time

- [ ] A file whose schema contains `POS_1` / `POS_2` (or recognized aliases) is classified as the canonical path.
- [ ] A file whose schema contains the legacy raw columns (e.g., `hg19_pos_1`) but not `POS_1` / `POS_2` is classified as the raw/legacy path.
- [ ] A file whose schema matches neither pattern raises `ValueError` with a message identifying the columns found.

## 5. Reader Init — Canonical Path

- [ ] Canonical files are opened with `pq.ParquetFile`, not `pyarrow.Dataset`.
- [ ] The reader validates that all six required columns (or their aliases) are present; missing columns raise `ValueError`.
- [ ] If more than one path is supplied for a single chromosome in canonical mode, `ValueError` is raised immediately. The message includes the number of paths received and the chromosome label.
- [ ] `_runtime_layout` is set to `"canonical"` on the reader instance.
- [ ] `_pf` holds the open `pq.ParquetFile` handle after init.
- [ ] No pair data is read from disk during init (only the footer is accessed to build `_rg_bounds`).

## 6. Reader Init — 3-Tier Build Validation (Canonical Path Only)

- [ ] **Tier 1:** `ldsc:sorted_by_build` present in schema metadata and matches the query `genome_build` → init proceeds silently, no warning emitted.
- [ ] **Tier 2:** `ldsc:sorted_by_build` present but does not match the query `genome_build` → `ValueError` raised. Message includes both the parquet's declared build and the analysis build.
- [ ] **Tier 3:** `ldsc:sorted_by_build` absent or `None` → build is inferred from the first row group's `(CHR, POS_1)` data; a `WARNING` is logged; if the inferred build matches the query build, init continues; if it does not match, `ValueError` is raised.
- [ ] Tier 3 inference reads only the first row group (not the entire file).

## 7. Reader Init — Coarse Row Group Warning

- [ ] If `num_rows / num_row_groups > 500_000`, a `WARNING` is logged after build validation.
- [ ] The warning includes the number of row groups and the average rows-per-group.
- [ ] The file is still used; a coarse row group count does not cause a failure.

## 8. Reader Init — Row-Group Bounds Index

- [ ] `_rg_bounds` is a list of `(min_pos1, max_pos1, rg_index)` tuples, one per row group that has valid min/max statistics for `POS_1`.
- [ ] Row groups without min/max footer statistics are excluded from `_rg_bounds`.
- [ ] `min_pos1 <= max_pos1` for every entry.
- [ ] `rg_index` is the zero-based integer index of the row group within the file.
- [ ] `_rg_bounds` is built from footer statistics only; no pair data rows are read.

## 9. Reader Init — Legacy Raw Path

- [ ] Raw-schema files are opened as `pyarrow.Dataset` (existing behavior, unchanged).
- [ ] A `WARNING` is logged at init time for raw-schema files, including the recommendation to convert via `ldsc normalize-r2-parquet` (or `ldsc build-ref-panel`).
- [ ] `_runtime_layout` is set to `"raw"` (or equivalent) so downstream query methods branch correctly.
- [ ] Numerical output of the raw path is identical to the pre-refactor behavior.

## 10. Window Query — Row-Group Pruning

- [ ] For a given `(pos_min, pos_max)` window, only row groups satisfying `mn <= pos_max AND mx >= pos_min` are passed to `read_row_groups`.
- [ ] Row groups entirely to the left (`mx < pos_min`) or entirely to the right (`mn > pos_max`) of the window are not read.
- [ ] If no row groups overlap the window, an empty result is returned immediately without reading any data.
- [ ] `read_row_groups` is called with column projection limited to the identifier and `R2` columns; no unnecessary columns are read.

## 11. Window Query — Fine-Grained Positional Mask

- [ ] After `read_row_groups`, a numpy boolean mask filters rows to those satisfying `POS_1 >= pos_min AND POS_2 <= pos_max`.
- [ ] Rows from overlapping row groups that fall outside the window bounds are excluded before index mapping.

## 12. Window Query — Identifier Mapping

- [ ] In `rsid` mode, `SNP_1` and `SNP_2` columns are used to look up matrix indices via `index_map`.
- [ ] In `chr_pos` mode, `POS_1` and `POS_2` integer values are used to look up matrix indices via `index_map`.
- [ ] Pairs whose identifiers are not found in `index_map` are silently dropped (not raised as errors).

## 13. Dense Matrix Construction

- [ ] The within-block matrix diagonal is filled with `1.0`.
- [ ] Off-diagonal R² values are assigned symmetrically: `matrix[i, j] = r2` and `matrix[j, i] = r2`.
- [ ] The matrix dtype is `float32`.
- [ ] The cross-block matrix is filled one-directionally (no symmetric fill).

## 14. SNP Universe A — Sidecar Authoritative

- [ ] `A` (the raw reference-panel SNP universe) is defined entirely by the metadata sidecar rows. No parquet pair-row scan is performed at any point in the normal workflow.
- [ ] `A' = A ∩ ref_panel_snps_path` when `--ref-panel-snps-path` is supplied; `A' = A` otherwise.
- [ ] `ld_reference_snps = B ∩ A'`.
- [ ] `.M` and `.M_5_50` are computed over `ld_reference_snps` using MAF from the sidecar.
- [ ] Sidecar SNPs that have no pair rows in the parquet (isolated SNPs) are included in `ld_reference_snps` and receive LD score `1.0` via the diagonal fill.

## 15. Removed Functions

- [ ] `read_sorted_r2_presence` does not exist anywhere in the runtime source (`src/`).
- [ ] `get_present_identifiers` does not exist anywhere in the runtime source (`src/`).
- [ ] `filter_reference_to_present_r2` either does not exist or has no callers in the runtime source.
- [ ] No call site in `compute_chrom_from_parquet` or any other runtime function invokes a parquet identifier scan.
- [ ] A grep for all three names across `src/` returns zero matches.

## 16. Alias-Tolerant Loading

- [ ] The parquet reader accepts lowercase or synonym column names (e.g., `chr`, `bp_1`, `bp_2`, `rsid_1`, `rsid_2`) and resolves them to the canonical logical fields at load time.
- [ ] The sidecar reader accepts recognized synonyms for `CHR`, `POS`, `SNP`, `CM`, `MAF` as listed in `REFERENCE_METADATA_SPEC_MAP`.
- [ ] A file with only alias names (no canonical-case columns) loads without error and produces numerically identical results to a file with canonical names.

## 17. Sidecar Hard-Fail

- [ ] `ParquetR2RefPanel.load_metadata()` raises `ImportError` immediately if no sidecar paths are supplied via `metadata_paths`.

## 18. Numerical Parity — Canonical vs. Legacy

- [ ] For the same reference panel data, the canonical parquet path produces LD scores identical (within float32 precision) to the legacy raw-schema Dataset path.
- [ ] `.M` and `.M_5_50` counts are identical between canonical and legacy paths for the same input data and sidecar.
