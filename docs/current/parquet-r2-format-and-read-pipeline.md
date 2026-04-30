# Parquet R2 Table: Format Specification and Read Pipeline

This document defines the canonical format for parquet-backed pairwise R² reference
panel files and describes how the read pipeline validates, indexes, and converts them
into the dense LD matrices required by the LD score computation algorithm.

---

## 1. Motivation

The parquet R² backend replaces repeated `pyarrow.Dataset.to_table(filter=...)` calls
with a row-group index strategy modelled after tabix: parquet footer statistics
(per-column min/max per row group) are read once at file open time and used to
identify the minimal set of row groups that overlap each genomic query window.
Only those row groups are read from disk per query.

Benchmark on a synthetic 25.6M-pair chromosome parquet:

| Approach | Setup | ms / query | Peak RAM |
|---|---|---|---|
| Current (`Dataset` + `to_pandas`) | 56 s | ~0 ms (setup pays all) | 1 358 MB |
| Pattern 2 fine-grained (50k rows/group) | 0.01 s | 1.2 ms | 3 MB |
| Pattern 2 coarse (3 row groups) | 0.01 s | 27 497 ms | 1 460 MB |

Fine-grained row groups are required for the performance benefit. Coarse files are
still accepted but trigger a startup warning (see §3.2).

---

## 2. Parquet Format Specification

### 2.1 Column Schema

The canonical parquet has exactly these six columns. No additional columns are written
or expected by the reader.

| Column | Type | Description |
|---|---|---|
| `CHR` | `string` | Chromosome label, normalized to bare integer string (`"1"` … `"22"`, `"X"`) |
| `POS_1` | `int64` | Genomic position of the left SNP in the pair (1-based, sort build) |
| `POS_2` | `int64` | Genomic position of the right SNP in the pair (1-based, sort build) |
| `R2` | `float32` | Squared Pearson correlation between allele dosages |
| `SNP_1` | `string` | dbSNP rsID of the left SNP |
| `SNP_2` | `string` | dbSNP rsID of the right SNP |

**Dropped from the legacy schema:** `hg19_pos_1`, `hg38_pos_1`, `hg19_pos_2`,
`hg38_pos_2`, `hg19_Uniq_ID_1/2`, `hg38_Uniq_ID_1/2`, `Dprime`, `+/-corr`.
`Dprime` was always written as NaN and is unused. `+/-corr` is irrelevant to LD
score regression because R² is unsigned. The build-specific position columns are
replaced by the single canonical pair `POS_1`/`POS_2` whose build is recorded in
the file metadata (see §2.4).

#### Alias-tolerant loading

The default on-disk column names are the uppercase canonical fields above:
`CHR`, `POS_1`, `POS_2`, `SNP_1`, `SNP_2`, `R2`. At load time, the parquet reader
resolves aliases through `src/ldsc/column_inference.py` before using the file.
For example, a file with columns `chr`, `bp_1`, `bp_2`, `rsid_1`, `rsid_2`, `R2`
is treated as the same logical schema as `CHR`, `POS_1`, `POS_2`, `SNP_1`,
`SNP_2`, `R2`.

This means the docs define the canonical names used by package-written parquet
artifacts, while the loader remains permissive for external parquet sources.

Writing canonical R2 parquet requires PyArrow. The package writer uses
`pyarrow.parquet.ParquetWriter` so it can set schema metadata and explicit row-group
sizes; pandas/fastparquet fallback writing is intentionally not part of the canonical
path.

### 2.2 Invariants

- `POS_1 < POS_2` for every row (canonical pair orientation, enforced at write time).
- Rows are sorted by `POS_1 ASC` (non-decreasing) across the entire chromosome. Ordering of `POS_2` within rows that share the same `POS_1` is not enforced; row-group pruning depends only on `POS_1` footer statistics.
- Each file contains exactly one chromosome. Multi-chromosome files are not supported.
- All pairs lie within the configured LD window at write time; no pairs outside the
  window are stored.

The writer **asserts** the sort invariant on every incoming pair. If a pair arrives
with `POS_1` less than the previous pair's `POS_1`, `write_r2_parquet`
raises immediately with a clear message:

> *"Pairs must arrive in non-decreasing POS_1 order. Received POS_1={x} after
> POS_1={prev}. Verify that the reference panel builder traverses SNPs in ascending
> positional order."*

This fail-fast behaviour is intentional: a violated sort invariant corrupts
row-group footer statistics, producing wrong window queries and silently incorrect
LD scores with no downstream diagnostic signal.

### 2.3 Row Group Size

Default: **50 000 rows per row group**.

This is tuned to the typical number of LD pairs per 1 Mb query window under HapMap3
SNP density (~50 000 pairs), so that most window queries touch 1–2 row groups.

Tradeoff guide for non-default reference panels:

| Row group size | Pairs/query ≪ size | Pairs/query ≈ size | Pairs/query ≫ size |
|---|---|---|---|
| Too small | Many groups per query; decompression overhead dominates | — | — |
| **Ideal** | — | **1–2 groups per query** | — |
| Too large | — | — | Reads excess data per query |

Recommended range: **25 000 – 200 000** rows per group. Choose closer to the
expected pairs-per-window for the target reference panel density.

### 2.4 Parquet Schema Metadata

Two key-value entries are written into the parquet `schema.metadata` at file creation:

| Key | Type | Description |
|---|---|---|
| `ldsc:sorted_by_build` | string | Genome build used to define `POS_1`/`POS_2`. One of `"hg19"` or `"hg38"`. |
| `ldsc:row_group_size` | string (int) | Intended row group size used when writing. Informational; reader uses actual footer stats. |

These are accessed via `pq.ParquetFile(path).schema_arrow.metadata` and are stored
as UTF-8 byte strings by PyArrow.

### 2.5 Canonical Example

**Schema metadata:**

```
ldsc:sorted_by_build = "hg19"
ldsc:row_group_size  = "50000"
```

**Row group footer statistics (first three groups of a chr1 file):**

| Row group | Rows | POS_1 min | POS_1 max |
|---|---|---|---|
| 0 | 50 000 | 752 566 | 1 823 421 |
| 1 | 50 000 | 1 823 422 | 2 910 884 |
| 2 | 50 000 | 2 910 885 | 4 011 203 |
| … | … | … | … |

**Data rows (first 10 rows of the file):**

| CHR | POS_1 | POS_2 | R2 | SNP_1 | SNP_2 |
|---|---|---|---|---|---|
| 1 | 752 566 | 776 546 | 0.8214 | rs3094315 | rs2905036 |
| 1 | 752 566 | 800 007 | 0.1342 | rs3094315 | rs11240777 |
| 1 | 752 566 | 817 186 | 0.0891 | rs3094315 | rs4040617 |
| 1 | 752 566 | 823 656 | 0.3124 | rs3094315 | rs2980319 |
| 1 | 776 546 | 800 007 | 0.4451 | rs2905036 | rs11240777 |
| 1 | 776 546 | 817 186 | 0.0674 | rs2905036 | rs4040617 |
| 1 | 776 546 | 823 656 | 0.1982 | rs2905036 | rs2980319 |
| 1 | 800 007 | 817 186 | 0.7231 | rs11240777 | rs4040617 |
| 1 | 800 007 | 823 656 | 0.1562 | rs11240777 | rs2980319 |
| 1 | 817 186 | 823 656 | 0.8913 | rs4040617 | rs2980319 |

Note: rows are sorted by non-decreasing `POS_1`. The left SNP (`POS_1 = 752 566`) appears
in multiple consecutive rows — one per right-side neighbor within the LD window.

---

### 2.6 SNP Metadata Sidecar

The R² parquet stores **only pairwise data** (`CHR`, `POS_1`, `POS_2`, `R2`,
`SNP_1`, `SNP_2`). Per-SNP metadata — chromosome, base position, rsID, genetic
distance, and minor allele frequency — is kept in a separate **metadata sidecar**
file. The sidecar is optional, but strongly recommended for production use.

Treat the parquet R² file and the metadata sidecar as a paired artifact when
possible. Together they define one logical reference panel: the parquet provides
pairwise LD entries, while the sidecar provides the complete per-SNP metadata
needed for MAF filtering, common counts, and cM windows.

#### Fallback when the sidecar is absent

R2 parquet files are required. In `ref_panel_dir` mode, chromosome `N` must have
`chrN_r2.parquet`; missing R2 is a hard error. If the matching
`chrN_meta.tsv.gz` sidecar is absent and `fail_on_missing_metadata=False`,
`ParquetR2RefPanel.load_metadata()` synthesizes a minimal metadata table by
scanning the R2 endpoint columns:

- `CHR`, `POS`, and `SNP` are derived from the union of `CHR/POS_1/SNP_1` and
  `CHR/POS_2/SNP_2`.
- `CM` is present but missing for every row.
- `MAF` is unavailable.

This fallback cannot recover SNPs that have no emitted pair rows, cannot apply
`--maf-min`, cannot compute common-SNP counts, and cannot support
`--ld-wind-cm` unless annotation metadata supplies cM values. Setting
`fail_on_missing_metadata=True` turns a missing sidecar into a hard error.

#### Sidecar format

A gzip-compressed or plain whitespace-delimited text file with one row per SNP.
Column names are resolved flexibly via `REFERENCE_METADATA_SPEC_MAP` (defined in
`src/ldsc/column_inference.py`); any recognized synonym is accepted.

| Column | Recognized synonyms | Required when | Description |
| --- | --- | --- | --- |
| `CHR` | `chr`, `chromosome` | `chr_pos` identifier mode | Chromosome label (bare integer or `chrN`) |
| `POS` | `BP`, `pos`, `position` | `chr_pos` identifier mode | 1-based genomic position |
| `SNP` | `rsID`, `rs`, `snpid` | `rsid` identifier mode | dbSNP rsID |
| `CM` | `cM`, `genetic_map_cm` | Never | Genetic map position (cM); filled with `NA` if absent |
| `MAF` | `maf`, `freq`, `AF`, `alt_freqs` | When `fail_on_missing_metadata=True` | Minor allele frequency; folded to [0, 0.5] on load |

In `rsid` mode the `SNP` column is mandatory. In `chr_pos` mode `CHR` and `POS` are
mandatory. A file providing all five columns works in either mode and is recommended
for production use.

#### Role in the downstream workflow

When present, the metadata sidecar feeds directly into three steps of the partitioned-LDSC workflow
(see `docs/design/partitioned-ldsc-workflow.md`, §6):

1. **Reference panel universe (A').** `load_metadata()` reads the sidecar, optionally
   applies `_apply_snp_restriction()` when `ref_panel_snps_file` is set, and returns
   the restricted per-SNP table A'. This is what the annotation bundle is aligned
   against (`B_chrom ∩ A'_chrom`) to materialize `ld_reference_snps`.

2. **All and common count vectors.** The `MAF` column propagates through the
   metadata DataFrame and is consumed by `compute_counts()` to compute
   `common_reference_snp_count` using `MAF >= common_maf_min`.
   If `MAF` is absent and `fail_on_missing_metadata` is `True`, `_read_metadata_table`
   raises `ValueError`; if `False`, common counts cannot be computed and regression
   falls back to all counts.

3. **Genetic-map window (`ld_wind_cm`).** For LD-score calculation, the annotation
   file's `CM` is the first source. The sidecar metadata only fills missing `CM`
   values. The final merged `CM` column is required when `--ld-wind-cm` is
   specified. If `CM` is absent from both the annotation metadata and sidecar, the
   pipeline hard-fails at window construction time. Provide annotation metadata or
   a sidecar that includes the `CM` column whenever using `--ld-wind-cm`.

#### Relationship to the R² parquet

`SortedR2BlockReader` receives the loaded sidecar table as a `metadata: pd.DataFrame`
argument at construction time. It uses this table to build `index_map` — the mapping
from SNP identifier to matrix row/column index — before any window query is issued.
The parquet's own identifier columns (`SNP_1`, `SNP_2`, `POS_1`, `POS_2`) serve only
as lookup keys at query time; all per-SNP analysis metadata comes from the sidecar.

This two-file design separates concerns: the parquet is a dense sorted pair table
optimised for row-group pruning; the sidecar is a lightweight per-SNP lookup table
that drives SNP universe construction and MAF/CM metadata.

In the workflow notation of `docs/design/partitioned-ldsc-workflow.md`, the raw
reference-panel SNP universe `A` comes from this paired artifact, but its authoritative
materialization is the sidecar SNP table. The parquet pair rows are not scanned at
runtime to define or validate `A`; they are used only to answer LD window queries once
`ld_reference_snps = B ∩ A'` has been established from annotations plus sidecar metadata.

---

## 3. Read Pipeline

The read pipeline is implemented in `SortedR2BlockReader`
(`src/ldsc/_kernel/ldscore.py`). The class is the single point of contact between
the LD score computation algorithm and the parquet file.

### 3.1 File Open and Schema Validation

At `__init__`, the reader:

1. Opens the file as `pq.ParquetFile(path)` (not `pyarrow.Dataset`).
2. Confirms that the six required columns are present in the schema; raises
   `ValueError` with a clear message if any are missing.
3. Reads `schema.metadata` and proceeds through the **3-tier build validation**
   described in §3.2.

### 3.2 Build Validation (3-tier)

The parquet is sorted by one genome build's positions. The query build must match
that sort build for row-group pruning to be correct.

| Tier | Condition | Action |
|---|---|---|
| 1 | `ldsc:sorted_by_build` present **and matches** query build | Proceed silently |
| 2 | `ldsc:sorted_by_build` present **and mismatches** query build | Raise `ValueError`: *"Parquet sorted for {parquet_build} but analysis uses {query_build}. Use the matching `{query_build}/r2` reference file or regenerate the panel with the correct source/liftover inputs."* |
| 3 | `ldsc:sorted_by_build` absent or `None` | Infer build from first row group (see below); log `WARNING`; proceed if inferred build matches query build, raise `ValueError` if not |

**Tier-3 inference:** read only the first row group, extract `(CHR, POS_1)` as a
`(CHR, POS)` table, and call the existing `infer_chr_pos_build()` function
(`src/ldsc/genome_build_inference.py`). This reads ~1 MB from disk and requires
≥ 200 matched HM3 reference SNPs. The WARNING message is:
*"No build metadata found in {path}; inferred {build} from first row group. To
silence this warning, regenerate the parquet with `ldsc build-ref-panel`."*

**Coarse row group warning:** after build validation, compute
`avg_rows_per_rg = metadata.num_rows / metadata.num_row_groups`. If
`avg_rows_per_rg > 500 000`, emit a `WARNING`:
*"{path} has {n} row group(s) (avg {avg:,.0f} rows/group). Query performance will
be severely degraded. Regenerate with `row_group_size=50000` for optimal speed."*
The file is still used; coarse row groups reduce pruning effectiveness but do not
affect numerical correctness.

### 3.3 Row-Group Index Construction

After validation, the reader builds an in-memory index from the parquet footer:

```python
col_idx = schema.get_field_index("POS_1")
self._rg_bounds = [
    (rg.column(col_idx).statistics.min,
     rg.column(col_idx).statistics.max,
     i)
    for i in range(metadata.num_row_groups)
    if rg.column(col_idx).statistics.has_min_max
]
```

This reads only the footer (a few KB) — no pair data is loaded at init time.
Row groups without valid `POS_1` min/max footer statistics are excluded from
`_rg_bounds` rather than treated as fatal. They are therefore not used by the
pruned canonical query path; regenerate such files with PyArrow statistics enabled
if those rows must be queryable.

### 3.4 Window Queries (`_query_union_rows`)

For each sliding-window block, `cross_block_matrix` and `within_block_matrix` call
`_query_union_rows(pos_min, pos_max)`, which:

1. Identifies overlapping row groups from `_rg_bounds`:
   ```python
   rg_idxs = [i for mn, mx, i in self._rg_bounds if mn <= pos_max and mx >= pos_min]
   ```
2. Reads only those row groups via:
   ```python
   tbl = self._pf.read_row_groups(rg_idxs, columns=["POS_1", "POS_2", "R2"])
   ```
3. Converts to numpy arrays using zero-copy `.to_numpy()` — no Python object
   creation, no intermediate pandas DataFrame.
4. Applies a positional mask: `(POS_1 >= pos_min) & (POS_2 <= pos_max)`.
5. Maps `POS_1`/`POS_2` (or `SNP_1`/`SNP_2`) to matrix indices via `self.index_map`.
6. Applies optional R² bias correction (`_transform_r2`).
7. Returns `(i_array, j_array, r2_array)` as numpy arrays.

**Coarse-file behaviour:** if the file has few, large row groups, step 1 may return
all row groups for every query — performance degrades to O(N) per query but
correctness is maintained.

### 3.5 Dense Matrix Construction

`cross_block_matrix` and `within_block_matrix` consume the output of
`_query_union_rows` and fill a dense `float32` numpy matrix via vectorized
fancy-index assignment:

```python
matrix[i_array, j_array] = r2_array   # all pairs in one numpy write
matrix[j_array, i_array] = r2_array   # symmetric fill (within-block only)
np.fill_diagonal(matrix, 1.0)
```

These matrices are the direct input to `ld_score_var_blocks_from_r2_reader`, which
accumulates annotation-weighted LD scores via the sliding-block algorithm.

---

## 4. Caveats and Constraints

**Build specificity.** The parquet encodes positions from one genome build only.
Using the same file for analyses under a different build produces silently incorrect
row-group pruning (wrong genomic windows are queried). The 3-tier validation in §3.2
is the safeguard. Users who need both hg19 and hg38 must maintain two separate
parquet files, one sorted per build.

**One chromosome per file.** The reader is instantiated per chromosome and expects
a single-chromosome file. Multi-chromosome files are not tested and may produce
incorrect results if the `CHR` filter is not applied at the row-group level.

**Coarse row groups.** Files with few large row groups (e.g., the legacy 2-group
files) are accepted but may degrade query performance by orders of magnitude
(benchmark: 27 497 ms/query with 3 row groups vs. 1.2 ms/query with 514 groups).
A startup warning is emitted; users should regenerate for production use.

**`POS_1 < POS_2` invariant.** The reader assumes canonical pair orientation is
enforced at write time by the sort assertion described in §2.2. Pairs where
`POS_1 >= POS_2` are not re-validated at read time but will cause incorrect matrix
index lookups if present.

**Identifier mode.** `SNP_1`/`SNP_2` are always written. `identifier_mode` is a
runtime parameter; the reader selects which columns to use for index mapping
(`POS_1`/`POS_2` for `chr_pos` mode, `SNP_1`/`SNP_2` for `rsid` mode) without
requiring different parquet files.

**Legacy raw-schema parquets.** Parquets in the old schema (containing
`hg19_pos_1`, `hg38_pos_1`, etc. but no canonical `POS_1`/`POS_2` columns) are
detected at `__init__` by the absence of `POS_1` and `POS_2` from the schema. The
reader emits a `WARNING` and falls back to the existing `pyarrow.Dataset` full-scan
path:

> *"{path} uses the legacy raw schema. Row-group pruning is disabled and query
> performance will be severely degraded. Convert to the canonical format with
> `ldsc normalize-r2-parquet --input {path} --genome-build hg19 --output {new_path}`
> to restore full performance."*

The full-scan fallback is numerically identical to the current implementation;
correctness is preserved. This compatibility path is planned for removal in a future
major version once the ecosystem has migrated to the canonical schema.

---

## 5. Affected Modules

| Module | Change |
|---|---|
| `_kernel/ref_panel_builder.py` | `write_r2_parquet`: require PyArrow; assert sort invariant on each incoming pair; write new 6-column schema; preserve canonical dtypes for empty outputs; write `ldsc:sorted_by_build` and `ldsc:row_group_size` metadata; default `row_group_size=50_000` |
| `_kernel/ldscore.py` — `SortedR2BlockReader.__init__` | Detect schema (canonical vs legacy raw); canonical path: open as `pq.ParquetFile`, build `_rg_bounds` index from footer, validate build (3-tier), warn if coarse row groups; legacy path: open as `pyarrow.Dataset`, emit deprecation warning, proceed with full-scan fallback |
| `_kernel/ldscore.py` — `SortedR2BlockReader._query_union_rows` | Canonical path: row-group index lookup + `read_row_groups` + `.to_numpy()`; legacy path: existing `dataset.to_table(filter=...)` behaviour unchanged |
| `_kernel/ldscore.py` — `read_sorted_r2_presence` | Remove standalone function; the parquet backend no longer performs a runtime SNP-presence scan |
| `_kernel/ldscore.py` — `compute_chrom_from_parquet` | Use the sidecar-defined `A'` directly when forming `ld_reference_snps`; no parquet presence scan is performed |
| Tests — canonical schema | Add fixtures writing new 6-column sorted parquets; assert row-group index is built; assert window queries return correct pairs |
| Tests — legacy schema | Retain existing raw-schema fixtures; assert deprecation warning is emitted; assert numerical output is unchanged |
