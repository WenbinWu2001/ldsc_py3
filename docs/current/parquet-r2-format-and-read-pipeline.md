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
Canonical row groups are decoded into numeric endpoint arrays and held in a
chromosome-local LRU cache sized automatically from the LD window and
`snp_batch_size`, so overlapping sliding windows reuse decoded data instead of
re-reading parquet row groups.

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

Package-written canonical parquet files have six base columns plus endpoint
allele columns when the builder has allele metadata. The reader requires the
six base logical fields in every mode, and requires endpoint alleles in
allele-aware modes. Aliases are accepted at load time for external base-mode
inputs; unrelated extra columns are ignored.

| Column | Type | Description |
|---|---|---|
| `CHR` | `string` | Chromosome label, normalized to bare integer string (`"1"` … `"22"`, `"X"`) |
| `POS_1` | `int64` | Genomic position of the left SNP in the pair (1-based, sort build) |
| `POS_2` | `int64` | Genomic position of the right SNP in the pair (1-based, sort build) |
| `R2` | `float32` | Squared Pearson correlation between allele dosages |
| `SNP_1` | `string` | dbSNP rsID of the left SNP |
| `SNP_2` | `string` | dbSNP rsID of the right SNP |
| `A1_1`, `A2_1` | `string` | left endpoint alleles in reference-panel orientation |
| `A1_2`, `A2_2` | `string` | right endpoint alleles in reference-panel orientation |

**Dropped from the legacy schema:** `hg19_pos_1`, `hg38_pos_1`, `hg19_pos_2`,
`hg38_pos_2`, `hg19_Uniq_ID_1/2`, `hg38_Uniq_ID_1/2`, `Dprime`, `+/-corr`.
`Dprime` was always written as NaN and is unused. `+/-corr` is irrelevant to LD
score regression because R² is unsigned. The build-specific position columns are
replaced by the single canonical pair `POS_1`/`POS_2` whose build is recorded in
the file metadata (see §2.4).

#### Alias-tolerant loading

The default on-disk column names are the uppercase canonical fields above:
`CHR`, `POS_1`, `POS_2`, `SNP_1`, `SNP_2`, `R2`, plus endpoint
`A1_1/A2_1/A1_2/A2_2` columns in allele-aware package-built files. At load time, the parquet reader
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

### 2.4 SNP Identity Modes And External R2 Limits

The public SNP identifier modes are exactly `rsid`, `rsid_allele_aware`,
`chr_pos`, and `chr_pos_allele_aware`; the package default is
`chr_pos_allele_aware`. Base modes are fully allele-blind: `rsid` uses only
`SNP_1/SNP_2`, and `chr_pos` uses only `CHR/POS_1/POS_2`. Endpoint allele
columns may be present in base modes, but they do not affect identity,
duplicate filtering, retention, or drop reasons.

Allele-aware modes require package-built canonical R2 parquet with endpoint
allele columns `A1_1/A2_1/A1_2/A2_2`. These modes use the unordered,
strand-aware endpoint allele sets only to make merge keys safer. External raw
R2 parquet inputs are supported only in `rsid` and `chr_pos`; run
`ldsc build-ref-panel` with the current package to produce allele-aware R2
artifacts.

### 2.5 Parquet Schema Metadata

Package-written parquet files include these technical metadata entries:

| Key | Type | Description |
|---|---|---|
| `ldsc:sorted_by_build` | string | Genome build used to define `POS_1`/`POS_2`. One of `"hg19"` or `"hg38"`. |
| `ldsc:row_group_size` | string (int) | Intended row group size used when writing. Informational; reader uses actual footer stats. |
| `ldsc:n_samples` | string (int) | Number of reference-panel individuals used when computing the stored R2 values. Package-built panels write `geno.n`. |
| `ldsc:r2_bias` | string | Bias state of stored R2 values. Package-built panels currently write `"unbiased"`; `"raw"` is reserved for raw sample R2 values that need downstream correction. |

They also include the minimal identity provenance required by current
reloaders:

| Key | Type | Description |
|---|---|---|
| `ldsc:schema_version` | string (int) | Current identity schema version, `1`. |
| `ldsc:artifact_type` | string | `ref_panel_r2`. |
| `ldsc:snp_identifier` | string | Exact mode used to build the artifact. |
| `ldsc:genome_build` | string | Genome build of the endpoint coordinates. |

These are accessed via `pq.ParquetFile(path).schema_arrow.metadata` and are stored
as UTF-8 byte strings by PyArrow.

Downstream readers inspect this metadata before applying R2-bias defaults. For
package-built panels, users do not need to pass `--r2-bias-mode` or
`--r2-sample-size`: unbiased panels resolve to `r2_bias_mode="unbiased"`, and a
future raw-R2 panel with `ldsc:r2_bias="raw"` would auto-fill its correction
sample size from `ldsc:n_samples`. Legacy parquet files without these keys keep
the historical default of treating omitted bias mode as `"unbiased"`; users can
still pass `--r2-bias-mode raw --r2-sample-size N` for external raw R2 files
when running in `rsid` or `chr_pos` mode.

### 2.6 Canonical Example

**Schema metadata:**

```
ldsc:sorted_by_build = "hg19"
ldsc:row_group_size  = "50000"
ldsc:n_samples       = "3202"
ldsc:r2_bias         = "unbiased"
ldsc:schema_version  = "1"
ldsc:artifact_type   = "ref_panel_r2"
ldsc:snp_identifier  = "chr_pos_allele_aware"
ldsc:genome_build    = "hg19"
```

**Row group footer statistics (first three groups of a chr1 file):**

| Row group | Rows | POS_1 min | POS_1 max |
|---|---|---|---|
| 0 | 50 000 | 752 566 | 1 823 421 |
| 1 | 50 000 | 1 823 422 | 2 910 884 |
| 2 | 50 000 | 2 910 885 | 4 011 203 |
| … | … | … | … |

**Data rows (first rows of the file):**

| CHR | POS_1 | POS_2 | SNP_1 | SNP_2 | A1_1 | A2_1 | A1_2 | A2_2 | R2 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 752 566 | 776 546 | rs3094315 | rs2905036 | A | C | A | G | 0.8214 |
| 1 | 752 566 | 800 007 | rs3094315 | rs11240777 | A | C | C | T | 0.1342 |
| 1 | 752 566 | 817 186 | rs3094315 | rs4040617 | A | C | A | G | 0.0891 |

Note: rows are sorted by non-decreasing `POS_1`. The left SNP (`POS_1 = 752 566`) appears
in multiple consecutive rows — one per right-side neighbor within the LD window.

---

### 2.7 SNP Metadata Sidecar

The R² parquet stores **only pairwise data** (`CHR`, `POS_1`, `POS_2`,
`SNP_1`, `SNP_2`, `R2`). Per-SNP metadata — chromosome, base position, rsID, genetic
distance, and minor allele frequency — is kept in a separate **metadata sidecar**
file. The sidecar is optional, but strongly recommended for production use.

Treat the parquet R² file and the metadata sidecar as a paired artifact when
possible. Together they define one logical reference panel: the parquet provides
pairwise LD entries, while the sidecar provides the complete per-SNP metadata
needed for MAF filtering, common counts, and cM windows.

#### Fallback when the sidecar is absent

R2 parquet files are required. In `r2_dir` mode, chromosome `N` must have
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
| `A1`, `A2` | standard allele aliases | Allele-aware modes | Reference-panel orientation alleles |

In `rsid`-family modes the `SNP` column is mandatory. In `chr_pos`-family modes
`CHR` and `POS` are mandatory. A file providing all identity columns works in
any mode and is recommended for production use. In base modes, `A1/A2` are
passive for identity; in allele-aware modes, missing or unusable `A1/A2` makes
the sidecar malformed.

#### Role in the downstream workflow

When present, the metadata sidecar feeds directly into three steps of the partitioned-LDSC workflow
(see `docs/current/partitioned-ldsc-workflow.md`, §6):

1. **Reference panel universe (A').** `load_metadata()` reads the sidecar,
   optionally applies `_apply_snp_restriction()` when `ref_panel_snps_file` or
   `use_hm3_ref_panel_snps` is set, and returns the restricted per-SNP table A'.
   This is what the annotation bundle is aligned against (`B_chrom ∩ A'_chrom`)
   to materialize `ld_reference_snps`.

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

In the workflow notation of `docs/current/partitioned-ldsc-workflow.md`, the raw
reference-panel SNP universe `A` comes from this paired artifact. When the
sidecar is present, its SNP table is authoritative and parquet pair rows are
used only to answer LD window queries once `ld_reference_snps = B ∩ A'` has been
established. When the sidecar is missing, the loader falls back to scanning
parquet endpoint columns to synthesize a partial SNP table with the limitations
listed above.

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

### 3.4 Auto Decoded Row-Group Cache

Before the parquet sliding-block loop starts, `ld_score_var_blocks_from_r2_reader`
calls `SortedR2BlockReader.configure_auto_row_group_cache(block_left, snp_batch_size)`.
The reader simulates the same index-window query sequence that the LD-score
sliding algorithm will issue for the chromosome. For every adjacent query pair it
computes the union of required row groups from the footer bounds, then sets:

```python
capacity = min(max_adjacent_union + 1, num_row_groups)
```

The cache unit is one decoded parquet row group. The cached payload is numeric
arrays only: endpoint matrix indices `i:int32`, `j:int32`, and `r2:float32`.
In base `rsid` mode, `SNP_1`/`SNP_2` strings are converted to retained-SNP
matrix indices during row-group decode. In base `chr_pos` mode, `POS_1`/`POS_2`
are mapped to retained-SNP matrix indices during decode. In allele-aware modes,
endpoint alleles are included in the effective key. Unmapped endpoints are
dropped at decode time.

The cache is per chromosome and is discarded when the next chromosome reader is
created. There is no public row-group cache-size option. Cache diagnostics
(capacity, hits, misses, evictions, and parquet row-group reads) are logged at
`DEBUG` level only.

Correctness does not depend on cache state. Every query recomputes its full
required row-group set from footer bounds. A cache miss reads and decodes the
row group from parquet; it never means that a pair is absent.

### 3.5 Window Queries (`_query_union_rows`)

For each sliding-window block, `cross_block_matrix` and `within_block_matrix` call
the canonical index-window query path, which:

1. Identifies overlapping row groups from `_rg_bounds`:
   ```python
   rg_idxs = [i for mn, mx, i in self._rg_bounds if mn <= pos_max and mx >= pos_min]
   ```
2. Retrieves each required row group from the decoded cache, or reads one row
   group from parquet on cache miss:
   ```python
   tbl = self._pf.read_row_group(rg_idx, columns=[...])
   ```
3. Converts Arrow columns to numpy arrays and maps pair endpoints to retained-SNP
   matrix indices during row-group decode.
4. Applies optional R² bias correction (`_transform_r2`) during decode.
5. Concatenates decoded row groups and filters numeric `i`/`j` endpoints to the
   current union/index window.
6. Returns numeric pair rows for the dense matrix builder. The dense matrix path
   still filters to the block-local matrix region before assignment.

**Coarse-file behaviour:** if the file has few, large row groups, step 1 may return
all row groups for every query — performance degrades to O(N) per query but
correctness is maintained.

### 3.6 Dense Matrix Construction

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
runtime parameter; the reader selects which columns to use for index mapping:
`SNP_1`/`SNP_2` for base `rsid`, `POS_1`/`POS_2` for base `chr_pos`, and the
corresponding endpoint allele-aware effective key for allele-aware modes,
without requiring different package-built parquet files.

**Legacy raw-schema parquets.** Parquets in the old schema (containing
`hg19_pos_1`, `hg38_pos_1`, etc. but no canonical `POS_1`/`POS_2` columns) are
detected at `__init__` by the absence of `POS_1` and `POS_2` from the schema. The
reader emits a `WARNING` and falls back to the existing `pyarrow.Dataset` full-scan
path:

> *"{path} uses the legacy raw schema. Row-group pruning is disabled and query
> performance will be severely degraded."*

The full-scan fallback is numerically identical to the legacy implementation;
correctness is preserved. There is no public `normalize-r2-parquet` subcommand
in the current CLI; regenerate canonical files with `ldsc build-ref-panel`.

---

## 5. Affected Modules

| Module | Change |
|---|---|
| `_kernel/ref_panel_builder.py` | `write_r2_parquet`: require PyArrow; assert sort invariant on each incoming pair; write new 6-column schema; preserve canonical dtypes for empty outputs; write `ldsc:sorted_by_build` and `ldsc:row_group_size` metadata; default `row_group_size=50_000` |
| `_kernel/ldscore.py` — `SortedR2BlockReader.__init__` | Detect schema (canonical vs legacy raw); canonical path: open as `pq.ParquetFile`, build `_rg_bounds` index from footer, validate build (3-tier), warn if coarse row groups; legacy path: open as `pyarrow.Dataset`, emit deprecation warning, proceed with full-scan fallback |
| `_kernel/ldscore.py` — `SortedR2BlockReader` decoded cache | Canonical path: auto-size a chromosome-local LRU cache from `block_left` and `snp_batch_size`; cache decoded row groups as numeric `i`, `j`, `r2` arrays; log DEBUG-only cache diagnostics |
| `_kernel/ldscore.py` — `SortedR2BlockReader._query_union_rows` | Canonical path: row-group index lookup + decoded cache lookup or `read_row_group` miss + numeric endpoint filtering; legacy path: existing `dataset.to_table(filter=...)` behaviour unchanged |
| `_kernel/ldscore.py` — `read_sorted_r2_presence` | Remove standalone function; the parquet backend no longer performs a runtime SNP-presence scan |
| `_kernel/ldscore.py` — `compute_chrom_from_parquet` | Use the sidecar-defined `A'` directly when forming `ld_reference_snps`; no parquet presence scan is performed |
| Tests — canonical schema | Add fixtures writing new 6-column sorted parquets; assert row-group index is built; assert window queries return correct pairs |
| Tests — legacy schema | Retain existing raw-schema fixtures; assert deprecation warning is emitted; assert numerical output is unchanged |
