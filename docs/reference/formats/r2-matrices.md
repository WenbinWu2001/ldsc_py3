# R² Matrix: Format Specification and Read Pipeline

Last updated on: 2026-10-28

This document defines the canonical **index** format for R² matrices and describes how the read pipeline validates, binds, and streams them into the LD-score accumulation.

The R² Parquet file is a compact pair table of **sidecar-row indices**; it carries no SNP identity of its own.
Identity lives once-per-SNP in the paired sidecar, and the Parquet file references it by integer index.
The two files are an inseparable pair, bound by a content hash.

---

## 1. Motivation

A pairwise R² matrix stores, for every SNP pair within an LD window, the squared
correlation between allele dosages. The builder computes each pair as a pair of
integer indices into the retained reference-SNP table; the read path ultimately
needs those same indices to accumulate LD scores. The index format stores the
indices directly instead of expanding them to positions, rsIDs, and alleles
(which the old 10-column format did, only for the reader to collapse them back).

Removing that expansion/collapse round-trip:

- shrinks the identity overhead to two `int32` columns,
- makes the per-pair read hot path a pure integer gather (mode-independent), and
- keeps the stored pairs independent of how sidecar SNPs are matched, so each reader applies the recorded mode in its remap (§2.4).

A `SIGN_R` bit (sign of the Pearson correlation r) is preserved for downstream
consumers, though LD-score calculation uses the unbiased R² values and ignores this sign bit.

---

## 2. Parquet Format Specification

### 2.1 Column Schema (exactly four columns)

| Column | Type | Description |
|---|---|---|
| `IDX_1` | `int32` | left endpoint: 0-based row index into the R² matrix sidecar |
| `IDX_2` | `int32` | right endpoint: 0-based row index into the R² matrix sidecar |
| `R2` | `int16` (quantized) | unbiased squared correlation, symmetric int16 quantization (scale 32767); reader dequantizes to `float32`. See §2.3. |
| `SIGN_R` | `bool` (bit-packed) | `True` ⇔ Pearson r ≥ 0 in the sidecar's A1/A2 orientation |

There are no `CHR`, `POS_*`, `SNP_*`, or allele columns. The chromosome and all
per-SNP identity come from the paired sidecar (§2.7). `SIGN_R` is stored as Arrow
`bool` (1 bit/value, RLE-friendly); it is the sign of the correlation r, not of
the (possibly negative) unbiased R² value, and is meaningful only together with
the sidecar's A1/A2 alleles.

### 2.2 Invariants

- `IDX_1 < IDX_2` for every row (canonical pair orientation).
- **Each unordered SNP pair appears exactly once.** With `IDX_1 < IDX_2`, every
  `(IDX_1, IDX_2)` row is unique — the builder never emits a pair twice. This is
  an index-level (identifier-mode-independent) property: `IDX_k` are sidecar row
  numbers, so pair identity does not depend on `snp_identifier`. This is
  **guaranteed by the package builder**; the §4 SHA-256 binding separately
  fixes the sidecar identity/order but does not authenticate or deduplicate pair rows.
  The read path trusts the builder contract and does not scan the parquet for duplicate rows.
  There is no per-window or whole-file deduplication.
  The only uniqueness enforcement at read time is the
  up-front, mode-aware collapse rejection at remap time (§3.1), which guards the
  one case that *can* induce a duplicate in the retained matrix: two distinct
  R² matrix SNPs sharing an identity key under the active mode.
- Rows are sorted by non-decreasing `IDX_1`. Because sidecar row order is
  position-sorted, this equals `POS_1` order, so row-group pruning works on
  `IDX_1` footer statistics.
- Each file contains exactly one chromosome (defined by the paired sidecar).
- `int32` indices suffice; the writer asserts `n_snps < 2³¹`.
- All pairs lie within the configured LD window at write time.

The writer **asserts** the sort invariant per batch on `IDX_1`. A pair arriving
with `IDX_1` below the previous row's raises immediately, because a violated sort
invariant corrupts row-group footer statistics and silently produces wrong window
queries.

### 2.3 Row Group Size, Compression, and Column Encodings

The default is **50,000 rows per row group**.
Row-group granularity affects pair-query pruning and read overhead; the groups touched depend on pair density and the requested endpoints.
Column chunks use **zstd level 9** when available, with a snappy fallback and a `UserWarning`.
Compression is recorded per column-chunk in the footer and auto-detected on read.

**Per-column Parquet encodings** are set explicitly to match each column's
structure:

| Column | Encoding | Rationale |
|---|---|---|
| `IDX_1` | `RLE_DICTIONARY` (PyArrow default) | Repeated left endpoints form constant runs suitable for dictionary and run-length encoding. |
| `IDX_2` | **`DELTA_BINARY_PACKED`** | Sorted right endpoints within each left-endpoint run have small differences suited to delta encoding. `use_dictionary` is restricted to `IDX_1` so it does not override this choice. |
| `R2` | **`int16` quantization + `BYTE_STREAM_SPLIT`** | Symmetric int16 (`round(R²·32767)`, clipped `[-32767, 32767]`) replaces float32; byte-stream splitting compresses the two byte planes separately. Quantization is lossy; error and storage ratios depend on inputs and downstream weighted sums. |
| `SIGN_R` | `RLE` + bit-packing (PyArrow bool default) | Arrow `bool` is natively bit-packed at 1 bit/value; RLE handles the encoding automatically. |

**R² quantization (lossy) — endpoint-exact.**
The scale `32767` is the int16 maximum, so the value `1.0` round-trips exactly.
The writer upper-clips the unbiased estimate at `1.0`, encodes it as `32767`, and the reader decodes `32767 / 32767 = 1.0`.
Negative unbiased values are preserved within the encoding's range, not floored at zero.
For values within the clipping range, ideal rounding contributes at most half a step, $1/(2\cdot32767)$, plus float32 decoding roundoff.
Clipping outside that range has additional error.
An LD-score sum combines many such errors with annotation weights, so historical per-SNP differences are workload measurements, not universal bounds or an assertion of losslessness.
This differs from [binary annotation packing](../../dev/design/memory-and-cpu.md#annotation-format-policy), which preserves exact 0/1 values.

**Adoption evidence (2026-06-01 workload measurement).**
Before int16 storage was adopted, a genome-wide evaluation of one reference panel (all 22 autosomes, 473.5 million pairs, $N=3202$ samples) compared symmetric int16 values with the source float32 values.
No value saturated the clipping range.
The per-pair absolute error had maximum 1.5e-5 and mean 8e-6, about 20 times below the unbiased-estimator noise floor $1/N \approx 3.1\times10^{-4}$.
The resulting LD scores changed by 2.1e-4 on average and by at most 2.1e-3 (chromosome 3).
Storage fell from 2,710.9 MiB (float32) to 1,014.2 MiB (int16 with `BYTE_STREAM_SPLIT`), a 62.6% reduction.
Lossless `BYTE_STREAM_SPLIT` on float32 alone saved 31.5%, and int8 was rejected because its maximum error of about 2e-3 is visible in the small-R² range that holds about 59% of pairs.
These numbers describe that R² matrix; they are not error bounds for other inputs.

All encodings are recorded in the parquet footer and PyArrow auto-decodes them.
The int16→float32 **dequantization is a reader-side step** (divide by
`ldsc:r2_scale`), not a Parquet codec; see §3.2.
Readers reject any other R² dtype, encoding, or scale; an R² matrix with a `float32` R² column must be rebuilt.

### 2.4 SNP Identifier Modes — handled entirely off the parquet

The four identifier modes are `rsid`, `rsid_allele_aware`, `chr_pos`, and
`chr_pos_allele_aware` (package default `chr_pos_allele_aware`). The index
parquet stores only indices; the mode affects only how sidecar SNPs are matched
to the retained reference SNPs when the reader builds its remap (§3).
`ldsc:snp_identifier` records the mode the R² matrix was built in, and readers enforce it: `ldscore` requires the active `--snp-identifier` to equal the recorded mode, and `query-r2-matrix` accepts the recorded mode or its base form (`rsid` or `chr_pos`).

### 2.5 Parquet Schema Metadata

| Key | Type | Description |
|---|---|---|
| `ldsc:artifact_type` | string | `r2_matrix` |
| `ldsc:snp_identifier` | string | mode the R² matrix was built in; `ldscore` requires an exact match, `query-r2-matrix` also accepts its base form |
| `ldsc:genome_build` | string | genome build of the sidecar coordinates |
| `ldsc:sorted_by_build` | string | build defining the `IDX`/position sort order (`"hg19"` or `"hg38"`); required |
| `ldsc:row_group_size` | string(int) | intended row-group size (informational) |
| `ldsc:n_samples` | string(int) | reference individuals used to compute R²; required, at least 3, and equal across an R² matrix's chromosome files for `query-r2-matrix` |
| `ldsc:r2_bias` | string | required; must be `"unbiased"` (stored values are already bias-corrected) |
| `ldsc:r2_min` | string(float) | required; unbiased-R² floor applied at build, at most 1; non-positive means complete. `ldscore` rejects a positive value; `query-r2-matrix` accepts it |
| `ldsc:n_snps` | string(int) | length of the R² matrix sidecar = size of the index space |
| `ldsc:sidecar_identity_sha256` | string(hex) | binding hash of the sidecar identity (§4) |
| `ldsc:r2_encoding` | string | required; must be `"int16_symmetric"` |
| `ldsc:r2_scale` | string(int) | required; must be `"32767"`. The reader divides the int16 values by this to recover float32 R². |
| `ldsc:ld_window_mode`, `ldsc:ld_window_value` | string, string(float) | required; build window unit (`snps`, `kb`, or `cm`) and positive width. `ldscore` rejects a window in another unit or wider than this one ([window rules](#37-ld-windows-and-sidecar-values-in-ldscore)) |

Stored as UTF-8 byte strings, accessed via
`pq.ParquetFile(path).schema_arrow.metadata`.
Both readers check every key above (except the informational `ldsc:row_group_size`) in `read_r2_parquet_contracts` before decoding any pair and report all defects of all files in one error; no value is defaulted.

### 2.6 Canonical Example

**Schema metadata:**
```
ldsc:artifact_type           = "r2_matrix"
ldsc:snp_identifier          = "chr_pos"
ldsc:genome_build            = "hg19"
ldsc:sorted_by_build         = "hg19"
ldsc:n_samples               = "3202"
ldsc:r2_bias                 = "unbiased"
ldsc:r2_min                  = "0.0"
ldsc:n_snps                  = "969730"
ldsc:sidecar_identity_sha256 = "a1b2…(64 hex chars)"
ldsc:r2_encoding             = "int16_symmetric"
ldsc:r2_scale                = "32767"
ldsc:ld_window_mode          = "cm"
ldsc:ld_window_value         = "1.0"
```

**Data rows (sorted by non-decreasing `IDX_1`):** R² is stored on disk as int16
(`round(R²·32767)`); the decoded float32 value the reader yields is shown in
parentheses.

| IDX_1 | IDX_2 | R2 (int16 → decoded) | SIGN_R |
|---|---|---|---|
| 0 | 1 | 26915 → 0.8214 | True |
| 0 | 2 | 4397 → 0.1342 | False |
| 1 | 2 | 2920 → 0.0891 | True |

`IDX_1 = 0` is sidecar row 0; the left SNP appears in multiple consecutive rows,
one per right-side neighbor within the LD window.

### 2.7 SNP Metadata Sidecar (mandatory, authoritative)

The sidecar (`chrN_meta.tsv.gz`) is the per-SNP table for the whole R² matrix and
**defines the index space**: parquet `IDX_k` is a 0-based row number into it. It
is a gzip TSV with `# ldsc:*` identity-provenance comments followed by columns
`CHR, POS, SNP, A1, A2, CM, MAF` (one row per SNP, in build/position order).
`CM` and `MAF` are written with 9 significant digits (`SIDECAR_FLOAT_FORMAT` in `_refpanel/r2_format.py`).
`build-r2-matrix` computes its cM window from `CM` rounded to that precision, so the stored pairs match the window `ldscore` rebuilds from the sidecar.

The sidecar is **required**. Because the parquet stores only indices, it is
meaningless without the exact sidecar; there is no synthesize-from-parquet
fallback. A missing sidecar is a hard error.

The sidecar feeds the partitioned-heritability workflow as before: the R² matrix SNPs `A` before filtering,
all-SNP and common-SNP counts via `MAF`, and the genetic-map (`CM`) window.

---

## 3. Read Pipeline

Implemented in `SortedR2BlockReader` (`src/ldsc/_ldscore/compute.py`), the single
point of contact between the LD-score algorithm and the parquet.

### 3.1 Open, Bind, Remap (`__init__` → `_init_index_path`)

The reader receives the **retained** metadata (the retained reference SNPs, in
matrix order) exactly as the rest of the pipeline produces it — the index format
does not change any filtering/restriction/intersection step. It then:

1. Confirms the schema is the index layout (`IDX_1/IDX_2/R2/SIGN_R` present); any other
   schema raises an actionable "regenerate with `ldsc build-r2-matrix`" error.
2. Reads `ldsc:sorted_by_build` (required) and checks it against the analysis
   build.
3. Loads the **full** R² matrix sidecar (`chrN_meta.tsv.gz`) — all `n_snps` rows in
   build order.
4. Validates the binding (§4): recomputes the identity hash and checks
   `len(sidecar) == n_snps`; mismatch is a hard error.
5. Builds `remap: build_idx → retained_idx` (int32, `-1` for SNPs that are not
   retained reference SNPs) and `retained_build_idx: matrix_idx → build_idx`
   (ascending, for pruning) via `build_index_remap`. Matching uses the same
   mode-dependent identity key (`effective_merge_key_series`) that the retained
   reference SNPs are keyed by, so the result is identical to the former per-pair
   identity lookup. **This is the only place identifier-mode matching happens —
   once per chromosome.** The remap **must be injective over retained indices**:
   if two distinct R² matrix sidecar rows share the same identity key under the
   active `snp_identifier` and both map onto the retained reference SNPs, they collapse
   onto one retained SNP and would silently double-count that SNP's pairs. Any
   such collapse is a **hard error** (rebuild the R² matrix or use an allele-aware
   `--snp-identifier`). Because the check is the injectivity of a remap built
   from mode-derived keys (`base_key [+ ":" + normalized allele set]`), it is
   automatically mode-aware. Two R² matrix rows collapse only when they share the
   full key for the active mode: `rsid` → same **rsID**; `rsid_allele_aware` →
   same **rsID + normalized allele pair**; `chr_pos` → same **CHR:POS**;
   `chr_pos_allele_aware` → same **CHR:POS + normalized allele pair**. So two
   rows at the same `CHR:POS` with different rsIDs stay distinct in both
   rsID-family modes, and a multi-allelic site splits into distinct SNPs in both
   allele-aware modes (the allele set is order-independent, so `A/G` ≡ `G/A`).
   This up-front rejection is what
   lets §3.4 accumulate LD scores with no pair deduplication; with each
   `(IDX_1, IDX_2)` pair unique (§2.2) and the remap injective, every decoded
   retained pair is necessarily unique.
6. No window index is built: LD-score computation streams every row group once
   (`iter_all_pairs`, §3.3). The coarse-row-group startup warning from the parquet
   footer statistics is retained.

<a id="32-decode-_decode_index_row_group"></a>

### 3.2 Decode (`_decode_index_batch`)

Per bounded Arrow record batch, the reader reads only `IDX_1, IDX_2, R2` and gathers:

```python
i = remap[IDX_1]; j = remap[IDX_2]
keep = (i >= 0) & (j >= 0)
```

Endpoints that are not retained reference SNPs map to `-1` and are dropped.
The reader **dequantizes** the int16 column (`r2 = int16_values / ldsc:r2_scale`, with the scale taken from the validated contract) to `float32` and applies no further correction.
`SIGN_R` is not read (unused by LD-score computation).
Each decoded batch is a numeric
`(i:int32, j:int32, r2:float32)` triple of retained-matrix index pairs.

The public chromosome workflow calls `R2MatrixRefPanel.prepare_chromosome`, which aligns annotation rows to the filtered sidecar, takes authoritative CM/MAF from that sidecar, validates LD windows, and calls `build_reader`. That single reader factory checks the parquet against the shared contract (`read_r2_parquet_contract`) and the active identity configuration; production and adapter tests share the same path. The numerical kernel receives this prepared state, and its context manager closes the reader after success or failure. Stored values are already unbiased, so ordinary and regression-weight LD scores apply no sample-size correction.

### 3.3 Pair streaming (`iter_all_pairs`)

`iter_all_pairs()` visits every row group once in `IDX_1` order and calls `iter_batches` with a maximum of 65,536 rows.
Large stored row groups are decoded in bounded record batches rather than materialized whole.
There is no row-group pruning, decoded-row-group cache, or dense LD-window matrix.
The driver consumes these arrays into a separate bounded pair buffer (§3.4).
This avoids retaining every pair in a chromosome, but it is not a total RSS bound (§3.6).

### 3.4 LD-score accumulation (`ld_score_streaming_from_r2_reader`)

Conceptually, the full LD-score matrix is $C=RA$, where $R$ is symmetric with unit diagonal and $A$ contains annotations plus the optional regression-weight column.
Production retains only requested output rows $O$, computing $C_{\mathrm{out}}=(RA)[O,:]$ while preserving all retained reference SNPs.
The corresponding pair updates are:

```text
cor_sum[k] ← annot[k]                         # diagonal once, k in output rows
for each stored pair (i, j, r2) with i < j and i ≥ block_left[j]:
    if i is an output row: cor_sum[i] += r2 · annot[j]
    if j is an output row: cor_sum[j] += r2 · annot[i]
```

`block_left[j]` applies the read-time LD window, which may be narrower than the build window.
The driver splits decoded batches to buffer at most `_CSR_CHUNK_PAIRS` pairs, currently 16 million, before projection.
[`ProjectionAccumulator.add_pairs()`](../../../src/ldsc/_ldscore/compute.py) filters the window and output targets, expands the required pair orientations, and constructs a rectangular CSR over the chunk's involved output rows and reference SNPs.
It multiplies in float64 using selected annotation tiles, reusing each reference-SNP tile across output blocks.
Shared baseline columns and weights are projected together, followed by focal column tiles for the active execution batch.
The next execution batch repeats traversal.
The production path does not require a full `m×m` upper-triangular CSR or a complete dense annotation matrix.

Scores are accumulated in float64 and published in the existing float32 format.
Changes in tile or pair grouping can reorder floating-point sums; numerical tolerances must be verified against independent expectations, not inferred from algebra alone.
The [mathematical specification](../methods/ld-score-computation.md) distinguishes the conceptual full product from the selected-output implementation.

### 3.5 Mode-agnosticism

Because the parquet carries no identity, the four identifier modes read the same
bytes; only `remap` construction (step 5) differs. The cross-mode parity test
builds one index parquet and asserts identical LD scores in all four modes.

### 3.6 Read-side memory (streaming)

The main float64 score array has payload $8\,|O|\,n_a$ bytes, where $n_a$ includes only the active execution batch, baseline columns, and optional weight.
Other live allocations include the decoded record batch, buffered/concatenated pairs, remap arrays and sidecar, bidirectional pair expansion, CSR construction/sorting, annotation tiles, sparse slices, gathered outputs, products, and public-format conversion.
The default 16 MiB projection read budget targets one numerical operand tile; it does not cover those other allocations.
Likewise, 16 million buffered pairs is a count limit, not a 0.45 GiB RSS guarantee.

There is no dense window matrix or all-chromosome pair cache.
A wider or denser R² matrix therefore increases streamed work without requiring all its pairs to be resident simultaneously.
The read-time window still filters pairs before sparse construction, so it can change the working set and arithmetic even though all stored pairs are read.
Each chromosome worker adds its own working set, while staged output and file cache have separate lifetimes.
See the [package resource model](../../dev/design/memory-and-cpu.md#precision-snp-sets-and-resource-accounting).

The builder has a different working set.
Unrestricted builds stream genotypes from disk; with `--ref-panel-snps-file`, the selected genotypes are read into RAM, with individual filtering applied during the read.
Pair computation keeps the current batch, window-spanning carry-over columns and pending pair batches, so memory depends on retained individuals, selected SNPs, window density and `--snp-batch-size`; wider windows also increase computation and output size.
A positive `--r2-min` reduces the stored pairs, and `ldscore` rejects such an R² matrix (`query-r2-matrix` still reads it).
See `R2MatrixBuilder._build_chromosome` in `src/ldsc/_refpanel/builder.py` and `yield_pairwise_r2_rows` in `src/ldsc/_refpanel/build_kernel.py`.

### 3.7 LD windows and sidecar values in `ldscore`

`ldscore` builds its LD window from the sidecar: `--ld-wind-snps` counts sidecar rows, `--ld-wind-kb` uses sidecar `POS`, and `--ld-wind-cm` uses sidecar `CM`.
Exactly one of the three options is required.

**Window coverage.**
An R² matrix stores only the SNP pairs inside the window it was built with, recorded as `ldsc:ld_window_mode` and `ldsc:ld_window_value`, so `ldscore` accepts only a window the R² matrix is known to cover:

| Requested window | Result |
| --- | --- |
| Same unit, no wider than the build window | Accepted. |
| Same unit, wider than the build window | `LDSCUsageError`: "is wider than the input R2 matrix window". |
| Different unit, for example a kb R² matrix with `--ld-wind-snps` | `LDSCUsageError`: "windows in different units cannot be checked for coverage". |

A window in a different unit is rejected even when the R² matrix happens to cover it, because coverage is not inferred from per-SNP coordinates.
The CLI checks every selected chromosome file in its header preflight, before pair scans, annotation preparation or LD-score computation.
`LDScoreCalculator.run` runs the same checks for every chromosome before computing any, and lists all issues in one `LDSCInputError`.
`R2MatrixRefPanel.prepare_chromosome` repeats the check for each chromosome it prepares, for example through `LDScoreCalculator.compute_chromosome`.
Because `build-r2-matrix` chooses the stored pairs from the `CM` values as written to the sidecar (section 2.7), a window no wider than the R² matrix's never needs a pair the R² matrix did not store, in any unit.
A PLINK `--ld-wind-cm` run matches an R² matrix cM run only when both use the same `CM` values: `.bim` column 3, or the genetic map the R² matrix was built with.
Code: `validate_ldscore_window_within_r2_matrix_window` in `src/ldsc/_ldscore/compute.py` and `inspect_r2_headers` in `src/ldsc/_runtime/preflight.py`.

**`--ld-wind-snps` counts R² matrix rows.**
For an R² matrix, `--ld-wind-snps K` spans `K` sidecar rows on each side, whatever SNPs the run keeps.
When the run keeps fewer SNPs than the R² matrix (`--maf-min`, `--ref-panel-snps-file`, or annotations that cover part of the R² matrix), its SNPs have fewer than `K` reference SNPs on each side.
The PLINK backend counts the SNPs the run keeps, as legacy LDSC did, so the two backends then differ for the same options.
Each such run logs one warning with the number of affected chromosomes, one example with its kept and R² matrix SNP counts, and the location of the per-chromosome counts.
Every R² matrix `--ld-wind-snps` run writes `diagnostics/ld_window_snp_counts.tsv` with `CHR`, `panel_snps`, `kept_snps`, `ld_wind_snps` and `counting_basis` (`panel_rows`), listed in the `files` entry of `metadata.json`.
SNPs that share a position are computed in sidecar row order and written sorted by `CHR`, `POS`, then `SNP`.
That output order matches an annotation-ordered run only because SNPs that share `CHR`, `POS` and `SNP` never reach computation: the identity policy drops them as duplicate or multi-allelic.
If that policy changes to keep such SNPs, `R2MatrixRefPanel.prepare_chromosome` would also need to return them to annotation order.
Code: `build_window_coordinates` in `src/ldsc/_ldscore/compute.py`, `R2MatrixRefPanel.prepare_chromosome` in `src/ldsc/_refpanel/adapters.py`, and `_ld_window_snp_counts` in `src/ldsc/_ldscore/calculator.py`.

**`CM` and `MAF` come from the reference panel.**
Annotation `CM` and `MAF` values are ignored.
For an R² matrix the sidecar supplies both for every retained SNP; the PLINK backend uses `.bim` `CM`, or a genetic map interpolated at the `.bim` positions, and genotype `MAF`.
`--ld-wind-cm` requires a non-missing `CM` for every retained SNP; otherwise chromosome computation raises "--ld-wind-cm requires non-missing CM values for all retained SNPs."
Sidecars from builds without a genetic map store `CM=NA`; those R² matrices were built with a kb or SNP-count window, so the window-unit rule rejects `--ld-wind-cm` first.
For PLINK input, `--ld-wind-cm` also rejects an unusable `.bim` `CM` (fewer than two distinct finite values on a chromosome) with an error that `--yes-really` does not bypass; supply `--genetic-map-hg19-sources` or `--genetic-map-hg38-sources` instead ([troubleshooting](../../user-guide/troubleshooting.md#ldscore-unusable-cm-for-ld-wind-cm)).
Any missing or non-numeric sidecar `MAF` is an error when the sidecar is read.
`--maf-min` filters on the sidecar `MAF`, and the common-SNP counts and the common-SNP overlap matrix use inclusive `MAF >= common_maf_min` from the sidecar; all-SNP counts come from the retained annotation rows.
Without baseline annotations, `ldscore` builds an all-ones `base` annotation from the sidecar.

---

## 4. Index ↔ Sidecar Binding

The binding breaks if the R² matrix's SNP **set or order** changes, so the hash
covers exactly that:

- **Basis:** SHA-256 over the newline-joined `CHR:POS:A1:A2` of each sidecar row,
  in order, computed in memory at build and load (never written as a column). The
  rsID is excluded (it is the most volatile identifier across dbSNP builds);
  `CHR:POS:A1:A2` is the stable physical identity and disambiguates multi-allelic
  sites. `CM`/`MAF` are excluded so a benign recompute does not trip the binding.
- **Storage:** hex digest in `ldsc:sidecar_identity_sha256`, length in
  `ldsc:n_snps`.
- **Validation (load):** recompute from the full sidecar; hard-error on digest
  mismatch (wrong/reordered/edited sidecar) or `n_snps` length mismatch.

---

## 5. Caveats and Constraints

**Current names (2026-10-28).** Use `ldsc build-r2-matrix` and `ldsc query-r2-matrix`, which replaced `build-r2-panel` and `query-r2`, and the required `SIGN_R` Parquet column.
R² matrices whose Parquet footer records `ldsc:artifact_type="ref_panel_r2"`, or that use the older `SIGN` schema, are unsupported and must be rebuilt.
`query-r2-matrix` returns nullable `sign_r` (+1/-1 in query allele orientation), alongside `r2`, `r`, and `status`.
There are no compatibility aliases.
Workflow logs use `build-r2-matrix[.chrN].log`; the chromosome artifact filenames `chrN_r2.parquet` and `chrN_meta.tsv.gz` are unchanged.

**Sidecar mandatory.** The parquet is opaque without its exact sidecar; ship them
together. There is no metadata-synthesis fallback.

**Build specificity.** The parquet encodes one build's sort order. The 3-tier
inference of the old format is gone — the index format always records
`ldsc:sorted_by_build`, and its absence is a hard error. Maintain separate hg19
and hg38 R² matrices.

**One chromosome per file.** The reader is instantiated per chromosome and
expects a single-chromosome parquet + sidecar.

**Clean break — index only.** The former 10-column canonical schema, the raw
external-input read path, and the synthesize-from-parquet metadata fallback are
removed. External R² formats (e.g. TOP-LD pairwise tables, SBayesR binary LD
matrices) are not read directly; a future converter would turn them into the
canonical index format offline. Existing 10-column R² matrices must be regenerated.

---

## 6. Affected Modules

| Module | Role |
|---|---|
| `_refpanel/r2_format.py` | `write_r2_parquet`: 4-column index schema, sort assertion on `IDX_1`, binding metadata; `sidecar_identity_sha256` (re-exported from `snp_identity`). Pairs flow as columnar `PairColumns` `(i, j, r2, sign)` batches from `yield_pairwise_r2_rows` straight into the writer — no per-pair dict — and int16 quantization/encodings are applied at table build; the writer emits the current `SIGN_R` schema. |
| `_core/snp_identity.py` | `sidecar_identity_sha256` binding hash |
| `_refpanel/plink.py` | `PlinkBEDFile` genotype source for the build: selective reads retain the selected packed payload; unrestricted R² matrix builds opt into disk streaming. Both decode bounded standardized blocks |
| `_refpanel/builder.py` | build loop: sidecar built first, hash + `n_snps` passed to the writer |
| `_refpanel/r2_format.py` | `SortedR2BlockReader` index path: full-sidecar load, binding validation, `build_index_remap`, gather decode, `iter_all_pairs` streaming; `ld_score_streaming_from_r2_reader` pair accumulation; raw/former-format and block/query/cache paths removed |
| `_ldscore/compute.py` | `ProjectionAccumulator`: selected-output float64 scores, rectangular sparse pair operators, bounded annotation operands, and shared weight projection |
| `_refpanel/adapters.py` | mandatory sidecar loading; shared `prepare_chromosome` alignment/window policy and `build_reader` contract and identity checks |
| `tests/` | index writer/reader/binding/remap tests; cross-mode parity gate; build→read parity |
