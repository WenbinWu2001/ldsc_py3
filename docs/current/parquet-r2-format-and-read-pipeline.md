# Parquet R2 Table: Format Specification and Read Pipeline

This document defines the canonical **index** format for parquet-backed pairwise
R² reference panels and describes how the read pipeline validates, binds, and
streams them into the LD-score accumulation.

The R² parquet is a compact pair table of **sidecar-row indices**; it carries no
SNP identity of its own. Identity lives once-per-SNP in the paired metadata
sidecar, and the parquet references it by integer index. The two files are an
inseparable pair, bound by a content hash.

---

## 1. Motivation

A pairwise R² panel stores, for every SNP pair within an LD window, the squared
correlation between allele dosages. The builder computes each pair as a pair of
integer indices into the retained reference-SNP table; the read path ultimately
needs those same indices to accumulate LD scores. The index format stores the
indices directly instead of expanding them to positions, rsIDs, and alleles
(which the old 10-column format did, only for the reader to collapse them back).

Removing that expansion/collapse round-trip:

- shrinks the identity overhead to two `int32` columns,
- makes the per-pair read hot path a pure integer gather (mode-independent), and
- lets a **single** parquet serve all four SNP identifier modes.

A `SIGN` bit (sign of the Pearson correlation r) is preserved for downstream
consumers, though LD-score regression — which uses unsigned R² — ignores it.

---

## 2. Parquet Format Specification

### 2.1 Column Schema (exactly four columns)

| Column | Type | Description |
|---|---|---|
| `IDX_1` | `int32` | left endpoint: 0-based row index into the panel sidecar |
| `IDX_2` | `int32` | right endpoint: 0-based row index into the panel sidecar |
| `R2` | `int16` (quantized) | unbiased squared correlation, symmetric int16 quantization (scale 32767); reader dequantizes to `float32`. See §2.3. |
| `SIGN` | `bool` (bit-packed) | `True` ⇔ Pearson r ≥ 0 in the sidecar's A1/A2 orientation |

There are no `CHR`, `POS_*`, `SNP_*`, or allele columns. The chromosome and all
per-SNP identity come from the paired sidecar (§2.7). `SIGN` is stored as Arrow
`bool` (1 bit/value, RLE-friendly); it is the sign of the correlation r, not of
the (possibly negative) unbiased R² value, and is meaningful only together with
the sidecar's A1/A2 alleles.

### 2.2 Invariants

- `IDX_1 < IDX_2` for every row (canonical pair orientation).
- **Each unordered SNP pair appears exactly once.** With `IDX_1 < IDX_2`, every
  `(IDX_1, IDX_2)` row is unique — the builder never emits a pair twice. This is
  an index-level (identifier-mode-independent) property: `IDX_k` are sidecar row
  numbers, so pair identity does not depend on `snp_identifier`. This is
  **guaranteed by construction** (builder + the §4 SHA-256 sidecar binding, which
  fixes the exact `(IDX_1, IDX_2)` set), so the read path **trusts it and does
  not scan the parquet for duplicate rows** — there is no per-window and no
  whole-file deduplication. The only uniqueness enforcement at read time is the
  up-front, mode-aware collapse rejection at remap time (§3.1), which guards the
  one case that *can* induce a duplicate in the retained matrix: two distinct
  panel SNPs sharing an identity key under the active mode.
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

Default: **50 000 rows per row group**, tuned so most 1 Mb window queries touch
1–2 row groups. Recommended range 25 000–200 000. Column chunks are compressed
with **zstd level 9** when available (snappy fallback, with a `UserWarning`).
Compression is recorded per column-chunk in the footer and auto-detected on read.

**Per-column Parquet encodings** are set explicitly to match each column's
structure:

| Column | Encoding | Rationale |
|---|---|---|
| `IDX_1` | `RLE_DICTIONARY` (PyArrow default) | Long constant runs (each left SNP appears ~3k times at 1cM/1KG density); dictionary + RLE reduces to ≈0.05 bits/pair. |
| `IDX_2` | **`DELTA_BINARY_PACKED`** | Within each `IDX_1` run, right-neighbor indices are sorted with a median forward-gap of 1. DELTA stores successive differences, collapsing ~650 MB to ~6 MB vs the default dictionary path. `use_dictionary` is restricted to `IDX_1` only so PyArrow cannot override this encoding. |
| `R2` | **`int16` quantization + `BYTE_STREAM_SPLIT`** | R² is stored as symmetric int16 (`round(R²·32767)`, clipped `[-32767, 32767]`) instead of float32. `BYTE_STREAM_SPLIT` separates the two int16 byte planes so zstd compresses the low-entropy high byte (R² clusters near 0) apart from the noisy low byte. Genome-wide this cuts the R² column ~63% vs float32 PLAIN. Lossy but effectively lossless for LD scores: per-pair `|error| ≤ 1.5e-5` (½ step, ~20× below the `1/N` unbiased-correction floor), per-SNP LD-score delta `≤ ~2e-3`. |
| `SIGN` | `RLE` + bit-packing (PyArrow bool default) | Arrow `bool` is natively bit-packed at 1 bit/value; RLE handles the encoding automatically. |

**R² quantization (lossy) — endpoint-exact.** The scale `32767` is the int16
maximum, so the maximum R² round-trips exactly: the writer upper-clips the
unbiased estimate at `1.0`, `round(1.0·32767) = 32767`, and the reader decodes
`32767 / 32767 = 1.0` exactly. Negative unbiased values (~15% of pairs) are
preserved, not floored. The `[-32767, 32767]` clip is defensive — with the `≤1.0`
clip and `|negatives| ~ 1/N` it never fires on real data.

All encodings are recorded in the parquet footer and PyArrow auto-decodes them.
The int16→float32 **dequantization is a reader-side step** (divide by
`ldsc:r2_scale`), not a Parquet codec; see §3.2. A `float32` R² column (legacy
panels) is detected by dtype and read unscaled, so old panels keep working.

### 2.4 SNP Identity Modes — handled entirely off the parquet

The four identifier modes are `rsid`, `rsid_allele_aware`, `chr_pos`, and
`chr_pos_allele_aware` (package default `chr_pos_allele_aware`). The index
parquet is **identifier-mode-agnostic**: it stores only indices, so the *same*
file is read in any mode. The mode affects only how sidecar SNPs are matched to
the analysis universe when the reader builds its remap (§3); it never changes
which parquet is required. `ldsc:snp_identifier` in the metadata records the mode
the panel was built for, as provenance only.

### 2.5 Parquet Schema Metadata

| Key | Type | Description |
|---|---|---|
| `ldsc:schema_version` | string(int) | identity-provenance contract version, `"1"` (shared across all LDSC artifact types; **not** a per-layout marker — the index layout is identified structurally by its columns and the keys below) |
| `ldsc:artifact_type` | string | `ref_panel_r2` |
| `ldsc:snp_identifier` | string | mode the panel was built for (provenance only) |
| `ldsc:genome_build` | string | genome build of the sidecar coordinates |
| `ldsc:sorted_by_build` | string | build defining the `IDX`/position sort order (`"hg19"` or `"hg38"`); required |
| `ldsc:row_group_size` | string(int) | intended row-group size (informational) |
| `ldsc:n_samples` | string(int) | reference individuals used to compute R² |
| `ldsc:r2_bias` | string | `"unbiased"` for package panels; `"raw"` reserved for external raw-R² panels needing correction |
| `ldsc:min_r2` | string(float) | unbiased-R² floor applied upstream; `"0.0"` means complete |
| `ldsc:n_snps` | string(int) | length of the panel sidecar = size of the index space |
| `ldsc:sidecar_identity_sha256` | string(hex) | binding hash of the sidecar identity (§4) |
| `ldsc:r2_encoding` | string | R² storage encoding; `"int16_symmetric"` for package panels. Absent on legacy float32 panels. |
| `ldsc:r2_scale` | string(int) | dequant divisor for an integer R² column (`"32767"`). The reader divides the int16 values by this to recover float32 R². |

Stored as UTF-8 byte strings, accessed via
`pq.ParquetFile(path).schema_arrow.metadata`.

### 2.6 Canonical Example

**Schema metadata:**
```
ldsc:schema_version          = "1"
ldsc:artifact_type           = "ref_panel_r2"
ldsc:snp_identifier          = "chr_pos"
ldsc:genome_build            = "hg19"
ldsc:sorted_by_build         = "hg19"
ldsc:n_samples               = "3202"
ldsc:r2_bias                 = "unbiased"
ldsc:min_r2                  = "0.0"
ldsc:n_snps                  = "969730"
ldsc:sidecar_identity_sha256 = "a1b2…(64 hex chars)"
ldsc:r2_encoding             = "int16_symmetric"
ldsc:r2_scale                = "32767"
```

**Data rows (sorted by non-decreasing `IDX_1`):** R² is stored on disk as int16
(`round(R²·32767)`); the decoded float32 value the reader yields is shown in
parentheses.

| IDX_1 | IDX_2 | R2 (int16 → decoded) | SIGN |
|---|---|---|---|
| 0 | 1 | 26915 → 0.8214 | True |
| 0 | 2 | 4397 → 0.1342 | False |
| 1 | 2 | 2920 → 0.0891 | True |

`IDX_1 = 0` is sidecar row 0; the left SNP appears in multiple consecutive rows,
one per right-side neighbor within the LD window.

### 2.7 SNP Metadata Sidecar (mandatory, authoritative)

The sidecar (`chrN_meta.tsv.gz`) is the per-SNP table for the whole panel and
**defines the index space**: parquet `IDX_k` is a 0-based row number into it. It
is a gzip TSV with `# ldsc:*` identity-provenance comments followed by columns
`CHR, POS, SNP, A1, A2, CM, MAF` (one row per SNP, in build/position order).

The sidecar is **required**. Because the parquet stores only indices, it is
meaningless without the exact sidecar; there is no synthesize-from-parquet
fallback. A missing sidecar is a hard error.

The sidecar feeds the partitioned-LDSC workflow as before: SNP universe `A`,
all/common counts via `MAF`, and the genetic-map (`CM`) window.

---

## 3. Read Pipeline

Implemented in `SortedR2BlockReader` (`src/ldsc/_kernel/ldscore.py`), the single
point of contact between the LD-score algorithm and the parquet.

### 3.1 Open, Bind, Remap (`__init__` → `_init_index_path`)

The reader receives the **retained** metadata (the analysis matrix universe, in
matrix order) exactly as the rest of the pipeline produces it — the index format
does not change any filtering/restriction/intersection step. It then:

1. Confirms the schema is the index layout (`IDX_1/IDX_2/R2` present); any other
   schema raises an actionable "regenerate with `ldsc build-ref-panel`" error.
2. Reads `ldsc:sorted_by_build` (required) and checks it against the analysis
   build.
3. Loads the **full** panel sidecar (`chrN_meta.tsv.gz`) — all `n_snps` rows in
   build order.
4. Validates the binding (§4): recomputes the identity hash and checks
   `len(sidecar) == n_snps`; mismatch is a hard error.
5. Builds `remap: build_idx → retained_idx` (int32, `-1` for SNPs not in the
   analysis universe) and `retained_build_idx: matrix_idx → build_idx`
   (ascending, for pruning) via `build_index_remap`. Matching uses the same
   mode-dependent identity key (`effective_merge_key_series`) that the matrix
   universe is keyed by, so the result is identical to the legacy per-pair
   identity lookup. **This is the only place identifier-mode matching happens —
   once per chromosome.** The remap **must be injective over retained indices**:
   if two distinct panel sidecar rows share the same identity key under the
   active `snp_identifier` and both map onto the analysis universe, they collapse
   onto one retained SNP and would silently double-count that SNP's pairs. Any
   such collapse is a **hard error** (rebuild the panel or use an allele-aware
   `--snp-identifier`). Because the check is the injectivity of a remap built
   from mode-derived keys (`base_key [+ ":" + normalized allele set]`), it is
   automatically mode-aware. Two panel rows collapse only when they share the
   full key for the active mode: `rsid` → same **rsID**; `rsid_allele_aware` →
   same **rsID + normalized allele pair**; `chr_pos` → same **CHR:POS**;
   `chr_pos_allele_aware` → same **CHR:POS + normalized allele pair**. So two
   rows at the same `CHR:POS` with different rsIDs stay distinct in both
   rsid-family modes, and a multi-allelic site splits into distinct SNPs in both
   allele-aware modes (the allele set is order-independent, so `A/G` ≡ `G/A`).
   This up-front rejection is what
   lets §3.4 accumulate LD scores with no pair deduplication; with each
   `(IDX_1, IDX_2)` pair unique (§2.2) and the remap injective, every decoded
   retained pair is necessarily unique.
6. No window index is built: LD-score computation streams every row group once
   (`iter_all_pairs`, §3.3). The coarse-row-group startup warning from the parquet
   footer statistics is retained.

### 3.2 Decode (`_decode_index_row_group`)

Per row group, the reader reads only `IDX_1, IDX_2, R2` and gathers:

```python
i = remap[IDX_1]; j = remap[IDX_2]
keep = (i >= 0) & (j >= 0)
```

Endpoints not in the analysis universe map to `-1` and are dropped. When the R²
column is an integer dtype (quantized panels), the reader **dequantizes**
(`r2 = int16_values / ldsc:r2_scale`, scale resolved once at open in
`_resolve_r2_scale`) to `float32` before the optional raw→unbiased R² correction
(`_transform_r2`); a float32 R² column (legacy) skips dequant. `SIGN` is not read
(unused by LD-score computation). Each decoded group is a numeric
`(i:int32, j:int32, r2:float32)` triple of retained-matrix index pairs.

### 3.3 Pair streaming (`iter_all_pairs`)

LD-score computation does **not** query windows or tile the panel into blocks. The
reader exposes `iter_all_pairs()`, which decodes **every row group exactly once**
in `IDX_1` order and yields its `(i, j, r2)` arrays. There is no row-group pruning,
no decoded-row-group cache, and no dense block matrix: each stored pair is produced
a single time and consumed immediately (§3.4). Read-side peak RSS is therefore
bounded by the accumulator (§3.6), not by the LD window.

### 3.4 LD-score accumulation (`ld_score_streaming_from_r2_reader`)

The LD score is `cor_sum = R · annot`, where `R` is the `m×m` symmetric
within-window R² matrix with unit diagonal. The streaming driver computes it
directly from the stored pairs:

```
cor_sum[k]  ← annot[k]                         # diagonal R²(k,k) = 1, every SNP
for each stored pair (i, j, r2) with i < j and i ≥ block_left[j]:
    cor_sum[i] ← cor_sum[i] + r2 · annot[j]    # off-diagonal, both orientations
    cor_sum[j] ← cor_sum[j] + r2 · annot[i]
```

`block_left[j]` applies the ldscore LD window (`--ld-wind-*`), which may be
narrower than the panel's build window: a stored pair contributes iff
`i ≥ block_left[j]`. The two scatter-adds use `np.add.at` (unbuffered) because a
row group repeats each left index `i` across its many right neighbours; a buffered
`+=` would silently collapse them. Cost is `O(nnz · n_a)` — each stored pair is
touched once per annotation. `annot` carries the partitioned annotation columns
plus the regression-weight column, so a single pass yields both the reference and
regression-universe LD scores.

`cor_sum` accumulates in `float64` from `float32` products, so results are *not
bit-identical* to the legacy genotype-block path (≈1e-6 drift from accumulation
order), well within LD-score tolerances.

The full step-by-step math, with a glossary defining every symbol, is in
`docs/current/ldscore-parquet-accumulation.md`.

### 3.5 Mode-agnosticism

Because the parquet carries no identity, the four identifier modes read the same
bytes; only `remap` construction (step 5) differs. The cross-mode parity test
builds one index parquet and asserts identical LD scores in all four modes.

### 3.6 Read-side memory (streaming)

The streaming reader holds no decoded-row-group cache and no dense block matrices.
Peak reader RSS is bounded by the accumulator, plus one decoded row group at a
time and the workflow/import floor:

```
RSS ≈ m · n_a · 8 bytes  (cor_sum, float64)  +  one decoded row group  +  fixed floor
```

For chr6 (`m ≈ 664k`, `n_a ≈ 57`): `cor_sum ≈ 303 MiB`; one row group adds a few
MiB. Two consequences vs. the prior block reader:

- **`--ld-wind-*` is not a read-side peak-RSS lever.** The window only filters
  which streamed pairs contribute (`i ≥ block_left[j]`); it does not change
  resident memory.
- **Dense wide-window regions (the chr6 MHC) no longer set a memory high-water
  mark** — they add pairs to stream, not memory to hold.

This matches the *builder*, whose RSS is also flat in window size (it streams a
bounded genotype window from disk — see
`docs/current/ld-window-parquet-r2-sidecar-behavior.md`). Under cross-chromosome
parallelism each worker holds its own `cor_sum`, so aggregate RSS ≈
workers × (largest concurrent chromosome's `cor_sum` + fixed floor).

---

## 4. Index ↔ Sidecar Binding

The binding breaks if the panel's SNP **set or order** changes, so the hash
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

**Sidecar mandatory.** The parquet is opaque without its exact sidecar; ship them
together. There is no metadata-synthesis fallback.

**Build specificity.** The parquet encodes one build's sort order. The 3-tier
inference of the old format is gone — the index format always records
`ldsc:sorted_by_build`, and its absence is a hard error. Maintain separate hg19
and hg38 panels.

**One chromosome per file.** The reader is instantiated per chromosome and
expects a single-chromosome parquet + sidecar.

**Clean break — index only.** The legacy 10-column canonical schema, the raw
external-input read path, and the synthesize-from-parquet metadata fallback are
removed. External R² formats (e.g. TOP-LD pairwise tables, SBayesR binary LD
matrices) are not read directly; a future converter would turn them into the
canonical index format offline. Existing 10-column panels must be regenerated.

---

## 6. Affected Modules

| Module | Role |
|---|---|
| `_kernel/ref_panel_builder.py` | `write_r2_parquet`: 4-column index schema, sort assertion on `IDX_1`, binding metadata; `sidecar_identity_sha256` (re-exported from `snp_identity`). Pairs flow as columnar `PairColumns` `(i, j, r2, sign)` batches from `yield_pairwise_r2_rows` straight into the writer — no per-pair dict — and int16 quantization/encodings are applied at table build; the on-disk format is unchanged. |
| `_kernel/snp_identity.py` | `sidecar_identity_sha256` binding hash |
| `_kernel/plink_bed.py` | `PlinkBEDFile` genotype source for the build: selective per-SNP read (restricted) or disk streaming (unrestricted), feeding standardized columns to pairwise-R2 emission without loading the whole chromosome |
| `ref_panel_builder.py` | build loop: sidecar built first, hash + `n_snps` passed to the writer |
| `_kernel/ldscore.py` | `SortedR2BlockReader` index path: full-sidecar load, binding validation, `build_index_remap`, gather decode, `iter_all_pairs` streaming; `ld_score_streaming_from_r2_reader` pair accumulation; raw/legacy and block/query/cache paths removed |
| `_kernel/ref_panel.py` | sidecar mandatory in `load_metadata` (synthesis fallback removed) |
| `tests/` | index writer/reader/binding/remap tests; cross-mode parity gate; build→read parity |
