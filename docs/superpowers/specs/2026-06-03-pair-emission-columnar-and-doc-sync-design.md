# Columnar pair emission (O6) + memory doc sync (O7) — design

Date: 2026-06-03
Status: proposed (awaiting approval)

Builds on the reader streaming work (`2026-06-02-build-ref-panel-reader-streaming-design.md`,
commits `f9d3c7b`…`ce341b6`). O4 (cM map) and O5 (`.bim` dtypes) are explicitly dropped.

## O6 — columnar pair emission

### Problem

The `build-ref-panel` pairwise path is already columnar at both ends but
round-trips through per-pair dictionaries in the middle:

1. `_emit_*` / `_finalize_pair_columns` produce `PairColumns` —
   `(i, j, r2:float32, sign:int8)` numpy arrays (`ref_panel_builder.py:283-334`).
2. `_stash_pair_rows` groups them by left index into `_PendingPairs` (lists of
   arrays) (`:337`).
3. `_pop_pair_rows_before` concatenates each bucket and **emits one `dict[str,
   float|int|str]` per pair**, encoding sign as the string `"+"`/`"-"` (`:351`).
4. `write_r2_parquet` accumulates those dicts into a 100k list, and
   `_standard_r2_index_table` rebuilds arrays with `np.fromiter((int(row["i"]) …))`
   and `row["sign"] == "+"` (`:561-584`).

So pairs go `arrays → dict-per-pair → arrays`. The dict layer costs ~31 MB for a
100k flush batch (vs ~1.2 MB columnar), allocates one dict per pair (millions to
billions per genome-wide build), and is weakly typed.

### Goal

Keep pairs columnar end-to-end: the generator yields `PairColumns` chunks and the
writer consumes them directly. Remove the dict carrier, the `"+"`/`"-"` string
round-trip, and the `np.fromiter`-from-dict reconstruction. **No change to the
output parquet** — bytes, schema, int16 quantization, encodings, and metadata are
identical.

Non-goals: changing the int16/`BYTE_STREAM_SPLIT`/`DELTA_BINARY_PACKED`/zstd
encoding, the `R2_QUANT_SCALE`/`R2_ENCODING` constants, `_quantize_r2`,
`_unbiased_r2_array`, the row-group/batch sizes, or the reader.

### Design

The pair carrier becomes `PairColumns` (already defined) at every hop:

- **`_pop_pair_rows_before` yields one `PairColumns` chunk per finalized left
  index** instead of per-pair dicts. Each chunk is `(i_arr, j_arr, r2_arr,
  sign_arr)` for that left index, with `j` sorted ascending (today's order).
  Because chunks are emitted in increasing left-index order, IDX_1 stays
  non-decreasing — the existing guard still holds.
- **`yield_pairwise_r2_rows` / `iter_pairwise_r2_rows`** return
  `Iterator[PairColumns]` / `list[PairColumns]`.
- **`write_r2_parquet(pair_chunks: Iterable[PairColumns], …)`** accumulates chunks
  until the buffered row count reaches `batch_size`, concatenates them, and flushes
  one table. The IDX_1-monotonic guard runs on the concatenated `i` array exactly
  as today.
- **`_standard_r2_index_table`** takes the four concatenated arrays directly:
  `i→int32`, `j→int32`, `_quantize_r2(r2)→int16`, `(sign == 1)→bool`. The `+1/-1`
  int8 sign maps to the same `SIGN` bool the dict path produced (`"+" → True`).

This is the minimal change that removes the round-trip; the flush/batching model,
ordering guard, and encoding are untouched.

### Contract change and test impact

`write_r2_parquet` / `_standard_r2_index_table` currently accept an iterable of
pair **dicts**, and that is a *tested contract*: `test_plink_io.py` builds dict
lists like `{"i":0,"j":1,"R2":0.4,"sign":"+"}` in ~15 cases, and
`test_ref_panel_builder.py` asserts `row["i"]/row["R2"]/row["sign"]` on
`iter_pairwise_r2_rows` output (`:266-289`, `:303…`). The clean refactor changes
that input contract to `PairColumns` chunks, so those tests must be rewritten to
build small columnar inputs (a one-line `_chunk(i, j, r2, sign)` test helper). No
dual dict/column input is kept — that would re-introduce the very layer we are
removing.

### Safety

- The **build-ref-panel parquet golden** from the streaming session
  (`tests/fixtures/golden/build_ref_panel_r2.npz`, asserted by
  `test_build_ref_panel_parquet_golden`) guards end-to-end byte-equivalence of the
  emitted R2 parquet through the refactor.
- A new unit test asserts the generator yields arrays (not dicts) and that a
  columnar `write_r2_parquet` reproduces the same table as the pre-refactor dict
  path for a hand-built input.
- `min_r2` filtering stays in `_finalize_pair_columns` (unchanged).

## O7 — memory documentation sync

Bring the design docs and tutorial in line with the implemented memory work
(O1/O2/O3 reader + O6 emission) and set correct user expectations about which
knobs affect peak memory.

### Content

- **Reader behavior:** `build-ref-panel` (and the shared `PlinkBEDFile`) no longer
  loads the whole chromosome. Restricted builds read only kept SNP blocks; the
  default unrestricted build streams a sliding window. Individual filtering no
  longer materializes a second bitarray.
- **Knob expectations:** `--snp-batch-size`, `--min-r2`, `--ld-wind-*`, and
  `--maf-min` affect speed and output size, not peak RSS — peak is governed by the
  genotype read (now bounded) plus the workflow/import floor.
- **Emission:** R2 pairs are produced and written as columnar batches (no per-pair
  dict); the parquet format, int16 quantization, and encodings are unchanged.

### Files to update (`docs/current/`)

- `parquet-r2-format-and-read-pipeline.md` — note columnar emission; confirm
  format/encoding unchanged.
- `ld-window-parquet-r2-sidecar-behavior.md` — reader read/stream behavior and the
  knob-vs-memory note for window options.
- `io-argument-inventory.md` — annotate the memory implication (or lack thereof) of
  the affected build-ref-panel arguments.
- `artifact-metadata-field-inventory.md` — confirm the R2 metadata fields are
  unchanged (no new fields from O6).
- `architecture.md`, `code-structure.md`, `data-flow.md`, `class-and-features.md` —
  reflect the new `_kernel/plink_bed.py` module, the selective/streaming read, and
  columnar emission where these docs describe the reader or pair path.

### Tutorial

- `tutorials/build-parquet-reference-panel-from-plink.md` — a short "memory" note:
  what drives peak, that restricted builds read selectively and unrestricted builds
  stream, and that the tuning knobs are not peak-memory levers.

### Also

- Update `design_map.md` to map this spec/plan and the reader docs to
  `_kernel/plink_bed.py` and the columnar emission functions.

## Implementation sequence

1. O6 columnar refactor (TDD; golden + new unit test green; full suite green).
2. O7 doc sync (no code).

Each step is its own commit.
