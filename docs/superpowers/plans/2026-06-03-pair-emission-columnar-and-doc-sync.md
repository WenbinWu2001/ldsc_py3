# Columnar pair emission (O6) + memory doc sync (O7) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Remove the per-pair dict round-trip in the `build-ref-panel` R2 path (pairs stay columnar from generator to parquet writer), with byte-identical parquet output; then sync the design docs and tutorial to the implemented memory work.

**Architecture:** The pair carrier becomes `PairColumns` `(i, j, r2:float32, sign:int8)` at every hop. `_pop_pair_rows_before` yields one chunk per finalized left index; `write_r2_parquet` accumulates chunks to `batch_size` rows, concatenates, and flushes. Encoding (int16/`BYTE_STREAM_SPLIT`/`DELTA_BINARY_PACKED`/zstd), the `IDX_1`-monotonic guard, and batch/row-group sizes are unchanged.

**Tech Stack:** Python 3.13, numpy, pyarrow; pytest.

**Reference spec:** `docs/superpowers/specs/2026-06-03-pair-emission-columnar-and-doc-sync-design.md`

**Environment for every command:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
```
Repo root: `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`.

**Safety gate:** `tests/test_plink_bed_reader.py::test_build_ref_panel_parquet_golden` must stay green — it asserts byte-equivalence of the emitted R2 parquet. If it fails, the refactor changed output; fix the code, not the golden.

---

## Task 1: Columnar production path

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py` (`_pop_pair_rows_before`, `yield_pairwise_r2_rows`, `iter_pairwise_r2_rows`, `_standard_r2_index_table`, `write_r2_parquet`)
- Modify: `src/ldsc/ref_panel_builder.py:870` (call-site kwarg rename)

- [ ] **Step 1: Baseline green**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: all pass (record count).

- [ ] **Step 2: `_pop_pair_rows_before` yields `PairColumns` chunks**

Replace the dict-yielding body (`ref_panel_builder.py:351-369`) with:
```python
def _pop_pair_rows_before(
    pending: dict[int, _PendingPairs],
    min_future_i: int,
) -> Iterator[PairColumns]:
    """Yield one column chunk per left index whose pairs are final.

    Each chunk is ``(i, j, r2, sign)`` for a single left index with ``j`` sorted
    ascending. Chunks are emitted in increasing left-index order, so concatenated
    ``i`` stays non-decreasing downstream.
    """
    flushable = sorted(i for i in pending if i < int(min_future_i))
    for i in flushable:
        bucket = pending.pop(i)
        j = np.concatenate(bucket.j)
        r2 = np.concatenate(bucket.r2)
        sign = np.concatenate(bucket.sign)
        order = np.argsort(j, kind="stable")
        j = j[order]
        i_arr = np.full(j.shape, i, dtype=np.int64)
        yield (i_arr, j, r2[order], sign[order])
```

- [ ] **Step 3: Update generator return-type annotations**

In `yield_pairwise_r2_rows` change the return annotation `-> Iterator[dict[str, float | int | str]]` to `-> Iterator[PairColumns]` (body unchanged — its three `yield from _pop_pair_rows_before(...)` now yield chunks). In `iter_pairwise_r2_rows` change `-> list[dict[str, float | int | str]]` to `-> list[PairColumns]`.

- [ ] **Step 4: `_standard_r2_index_table` takes arrays**

Replace (`ref_panel_builder.py:561-584`) with:
```python
def _standard_r2_index_table(pa, schema, *, i, j, r2, sign):
    """Build one 4-column index ``pyarrow.Table`` from concatenated pair columns.

    ``i``/``j`` are sidecar-row indices (int32), ``r2`` is unbiased float32 stored
    as symmetric int16 (scale ``R2_QUANT_SCALE``), and ``sign`` is int8 +1/-1
    stored as the ``SIGN`` bool (``True`` when r >= 0).
    """
    if i.size == 0:
        return schema.empty_table()
    return pa.Table.from_arrays(
        [
            pa.array(i.astype(np.int32), type=pa.int32()),
            pa.array(j.astype(np.int32), type=pa.int32()),
            pa.array(_quantize_r2(r2), type=pa.int16()),
            pa.array(sign == 1, type=pa.bool_()),
        ],
        schema=schema,
    )
```

- [ ] **Step 5: `write_r2_parquet` consumes chunks**

Change the signature `pair_rows: Iterable[dict[str, float | int | str]]` to
`pair_chunks: Iterable[PairColumns]`. Replace the buffer/`_flush`/consume block
(`ref_panel_builder.py:664-721`) with:
```python
    writer = None
    buf_i: list[np.ndarray] = []
    buf_j: list[np.ndarray] = []
    buf_r2: list[np.ndarray] = []
    buf_sign: list[np.ndarray] = []
    buffered = 0
    prev_idx1: int | None = None

    def _flush() -> None:
        nonlocal writer, prev_idx1, buf_i, buf_j, buf_r2, buf_sign, buffered
        i = np.concatenate(buf_i) if buf_i else np.empty(0, dtype=np.int64)
        j = np.concatenate(buf_j) if buf_j else np.empty(0, dtype=np.int64)
        r2 = np.concatenate(buf_r2) if buf_r2 else np.empty(0, dtype=np.float32)
        sign = np.concatenate(buf_sign) if buf_sign else np.empty(0, dtype=np.int8)
        table = _standard_r2_index_table(pa, schema, i=i, j=j, r2=r2, sign=sign)
        if i.size:
            sequence = i if prev_idx1 is None else np.concatenate(([prev_idx1], i))
            bad = np.flatnonzero(np.diff(sequence) < 0)
            if bad.size:
                k = int(bad[0])
                raise ValueError(
                    "Pairs must arrive in non-decreasing IDX_1 order. "
                    f"Received IDX_1={int(sequence[k + 1])} after IDX_1={int(sequence[k])}. "
                    "Verify that the reference panel builder traverses SNPs in ascending index order."
                )
            prev_idx1 = int(i[-1])
        if writer is None:
            # IDX_2: DELTA_BINARY_PACKED exploits sorted right-neighbors. R2:
            # BYTE_STREAM_SPLIT splits the int16 byte planes for zstd.
            # use_dictionary restricts auto-dictionary to IDX_1 only.
            enc_kwargs = dict(
                column_encoding={"IDX_2": "DELTA_BINARY_PACKED", "R2": "BYTE_STREAM_SPLIT"},
                use_dictionary=["IDX_1"],
            )
            if pa.Codec.is_available("zstd"):
                writer = pq.ParquetWriter(
                    str(path), schema, compression="zstd", compression_level=9, **enc_kwargs
                )
            else:
                warnings.warn(
                    "zstd codec is not available in this pyarrow build; "
                    "falling back to snappy compression. "
                    "Install pyarrow from conda-forge or PyPI to enable zstd.",
                    UserWarning,
                    stacklevel=2,
                )
                writer = pq.ParquetWriter(str(path), schema, compression="snappy", **enc_kwargs)
        writer.write_table(table, row_group_size=row_group_size)
        buf_i, buf_j, buf_r2, buf_sign = [], [], [], []
        buffered = 0

    try:
        for chunk in pair_chunks:
            ci, cj, cr2, csign = chunk
            if ci.size == 0:
                continue
            buf_i.append(ci)
            buf_j.append(cj)
            buf_r2.append(cr2)
            buf_sign.append(csign)
            buffered += int(ci.size)
            if buffered >= batch_size:
                _flush()
        if buffered or writer is None:
            _flush()
    finally:
        if writer is not None:
            writer.close()
    return str(path)
```
Update the function docstring's first input line to describe `pair_chunks` (an
iterable of `(i, j, r2, sign)` arrays) instead of pair dicts. Leave all metadata,
schema, and encoding code above this block unchanged.

- [ ] **Step 6: Rename the workflow call-site kwarg**

In `src/ldsc/ref_panel_builder.py` (~line 870), change
`pair_rows=kernel_builder.yield_pairwise_r2_rows(` to
`pair_chunks=kernel_builder.yield_pairwise_r2_rows(`.

- [ ] **Step 7: Run the golden (expect failure only in dict-based tests, not the golden)**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest tests/test_plink_bed_reader.py::test_build_ref_panel_parquet_golden -v
```
Expected: PASS (production path is internally consistent; output unchanged). Dict-based unit tests in `test_plink_io.py`/`test_ref_panel_builder.py` will fail until Task 2 — that is expected.

- [ ] **Step 8: Commit**
```bash
git add src/ldsc/_kernel/ref_panel_builder.py src/ldsc/ref_panel_builder.py
git commit -m "refactor(ref-panel): columnar pair emission, drop dict round-trip"
```

---

## Task 2: Migrate dict-based tests to columnar inputs

**Files:**
- Modify: `tests/test_plink_io.py`, `tests/test_ref_panel_builder.py`

- [ ] **Step 1: Add test helpers**

At module scope in `tests/test_plink_io.py` (and import/duplicate the same two
helpers in `tests/test_ref_panel_builder.py`):
```python
import numpy as np

def pair_chunk(rows):
    """Build one PairColumns chunk from legacy {i,j,R2,sign} dict rows."""
    i = np.array([r["i"] for r in rows], dtype=np.int64)
    j = np.array([r["j"] for r in rows], dtype=np.int64)
    r2 = np.array([r["R2"] for r in rows], dtype=np.float32)
    sign = np.array([1 if r["sign"] == "+" else -1 for r in rows], dtype=np.int8)
    return (i, j, r2, sign)

def flatten_pairs(chunks):
    """Flatten PairColumns chunks back to (i, j, r2, sign) dict rows for assertions."""
    out = []
    for ci, cj, cr2, csign in chunks:
        for k in range(int(ci.size)):
            out.append({"i": int(ci[k]), "j": int(cj[k]),
                        "R2": float(cr2[k]), "sign": "+" if csign[k] == 1 else "-"})
    return out
```

- [ ] **Step 2: Convert `write_r2_parquet`/`_standard_r2_index_table` call sites**

Transformation rule (apply to every `pair_rows=<dict-list>` call in
`test_plink_io.py` — the literals stay, only the wrapping changes):
- `write_r2_parquet(..., pair_rows=ROWS, ...)` → `write_r2_parquet(..., pair_chunks=[pair_chunk(ROWS)], ...)`
- `_write_index_panel(..., ROWS, ...)` helper: change its internal `write_r2_parquet` call the same way (one edit covers all tests that route through it).
- `_standard_r2_index_table(pa, schema, pair_rows=ROWS)` (`test_ref_panel_builder.py:3083`) → build arrays first:
  ```python
  i, j, r2, sign = pair_chunk(ROWS)
  table = kb._standard_r2_index_table(pa, schema, i=i, j=j, r2=r2, sign=sign)
  ```

- [ ] **Step 3: Convert `iter_pairwise_r2_rows` assertions**

In `test_ref_panel_builder.py` (`:266-289`, `:303…`, `:326…`, `:349…`) wrap the
result: `rows = flatten_pairs(kernel_builder.iter_pairwise_r2_rows(...))`. The
existing `row["i"]/row["j"]/row["R2"]/row["sign"]` assertions then work unchanged.

- [ ] **Step 4: Run the affected suites**

Run:
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev \
  && pytest tests/test_plink_io.py tests/test_ref_panel_builder.py tests/test_plink_bed_reader.py -q
```
Expected: all pass.

- [ ] **Step 5: Full suite**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q
```
Expected: all pass (same count as Task 1 Step 1).

- [ ] **Step 6: Commit**
```bash
git add tests/test_plink_io.py tests/test_ref_panel_builder.py
git commit -m "test(ref-panel): columnar pair_chunks inputs"
```

---

## Task 3: O7 — sync design docs and tutorial

**Files (docs only; no code):**
- Modify: `docs/current/parquet-r2-format-and-read-pipeline.md`, `docs/current/ld-window-parquet-r2-sidecar-behavior.md`, `docs/current/io-argument-inventory.md`, `docs/current/artifact-metadata-field-inventory.md`, `docs/current/architecture.md`, `docs/current/code-structure.md`, `docs/current/data-flow.md`, `docs/current/class-and-features.md`
- Modify: `tutorials/build-parquet-reference-panel-from-plink.md`
- Modify: `design_map.md`

- [ ] **Step 1: Reader behavior in the design docs**

In `architecture.md`, `code-structure.md`, `data-flow.md`, `class-and-features.md`:
add/adjust the PLINK genotype reader description to state it lives in
`src/ldsc/_kernel/plink_bed.py`, reads selectively (kept SNP blocks) for restricted
builds, streams a sliding window for unrestricted builds, and never materializes
the whole-chromosome bitarray; individual filtering is fused into the per-SNP read.
Cite `plink_bed.py` (`__read__`, `_source_snp_bits`, `__filter_snps_maf__`,
`nextSNPs`).

- [ ] **Step 2: Emission + format docs**

In `parquet-r2-format-and-read-pipeline.md` add a sentence that R2 pairs are
produced and written as columnar batches (`PairColumns`), with the parquet schema,
int16 quantization, `BYTE_STREAM_SPLIT`/`DELTA_BINARY_PACKED` encodings, and
metadata fields unchanged. In `artifact-metadata-field-inventory.md` confirm no R2
metadata fields changed.

- [ ] **Step 3: Knob-vs-memory note**

In `io-argument-inventory.md` and `ld-window-parquet-r2-sidecar-behavior.md`, note
that `--snp-batch-size`, `--min-r2`, `--ld-wind-*`, and `--maf-min` affect speed
and output size, not peak RSS; peak is governed by the (now bounded) genotype read
plus the workflow/import floor.

- [ ] **Step 4: Tutorial memory note**

In `tutorials/build-parquet-reference-panel-from-plink.md` add a short "Memory"
subsection: restricted builds read only kept SNP blocks; the default unrestricted
build streams a window; the tuning knobs above are not peak-memory levers.

- [ ] **Step 5: design_map.md**

Add rows mapping this spec/plan and the reader docs to `_kernel/plink_bed.py` and
the columnar emission functions in `_kernel/ref_panel_builder.py`.

- [ ] **Step 6: Sanity-check docs build/links and commit**

Verify no broken intra-doc references were introduced (grep for renamed symbols).
```bash
git add docs/current/ tutorials/build-parquet-reference-panel-from-plink.md design_map.md \
        docs/superpowers/specs/2026-06-03-pair-emission-columnar-and-doc-sync-design.md \
        docs/superpowers/plans/2026-06-03-pair-emission-columnar-and-doc-sync.md
git commit -m "docs: sync reader/emission memory behavior; knob expectations"
```

---

## Self-review notes (author)

- **Spec coverage:** O6 production → Task 1; test contract migration → Task 2;
  O7 doc sync → Task 3. Golden equivalence gate referenced in Task 1 Step 7 and
  Task 2 Steps 4-5.
- **Encoding preserved:** int16 `_quantize_r2`, `BYTE_STREAM_SPLIT`,
  `DELTA_BINARY_PACKED`, zstd-9, `R2_QUANT_SCALE`/`R2_ENCODING`, batch/row-group
  sizes — all untouched (only the row carrier changes).
- **Ordering invariant:** `IDX_1`-monotonic guard runs on concatenated `i` in
  `_flush`, same as before; chunks emitted in increasing left-index order.
- **Empty-output parity:** final `if buffered or writer is None: _flush()`
  reproduces today's empty-parquet-with-schema behavior.
- **Naming consistency:** `pair_chunks`, `PairColumns`, `pair_chunk`,
  `flatten_pairs` used consistently across tasks.
