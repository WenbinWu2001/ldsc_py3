# int16 R² Quantization + BYTE_STREAM_SPLIT Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store the index-format R² parquet column as symmetric int16 (scale 32767) with `BYTE_STREAM_SPLIT`, dequantizing transparently on read, cutting the R² parquet ~63% smaller with no material effect on LD scores.

**Architecture:** The writer clips unbiased R² at `1.0`, quantizes to int16 (`round(R²·32767)`), writes an `int16` `R2` column with `BYTE_STREAM_SPLIT` plus `ldsc:r2_encoding`/`ldsc:r2_scale` provenance metadata. The reader auto-detects an integer `R2` column, reads the scale from metadata, and dequantizes (`q/scale`) to float32 before the existing raw→unbiased transform. float32 panels keep their current path unchanged.

**Tech Stack:** Python, NumPy, PyArrow (≥21 already pinned — int16 BSS supported), pytest/unittest.

**This is refactor #11** of the build-ref-panel optimization spec. Design:
`docs/superpowers/specs/2026-06-01-build-ref-panel-memory-optimization-design.md` §`#11`.

**Worktree (all work here):**
`/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/.worktrees/int16-r2-quantization`
Branch: `feat/int16-r2-quantization` (off `restructure`).

**Test runner (run from the worktree root):**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
cd .../.worktrees/int16-r2-quantization && python -m pytest -q
```
Baseline before this work: **870 passed, 1 skipped.**

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `src/ldsc/_kernel/ref_panel_builder.py` | Writer: clip, quantize, schema, encoding, metadata | Modify `_unbiased_r2_array`, `_standard_r2_index_table`, `write_r2_parquet`; add `_quantize_r2`, constants |
| `src/ldsc/_kernel/ldscore.py` | Reader: detect int16, dequantize, raw-path clip | Modify `__init__`, `_init_index_path`, `_decode_index_row_group`, `_transform_r2`; add `_resolve_r2_scale` |
| `tests/test_ref_panel_builder.py` | Writer/table tests | Modify `IndexArrowTableTest`; add quantization tests |
| `tests/test_ldscore_workflow.py` | Reader/parity tests | Modify `IndexCrossModeParityTest`; add dequant + backward-compat tests |
| `docs/current/parquet-r2-format-and-read-pipeline.md` | Format spec §2.3 | Update encoding table + metadata keys |

**No `environment.yml` change** — `pyarrow>=21,<24` already covers int16 BYTE_STREAM_SPLIT.

---

## Task 1: Upper-clip unbiased R² at 1.0 (writer + reader-raw path)

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py:239-243` (`_unbiased_r2_array`)
- Modify: `src/ldsc/_kernel/ldscore.py:1379-1389` (`_transform_r2` raw branch)
- Test: `tests/test_ref_panel_builder.py` (new test in `PairwiseEmissionTest`), `tests/test_ldscore_workflow.py` (new `R2TransformClipTest`)

- [ ] **Step 1: Write the failing writer-clip test**

Add to `tests/test_ref_panel_builder.py` inside `PairwiseEmissionTest` (after line 338):

```python
    def test_unbiased_r2_array_upper_clips_at_one_keeps_negatives(self):
        import numpy as np
        from ldsc._kernel import ref_panel_builder as kb
        # corr=1.0 -> raw sq=1.0; float roundoff can give sq slightly >1.
        corr = np.array([1.0, np.nextafter(np.float64(1.0), 2.0), 0.0], dtype=np.float64)
        out = kb._unbiased_r2_array(corr, n_samples=100)
        self.assertLessEqual(float(out.max()), 1.0)            # never exceeds 1.0
        self.assertEqual(float(out[0]), 1.0)                   # perfect LD -> exactly 1.0
        self.assertLess(float(out[2]), 0.0)                    # corr=0 -> negative unbiased kept
```

- [ ] **Step 2: Run it to verify it fails**

Run: `python -m pytest tests/test_ref_panel_builder.py -k unbiased_r2_array_upper_clips -q`
Expected: FAIL — `out.max()` exceeds 1.0 (no clip yet) or `out[0] != 1.0`.

- [ ] **Step 3: Implement the writer clip**

In `src/ldsc/_kernel/ref_panel_builder.py`, replace `_unbiased_r2_array` (lines 239-243):

```python
def _unbiased_r2_array(correlation: np.ndarray, n_samples: int) -> np.ndarray:
    """Convert correlation coefficients to the unbiased :math:`R^2` estimate.

    Upper-clipped at ``1.0``: float roundoff at perfect LD can push raw
    :math:`r^2` a hair above 1, and the int16 quantization endpoint must map a
    true maximum of exactly ``1.0`` to ``32767``. Negative unbiased values are
    kept (not floored) so the bias correction stays visible downstream.
    """
    sq = correlation * correlation
    denom = n_samples - 2 if n_samples > 2 else n_samples
    return np.minimum(sq - (1.0 - sq) / denom, 1.0)
```

- [ ] **Step 4: Write the failing reader-raw-clip test**

Add to `tests/test_ldscore_workflow.py` (top-level class):

```python
class R2TransformClipTest(unittest.TestCase):
    def test_raw_transform_upper_clips_at_one(self):
        import numpy as np
        from ldsc._kernel.ldscore import SortedR2BlockReader
        reader = SortedR2BlockReader.__new__(SortedR2BlockReader)
        reader.r2_bias_mode = "raw"
        reader.r2_sample_size = 100.0
        # raw r2 = 1.0 -> corrected 1.0 - 0/98 = 1.0; values >1 must clip to 1.0.
        out = reader._transform_r2(np.array([1.0, 1.05, 0.0], dtype=np.float32))
        self.assertLessEqual(float(out.max()), 1.0)
        self.assertEqual(float(out[0]), 1.0)
```

- [ ] **Step 5: Run it to verify it fails**

Run: `python -m pytest tests/test_ldscore_workflow.py -k raw_transform_upper_clips -q`
Expected: FAIL — `out.max()` is 1.05-ish (raw correction does not clip yet).

- [ ] **Step 6: Implement the reader-raw clip**

In `src/ldsc/_kernel/ldscore.py`, `_transform_r2` (lines 1379-1389), add the clip after the correction:

```python
    def _transform_r2(self, values: np.ndarray) -> np.ndarray:
        """Apply the configured raw-to-unbiased R2 correction when required."""
        values = values.astype(np.float32, copy=False)
        if self.r2_bias_mode == "raw":
            if self.r2_sample_size is None:
                raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
            denom = self.r2_sample_size - 2
            if denom <= 0:
                raise ValueError("--r2-sample-size must be greater than 2 for raw R2 correction.")
            values = values - (1.0 - values) / denom
            # Share the writer's R2<=1 invariant: roundoff/raw inputs can exceed 1.
            values = np.minimum(values, np.float32(1.0))
        return values
```

- [ ] **Step 7: Run both tests to verify they pass**

Run: `python -m pytest tests/test_ref_panel_builder.py -k unbiased_r2_array_upper_clips tests/test_ldscore_workflow.py -k raw_transform_upper_clips -q`
Expected: PASS (2 passed).

- [ ] **Step 8: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py src/ldsc/_kernel/ldscore.py tests/test_ref_panel_builder.py tests/test_ldscore_workflow.py
git commit -m "feat(ref-panel): upper-clip unbiased R2 at 1.0

Clip the unbiased R2 estimate at 1.0 in the writer and in the reader's
raw-correction path so the int16 quantization endpoint maps a true max
of exactly 1.0. Negatives are kept, not floored."
```

---

## Task 2: int16 quantization — table builder, schema, encoding, metadata

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py` — add constants + `_quantize_r2` (near line 239); `_standard_r2_index_table` (lines 538-560); `write_r2_parquet` schema (633), metadata (612-628), encoding (661-663)
- Test: `tests/test_ref_panel_builder.py` — modify `IndexArrowTableTest` (line 3017); add tests in `StandardTableFormattingTest`/`IndexWriterTest`

- [ ] **Step 1: Write the failing endpoint-exactness + quantization tests**

Add a new class to `tests/test_ref_panel_builder.py`:

```python
class R2QuantizationTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_quantize_r2_endpoint_is_exact(self):
        import numpy as np
        from ldsc._kernel import ref_panel_builder as kb
        q = kb._quantize_r2(np.array([1.0, 1.000010, 0.0, -0.0003], dtype=np.float32))
        self.assertEqual(q.dtype, np.int16)
        self.assertEqual(int(q[0]), 32767)                 # 1.0 -> 32767
        self.assertEqual(int(q[1]), 32767)                 # >1 clips to 32767
        self.assertEqual(int(q[2]), 0)
        self.assertLess(int(q[3]), 0)                      # negative preserved
        # endpoint round-trips to EXACTLY 1.0
        self.assertEqual(float(np.float32(32767) / np.float32(kb.R2_QUANT_SCALE)), 1.0)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_quantize_r2_precision_within_half_step(self):
        import numpy as np
        from ldsc._kernel import ref_panel_builder as kb
        r2 = np.linspace(-0.0003, 1.0, 5000).astype(np.float32)
        decoded = kb._quantize_r2(r2).astype(np.float32) / np.float32(kb.R2_QUANT_SCALE)
        self.assertLessEqual(float(np.max(np.abs(decoded - r2))), 2e-5)
```

- [ ] **Step 2: Run to verify it fails**

Run: `python -m pytest tests/test_ref_panel_builder.py -k R2QuantizationTest -q`
Expected: FAIL — `kb._quantize_r2` / `kb.R2_QUANT_SCALE` do not exist.

- [ ] **Step 3: Add constants and the `_quantize_r2` helper**

In `src/ldsc/_kernel/ref_panel_builder.py`, add directly above `_unbiased_r2_array` (line 239):

```python
# int16 symmetric R2 quantization (refactor #11). Scale 32767 = int16 max, so
# the endpoint is exact: round(1.0 * 32767) = 32767, decoded 32767/32767 = 1.0.
R2_QUANT_SCALE = 32767
R2_ENCODING = "int16_symmetric"


def _quantize_r2(r2: np.ndarray) -> np.ndarray:
    """Symmetric int16 quantization of unbiased R2 (scale ``R2_QUANT_SCALE``).

    ``stored = clip(round(R2 * 32767), -32767, 32767)``. The clip is defensive:
    with the upstream ``min(R2, 1.0)`` and ``|negatives| ~ 1/N`` it never fires
    on real data. float64 multiply makes rounding deterministic at half-steps.
    """
    scaled = np.rint(r2.astype(np.float64) * R2_QUANT_SCALE)
    return np.clip(scaled, -R2_QUANT_SCALE, R2_QUANT_SCALE).astype(np.int16)
```

- [ ] **Step 4: Run to verify the helper tests pass**

Run: `python -m pytest tests/test_ref_panel_builder.py -k R2QuantizationTest -q`
Expected: PASS (2 passed).

- [ ] **Step 5: Update the failing `IndexArrowTableTest` for int16 R2**

The table builder will emit int16 R2, so its test's schema and assertions must change. Replace `IndexArrowTableTest.test_index_arrow_table_has_four_columns` (line ~3017) with:

```python
    def test_index_arrow_table_has_four_columns(self):
        import pyarrow as pa
        from ldsc._kernel import ref_panel_builder as kb

        schema = pa.schema([("IDX_1", pa.int32()), ("IDX_2", pa.int32()),
                            ("R2", pa.int16()), ("SIGN", pa.bool_())])
        rows = [
            {"i": 0, "j": 2, "R2": 0.5, "sign": "+"},
            {"i": 0, "j": 3, "R2": -0.01, "sign": "-"},
        ]
        table = kb._standard_r2_index_table(pa, schema, pair_rows=rows)
        self.assertEqual(table.schema.names, ["IDX_1", "IDX_2", "R2", "SIGN"])
        self.assertEqual(table.column("IDX_1").to_pylist(), [0, 0])
        self.assertEqual(table.column("IDX_2").to_pylist(), [2, 3])
        self.assertEqual(table.column("SIGN").to_pylist(), [True, False])
        self.assertEqual(table.schema.field("R2").type, pa.int16())
        # 0.5 -> round(0.5*32767)=16384 (banker's rounding to even); -0.01 -> -328
        self.assertEqual(table.column("R2").to_pylist(), [16384, -328])
```

- [ ] **Step 6: Run to verify it now fails for the right reason**

Run: `python -m pytest tests/test_ref_panel_builder.py -k test_index_arrow_table_has_four_columns -q`
Expected: FAIL — builder still emits float32 R2 (`pa.array(r2, type=pa.float32())`), so `from_arrays` against the int16 schema raises a type error.

- [ ] **Step 7: Quantize inside `_standard_r2_index_table`**

In `src/ldsc/_kernel/ref_panel_builder.py`, edit `_standard_r2_index_table` (lines 547-560). Replace the R2 array build and the `pa.array` for R2:

```python
    count = len(pair_rows)
    i = np.fromiter((int(row["i"]) for row in pair_rows), dtype=np.int32, count=count)
    j = np.fromiter((int(row["j"]) for row in pair_rows), dtype=np.int32, count=count)
    r2 = np.fromiter((float(row["R2"]) for row in pair_rows), dtype=np.float32, count=count)
    sign = np.fromiter((row["sign"] == "+" for row in pair_rows), dtype=bool, count=count)
    return pa.Table.from_arrays(
        [
            pa.array(i, type=pa.int32()),
            pa.array(j, type=pa.int32()),
            pa.array(_quantize_r2(r2), type=pa.int16()),
            pa.array(sign, type=pa.bool_()),
        ],
        schema=schema,
    )
```

Also update the function's docstring R2 line (line 542) to say "the unbiased ``R2`` (stored as symmetric int16, scale 32767)".

- [ ] **Step 8: Write the failing writer schema/metadata/encoding test**

Add to `tests/test_ref_panel_builder.py` `R2QuantizationTest`:

```python
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_write_r2_parquet_stores_int16_with_bss_and_scale_metadata(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        from ldsc._kernel import ref_panel_builder as kb
        if not pa.Codec.is_available("zstd"):
            self.skipTest("zstd unavailable")
        rows = [{"i": 0, "j": 1, "R2": 1.0, "sign": "+"},
                {"i": 0, "j": 2, "R2": -0.0003, "sign": "-"}]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chr1_r2.parquet"
            kb.write_r2_parquet(pair_rows=iter(rows), path=path, genome_build="hg19",
                                n_samples=100, snp_identifier="chr_pos", n_snps=3,
                                sidecar_identity_sha256="0" * 64)
            pf = pq.ParquetFile(path)
            self.assertEqual(pf.schema_arrow.field("R2").type, pa.int16())
            meta = pf.schema_arrow.metadata
            self.assertEqual(meta[b"ldsc:r2_encoding"], b"int16_symmetric")
            self.assertEqual(meta[b"ldsc:r2_scale"], b"32767")
            self.assertEqual(meta[b"ldsc:schema_version"], b"1")
            rg = pf.metadata.row_group(0)
            encs = {rg.column(c).path_in_schema: rg.column(c).encodings
                    for c in range(rg.num_columns)}
            self.assertIn("BYTE_STREAM_SPLIT", encs.get("R2", ()))
```

- [ ] **Step 9: Run to verify it fails**

Run: `python -m pytest tests/test_ref_panel_builder.py -k test_write_r2_parquet_stores_int16 -q`
Expected: FAIL — schema R2 is float32; metadata keys / BSS encoding absent.

- [ ] **Step 10: Change schema, metadata, and encoding in `write_r2_parquet`**

In `src/ldsc/_kernel/ref_panel_builder.py`:

(a) Schema (line 633): `("R2", pa.float32()),` → `("R2", pa.int16()),`

(b) Metadata dict (after line 627, inside `pa_meta`):

```python
        b"ldsc:sidecar_identity_sha256": str(sidecar_identity_sha256).encode("utf-8"),
        b"ldsc:r2_encoding": R2_ENCODING.encode("utf-8"),
        b"ldsc:r2_scale": str(R2_QUANT_SCALE).encode("utf-8"),
```

(c) Encoding kwargs (lines 661-663):

```python
            enc_kwargs = dict(
                column_encoding={"IDX_2": "DELTA_BINARY_PACKED", "R2": "BYTE_STREAM_SPLIT"},
                use_dictionary=["IDX_1"],
            )
```

Update the comment at line 657-660 to note R2 now uses BYTE_STREAM_SPLIT (int16 byte-plane split) alongside IDX_2's DELTA encoding.

- [ ] **Step 11: Run the writer + table tests**

Run: `python -m pytest tests/test_ref_panel_builder.py -k "R2QuantizationTest or test_index_arrow_table_has_four_columns or test_write_r2_parquet_uses_delta_encoding_for_idx2 or test_write_r2_parquet_uses_zstd_level_9" -q`
Expected: PASS. (`test_write_r2_parquet_uses_delta_encoding_for_idx2` still passes — it only asserts IDX_2 has DELTA and R2 lacks DELTA, both still true.)

- [ ] **Step 12: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): quantize R2 to int16 + BYTE_STREAM_SPLIT

Store the R2 column as symmetric int16 (scale 32767) with
BYTE_STREAM_SPLIT and ldsc:r2_encoding/ldsc:r2_scale provenance.
~63% smaller R2 parquet; endpoint 1.0 round-trips exactly."
```

---

## Task 3: Reader dequantization (auto-detect by column dtype)

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` — `__init__` (line ~1283), `_init_index_path` (line 1329), add `_resolve_r2_scale`, `_decode_index_row_group` (lines 1520-1523)
- Test: `tests/test_ldscore_workflow.py` — add `R2DequantizationTest`

- [ ] **Step 1: Write the failing round-trip + backward-compat tests**

Add to `tests/test_ldscore_workflow.py`. (`_write_index_sidecar` already exists in this file and is used by `IndexCrossModeParityTest`.)

```python
class R2DequantizationTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_reader_dequantizes_int16_panel_to_float32(self):
        import numpy as np, pandas as pd, tempfile
        from pathlib import Path
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel import ldscore as kernel_ldscore
        from ldsc._kernel.snp_identity import sidecar_identity_sha256
        panel = pd.DataFrame({"CHR": ["1"] * 3, "POS": [100, 120, 140],
                              "SNP": ["rs1", "rs2", "rs3"], "A1": ["A"] * 3, "A2": ["C"] * 3,
                              "CM": [0.0] * 3, "MAF": [0.3] * 3})
        pairs = [{"i": 0, "j": 1, "R2": 1.0, "sign": "+"},
                 {"i": 0, "j": 2, "R2": 0.2, "sign": "-"}]
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            kb.write_r2_parquet(pair_rows=pairs, path=r2, genome_build="hg19", n_samples=100,
                                snp_identifier="chr_pos", min_r2=0.0, n_snps=3,
                                sidecar_identity_sha256=sidecar_identity_sha256(panel))
            reader = kernel_ldscore.SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=panel.copy(),
                identifier_mode="chr_pos", r2_bias_mode="unbiased",
                r2_sample_size=None, genome_build="hg19")
            self.assertEqual(reader._r2_scale, 32767.0)
            decoded = reader._decode_index_row_group(0)
            self.assertEqual(decoded.r2.dtype, np.float32)
            # endpoint exact; 0.2 within half-step
            got = dict(zip(zip(decoded.i.tolist(), decoded.j.tolist()), decoded.r2.tolist()))
            self.assertEqual(got[(0, 1)], 1.0)
            self.assertAlmostEqual(got[(0, 2)], 0.2, delta=2e-5)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_reader_leaves_float32_panel_unscaled(self):
        # A hand-built float32 R2 index parquet must read with _r2_scale is None.
        import pyarrow as pa, pyarrow.parquet as pq, pandas as pd, tempfile
        from pathlib import Path
        from ldsc._kernel import ldscore as kernel_ldscore
        from ldsc._kernel.snp_identity import sidecar_identity_sha256
        panel = pd.DataFrame({"CHR": ["1", "1"], "POS": [100, 120], "SNP": ["rs1", "rs2"],
                              "A1": ["A", "A"], "A2": ["C", "G"], "CM": [0.0, 0.0], "MAF": [0.3, 0.3]})
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write_index_sidecar(tmp, panel)
            r2 = meta.with_name("chr1_r2.parquet")
            schema = pa.schema(
                [("IDX_1", pa.int32()), ("IDX_2", pa.int32()),
                 ("R2", pa.float32()), ("SIGN", pa.bool_())]
            ).with_metadata({
                b"ldsc:schema_version": b"1", b"ldsc:artifact_type": b"ref_panel_r2",
                b"ldsc:snp_identifier": b"chr_pos", b"ldsc:genome_build": b"hg19",
                b"ldsc:sorted_by_build": b"hg19", b"ldsc:n_samples": b"100",
                b"ldsc:r2_bias": b"unbiased", b"ldsc:min_r2": b"0.0", b"ldsc:n_snps": b"2",
                b"ldsc:sidecar_identity_sha256": sidecar_identity_sha256(panel).encode(),
                b"ldsc:row_group_size": b"50000",
            })
            tbl = pa.table({"IDX_1": pa.array([0], pa.int32()), "IDX_2": pa.array([1], pa.int32()),
                            "R2": pa.array([0.5], pa.float32()), "SIGN": pa.array([True], pa.bool_())},
                           schema=schema)
            pq.write_table(tbl, str(r2))
            reader = kernel_ldscore.SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=panel.copy(),
                identifier_mode="chr_pos", r2_bias_mode="unbiased",
                r2_sample_size=None, genome_build="hg19")
            self.assertIsNone(reader._r2_scale)
            decoded = reader._decode_index_row_group(0)
            self.assertAlmostEqual(float(decoded.r2[0]), 0.5, delta=1e-7)
```

- [ ] **Step 2: Run to verify they fail**

Run: `python -m pytest tests/test_ldscore_workflow.py -k R2DequantizationTest -q`
Expected: FAIL — `reader._r2_scale` attribute does not exist; int16 panel reads raw integer values (1.0 endpoint becomes 32767.0).

- [ ] **Step 3: Initialize `self._r2_scale` in `__init__`**

In `src/ldsc/_kernel/ldscore.py`, after line 1283 (`self._pf = None`):

```python
        self._pf = None
        self._r2_scale: float | None = None
```

- [ ] **Step 4: Detect scale in `_init_index_path` and add `_resolve_r2_scale`**

In `_init_index_path`, right after `schema_meta = self._pf.schema_arrow.metadata or {}` (line 1329):

```python
        schema_meta = self._pf.schema_arrow.metadata or {}
        self._r2_scale = self._resolve_r2_scale(schema_meta, path)
```

Add this method to the class (e.g. immediately before `_transform_r2`, line 1379):

```python
    def _resolve_r2_scale(self, schema_meta: dict, path: str) -> float | None:
        """Return the dequant scale when R2 is stored as quantized integers.

        Detection is by on-disk column dtype: an integer ``R2`` column is
        quantized and the scale comes from ``ldsc:r2_scale`` (defaulting to
        32767 with a warning if the key is absent). A float ``R2`` column is the
        legacy/unquantized path and returns ``None``.
        """
        import pyarrow as pa
        if self._pf is None or not pa.types.is_integer(self._pf.schema_arrow.field("R2").type):
            return None
        scale_raw = schema_meta.get(b"ldsc:r2_scale")
        if scale_raw is None:
            LOGGER.warning(
                f"'{path}' has an integer R2 column but no ldsc:r2_scale; defaulting to 32767."
            )
            return 32767.0
        return float(scale_raw.decode("utf-8"))
```

- [ ] **Step 5: Dequantize in `_decode_index_row_group`**

In `src/ldsc/_kernel/ldscore.py`, replace line 1523:

```python
        r2_raw = _arrow_column_to_numpy(table.column("R2")).astype(np.float32, copy=False)
        if self._r2_scale is not None:
            r2_raw = r2_raw / np.float32(self._r2_scale)
        r2 = self._transform_r2(r2_raw)
```

- [ ] **Step 6: Run to verify they pass**

Run: `python -m pytest tests/test_ldscore_workflow.py -k R2DequantizationTest -q`
Expected: PASS (2 passed).

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git commit -m "feat(ldscore): dequantize int16 R2 panels on read

Auto-detect an integer R2 column, read ldsc:r2_scale, and divide to
float32 before the raw->unbiased transform. float32 panels keep the
existing unscaled path."
```

---

## Task 4: Relax parity tests for the quantized path

**Files:**
- Modify: `tests/test_ldscore_workflow.py` — `IndexCrossModeParityTest` (line 3033)
- Verify: `tests/test_ref_panel_builder.py` — `ReferencePanelBuilderParityTest.test_hm3_chr22_subset_runs_direct_and_parquet_ldscore_paths` (line 2898) needs no change (smoke test asserts finiteness/order, not equality)

- [ ] **Step 1: Update `IndexCrossModeParityTest` — keep modes mutually exact, relax the oracle**

The four modes read the *same* quantized parquet, so they must still be **bit-identical to each other**; only the dense oracle (built from exact 0.4/0.2/0.6/0.5) needs the quantization tolerance. Replace the per-mode loop + assertion (lines ~3063-3075) with:

```python
            mode_scores = {}
            for mode in ("rsid", "chr_pos", "rsid_allele_aware", "chr_pos_allele_aware"):
                reader = kernel_ldscore.SortedR2BlockReader(
                    paths=[str(r2)], chrom="1",
                    metadata=panel[["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]].copy(),
                    identifier_mode=mode, r2_bias_mode="unbiased",
                    r2_sample_size=None, genome_build="hg19",
                )
                scores = kernel_ldscore.ld_score_var_blocks_from_r2_reader(
                    block_left=np.zeros(4, dtype=np.int64), snp_batch_size=2,
                    annot=np.ones((4, 1), dtype=np.float32), block_reader=reader,
                )
                mode_scores[mode] = scores[:, 0]
            # All four modes read one quantized parquet -> bit-identical to each other.
            first = mode_scores["rsid"]
            for mode, s in mode_scores.items():
                np.testing.assert_array_equal(s, first, err_msg=f"mode {mode} differs from rsid")
            # Dense oracle uses exact R2; int16 quantization adds <= a few half-steps.
            np.testing.assert_allclose(first, expected, rtol=0, atol=5e-5,
                                       err_msg="quantized LD scores drifted from dense oracle")
```

- [ ] **Step 2: Run the parity tests**

Run: `python -m pytest tests/test_ldscore_workflow.py -k IndexCrossModeParityTest tests/test_ref_panel_builder.py -k test_hm3_chr22_subset_runs_direct_and_parquet -q`
Expected: PASS. If the hm3 smoke test fails, do NOT loosen blindly — inspect with systematic-debugging (it should still produce finite scores with max>1 after dequant).

- [ ] **Step 3: Commit**

```bash
git add tests/test_ldscore_workflow.py
git commit -m "test(ldscore): relax cross-mode oracle to int16 tolerance

Modes still read one parquet and must be bit-identical to each other;
only the exact-R2 dense oracle gets the int16 half-step tolerance."
```

---

## Task 5: Documentation — format spec §2.3

**Files:**
- Modify: `docs/current/parquet-r2-format-and-read-pipeline.md` (§2.3 encoding table + metadata keys)

- [ ] **Step 1: Locate the encoding/metadata section**

Run: `grep -n "R2\|float32\|encoding\|BYTE_STREAM\|DELTA\|ldsc:r2\|2.3" docs/current/parquet-r2-format-and-read-pipeline.md | head -40`

- [ ] **Step 2: Update the R2 row and metadata list**

Change the R2 column description from `float32` to: `int16 (symmetric quantization, scale 32767, BYTE_STREAM_SPLIT)`. Add the two metadata keys to the `ldsc:*` table: `ldsc:r2_encoding = int16_symmetric`, `ldsc:r2_scale = 32767`. Document the read-path dequant rule (integer R2 column → divide by `ldsc:r2_scale`; float column → legacy path) and the endpoint-exactness invariant (`1.0 ↔ 32767`). State the lossy bound: per-pair `|error| ≤ 1.5e-5`, LD-score `≤ ~2e-3`, effectively lossless for LD-score regression.

- [ ] **Step 3: Commit**

```bash
git add docs/current/parquet-r2-format-and-read-pipeline.md
git commit -m "docs(ref-panel): document int16 R2 encoding and dequant"
```

---

## Task 6: Full-suite verification

- [ ] **Step 1: Run the entire suite**

Run: `python -m pytest -q`
Expected: all green. Baseline was 870 passed, 1 skipped; this work **adds** tests (new `R2QuantizationTest`, `R2DequantizationTest`, `R2TransformClipTest`, two clip tests) and modifies a few. Expect ~877+ passed, 1 skipped, **0 failed**. Paste the actual summary line.

- [ ] **Step 2: Manual smoke — build a tiny panel and confirm on-disk encoding**

Run (adjust the fixture path if needed):

```bash
python - <<'PY'
import pyarrow.parquet as pq, tempfile
from pathlib import Path
from ldsc._kernel import ref_panel_builder as kb
rows = [{"i":0,"j":1,"R2":1.0,"sign":"+"},{"i":0,"j":2,"R2":0.0173,"sign":"-"}]
p = Path(tempfile.mkdtemp())/"chr1_r2.parquet"
kb.write_r2_parquet(pair_rows=iter(rows), path=p, genome_build="hg19", n_samples=3202,
                    snp_identifier="chr_pos", n_snps=3, sidecar_identity_sha256="0"*64)
pf = pq.ParquetFile(p)
print("R2 dtype:", pf.schema_arrow.field("R2").type)
print("encoding:", pf.metadata.row_group(0).column(2).encodings)
print("scale meta:", pf.schema_arrow.metadata[b"ldsc:r2_scale"])
PY
```

Expected: `R2 dtype: int16`, encodings include `BYTE_STREAM_SPLIT`, `scale meta: b'32767'`.

- [ ] **Step 3: Update the design status table**

In `docs/superpowers/specs/2026-06-01-build-ref-panel-memory-optimization-design.md`, change refactor #11 status from `planned` to `implemented` with the merge commit once merged. (Do this on `restructure` after merge, not on the feature branch.)

- [ ] **Step 4: Finish the branch**

Use `superpowers:finishing-a-development-branch` to choose merge/PR. Suggested PR title: `feat(ref-panel): int16 R² quantization + BYTE_STREAM_SPLIT (#11)`.

---

## Self-Review Notes

- **Spec coverage:** clip (Task 1) ↔ §#11 "Endpoint invariant" + clip decision; quantize/schema/encoding/metadata (Task 2) ↔ §#11 "Writer"; reader auto-detect + dequant (Task 3) ↔ §#11 "Reader"; relaxed parity (Task 4) ↔ tests #14; docs (Task 5) ↔ §#11 "Pin"/format doc; endpoint exactness pinned by `==` in Task 2 Step 1 ↔ user requirement.
- **Type consistency:** `R2_QUANT_SCALE` (int 32767), `R2_ENCODING` ("int16_symmetric"), `_quantize_r2(np.ndarray)->int16 ndarray`, `_resolve_r2_scale(dict,str)->float|None`, `self._r2_scale: float|None` used consistently across writer/reader/tests.
- **No env change:** `pyarrow>=21,<24` already supports int16 BYTE_STREAM_SPLIT (verified on 23.0.1).
- **Backward compat:** float32 panels detected by non-integer R2 dtype → `_r2_scale=None` → unchanged path (Task 3 Step 1 backward-compat test).
