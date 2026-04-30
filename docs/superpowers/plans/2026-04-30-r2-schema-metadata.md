# R2 Schema Metadata Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store `ldsc:n_samples` and `ldsc:r2_bias` in every R2 parquet file's Arrow schema metadata, and auto-load them at both downstream read sites so users never need to manually supply `--r2-bias-mode` or `--r2-sample-size` for panels built by this codebase.

**Architecture:** `write_r2_parquet` gains an `n_samples: int` parameter; `_build_chromosome` passes `geno.n`. A new helper `_read_r2_schema_meta()` and resolver `_resolve_r2_bias_from_meta()` in `_kernel/ref_panel.py` are called at two sites: the CLI validation block in `_kernel/ldscore.py` and `ParquetR2RefPanel.load_r2()`. `SortedR2BlockReader` is unchanged.

**Tech Stack:** Python 3.11, pyarrow, pandas, unittest

---

## File Map

| File | Change |
|---|---|
| `src/ldsc/_kernel/ref_panel_builder.py` | Add `n_samples: int` to `write_r2_parquet`; write two new schema keys |
| `src/ldsc/ref_panel_builder.py` | Pass `geno.n` as `n_samples` to `write_r2_parquet` in `_build_chromosome` |
| `src/ldsc/_kernel/ref_panel.py` | Add `_R2SchemaMeta`, `_read_r2_schema_meta()`, `_resolve_r2_bias_from_meta()` |
| `src/ldsc/_kernel/ldscore.py` | Apply auto-load in the CLI `r2_table` validation block (lines ~1987–1993) |
| `tests/test_ref_panel_builder.py` | Test that `write_r2_parquet` stores both keys correctly |
| `tests/test_ldscore_workflow.py` | Unit tests for `_resolve_r2_bias_from_meta`; integration tests for both auto-load sites |

---

## Task 1: Store `ldsc:n_samples` and `ldsc:r2_bias` in `write_r2_parquet`

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py` — `write_r2_parquet()` signature and `pa_meta`
- Modify: `src/ldsc/ref_panel_builder.py` — `_build_chromosome()` call to `write_r2_parquet`
- Test: `tests/test_ref_panel_builder.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_ref_panel_builder.py` in a new class `R2SchemaMetadataTest`:

```python
class R2SchemaMetadataTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_write_r2_parquet_stores_n_samples_and_r2_bias(self):
        import pyarrow.parquet as pq
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            ref_table = kernel_builder.build_reference_snp_table(
                metadata=pd.DataFrame({
                    "CHR": ["1", "1"],
                    "SNP": ["rs1", "rs2"],
                    "A1": ["A", "T"],
                    "A2": ["C", "G"],
                    "MAF": [0.3, 0.4],
                }),
                hg19_positions=np.array([100, 200], dtype=int),
                hg38_positions=None,
            )
            kernel_builder.write_r2_parquet(
                pair_rows=iter([]),
                reference_snp_table=ref_table,
                path=path,
                genome_build="hg19",
                n_samples=42,
            )
            meta = pq.read_schema(str(path)).metadata
            self.assertEqual(meta[b"ldsc:n_samples"], b"42")
            self.assertEqual(meta[b"ldsc:r2_bias"], b"unbiased")
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::R2SchemaMetadataTest -v
```

Expected: `FAILED` — `write_r2_parquet` has no `n_samples` parameter.

- [ ] **Step 3: Add `n_samples` to `write_r2_parquet` in `src/ldsc/_kernel/ref_panel_builder.py`**

Change the function signature from:
```python
def write_r2_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    path: str | PathLike[str],
    genome_build: str,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
) -> str:
```
to:
```python
def write_r2_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    reference_snp_table: pd.DataFrame,
    path: str | PathLike[str],
    genome_build: str,
    n_samples: int,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
) -> str:
```

Change `pa_meta` (currently at line ~644) from:
```python
    pa_meta = {
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
    }
```
to:
```python
    pa_meta = {
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size":  str(row_group_size).encode("utf-8"),
        b"ldsc:n_samples":       str(n_samples).encode("utf-8"),
        b"ldsc:r2_bias":         b"unbiased",
    }
```

- [ ] **Step 4: Pass `geno.n` to `write_r2_parquet` in `src/ldsc/ref_panel_builder.py`**

In `_build_chromosome`, the call to `kernel_builder.write_r2_parquet(...)` (around line 468) does not yet pass `n_samples`. Add `n_samples=geno.n`:

```python
            kernel_builder.write_r2_parquet(
                pair_rows=kernel_builder.yield_pairwise_r2_rows(
                    block_left=block_left,
                    chunk_size=config.chunk_size,
                    standardized_snp_getter=geno.nextSNPs,
                    m=geno.m,
                    n=geno.n,
                ),
                reference_snp_table=reference_snp_table,
                path=r2_path,
                genome_build=build,
                n_samples=geno.n,
            )
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py::R2SchemaMetadataTest -v
```

Expected: PASS.

- [ ] **Step 6: Run full suite to check for regressions**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ref_panel_builder.py -v
```

Expected: all tests PASS. Any test that calls `write_r2_parquet` without `n_samples` will fail with a `TypeError` — fix those call sites in the test fixtures by adding `n_samples=10` (a plausible dummy value).

- [ ] **Step 7: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/_kernel/ref_panel_builder.py src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py && git commit -m "feat: store ldsc:n_samples and ldsc:r2_bias in R2 parquet schema metadata"
```

---

## Task 2: Add `_read_r2_schema_meta()` and `_resolve_r2_bias_from_meta()` helpers

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel.py` — add dataclass + two helper functions
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing tests**

Add a new class `R2SchemaMetaReaderTest` in `tests/test_ldscore_workflow.py`:

```python
class R2SchemaMetaReaderTest(unittest.TestCase):
    """Tests for _read_r2_schema_meta and _resolve_r2_bias_from_meta."""

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def _write_minimal_parquet(self, path, extra_meta=None):
        """Write a zero-row canonical parquet with optional schema metadata."""
        import pyarrow as pa
        import pyarrow.parquet as pq
        schema = pa.schema([
            pa.field("CHR", pa.string()),
            pa.field("POS_1", pa.int64()),
            pa.field("POS_2", pa.int64()),
            pa.field("R2", pa.float32()),
            pa.field("SNP_1", pa.string()),
            pa.field("SNP_2", pa.string()),
        ])
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"50000"}
        if extra_meta:
            meta.update(extra_meta)
        table = pa.table({col: pa.array([], type=schema.field(col).type) for col in schema.names})
        table = table.replace_schema_metadata(meta)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, str(path))

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_read_r2_schema_meta_returns_both_fields(self):
        from ldsc._kernel.ref_panel import _read_r2_schema_meta
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1_r2.parquet"
            self._write_minimal_parquet(path, {b"ldsc:n_samples": b"200", b"ldsc:r2_bias": b"unbiased"})
            meta = _read_r2_schema_meta(str(path))
        self.assertEqual(meta.n_samples, 200)
        self.assertEqual(meta.r2_bias, "unbiased")

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_read_r2_schema_meta_legacy_file_returns_none_fields(self):
        from ldsc._kernel.ref_panel import _read_r2_schema_meta
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1_r2.parquet"
            self._write_minimal_parquet(path)  # no ldsc:n_samples or ldsc:r2_bias
            meta = _read_r2_schema_meta(str(path))
        self.assertIsNone(meta.n_samples)
        self.assertIsNone(meta.r2_bias)

    def test_resolve_unbiased_no_user_n_returns_unbiased(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta
        meta = _R2SchemaMeta(n_samples=100, r2_bias="unbiased")
        mode, n = _resolve_r2_bias_from_meta(None, None, meta)
        self.assertEqual(mode, "unbiased")
        self.assertIsNone(n)

    def test_resolve_unbiased_with_user_n_warns_and_returns_unbiased(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta
        meta = _R2SchemaMeta(n_samples=100, r2_bias="unbiased")
        import logging
        with self.assertLogs("LDSC", level="WARNING"):
            mode, n = _resolve_r2_bias_from_meta(None, 999.0, meta)
        self.assertEqual(mode, "unbiased")
        self.assertIsNone(n)

    def test_resolve_raw_no_user_n_autofills_stored_n(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta
        meta = _R2SchemaMeta(n_samples=150, r2_bias="raw")
        mode, n = _resolve_r2_bias_from_meta(None, None, meta)
        self.assertEqual(mode, "raw")
        self.assertEqual(n, 150.0)

    def test_resolve_raw_user_n_overrides_stored_n(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta
        meta = _R2SchemaMeta(n_samples=150, r2_bias="raw")
        mode, n = _resolve_r2_bias_from_meta(None, 999.0, meta)
        self.assertEqual(mode, "raw")
        self.assertEqual(n, 999.0)

    def test_resolve_absent_metadata_defaults_to_unbiased(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta
        meta = _R2SchemaMeta(n_samples=None, r2_bias=None)
        mode, n = _resolve_r2_bias_from_meta(None, None, meta)
        self.assertEqual(mode, "unbiased")
        self.assertIsNone(n)

    def test_resolve_absent_metadata_raw_mode_with_user_n(self):
        from ldsc._kernel.ref_panel import _R2SchemaMeta, _resolve_r2_bias_from_meta
        meta = _R2SchemaMeta(n_samples=None, r2_bias=None)
        mode, n = _resolve_r2_bias_from_meta("raw", 300.0, meta)
        self.assertEqual(mode, "raw")
        self.assertEqual(n, 300.0)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ldscore_workflow.py::R2SchemaMetaReaderTest -v
```

Expected: `ImportError` — `_read_r2_schema_meta` does not exist yet.

- [ ] **Step 3: Implement the helpers in `src/ldsc/_kernel/ref_panel.py`**

Add near the top of the module (after existing imports):

```python
from dataclasses import dataclass as _dataclass

@_dataclass
class _R2SchemaMeta:
    n_samples: int | None
    r2_bias: str | None   # "unbiased" | "raw" | None
```

Add as module-level functions (before `ParquetR2RefPanel`):

```python
def _read_r2_schema_meta(path: str) -> "_R2SchemaMeta":
    """Read ldsc:n_samples and ldsc:r2_bias from a canonical R2 parquet file."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return _R2SchemaMeta(n_samples=None, r2_bias=None)

    raw_meta = pq.read_schema(path).metadata or {}
    n_raw = raw_meta.get(b"ldsc:n_samples")
    bias_raw = raw_meta.get(b"ldsc:r2_bias")

    n_samples = int(n_raw.decode("utf-8")) if n_raw is not None else None
    r2_bias = bias_raw.decode("utf-8") if bias_raw is not None else None

    if n_samples is not None and r2_bias is None:
        import logging
        logging.getLogger("LDSC").warning(
            f"'{path}' has ldsc:n_samples but no ldsc:r2_bias — treating as 'raw'."
        )
        r2_bias = "raw"

    return _R2SchemaMeta(n_samples=n_samples, r2_bias=r2_bias)


def _resolve_r2_bias_from_meta(
    r2_bias_mode: str | None,
    r2_sample_size: float | None,
    meta: "_R2SchemaMeta",
) -> tuple[str, float | None]:
    """Resolve r2_bias_mode and r2_sample_size using stored schema metadata.

    Returns (resolved_bias_mode, resolved_sample_size).
    """
    import logging
    logger = logging.getLogger("LDSC")

    stored_bias = meta.r2_bias
    stored_n = float(meta.n_samples) if meta.n_samples is not None else None

    # Determine effective bias mode
    effective_bias = r2_bias_mode if r2_bias_mode is not None else (stored_bias or "unbiased")

    if effective_bias == "unbiased":
        if r2_sample_size is not None:
            logger.warning(
                "r2_sample_size is ignored because R2 values are already unbiased "
                "(ldsc:r2_bias=unbiased in parquet schema metadata)."
            )
        return "unbiased", None

    # effective_bias == "raw"
    resolved_n = r2_sample_size if r2_sample_size is not None else stored_n
    return "raw", resolved_n
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ldscore_workflow.py::R2SchemaMetaReaderTest -v
```

Expected: all tests PASS.

- [ ] **Step 5: Run full suite**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/ -v
```

Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/_kernel/ref_panel.py tests/test_ldscore_workflow.py && git commit -m "feat: add _read_r2_schema_meta and _resolve_r2_bias_from_meta helpers"
```

---

## Task 3: Apply auto-load at the CLI validation site in `_kernel/ldscore.py`

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` — CLI `r2_table` validation block (~lines 1987–1993)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing tests**

Add to class `LDScoreParserTest` (or a new `R2AutoLoadCLITest`) in `tests/test_ldscore_workflow.py`.
These tests use the `_write_minimal_parquet` helper from Task 2's test class
(extract it to a module-level helper so both classes can use it):

```python
class R2AutoLoadCLITest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def _write_minimal_parquet(self, path, extra_meta=None):
        import pyarrow as pa, pyarrow.parquet as pq
        schema = pa.schema([
            pa.field("CHR", pa.string()), pa.field("POS_1", pa.int64()),
            pa.field("POS_2", pa.int64()), pa.field("R2", pa.float32()),
            pa.field("SNP_1", pa.string()), pa.field("SNP_2", pa.string()),
        ])
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"50000"}
        if extra_meta:
            meta.update(extra_meta)
        table = pa.table({col: pa.array([], type=schema.field(col).type) for col in schema.names})
        table = table.replace_schema_metadata(meta)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, str(path))

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_cli_autofills_unbiased_from_schema_when_mode_is_none(self):
        from ldsc._kernel import ldscore as kernel_ldscore
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            self._write_minimal_parquet(path, {
                b"ldsc:n_samples": b"200", b"ldsc:r2_bias": b"unbiased",
            })
            args = kernel_ldscore.build_parser().parse_args([
                "--r2-table", str(Path(tmpdir) / "hg19"),
                "--snp-identifier", "rsid",
                "--baseline-annot", "fake",
                "--out", "fake",
                "--ld-wind-kb", "1",
            ])
            args.r2_bias_mode = None
            args.r2_sample_size = None
            kernel_ldscore.validate_run_args(args, keep=False)
        self.assertEqual(args.r2_bias_mode, "unbiased")
        self.assertIsNone(args.r2_sample_size)

    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_cli_autofills_raw_and_n_from_schema(self):
        from ldsc._kernel import ldscore as kernel_ldscore
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            self._write_minimal_parquet(path, {
                b"ldsc:n_samples": b"150", b"ldsc:r2_bias": b"raw",
            })
            args = kernel_ldscore.build_parser().parse_args([
                "--r2-table", str(Path(tmpdir) / "hg19"),
                "--snp-identifier", "rsid",
                "--baseline-annot", "fake",
                "--out", "fake",
                "--ld-wind-kb", "1",
            ])
            args.r2_bias_mode = None
            args.r2_sample_size = None
            kernel_ldscore.validate_run_args(args, keep=False)
        self.assertEqual(args.r2_bias_mode, "raw")
        self.assertAlmostEqual(args.r2_sample_size, 150.0)
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ldscore_workflow.py::R2AutoLoadCLITest -v
```

Expected: `FAILED` — auto-load not yet applied.

- [ ] **Step 3: Apply auto-load in the `r2_table` validation block in `_kernel/ldscore.py`**

Add an import at the top of `_kernel/ldscore.py` (with other kernel imports):

```python
from . import ref_panel as _ref_panel_mod
```

Replace the current validation block (lines ~1987–1993):

```python
    if args.r2_table:
        if keep:
            raise ValueError("--keep-indivs-file is only supported in PLINK mode.")
        if args.r2_bias_mode is None:
            args.r2_bias_mode = "unbiased"
        if args.r2_bias_mode == "raw" and args.r2_sample_size is None:
            raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
```

with:

```python
    if args.r2_table:
        if keep:
            raise ValueError("--keep-indivs-file is only supported in PLINK mode.")
        _first_r2 = _first_r2_parquet_in_dir(args.r2_table)
        if _first_r2 is not None:
            _stored = _ref_panel_mod._read_r2_schema_meta(_first_r2)
            args.r2_bias_mode, args.r2_sample_size = _ref_panel_mod._resolve_r2_bias_from_meta(
                args.r2_bias_mode,
                getattr(args, "r2_sample_size", None),
                _stored,
            )
        if args.r2_bias_mode is None:
            args.r2_bias_mode = "unbiased"
        if args.r2_bias_mode == "raw" and args.r2_sample_size is None:
            raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
```

Also add the helper (before `validate_run_args`):

```python
def _first_r2_parquet_in_dir(r2_table: str) -> str | None:
    """Return the path to the first chr*_r2.parquet found in r2_table, or None."""
    import glob
    pattern = str(Path(r2_table) / "chr*_r2.parquet")
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ldscore_workflow.py::R2AutoLoadCLITest -v
```

Expected: all tests PASS.

- [ ] **Step 5: Run full suite**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/ -v
```

Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py && git commit -m "feat: auto-load r2_bias_mode and r2_sample_size from parquet schema at CLI"
```

---

## Task 4: Apply auto-load in `ParquetR2RefPanel.load_r2()`

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel.py` — `ParquetR2RefPanel.load_r2()`
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

Add to `R2AutoLoadCLITest` (or a new `R2AutoLoadPanelTest`):

```python
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_parquet_panel_autofills_raw_and_n_from_schema(self):
        from ldsc._kernel.ref_panel import ParquetR2RefPanel
        from ldsc.config import GlobalConfig, RefPanelSpec
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "hg19" / "chr1_r2.parquet"
            self._write_minimal_parquet(path, {
                b"ldsc:n_samples": b"200", b"ldsc:r2_bias": b"raw",
            })
            spec = RefPanelSpec(
                r2_dir=str(Path(tmpdir) / "hg19"),
                r2_bias_mode=None,
                r2_sample_size=None,
            )
            panel = ParquetR2RefPanel(
                spec=spec,
                global_config=GlobalConfig(snp_identifier="rsid", genome_build="hg19"),
            )
            metadata = pd.DataFrame({
                "CHR": ["1"], "POS": [100], "SNP": ["rs1"], "CM": [0.0], "MAF": [0.3],
            })
            reader = panel.load_r2("1", metadata=metadata)
        self.assertEqual(reader.r2_bias_mode, "raw")
        self.assertAlmostEqual(reader.r2_sample_size, 200.0)
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ldscore_workflow.py::R2AutoLoadCLITest::test_parquet_panel_autofills_raw_and_n_from_schema -v
```

Expected: `FAILED` — `reader.r2_bias_mode` is `None`, not `"raw"`.

- [ ] **Step 3: Apply auto-load in `ParquetR2RefPanel.load_r2()` in `src/ldsc/_kernel/ref_panel.py`**

Current `load_r2` method (around line 315):

```python
    def load_r2(
        self,
        chrom: str,
        *,
        metadata: pd.DataFrame | None = None,
        r2_bias_mode: str | None = None,
        r2_sample_size: float | None = None,
    ):
        metadata = metadata if metadata is not None else self.load_metadata(chrom)
        return kernel_ldscore.SortedR2BlockReader(
            paths=self.resolve_r2_paths(chrom),
            chrom=chrom,
            metadata=metadata,
            identifier_mode=self.global_config.snp_identifier,
            r2_bias_mode=self.spec.r2_bias_mode if r2_bias_mode is None else r2_bias_mode,
            r2_sample_size=r2_sample_size if r2_sample_size is not None else self.spec.sample_size,
            genome_build=self.global_config.genome_build,
        )
```

Replace with:

```python
    def load_r2(
        self,
        chrom: str,
        *,
        metadata: pd.DataFrame | None = None,
        r2_bias_mode: str | None = None,
        r2_sample_size: float | None = None,
    ):
        metadata = metadata if metadata is not None else self.load_metadata(chrom)
        paths = self.resolve_r2_paths(chrom)

        effective_bias = self.spec.r2_bias_mode if r2_bias_mode is None else r2_bias_mode
        effective_n = r2_sample_size if r2_sample_size is not None else self.spec.sample_size

        if paths:
            stored = _read_r2_schema_meta(paths[0])
            effective_bias, effective_n = _resolve_r2_bias_from_meta(
                effective_bias, effective_n, stored
            )

        return kernel_ldscore.SortedR2BlockReader(
            paths=paths,
            chrom=chrom,
            metadata=metadata,
            identifier_mode=self.global_config.snp_identifier,
            r2_bias_mode=effective_bias,
            r2_sample_size=effective_n,
            genome_build=self.global_config.genome_build,
        )
```

- [ ] **Step 4: Run test to verify it passes**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/test_ldscore_workflow.py::R2AutoLoadCLITest -v
```

Expected: all tests PASS.

- [ ] **Step 5: Run full suite**

```bash
cd ldsc_py3_Jerry && source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3 && pytest tests/ -v
```

Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
cd ldsc_py3_Jerry && git add src/ldsc/_kernel/ref_panel.py tests/test_ldscore_workflow.py && git commit -m "feat: auto-load r2_bias_mode and r2_sample_size from parquet schema in ParquetR2RefPanel"
```
