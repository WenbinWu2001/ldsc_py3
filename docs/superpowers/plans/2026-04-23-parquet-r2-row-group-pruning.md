# Parquet R2 Row-Group Pruning — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current `pyarrow.Dataset` full-chromosome scan in `SortedR2BlockReader` with a tabix-style row-group index (`pq.ParquetFile` + footer statistics), eliminating ~46 000 repeated full-chromosome reads during LD score computation and reducing per-chromosome runtime from ~minutes to ~seconds.

**Architecture:** On file open, read parquet footer statistics (per-row-group min/max of `pos_1`) into an in-memory list `_rg_bounds`. Per window query, intersect the requested genomic range against `_rg_bounds` and call `pf.read_row_groups(matched_idxs)` — only overlapping row groups are read. The writer is updated to produce a compact 6-column canonical schema with sort-invariant assertion and parquet schema metadata (`ldsc:sorted_by_build`, `ldsc:row_group_size`). Legacy raw-schema files continue to work via a `pyarrow.Dataset` fallback path with a deprecation warning.

**Tech Stack:** `pyarrow` ≥ 7 (`pq.ParquetFile`, `read_row_groups`, schema metadata), `numpy`, `pandas`. No new external dependencies.

**Design reference:** `docs/design/parquet-r2-format-and-read-pipeline.md` — read §2 (format) and §3 (read pipeline) before implementing. Terminology: "canonical schema" = 6-column format with `pos_1`/`pos_2`; "raw schema" = legacy 14-column format with `hg19_pos_1` etc.

---

## File Map

| File | Change |
|---|---|
| `src/ldsc/_kernel/ref_panel_builder.py` | Shrink `_STANDARD_LD_COLUMNS` to 6 cols; rewrite `build_standard_ld_table`; add `genome_build` + `row_group_size` + sort assertion + schema metadata to `write_standard_ld_parquet` |
| `src/ldsc/_kernel/ldscore.py` | Refactor `SortedR2BlockReader.__init__` (canonical: `pq.ParquetFile` + build validation + `_rg_bounds`; raw: `Dataset` + deprecation warning); rewrite `_query_union_rows` canonical path; add `get_present_identifiers()`; update `compute_chrom_from_parquet` to open reader first; remove `read_sorted_r2_presence` |
| `tests/test_ref_panel_builder.py` | Update `test_build_standard_ld_table_uses_exact_schema` for 6-column schema; add sort assertion + schema metadata tests |
| `tests/test_ldscore_legacy.py` | Update `RawParquetRuntimeTest` to remove `read_sorted_r2_presence` call; add `CanonicalParquetRuntimeTest` class for canonical schema |

---

## Task 1: Update parquet writer to 6-column canonical schema

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py:32-47` (`_STANDARD_LD_COLUMNS`)
- Modify: `src/ldsc/_kernel/ref_panel_builder.py:453-487` (`build_standard_ld_table`)
- Modify: `src/ldsc/_kernel/ref_panel_builder.py:528-572` (`write_standard_ld_parquet`)
- Test: `tests/test_ref_panel_builder.py`

### Context

Currently `_STANDARD_LD_COLUMNS` has 14 columns including `hg19_pos_1`, `hg38_pos_1`, `Dprime`, `+/-corr`. The new canonical schema has 6 columns: `chr`, `pos_1`, `pos_2`, `R2`, `rsID_1`, `rsID_2`. `Dprime` was always NaN. `+/-corr` is irrelevant to LD score regression. The build-specific position columns collapse to the single canonical `pos_1`/`pos_2` pair whose build is declared in schema metadata.

`build_standard_ld_table` currently takes `pair_rows` (list of `{i, j, R2, sign}` dicts) and `annotation_table`. It needs a `genome_build` parameter to select the right position column (`hg19_pos` vs `hg38_pos`).

`write_standard_ld_parquet` needs: `genome_build` (required), `row_group_size=50_000`, sort assertion on each incoming pair, and parquet schema metadata.

- [ ] **Step 1.1: Write the failing test for 6-column schema**

In `tests/test_ref_panel_builder.py`, replace the body of `test_build_standard_ld_table_uses_exact_schema`:

```python
def test_build_standard_ld_table_uses_exact_schema(self):
    annotation_table = pd.DataFrame(
        {
            "chr": ["1", "1"],
            "hg19_pos": [100, 120],
            "hg38_pos": [110, 130],
            "hg19_Uniq_ID": ["1:100:A:G", "1:120:C:T"],
            "hg38_Uniq_ID": ["1:110:A:G", "1:130:C:T"],
            "rsID": ["rs1", "rs2"],
            "MAF": [0.2, 0.3],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
        }
    )
    pair_rows = [{"i": 0, "j": 1, "R2": 0.75, "sign": "-"}]

    table = kernel_builder.build_standard_ld_table(
        pair_rows=pair_rows,
        annotation_table=annotation_table,
        genome_build="hg19",
    )

    self.assertEqual(
        table.columns.tolist(),
        ["chr", "pos_1", "pos_2", "R2", "rsID_1", "rsID_2"],
    )
    self.assertEqual(table.loc[0, "chr"], "1")
    self.assertEqual(table.loc[0, "pos_1"], 100)
    self.assertEqual(table.loc[0, "pos_2"], 120)
    self.assertAlmostEqual(table.loc[0, "R2"], 0.75, places=4)
    self.assertEqual(table.loc[0, "rsID_1"], "rs1")
    self.assertEqual(table.loc[0, "rsID_2"], "rs2")
```

- [ ] **Step 1.2: Run to verify it fails**

```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
python -m pytest tests/test_ref_panel_builder.py::RefPanelBuilderTest::test_build_standard_ld_table_uses_exact_schema -v
```

Expected: FAIL — `build_standard_ld_table()` takes no `genome_build` argument.

- [ ] **Step 1.3: Write the failing test for sort assertion**

Add this test method to the existing test class in `tests/test_ref_panel_builder.py`:

```python
@unittest.skipUnless(_has_module("pyarrow"), "pyarrow dependency is not installed")
def test_write_standard_ld_parquet_asserts_sort_invariant(self):
    import tempfile, pathlib
    annotation_table = pd.DataFrame(
        {
            "chr": ["1", "1", "1"],
            "hg19_pos": [100, 80, 120],   # 80 comes after 100 → sort violation
            "hg38_pos": [110, 90, 130],
            "hg19_Uniq_ID": ["1:100:A:G", "1:80:C:T", "1:120:G:A"],
            "hg38_Uniq_ID": ["1:110:A:G", "1:90:C:T", "1:130:G:A"],
            "rsID": ["rs1", "rs2", "rs3"],
            "MAF": [0.2, 0.3, 0.1],
            "REF": ["A", "C", "G"],
            "ALT": ["G", "T", "A"],
        }
    )
    # pair (0→1) has pos_1=100, then pair (1→2) has pos_1=80 → violates sort
    pair_rows = [{"i": 0, "j": 2, "R2": 0.5, "sign": "+"}, {"i": 1, "j": 2, "R2": 0.3, "sign": "+"}]
    with tempfile.TemporaryDirectory() as tmpdir:
        path = pathlib.Path(tmpdir) / "chr1.parquet"
        with self.assertRaises(ValueError, msg="Expected sort invariant violation error"):
            kernel_builder.write_standard_ld_parquet(
                pair_rows=iter(pair_rows),
                annotation_table=annotation_table,
                path=path,
                genome_build="hg19",
            )
```

- [ ] **Step 1.4: Write the failing test for schema metadata**

Add this test method:

```python
@unittest.skipUnless(_has_module("pyarrow"), "pyarrow dependency is not installed")
def test_write_standard_ld_parquet_writes_schema_metadata(self):
    import tempfile, pathlib
    import pyarrow.parquet as pq
    annotation_table = pd.DataFrame(
        {
            "chr": ["1", "1"],
            "hg19_pos": [100, 120],
            "hg38_pos": [110, 130],
            "hg19_Uniq_ID": ["1:100:A:G", "1:120:C:T"],
            "hg38_Uniq_ID": ["1:110:A:G", "1:130:C:T"],
            "rsID": ["rs1", "rs2"],
            "MAF": [0.2, 0.3],
            "REF": ["A", "C"],
            "ALT": ["G", "T"],
        }
    )
    pair_rows = [{"i": 0, "j": 1, "R2": 0.75, "sign": "+"}]
    with tempfile.TemporaryDirectory() as tmpdir:
        path = pathlib.Path(tmpdir) / "chr1.parquet"
        kernel_builder.write_standard_ld_parquet(
            pair_rows=iter(pair_rows),
            annotation_table=annotation_table,
            path=path,
            genome_build="hg19",
            row_group_size=50_000,
        )
        pf = pq.ParquetFile(str(path))
        meta = pf.schema_arrow.metadata
        self.assertIn(b"ldsc:sorted_by_build", meta)
        self.assertEqual(meta[b"ldsc:sorted_by_build"].decode("utf-8"), "hg19")
        self.assertIn(b"ldsc:row_group_size", meta)
        self.assertEqual(meta[b"ldsc:row_group_size"].decode("utf-8"), "50000")
```

- [ ] **Step 1.5: Run to verify these tests fail**

```bash
python -m pytest tests/test_ref_panel_builder.py -k "standard_ld" -v
```

Expected: all three new/updated tests FAIL.

- [ ] **Step 1.6: Implement the 6-column schema**

In `src/ldsc/_kernel/ref_panel_builder.py`, replace `_STANDARD_LD_COLUMNS` (lines 32-47):

```python
_STANDARD_LD_COLUMNS = [
    "chr",
    "pos_1",
    "pos_2",
    "R2",
    "rsID_1",
    "rsID_2",
]
```

Replace `build_standard_ld_table` (lines 453-487) with:

```python
def build_standard_ld_table(
    *,
    pair_rows: list[dict[str, float | int | str]],
    annotation_table: pd.DataFrame,
    genome_build: str,
) -> pd.DataFrame:
    """Build the canonical 6-column LD parquet table for one chromosome."""
    if not pair_rows:
        return pd.DataFrame(columns=_STANDARD_LD_COLUMNS)

    i = np.asarray([int(row["i"]) for row in pair_rows], dtype=np.int64)
    j = np.asarray([int(row["j"]) for row in pair_rows], dtype=np.int64)
    r2 = np.asarray([float(row["R2"]) for row in pair_rows], dtype=np.float32)
    left = annotation_table.iloc[i].reset_index(drop=True)
    right = annotation_table.iloc[j].reset_index(drop=True)
    pos_col = f"{genome_build}_pos"
    return pd.DataFrame(
        {
            "chr": left["chr"].astype(str),
            "pos_1": left[pos_col].to_numpy(dtype=np.int64),
            "pos_2": right[pos_col].to_numpy(dtype=np.int64),
            "R2": r2,
            "rsID_1": left["rsID"].astype(str),
            "rsID_2": right["rsID"].astype(str),
        },
        columns=_STANDARD_LD_COLUMNS,
    )
```

- [ ] **Step 1.7: Implement `write_standard_ld_parquet` with sort assertion and metadata**

Replace `write_standard_ld_parquet` (lines 528-572) with:

```python
def write_standard_ld_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    annotation_table: pd.DataFrame,
    path: str | PathLike[str],
    genome_build: str,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
) -> str:
    """Write the canonical LD parquet table, streaming batches when pyarrow is available."""
    _ensure_parent_dir(path)
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise ImportError(
            "Writing reference-panel LD parquet artifacts requires pyarrow."
        ) from exc

    pos_col = f"{genome_build}_pos"
    pa_meta = {
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
    }
    writer = None
    batch: list[dict[str, float | int | str]] = []
    prev_pos1 = -1
    try:
        for row in pair_rows:
            current_pos1 = int(annotation_table.iloc[int(row["i"])][pos_col])
            if current_pos1 < prev_pos1:
                raise ValueError(
                    f"Pairs must arrive in non-decreasing pos_1 order. "
                    f"Received pos_1={current_pos1} after pos_1={prev_pos1}. "
                    "Verify that the reference panel builder traverses SNPs in ascending positional order."
                )
            prev_pos1 = current_pos1
            batch.append(row)
            if len(batch) < batch_size:
                continue
            frame = build_standard_ld_table(pair_rows=batch, annotation_table=annotation_table, genome_build=genome_build)
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                enriched_schema = table.schema.with_metadata({**(table.schema.metadata or {}), **pa_meta})
                writer = pq.ParquetWriter(str(path), enriched_schema)
            writer.write_table(table, row_group_size=row_group_size)
            batch = []

        frame = build_standard_ld_table(pair_rows=batch, annotation_table=annotation_table, genome_build=genome_build)
        table = pa.Table.from_pandas(frame, preserve_index=False)
        if writer is None:
            enriched_schema = table.schema.with_metadata({**(table.schema.metadata or {}), **pa_meta})
            writer = pq.ParquetWriter(str(path), enriched_schema)
        writer.write_table(table, row_group_size=row_group_size)
    finally:
        if writer is not None:
            writer.close()
    return str(path)
```

- [ ] **Step 1.8: Run all three tests to verify they pass**

```bash
python -m pytest tests/test_ref_panel_builder.py -k "standard_ld" -v
```

Expected: all three PASS.

- [ ] **Step 1.9: Run the full ref_panel_builder test suite**

```bash
python -m pytest tests/test_ref_panel_builder.py -v
```

Expected: all tests PASS. If any test still calls `build_standard_ld_table` without `genome_build`, update that call to pass `genome_build="hg19"`.

- [ ] **Step 1.10: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "refactor: collapse parquet writer to 6-column canonical schema with sort assertion and metadata"
```

---

## Task 2: Refactor `SortedR2BlockReader.__init__` — canonical path

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py:1262-1345` (`SortedR2BlockReader.__init__`)
- Test: `tests/test_ldscore_legacy.py` (add `CanonicalParquetRuntimeTest` class)

### Context

Currently `__init__` opens all files as `pyarrow.Dataset` regardless of schema. After this task, the canonical path (schema has `pos_1`/`pos_2`) opens the file with `pq.ParquetFile`, validates the build (3-tier), warns if row groups are coarse, and builds `self._rg_bounds` — the list of `(min_pos1, max_pos1, rg_index)` tuples that powers the fast window query.

The existing `_parquet_schema_layout()` function already returns `"normalized"` for files with `pos_1`/`pos_2` (these are our canonical files). We rename this case internally to `"canonical"` in `self._runtime_layout` so downstream code branches cleanly.

**New attributes for canonical path:**
- `self._pf: pq.ParquetFile` — the open file handle
- `self._rg_bounds: list[tuple[int, int, int]]` — `(min_pos1, max_pos1, rg_index)` per row group
- `self._runtime_layout = "canonical"`

**Existing attributes preserved for both paths:**
- `self.chrom`, `self.identifier_mode`, `self.r2_bias_mode`, `self.r2_sample_size`, `self.genome_build`
- `self.pos`, `self.m`, `self.index_map`
- `self._last_query_key`, `self._last_query_rows`

**Legacy raw path** keeps `self.dataset` (pyarrow Dataset) and emits a `WARNING`. The existing `self.ds` (dataset module) is only needed for the raw path.

**3-tier build validation** (for canonical path only):
- Tier 1: `ldsc:sorted_by_build` metadata present and matches `self.genome_build` → proceed silently.
- Tier 2: metadata present but mismatches → raise `ValueError`.
- Tier 3: metadata absent → call `resolve_chr_pos_table` on the first row group's `(chr, pos_1)` pairs to infer the build; if inferred build matches, log WARNING; if not, raise `ValueError`.

**Coarse row group warning:** after build validation, if `metadata.num_rows / metadata.num_row_groups > 500_000`, emit a `WARNING` with the message from §3.2 of the design doc.

**Single-file requirement for canonical path:** canonical mode requires exactly one file per chromosome. If `len(paths) != 1`, raise `ValueError` immediately — no silent fallback. The error message should be:
`"canonical parquet_r2 backend requires exactly one file per chromosome; got {len(paths)} paths for chromosome {chrom}"`

- [ ] **Step 2.1: Write the failing test for canonical __init__**

Add a new test class to `tests/test_ldscore_legacy.py` below the existing `RawParquetRuntimeTest`:

```python
@unittest.skipIf(ld is None, "ldscore kernel is not available")
@unittest.skipUnless(_has_module("pyarrow"), "pyarrow dependency is not installed")
class CanonicalParquetRuntimeTest(unittest.TestCase):
    """Tests for the canonical 6-column parquet schema path."""

    def _write_canonical_parquet(self, path, *, genome_build="hg19", n_pairs=5, row_group_size=3):
        """Write a tiny canonical parquet for testing."""
        import pyarrow as pa
        import pyarrow.parquet as pq

        pairs = [
            {"chr": "1", "pos_1": 100 + i * 10, "pos_2": 120 + i * 10, "R2": 0.5, "rsID_1": f"rs{i+1}", "rsID_2": f"rs{i+2}"}
            for i in range(n_pairs)
        ]
        df = pd.DataFrame(pairs)
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {
            b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
            b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
        }
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched), row_group_size=row_group_size)

    def test_canonical_init_builds_rg_bounds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "chr1.parquet"
            self._write_canonical_parquet(path, n_pairs=6, row_group_size=2)
            metadata = pd.DataFrame(
                {"CHR": ["1"] * 7, "SNP": [f"rs{i+1}" for i in range(7)], "BP": [100 + i * 10 for i in range(7)], "CM": [0.0] * 7}
            )
            reader = ld.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )
            self.assertEqual(reader._runtime_layout, "canonical")
            self.assertIsNotNone(reader._rg_bounds)
            self.assertGreater(len(reader._rg_bounds), 1, "expected multiple row groups")
            for mn, mx, idx in reader._rg_bounds:
                self.assertLessEqual(mn, mx)
                self.assertIsInstance(idx, int)
```

- [ ] **Step 2.2: Run to verify it fails**

```bash
python -m pytest tests/test_ldscore_legacy.py::CanonicalParquetRuntimeTest::test_canonical_init_builds_rg_bounds -v
```

Expected: FAIL — `reader._runtime_layout` is `"normalized"`, not `"canonical"`; `_rg_bounds` attribute does not exist.

- [ ] **Step 2.3: Write the failing test for Tier 2 build validation (mismatch)**

Add to `CanonicalParquetRuntimeTest`:

```python
def test_canonical_init_raises_on_build_mismatch(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "chr1.parquet"
        self._write_canonical_parquet(path, genome_build="hg19")
        metadata = pd.DataFrame(
            {"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]}
        )
        with self.assertRaises(ValueError) as ctx:
            ld.SortedR2BlockReader(
                paths=[str(path)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg38",   # mismatch: file is hg19, analysis wants hg38
            )
        self.assertIn("hg19", str(ctx.exception))
        self.assertIn("hg38", str(ctx.exception))
```

- [ ] **Step 2.4: Write the failing test for coarse row group warning and multi-file hard-fail**

Add to `CanonicalParquetRuntimeTest`. The coarse warning uses `assertLogs` (consistent with the codebase's `LOGGER.warning` convention). The multi-file test verifies the hard-fail:

```python
def test_canonical_init_warns_on_coarse_row_groups(self):
    import pyarrow as pa
    import pyarrow.parquet as pq
    import unittest.mock

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "chr1.parquet"
        df = pd.DataFrame({"chr": ["1"], "pos_1": [100], "pos_2": [120], "R2": [0.5], "rsID_1": ["rs1"], "rsID_2": ["rs2"]})
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"50000"}
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched))

        # Patch metadata so num_rows / num_row_groups > 500_000
        real_pf_cls = pq.ParquetFile
        class FakeMeta:
            num_rows = 600_000
            num_row_groups = 1
            def row_group(self, i):
                return real_pf_cls(str(path)).metadata.row_group(0)
        class FakePF:
            def __init__(self, p):
                self._real = real_pf_cls(p)
                self.metadata = FakeMeta()
                self.schema_arrow = self._real.schema_arrow
            def read_row_groups(self, *a, **kw):
                return self._real.read_row_groups(*a, **kw)
            def read_row_group(self, *a, **kw):
                return self._real.read_row_group(*a, **kw)

        with unittest.mock.patch("pyarrow.parquet.ParquetFile", FakePF):
            metadata = pd.DataFrame({"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]})
            with self.assertLogs("LDSC.ldscore", level="WARNING") as log_ctx:
                ld.SortedR2BlockReader(
                    paths=[str(path)],
                    chrom="1",
                    metadata=metadata,
                    identifier_mode="rsID",
                    r2_bias_mode="unbiased",
                    r2_sample_size=None,
                    genome_build="hg19",
                )
            self.assertTrue(any("row group" in msg.lower() for msg in log_ctx.output))

def test_canonical_init_hard_fails_on_multiple_paths(self):
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as tmpdir:
        path1 = Path(tmpdir) / "chr1_a.parquet"
        path2 = Path(tmpdir) / "chr1_b.parquet"
        for path in (path1, path2):
            df = pd.DataFrame({"chr": ["1"], "pos_1": [100], "pos_2": [120], "R2": [0.5], "rsID_1": ["rs1"], "rsID_2": ["rs2"]})
            table = pa.Table.from_pandas(df, preserve_index=False)
            meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"50000"}
            enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
            with pq.ParquetWriter(str(path), enriched) as writer:
                writer.write_table(table.cast(enriched))

        metadata = pd.DataFrame({"CHR": ["1", "1"], "SNP": ["rs1", "rs2"], "BP": [100, 120], "CM": [0.0, 0.0]})
        with self.assertRaises(ValueError) as ctx:
            ld.SortedR2BlockReader(
                paths=[str(path1), str(path2)],
                chrom="1",
                metadata=metadata,
                identifier_mode="rsID",
                r2_bias_mode="unbiased",
                r2_sample_size=None,
                genome_build="hg19",
            )
        self.assertIn("exactly one file", str(ctx.exception))
```

Note: the logger name `"LDSC.ldscore"` must match the module-level `LOGGER = logging.getLogger("LDSC.ldscore")` in `ldscore.py`. Verify this before running the test.

- [ ] **Step 2.5: Run to verify new tests fail**

```bash
python -m pytest tests/test_ldscore_legacy.py::CanonicalParquetRuntimeTest -v
```

Expected: all three new tests FAIL.

- [ ] **Step 2.6: Implement the canonical `__init__` path**

In `src/ldsc/_kernel/ldscore.py`, replace the `__init__` body of `SortedR2BlockReader` (lines 1272-1345). The key structural change: detect layout first, then branch. Keep the entire existing metadata/index_map logic unchanged — it runs for both paths.

```python
def __init__(
    self,
    paths: Sequence[str],
    chrom: str,
    metadata: pd.DataFrame,
    identifier_mode: str,
    r2_bias_mode: str,
    r2_sample_size: float | None,
    genome_build: str | None = None,
) -> None:
    """Open one chromosome's sorted parquet R2 tables and build index maps."""
    if not paths:
        raise FileNotFoundError(f"No sorted parquet R2 files resolved for chromosome {chrom}.")

    self.chrom = normalize_chromosome(chrom)
    self.identifier_mode = identifier_mode
    self.r2_bias_mode = r2_bias_mode
    self.r2_sample_size = r2_sample_size
    self.genome_build = genome_build

    # --- Normalize metadata sidecar columns ---
    metadata = metadata.copy()
    metadata_context = f"SortedR2BlockReader[{self.chrom}] metadata"
    renamed = {
        resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=metadata_context): "CHR",
        resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=metadata_context): "POS",
        resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=metadata_context): "SNP",
    }
    optional_cm = resolve_optional_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=metadata_context)
    if optional_cm is not None:
        renamed[optional_cm] = "CM"
    metadata = metadata.rename(columns=renamed)
    if self.genome_build == "auto":
        validate_auto_genome_build_mode(self.identifier_mode, self.genome_build)
        metadata, inference = resolve_chr_pos_table(
            metadata,
            context=f"SortedR2BlockReader[{self.chrom}]",
            reference_table=load_packaged_reference_table(),
            logger=LOGGER,
        )
        self.genome_build = inference.genome_build
    validate_retained_identifier_uniqueness(metadata, identifier_mode, chrom)
    self.pos = metadata["POS"].to_numpy(dtype=np.int64)
    self.m = len(metadata)
    if self.identifier_mode == "rsid":
        self.index_map = {str(snp): idx for idx, snp in enumerate(metadata["SNP"].astype(str))}
    else:
        self.index_map = {int(pos): idx for idx, pos in enumerate(metadata["POS"].astype(np.int64))}

    self._last_query_key: tuple[int, int] | None = None
    self._last_query_rows: pd.DataFrame | None = None

    # --- Detect schema layout and open file handle ---
    # Peek at the first file's schema to classify as canonical or raw.
    try:
        import pyarrow.parquet as pq as _pq_mod
    except ImportError:
        _pq_mod = None
    try:
        import pyarrow.dataset as _ds_mod
    except ImportError:
        _ds_mod = None

    probe_schema_names: list[str] = []
    if _pq_mod is not None:
        probe_schema_names = list(_pq_mod.ParquetFile(paths[0]).schema_arrow.names)
    elif _ds_mod is not None:
        probe_schema_names = list(_ds_mod.dataset(list(paths), format="parquet").schema.names)

    layout = _parquet_schema_layout(probe_schema_names)

    if layout == "normalized":
        # Canonical path requires exactly one file — hard-fail if multiple paths given
        if len(paths) != 1:
            raise ValueError(
                f"canonical parquet_r2 backend requires exactly one file per chromosome; "
                f"got {len(paths)} paths for chromosome {chrom}"
            )
        if _pq_mod is None:
            raise ImportError("pyarrow is required for canonical parquet R2 input.")
        self._runtime_layout = "canonical"
        self._pf = _pq_mod.ParquetFile(paths[0])
        self.dataset = None
        self.ds = None
        self._raw_pos_columns = None
        self._raw_query_columns = None
        self._init_canonical_path(paths[0], _pq_mod)
    elif layout == "raw":
        # Legacy raw schema: fall back to Dataset full-scan with deprecation warning
        self._runtime_layout = "raw"
        self._pf = None
        self._rg_bounds = []
        if _ds_mod is None:
            raise ImportError("pyarrow is required for parquet R2 input.")
        self.ds = _ds_mod
        self.dataset = _ds_mod.dataset(list(paths), format="parquet")
        self.genome_build = _require_runtime_genome_build(self.genome_build)
        raw_mapping = _resolve_r2_source_columns(self.dataset.schema.names, context="raw parquet R2 schema")
        self._raw_pos_columns = get_r2_build_columns(self.genome_build, self.dataset.schema.names)
        self._raw_query_columns = [raw_mapping[canonical] for canonical in R2_CANONICAL_SOURCE_COLUMNS]
        LOGGER.warning(
            "%s uses the legacy raw schema. Row-group pruning is disabled and "
            "query performance will be severely degraded. Convert to the canonical "
            "format with `ldsc normalize-r2-parquet` to restore full performance.",
            paths[0],
        )
    else:
        raise ValueError(
            "Parquet R2 input must contain either normalized canonical columns "
            "`chr`, `pos_1`, `pos_2` or the raw legacy pairwise columns."
        )
```

Note: the `import pyarrow.parquet as pq as _pq_mod` syntax is invalid Python. Write it as two separate lines:

```python
    try:
        import pyarrow.parquet as _pq_mod
    except ImportError:
        _pq_mod = None
    try:
        import pyarrow.dataset as _ds_mod
    except ImportError:
        _ds_mod = None
```

- [ ] **Step 2.7: Implement `_init_canonical_path` helper**

Add this method to `SortedR2BlockReader` (place it after `__init__`):

```python
def _init_canonical_path(self, path: str, pq_mod) -> None:
    """Validate build, warn if coarse, and build the row-group bounds index."""
    pf = self._pf
    pf_meta = pf.metadata
    schema_arrow = pf.schema_arrow
    raw_schema_meta = schema_arrow.metadata or {}

    # Validate required 6 columns are present
    required = {"chr", "pos_1", "pos_2", "R2", "rsID_1", "rsID_2"}
    missing = required - set(schema_arrow.names)
    if missing:
        raise ValueError(
            f"Canonical parquet {path} is missing required columns: {sorted(missing)}"
        )

    # 3-tier build validation
    parquet_build_raw = raw_schema_meta.get(b"ldsc:sorted_by_build")
    if parquet_build_raw is not None:
        parquet_build = parquet_build_raw.decode("utf-8")
        if self.genome_build is not None and self.genome_build not in {"auto", parquet_build}:
            raise ValueError(
                f"Parquet sorted for {parquet_build} but analysis uses {self.genome_build}. "
                f"Use the correct reference file or regenerate with "
                f"`--genome-build {self.genome_build}`."
            )
        self.genome_build = parquet_build   # adopt file's declared build
    else:
        # Tier 3: infer from first row group
        first_tbl = pf.read_row_group(0, columns=["chr", "pos_1"])
        infer_df = pd.DataFrame({
            "CHR": first_tbl.column("chr").to_numpy(zero_copy_only=False),
            "POS": first_tbl.column("pos_1").to_numpy(zero_copy_only=False).astype(np.int64),
        })
        _, inference = resolve_chr_pos_table(
            infer_df,
            context=path,
            reference_table=load_packaged_reference_table(),
            logger=LOGGER,
        )
        inferred_build = inference.genome_build
        LOGGER.warning(
            "No build metadata found in %s; inferred %s from first row group. "
            "To silence this warning, regenerate the parquet with `ldsc build-ref-panel`.",
            path,
            inferred_build,
        )
        if self.genome_build is not None and self.genome_build not in {"auto", inferred_build}:
            raise ValueError(
                f"Parquet inferred as {inferred_build} but analysis uses {self.genome_build}."
            )
        self.genome_build = inferred_build

    # Warn if row groups are coarse
    if pf_meta.num_row_groups > 0:
        avg_rows_per_rg = pf_meta.num_rows / pf_meta.num_row_groups
        if avg_rows_per_rg > 500_000:
            LOGGER.warning(
                "%s has %d row group(s) (avg %,.0f rows/group). Query performance will be "
                "severely degraded. Regenerate with `row_group_size=50000` for optimal speed.",
                path,
                pf_meta.num_row_groups,
                avg_rows_per_rg,
            )

    # Build row-group bounds index from footer statistics
    col_idx = schema_arrow.get_field_index("pos_1")
    self._rg_bounds: list[tuple[int, int, int]] = []
    for i in range(pf_meta.num_row_groups):
        rg = pf_meta.row_group(i)
        stats = rg.column(col_idx).statistics
        if stats is not None and stats.has_min_max:
            self._rg_bounds.append((int(stats.min), int(stats.max), i))
```

- [ ] **Step 2.8: Run canonical tests to verify they pass**

```bash
python -m pytest tests/test_ldscore_legacy.py::CanonicalParquetRuntimeTest -v
```

Expected: all three PASS. If the coarse warning test fails because `warnings.warn` doesn't integrate with `assertWarns`, switch to `LOGGER.warning(...)` and use `assertLogs` in the test.

- [ ] **Step 2.9: Verify legacy test still passes**

```bash
python -m pytest tests/test_ldscore_legacy.py::RawParquetRuntimeTest -v
```

Expected: PASS (legacy Dataset path is unchanged).

- [ ] **Step 2.10: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_legacy.py
git commit -m "feat: SortedR2BlockReader canonical path — ParquetFile open, 3-tier build validation, rg_bounds index"
```

---

## Task 3: Canonical `_query_union_rows` with row-group pruning

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py:1359-1415` (`SortedR2BlockReader._query_union_rows`)
- Test: `tests/test_ldscore_legacy.py` (`CanonicalParquetRuntimeTest`)

### Context

The current `_query_union_rows` always calls `dataset.to_table(filter=...)` — a full-chromosome scan for every window. The new canonical path:
1. Looks up overlapping row groups from `self._rg_bounds`: `[i for mn, mx, i in self._rg_bounds if mn <= pos_max and mx >= pos_min]`
2. Reads only those groups: `pf.read_row_groups(rg_idxs, columns=[...])`
3. Applies fine-grained positional mask via numpy (row groups overlap the window but may contain rows outside it)
4. Maps identifiers to matrix indices via `self.index_map`
5. Returns a `pd.DataFrame` with columns `i`, `j`, `R2` (same contract as before — callers are unchanged)

The legacy raw path (`self._runtime_layout == "raw"`) is unchanged.

If `self._runtime_layout == "canonical"` but `self._pf is None` (multi-file fallback from Task 2), fall back to the existing Dataset query using `self.dataset`.

- [ ] **Step 3.1: Write the failing test for canonical window query correctness**

Add to `CanonicalParquetRuntimeTest` in `tests/test_ldscore_legacy.py`:

```python
def test_canonical_within_block_matrix_correctness(self):
    """Canonical path should produce the same within-block matrix as the raw path."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "chr1.parquet"
        # Three SNPs: rs1@100, rs2@120, rs3@140. Pairs: (rs1,rs2)=0.4, (rs1,rs3)=0.2, (rs2,rs3)=0.6.
        df = pd.DataFrame({
            "chr": ["1", "1", "1"],
            "pos_1": [100, 100, 120],
            "pos_2": [120, 140, 140],
            "R2": [0.4, 0.2, 0.6],
            "rsID_1": ["rs1", "rs1", "rs2"],
            "rsID_2": ["rs2", "rs3", "rs3"],
        })
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"2"}
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched), row_group_size=2)

        metadata = pd.DataFrame({
            "CHR": ["1", "1", "1"],
            "SNP": ["rs1", "rs2", "rs3"],
            "BP": [100, 120, 140],
            "CM": [0.0, 0.0, 0.0],
        })
        reader = ld.SortedR2BlockReader(
            paths=[str(path)],
            chrom="1",
            metadata=metadata,
            identifier_mode="rsID",
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            genome_build="hg19",
        )
        self.assertEqual(reader._runtime_layout, "canonical")
        matrix = reader.within_block_matrix(l_B=0, c=3)
        expected = np.array(
            [[1.0, 0.4, 0.2],
             [0.4, 1.0, 0.6],
             [0.2, 0.6, 1.0]],
            dtype=np.float32,
        )
        assert_array_almost_equal(matrix, expected)
```

- [ ] **Step 3.2: Write the failing test for row-group pruning (only correct groups are read)**

Add to `CanonicalParquetRuntimeTest`:

```python
def test_canonical_query_reads_only_overlapping_row_groups(self):
    """_query_union_rows must touch only row groups that overlap the query window."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "chr1.parquet"
        # 6 pairs sorted by pos_1; row_group_size=2 → 3 row groups.
        # Window query [100, 120] should touch only rg0.
        df = pd.DataFrame({
            "chr": ["1"] * 6,
            "pos_1": [100, 100, 200, 200, 300, 300],
            "pos_2": [120, 130, 220, 230, 320, 330],
            "R2": [0.4, 0.3, 0.5, 0.6, 0.7, 0.8],
            "rsID_1": ["rs1", "rs1", "rs3", "rs3", "rs5", "rs5"],
            "rsID_2": ["rs2", "rs2b", "rs4", "rs4b", "rs6", "rs6b"],
        })
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"2"}
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched), row_group_size=2)

        metadata = pd.DataFrame({
            "CHR": ["1"] * 7,
            "SNP": [f"rs{i+1}" for i in range(7)],
            "BP": [100, 120, 130, 200, 220, 300, 320],
            "CM": [0.0] * 7,
        })
        reader = ld.SortedR2BlockReader(
            paths=[str(path)],
            chrom="1",
            metadata=metadata,
            identifier_mode="rsID",
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            genome_build="hg19",
        )
        # Row groups: rg0=[100,100], rg1=[200,200], rg2=[300,300]
        # Query [100, 130] should match only rg0
        rg_idxs_used = []
        original_read_row_groups = reader._pf.read_row_groups
        def mock_read_rgs(idxs, **kwargs):
            rg_idxs_used.extend(idxs)
            return original_read_row_groups(idxs, **kwargs)
        reader._pf.read_row_groups = mock_read_rgs

        rows = reader._query_union_rows(100, 130)
        self.assertEqual(sorted(set(rg_idxs_used)), [0], "only rg0 should have been read")
        self.assertEqual(len(rows), 2, "expected 2 pairs in window [100,130]")
```

- [ ] **Step 3.3: Run to verify tests fail**

```bash
python -m pytest tests/test_ldscore_legacy.py::CanonicalParquetRuntimeTest::test_canonical_within_block_matrix_correctness tests/test_ldscore_legacy.py::CanonicalParquetRuntimeTest::test_canonical_query_reads_only_overlapping_row_groups -v
```

Expected: FAIL — canonical `_query_union_rows` doesn't exist yet (falls through to Dataset path or errors).

- [ ] **Step 3.4: Implement the canonical `_query_union_rows`**

Replace `_query_union_rows` in `SortedR2BlockReader` (lines 1359-1415):

```python
def _query_union_rows(self, pos_min: int, pos_max: int) -> pd.DataFrame:
    """Query cached or on-disk pair rows spanning a union genomic window."""
    key = (int(pos_min), int(pos_max))
    if self._last_query_key == key and self._last_query_rows is not None:
        return self._last_query_rows.copy()

    empty_cols = {"i": pd.Series([], dtype=np.int64), "j": pd.Series([], dtype=np.int64), "R2": pd.Series([], dtype=np.float32)}

    if self._runtime_layout == "canonical" and self._pf is not None:
        rows = self._query_union_rows_canonical(pos_min, pos_max, empty_cols)
    elif self._runtime_layout == "canonical" and self._pf is None:
        # Multi-file fallback: use Dataset (no pruning)
        rows = self._query_union_rows_dataset_canonical(pos_min, pos_max, empty_cols)
    else:
        # Legacy raw schema path (unchanged logic)
        rows = self._query_union_rows_raw(pos_min, pos_max, empty_cols)

    self._last_query_key = key
    self._last_query_rows = rows.copy()
    return rows

def _query_union_rows_canonical(self, pos_min: int, pos_max: int, empty_cols: dict) -> pd.DataFrame:
    """Fast path: row-group index + read_row_groups + to_numpy."""
    rg_idxs = [i for mn, mx, i in self._rg_bounds if mn <= pos_max and mx >= pos_min]
    if not rg_idxs:
        return pd.DataFrame(empty_cols)

    if self.identifier_mode == "rsid":
        read_cols = ["pos_1", "pos_2", "R2", "rsID_1", "rsID_2"]
    else:
        read_cols = ["pos_1", "pos_2", "R2"]

    tbl = self._pf.read_row_groups(rg_idxs, columns=read_cols)
    pos_1 = tbl.column("pos_1").to_numpy(zero_copy_only=False).astype(np.int64)
    pos_2 = tbl.column("pos_2").to_numpy(zero_copy_only=False).astype(np.int64)
    r2_raw = tbl.column("R2").to_numpy(zero_copy_only=False).astype(np.float32)

    # Fine-grained mask: row groups overlap window but may have rows outside
    mask = (pos_1 >= int(pos_min)) & (pos_2 <= int(pos_max))
    pos_1 = pos_1[mask]
    pos_2 = pos_2[mask]
    r2_raw = r2_raw[mask]

    if len(pos_1) == 0:
        return pd.DataFrame(empty_cols)

    r2 = self._transform_r2(r2_raw)

    if self.identifier_mode == "rsid":
        id1 = tbl.column("rsID_1").to_numpy(zero_copy_only=False)[mask].astype(str)
        id2 = tbl.column("rsID_2").to_numpy(zero_copy_only=False)[mask].astype(str)
        i_raw = [self.index_map.get(k) for k in id1]
        j_raw = [self.index_map.get(k) for k in id2]
    else:
        i_raw = [self.index_map.get(int(k)) for k in pos_1]
        j_raw = [self.index_map.get(int(k)) for k in pos_2]

    rows = pd.DataFrame({
        "i": pd.array(i_raw, dtype=pd.Int64Dtype()),
        "j": pd.array(j_raw, dtype=pd.Int64Dtype()),
        "R2": r2,
    })
    rows = rows.dropna(subset=["i", "j"]).copy()
    rows["i"] = rows["i"].astype(np.int64)
    rows["j"] = rows["j"].astype(np.int64)
    return rows

def _query_union_rows_dataset_canonical(self, pos_min: int, pos_max: int, empty_cols: dict) -> pd.DataFrame:
    """Multi-file canonical fallback: Dataset full-scan (no row-group pruning)."""
    filter_expr = (
        (self.ds.field("chr") == self.chrom)
        & (self.ds.field("pos_1") >= int(pos_min))
        & (self.ds.field("pos_2") <= int(pos_max))
    )
    cols = ["pos_1", "pos_2", "R2"]
    if self.identifier_mode == "rsid":
        cols += ["rsID_1", "rsID_2"]
    table = self.dataset.to_table(columns=cols, filter=filter_expr)
    rows = table.to_pandas()
    if len(rows) == 0:
        return pd.DataFrame(empty_cols)
    rows["R2"] = self._transform_r2(rows["R2"].to_numpy(dtype=np.float32))
    if self.identifier_mode == "rsid":
        rows["i"] = rows["rsID_1"].astype(str).map(self.index_map)
        rows["j"] = rows["rsID_2"].astype(str).map(self.index_map)
    else:
        rows["i"] = rows["pos_1"].astype(np.int64).map(self.index_map)
        rows["j"] = rows["pos_2"].astype(np.int64).map(self.index_map)
    rows = rows.dropna(subset=["i", "j"]).copy()
    rows["i"] = rows["i"].astype(np.int64)
    rows["j"] = rows["j"].astype(np.int64)
    return rows

def _query_union_rows_raw(self, pos_min: int, pos_max: int, empty_cols: dict) -> pd.DataFrame:
    """Legacy raw schema path: Dataset full-scan with canonicalization."""
    left_pos_col, right_pos_col = self._raw_pos_columns
    filter_expr = (
        (self.ds.field("chr") == self.chrom)
        & (self.ds.field(left_pos_col) >= int(pos_min))
        & (self.ds.field(left_pos_col) <= int(pos_max))
        & (self.ds.field(right_pos_col) >= int(pos_min))
        & (self.ds.field(right_pos_col) <= int(pos_max))
    )
    table = self.dataset.to_table(columns=self._raw_query_columns, filter=filter_expr)
    rows = table.to_pandas()
    if len(rows) == 0:
        base_columns = ["chr", "pos_1", "pos_2", "R2"]
        if self.identifier_mode == "rsid":
            base_columns.extend(["rsID_1", "rsID_2"])
        return pd.DataFrame(columns=base_columns + ["i", "j"])
    rows = rows.loc[
        rows["chr"].map(lambda value: normalize_chromosome(value, context="raw parquet R2 query")) == self.chrom
    ].copy()
    keep = (
        pd.to_numeric(rows[left_pos_col], errors="raise").between(int(pos_min), int(pos_max))
        & pd.to_numeric(rows[right_pos_col], errors="raise").between(int(pos_min), int(pos_max))
    )
    rows = rows.loc[keep].reset_index(drop=True)
    if len(rows) > 0:
        rows = canonicalize_r2_pairs(rows, self.genome_build)
    rows["R2"] = self._transform_r2(pd.to_numeric(rows["R2"], errors="raise").to_numpy(dtype=np.float32))
    if self.identifier_mode == "rsid":
        rows["i"] = rows["rsID_1"].astype(str).map(self.index_map)
        rows["j"] = rows["rsID_2"].astype(str).map(self.index_map)
    else:
        rows["i"] = pd.to_numeric(rows["pos_1"], errors="raise").astype(np.int64).map(self.index_map)
        rows["j"] = pd.to_numeric(rows["pos_2"], errors="raise").astype(np.int64).map(self.index_map)
    rows = rows.dropna(subset=["i", "j"]).copy()
    rows["i"] = rows["i"].astype(np.int64)
    rows["j"] = rows["j"].astype(np.int64)
    return rows
```

Also remove `self.query_columns` and `self._runtime_layout = _parquet_schema_layout(...)` assignments from the old `__init__` (they are no longer needed as the layout is now stored directly as `"canonical"` or `"raw"`).

- [ ] **Step 3.5: Run the new canonical tests**

```bash
python -m pytest tests/test_ldscore_legacy.py::CanonicalParquetRuntimeTest -v
```

Expected: all PASS.

- [ ] **Step 3.6: Run the full legacy ldscore test suite**

```bash
python -m pytest tests/test_ldscore_legacy.py -v
```

Expected: all PASS. The raw schema test at `RawParquetRuntimeTest::test_raw_parquet_presence_and_reader_support_hg19` should still pass (legacy path is unchanged). If it fails because of the `read_sorted_r2_presence` call, note that we will fix that in Task 4.

- [ ] **Step 3.7: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_legacy.py
git commit -m "feat: SortedR2BlockReader canonical _query_union_rows — row-group pruning with ParquetFile.read_row_groups"
```

---

## Task 4: `get_present_identifiers()` + `compute_chrom_from_parquet` migration + remove `read_sorted_r2_presence`

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` — add `get_present_identifiers()`; update `compute_chrom_from_parquet`; remove `read_sorted_r2_presence`
- Modify: `tests/test_ldscore_legacy.py` — update `RawParquetRuntimeTest` to remove `read_sorted_r2_presence` call; add `get_present_identifiers` tests

### Context

`read_sorted_r2_presence` (lines 888-938) is a standalone function that opens the parquet via `pyarrow.Dataset`, reads the full chromosome, and returns a set of all SNP identifiers present. It causes a ~6 GB RAM spike for large files and forces the file to be opened twice per chromosome.

The replacement is `SortedR2BlockReader.get_present_identifiers()`: streams through row groups one at a time via `pf.read_row_group()`, accumulating identifiers into a Python set. Peak RAM is ~800 KB (one row group at a time) regardless of file size.

`compute_chrom_from_parquet` currently calls `read_sorted_r2_presence(args, chrom)` BEFORE constructing `SortedR2BlockReader`, opening the file twice. After this task, `SortedR2BlockReader` is constructed first, then `block_reader.get_present_identifiers()` is called — one file open per chromosome.

- [ ] **Step 4.1: Write the failing test for `get_present_identifiers` (canonical path)**

Add to `CanonicalParquetRuntimeTest` in `tests/test_ldscore_legacy.py`:

```python
def test_get_present_identifiers_canonical(self):
    """get_present_identifiers must return all unique rsIDs present in the parquet."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "chr1.parquet"
        df = pd.DataFrame({
            "chr": ["1", "1", "1"],
            "pos_1": [100, 100, 200],
            "pos_2": [120, 200, 300],
            "R2": [0.4, 0.2, 0.6],
            "rsID_1": ["rs1", "rs1", "rs2"],
            "rsID_2": ["rs2", "rs3", "rs3"],
        })
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"2"}
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched), row_group_size=2)

        metadata = pd.DataFrame({
            "CHR": ["1", "1", "1"],
            "SNP": ["rs1", "rs2", "rs3"],
            "BP": [100, 200, 300],
            "CM": [0.0, 0.0, 0.0],
        })
        reader = ld.SortedR2BlockReader(
            paths=[str(path)],
            chrom="1",
            metadata=metadata,
            identifier_mode="rsID",
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            genome_build="hg19",
        )
        present = reader.get_present_identifiers()
        self.assertEqual(present, {"rs1", "rs2", "rs3"})
```

- [ ] **Step 4.2: Write the failing test for `get_present_identifiers` (chr_pos mode)**

Add to `CanonicalParquetRuntimeTest`:

```python
def test_get_present_identifiers_chr_pos_mode(self):
    import pyarrow as pa
    import pyarrow.parquet as pq

    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "chr1.parquet"
        df = pd.DataFrame({
            "chr": ["1", "1"],
            "pos_1": [100, 200],
            "pos_2": [200, 300],
            "R2": [0.4, 0.6],
            "rsID_1": ["rs1", "rs2"],
            "rsID_2": ["rs2", "rs3"],
        })
        table = pa.Table.from_pandas(df, preserve_index=False)
        meta = {b"ldsc:sorted_by_build": b"hg19", b"ldsc:row_group_size": b"2"}
        enriched = table.schema.with_metadata({**(table.schema.metadata or {}), **meta})
        with pq.ParquetWriter(str(path), enriched) as writer:
            writer.write_table(table.cast(enriched), row_group_size=2)

        metadata = pd.DataFrame({
            "CHR": ["1", "1", "1"],
            "SNP": ["rs1", "rs2", "rs3"],
            "BP": [100, 200, 300],
            "CM": [0.0, 0.0, 0.0],
        })
        reader = ld.SortedR2BlockReader(
            paths=[str(path)],
            chrom="1",
            metadata=metadata,
            identifier_mode="chr_pos",
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            genome_build="hg19",
        )
        present = reader.get_present_identifiers()
        self.assertEqual(present, {100, 200, 300})
```

- [ ] **Step 4.3: Update `RawParquetRuntimeTest` to remove `read_sorted_r2_presence` call**

In `tests/test_ldscore_legacy.py`, update `test_raw_parquet_presence_and_reader_support_hg19` to remove the `read_sorted_r2_presence` assertion and replace it with a `get_present_identifiers()` call:

```python
def test_raw_parquet_presence_and_reader_support_hg19(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "raw_chr1.parquet"
        pd.DataFrame(
            {
                "chr": ["1"],
                "rsID_1": ["rs2"],
                "rsID_2": ["rs1"],
                "hg38_bp1": [120],
                "hg38_bp2": [100],
                "hg19_bp_1": [20],
                "hg19_bp_2": [10],
                "hg38_Uniq_ID_1": ["1:120"],
                "hg38_Uniq_ID_2": ["1:100"],
                "hg19_Uniq_ID_1": ["1:20"],
                "hg19_Uniq_ID_2": ["1:10"],
                "R2": [0.4],
                "Dprime": [0.5],
                "+/-corr": ["+"],
            }
        ).to_parquet(path, index=False)

        metadata = pd.DataFrame(
            {
                "CHR": ["1", "1"],
                "SNP": ["rs1", "rs2"],
                "BP": [10, 20],
                "CM": [0.1, 0.2],
            }
        )
        reader = ld.SortedR2BlockReader(
            paths=[str(path)],
            chrom="1",
            metadata=metadata,
            identifier_mode="rsID",
            r2_bias_mode="unbiased",
            r2_sample_size=None,
            genome_build="hg19",
        )
        # Legacy raw path still exposes get_present_identifiers
        self.assertEqual(reader.get_present_identifiers(), {"rs1", "rs2"})
        matrix = reader.within_block_matrix(l_B=0, c=2)
        expected = np.array([[1.0, 0.4], [0.4, 1.0]], dtype=np.float32)
        assert_array_almost_equal(matrix, expected)
```

- [ ] **Step 4.4: Run to verify tests fail**

```bash
python -m pytest tests/test_ldscore_legacy.py -k "present_identifiers or raw_parquet" -v
```

Expected: `get_present_identifiers` tests FAIL (method doesn't exist); raw parquet test FAILS (method doesn't exist).

- [ ] **Step 4.5: Implement `get_present_identifiers()`**

Add this method to `SortedR2BlockReader` (after `_query_union_rows_raw`):

```python
def get_present_identifiers(self) -> set:
    """Return the set of all SNP identifiers present in the parquet as left or right SNPs.

    Uses streaming row-group reads to keep peak RAM constant (one row group at a time)
    regardless of file size. For canonical files this replaces the old
    read_sorted_r2_presence standalone function.
    """
    if self.identifier_mode == "rsid":
        cols = ["rsID_1", "rsID_2"]
    else:
        cols = ["pos_1", "pos_2"]

    present: set = set()

    if self._runtime_layout == "canonical" and self._pf is not None:
        for i in range(self._pf.metadata.num_row_groups):
            tbl = self._pf.read_row_group(i, columns=cols)
            for col in cols:
                arr = tbl.column(col).to_numpy(zero_copy_only=False)
                if self.identifier_mode == "rsid":
                    present.update(arr.astype(str))
                else:
                    present.update(arr.astype(np.int64))
    else:
        # Legacy raw schema or multi-file canonical: full Dataset scan
        if self.identifier_mode == "rsid":
            read_cols = ["rsID_1", "rsID_2"]
        else:
            # For raw schema, positions are in build-specific columns
            if self._raw_pos_columns is not None:
                read_cols = list(self._raw_pos_columns)
            else:
                read_cols = ["pos_1", "pos_2"]
        tbl = self.dataset.to_table(columns=read_cols)
        for col in read_cols:
            arr = tbl.column(col).to_numpy(zero_copy_only=False)
            if self.identifier_mode == "rsid":
                present.update(arr.astype(str))
            else:
                present.update(arr.astype(np.int64))

    return present
```

- [ ] **Step 4.6: Run presence tests to verify they pass**

```bash
python -m pytest tests/test_ldscore_legacy.py -k "present_identifiers or raw_parquet" -v
```

Expected: all PASS.

- [ ] **Step 4.7: Update `compute_chrom_from_parquet` to open reader first**

In `src/ldsc/_kernel/ldscore.py`, replace lines 1615-1640 of `compute_chrom_from_parquet`:

```python
    # Old order: read_sorted_r2_presence → filter → SortedR2BlockReader
    # New order: SortedR2BlockReader (file open once) → get_present_identifiers → filter
    block_reader = SortedR2BlockReader(
        paths=resolve_parquet_files(args, chrom=chrom),
        chrom=chrom,
        metadata=metadata,
        identifier_mode=args.snp_identifier,
        r2_bias_mode=args.r2_bias_mode,
        r2_sample_size=args.r2_sample_size,
        genome_build=args.genome_build,
    )
    present_values = block_reader.get_present_identifiers()
    metadata, annotations, regression_keys = filter_reference_to_present_r2(
        metadata, annotations, regression_keys, present_values, args.snp_identifier, chrom
    )
```

Remove the old lines:
```python
    present_values = read_sorted_r2_presence(args, chrom=chrom)
    metadata, annotations, regression_keys = filter_reference_to_present_r2(
        metadata, annotations, regression_keys, present_values, args.snp_identifier, chrom
    )
```

and the existing `block_reader = SortedR2BlockReader(...)` block that comes later (now consolidated above).

- [ ] **Step 4.8: Remove `read_sorted_r2_presence` standalone function**

Delete lines 888-938 in `src/ldsc/_kernel/ldscore.py` (the `read_sorted_r2_presence` function body). Verify no other callers exist:

```bash
grep -rn "read_sorted_r2_presence" /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured/src/
```

Expected: zero matches. If any remain, update them to use `block_reader.get_present_identifiers()`.

- [ ] **Step 4.9: Run the full test suite**

```bash
python -m pytest tests/ -v 2>&1 | tail -30
```

Expected: all tests PASS. Key tests to watch:
- `test_ldscore_legacy.py::CanonicalParquetRuntimeTest` — all canonical tests
- `test_ldscore_legacy.py::RawParquetRuntimeTest` — legacy schema still works
- `test_ref_panel_builder.py` — writer tests

If `test_ldscore_workflow.py` contains integration tests that call `compute_chrom_from_parquet`, verify they pass.

- [ ] **Step 4.10: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_legacy.py
git commit -m "feat: add SortedR2BlockReader.get_present_identifiers, migrate compute_chrom_from_parquet to single file open, remove read_sorted_r2_presence"
```

---

## Self-Review Checklist

- [x] **Spec coverage — §2.1 Column Schema:** Task 1 implements 6-column schema in writer.
- [x] **Spec coverage — §2.2 Sort invariant:** Task 1 adds sort assertion in `write_standard_ld_parquet`.
- [x] **Spec coverage — §2.4 Parquet metadata:** Task 1 writes `ldsc:sorted_by_build` and `ldsc:row_group_size`.
- [x] **Spec coverage — §3.1 File open:** Task 2 opens canonical files with `pq.ParquetFile`.
- [x] **Spec coverage — §3.2 Build validation (3 tiers):** Task 2 implements all three tiers in `_init_canonical_path`.
- [x] **Spec coverage — §3.2 Coarse warning:** Task 2 emits `UserWarning` when `avg_rows_per_rg > 500_000`.
- [x] **Spec coverage — §3.3 Row-group index:** Task 2 builds `_rg_bounds` from footer statistics.
- [x] **Spec coverage — §3.4 `get_present_identifiers`:** Task 4 adds streaming scan method.
- [x] **Spec coverage — §3.5 `_query_union_rows`:** Task 3 implements row-group pruning + `to_numpy()`.
- [x] **Spec coverage — §4 Legacy shim:** Task 2 detects raw schema, opens as Dataset, emits deprecation warning.
- [x] **Spec coverage — §5 Affected modules table:** All listed modules are covered by tasks 1-4.
- [x] **`read_sorted_r2_presence` removal:** Task 4 removes it and the `RawParquetRuntimeTest` test is updated.
- [x] **Type consistency:** `_rg_bounds` defined as `list[tuple[int, int, int]]` in Task 2; used in Task 3 as `[i for mn, mx, i in self._rg_bounds if ...]`. Consistent.
- [x] **`build_standard_ld_table` signature:** `genome_build` parameter added in Task 1; used with keyword in all call sites.
- [x] **No placeholders:** all test code and implementation code is complete.
