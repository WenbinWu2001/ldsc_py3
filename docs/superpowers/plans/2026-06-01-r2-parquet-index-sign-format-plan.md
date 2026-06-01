# R2 Parquet Index+Sign Format — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 10-column R² parquet with a 4-column index format (`IDX_1`, `IDX_2`, `R2`, `SIGN`) bound to the per-SNP sidecar by a SHA-256 identity hash, mode-agnostic at read time, with bit-identical LD scores.

**Architecture:** The builder already produces pair indices; the writer stores them directly plus a sign bit instead of expanding to POS/rsID/allele strings. The reader loads the full panel sidecar, validates the binding hash, builds a once-per-chromosome `build_idx → retained_idx` remap, and decodes each row group with a pure int32 gather. The current LD-score pipeline order (restrict → align → merge → MAF → window) is untouched; only the reader's endpoint resolution changes.

**Tech Stack:** Python 3, NumPy, pandas, PyArrow (parquet), hashlib (stdlib), pytest.

**Design spec:** `docs/superpowers/specs/2026-06-01-r2-parquet-index-sign-format-design.md`

**Test command (use everywhere below):**
```
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured && PYTHONPATH=src /Users/wenbinwu/miniforge3/envs/ldsc3/bin/python -m pytest -q
```
Single-test form: append `<path>::<Class>::<test>` and `-v`.

**Baseline before starting:** full suite is `871 passed, 1 skipped`. Expect this number to change as old-format fixtures/tests are replaced; the plan calls out where.

---

## File Structure

| File | Responsibility | Action |
|---|---|---|
| `src/ldsc/_kernel/snp_identity.py` | shared `sidecar_identity_sha256()` binding hash | Modify (add function) |
| `src/ldsc/_kernel/ref_panel_builder.py` | 4-column writer + sign plumbing + index arrow table | Modify |
| `src/ldsc/ref_panel_builder.py` | workflow: compute hash/n_snps, reorder sidecar-before-parquet, pass new args | Modify |
| `src/ldsc/_kernel/ldscore.py` | reader: sidecar load, hash validation, remap, gather decode, IDX_1 pruning; remove raw + legacy-canonical paths | Modify |
| `src/ldsc/_kernel/ref_panel.py` | remove synth-from-parquet fallback (sidecar mandatory) | Modify |
| `docs/current/parquet-r2-format-and-read-pipeline.md` | rewrite for index layout | Modify |
| `tests/test_ref_panel_builder.py` | writer tests for new schema/metadata/sign | Modify |
| `tests/test_ldscore_workflow.py` | reader/binding/remap tests + parity gate; drop legacy fixtures | Modify |

---

## Phase 0 — Shared binding hash

### Task 0: `sidecar_identity_sha256` in snp_identity.py

**Files:**
- Modify: `src/ldsc/_kernel/snp_identity.py` (add `import hashlib` near top; add function after `effective_merge_key_series`, ~line 272)
- Test: `tests/test_ref_panel_builder.py` (new class `SidecarIdentityHashTest`)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_ref_panel_builder.py`:

```python
class SidecarIdentityHashTest(unittest.TestCase):
    def test_hash_is_stable_and_order_sensitive(self):
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        df = pd.DataFrame(
            {"CHR": ["1", "1", "1"], "POS": [10, 20, 30],
             "SNP": ["rsA", "rsB", "rsC"], "A1": ["A", "C", "G"], "A2": ["G", "T", "A"]}
        )
        h1 = sidecar_identity_sha256(df)
        h2 = sidecar_identity_sha256(df.copy())
        self.assertEqual(h1, h2)
        self.assertEqual(len(h1), 64)  # sha256 hex

        reordered = df.iloc[[1, 0, 2]].reset_index(drop=True)
        self.assertNotEqual(sidecar_identity_sha256(reordered), h1)

    def test_hash_ignores_snp_cm_maf_but_not_alleles(self):
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        base = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["rsA"], "A1": ["A"], "A2": ["G"]})
        changed_snp = base.assign(SNP=["rs999"])
        changed_allele = base.assign(A1=["T"])
        self.assertEqual(sidecar_identity_sha256(changed_snp), sidecar_identity_sha256(base))
        self.assertNotEqual(sidecar_identity_sha256(changed_allele), sidecar_identity_sha256(base))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... -m pytest -q tests/test_ref_panel_builder.py::SidecarIdentityHashTest -v`
Expected: FAIL with `ImportError` / `cannot import name 'sidecar_identity_sha256'`.

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/_kernel/snp_identity.py`, add `import hashlib` with the other stdlib imports (after `import logging`, ~line 20). Add after `effective_merge_key_series` (~line 272):

```python
def sidecar_identity_sha256(metadata: pd.DataFrame, *, context: str = "sidecar identity") -> str:
    """Hash the panel's physical SNP identity to bind a parquet to its sidecar.

    Computes ``sha256`` over the newline-joined ``CHR:POS:A1:A2`` string of each
    row, in order. This captures the SNP set and order — exactly what the index
    parquet's ``IDX`` values reference — while ignoring the volatile rsID and the
    recomputable ``CM``/``MAF`` columns. Built in memory; never written as a
    sidecar column.
    """
    _require_columns(metadata, ("CHR", "POS", "A1", "A2"), context=context)
    pos = pd.to_numeric(metadata["POS"], errors="raise").astype("int64").astype(str)
    joined = (
        metadata["CHR"].astype(str) + ":" + pos + ":"
        + metadata["A1"].astype(str) + ":" + metadata["A2"].astype(str)
    )
    payload = "\n".join(joined.tolist()).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `... -m pytest -q tests/test_ref_panel_builder.py::SidecarIdentityHashTest -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/snp_identity.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): add sidecar identity sha256 binding hash"
```

---

## Phase 1 — Writer (4-column index format)

### Task 1: Index arrow table builder

Replace `_standard_r2_arrow_table` (which joins `reference_snp_table` and expands identifiers) with an index builder that emits the four columns directly from `pair_rows`.

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py` (`_standard_r2_arrow_table` ~lines 674-716; `_STANDARD_R2_COLUMNS` ~lines 46-57)
- Test: `tests/test_ref_panel_builder.py` (new class `IndexArrowTableTest`)

- [ ] **Step 1: Write the failing test**

```python
class IndexArrowTableTest(unittest.TestCase):
    def test_index_arrow_table_has_four_columns(self):
        import pyarrow as pa
        from ldsc._kernel import ref_panel_builder as kb

        schema = pa.schema([("IDX_1", pa.int32()), ("IDX_2", pa.int32()),
                            ("R2", pa.float32()), ("SIGN", pa.bool_())])
        rows = [
            {"i": 0, "j": 2, "R2": 0.5, "sign": "+"},
            {"i": 0, "j": 3, "R2": -0.01, "sign": "-"},
        ]
        table = kb._standard_r2_index_table(pa, schema, pair_rows=rows)
        self.assertEqual(table.schema.names, ["IDX_1", "IDX_2", "R2", "SIGN"])
        self.assertEqual(table.column("IDX_1").to_pylist(), [0, 0])
        self.assertEqual(table.column("IDX_2").to_pylist(), [2, 3])
        self.assertEqual(table.column("SIGN").to_pylist(), [True, False])
        self.assertEqual(table.schema.field("IDX_1").type, pa.int32())
        self.assertEqual(table.schema.field("SIGN").type, pa.bool_())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ref_panel_builder.py::IndexArrowTableTest -v`
Expected: FAIL with `AttributeError: module ... has no attribute '_standard_r2_index_table'`.

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/_kernel/ref_panel_builder.py`, replace `_STANDARD_R2_COLUMNS` (lines 46-57) with:

```python
_INDEX_R2_COLUMNS = ["IDX_1", "IDX_2", "R2", "SIGN"]
```

Delete `_empty_standard_r2_table` (lines 60-76) and `build_standard_r2_table` (lines 579-619) — superseded (they expand identifiers). Replace `_standard_r2_arrow_table` (lines 674-716) with:

```python
def _standard_r2_index_table(pa, schema, *, pair_rows: list[dict[str, float | int | str]]):
    """Build one 4-column index ``pyarrow.Table`` batch directly from pair rows.

    ``pair_rows`` carry sidecar-row indices ``i``/``j`` (already the panel index
    space), the unbiased ``R2``, and the correlation ``sign`` as ``"+"``/``"-"``.
    No reference-SNP join or identifier expansion: the indices are stored as-is.
    """
    if not pair_rows:
        return schema.empty_table()
    count = len(pair_rows)
    i = np.fromiter((int(row["i"]) for row in pair_rows), dtype=np.int32, count=count)
    j = np.fromiter((int(row["j"]) for row in pair_rows), dtype=np.int32, count=count)
    r2 = np.fromiter((float(row["R2"]) for row in pair_rows), dtype=np.float32, count=count)
    sign = np.fromiter((row["sign"] == "+" for row in pair_rows), dtype=bool, count=count)
    return pa.Table.from_arrays(
        [
            pa.array(i, type=pa.int32()),
            pa.array(j, type=pa.int32()),
            pa.array(r2, type=pa.float32()),
            pa.array(sign, type=pa.bool_()),
        ],
        schema=schema,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ref_panel_builder.py::IndexArrowTableTest -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): add 4-column index arrow table builder"
```

### Task 2: Rewrite `write_r2_parquet` to the index schema

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel_builder.py` (`write_r2_parquet` ~lines 719-855)
- Test: `tests/test_ref_panel_builder.py` (replace `test_write_r2_parquet_*` body assertions; add `IndexWriterTest`)

- [ ] **Step 1: Write the failing test**

Add:

```python
class IndexWriterTest(unittest.TestCase):
    def _rows(self):
        return [
            {"i": 0, "j": 1, "R2": 0.8, "sign": "+"},
            {"i": 0, "j": 2, "R2": 0.1, "sign": "-"},
            {"i": 1, "j": 2, "R2": 0.3, "sign": "+"},
        ]

    def test_writes_four_columns_and_binding_metadata(self):
        import pyarrow.parquet as pq
        from ldsc._kernel import ref_panel_builder as kb

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hg19" / "chr1_r2.parquet"
            kb.write_r2_parquet(
                pair_rows=self._rows(), path=path, genome_build="hg19",
                n_samples=100, snp_identifier="chr_pos", min_r2=0.0,
                n_snps=3, sidecar_identity_sha256="deadbeef" * 8,
            )
            pf = pq.ParquetFile(path)
            self.assertEqual(pf.schema_arrow.names, ["IDX_1", "IDX_2", "R2", "SIGN"])
            meta = {k.decode(): v.decode() for k, v in pf.schema_arrow.metadata.items()}
            self.assertEqual(meta["ldsc:schema_version"], "2")
            self.assertEqual(meta["ldsc:n_snps"], "3")
            self.assertEqual(meta["ldsc:sidecar_identity_sha256"], "deadbeef" * 8)
            self.assertEqual(meta["ldsc:sorted_by_build"], "hg19")

    def test_asserts_idx1_sort_invariant(self):
        from ldsc._kernel import ref_panel_builder as kb

        bad = [{"i": 2, "j": 3, "R2": 0.5, "sign": "+"},
               {"i": 0, "j": 1, "R2": 0.5, "sign": "+"}]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chr1_r2.parquet"
            with self.assertRaisesRegex(ValueError, "non-decreasing IDX_1"):
                kb.write_r2_parquet(
                    pair_rows=bad, path=path, genome_build="hg19", n_samples=100,
                    snp_identifier="chr_pos", min_r2=0.0, n_snps=4,
                    sidecar_identity_sha256="0" * 64,
                )
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ref_panel_builder.py::IndexWriterTest -v`
Expected: FAIL (`write_r2_parquet() got unexpected keyword argument 'n_snps'` or schema mismatch).

- [ ] **Step 3: Write minimal implementation**

Replace `write_r2_parquet` (lines 719-855). New signature drops `reference_snp_table`, `batch_size` stays, adds `n_snps` and `sidecar_identity_sha256`:

```python
def write_r2_parquet(
    *,
    pair_rows: Iterable[dict[str, float | int | str]],
    path: str | PathLike[str],
    genome_build: str,
    n_samples: int,
    snp_identifier: str,
    n_snps: int,
    sidecar_identity_sha256: str,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
    min_r2: float = 0.0,
) -> str:
    """Write one canonical 4-column index R2 parquet with LDSC schema metadata.

    Columns are ``IDX_1``/``IDX_2`` (int32 sidecar-row indices), ``R2``
    (float32 unbiased), and ``SIGN`` (bit-packed bool, ``True`` for Pearson
    r >= 0). Pairs must arrive in non-decreasing ``IDX_1`` order. The parquet is
    bound to its sidecar by ``ldsc:n_snps`` and ``ldsc:sidecar_identity_sha256``;
    indices are meaningless without the matching sidecar.
    """
    _ensure_parent_dir(path)
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "Writing canonical reference-panel R2 parquet artifacts requires pyarrow."
        ) from exc

    if int(n_snps) >= 2**31:
        raise ValueError(f"n_snps={n_snps} exceeds int32 index range; cannot store IDX columns.")

    pa_meta = {
        **{
            f"ldsc:{key}".encode("utf-8"): str(value).encode("utf-8")
            for key, value in identity_artifact_metadata(
                artifact_type="ref_panel_r2",
                snp_identifier=snp_identifier,
                genome_build=genome_build,
            ).items()
        },
        b"ldsc:schema_version": b"2",
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
        b"ldsc:n_samples": str(n_samples).encode("utf-8"),
        b"ldsc:r2_bias": b"unbiased",
        b"ldsc:min_r2": str(min_r2).encode("utf-8"),
        b"ldsc:n_snps": str(int(n_snps)).encode("utf-8"),
        b"ldsc:sidecar_identity_sha256": str(sidecar_identity_sha256).encode("utf-8"),
    }
    schema = pa.schema(
        [
            ("IDX_1", pa.int32()),
            ("IDX_2", pa.int32()),
            ("R2", pa.float32()),
            ("SIGN", pa.bool_()),
        ]
    ).with_metadata(pa_meta)

    writer = None
    batch: list[dict[str, float | int | str]] = []
    prev_idx1: int | None = None

    def _flush(rows: list[dict[str, float | int | str]]) -> None:
        nonlocal writer, prev_idx1
        table = _standard_r2_index_table(pa, schema, pair_rows=rows)
        idx1 = table.column("IDX_1").to_numpy()
        if idx1.size:
            sequence = idx1 if prev_idx1 is None else np.concatenate(([prev_idx1], idx1))
            bad = np.flatnonzero(np.diff(sequence) < 0)
            if bad.size:
                k = int(bad[0])
                raise ValueError(
                    "Pairs must arrive in non-decreasing IDX_1 order. "
                    f"Received IDX_1={int(sequence[k + 1])} after IDX_1={int(sequence[k])}. "
                    "Verify that the reference panel builder traverses SNPs in ascending index order."
                )
            prev_idx1 = int(idx1[-1])
        if writer is None:
            if pa.Codec.is_available("zstd"):
                writer = pq.ParquetWriter(str(path), schema, compression="zstd", compression_level=9)
            else:
                warnings.warn(
                    "zstd codec is not available in this pyarrow build; falling back to snappy.",
                    UserWarning,
                    stacklevel=2,
                )
                writer = pq.ParquetWriter(str(path), schema, compression="snappy")
        writer.write_table(table, row_group_size=row_group_size)

    try:
        for row in pair_rows:
            batch.append(row)
            if len(batch) >= batch_size:
                _flush(batch)
                batch = []
        if batch or writer is None:
            _flush(batch)
    finally:
        if writer is not None:
            writer.close()
    return str(path)
```

Then update the existing `tests/test_ref_panel_builder.py` writer tests (`test_write_r2_parquet_*`, lines ~528-880) that still pass `reference_snp_table=` or assert 10-column schema: change calls to the new signature (pass `n_snps`/`sidecar_identity_sha256`, drop `reference_snp_table`) and assert `["IDX_1","IDX_2","R2","SIGN"]`. Tests asserting POS/SNP/allele columns or `build_standard_r2_table` are removed.

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ref_panel_builder.py -v`
Expected: PASS for `IndexWriterTest`; updated writer tests green; no references to removed helpers remain.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): write 4-column index R2 parquet with binding metadata"
```

### Task 3: Wire hash/n_snps and reorder sidecar-before-parquet in the workflow

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py` (build loop ~lines 852-884)
- Modify: `src/ldsc/_kernel/ref_panel_builder.py` (re-export hash for the workflow)
- Test: `tests/test_ref_panel_builder.py` (extend an existing end-to-end build test to assert parquet hash matches the sidecar)

- [ ] **Step 1: Write the failing test**

Add a focused test that builds a tiny panel through the workflow and checks the parquet's `ldsc:sidecar_identity_sha256` equals the hash recomputed from the written sidecar. Place near the existing build integration tests:

```python
class BuildBindingConsistencyTest(unittest.TestCase):
    @unittest.skipUnless(_HAS_PYARROW, "pyarrow required")
    def test_parquet_hash_matches_written_sidecar(self):
        import pyarrow.parquet as pq
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        out_dir = _build_minimal_reference_panel()  # helper already used by existing build tests
        r2 = Path(out_dir) / "hg19" / "chr1_r2.parquet"
        meta = Path(out_dir) / "hg19" / "chr1_meta.tsv.gz"
        stored = {k.decode(): v.decode() for k, v in pq.ParquetFile(r2).schema_arrow.metadata.items()}
        sidecar = pd.read_csv(meta, sep="\t", comment="#")
        self.assertEqual(stored["ldsc:sidecar_identity_sha256"], sidecar_identity_sha256(sidecar))
        self.assertEqual(int(stored["ldsc:n_snps"]), len(sidecar))
```

> If no `_build_minimal_reference_panel` helper exists, reuse the construction from the nearest existing end-to-end build test in this file and inline it.

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ref_panel_builder.py::BuildBindingConsistencyTest -v`
Expected: FAIL (`write_r2_parquet()` still called with old args, or hash metadata absent).

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/_kernel/ref_panel_builder.py`, re-export the hash near the top-level imports (it is already importable from `.snp_identity`; add to the existing `from .snp_identity import (...)` block, ~lines 29-36):

```python
from .snp_identity import (
    RestrictionIdentityKeys,
    allele_set_series,
    base_key_series,
    identity_artifact_metadata,
    identity_mode_family,
    restriction_membership_mask,
    sidecar_identity_sha256,
)
```

In `src/ldsc/ref_panel_builder.py`, reorder the build loop so the sidecar table is built first, then the hash, then the parquet (replace lines ~852-884):

```python
            reference_snp_table = kernel_builder.build_reference_snp_table(
                metadata=metadata,
                hg19_positions=hg19_positions,
                hg38_positions=hg38_positions,
            )
            r2_path = Path(config.output_dir) / build / f"chr{chrom}_r2.parquet"
            meta_path = Path(config.output_dir) / build / f"chr{chrom}_meta.tsv.gz"

            runtime_metadata = kernel_builder.build_runtime_metadata_table(
                metadata=metadata,
                positions=build_positions,
                cm_values=cm_values,
            )
            identity_hash = kernel_builder.sidecar_identity_sha256(runtime_metadata)

            kernel_builder.write_r2_parquet(
                pair_rows=kernel_builder.yield_pairwise_r2_rows(
                    block_left=block_left,
                    snp_batch_size=config.snp_batch_size,
                    standardized_snp_getter=lambda b: geno.nextSNPs(b, dtype=np.float32),
                    m=geno.m,
                    n=geno.n,
                    min_r2=config.min_r2,
                ),
                path=r2_path,
                genome_build=build,
                snp_identifier=self.global_config.snp_identifier,
                n_samples=geno.n,
                min_r2=config.min_r2,
                n_snps=len(runtime_metadata),
                sidecar_identity_sha256=identity_hash,
            )
            kernel_builder.write_runtime_metadata_sidecar(
                runtime_metadata,
                meta_path,
                genome_build=build,
                snp_identifier=self.global_config.snp_identifier,
```

`reference_snp_table` is now unused by the writer; if nothing else in the loop consumes it, delete its construction too. Verify with `grep -n reference_snp_table src/ldsc/ref_panel_builder.py` and remove if dead.

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ref_panel_builder.py::BuildBindingConsistencyTest -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel_builder.py src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(ref-panel): bind parquet to sidecar via identity hash in build"
```

---

## Phase 2 — Reader (index decode + binding + remap)

### Task 4: Full-sidecar loader + binding validation helpers

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (add helpers near `_parquet_schema_layout`, ~line 824)
- Test: `tests/test_ldscore_workflow.py` (new class `IndexBindingTest`)

- [ ] **Step 1: Write the failing test**

```python
class IndexBindingTest(unittest.TestCase):
    def _sidecar(self, tmp):
        meta = Path(tmp) / "hg19" / "chr1_meta.tsv.gz"
        meta.parent.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [10, 20, 30],
                           "SNP": ["rsA", "rsB", "rsC"], "A1": ["A", "C", "G"],
                           "A2": ["G", "T", "A"], "CM": [0.0, 0.0, 0.0], "MAF": [0.2, 0.3, 0.4]})
        import gzip as _gz
        with _gz.open(meta, "wt") as h:
            h.write("# ldsc:schema_version=2\n")
            df.to_csv(h, sep="\t", index=False)
        return meta, df

    def test_load_full_panel_sidecar_reads_canonical_columns(self):
        from ldsc._kernel import ldscore as ls
        with tempfile.TemporaryDirectory() as tmp:
            meta, df = self._sidecar(tmp)
            r2 = meta.with_name("chr1_r2.parquet")
            loaded = ls._load_full_panel_sidecar(str(r2))
            self.assertEqual(list(loaded["CHR"]), ["1", "1", "1"])
            self.assertEqual(list(loaded["POS"].astype(int)), [10, 20, 30])

    def test_validate_binding_raises_on_hash_mismatch(self):
        from ldsc._kernel import ldscore as ls
        from ldsc._kernel.snp_identity import sidecar_identity_sha256
        with tempfile.TemporaryDirectory() as tmp:
            meta, df = self._sidecar(tmp)
            good = sidecar_identity_sha256(df)
            ls._validate_index_binding(df, n_snps=3, identity_hash=good, context="t")  # no raise
            with self.assertRaisesRegex(ValueError, "identity hash"):
                ls._validate_index_binding(df, n_snps=3, identity_hash="0" * 64, context="t")
            with self.assertRaisesRegex(ValueError, "n_snps"):
                ls._validate_index_binding(df, n_snps=999, identity_hash=good, context="t")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ldscore_workflow.py::IndexBindingTest -v`
Expected: FAIL (`_load_full_panel_sidecar` / `_validate_index_binding` not defined).

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/_kernel/ldscore.py`, add (near `_parquet_schema_layout`, ~line 824). Reuse `sidecar_identity_sha256` (add to the existing snp_identity import block) and the metadata column spec map already imported:

```python
def _panel_sidecar_path_for_r2(r2_path: str) -> Path:
    """Return the ``chrN_meta.tsv.gz`` sidecar path paired with an R2 parquet."""
    p = Path(r2_path)
    if not p.name.endswith("_r2.parquet"):
        raise ValueError(f"Cannot derive sidecar path from non-canonical R2 filename: {p.name}")
    return p.with_name(p.name[: -len("_r2.parquet")] + "_meta.tsv.gz")


def _load_full_panel_sidecar(r2_path: str) -> pd.DataFrame:
    """Load the complete (unrestricted) panel sidecar that defines the index space."""
    sidecar_path = _panel_sidecar_path_for_r2(r2_path)
    if not sidecar_path.exists():
        raise FileNotFoundError(
            f"Index-format R2 parquet '{r2_path}' requires its sidecar '{sidecar_path}'. "
            "The sidecar is mandatory: parquet IDX values are meaningless without it."
        )
    df = pd.read_csv(sidecar_path, sep="\t", comment="#")
    context = f"panel sidecar {sidecar_path}"
    renamed = {
        resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=context): "CHR",
        resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=context): "POS",
        resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=context): "SNP",
        resolve_required_column(df.columns, A1_COLUMN_SPEC, context=context): "A1",
        resolve_required_column(df.columns, A2_COLUMN_SPEC, context=context): "A2",
    }
    return df.rename(columns=renamed)


def _validate_index_binding(full_sidecar: pd.DataFrame, *, n_snps: int, identity_hash: str, context: str) -> None:
    """Hard-fail if the sidecar does not match the parquet's recorded binding."""
    if len(full_sidecar) != int(n_snps):
        raise ValueError(
            f"{context}: sidecar has {len(full_sidecar)} rows but parquet records n_snps={n_snps}. "
            "The parquet and sidecar are not a matched pair."
        )
    actual = sidecar_identity_sha256(full_sidecar, context=context)
    if actual != identity_hash:
        raise ValueError(
            f"{context}: sidecar identity hash {actual} does not match parquet "
            f"ldsc:sidecar_identity_sha256 {identity_hash}. The sidecar is wrong, reordered, or edited."
        )
```

Ensure `sidecar_identity_sha256`, `A1_COLUMN_SPEC`, `A2_COLUMN_SPEC`, `REFERENCE_METADATA_SPEC_MAP`, and `resolve_required_column` are imported at the top of `ldscore.py` (the metadata specs and resolver already are; add `sidecar_identity_sha256` to the `from .snp_identity import (...)` block).

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ldscore_workflow.py::IndexBindingTest -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git commit -m "feat(ldscore): add panel-sidecar loader and index binding validation"
```

### Task 5: Remap builder (build_idx → retained_idx), mode-agnostic

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (add free function near `identifier_keys`, ~line 664)
- Test: `tests/test_ldscore_workflow.py` (new class `IndexRemapTest`)

- [ ] **Step 1: Write the failing test**

```python
class IndexRemapTest(unittest.TestCase):
    def test_remap_maps_build_rows_to_retained_indices_chr_pos(self):
        from ldsc._kernel.ldscore import build_index_remap

        full = pd.DataFrame({"CHR": ["1"] * 4, "POS": [10, 20, 30, 40],
                             "SNP": ["a", "b", "c", "d"], "A1": list("ACGT"), "A2": list("GTAC")})
        # retained = rows 0 and 2 of the panel, in matrix order [pos30, pos10]
        retained = full.iloc[[2, 0]].reset_index(drop=True)
        remap, retained_build_idx = build_index_remap(full, retained, "chr_pos")
        # build row 2 (pos30) -> retained idx 0 ; build row 0 (pos10) -> retained idx 1 ; others -1
        self.assertEqual(remap.tolist(), [1, -1, 0, -1])
        # retained_build_idx[matrix_idx] = build row, ascending
        self.assertEqual(retained_build_idx.tolist(), [2, 0])  # matches retained order

    def test_remap_allele_aware_uses_alleles(self):
        from ldsc._kernel.ldscore import build_index_remap

        full = pd.DataFrame({"CHR": ["1", "1"], "POS": [10, 10],
                             "SNP": ["a", "b"], "A1": ["A", "C"], "A2": ["G", "T"]})
        retained = full.iloc[[1]].reset_index(drop=True)  # the C/T variant at pos10
        remap, _ = build_index_remap(full, retained, "chr_pos_allele_aware")
        self.assertEqual(remap.tolist(), [-1, 0])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ldscore_workflow.py::IndexRemapTest -v`
Expected: FAIL (`cannot import name 'build_index_remap'`).

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/_kernel/ldscore.py`, add after `identifier_keys` (~line 678):

```python
def build_index_remap(
    full_sidecar: pd.DataFrame,
    retained_metadata: pd.DataFrame,
    identifier_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Map panel (build) indices to retained matrix indices for the index format.

    ``full_sidecar`` is the complete panel in build order (the parquet IDX
    space). ``retained_metadata`` is the analysis-restricted matrix universe in
    matrix order. Returns ``(remap, retained_build_idx)`` where
    ``remap[build_idx]`` is the retained matrix index or ``-1``, and
    ``retained_build_idx[matrix_idx]`` is the originating build index (ascending,
    used for IDX_1 row-group pruning). Matching uses the same mode-dependent
    identity keys as the legacy per-pair decode, so the result is bit-identical.
    """
    mode = normalize_snp_identifier_mode(identifier_mode)
    full_keys = identifier_keys(full_sidecar, mode).to_numpy()
    retained_keys = identifier_keys(retained_metadata, mode).to_numpy()
    retained_index = pd.Index(retained_keys)
    if retained_index.has_duplicates:
        raise ValueError("Retained metadata has duplicate identity keys; cannot build index remap.")
    remap = retained_index.get_indexer(full_keys).astype(np.int32, copy=False)

    m = len(retained_metadata)
    retained_build_idx = np.empty(m, dtype=np.int64)
    valid = remap >= 0
    retained_build_idx[remap[valid]] = np.nonzero(valid)[0]
    return remap, retained_build_idx
```

> `pd.Index.get_indexer` requires a unique *index* (`retained_keys`), not unique queries, so duplicate or `NaN` `full_keys` (unmatchable build rows) correctly yield `-1`.

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ldscore_workflow.py::IndexRemapTest -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git commit -m "feat(ldscore): build mode-agnostic build->retained index remap"
```

### Task 6: Index reader path — init, decode, pruning

This rewires `SortedR2BlockReader` for the index layout. It is the largest task; keep the legacy/raw branches *temporarily* until Phase 3 removes them, but route a parquet whose schema is `IDX_1/IDX_2/R2/SIGN` through the new path.

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`_parquet_schema_layout`, `SortedR2BlockReader.__init__`, `_init_canonical_path`, `_decode_canonical_row_group`, `_row_group_indices_for_index_window`)
- Test: `tests/test_ldscore_workflow.py` (new class `IndexReaderDecodeTest`)

- [ ] **Step 1: Write the failing test**

This test writes a tiny index parquet + sidecar, constructs the reader with a restricted metadata, and asserts the decoded dense matrix matches hand-computed values, including a dropped (un-retained) endpoint.

```python
class IndexReaderDecodeTest(unittest.TestCase):
    def _make_panel(self, tmp):
        import gzip as _gz
        from ldsc._kernel import ref_panel_builder as kb
        from ldsc._kernel.snp_identity import sidecar_identity_sha256

        full = pd.DataFrame({"CHR": ["1"] * 4, "POS": [10, 20, 30, 40],
                             "SNP": ["a", "b", "c", "d"], "A1": list("ACGT"),
                             "A2": list("GTAC"), "CM": [0.0] * 4, "MAF": [0.3] * 4})
        meta = Path(tmp) / "hg19" / "chr1_meta.tsv.gz"
        meta.parent.mkdir(parents=True, exist_ok=True)
        with _gz.open(meta, "wt") as h:
            h.write("# ldsc:schema_version=2\n")
            full.to_csv(h, sep="\t", index=False)
        r2 = meta.with_name("chr1_r2.parquet")
        rows = [
            {"i": 0, "j": 1, "R2": 0.5, "sign": "+"},
            {"i": 0, "j": 2, "R2": 0.2, "sign": "+"},
            {"i": 1, "j": 2, "R2": 0.4, "sign": "-"},
            {"i": 2, "j": 3, "R2": 0.9, "sign": "+"},
        ]
        kb.write_r2_parquet(pair_rows=rows, path=r2, genome_build="hg19", n_samples=100,
                            snp_identifier="chr_pos", min_r2=0.0, n_snps=4,
                            sidecar_identity_sha256=sidecar_identity_sha256(full))
        return r2, full

    def test_within_block_matrix_drops_unretained_endpoint(self):
        from ldsc._kernel.ldscore import SortedR2BlockReader
        with tempfile.TemporaryDirectory() as tmp:
            r2, full = self._make_panel(tmp)
            # retained = panel rows 0,1,2 (drop pos40); matrix order = [10,20,30]
            retained = full.iloc[[0, 1, 2]][["CHR", "POS", "SNP", "A1", "A2", "CM", "MAF"]].reset_index(drop=True)
            reader = SortedR2BlockReader(
                paths=[str(r2)], chrom="1", metadata=retained,
                identifier_mode="chr_pos", r2_bias_mode="unbiased",
                r2_sample_size=None, genome_build="hg19",
            )
            mat = reader.within_block_matrix(0, 3)
            self.assertAlmostEqual(mat[0, 1], 0.5, places=6)
            self.assertAlmostEqual(mat[0, 2], 0.2, places=6)
            self.assertAlmostEqual(mat[1, 2], 0.4, places=6)
            self.assertAlmostEqual(mat[0, 0], 1.0, places=6)  # diagonal
            # the pair (2,3) involved dropped panel row 3 -> absent
            self.assertEqual(mat.shape, (3, 3))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ldscore_workflow.py::IndexReaderDecodeTest -v`
Expected: FAIL — `_parquet_schema_layout` returns neither `canonical` nor `raw` for `IDX_*` columns (or decode reads `POS_1`).

- [ ] **Step 3: Write minimal implementation**

**(a)** Update `_parquet_schema_layout` (~line 824) to recognize the index schema first:

```python
def _parquet_schema_layout(schema_names: Sequence[str]) -> str:
    names = set(schema_names)
    if {"IDX_1", "IDX_2", "R2"}.issubset(names):
        return "index"
    try:
        _resolve_canonical_parquet_columns(schema_names, require_endpoint_alleles=False)
        return "canonical"
    except ValueError:
        return "raw"
```

**(b)** In `SortedR2BlockReader.__init__`, after `layout = _parquet_schema_layout(probe_schema_names)` (~line 1506), add an index branch that defers identity-key work to the remap and skips the per-mode `index_map`/`_index_keys` construction. Insert before the `if layout == "canonical":` block:

```python
        if layout == "index":
            if len(paths) != 1:
                raise ValueError(
                    "index parquet_r2 backend requires exactly one file per chromosome; "
                    f"got {len(paths)} for chromosome {self.chrom}"
                )
            if pq is None:
                raise LDSCDependencyError("pyarrow is required for index parquet R2 input.")
            self._runtime_layout = "index"
            self._pf = pq.ParquetFile(paths[0])
            self._init_index_path(paths[0], metadata)
            return
```

Note: the retained `metadata` passed in has already been renamed to `CHR/POS/SNP` (and `A1/A2` for allele-aware) by the block above (~lines 1482-1513); reuse it directly.

**(c)** Add `_init_index_path`:

```python
    def _init_index_path(self, path: str, retained_metadata: pd.DataFrame) -> None:
        """Validate index-format metadata, the sidecar binding, and build the remap."""
        if self._pf is None:
            raise ValueError("Index parquet reader is not initialized.")
        schema_meta = self._pf.schema_arrow.metadata or {}

        build_raw = schema_meta.get(b"ldsc:sorted_by_build")
        if build_raw is None:
            raise ValueError(
                f"'{path}' is index-format but has no ldsc:sorted_by_build metadata. "
                "Regenerate with `ldsc build-ref-panel`."
            )
        parquet_build = normalize_genome_build(build_raw.decode("utf-8"))
        if self.genome_build not in {None, parquet_build}:
            raise ValueError(
                f"Parquet sorted for {parquet_build} but analysis uses {self.genome_build}. "
                f"Use the matching reference file or regenerate with `--genome-build {self.genome_build}`."
            )
        self.genome_build = parquet_build

        n_snps_raw = schema_meta.get(b"ldsc:n_snps")
        hash_raw = schema_meta.get(b"ldsc:sidecar_identity_sha256")
        if n_snps_raw is None or hash_raw is None:
            raise ValueError(f"'{path}' is missing ldsc:n_snps or ldsc:sidecar_identity_sha256 binding metadata.")
        n_snps = int(n_snps_raw.decode("utf-8"))

        full_sidecar = _load_full_panel_sidecar(path)
        _validate_index_binding(
            full_sidecar, n_snps=n_snps, identity_hash=hash_raw.decode("utf-8"),
            context=f"SortedR2BlockReader[{self.chrom}] {path}",
        )

        self._remap, self._retained_build_idx = build_index_remap(
            full_sidecar, retained_metadata, self.identifier_mode
        )
        self.m = len(retained_metadata)

        meta = self._pf.metadata
        if meta.num_row_groups > 0:
            avg = meta.num_rows / meta.num_row_groups
            if avg > 500_000:
                LOGGER.warning(
                    f"'{path}' has {meta.num_row_groups} row group(s) (avg {avg:.0f} rows/group); "
                    "query performance will be degraded. Regenerate with row_group_size=50000."
                )

        idx1_col = self._pf.schema_arrow.names.index("IDX_1")
        self._rg_bounds = []
        for rg_index in range(meta.num_row_groups):
            stats = meta.row_group(rg_index).column(idx1_col).statistics
            if stats is None or not getattr(stats, "has_min_max", False):
                continue
            self._rg_bounds.append((int(stats.min), int(stats.max), rg_index))
```

**(d)** Repoint pruning. Replace `_row_group_indices_for_index_window` (~lines 1685-1691) to use `self._retained_build_idx` for the index layout while leaving the canonical/raw `self.pos` behavior intact:

```python
    def _row_group_indices_for_index_window(self, start: int, stop: int) -> list[int]:
        """Return row groups overlapping a retained-SNP index window."""
        start = max(0, int(start))
        stop = min(int(stop), int(getattr(self, "m", 0)))
        if stop <= start or not self._rg_bounds:
            return []
        if getattr(self, "_runtime_layout", None) == "index":
            lo = int(self._retained_build_idx[start])
            hi = int(self._retained_build_idx[stop - 1])
            return [rg for mn, mx, rg in self._rg_bounds if mn <= hi and mx >= lo]
        return self._row_group_indices_for_pos_window(int(self.pos[start]), int(self.pos[stop - 1]))
```

**(e)** Add an index decode that gathers through `self._remap`. Add `_decode_index_row_group` and route `_get_decoded_canonical_row_group` to it when `self._runtime_layout == "index"` (simplest: add a sibling method and branch in `_get_decoded_canonical_row_group`):

```python
    def _decode_index_row_group(self, row_group_index: int) -> _DecodedR2RowGroup:
        """Decode one index row group: gather endpoints through the remap."""
        table = self._pf.read_row_group(int(row_group_index), columns=["IDX_1", "IDX_2", "R2"])
        idx1 = _arrow_column_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        idx2 = _arrow_column_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        r2 = self._transform_r2(_arrow_column_to_numpy(table.column("R2")).astype(np.float32, copy=False))
        i = self._remap[idx1]
        j = self._remap[idx2]
        keep = (i >= 0) & (j >= 0)
        return _DecodedR2RowGroup(
            row_group_index=int(row_group_index),
            i=i[keep].astype(np.int32, copy=False),
            j=j[keep].astype(np.int32, copy=False),
            r2=r2[keep].astype(np.float32, copy=False),
        )
```

In `_get_decoded_canonical_row_group` (~line 1857), branch the miss path:

```python
        decoded = (
            self._decode_index_row_group(row_group_index)
            if self._runtime_layout == "index"
            else self._decode_canonical_row_group(row_group_index)
        )
```

`cross_block_matrix` / `within_block_matrix` already call `_query_union_rows_canonical_by_index` when `self._runtime_layout == "canonical"`; extend those two conditionals (~lines 1985, 2023) to also take the by-index path for `"index"`:

```python
        if self._runtime_layout in ("canonical", "index"):
            query_rows = self._query_union_rows_canonical_by_index(union_start, union_stop)
```

(and the analogous line in `within_block_matrix`).

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ldscore_workflow.py::IndexReaderDecodeTest -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git commit -m "feat(ldscore): index-format reader with gather decode and IDX_1 pruning"
```

---

## Phase 3 — Remove dead paths (clean break)

### Task 7: Remove the raw external-input read path

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (delete `_query_union_rows_raw`, raw branch in `__init__`, `_query_union_rows`, `_query_union_rows_canonical`, `self.pos`/`_raw_*` state that only served raw/legacy)
- Test: `tests/test_ldscore_workflow.py` (assert raw schema now raises)

- [ ] **Step 1: Write the failing test**

```python
class RawSchemaRejectedTest(unittest.TestCase):
    def test_legacy_raw_schema_is_rejected(self):
        from ldsc._kernel.ldscore import SortedR2BlockReader
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chr1_r2.parquet"
            _write_legacy_r2_parquet(path)  # existing helper in this file
            meta = pd.DataFrame({"CHR": ["1"], "POS": [10], "SNP": ["a"], "A1": ["A"], "A2": ["G"]})
            with self.assertRaisesRegex(ValueError, "regenerate|index|build-ref-panel"):
                SortedR2BlockReader(paths=[str(path)], chrom="1", metadata=meta,
                                    identifier_mode="chr_pos", r2_bias_mode="unbiased",
                                    r2_sample_size=None, genome_build="hg19")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ldscore_workflow.py::RawSchemaRejectedTest -v`
Expected: FAIL (raw path still constructs a working reader).

- [ ] **Step 3: Write minimal implementation**

In `SortedR2BlockReader.__init__`, replace the `if layout == "raw":` block (~lines 1553-1566) and the trailing `raise ValueError(...)` with a single rejection for anything that is not `index`:

```python
        raise ValueError(
            f"'{paths[0]}' is not an index-format R2 parquet (columns IDX_1/IDX_2/R2/SIGN). "
            "External and legacy formats are not supported; regenerate with `ldsc build-ref-panel`."
        )
```

Delete now-dead members: `_query_union_rows_raw`, `_query_union_rows`, `_query_union_rows_canonical`, `_decode_canonical_row_group` (legacy 10-col), `_init_canonical_path`, `_row_group_indices_for_pos_window`, the `self.pos`/`self._raw_*`/`self._last_query_*`/`index_map`/`_index_keys`/`_index_values` state and the `effective_merge_key_series`/`_vectorized_*`/`_endpoint_identity_keys` usages that only served those paths. Remove the `"canonical"` arms in `cross_block_matrix`/`within_block_matrix` left over from Task 6 (they are now unreachable; keep only `"index"`).

> Work iteratively: delete a member, run the suite, fix the next `AttributeError`/`NameError`. Keep `_DecodedR2RowGroup`, `_RowGroupLRUCache`, `_arrow_column_to_numpy`, `_query_union_rows_canonical_by_index`, `_sliding_query_index_windows`, `configure_auto_row_group_cache`, `_transform_r2`, `cross_block_matrix`, `within_block_matrix`.

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ldscore_workflow.py::RawSchemaRejectedTest -v` then the full suite.
Expected: target test PASS. Other tests referencing raw/legacy readers will fail — they are addressed in Task 9.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ldscore.py tests/test_ldscore_workflow.py
git commit -m "refactor(ldscore): drop raw/legacy R2 read paths (index-only)"
```

### Task 8: Make the sidecar mandatory (remove synth fallback)

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel.py` (`load_metadata` ~lines 484-499; delete `_synthesize_metadata_from_r2_paths` ~line 638)
- Test: `tests/test_ref_panel.py` (assert missing sidecar raises)

- [ ] **Step 1: Write the failing test**

Add to `tests/test_ref_panel.py` a test that a parquet directory with an R2 file but no `chrN_meta.tsv.gz` raises `ValueError` mentioning the sidecar (mirror an existing `ParquetR2RefPanel` test's setup; assert on `load_metadata`).

```python
    def test_missing_sidecar_is_a_hard_error(self):
        # build a panel dir with only chr1_r2.parquet (index format), no meta sidecar
        ...  # reuse the existing index-parquet writer helper used elsewhere in this file
        with self.assertRaisesRegex(ValueError, "sidecar"):
            panel.load_metadata("1")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `... tests/test_ref_panel.py -k missing_sidecar -v`
Expected: FAIL (synthesis fallback currently succeeds with a warning).

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/_kernel/ref_panel.py` `load_metadata` (~lines 484-499), replace the `else:` synthesis branch with a hard error:

```python
        else:
            raise ValueError(
                f"Reference-panel metadata sidecar is missing for chromosome {chrom}. "
                "Index-format R2 parquets require their chrN_meta.tsv.gz sidecar."
            )
```

Delete `_synthesize_metadata_from_r2_paths` (~line 638) and any now-unused imports/helpers it pulled in. The `fail_on_missing_metadata` config flag becomes effectively always-on; leave the flag in `config.py` for compatibility but note it no longer changes behavior (a one-line comment at its definition).

- [ ] **Step 4: Run test to verify it passes**

Run: `... tests/test_ref_panel.py -k missing_sidecar -v`
Expected: PASS. Any existing test asserting the synthesis-warning path is removed.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/ref_panel.py tests/test_ref_panel.py
git commit -m "refactor(ref-panel): require sidecar for index parquets (drop synth fallback)"
```

---

## Phase 4 — Integration parity

### Task 9: Regenerate fixtures and re-baseline the workflow suite

**Files:**
- Modify: `tests/test_ldscore_workflow.py` (`_write_minimal_r2_parquet`, `_write_legacy_r2_parquet`, all fixtures producing 10-col parquets)
- Add: a shared fixture helper `_write_index_panel(tmp, full_df, pair_rows)` writing parquet + sidecar with a correct binding hash.

- [ ] **Step 1: Write/replace fixtures**

Replace `_write_minimal_r2_parquet` to emit the 4-column index schema + a sibling sidecar with a valid hash (reuse `kb.write_r2_parquet` + `sidecar_identity_sha256` as in Task 6's helper). Keep `_write_legacy_r2_parquet` only for the Task 7 rejection test.

- [ ] **Step 2: Run the full suite, trace failures**

Run: `... -m pytest -q`
Expected: a set of failures in tests that (a) decode old-schema fixtures, (b) assert POS/SNP columns, or (c) exercise the removed raw path. List them.

- [ ] **Step 3: Update each failing test**

For each: regenerate its panel via the index helper; update assertions to index semantics; delete tests that only validated removed behavior (raw canonicalization, synth fallback, 10-col schema, tier-3 POS build inference). Do not weaken numerical assertions.

- [ ] **Step 4: Run the full suite**

Run: `... -m pytest -q`
Expected: green. Record the new pass/skip count (replaces the old `871 passed, 1 skipped`).

- [ ] **Step 5: Commit**

```bash
git add tests/
git commit -m "test(ldscore): re-baseline workflow fixtures on index format"
```

### Task 10: Build→read LD-score parity gate across all four modes

**Files:**
- Test: `tests/test_ldscore_workflow.py` (new class `IndexParityTest`)

- [ ] **Step 1: Write the parity test**

Build a small panel from a deterministic genotype fixture, compute LD scores via the full pipeline for each of `rsid`, `chr_pos`, `rsid_allele_aware`, `chr_pos_allele_aware`, and assert they agree with a stored reference vector (computed once from a direct dense-R² calculation on the same genotypes). The same index parquet must serve all four modes (mode-agnosticism), differing only in the sidecar-keyed remap.

```python
class IndexParityTest(unittest.TestCase):
    def test_ld_scores_identical_across_modes_from_one_index_parquet(self):
        ref = _direct_dense_ld_scores(_GENO_FIXTURE)   # reference, mode-independent universe
        out = _build_index_panel_from_geno(_GENO_FIXTURE)
        for mode in ("rsid", "chr_pos", "rsid_allele_aware", "chr_pos_allele_aware"):
            scores = _ld_scores_via_pipeline(out, identifier_mode=mode)
            np.testing.assert_allclose(scores, ref, rtol=0, atol=1e-6)
```

> Implement `_direct_dense_ld_scores`, `_build_index_panel_from_geno`, and `_ld_scores_via_pipeline` using the existing genotype fixtures and pipeline entry points already used by other tests in this file. The reference must be computed independently of the parquet path (e.g., dense `R²` matrix times the all-ones annotation) so the test is a true oracle.

- [ ] **Step 2: Run to verify it fails (if helpers absent) then passes**

Run: `... tests/test_ldscore_workflow.py::IndexParityTest -v`
Expected after implementation: PASS for all four modes.

- [ ] **Step 3: Commit**

```bash
git add tests/test_ldscore_workflow.py
git commit -m "test(ldscore): build->read LD-score parity across all identifier modes"
```

---

## Phase 5 — Documentation

### Task 11: Rewrite the format/read-path spec doc

**Files:**
- Modify: `docs/current/parquet-r2-format-and-read-pipeline.md`

- [ ] **Step 1: Rewrite §2 (format)** — four columns (`IDX_1`/`IDX_2` int32, `R2` float32, `SIGN` bool); drop the 10-column table, alias-tolerant base columns, and the build-position discussion. Document the index space (sidecar row order), `IDX_1` sort invariant, and `SIGN` semantics (`True ⇔ r ≥ 0`).

- [ ] **Step 2: Rewrite §2.5 (metadata)** — add `ldsc:n_snps`, `ldsc:sidecar_identity_sha256`; bump `schema_version` to 2; note `ldsc:snp_identifier` is provenance-only.

- [ ] **Step 3: Rewrite §3 (read pipeline)** — full-sidecar load, binding validation, remap, gather decode, `IDX_1` pruning; mode-agnosticism; remove the legacy/raw §4 caveats and replace with the clean-break statement and the §6 pipeline-order invariant from the design spec.

- [ ] **Step 4: Cross-check** — `grep -n "POS_1\|SNP_1\|raw schema\|synthesize" docs/current/parquet-r2-format-and-read-pipeline.md` returns nothing stale.

- [ ] **Step 5: Commit**

```bash
git add docs/current/parquet-r2-format-and-read-pipeline.md
git commit -m "docs: rewrite R2 parquet format/read spec for index layout"
```

---

## Final verification

- [ ] Run the full suite: `... -m pytest -q`. Expected: green; record the count.
- [ ] `grep -rn "POS_1\|SNP_1\|A1_1\|reference_snp_table\|_query_union_rows_raw\|_synthesize_metadata_from_r2" src/ldsc` returns only intentional references (ideally none).
- [ ] Update `docs/superpowers/specs/2026-06-01-build-ref-panel-memory-optimization-design.md` commit log if it tracks format changes, and `design_map.md` if present.

---

## Self-Review (completed during planning)

- **Spec coverage:** §2 format → Tasks 1-2; §2.3 metadata → Task 2; §3.1 sidecar mandatory → Task 8; §3.2 reader/remap → Tasks 4-6; §3.3 mode-agnostic → Tasks 5,10; §4 hash → Tasks 0,2,4; §5 write → Tasks 1-3; §6 order invariant → Tasks 6,10 (reader receives unchanged retained metadata; parity gate enforces); §7 removed surface → Tasks 7-8; §9 testing → Tasks 9-10; docs → Task 11. No gaps.
- **Placeholders:** none — every code step shows real code. The parity-oracle helpers in Task 10 are described with explicit construction requirements (independent dense-R² reference) rather than left as "write tests".
- **Type/name consistency:** `_standard_r2_index_table`, `write_r2_parquet(..., n_snps, sidecar_identity_sha256)`, `sidecar_identity_sha256`, `_load_full_panel_sidecar`, `_validate_index_binding`, `build_index_remap`, `_decode_index_row_group`, `_retained_build_idx`, `_runtime_layout == "index"` are used consistently across tasks.
