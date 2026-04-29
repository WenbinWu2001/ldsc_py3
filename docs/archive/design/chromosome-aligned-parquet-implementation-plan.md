# Implementation Plan: Chromosome-Aligned Parquet Row Groups

Design reference: `docs/design/chromosome-aligned-parquet-row-groups.md`

Status: implemented.

Only two files change: `src/ldsc/outputs.py` and `tests/test_output.py`.
No other module touches the parquet write path.

---

## Step 1 — Add `_write_chromosome_aligned_parquet` to `outputs.py`

**File:** `src/ldsc/outputs.py`

### 1a — Add `pyarrow` imports

After the existing imports block (currently ends at line 24 with `from ._row_alignment import assert_same_snp_rows`), add:

```python
import pyarrow as pa
import pyarrow.parquet as pq
```

### 1b — Add the helper function

Insert after `_cast_parquet_floats` (currently ends around line 29):

```python
def _write_chromosome_aligned_parquet(
    df: pd.DataFrame, path: Path, compression: str | None
) -> list[dict]:
    """Write one pyarrow row group per chromosome; return row-group metadata list.

    Assumes df is already sorted by (CHR, POS), which _aggregate_chromosome_results
    guarantees via sort_frame_by_genomic_position. groupby(sort=False) preserves
    that order without re-sorting.
    """
    df = _cast_parquet_floats(df)
    schema = pa.Schema.from_pandas(df, preserve_index=False)
    row_group_meta: list[dict] = []
    offset = 0
    with pq.ParquetWriter(path, schema, compression=compression) as writer:
        for chrom, chrom_df in df.groupby("CHR", sort=False):
            writer.write_table(pa.Table.from_pandas(chrom_df, preserve_index=False))
            row_group_meta.append({
                "chrom": str(chrom),
                "row_group_index": len(row_group_meta),
                "row_offset": offset,
                "n_rows": len(chrom_df),
            })
            offset += len(chrom_df)
    return row_group_meta
```

---

## Step 2 — Replace `to_parquet` calls in `LDScoreDirectoryWriter.write`

**File:** `src/ldsc/outputs.py`, currently lines 99–101.

### Current code

```python
compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
_cast_parquet_floats(baseline_table).to_parquet(paths["baseline"], index=False, compression=compression)
if query_table is not None:
    _cast_parquet_floats(query_table).to_parquet(paths["query"], index=False, compression=compression)
manifest = self.build_manifest(result, files={name: path.name for name, path in paths.items() if name != "manifest"})
```

### Replacement

```python
compression = None if output_config.parquet_compression in {None, "none"} else output_config.parquet_compression
baseline_rg = _write_chromosome_aligned_parquet(baseline_table, paths["baseline"], compression)
query_rg = None
if query_table is not None:
    query_rg = _write_chromosome_aligned_parquet(query_table, paths["query"], compression)
manifest = self.build_manifest(
    result,
    files={name: path.name for name, path in paths.items() if name != "manifest"},
    baseline_rg=baseline_rg,
    query_rg=query_rg,
)
```

---

## Step 3 — Update `build_manifest` signature and body

**File:** `src/ldsc/outputs.py`, `build_manifest` method.

### 3a — New signature

```python
def build_manifest(
    self,
    result: Any,
    files: dict[str, str],
    baseline_rg: list[dict] | None = None,
    query_rg: list[dict] | None = None,
) -> dict[str, Any]:
```

### 3b — Add three new keys to the returned dict

After the existing `"n_query_rows"` entry, add:

```python
"row_group_layout": "one_per_chromosome",
"baseline_row_groups": baseline_rg or [],
"query_row_groups": query_rg,
```

The `baseline_rg or []` guard keeps the manifest valid if `build_manifest` is
called directly in tests without providing `baseline_rg`.

---

## Step 4 — Add multi-chromosome test fixture to `test_output.py`

**File:** `tests/test_output.py`

Add after `make_split_ldscore_result` (currently ends at line 71):

```python
def make_multi_chrom_result(chromosomes: list[str] | None = None) -> LDScoreResult:
    """Three-chromosome result for row-group layout tests."""
    if chromosomes is None:
        chromosomes = ["1", "2", "22"]
    rows_per_chrom = {"1": 3, "2": 2, "22": 1}
    baseline_rows = []
    query_rows = []
    snp_idx = 1
    for chrom in chromosomes:
        n = rows_per_chrom.get(chrom, 2)
        for i in range(n):
            baseline_rows.append({
                "CHR": chrom, "SNP": f"rs{snp_idx}", "POS": (i + 1) * 100,
                "regr_weight": 1.0, "base": float(snp_idx),
            })
            query_rows.append({
                "CHR": chrom, "SNP": f"rs{snp_idx}", "POS": (i + 1) * 100,
                "query": float(snp_idx) * 0.5,
            })
            snp_idx += 1
    baseline_table = pd.DataFrame(baseline_rows)
    query_table = pd.DataFrame(query_rows)
    return LDScoreResult(
        baseline_table=baseline_table,
        query_table=query_table,
        count_records=[],
        baseline_columns=["base"],
        query_columns=["query"],
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(),
        chromosome_results=[],
        config_snapshot=GlobalConfig(genome_build="hg38", snp_identifier="rsid"),
    )
```

---

## Step 5 — Add new test cases to `LDScoreDirectoryWriterTest`

**File:** `tests/test_output.py`

Add inside `LDScoreDirectoryWriterTest`:

```python
def test_one_row_group_per_chromosome(self):
    import pyarrow.parquet as pq

    result = make_multi_chrom_result()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "out"
        LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

        for fname in ("baseline.parquet", "query.parquet"):
            pf = pq.ParquetFile(output_dir / fname)
            self.assertEqual(pf.metadata.num_row_groups, 3, f"{fname}: expected 3 row groups")
            for i in range(pf.metadata.num_row_groups):
                chroms = pf.read_row_group(i)["CHR"].unique().to_pylist()
                self.assertEqual(len(chroms), 1, f"{fname} row group {i} mixes chromosomes: {chroms}")

def test_chromosome_read_via_row_group_index_excludes_other_chromosomes(self):
    import pyarrow.parquet as pq

    result = make_multi_chrom_result()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "out"
        LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

        manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
        rg_by_chrom = {e["chrom"]: e["row_group_index"]
                       for e in manifest["baseline_row_groups"]}

        pf = pq.ParquetFile(output_dir / "baseline.parquet")
        df = pf.read_row_group(rg_by_chrom["1"]).to_pandas()
        self.assertTrue((df["CHR"] == "1").all())
        self.assertNotIn("2", df["CHR"].values)
        self.assertNotIn("22", df["CHR"].values)

def test_manifest_row_group_metadata_is_consistent(self):
    result = make_multi_chrom_result()
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "out"
        LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

        manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["row_group_layout"], "one_per_chromosome")

        for field_name, total_field in [
            ("baseline_row_groups", "n_baseline_rows"),
            ("query_row_groups", "n_query_rows"),
        ]:
            entries = manifest[field_name]
            self.assertIsNotNone(entries)
            expected_offset = 0
            for e in entries:
                self.assertEqual(e["row_offset"], expected_offset)
                expected_offset += e["n_rows"]
            self.assertEqual(expected_offset, manifest[total_field])

def test_query_row_groups_null_when_no_query_table(self):
    result = make_split_ldscore_result(query=False)
    with tempfile.TemporaryDirectory() as tmpdir:
        output_dir = Path(tmpdir) / "out"
        LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=output_dir))

        manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
        self.assertIsNone(manifest["query_row_groups"])
        self.assertIsNotNone(manifest["baseline_row_groups"])
        self.assertEqual(len(manifest["baseline_row_groups"]), 1)
```

---

## Step 6 — Verify

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
pytest tests/test_output.py -v
```

All existing tests must still pass. The four new tests must also pass.

After the run, spot-check the row group count on a real output if available:

```python
import pyarrow.parquet as pq
pf = pq.ParquetFile("path/to/baseline.parquet")
print(pf.metadata.num_row_groups)   # expect: number of chromosomes in the run
```

---

## Checklist

- [x] Step 1a — add `pyarrow` imports to `outputs.py`
- [x] Step 1b — add `_write_chromosome_aligned_parquet` helper
- [x] Step 2  — replace `to_parquet` calls in `write()`
- [x] Step 3  — update `build_manifest` signature and body
- [x] Step 4  — add `make_multi_chrom_result` fixture to `test_output.py`
- [x] Step 5  — add four new test cases
- [x] Step 6  — run `pytest tests/test_output.py -v` and confirm all pass

---

## Non-Goals

- No changes to `ldscore_calculator.py` — the sort order it already produces is
  what the writer relies on.
- No changes to regression code (`h2`, `rg`, `partitioned-h2`) — they read the
  full parquet file and are unaffected.
- No changes to the directory structure or filenames.
- No new public API functions beyond what the manifest enables.
