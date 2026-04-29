# Design: Chromosome-Aligned Parquet Row Groups

## Problem

`baseline.parquet` and `query.parquet` were previously written by
`LDScoreDirectoryWriter` through `pandas.DataFrame.to_parquet()`. Pandas
delegates row group layout to pyarrow, which produces a single large row group
by default. Downstream code that wants only one chromosome must then read the
entire file.

For a run with 50 annotation columns and ~800 k SNPs at float32, `baseline.parquet`
is roughly 160 MB. Loading the full file to use 4 % of it (chromosome 22,
~34 k SNPs) wastes I/O and memory.

---

## Goal

Allow callers to read all SNPs for one chromosome without reading any other
chromosome's data. The solution must not change the public output directory
structure or break existing readers that load the full file.

---

## Solution: One Row Group Per Chromosome

Write each chromosome as a separate pyarrow row group inside the existing flat
`baseline.parquet` and `query.parquet` files. A parquet file with
chromosome-aligned row groups lets pyarrow skip entire row groups at the I/O
level based on column statistics — no full-file scan required.

### Invariant

> Every row group in `baseline.parquet` and `query.parquet` contains rows from
> exactly one chromosome. No row group mixes chromosomes.

The simplest implementation is **one row group per chromosome**. LDSC chromosome
sizes are bounded (~34 k–80 k SNPs), which are well below pyarrow's 64 MB
default row group target even at 50+ annotation columns. Splitting large
chromosomes into multiple row groups is not needed and is out of scope.

---

## Why Not a Partitioned Dataset Layout?

A `baseline/CHR=1/part-0.parquet` Hive-style layout achieves the same I/O
savings but:

- changes the public directory contract (`baseline.parquet` → a subdirectory),
- requires the pyarrow Dataset API on the read side,
- complicates manual inspection, rsync, and downstream tooling.

Chromosome-aligned row groups in a single flat file provide identical read
efficiency with no structural change.

---

## Data Contract Changes

### Parquet Files (no public API change)

`baseline.parquet` and `query.parquet` continue to be single flat files. The
only change is internal: the row group layout is now guaranteed to be
chromosome-aligned. Existing readers that call `pd.read_parquet()` or
`pyarrow.parquet.read_table()` on the full file are unaffected.

### `manifest.json` (additive change)

Three new fields are added:

```json
{
  "row_group_layout": "one_per_chromosome",
  "baseline_row_groups": [
    {"chrom": "1",  "row_group_index": 0, "row_offset": 0,      "n_rows": 78423},
    {"chrom": "2",  "row_group_index": 1, "row_offset": 78423,  "n_rows": 72841},
    {"chrom": "22", "row_group_index": 21, "row_offset": 756190, "n_rows": 33810}
  ],
  "query_row_groups": null
}
```

- `row_group_layout` — version tag. Readers must assert `"one_per_chromosome"`;
  unknown values should raise a clear error.
- `baseline_row_groups` / `query_row_groups` — ordered list of per-row-group
  metadata. `query_row_groups` is `null` when no `query.parquet` was written.
  Each entry: `chrom` (str), `row_group_index` (int, 0-based), `row_offset`
  (int, cumulative SNP offset to the first row of this group),
  `n_rows` (int, SNP count in this group).
- `row_offset + n_rows` for the last entry must equal `n_baseline_rows` (already
  recorded). This is a verifiable consistency check.

These fields are purely additive. Old readers that parse only `files`,
`chromosomes`, `baseline_columns`, etc. are unaffected.

---

## Writer Implementation

`LDScoreDirectoryWriter.write()` in `outputs.py` replaces the two
`pandas.DataFrame.to_parquet()` calls with `_write_chromosome_aligned_parquet`.

```python
import pyarrow as pa
import pyarrow.parquet as pq

def _write_chromosome_aligned_parquet(
    df: pd.DataFrame, path: Path, compression: str | None
) -> list[dict]:
    """Write one row group per chromosome; return row-group metadata for manifest."""
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

`groupby("CHR", sort=False)` preserves the (CHR, POS) order that
`_aggregate_chromosome_results` already establishes via
`sort_frame_by_genomic_position`. It does not re-sort.

`build_manifest` receives `baseline_rg` and `query_rg` as new arguments and
adds the three manifest fields before returning.

---

## Reader-Side API

### Tier 1 — Full chromosome load via pyarrow Dataset API

No manifest needed. Pyarrow uses column statistics for predicate pushdown and
skips non-matching row groups automatically.

```python
import pyarrow.dataset as ds

table = ds.dataset("path/to/baseline.parquet").to_table(
    filter=ds.field("CHR") == "1"
)
df = table.to_pandas()
```

### Tier 2 — Exact row group read via manifest lookup

Zero per-row predicate evaluation. Reads exactly one row group by index.

```python
import json
import pyarrow.parquet as pq
from pathlib import Path

ldscore_dir = Path("path/to/ldscore_dir")
manifest = json.loads((ldscore_dir / "manifest.json").read_text())
assert manifest["row_group_layout"] == "one_per_chromosome"

rg_by_chrom = {e["chrom"]: e["row_group_index"]
               for e in manifest["baseline_row_groups"]}

pf = pq.ParquetFile(ldscore_dir / "baseline.parquet")
df = pf.read_row_group(rg_by_chrom["1"]).to_pandas()
```

### Tier 3 — SNP-list filtering within a chromosome

Load the chromosome first (Tier 1 or 2), then filter in memory. The
chromosome-level read already skips ~96 % of the file; in-memory SNP filtering
is negligible.

```python
target_snps = {"rs1234", "rs5678"}
df = df[df["SNP"].isin(target_snps)]
```

---

## Tests

All new tests live in `tests/test_output.py`, extending the existing
`LDScoreDirectoryWriterTest` class. A new multi-chromosome fixture
`make_multi_chrom_result` is added.

### Test 1 — `test_one_row_group_per_chromosome`

Assert `pf.metadata.num_row_groups == len(chromosomes)` and that each row group
contains rows from exactly one chromosome.

### Test 2 — `test_chromosome_read_via_row_group_index_excludes_other_chromosomes`

Read row group 0 by index; assert only chromosome 1 rows are returned and
chromosomes 2–22 are absent. This proves the read is truly row-group-scoped.

### Test 3 — `test_manifest_row_group_metadata_is_consistent`

Parse the manifest; verify cumulative `row_offset + n_rows` matches
`n_baseline_rows`. Verify `query_row_groups` matches `n_query_rows`. Verify
`row_group_layout == "one_per_chromosome"`.

### Test 4 — `test_query_row_group_null_when_no_query_table`

Baseline-only result: assert `manifest["query_row_groups"] is None`.

### Test 5 — `test_existing_tests_still_pass`

No new assertions needed; the existing test suite is the regression guard.
The existing `make_split_ldscore_result` fixture (single chromosome) must
continue to pass — it verifies that the new writer produces one row group for
the single-chromosome case.

---

## Out of Scope

- Partitioned dataset layout.
- Splitting large chromosomes into multiple row groups.
- A public Python reader helper function (the manifest + pyarrow API is
  sufficient; a helper can be added later without a design change).
- Changes to any file other than `outputs.py` and `tests/test_output.py`.
