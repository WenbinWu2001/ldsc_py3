# Implementation Plan: Enforce Column Schema

Enforce the column schema defined in `docs/design/column-schema.md` across the
codebase. All changes are mechanical renames or reorderings; no algorithm or
data-model changes are involved.

---

## Changes Required

### File 1: `src/ldsc/column_inference.py`

**Change 1a — `INTERNAL_SUMSTATS_ARTIFACT_SPECS` tuple order**

Current (lines 179–188):
```python
INTERNAL_SUMSTATS_ARTIFACT_SPECS = (
    _canonical_only_spec("SNP", "munged sumstats SNP identifier"),
    _canonical_only_spec("CHR", "munged sumstats chromosome"),
    _canonical_only_spec("POS", "munged sumstats position"),
    _canonical_only_spec("N", ...),
    ...
)
```

Target: move `CHR` and `POS` before `SNP`:
```python
INTERNAL_SUMSTATS_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", "munged sumstats chromosome"),
    _canonical_only_spec("POS", "munged sumstats position"),
    _canonical_only_spec("SNP", "munged sumstats SNP identifier"),
    _canonical_only_spec("N", ...),
    ...
)
```

**Change 1b — `INTERNAL_LDSCORE_ARTIFACT_SPECS` tuple order**

Current (lines 204–210):
```python
INTERNAL_LDSCORE_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", ...),
    _canonical_only_spec("SNP", ...),
    _canonical_only_spec("POS", ...),
    _canonical_only_spec("CM", ...),
    _canonical_only_spec("MAF", ...),
)
```

Target: swap `POS` before `SNP`:
```python
INTERNAL_LDSCORE_ARTIFACT_SPECS = (
    _canonical_only_spec("CHR", ...),
    _canonical_only_spec("POS", ...),
    _canonical_only_spec("SNP", ...),
    _canonical_only_spec("CM", ...),
    _canonical_only_spec("MAF", ...),
)
```

Note: `INTERNAL_ANNOT_ARTIFACT_SPECS` is already correct (`CHR, POS, SNP, CM,
MAF`) — do not change it.

---

### File 2: `src/ldsc/_kernel/ldscore.py`

**Change 2a — `ANNOT_META_COLUMNS` constant (line 238)**

Current:
```python
ANNOT_META_COLUMNS = ("CHR", "SNP", "POS", "CM", "MAF")
```
Target:
```python
ANNOT_META_COLUMNS = ("CHR", "POS", "SNP", "CM", "MAF")
```

**Change 2b — Explicit column reorder list (line ~1904)**

Current:
```python
out = out.loc[:, [col for col in ["CHR", "SNP", "POS", "CM", "MAF"] if col in out.columns] + ...]
```
Target:
```python
out = out.loc[:, [col for col in ["CHR", "POS", "SNP", "CM", "MAF"] if col in out.columns] + ...]
```

**Change 2c — Module-level docstring (lines ~67, ~69)**

Current:
```
- ``baseline.parquet``, containing ``CHR``, ``SNP``, ``POS``, ``regression_ld_scores``, ...
- optional ``query.parquet``, containing ``CHR``, ``SNP``, ``POS``, and query ...
```
Target:
```
- ``baseline.parquet``, containing ``CHR``, ``POS``, ``SNP``, ``regression_ld_scores``, ...
- optional ``query.parquet``, containing ``CHR``, ``POS``, ``SNP``, and query ...
```

---

### File 3: `src/ldsc/outputs.py`

**Change 3a — `required_baseline` validation list (line ~133)**

Current:
```python
required_baseline = ["CHR", "SNP", "POS", "regression_ld_scores", *baseline_columns]
```
Target:
```python
required_baseline = ["CHR", "POS", "SNP", "regression_ld_scores", *baseline_columns]
```

**Change 3b — `required_query` validation list (line ~143)**

Current:
```python
required_query = ["CHR", "SNP", "POS", *query_columns]
```
Target:
```python
required_query = ["CHR", "POS", "SNP", *query_columns]
```

Note: These lists are used only for presence-checking, so the reorder has no
behavioral effect. Update them anyway for consistency with the schema doc.

---

### File 4: `src/ldsc/_kernel/annotation.py`

**Verify only** — `REQUIRED_ANNOT_COLUMNS` at line 102 is already
`("CHR", "POS", "SNP", "CM")` and the write at line 1072 uses it directly.
Confirm no other reorder lists in this file use `SNP` before `POS`.

---

### File 5: `src/ldsc/outputs.py` — float32 narrowing before parquet writes

Add a private helper and call it at both write sites.

**Change 5a — add helper (insert near top of file, after imports):**
```python
def _cast_parquet_floats(df: pd.DataFrame) -> pd.DataFrame:
    """Narrow float64 columns to float32 before parquet write."""
    float64_cols = [c for c in df.columns if df[c].dtype == np.float64]
    return df.astype({c: np.float32 for c in float64_cols}) if float64_cols else df
```

**Change 5b — apply at the baseline/query write site (line ~94):**

Current:
```python
baseline_table.to_parquet(paths["baseline"], index=False, compression=compression)
if query_table is not None:
    query_table.to_parquet(paths["query"], index=False, compression=compression)
```
Target:
```python
_cast_parquet_floats(baseline_table).to_parquet(paths["baseline"], index=False, compression=compression)
if query_table is not None:
    _cast_parquet_floats(query_table).to_parquet(paths["query"], index=False, compression=compression)
```

---

### File 6: Pairwise R² parquet write site

Locate where the pairwise R² parquet files are written (search for
`to_parquet` in `src/ldsc/_kernel/ref_panel_builder.py` and
`src/ldsc/ref_panel_builder.py`). Before each write, cast `R2` to float32:

```python
df["R2"] = df["R2"].astype(np.float32)
df.to_parquet(path, index=False, ...)
```

If a `_cast_parquet_floats` helper is accessible there (or can be imported from
`outputs.py` after it is added in Change 5a), use it instead of the inline cast.
If the helper is not easily importable, use the inline cast to keep the change
self-contained.

---

## Verification Steps

After each file change, run:

```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3
pytest tests/ -x -q
```

If tests pass, inspect an actual output file to confirm column order:

```python
import pandas as pd
df = pd.read_parquet("path/to/baseline.parquet")
print(list(df.columns[:6]))  # expect: ['CHR', 'POS', 'SNP', ...]
```

---

## Out of Scope

- The Jerry workspace (`ldsc_py3_Jerry/`) is a legacy reference; do not modify it.
- The hm3 curated map file on disk uses lowercase column names (`rsID`, `chr`,
  `hg19_pos`). Leave the file as-is; normalization happens on read.
- In-memory dtypes do not change — `str`, `np.int64`, and `float64` remain the
  canonical in-memory types. The float32 narrowing applies only at the parquet
  write boundary.
- `N`, `N_CAS`, `N_CON`, `NSTUDY` are stored in `.sumstats.gz` (text), never
  in parquet. They are not subject to float32 narrowing anywhere.

---
