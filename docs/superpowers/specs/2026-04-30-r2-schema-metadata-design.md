# R2 Schema Metadata: Sample Size and Bias Mode

**Date:** 2026-04-30
**Scope:** `src/ldsc/_kernel/ref_panel_builder.py`, `src/ldsc/ref_panel_builder.py`,
`src/ldsc/_kernel/ldscore.py`, `src/ldsc/_kernel/ref_panel.py`

---

## Problem

The R2 parquet files produced by `build-ref-panel` always store **unbiased** R2
values (via `_unbiased_r2_from_correlation`), but they carry no record of this
fact or of the sample size `N` used. Downstream LD-score calculation must either
receive this information from the user explicitly or rely on an implicit default
(`r2_bias_mode=None → "unbiased"` at ldscore.py:1990) whose correctness cannot
be verified from the files alone.

Consequences:
- Users sharing a panel must communicate N and bias mode out-of-band.
- A panel built with a future `--store-raw-r2` flag would be silently
  misinterpreted as unbiased by any reader that defaults to `"unbiased"`.
- N is lost after the build run even though it is available as `geno.n`.

---

## Design Decisions

| Question | Decision |
|---|---|
| Storage location | Arrow schema metadata on each R2 parquet file (same pattern as `ldsc:sorted_by_build` and `ldsc:row_group_size`) |
| Keys added | `ldsc:n_samples` (integer, UTF-8 encoded) and `ldsc:r2_bias` (`"unbiased"` or `"raw"`) |
| Current builder always writes | `ldsc:r2_bias = "unbiased"` — no new flag needed now |
| Reader helper location | `_kernel/ref_panel.py` — `_read_r2_schema_meta(path)` |
| Auto-load call sites | (1) CLI validation block in `_kernel/ldscore.py`; (2) `ParquetR2RefPanel.load_r2()` in `_kernel/ref_panel.py` |
| `SortedR2BlockReader` | **Unchanged** — receives fully resolved `r2_bias_mode` and `r2_sample_size` as before |

---

## New Schema Metadata Keys

```
ldsc:n_samples   integer encoded as UTF-8 string, e.g. b"3202"
ldsc:r2_bias     b"unbiased" | b"raw"
```

Added in `write_r2_parquet()` alongside the existing keys:

```python
pa_meta = {
    b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
    b"ldsc:row_group_size":  str(row_group_size).encode("utf-8"),
    b"ldsc:n_samples":       str(n_samples).encode("utf-8"),   # NEW
    b"ldsc:r2_bias":         b"unbiased",                       # NEW (always for now)
}
```

`write_r2_parquet` gains one new required parameter: `n_samples: int`.
`_build_chromosome` passes `geno.n` as `n_samples`.

---

## Reader Helper: `_read_r2_schema_meta()`

New function in `_kernel/ref_panel.py`:

```python
from dataclasses import dataclass

@dataclass
class _R2SchemaMeta:
    n_samples: int | None   # None if key absent
    r2_bias: str | None     # "unbiased" | "raw" | None if key absent

def _read_r2_schema_meta(path: str) -> _R2SchemaMeta:
    """Read ldsc:n_samples and ldsc:r2_bias from one R2 parquet file's schema."""
```

Uses `pyarrow.parquet.read_schema(path).metadata`. Returns `_R2SchemaMeta` with
`None` fields when a key is absent (legacy file).

---

## Resolution Rule

Applied at **both** auto-load call sites using a shared helper
`_resolve_r2_bias_from_meta(r2_bias_mode, r2_sample_size, meta)`:

| stored `ldsc:r2_bias` | user `r2_bias_mode` | user `r2_sample_size` | result |
|---|---|---|---|
| `"unbiased"` | `None` | not set | use stored → `"unbiased"`, N not needed |
| `"unbiased"` | `None` | set | use stored → `"unbiased"`, WARN N unused |
| `"unbiased"` | `"unbiased"` | any | consistent, N not needed |
| `"unbiased"` | `"raw"` | any | user override wins; apply correction |
| `"raw"` | `None` | not set | use stored → `"raw"`, auto-fill N from `ldsc:n_samples` |
| `"raw"` | `None` | set | use stored → `"raw"`, use user N (override) |
| `"raw"` | `"raw"` | not set | auto-fill N from `ldsc:n_samples` |
| `"raw"` | `"raw"` | set | use user N (override) |
| absent | `None` | not set | default to `"unbiased"` (existing behavior) |
| absent | `None` | set | WARN N unused (existing: no stored bias to apply it to) |
| absent (legacy) | `"raw"` | set | apply correction with user N (existing behavior) |
| absent (legacy) | `"raw"` | not set | raise: no N available (existing behavior) |

**One additional rule:** if `ldsc:r2_bias` is absent but `ldsc:n_samples` is
present (malformed file), treat as `"raw"` with a WARNING.

---

## Auto-Load Call Sites

### Site 1 — CLI validation in `_kernel/ldscore.py`

Current code (lines 1987–1993):
```python
if args.r2_table:
    if args.r2_bias_mode is None:
        args.r2_bias_mode = "unbiased"
    if args.r2_bias_mode == "raw" and args.r2_sample_size is None:
        raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
```

Replacement: before the `None → "unbiased"` promotion, read schema metadata from
the first resolved parquet file and apply the resolution rule. The promotion and
the `"raw"` raise remain as a fallback for legacy files.

### Site 2 — `ParquetR2RefPanel.load_r2()` in `_kernel/ref_panel.py`

Read schema metadata from the first resolved parquet path for the chromosome.
Apply the resolution rule to determine `r2_bias_mode` and `r2_sample_size` before
passing them to `SortedR2BlockReader`.

---

## Testing Strategy

### Unit tests — `write_r2_parquet` stores the keys

In `tests/test_ref_panel_builder.py`:
- After calling `write_r2_parquet(..., n_samples=42)`, read back
  `pyarrow.parquet.read_schema(path).metadata` and assert both keys are present
  with correct values.

### Unit tests — `_read_r2_schema_meta`

In `tests/test_ref_panel_builder.py` or a new `tests/test_ref_panel.py`:
- File with both keys → returns correct `_R2SchemaMeta`
- File with neither key (legacy) → returns `_R2SchemaMeta(None, None)`
- File with `n_samples` but no `r2_bias` → returns `_R2SchemaMeta(n, None)` with
  warning logged

### Unit tests — `_resolve_r2_bias_from_meta`

In `tests/test_ldscore_workflow.py`:
- Cover each row of the resolution rule table

### Integration test — CLI auto-load

In `tests/test_ldscore_workflow.py`:
- Build a minimal canonical parquet with `ldsc:r2_bias = "unbiased"` and
  `ldsc:n_samples = 100`; run validation with `r2_bias_mode=None` and verify
  no error is raised and `r2_bias_mode` is resolved to `"unbiased"`
- Build with `ldsc:r2_bias = "raw"` and `ldsc:n_samples = 100`; run with
  `r2_bias_mode=None`; verify `r2_bias_mode = "raw"` and `r2_sample_size = 100`
  are auto-populated without error

### Integration test — `ParquetR2RefPanel.load_r2` auto-load

In `tests/test_ldscore_workflow.py` or `tests/test_ref_panel.py`:
- Panel parquet has `ldsc:n_samples = 200`, `ldsc:r2_bias = "raw"`;
  `RefPanelSpec(r2_dir=..., r2_bias_mode=None, r2_sample_size=None)` → verify
  `SortedR2BlockReader` is constructed with `r2_bias_mode="raw"` and
  `r2_sample_size=200`
