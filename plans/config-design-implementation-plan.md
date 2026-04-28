# Implementation Plan: Immutable Config + Provenance-Carrying Results

## Overview

This plan refactors the LDSC Python package to make genome-build and SNP-identifier
assumptions explicit, immutable, and validated at the moment two result objects are
combined. `GlobalConfig` is already `frozen=True`; the remaining work is:

1. Add a new exception class and a compatibility-validation function to `config.py`
2. Add `config_snapshot: GlobalConfig` to result dataclasses that lack it
3. Propagate the snapshot through all workflow paths that produce those dataclasses
4. Call `validate_config_compatibility()` at every combination point
5. Update `__init__.py` to export `ConfigMismatchError`

**Estimated scope:** ~150–200 lines changed or added across 7 files.
No changes to public API signatures (all changes are additive).
No changes to CLI behavior.

---

## Step 0 — Baseline understanding (read before writing any code)

Read the following files in full before starting:

- `src/ldsc/config.py` — complete file, understand `GlobalConfig`, `_GLOBAL_CONFIG`,
  `get_global_config()`, `set_global_config()`
- `src/ldsc/ldscore_calculator.py` — focus on `ChromLDScoreResult` (lines 62–106),
  `LDScoreResult` (lines 108–169), `LDScoreCalculator.run()` (lines 171–354),
  `run_ldscore()` (lines 558–576)
- `src/ldsc/_kernel/annotation.py` — focus on `AnnotationBundle` (lines 131–189),
  `AnnotationBuilder.run()` method
- `src/ldsc/regression_runner.py` — focus on `RegressionDataset` (lines 44–63),
  `RegressionRunner.__init__()` (line 74), `RegressionRunner.build_dataset()`
- `src/ldsc/sumstats_munger.py` — focus on `SumstatsTable` (lines 89–156),
  `SumstatsMunger.run()`
- `src/ldsc/_kernel/ref_panel.py` — focus on `RefPanelConfig` (lines 53–79),
  `RefPanel.__init__()` (lines 81–142)
- `src/ldsc/outputs.py` — focus on `RunSummary` (lines 65–73) which already has
  `config_snapshot`

---

## Step 1 — Add `ConfigMismatchError` and `validate_config_compatibility()` to `config.py`

**File:** `src/ldsc/config.py`

### 1a. Add `ConfigMismatchError`

Insert after the existing imports and before the `GlobalConfig` class definition
(approximately line 40, before the type aliases or after the `__all__` if present):

```python
class ConfigMismatchError(ValueError):
    """Raised when two result objects carry incompatible GlobalConfig snapshots."""
```

### 1b. Add `validate_config_compatibility()`

Insert after the `reset_global_config()` function (approximately line 136).
The function signature and logic:

```python
def validate_config_compatibility(
    a: "GlobalConfig",
    b: "GlobalConfig",
    context: str = "",
) -> None:
    """
    Raise ConfigMismatchError if a and b differ on critical fields.

    Critical fields: genome_build, snp_identifier.
    Advisory fields: log_level, fail_on_missing_metadata, restrict_snps_path.

    Parameters
    ----------
    a, b : GlobalConfig
        Snapshots to compare (typically from two result objects being combined).
    context : str
        Human-readable description of what is being combined, e.g.
        "LDScoreResult and SumstatsTable". Included in the error message.
    """
    prefix = f" when combining {context}" if context else ""
    if a.genome_build != b.genome_build:
        raise ConfigMismatchError(
            f"genome_build mismatch{prefix}: {a.genome_build!r} vs "
            f"{b.genome_build!r}. These objects were computed under different "
            "genome-build assumptions and cannot be safely merged."
        )
    if a.snp_identifier != b.snp_identifier:
        raise ConfigMismatchError(
            f"snp_identifier mismatch{prefix}: {a.snp_identifier!r} vs "
            f"{b.snp_identifier!r}. These objects were computed under different "
            "SNP-identifier modes and cannot be safely merged."
        )
    if a.restrict_snps_path != b.restrict_snps_path:
        import warnings
        warnings.warn(
            f"restrict_snps_path differs{prefix}: {a.restrict_snps_path!r} vs "
            f"{b.restrict_snps_path!r}. Results are combinable but SNP sets may differ.",
            stacklevel=3,
        )
```

**Caveat:** `genome_build` may be `None` (when left at default) or a normalised string.
The comparison `a.genome_build != b.genome_build` handles this correctly because
`None != "hg38"` is `True`. Verify the default value in `GlobalConfig.__post_init__`
before finalising — if `genome_build="auto"` is normalised to a concrete value in
`__post_init__`, the comparison will work; if it remains `"auto"`, two objects both
set to `"auto"` will compare equal (safe) while `"auto"` vs `"hg38"` will mismatch
(also correct behaviour).

### 1c. Add to `__all__` in `config.py` if one exists

Add `"ConfigMismatchError"` and `"validate_config_compatibility"` to `__all__`.

---

## Step 2 — Add `config_snapshot` to result dataclasses

The following result classes currently **lack** a `config_snapshot` field.
`ReferencePanelBuildResult` and `RunSummary` already have it and do not need changes.

### 2a. `ChromLDScoreResult` — `src/ldsc/ldscore_calculator.py` lines 62–106

Add field at the end of the dataclass body (before any methods):

```python
config_snapshot: GlobalConfig | None = None
```

Use `None` as default so that existing code that constructs this object without a
snapshot does not break immediately. Downstream validation will check for `None`.

**Caveat:** `frozen=True` means you cannot set this field after construction. Every
call site that constructs `ChromLDScoreResult(...)` must pass `config_snapshot=cfg`.
Find all construction sites with `grep -n "ChromLDScoreResult(" src/`.

### 2b. `LDScoreResult` — `src/ldsc/ldscore_calculator.py` lines 108–169

Same pattern:

```python
config_snapshot: GlobalConfig | None = None
```

Find all construction sites: `grep -n "LDScoreResult(" src/`.

### 2c. `AnnotationBundle` — `src/ldsc/_kernel/annotation.py` lines 131–189

Same pattern:

```python
config_snapshot: GlobalConfig | None = None
```

Find all construction sites: `grep -n "AnnotationBundle(" src/`.

### 2d. `SumstatsTable` — `src/ldsc/sumstats_munger.py` lines 89–156

Same pattern:

```python
config_snapshot: GlobalConfig | None = None
```

Find all construction sites: `grep -n "SumstatsTable(" src/`.

### 2e. `RegressionDataset` — `src/ldsc/regression_runner.py` lines 44–63

Same pattern:

```python
config_snapshot: GlobalConfig | None = None
```

Find all construction sites: `grep -n "RegressionDataset(" src/`.

**General caveat for all Step 2 changes:** Because the field has a default of `None`,
adding it to a `frozen=True` dataclass with existing fields that have no defaults will
cause a `TypeError` at class definition time if any earlier field lacks a default.
Place `config_snapshot` as the **last field** in each dataclass, after all required
(no-default) fields. Verify field ordering by reading each class body before editing.

---

## Step 3 — Propagate snapshots through workflow construction sites

For each dataclass changed in Step 2, find every place the object is constructed and
pass the active `GlobalConfig` as `config_snapshot=`.

### 3a. `ChromLDScoreResult` construction

**Location:** `LDScoreCalculator.compute_chromosome()` in `ldscore_calculator.py`.
This method builds a `ChromLDScoreResult`. It already receives `global_config` as a
parameter. Pass it through:

```python
return ChromLDScoreResult(
    ...,
    config_snapshot=global_config,
)
```

### 3b. `LDScoreResult` construction

**Location:** `LDScoreCalculator.run()` in `ldscore_calculator.py`.
This method aggregates `ChromLDScoreResult` objects. The active `global_config` is
available. Pass it:

```python
return LDScoreResult(
    ...,
    config_snapshot=global_config,
)
```

Also update `LDScoreResult.validate()` to check that `config_snapshot` is not `None`
and that all `chromosome_results` carry the same snapshot (optional but recommended):

```python
def validate(self) -> None:
    ...
    if self.config_snapshot is None:
        raise ValueError("LDScoreResult.config_snapshot must not be None")
```

### 3c. `AnnotationBundle` construction

**Location:** `AnnotationBuilder.run()` in `_kernel/annotation.py`.
`AnnotationBuilder.__init__` already stores `self.global_config`. Pass it when
constructing the bundle:

```python
return AnnotationBundle(
    ...,
    config_snapshot=self.global_config,
)
```

### 3d. `SumstatsTable` construction

**Location:** `SumstatsMunger.run()` in `sumstats_munger.py` and `load_sumstats()`.

For `SumstatsMunger.run()`, `global_config` is passed in or resolved via
`get_global_config()`. Capture it before constructing:

```python
cfg = global_config or get_global_config()
...
return SumstatsTable(
    ...,
    config_snapshot=cfg,
)
```

For `load_sumstats()` (lines 169–207), this function loads a pre-munged `.sumstats`
file and does not go through a workflow config. Use `get_global_config()` as the
snapshot here, but emit a warning that the loaded file's original config is unknown:

```python
import warnings
warnings.warn(
    "load_sumstats() cannot recover the GlobalConfig that was active when this "
    "file was originally munged. Using the current global config as a proxy. "
    "Validate manually if genome_build or snp_identifier matters here.",
    UserWarning,
    stacklevel=2,
)
```

**Caveat:** `load_sumstats()` is an important edge case. Files on disk have no
embedded config provenance. The snapshot attached to the resulting `SumstatsTable`
reflects the caller's current global config, not the original munge config. This is
documented in `docs/config-design.md`. Do not silently attach the wrong config;
the warning above is required.

### 3e. `RegressionDataset` construction

**Location:** `RegressionRunner.build_dataset()` in `regression_runner.py`.
This method receives a `SumstatsTable` and an `LDScoreResult`. By the time
`build_dataset()` is called, Step 4 will have already validated their snapshots.
Capture the resolved config from the `LDScoreResult` snapshot:

```python
config_snapshot = ldscore_result.config_snapshot or self.global_config
return RegressionDataset(
    ...,
    config_snapshot=config_snapshot,
)
```

---

## Step 4 — Add `validate_config_compatibility()` calls at combination points

These are the four places where two independently-computed objects are first merged.

### 4a. `RegressionRunner.build_dataset()` — `regression_runner.py`

This is the highest-value check: it is where a `SumstatsTable` (from munging) meets
an `LDScoreResult` (from LD-score computation). Add at the top of `build_dataset()`,
before any merge logic:

```python
from ldsc.config import validate_config_compatibility
if sumstats_table.config_snapshot is not None and ldscore_result.config_snapshot is not None:
    validate_config_compatibility(
        sumstats_table.config_snapshot,
        ldscore_result.config_snapshot,
        context="SumstatsTable and LDScoreResult",
    )
```

**Note:** The `if` guard handles `None` snapshots gracefully (legacy objects or
objects constructed without a snapshot will skip validation rather than crash).

### 4b. `LDScoreCalculator.run()` — `ldscore_calculator.py`

This is where an `AnnotationBundle` is paired with a `RefPanel` to compute LD scores.
Add at the entry of `run()`, before the chromosome loop:

```python
if annotation_bundle.config_snapshot is not None:
    validate_config_compatibility(
        annotation_bundle.config_snapshot,
        global_config,
        context="AnnotationBundle and LDScoreCalculator runtime config",
    )
```

**Note:** `RefPanel` does not currently carry a full `GlobalConfig` snapshot —
it carries only a `genome_build` string in `RefPanelConfig`. Add a separate check:

```python
if ref_panel.spec.genome_build is not None and global_config.genome_build not in (None, "auto"):
    if ref_panel.spec.genome_build != global_config.genome_build:
        from ldsc.config import ConfigMismatchError
        raise ConfigMismatchError(
            f"genome_build mismatch between RefPanelConfig ({ref_panel.spec.genome_build!r}) "
            f"and active GlobalConfig ({global_config.genome_build!r})."
        )
```

**Caveat:** This check requires that `RefPanel` exposes its `spec` publicly, or that
`RefPanelConfig.genome_build` is accessible. Verify the attribute path in
`_kernel/ref_panel.py` before writing this code.

### 4c. `LDScoreCalculator._aggregate_chromosome_results()` or equivalent

When aggregating multiple `ChromLDScoreResult` objects into one `LDScoreResult`,
verify all chromosome results carry consistent snapshots:

```python
snapshots = [r.config_snapshot for r in chromosome_results if r.config_snapshot is not None]
for s in snapshots[1:]:
    validate_config_compatibility(snapshots[0], s, context="ChromLDScoreResult aggregation")
```

**Location:** Find the aggregation loop in `LDScoreCalculator.run()`.

### 4d. `SumstatsMunger.run()` — `sumstats_munger.py`

No combination of two result objects occurs here (single-input workflow), so no
cross-object check is needed. The config is simply captured into the output
`SumstatsTable` as described in Step 3d.

---

## Step 5 — Export `ConfigMismatchError` from the public API

**File:** `src/ldsc/__init__.py`

Add to the `from .config import ...` import block:

```python
from .config import (
    ...,
    ConfigMismatchError,
    validate_config_compatibility,
)
```

Add both names to `__all__`.

---

## Step 6 — Tests

Add or extend tests in the test suite (wherever existing config tests live —
search for `test_config` or `test_global_config` in the `tests/` directory).

### Test cases to add:

**T1.** `ConfigMismatchError` is raised when `genome_build` differs:
```python
a = GlobalConfig(genome_build="hg19", snp_identifier="chr_pos")
b = GlobalConfig(genome_build="hg38", snp_identifier="chr_pos")
with pytest.raises(ConfigMismatchError, match="genome_build"):
    validate_config_compatibility(a, b, context="test")
```

**T2.** `ConfigMismatchError` is raised when `snp_identifier` differs.

**T3.** No error raised when configs match.

**T4.** `restrict_snps_path` difference emits a `UserWarning`, not an error.

**T5.** `RegressionRunner.build_dataset()` raises `ConfigMismatchError` when
`SumstatsTable` and `LDScoreResult` have mismatched `genome_build` snapshots.

**T6.** `LDScoreCalculator.run()` raises `ConfigMismatchError` when `AnnotationBundle`
snapshot mismatches the active `global_config`.

**T7.** `None` snapshots do not cause crashes (backward-compat guard).

**T8.** `load_sumstats()` emits a `UserWarning` about unknown provenance.

---

## Step 7 — Update `docs/config-design.md` if needed

After implementation, verify that the caveats section in `config-design.md`
accurately reflects any surprises discovered during coding. In particular:

- If `genome_build="auto"` is normalised in `GlobalConfig.__post_init__`, document
  whether it normalises to a concrete value or stays as the string `"auto"`.
- If `load_sumstats()` behaviour was changed from what is described, update the caveat.

---

## File Change Summary

| File | Change type | What changes |
|---|---|---|
| `src/ldsc/config.py` | Add | `ConfigMismatchError` class, `validate_config_compatibility()`, export in `__all__` |
| `src/ldsc/ldscore_calculator.py` | Add field + propagate | `config_snapshot` in `ChromLDScoreResult` and `LDScoreResult`; pass through in `compute_chromosome()` and `run()`; validate in `run()` |
| `src/ldsc/_kernel/annotation.py` | Add field + propagate | `config_snapshot` in `AnnotationBundle`; pass through in `AnnotationBuilder.run()` |
| `src/ldsc/regression_runner.py` | Add field + validate | `config_snapshot` in `RegressionDataset`; validate in `build_dataset()` |
| `src/ldsc/sumstats_munger.py` | Add field + propagate | `config_snapshot` in `SumstatsTable`; pass through in `SumstatsMunger.run()` and `load_sumstats()` |
| `src/ldsc/__init__.py` | Add exports | `ConfigMismatchError`, `validate_config_compatibility` |
| `tests/` (existing test files) | Add tests | T1–T8 as described in Step 6 |

**Files that do NOT need changes:**
- `src/ldsc/outputs.py` — `RunSummary` already has `config_snapshot`
- `src/ldsc/ref_panel_builder.py` — `ReferencePanelBuildResult` already has `config_snapshot`
- `src/ldsc/_kernel/ref_panel.py` — `RefPanelConfig.genome_build` is sufficient for now
- `src/ldsc/config.py` registry functions (`get_global_config`, `set_global_config`,
  `reset_global_config`) — no behavioural changes needed

---

## Known Caveats and Edge Cases

**`RefPanelConfig.genome_build` is a `str`, not a `GlobalConfig`.** The check in Step 4b
compares it directly to `global_config.genome_build`. If `GlobalConfig.__post_init__`
normalises `"hg37"` → `"hg19"` and `"GRCh38"` → `"hg38"`, the comparison will work.
If normalisation is incomplete, the check may produce false mismatches. Verify in
`genome_build_inference.py` and `config.py` before implementing.

**`AnnotationBundle.metadata` genome-build inference.** The `AnnotationBundle` infers
genome build from coordinate metadata, not from the `GlobalConfig`. If a user loads
annotations whose on-disk genome build differs from what `GlobalConfig.genome_build`
says, the snapshot in `AnnotationBundle` will record the config setting, not the
inferred reality. This is a pre-existing issue outside the scope of this plan;
flagging it here so it is not lost.

**`load_sumstats()` provenance is fundamentally unknowable.** Files on disk have no
embedded config metadata. The `UserWarning` in Step 3d is the best possible behaviour.
A future improvement would be to write a sidecar metadata file during `munge_sumstats`
that records the `GlobalConfig` used, and have `load_sumstats()` read it back.
This is out of scope here.

**Backward compatibility of result constructors.** All `config_snapshot` fields are
added with `default=None`, so existing code that constructs result objects without the
field will continue to run. Validation guards on `None` (Step 4) mean old-style
objects skip compatibility checks silently. This is intentional for a gradual rollout
but means the safety guarantee is only fully active once all construction sites pass
a real snapshot (Step 3).

**Thread safety of the global registry.** `_GLOBAL_CONFIG` is a module-level variable
with no locking. This is unchanged by this plan. Parallel notebook cells or
multithreaded workflows that call `set_global_config()` concurrently are still unsafe.
This is a pre-existing limitation.
