# Genome-Build Configuration Refactor: Implementation Plan

All design decisions are finalized in `design.md §6`. No open questions remain.
Read `design.md` in full before starting. Steps must be executed in order;
each step depends on the previous one unless noted otherwise.

Use TDD: write tests first (red), implement (green), verify `pytest` passes
before moving to the next step.

---

## Step 0 — Pre-cleanup (separate commit, no behavior change)

**Files touched:** `tests/test_config_identifiers.py`, any file with stale
`genome_build` TODOs.

- [ ] Delete `test_config_identifiers.py` lines 263 and 267. These tests
  construct `GlobalConfig(ref_panel_snps_file=...)` and
  `GlobalConfig(regression_snps_file=...)`, which are fields that no longer
  exist on `GlobalConfig`. Confirm they fail before deletion to ensure they
  are not testing something real.
- [ ] Grep for `# TODO.*genome_build` and `# FIXME.*genome_build` across
  `src/ldsc/` and remove stale comments. Do not change any logic.
- [ ] Confirm existing `validate_auto_genome_build_mode()` call sites. Keep the
  call in `_kernel/ref_panel.py RefPanel.__init__` for now. The scattered
  kernel-internal calls will be removed in Step 4.

Run `pytest` — all tests must pass before proceeding.

Commit message: `chore: pre-cleanup before genome-build config refactor`

---

## Step 1 — Add `resolve_genome_build()` to `genome_build_inference.py`

**Files touched:**
- `src/ldsc/genome_build_inference.py`
- `src/ldsc/__init__.py`
- `tests/test_genome_build_inference.py`

### 1a — Implementation

Add the following function to `src/ldsc/genome_build_inference.py` after the
existing `validate_auto_genome_build_mode()` function:

```python
def resolve_genome_build(
    hint: str | None,
    snp_identifier: str,
    sample_frame: "pd.DataFrame | None",
    *,
    context: str,
    logger=None,
) -> str | None:
    """Resolve a genome-build hint to a concrete build or None.

    This is the single authoritative entry point for build resolution.
    After this call returns, the value is "hg19", "hg38", or None; it is
    never "auto". Internal kernel components must not call infer_chr_pos_build()
    directly for build resolution.

    Parameters
    ----------
    hint : str or None
        Raw value from CLI or config. Accepted values: "auto", "hg19", "hg38",
        alias strings accepted by normalize_genome_build() (e.g. "hg37",
        "GRCh38"), or None.
    snp_identifier : str
        Already-normalized SNP identifier mode ("rsid" or "chr_pos").
    sample_frame : pd.DataFrame or None
        DataFrame with "CHR" and "POS" columns. Required when hint="auto"
        and snp_identifier="chr_pos"; pass None in all other cases.
    context : str
        Human-readable label for the data source used in log and error
        messages (e.g. "annotation chr1 of baseline.@.annot.gz").
    logger : logging.Logger-like, optional
        Receives inference summary at INFO or WARNING level.

    Returns
    -------
    str or None
        "hg19", "hg38", or None. Returns None when snp_identifier="rsid",
        regardless of hint.

    Raises
    ------
    ValueError
        - hint="auto" and snp_identifier="rsid".
        - hint="auto" and sample_frame is None.
        - hint="auto" and sample_frame has insufficient HM3 overlap.
    """
    snp_identifier = normalize_snp_identifier_mode(snp_identifier)
    hint = normalize_genome_build(hint)

    if snp_identifier == "rsid":
        return None

    # chr_pos path
    if hint != AUTO_GENOME_BUILD:
        return hint  # concrete build: return as-is

    # hint == "auto"
    if snp_identifier == "rsid":
        raise ValueError(
            "genome_build='auto' is not valid for snp_identifier='rsid'."
        )
    if sample_frame is None:
        raise ValueError(
            f"sample_frame is required for genome-build inference of {context!r} "
            "but None was passed. Pass --genome-build hg19 or --genome-build hg38."
        )
    inference = infer_chr_pos_build(
        sample_frame.loc[:, ["CHR", "POS"]],
        context=context,
    )
    msg = f"Inferred genome build for {context}: {inference.summary_message}"
    if logger is not None:
        if inference.coordinate_basis == "0-based":
            logger.warning(msg)
        else:
            logger.info(msg)
    return inference.genome_build
```

Add `"resolve_genome_build"` to `__all__` in `genome_build_inference.py`.

Add to the public re-export block in `src/ldsc/__init__.py`:

```python
from .genome_build_inference import resolve_genome_build
```

### 1b — Tests

Add to `tests/test_genome_build_inference.py` (new test class or extend existing):

```python
class TestResolveGenomeBuild(unittest.TestCase):

    def _make_hg38_frame(self, n=500):
        """Return a synthetic CHR/POS frame in hg38 coordinates."""
        ref = load_packaged_reference_table()
        sample = ref.head(n).copy()
        return pd.DataFrame({"CHR": sample["CHR"], "POS": sample["hg38_POS"]})

    def _make_hg19_frame(self, n=500):
        ref = load_packaged_reference_table()
        sample = ref.head(n).copy()
        return pd.DataFrame({"CHR": sample["CHR"], "POS": sample["hg19_POS"]})

    def _make_tiny_frame(self, n=5):
        """Frame too small to reach MIN_INSPECTED_REFERENCE_SNPS."""
        return pd.DataFrame({"CHR": ["1"] * n, "POS": list(range(n))})

    # --- rsid mode ---
    def test_rsid_auto_hint_returns_none(self):
        self.assertIsNone(resolve_genome_build("auto", "rsid", None, context="test"))

    def test_rsid_hg38_hint_returns_none(self):
        self.assertIsNone(resolve_genome_build("hg38", "rsid", None, context="test"))

    def test_rsid_none_hint_returns_none(self):
        self.assertIsNone(resolve_genome_build(None, "rsid", None, context="test"))

    # --- chr_pos + concrete hint ---
    def test_chr_pos_hg38_hint_returns_hg38(self):
        self.assertEqual(
            resolve_genome_build("hg38", "chr_pos", None, context="test"), "hg38"
        )

    def test_chr_pos_hg19_hint_returns_hg19(self):
        self.assertEqual(
            resolve_genome_build("hg19", "chr_pos", None, context="test"), "hg19"
        )

    def test_chr_pos_alias_hg37_returns_hg19(self):
        self.assertEqual(
            resolve_genome_build("hg37", "chr_pos", None, context="test"), "hg19"
        )

    # --- chr_pos + auto ---
    def test_auto_hg38_sample_infers_hg38(self):
        frame = self._make_hg38_frame()
        result = resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertEqual(result, "hg38")

    def test_auto_hg19_sample_infers_hg19(self):
        frame = self._make_hg19_frame()
        result = resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertEqual(result, "hg19")

    def test_auto_none_sample_frame_raises(self):
        with self.assertRaises(ValueError) as ctx:
            resolve_genome_build("auto", "chr_pos", None, context="test")
        self.assertIn("--genome-build", str(ctx.exception))

    def test_auto_insufficient_overlap_raises(self):
        frame = self._make_tiny_frame()
        with self.assertRaises(ValueError) as ctx:
            resolve_genome_build("auto", "chr_pos", frame, context="test")
        self.assertIn("--genome-build", str(ctx.exception))

    def test_auto_logs_inference_summary(self):
        import logging
        frame = self._make_hg38_frame()
        with self.assertLogs(level="INFO") as log:
            resolve_genome_build(
                "auto", "chr_pos", frame, context="test",
                logger=logging.getLogger("test")
            )
        self.assertTrue(any("Inferred genome build" in m for m in log.output))
```

Run `pytest tests/test_genome_build_inference.py` — all new tests must pass.

Commit message: `feat(genome-build): add resolve_genome_build() to centralize auto-inference`

---

## Step 2 — Tighten `GlobalConfig` validation; update all callsites

**Files touched:**
- `src/ldsc/config.py`
- `src/ldsc/__init__.py`
- `src/ldsc/ldscore_calculator.py`
- `src/ldsc/_kernel/annotation.py`
- `tests/test_config_identifiers.py`
- `tests/test_ref_panel_builder.py`
- `tests/test_annotation.py`
- `tests/test_global_config_registry.py`
- `tutorials/build-parquet-reference-panel-from-plink.md`

### 2a — `config.py` field and `__post_init__`

- [ ] Change line 103: `genome_build: GenomeBuildInput | None = "hg38"`
  → `genome_build: GenomeBuildInput | None = None`

- [ ] Replace the `__post_init__` body with the version in `design.md §3.1`.
  The new `__post_init__` must:
  1. Validate `snp_identifier` is `"rsid"` or `"chr_pos"`.
  2. Call `normalize_genome_build` and write back with `object.__setattr__`.
  3. Call `_normalize_log_level` and write back.
  4. Raise `ValueError` if chr_pos + `genome_build is None`.
  5. Raise `ValueError` if rsid + `genome_build == "auto"`.
  6. Emit `UserWarning` and force `genome_build = None` if rsid + concrete build.

  Add `import warnings` at the top of `config.py` if not already present.

- [ ] Line 115: `_GLOBAL_CONFIG: GlobalConfig = GlobalConfig()`
  → `_GLOBAL_CONFIG: GlobalConfig = GlobalConfig(snp_identifier="rsid")`

- [ ] Line 138 inside `reset_global_config()`: `_GLOBAL_CONFIG = GlobalConfig()`
  → `_GLOBAL_CONFIG = GlobalConfig(snp_identifier="rsid")`

### 2b — Remove CLI fallback constructions

- [ ] `ldscore_calculator.py:814–817` — delete the block:
  ```python
  default_config = GlobalConfig()
  global_config = GlobalConfig(
      snp_identifier=normalized_mode,
      genome_build=default_config.genome_build if getattr(args, "genome_build", None) is None else getattr(args, "genome_build"),
      log_level=getattr(args, "log_level", "INFO"),
  )
  ```
  The replacement `GlobalConfig(...)` construction is added in Step 3b.

- [ ] `_kernel/annotation.py:719` — delete `default_config = GlobalConfig()` and
  any `default_config.genome_build` references in the same block.

- [ ] `ldscore_calculator.py:791` — remove the `validate_auto_genome_build_mode`
  call (now redundant; `GlobalConfig.__post_init__` enforces it).

### 2c — Update `__init__.py` docstring

- [ ] `src/ldsc/__init__.py:26` — update the docstring example:
  `GlobalConfig().snp_identifier` → `GlobalConfig(snp_identifier="rsid").snp_identifier`

### 2d — Update tests

**`tests/test_config_identifiers.py`**

- [ ] Lines 41, 49, 65, 66 — replace bare `GlobalConfig()` with
  `GlobalConfig(snp_identifier="rsid")`.
- [ ] Line 51 — verify `GlobalConfig(snp_identifier=..., ...)` construction
  passes the now-required `genome_build` for chr_pos.
- [ ] Line 72 — `GlobalConfig(genome_build="hg18")` was testing an invalid build
  alias with an implicit chr_pos default. Change to:
  `GlobalConfig(snp_identifier="chr_pos", genome_build="hg18")`.
- [ ] Line 74 — `GlobalConfig(log_level="trace")` had implicit chr_pos default.
  Change to: `GlobalConfig(snp_identifier="rsid", log_level="trace")`.
- [ ] Lines 77–79, 83 — `GlobalConfig(genome_build="hg37")` etc. — add
  `snp_identifier="chr_pos"` to each.
- [ ] Lines 263, 267 — already deleted in Step 0.

**`tests/test_ref_panel_builder.py`** (10 sites)

- [ ] Lines 739, 789, 819, 839, 865, 901, 914, 937, 975, 986 — replace bare
  `GlobalConfig()` with
  `GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")`.

**`tests/test_annotation.py`**

- [ ] Line 127 — `GlobalConfig(snp_identifier="chr_pos")` missing `genome_build`.
  Change to `GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")`.

**`tests/test_global_config_registry.py`**

- [ ] Line 234 — `GlobalConfig(genome_build="hg38", snp_identifier="rsid")` now
  emits `UserWarning` and normalizes `genome_build` to `None`. Update:
  ```python
  with self.assertWarns(UserWarning):
      cfg = GlobalConfig(genome_build="hg38", snp_identifier="rsid")
  self.assertIsNone(cfg.genome_build)
  ```

### 2e — Update tutorial

- [ ] `tutorials/build-parquet-reference-panel-from-plink.md:247` —
  change `GlobalConfig(log_level="INFO")` to
  `GlobalConfig(snp_identifier="chr_pos", genome_build="hg38", log_level="INFO")`.

### 2f — New validation tests

Add to `tests/test_config_identifiers.py` (new test class `TestGlobalConfigValidation`
or extend existing validation class):

```python
def test_chr_pos_no_genome_build_raises(self):
    with self.assertRaises(ValueError) as ctx:
        GlobalConfig(snp_identifier="chr_pos")
    self.assertIn("genome_build is required", str(ctx.exception))

def test_chr_pos_hg38_ok(self):
    cfg = GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
    self.assertEqual(cfg.genome_build, "hg38")

def test_chr_pos_auto_ok(self):
    cfg = GlobalConfig(snp_identifier="chr_pos", genome_build="auto")
    self.assertEqual(cfg.genome_build, "auto")

def test_rsid_auto_raises(self):
    with self.assertRaises(ValueError):
        GlobalConfig(snp_identifier="rsid", genome_build="auto")

def test_rsid_no_genome_build_ok(self):
    cfg = GlobalConfig(snp_identifier="rsid")
    self.assertIsNone(cfg.genome_build)

def test_rsid_with_concrete_build_warns_and_nulls(self):
    with self.assertWarns(UserWarning) as ctx:
        cfg = GlobalConfig(snp_identifier="rsid", genome_build="hg38")
    self.assertIn("ignored", str(ctx.warning))
    self.assertIsNone(cfg.genome_build)

def test_singleton_is_rsid(self):
    from ldsc.config import get_global_config
    cfg = get_global_config()
    self.assertEqual(cfg.snp_identifier, "rsid")
    self.assertIsNone(cfg.genome_build)

def test_reset_returns_rsid(self):
    from ldsc.config import reset_global_config
    cfg = reset_global_config()
    self.assertEqual(cfg.snp_identifier, "rsid")
    self.assertIsNone(cfg.genome_build)
```

Run `pytest` — all tests must pass before proceeding.

Commit message: `feat(config): GlobalConfig.genome_build defaults to None; chr_pos requires explicit build`

---

## Step 3 — Add `_chr_sampler.py` and wire workflow entry points

**Files touched:**
- `src/ldsc/_chr_sampler.py` (new)
- `src/ldsc/ldscore_calculator.py`
- `src/ldsc/_kernel/annotation.py`
- `src/ldsc/sumstats_munger.py`
- `src/ldsc/ref_panel_builder.py`
- `tests/test_chr_sampler.py` (new)
- `tests/test_ldscore_calculator.py` or integration test file

**Depends on:** Steps 1 and 2.

### 3a — Create `src/ldsc/_chr_sampler.py`

Create a new file with the following content. It is private to the `ldsc`
package (not re-exported from `__init__.py`).

```python
"""Private helper for sampling CHR/POS data from per-chromosome file patterns.

Used by workflow entry points to obtain a sample_frame for
resolve_genome_build() when genome_build='auto'.
"""
from __future__ import annotations

import os
from collections.abc import Sequence

import pandas as pd

from .column_inference import infer_chr_pos_columns
from .path_resolution import substitute_chromosome


def sample_frame_from_chr_pattern(
    tokens: Sequence[str],
    chromosomes: Sequence[str] | None,
    *,
    logger=None,
    nrows: int = 5000,
) -> tuple[pd.DataFrame, str]:
    """Load a small CHR/POS sample from the first resolvable chromosome file.

    Parameters
    ----------
    tokens : sequence of str
        Path tokens already split by split_cli_path_tokens(). Each token may
        contain "@" as a chromosome placeholder.
    chromosomes : sequence of str or None
        Chromosome list from args.chromosomes. Used as fallback when chr1
        is not found.
    logger : logging.Logger-like, optional
        Receives an INFO message before reading.
    nrows : int
        Maximum rows to read (default 5000).

    Returns
    -------
    frame : pd.DataFrame
        Columns ["CHR", "POS"] with at most nrows rows.
    resolved_path : str
        The file path actually read.

    Raises
    ------
    ValueError
        No token with "@" found, or no chromosome substitution resolves to
        an existing file.
    """
    pattern_token = next((t for t in tokens if "@" in t), None)
    if pattern_token is None:
        raise ValueError(
            f"No per-chromosome pattern (@) found in tokens {list(tokens)!r}. "
            "Cannot sample for genome-build inference. "
            "Pass --genome-build hg19 or --genome-build hg38 explicitly."
        )

    candidates = ["1"] + [c for c in (chromosomes or []) if c != "1"]
    resolved_path: str | None = None
    used_chrom: str | None = None
    for chrom in candidates:
        path = substitute_chromosome(pattern_token, chrom)
        if os.path.exists(path):
            resolved_path = path
            used_chrom = chrom
            break

    if resolved_path is None:
        tried = ", ".join(candidates)
        raise ValueError(
            f"Could not find a readable chromosome file for pattern "
            f"{pattern_token!r}. Tried chromosomes: {tried}. "
            "Pass --genome-build hg19 or --genome-build hg38 explicitly."
        )

    if logger is not None:
        logger.info(
            "Auto-inferring genome build from first %d rows of %s "
            "(chr%s of %s). Pass --genome-build hg19 or --genome-build hg38 "
            "to skip inference.",
            nrows,
            resolved_path,
            used_chrom,
            pattern_token,
        )

    df = _read_head(resolved_path, nrows)
    chr_col, pos_col = infer_chr_pos_columns(df.columns, context=resolved_path)
    return df.rename(columns={chr_col: "CHR", pos_col: "POS"})[["CHR", "POS"]], resolved_path


def _read_head(path: str, nrows: int) -> pd.DataFrame:
    """Read the first nrows rows of a tabular file (annot.gz or parquet)."""
    if path.endswith(".parquet"):
        import pyarrow.parquet as pq
        table = pq.read_table(path)
        return table.to_pandas().head(nrows)
    return pd.read_csv(path, sep="\t", compression="infer", nrows=nrows)
```

### 3b — Tests for `_chr_sampler.py`

Create `tests/test_chr_sampler.py`:

```python
"""Tests for _chr_sampler.sample_frame_from_chr_pattern."""
import os
import tempfile
import unittest

import pandas as pd

from ldsc._chr_sampler import sample_frame_from_chr_pattern


class TestSampleFrameFromChrPattern(unittest.TestCase):

    def _write_annot(self, directory, chrom, rows=100):
        """Write a minimal .annot.gz file for the given chromosome."""
        path = os.path.join(directory, f"baseline.{chrom}.annot.gz")
        df = pd.DataFrame({
            "CHR": [chrom] * rows,
            "BP": list(range(1, rows + 1)),
            "SNP": [f"rs{i}" for i in range(rows)],
            "CM": [0.0] * rows,
            "base": [1] * rows,
        })
        df.to_csv(path, sep="\t", index=False, compression="gzip")
        return path

    def test_chr1_file_found(self):
        with tempfile.TemporaryDirectory() as d:
            self._write_annot(d, "1")
            token = os.path.join(d, "baseline.@.annot.gz")
            frame, path = sample_frame_from_chr_pattern([token], None)
            self.assertIn("CHR", frame.columns)
            self.assertIn("POS", frame.columns)
            self.assertIn("1", path)

    def test_fallback_to_chromosomes_list(self):
        with tempfile.TemporaryDirectory() as d:
            self._write_annot(d, "22")  # only chr22 exists
            token = os.path.join(d, "baseline.@.annot.gz")
            frame, path = sample_frame_from_chr_pattern(
                [token], chromosomes=["22"]
            )
            self.assertIn("22", path)

    def test_no_at_token_raises(self):
        with self.assertRaises(ValueError) as ctx:
            sample_frame_from_chr_pattern(["/some/literal/path.gz"], None)
        self.assertIn("--genome-build", str(ctx.exception))

    def test_no_existing_file_raises(self):
        with self.assertRaises(ValueError) as ctx:
            sample_frame_from_chr_pattern(
                ["/nonexistent/baseline.@.annot.gz"], ["1", "22"]
            )
        self.assertIn("--genome-build", str(ctx.exception))

    def test_nrows_respected(self):
        with tempfile.TemporaryDirectory() as d:
            self._write_annot(d, "1", rows=200)
            token = os.path.join(d, "baseline.@.annot.gz")
            frame, _ = sample_frame_from_chr_pattern([token], None, nrows=50)
            self.assertLessEqual(len(frame), 50)

    def test_logs_before_reading(self):
        import logging
        with tempfile.TemporaryDirectory() as d:
            self._write_annot(d, "1")
            token = os.path.join(d, "baseline.@.annot.gz")
            logger = logging.getLogger("test_sampler")
            with self.assertLogs(logger, level="INFO") as log:
                sample_frame_from_chr_pattern([token], None, logger=logger)
            self.assertTrue(
                any("Auto-inferring genome build" in m for m in log.output)
            )
```

### 3c — Wire `ldscore_calculator.py`

In `_normalize_run_args()`, replace the removed `default_config` block (Step 2b)
with the logic from `design.md §3.4` (ldscore workflow). The exact insertion
point is after `normalized_args.snp_identifier = normalized_mode` and before
the `return` statement. Import `sample_frame_from_chr_pattern` at the top of the
function or at the module level:

```python
from ._chr_sampler import sample_frame_from_chr_pattern
from .genome_build_inference import resolve_genome_build
```

Key behaviour to implement:

- When `normalized_mode == "chr_pos"` and `hint == "auto"`:
  - Sample annotation (`baseline_annot_sources` tokens) → `annot_build`.
  - Sample reference panel (`r2_sources` tokens) → `ref_build`.
  - If both sources present and `annot_build != ref_build`: raise `ValueError`
    with the mismatch message from `design.md §3.4`.
  - If only one source present: use that build, skip cross-check.
- When `normalized_mode == "chr_pos"` and `hint` is a concrete build: use it
  directly; skip sampling.
- When `normalized_mode == "chr_pos"` and `hint is None`: raise `ValueError`
  ("genome_build is required …").
- When `normalized_mode == "rsid"`: construct
  `GlobalConfig(snp_identifier="rsid", log_level=...)`.

### 3d — Wire `annotation_builder.py`

Same dual-sampling + cross-check as §3c. Annotation tokens from
`--baseline-annot`; reference panel tokens from `--r2-table` (if present).

### 3e — Wire `sumstats_munger.py`

For chr_pos + auto: read first 5000 rows of the sumstats file directly
(no `@` pattern). Use `infer_chr_pos_columns` to find CHR and position columns.
Log before reading. Call `resolve_genome_build()`. For rsid: skip sampling.

### 3f — Wire `ref_panel_builder.py`

`source_genome_build` is always concrete. Call:
```python
resolved_build = resolve_genome_build(
    source_genome_build, "chr_pos", None,
    context="build-ref-panel source", logger=logger,
)
```
`resolve_genome_build` will return the concrete hint unchanged (see behaviour
table in `design.md §3.2`).

### 3g — Integration tests

Add to the relevant test file (or a new `tests/test_workflow_genome_build.py`):

```python
# Test: ldscore auto-inference, matching builds → succeeds
def test_ldscore_auto_infer_matching_builds(self):
    # Fixture: chr1 annotation (hg38) + chr1 ref panel (hg38)
    # Assert: GlobalConfig.genome_build == "hg38" after _normalize_run_args()

# Test: ldscore auto-inference, mismatched builds → raises ValueError
def test_ldscore_auto_infer_build_mismatch_raises(self):
    # Fixture: chr1 annotation (hg19) + chr1 ref panel (hg38)
    # Assert: ValueError containing both build names raised before any kernel call

# Test: ldscore chr1 missing, fallback to chromosomes list
def test_ldscore_auto_infer_fallback_chr(self):
    # Fixture: only chr22 annotation exists; args.chromosomes = ["22"]
    # Assert: chr22 file sampled; correct build inferred

# Test: chr_pos workflow with no --genome-build
def test_chr_pos_no_genome_build_arg_raises(self):
    # For each workflow: omit --genome-build; assert ValueError before computation
```

Run `pytest` — all tests must pass before proceeding.

Commit message: `feat(workflows): call resolve_genome_build() eagerly at each workflow entry`

---

## Step 4 — Remove scattered lazy `"auto"` blocks from kernel

**Files touched:**
- `src/ldsc/_kernel/ref_panel.py`
- `src/ldsc/_kernel/ldscore.py`
- `src/ldsc/_kernel/identifiers.py`
- `src/ldsc/_kernel/annotation.py`

**Depends on:** Step 3. Do not proceed until Step 3's `pytest` run is clean.

- [ ] `_kernel/ref_panel.py:341–347` — delete the block:
  ```python
  if snp_identifier == "chr_pos" and global_config.genome_build == "auto":
      out, _inference = resolve_chr_pos_table(...)
  ```
  Replace with a defensive assertion:
  ```python
  assert global_config.genome_build in {"hg19", "hg38", None}, (
      f"genome_build reached kernel as {global_config.genome_build!r}; "
      "should have been resolved at workflow entry."
  )
  ```

- [ ] `_kernel/ldscore.py:1257–1265` — delete `SortedR2BlockReader.__init__()`'s
  `"auto"` branch. Replace with:
  ```python
  assert genome_build in {"hg19", "hg38", None}, (
      f"genome_build={genome_build!r} must be concrete by this point."
  )
  ```

- [ ] `_kernel/identifiers.py:230` — delete `infer_build = (genome_build == "auto")`
  and the `if infer_build: ... resolve_chr_pos_table(...)` branch in
  `_finalize_chr_pos_restriction_frame`. The function now receives a concrete
  build or `None` and proceeds directly.

- [ ] `_kernel/annotation.py:497` — delete the auto-inference trigger block.

- [ ] Remove `validate_auto_genome_build_mode()` import and calls from all
  kernel-internal sites **except** `_kernel/ref_panel.py RefPanel.__init__`
  (keep that one as a temporary defensive assert; remove when confident).

Run `pytest`. Any failure here means a workflow in Step 3 still passes `"auto"`
to kernel code — investigate and fix before proceeding.

Commit message: `refactor(kernel): remove scattered lazy auto-inference; build resolved at entry`

---

## Step 5 — Remove `RefPanelConfig.genome_build`

**Files touched:**
- `src/ldsc/config.py`
- `src/ldsc/ldscore_calculator.py`
- `src/ldsc/_kernel/ref_panel.py`
- any test constructing `RefPanelConfig(genome_build=...)`

**Depends on:** Step 4.

- [ ] In `config.py` `RefPanelConfig`: delete the field
  `genome_build: GenomeBuildInput | None = None`.
- [ ] In `RefPanelConfig.__post_init__`: delete
  `object.__setattr__(self, "genome_build", normalize_genome_build(self.genome_build))`.
- [ ] In `ldscore_calculator.py:_ref_panel_from_args()`:
  - Line 848: remove `genome_build=global_config.genome_build` from the
    parquet R2 `RefPanelConfig(...)` construction.
  - Line 856: same removal from the plink `RefPanelConfig(...)` construction.
- [ ] In `_kernel/ref_panel.py:269`: replace
  `self.spec.genome_build or self.global_config.genome_build`
  with `self.global_config.genome_build`.
- [ ] Search for `RefPanelConfig(genome_build=` across `tests/` and update each
  to omit the argument.

Run `pytest` — all tests must pass.

Commit message: `feat(config): remove RefPanelConfig.genome_build; build flows through GlobalConfig`

---

## Step 6 — Update CLI help text and documentation

**Files touched:**
- `src/ldsc/ldscore_calculator.py` (parser)
- `src/ldsc/cli.py`
- `docs/class-and-features.md`
- `docs/data_flow.md`
- `docs/genome_build_refactor/design.md`

**Depends on:** Step 2 (can be done in parallel with Steps 3–5).

- [ ] In the ldscore parser and the annotate parser, update the `--genome-build`
  help string to:
  ```
  Genome build for chr_pos inputs. Required when --snp-identifier chr_pos
  (the default). Use 'auto' to infer hg19/hg38 and 0-based/1-based
  coordinates from data. Not used when --snp-identifier rsid.
  ```
- [ ] Confirm the `build-ref-panel` parser has no `"auto"` in the choices for
  `--source-genome-build`. If it does, remove it.
- [ ] Grep `docs/class-and-features.md` and `docs/data_flow.md` for
  `genome_build="hg38"` default references; update to reflect `None`.
- [ ] Update `docs/genome_build_refactor/design.md` status line to
  `Implemented`.

Commit message: `docs: update genome-build CLI help and design docs after refactor`

---

## Step 7 — Final verification

- [ ] Run full test suite: `pytest` — zero failures.
- [ ] Smoke-test with six CLI invocations:
  1. `ldsc ldscore --genome-build hg38 ...` (hg38 chr_pos) → completes.
  2. `ldsc ldscore --genome-build auto ...` (hg38 chr_pos) → infers hg38, completes.
  3. `ldsc ldscore --genome-build hg19 ...` (hg19 chr_pos) → completes.
  4. `ldsc ldscore ...` (chr_pos, no `--genome-build`) → `ValueError` before
     any file reading.
  5. `ldsc ldscore --snp-identifier rsid ...` (no `--genome-build`) → completes;
     `genome_build=null` in manifest.
  6. `ldsc ldscore --snp-identifier rsid --genome-build hg38 ...` → `UserWarning`
     logged; `genome_build=null` in manifest.
- [ ] Open a written manifest (`manifest.json`) from case 2 and confirm
  `config_snapshot.genome_build` is `"hg38"`, not `"auto"`.
- [ ] Open a written manifest from case 5 and confirm
  `config_snapshot.genome_build` is `null`.

---

## Risk Register

| Risk | Mitigation |
|---|---|
| Auto-inference fails on sparse chr_pos input (e.g. annotation file with very few HM3 SNPs) | Q5: fail immediately with error directing user to `--genome-build`; assertions in Step 4 surface any remaining `"auto"` leakage |
| Existing Python-API notebooks call `GlobalConfig()` expecting chr_pos/hg38 | Full callsite inventory in `design.md §6`; each site must be updated explicitly |
| Hidden kernel callsite still receives `"auto"` after Step 4 | Defensive `assert` statements added in Step 4 catch this at first test/smoke run |
| rsid manifests that previously recorded `genome_build="hg38"` lose provenance | After refactor, rsid manifests record `genome_build=null`; old manifests loaded with `config_snapshot=None` per existing compatibility rules in `config-design.md` |
| `substitute_chromosome` raises `FileNotFoundError` for tokens without `@` | `sample_frame_from_chr_pattern` pre-filters for tokens containing `"@"` before calling `substitute_chromosome` |
