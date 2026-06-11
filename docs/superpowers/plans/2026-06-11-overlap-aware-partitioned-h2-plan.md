# Overlap-Aware Partitioned Heritability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `partitioned-h2` emit legacy-compatible overlap-aware category summaries by computing an annotation overlap matrix during `ldscore`, storing it in a parquet sidecar, and feeding it to the already-ported `Hsq._overlap_output`.

**Architecture:** `ldscore` computes the baseline-rows overlap block `A_Bᵀ·A` plus query self-overlaps for the all-SNP and common-SNP universes, sums them across chromosomes, and persists them in `ldscore.overlap.parquet`. `partitioned-h2` loads the sidecar, assembles each fitted model's overlap matrix `O_R` by name, and calls the ported `Hsq._overlap_output` augmented with a one-sided `Coefficient_p`, conditional `Category_h2`, and an `overlap_aware` flag. Two regimes are auto-detected: a joint baseline fit when no query annotations exist (functional enrichment), and the existing baseline+one-query loop when they do (cell-type-specific).

**Tech Stack:** Python 3, NumPy, pandas, pyarrow, scipy.stats, pytest.

**Spec:** `docs/superpowers/specs/2026-06-11-overlap-aware-partitioned-h2-design.md`. Read it before starting.

**Conventions:**
- Run everything in the `restructure` worktree: `/Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured`.
- Activate env for every command: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && <command>`.
- Commit messages: Conventional Commits, subject ≤50 chars; footer `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`. No AI tool names elsewhere.
- TDD: write the failing test, run it red, implement minimally, run it green, commit.

---

## File Structure

**New files:**
- `src/ldsc/_kernel/overlap.py` — `OverlapContribution` dataclass, `compute_overlap()`, `sum_overlap_contributions()`. Pure-NumPy overlap-block math, no I/O.
- `src/ldsc/overlap_matrix.py` — public-layer `LDScoreOverlap` dataclass (labeled in-memory blocks), `overlap_to_long_frame()` / `overlap_from_long_frame()` (parquet (de)serialization), `assemble_model_overlap()` (build `O_R` for one model), `overlap_aware_category_table()` (call ported `_overlap_output` + augment), `model_collinearity_warning()`.
- `tests/test_overlap_matrix.py` — unit tests for the two modules above.

**Modified files:**
- `src/ldsc/_kernel/ldscore.py` — strict-`>` common mask; populate `ChromComputationResult.overlap`.
- `src/ldsc/ldscore_calculator.py` — carry/sum overlap through `_LegacyChromResult` → `ChromLDScoreResult` → `LDScoreResult`; strict-`>` in `count_config`.
- `src/ldsc/outputs.py` — write/read `ldscore.overlap.parquet`; `overlap_config` in metadata; single partitioned-h2 schema for the primary table in both regimes.
- `src/ldsc/regression_runner.py` — load overlap; two regimes; rewrite `summarize_partitioned_h2`; one-sided `Coefficient_p`; `auto` sort; log banner + metadata fields; collinearity warning; column constants.
- `tests/test_ldscore_workflow.py`, `tests/test_regression_workflow.py`, `tests/test_output.py` — extend/replace.
- `docs/current/partitioned-ldsc-workflow.md`, `docs/troubleshooting.md`, `design_map.md` — update.

---

## Milestone 0 — Strict common-MAF filter

### Task 1: Switch the common-SNP mask to strict `>`

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`compute_counts`, ~line 1553-1571)
- Modify: `src/ldsc/ldscore_calculator.py` (`_count_config_from_ldscore_config`, ~line 815-820)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_ldscore_workflow.py`:

```python
def test_compute_counts_common_mask_is_strict_greater_than():
    from ldsc._kernel.ldscore import compute_counts
    metadata = pd.DataFrame({"MAF": [0.05, 0.0500001, 0.2, 0.04]})
    annotations = pd.DataFrame({"base": [1.0, 1.0, 1.0, 1.0], "cat": [1.0, 0.0, 1.0, 1.0]})
    M, M_5_50 = compute_counts(metadata, annotations, common_maf_min=0.05)
    # MAF == 0.05 is excluded (strict >), so only rows 1 and 2 are common.
    assert list(M) == [4.0, 3.0]
    assert list(M_5_50) == [2.0, 1.0]


def test_count_config_records_strict_operator():
    from ldsc.ldscore_calculator import _count_config_from_ldscore_config
    from ldsc.config import LDScoreConfig
    cfg = _count_config_from_ldscore_config(LDScoreConfig(common_maf_min=0.05))
    assert cfg["common_reference_snp_maf_operator"] == ">"
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ldscore_workflow.py -k "strict" -v`
Expected: FAIL — current mask is `>=`, so `M_5_50` is `[3.0, 2.0]` and operator is `">="`.

- [ ] **Step 3: Implement**

In `src/ldsc/_kernel/ldscore.py` `compute_counts`, change the mask and docstring:

```python
    # docstring: "...restricted to rows with ``MAF > common_maf_min``..."
    common = metadata["MAF"] > common_maf_min
```

In `src/ldsc/ldscore_calculator.py` `_count_config_from_ldscore_config`:

```python
    return {
        "common_reference_snp_maf_min": float(ldscore_config.common_maf_min),
        "common_reference_snp_maf_operator": ">",
    }
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_ldscore_workflow.py -k "strict" -v`
Expected: PASS

- [ ] **Step 5: Fix any existing tests asserting `>=`**

Run: `pytest tests/test_ldscore_workflow.py tests/test_output.py -k "common or maf or M_5_50" -v`
Update any test that hard-coded the inclusive boundary or the `">="` operator string to the strict form. Show each change.

- [ ] **Step 6: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(ldscore): use strict MAF > common-maf-min for common counts"
```

---

## Milestone 1 — Overlap computation at `ldscore` time

### Task 2: `OverlapContribution` + `compute_overlap()` kernel function

**Files:**
- Create: `src/ldsc/_kernel/overlap.py`
- Test: `tests/test_overlap_matrix.py`

- [ ] **Step 1: Write the failing test**

Create `tests/test_overlap_matrix.py`:

```python
import numpy as np
import pandas as pd
import pytest


def test_compute_overlap_blocks_binary_and_continuous():
    from ldsc._kernel.overlap import compute_overlap
    # 4 SNPs, baseline = [base(all ones), catA], query = [q, cont]
    annotations = pd.DataFrame({
        "base": [1.0, 1.0, 1.0, 1.0],
        "catA": [1.0, 1.0, 0.0, 0.0],
        "q":    [1.0, 0.0, 1.0, 0.0],
        "cont": [0.5, 2.0, 1.0, 0.0],
    })
    metadata = pd.DataFrame({"MAF": [0.05, 0.06, 0.2, 0.3]})  # row 0 (MAF==0.05) is NOT common
    A = annotations.to_numpy(dtype=np.float64)

    contribution = compute_overlap(metadata, annotations, n_baseline=2, common_maf_min=0.05)

    # baseline-rows block is A_B.T @ A  (2 x 4)
    np.testing.assert_allclose(contribution.baseline_block_all, A[:, :2].T @ A)
    # query self-overlaps: sum of squares for q and cont
    np.testing.assert_allclose(contribution.query_diagonal_all, (A[:, 2:] ** 2).sum(axis=0))
    assert contribution.n_all == 4

    common = np.array([False, True, True, True])
    Ac = A[common]
    np.testing.assert_allclose(contribution.baseline_block_common, Ac[:, :2].T @ Ac)
    np.testing.assert_allclose(contribution.query_diagonal_common, (Ac[:, 2:] ** 2).sum(axis=0))
    assert contribution.n_common == 3


def test_compute_overlap_without_maf_has_no_common_universe():
    from ldsc._kernel.overlap import compute_overlap
    annotations = pd.DataFrame({"base": [1.0, 1.0], "q": [1.0, 0.0]})
    metadata = pd.DataFrame({"CHR": ["1", "1"]})  # no MAF
    contribution = compute_overlap(metadata, annotations, n_baseline=1, common_maf_min=0.05)
    assert contribution.baseline_block_common is None
    assert contribution.query_diagonal_common is None
    assert contribution.n_common is None
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_overlap_matrix.py -k compute_overlap -v`
Expected: FAIL — `ModuleNotFoundError: ldsc._kernel.overlap`.

- [ ] **Step 3: Implement `src/ldsc/_kernel/overlap.py`**

```python
"""Annotation overlap-matrix blocks computed during LD-score estimation.

The overlap matrix ``O = AᵀA`` (``A`` the SNP-by-annotation matrix) drives
legacy overlap-aware partitioned-heritability summaries. Because partitioned-h2
fits ``baseline + one query`` per model, only the baseline-rows block
``A_Bᵀ·A`` (shape ``B × (B+Q)``) and each query's self-overlap ``Σ_s A[s,q]²``
are ever needed; cross-query overlaps are not. This module computes those blocks
for the all-SNP and common-SNP universes, with no I/O.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class OverlapContribution:
    """Per-chromosome (or aggregated) overlap blocks for both SNP universes.

    Attributes
    ----------
    baseline_block_all, baseline_block_common : np.ndarray or None
        ``A_Bᵀ·A`` over the all/common universe, shape ``(B, B+Q)`` with columns
        ordered ``baseline + query``. The common block is ``None`` when MAF is
        unavailable.
    query_diagonal_all, query_diagonal_common : np.ndarray or None
        Query self-overlaps ``Σ_s A[s,q]²``, shape ``(Q,)``.
    n_all, n_common : int or None
        SNP-universe sizes (``M_tot``); ``n_common`` is ``None`` without MAF.
    """
    baseline_block_all: np.ndarray
    baseline_block_common: np.ndarray | None
    query_diagonal_all: np.ndarray
    query_diagonal_common: np.ndarray | None
    n_all: int
    n_common: int | None


def compute_overlap(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    n_baseline: int,
    common_maf_min: float = 0.05,
) -> OverlapContribution:
    """Compute baseline-rows overlap blocks and query self-overlaps.

    ``annotations`` columns are ordered ``baseline + query`` with ``n_baseline``
    leading baseline columns. The common universe is ``MAF > common_maf_min``
    (strict, matching the common-count mask); it is omitted when MAF metadata is
    absent or all-missing.
    """
    A = annotations.to_numpy(dtype=np.float64, copy=False)
    A_baseline = A[:, :n_baseline]
    block_all = A_baseline.T @ A
    query_diag_all = np.einsum("sq,sq->q", A[:, n_baseline:], A[:, n_baseline:])
    n_all = int(A.shape[0])
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        return OverlapContribution(block_all, None, query_diag_all, None, n_all, None)
    common = (metadata["MAF"] > common_maf_min).to_numpy()
    A_common = A[common]
    block_common = A_common[:, :n_baseline].T @ A_common
    query_diag_common = np.einsum(
        "sq,sq->q", A_common[:, n_baseline:], A_common[:, n_baseline:]
    )
    return OverlapContribution(
        block_all, block_common, query_diag_all, query_diag_common, n_all, int(common.sum())
    )
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_overlap_matrix.py -k compute_overlap -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(ldscore): compute baseline-rows annotation overlap blocks"
```

### Task 3: `sum_overlap_contributions()` aggregation

**Files:**
- Modify: `src/ldsc/_kernel/overlap.py`
- Test: `tests/test_overlap_matrix.py`

- [ ] **Step 1: Write the failing test**

```python
def test_sum_overlap_contributions_adds_blocks_and_sizes():
    from ldsc._kernel.overlap import OverlapContribution, sum_overlap_contributions
    a = OverlapContribution(np.ones((1, 2)), np.ones((1, 2)), np.array([2.0]), np.array([1.0]), 5, 3)
    b = OverlapContribution(np.full((1, 2), 3.0), np.full((1, 2), 2.0), np.array([4.0]), np.array([2.0]), 7, 4)
    total = sum_overlap_contributions([a, b])
    np.testing.assert_allclose(total.baseline_block_all, [[4.0, 4.0]])
    np.testing.assert_allclose(total.baseline_block_common, [[3.0, 3.0]])
    np.testing.assert_allclose(total.query_diagonal_all, [6.0])
    np.testing.assert_allclose(total.query_diagonal_common, [3.0])
    assert total.n_all == 12 and total.n_common == 7


def test_sum_overlap_contributions_drops_common_if_any_chrom_lacks_it():
    from ldsc._kernel.overlap import OverlapContribution, sum_overlap_contributions
    a = OverlapContribution(np.ones((1, 1)), np.ones((1, 1)), np.array([1.0]), np.array([1.0]), 3, 2)
    b = OverlapContribution(np.ones((1, 1)), None, np.array([1.0]), None, 2, None)
    total = sum_overlap_contributions([a, b])
    assert total.baseline_block_common is None and total.n_common is None
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_overlap_matrix.py -k sum_overlap -v`
Expected: FAIL — `sum_overlap_contributions` undefined.

- [ ] **Step 3: Implement (append to `src/ldsc/_kernel/overlap.py`)**

```python
def sum_overlap_contributions(contributions: list[OverlapContribution]) -> OverlapContribution:
    """Sum per-chromosome overlap contributions into one aggregate.

    The common universe survives only if every contribution carries it; if any
    chromosome lacks MAF, the aggregate common blocks are ``None`` (mirroring how
    aggregated common counts are dropped when unavailable).
    """
    block_all = np.sum([c.baseline_block_all for c in contributions], axis=0)
    query_diag_all = np.sum([c.query_diagonal_all for c in contributions], axis=0)
    n_all = int(sum(c.n_all for c in contributions))
    has_common = all(c.baseline_block_common is not None for c in contributions)
    if not has_common:
        return OverlapContribution(block_all, None, query_diag_all, None, n_all, None)
    block_common = np.sum([c.baseline_block_common for c in contributions], axis=0)
    query_diag_common = np.sum([c.query_diagonal_common for c in contributions], axis=0)
    n_common = int(sum(c.n_common for c in contributions))
    return OverlapContribution(
        block_all, block_common, query_diag_all, query_diag_common, n_all, n_common
    )
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_overlap_matrix.py -k sum_overlap -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(ldscore): aggregate overlap blocks across chromosomes"
```

### Task 4: Populate `ChromComputationResult.overlap` in both backends

**Files:**
- Modify: `src/ldsc/_kernel/ldscore.py` (`ChromComputationResult` ~line 281; `compute_chrom_from_parquet` ~line 1665; `compute_chrom_from_plink` ~line 1791)

- [ ] **Step 1: Add the field**

In `ChromComputationResult` (`@dataclass`), add after `query_columns`:

```python
    overlap: "OverlapContribution | None" = None
```

Add the import near the top of `ldscore.py`:

```python
from .overlap import OverlapContribution, compute_overlap
```

- [ ] **Step 2: Populate it in `compute_chrom_from_parquet`**

Immediately after the existing `M, M_5_50 = compute_counts(out_metadata, annotations, ...)` call, add:

```python
    overlap = compute_overlap(
        out_metadata,
        annotations,
        n_baseline=len(bundle.baseline_columns),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
```

and pass `overlap=overlap` into the returned `ChromComputationResult(...)`.

- [ ] **Step 3: Populate it in `compute_chrom_from_plink`**

After its `M, M_5_50 = compute_counts(out_metadata, annotation_matrix.reset_index(drop=True), ...)`, add the analogous call (note PLINK uses `annotation_matrix.reset_index(drop=True)` and `out_metadata`):

```python
    overlap = compute_overlap(
        out_metadata,
        annotation_matrix.reset_index(drop=True),
        n_baseline=len(bundle.baseline_columns),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
```

and pass `overlap=overlap` into the returned `ChromComputationResult(...)`.

- [ ] **Step 4: Verify nothing broke**

Run: `pytest tests/test_ldscore_workflow.py tests/test_ldscore_parallelism.py -x -q`
Expected: PASS (field is additive; existing assertions unaffected).

- [ ] **Step 5: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(ldscore): attach overlap blocks to chromosome results"
```

### Task 5: Carry and sum overlap through the calculator onto `LDScoreResult`

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (`_LegacyChromResult` ~93; `ChromLDScoreResult` ~107; `LDScoreResult` ~188; `_wrap_legacy_chrom_result` ~522; `_aggregate_chromosome_results` ~591)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_ldscore_workflow.py` (reuse whatever existing helper builds an in-memory ldscore run with PLINK or parquet fixtures; assert overlap survives aggregation). If an end-to-end helper named `run_small_ldscore(...)` exists, use it; otherwise model this on the nearest existing integration test:

```python
def test_ldscore_result_carries_aggregated_overlap(tmp_path):
    result = run_small_ldscore_with_query(tmp_path)  # existing-style helper; baseline+query
    assert result.overlap is not None
    B = len(result.baseline_columns)
    Q = len(result.query_columns)
    assert result.overlap.baseline_block_all.shape == (B, B + Q)
    assert result.overlap.query_diagonal_all.shape == (Q,)
    assert result.overlap.n_all == len(result.baseline_table)
```

If no such helper exists, build the minimal one in the test using the existing parquet/PLINK fixtures already used by `tests/test_ldscore_workflow.py` (copy the setup from the closest passing test in that file).

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ldscore_workflow.py -k carries_aggregated_overlap -v`
Expected: FAIL — `LDScoreResult` has no `overlap`.

- [ ] **Step 3: Add fields**

Add to `_LegacyChromResult` (frozen dataclass): `overlap: "OverlapContribution | None" = None` and import `from ._kernel.overlap import OverlapContribution, sum_overlap_contributions`.

Add to `ChromLDScoreResult`: `overlap: OverlapContribution | None = field(default=None, repr=False)`.

Add to `LDScoreResult`: `overlap: OverlapContribution | None = field(default=None, repr=False)`.

- [ ] **Step 4: Thread it through `_wrap_legacy_chrom_result`**

When constructing `ChromLDScoreResult(...)`, pass `overlap=getattr(legacy_result, "overlap", None)`.

- [ ] **Step 5: Sum it in `_aggregate_chromosome_results`**

After the `count_totals` aggregation block, add:

```python
        overlaps = [r.overlap for r in chromosome_results if r.overlap is not None]
        aggregated_overlap = (
            sum_overlap_contributions(overlaps) if len(overlaps) == len(chromosome_results) else None
        )
```

and pass `overlap=aggregated_overlap` into the `LDScoreResult(...)` constructor.

- [ ] **Step 6: Run to verify pass**

Run: `pytest tests/test_ldscore_workflow.py -k carries_aggregated_overlap -v`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(ldscore): aggregate overlap blocks onto LDScoreResult"
```

---

## Milestone 2 — Persist and load the overlap sidecar

### Task 6: `LDScoreOverlap` + long-format (de)serialization

**Files:**
- Create: `src/ldsc/overlap_matrix.py`
- Test: `tests/test_overlap_matrix.py`

- [ ] **Step 1: Write the failing test**

```python
def test_overlap_long_frame_round_trip():
    from ldsc._kernel.overlap import OverlapContribution
    from ldsc.overlap_matrix import (
        LDScoreOverlap, overlap_to_long_frame, overlap_from_long_frame,
    )
    contribution = OverlapContribution(
        baseline_block_all=np.array([[10.0, 4.0, 3.0]]),       # base x [base, q, cont]
        baseline_block_common=np.array([[8.0, 3.0, 2.0]]),
        query_diagonal_all=np.array([6.0, 9.0]),                # q, cont
        query_diagonal_common=np.array([5.0, 7.0]),
        n_all=10, n_common=8,
    )
    overlap = LDScoreOverlap.from_contribution(
        contribution, baseline_columns=["base"], query_columns=["q", "cont"]
    )
    frame = overlap_to_long_frame(overlap)
    assert set(frame.columns) == {"row_annotation", "col_annotation", "overlap_all_snps", "overlap_common_snps"}
    restored = overlap_from_long_frame(
        frame, baseline_columns=["base"], query_columns=["q", "cont"],
        total_all_reference_snps=10.0, total_common_reference_snps=8.0,
    )
    pd.testing.assert_frame_equal(restored.baseline_block_all, overlap.baseline_block_all)
    pd.testing.assert_series_equal(restored.query_diagonal_all, overlap.query_diagonal_all)
    assert restored.total_common_reference_snps == 8.0
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_overlap_matrix.py -k long_frame -v`
Expected: FAIL — module/symbols missing.

- [ ] **Step 3: Implement `src/ldsc/overlap_matrix.py` (this step)**

```python
"""Public-layer overlap-matrix container, (de)serialization, and assembly.

The labeled in-memory form keeps the baseline-rows block as a DataFrame indexed
by baseline annotation name with all annotation columns, plus a query
self-overlap Series. Persistence is a long/tidy parquet so the artifact scales
linearly with the query-annotation count and aligns by name.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ._kernel.overlap import OverlapContribution


@dataclass(frozen=True)
class LDScoreOverlap:
    """Labeled overlap blocks for both SNP universes.

    ``baseline_block_*`` is indexed by baseline annotation name with columns over
    all annotations (``baseline + query``). ``query_diagonal_*`` is indexed by
    query annotation name. ``total_*_reference_snps`` are the ``M_tot`` scalars.
    """
    baseline_block_all: pd.DataFrame
    baseline_block_common: pd.DataFrame | None
    query_diagonal_all: pd.Series
    query_diagonal_common: pd.Series | None
    total_all_reference_snps: float
    total_common_reference_snps: float | None

    @classmethod
    def from_contribution(
        cls,
        contribution: OverlapContribution,
        baseline_columns: list[str],
        query_columns: list[str],
    ) -> "LDScoreOverlap":
        all_columns = [*baseline_columns, *query_columns]
        block_all = pd.DataFrame(contribution.baseline_block_all, index=baseline_columns, columns=all_columns)
        qdiag_all = pd.Series(contribution.query_diagonal_all, index=query_columns)
        block_common = qdiag_common = total_common = None
        if contribution.baseline_block_common is not None:
            block_common = pd.DataFrame(contribution.baseline_block_common, index=baseline_columns, columns=all_columns)
            qdiag_common = pd.Series(contribution.query_diagonal_common, index=query_columns)
            total_common = float(contribution.n_common)
        return cls(block_all, block_common, qdiag_all, qdiag_common, float(contribution.n_all), total_common)


def overlap_to_long_frame(overlap: LDScoreOverlap) -> pd.DataFrame:
    """Flatten the baseline-rows block plus query diagonal into long/tidy rows."""
    rows = []
    block_all = overlap.baseline_block_all
    block_common = overlap.baseline_block_common
    for row_name in block_all.index:
        for col_name in block_all.columns:
            rows.append({
                "row_annotation": row_name,
                "col_annotation": col_name,
                "overlap_all_snps": float(block_all.at[row_name, col_name]),
                "overlap_common_snps": (
                    None if block_common is None else float(block_common.at[row_name, col_name])
                ),
            })
    for q in overlap.query_diagonal_all.index:
        rows.append({
            "row_annotation": q,
            "col_annotation": q,
            "overlap_all_snps": float(overlap.query_diagonal_all.at[q]),
            "overlap_common_snps": (
                None if overlap.query_diagonal_common is None else float(overlap.query_diagonal_common.at[q])
            ),
        })
    return pd.DataFrame(rows, columns=["row_annotation", "col_annotation", "overlap_all_snps", "overlap_common_snps"])


def overlap_from_long_frame(
    frame: pd.DataFrame,
    baseline_columns: list[str],
    query_columns: list[str],
    total_all_reference_snps: float,
    total_common_reference_snps: float | None,
) -> LDScoreOverlap:
    """Rebuild labeled blocks from the long/tidy parquet frame."""
    all_columns = [*baseline_columns, *query_columns]
    baseline_set = set(baseline_columns)
    block_rows = frame[frame["row_annotation"].isin(baseline_set)]
    block_all = (
        block_rows.pivot(index="row_annotation", columns="col_annotation", values="overlap_all_snps")
        .reindex(index=baseline_columns, columns=all_columns)
    )
    diag_rows = frame[~frame["row_annotation"].isin(baseline_set)].set_index("row_annotation")
    qdiag_all = diag_rows["overlap_all_snps"].reindex(query_columns)
    has_common = total_common_reference_snps is not None and frame["overlap_common_snps"].notna().any()
    block_common = qdiag_common = None
    if has_common:
        block_common = (
            block_rows.pivot(index="row_annotation", columns="col_annotation", values="overlap_common_snps")
            .reindex(index=baseline_columns, columns=all_columns)
        )
        qdiag_common = diag_rows["overlap_common_snps"].reindex(query_columns)
    return LDScoreOverlap(
        block_all, block_common, qdiag_all, qdiag_common,
        float(total_all_reference_snps),
        None if total_common_reference_snps is None else float(total_common_reference_snps),
    )
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_overlap_matrix.py -k long_frame -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(regression): labeled overlap container and long-frame serde"
```

### Task 7: Write `ldscore.overlap.parquet` + metadata; load it back

**Files:**
- Modify: `src/ldsc/outputs.py` (`LDScoreDirectoryWriter.write` ~267, `build_metadata` ~315)
- Modify: `src/ldsc/regression_runner.py` (`load_ldscore_from_dir` ~2142)
- Test: `tests/test_output.py`, `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing test (writer)**

Add to `tests/test_output.py`, modeled on the existing `LDScoreDirectoryWriter` round-trip test. Construct an `LDScoreResult` whose `overlap` is built via `LDScoreOverlap.from_contribution`, write it, and assert the sidecar + metadata:

```python
def test_ldscore_writer_emits_overlap_sidecar(tmp_path):
    result = make_ldscore_result_with_overlap()  # local helper; baseline=["base"], query=["q"]
    LDScoreDirectoryWriter().write(result, LDScoreOutputConfig(output_dir=str(tmp_path)))
    overlap_path = tmp_path / "ldscore.overlap.parquet"
    assert overlap_path.exists()
    meta = json.loads((tmp_path / "metadata.json").read_text())
    assert meta["files"]["overlap"] == "ldscore.overlap.parquet"
    assert meta["overlap_config"]["common_reference_snp_maf_operator"] == ">"
    assert meta["overlap_config"]["total_all_reference_snps"] == result.overlap.total_all_reference_snps
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_output.py -k overlap_sidecar -v`
Expected: FAIL — no overlap file/metadata key.

- [ ] **Step 3: Implement writer**

In `LDScoreDirectoryWriter.write`, after writing the baseline/query parquet and before `build_metadata`, when `getattr(result, "overlap", None) is not None`:

```python
        overlap = getattr(result, "overlap", None)
        if overlap is not None:
            from .overlap_matrix import overlap_to_long_frame
            overlap_path = output_dir / "ldscore.overlap.parquet"
            overlap_to_long_frame(overlap).to_parquet(overlap_path, index=False)
            paths["overlap"] = overlap_path
```

Add `"overlap"` to the family check list (`paths.values()` already includes it once added before the preflight; place the block before `preflight_output_artifact_family` so the new path is covered, and add `output_dir / "ldscore.overlap.parquet"` to `_ldscore_output_family`).

In `build_metadata`, add the overlap pointer + config:

```python
        overlap = getattr(result, "overlap", None)
        overlap_config = None
        if overlap is not None:
            overlap_config = {
                "total_all_reference_snps": float(overlap.total_all_reference_snps),
                "total_common_reference_snps": (
                    None if overlap.total_common_reference_snps is None
                    else float(overlap.total_common_reference_snps)
                ),
                "common_maf_min": float(dict(getattr(result, "count_config", {})).get("common_reference_snp_maf_min", 0.05)),
                "common_maf_operator": ">",
                "stored_block": "baseline_rows_plus_query_diagonal",
            }
```

and include `"overlap_config": overlap_config` in the returned dict (the `files` dict already receives `overlap` via the writer's `files=` argument because the writer maps every non-`metadata` path name).

- [ ] **Step 4: Run writer test to pass**

Run: `pytest tests/test_output.py -k overlap_sidecar -v`
Expected: PASS

- [ ] **Step 5: Write the failing test (loader)**

In `tests/test_regression_workflow.py`, extend the existing `test_load_ldscore_from_dir_reads_metadata_and_parquet_files` family with a directory that includes `ldscore.overlap.parquet`; assert `result.overlap` is a populated `LDScoreOverlap`.

- [ ] **Step 6: Implement loader**

In `load_ldscore_from_dir`, after reading `count_records` and before constructing `LDScoreResult`:

```python
    overlap = None
    overlap_rel = files.get("overlap")
    if overlap_rel:
        from .overlap_matrix import overlap_from_long_frame
        overlap_frame = pd.read_parquet(root / overlap_rel)
        oc = metadata.get("overlap_config", {})
        overlap = overlap_from_long_frame(
            overlap_frame,
            baseline_columns=baseline_columns,
            query_columns=query_columns,
            total_all_reference_snps=oc.get("total_all_reference_snps"),
            total_common_reference_snps=oc.get("total_common_reference_snps"),
        )
```

and pass `overlap=overlap` into the `LDScoreResult(...)` constructor.

- [ ] **Step 7: Run loader test + full output/regression suites**

Run: `pytest tests/test_output.py tests/test_regression_workflow.py -k "overlap or load_ldscore" -v`
Expected: PASS

- [ ] **Step 8: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(ldscore): persist and load ldscore.overlap.parquet sidecar"
```

---

## Milestone 3 — Overlap-aware summary math

### Task 8: `overlap_aware_category_table()` (ported `_overlap_output` + augmentation)

**Files:**
- Modify: `src/ldsc/overlap_matrix.py`
- Test: `tests/test_overlap_matrix.py`

- [ ] **Step 1: Write the failing test (legacy oracle + augmentation)**

```python
def _legacy_overlap_reference(prop, prop_cov, coef, coef_cov, n_blocks, O, M, M_tot):
    """Transcription of legacy Hsq._overlap_output math (pure NumPy oracle)."""
    from scipy.stats import t as tdist
    K = O.shape[0]
    P = np.zeros((K, K))
    for i in range(K):
        P[i, :] = O[i, :] / M
    prop_h2 = (P @ prop.reshape(-1)).reshape(K)
    prop_h2_se = np.sqrt(np.maximum(0, np.diag(P @ prop_cov @ P.T)))
    prop_snps = (M / M_tot).reshape(K)
    enrich = prop_h2 / prop_snps
    enrich_se = prop_h2_se / prop_snps
    D = np.zeros((K, K))
    for i in range(K):
        if M_tot != M[i]:
            D[i, :] = O[i, :] / M[i] - (M - O[i, :]) / (M_tot - M[i])
    diff_est = D @ coef
    diff_se = np.sqrt(np.diag(D @ coef_cov @ D.T))
    enr_p = [np.nan if diff_se[i] == 0 else 2 * tdist.sf(abs(diff_est[i] / diff_se[i]), n_blocks) for i in range(K)]
    return prop_snps, prop_h2, prop_h2_se, enrich, enrich_se, enr_p


def test_overlap_aware_category_table_matches_legacy_and_augments():
    from types import SimpleNamespace
    from scipy import stats
    from ldsc.overlap_matrix import overlap_aware_category_table

    names = ["base", "catA", "q"]
    K = 3
    M = np.array([10.0, 4.0, 6.0])
    M_tot = 10.0
    # overlap matrix: base overlaps everything fully; catA & q share 2 SNPs.
    O = np.array([[10.0, 4.0, 6.0],
                  [4.0, 4.0, 2.0],
                  [6.0, 2.0, 6.0]])
    rng = np.random.default_rng(0)
    coef = np.array([1e-7, 3e-7, 5e-7])
    coef_cov = np.diag([1e-15, 2e-15, 3e-15])
    prop = np.array([[0.5, 0.2, 0.3]])
    prop_cov = np.diag([1e-3, 2e-3, 3e-3])
    cat = M * coef
    cat_se = np.array([1e-7, 1e-7, 1e-7])
    hsq = SimpleNamespace(
        n_annot=K, prop=prop, prop_cov=prop_cov, coef=coef, coef_cov=coef_cov,
        coef_se=np.sqrt(np.diag(coef_cov)), n_blocks=200, cat=cat, cat_se=cat_se,
    )
    table = overlap_aware_category_table(hsq, O, M, M_tot, names)

    ps, ph, phse, en, ense, enp = _legacy_overlap_reference(
        prop, prop_cov, coef, coef_cov, 200, O, M, M_tot
    )
    np.testing.assert_allclose(table["Prop._SNPs"], ps)
    np.testing.assert_allclose(table["Prop._h2"], ph)
    np.testing.assert_allclose(table["Prop._h2_std_error"], phse)
    np.testing.assert_allclose(table["Enrichment"], en)
    np.testing.assert_allclose(table["Enrichment_std_error"], ense)
    np.testing.assert_allclose(np.asarray(table["Enrichment_p"], float), np.asarray(enp, float))
    # base contains all SNPs -> Enrichment_p is NaN
    assert np.isnan(float(table.loc[0, "Enrichment_p"]))
    # augmentation
    np.testing.assert_allclose(table["Coefficient_p"], stats.norm.sf(coef / np.sqrt(np.diag(coef_cov))))
    np.testing.assert_allclose(table["Category_h2"], cat)
    assert table["overlap_aware"].all()  # off-diagonals > 0
    assert "Coefficient_z" in table.columns and "Coefficient_z-score" not in table.columns
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_overlap_matrix.py -k overlap_aware_category_table -v`
Expected: FAIL — function missing.

- [ ] **Step 3: Implement (append to `src/ldsc/overlap_matrix.py`)**

```python
from scipy import stats

from ._kernel import regression as reg


def overlap_aware_category_table(hsq, overlap_matrix, m_annot, m_tot, category_names) -> pd.DataFrame:
    """Legacy overlap-aware category table plus LDSC3 augmentation columns.

    The overlap-aware columns are produced verbatim by the ported
    ``Hsq._overlap_output``. This function then renames the z-score column,
    appends the one-sided ``Coefficient_p`` (``H1: coefficient > 0``), the
    conditional ``Category_h2`` / ``Category_h2_std_error`` (from ``hsq.cat`` /
    ``hsq.cat_se``), and a per-model ``overlap_aware`` flag.
    """
    overlap_matrix = np.asarray(overlap_matrix, dtype=np.float64)
    table = reg.Hsq._overlap_output(
        hsq,
        list(category_names),
        overlap_matrix,
        np.asarray(m_annot, dtype=np.float64).reshape(1, -1),
        float(m_tot),
        True,
    ).reset_index(drop=True)
    table = table.rename(columns={"Coefficient_z-score": "Coefficient_z"})
    coef = np.ravel(hsq.coef).astype(np.float64)
    coef_se = np.ravel(hsq.coef_se).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(coef_se > 0, coef / coef_se, np.nan)
    table["Coefficient_p"] = stats.norm.sf(z)
    table["Category_h2"] = np.ravel(hsq.cat).astype(np.float64)
    table["Category_h2_std_error"] = np.ravel(hsq.cat_se).astype(np.float64)
    off_diagonal = overlap_matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    table["overlap_aware"] = bool(np.any(off_diagonal > 0))
    return table
```

- [ ] **Step 4: Run to verify pass**

Run: `pytest tests/test_overlap_matrix.py -k overlap_aware_category_table -v`
Expected: PASS

- [ ] **Step 5: Add the hand-computed ground-truth test**

```python
def test_overlap_aware_hand_example_disjoint_plus_base():
    from types import SimpleNamespace
    from ldsc.overlap_matrix import overlap_aware_category_table
    # base = all 6 SNPs; A = first 3; Bcat = last 3 (disjoint A,B but both overlap base)
    names = ["base", "A", "Bcat"]
    M = np.array([6.0, 3.0, 3.0]); M_tot = 6.0
    O = np.array([[6.0, 3.0, 3.0],
                  [3.0, 3.0, 0.0],
                  [3.0, 0.0, 3.0]])
    coef = np.array([0.0, 1.0, 2.0]); cat = M * coef
    hsq = SimpleNamespace(
        n_annot=3, prop=np.array([[0.0, 1/3, 2/3]]), prop_cov=np.eye(3) * 1e-6,
        coef=coef, coef_cov=np.eye(3) * 1e-6, coef_se=np.full(3, 1e-3),
        n_blocks=200, cat=cat, cat_se=np.full(3, 1e-3),
    )
    table = overlap_aware_category_table(hsq, O, M, M_tot, names).set_index("Category")
    # Prop._h2[A] = sum_j O[A,j] coef[j] / h2_tot ; h2_tot = sum M*coef = 3*1+3*2 = 9
    # = (3*0 + 3*1 + 0*2)/9 = 3/9
    assert abs(float(table.loc["A", "Prop._h2"]) - (3/9)) < 1e-9
    assert abs(float(table.loc["Bcat", "Prop._h2"]) - (6/9)) < 1e-9
    assert abs(float(table.loc["base", "Prop._h2"]) - 1.0) < 1e-9
    assert abs(float(table.loc["A", "Prop._SNPs"]) - 0.5) < 1e-12
    assert abs(float(table.loc["A", "Enrichment"]) - (3/9)/0.5) < 1e-9
```

- [ ] **Step 6: Run and commit**

Run: `pytest tests/test_overlap_matrix.py -k "overlap_aware or hand_example" -v`
Expected: PASS

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(regression): overlap-aware category table with augmentation"
```

### Task 9: `assemble_model_overlap()`

**Files:**
- Modify: `src/ldsc/overlap_matrix.py`
- Test: `tests/test_overlap_matrix.py`

- [ ] **Step 1: Write the failing test**

```python
def test_assemble_model_overlap_reconstructs_submatrix():
    from ldsc._kernel.overlap import OverlapContribution
    from ldsc.overlap_matrix import LDScoreOverlap, assemble_model_overlap
    # Full A: baseline=[base, catA], query=[q1, q2]
    A = np.array([
        [1, 1, 1, 0],
        [1, 1, 0, 1],
        [1, 0, 1, 0],
        [1, 0, 0, 1],
    ], dtype=float)
    O_full = A.T @ A
    contribution = OverlapContribution(
        baseline_block_all=A[:, :2].T @ A,
        baseline_block_common=None,
        query_diagonal_all=(A[:, 2:] ** 2).sum(axis=0),
        query_diagonal_common=None,
        n_all=4, n_common=None,
    )
    overlap = LDScoreOverlap.from_contribution(contribution, ["base", "catA"], ["q1", "q2"])
    # model for q1: retained = [base, catA, q1] -> rows/cols 0,1,2 of O_full
    O_R = assemble_model_overlap(overlap, ["base", "catA", "q1"], use_common=False)
    np.testing.assert_allclose(O_R, O_full[np.ix_([0, 1, 2], [0, 1, 2])])
    # symmetry of the assembled query row/col
    np.testing.assert_allclose(O_R[2, :], O_R[:, 2])
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_overlap_matrix.py -k assemble_model_overlap -v`
Expected: FAIL.

- [ ] **Step 3: Implement (append to `src/ldsc/overlap_matrix.py`)**

```python
def assemble_model_overlap(
    overlap: LDScoreOverlap, retained_columns: list[str], use_common: bool
) -> np.ndarray:
    """Build the ``K×K`` model overlap matrix for ``retained_columns`` by name.

    Baseline-baseline and baseline-query entries come from the baseline-rows
    block (using symmetry for query-baseline); the query-query diagonal comes
    from the stored self-overlap. ``use_common`` selects the common-SNP universe.
    """
    block = overlap.baseline_block_common if use_common else overlap.baseline_block_all
    qdiag = overlap.query_diagonal_common if use_common else overlap.query_diagonal_all
    if block is None:
        raise ValueError("common-universe overlap requested but not stored")
    baseline_index = set(block.index)
    k = len(retained_columns)
    out = np.zeros((k, k), dtype=np.float64)
    for a, ca in enumerate(retained_columns):
        for b, cb in enumerate(retained_columns):
            if ca in baseline_index:
                out[a, b] = float(block.at[ca, cb])
            elif cb in baseline_index:
                out[a, b] = float(block.at[cb, ca])
            else:  # both query -> only the diagonal ca == cb is reachable per model
                out[a, b] = float(qdiag.at[ca])
    return out
```

- [ ] **Step 4: Run and commit**

Run: `pytest tests/test_overlap_matrix.py -k assemble_model_overlap -v`
Expected: PASS

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(regression): assemble per-model overlap matrix by name"
```

---

## Milestone 4 — Regimes and wiring into partitioned-h2

### Task 10: Rewrite `summarize_partitioned_h2` to use the overlap matrix

**Files:**
- Modify: `src/ldsc/regression_runner.py` (`summarize_partitioned_h2` ~1213; column constants ~104-139; `_coefficient_p_value` ~1591; delete now-unused `_enrichment_p_value`, `_select_count_key` stays)
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Replace the two column constants with one**

Both regimes use a single schema. Replace `PARTITIONED_H2_AGGREGATE_COLUMNS` and `PARTITIONED_H2_FULL_COLUMNS` with:

```python
PARTITIONED_H2_COLUMNS = [
    "Category", "Prop._SNPs",
    "Category_h2", "Category_h2_std_error",
    "Prop._h2", "Prop._h2_std_error",
    "Enrichment", "Enrichment_std_error", "Enrichment_p",
    "Coefficient", "Coefficient_std_error", "Coefficient_z", "Coefficient_p",
    "overlap_aware",
]
```

Update every reference to the two old names in `regression_runner.py` and `outputs.py` (writer `_select_columns(...)` calls, any `pd.DataFrame(columns=...)` empty-frame construction) to `PARTITIONED_H2_COLUMNS`. Leave `PARTITIONED_H2_SUMMARY_SORT_COLUMNS` unchanged (its columns all still exist).

- [ ] **Step 2: Write the failing test**

In `tests/test_regression_workflow.py`, add a test that fits a baseline+query dataset (reuse `make_ldscore_result` which now must carry `overlap`; build it via `LDScoreOverlap.from_contribution`) and asserts the query row has the single-schema columns and `Prop._SNPs == M_query / M_tot` (not `M_query/ΣM`):

```python
def test_summarize_partitioned_h2_is_overlap_aware(self):
    # dataset with baseline=["base"], one query, M_tot from overlap totals
    summary = summarize_partitioned_h2(hsq, dataset, ["query1"])
    assert list(summary.columns) == PARTITIONED_H2_COLUMNS
    row = summary.iloc[0]
    # Prop._SNPs uses M_tot (universe size), not the sum of M
    assert row["Prop._SNPs"] == pytest.approx(M_query / M_tot)
    assert "overlap_aware" in summary.columns
```

(Provide the full dataset/hsq construction inline, mirroring the existing partitioned tests in this file.)

- [ ] **Step 3: Run to verify failure**

Run: `pytest tests/test_regression_workflow.py -k overlap_aware -v`
Expected: FAIL — current summary uses disjoint `Prop._SNPs` and lacks `overlap_aware`.

- [ ] **Step 4: Implement the rewrite**

Replace the body of `summarize_partitioned_h2` so it:
1. resolves the universe (`use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY`),
2. assembles `O_R = assemble_model_overlap(dataset.ldscore_overlap, dataset.retained_ld_columns, use_common)` (the dataset now carries the loaded `LDScoreOverlap` — see Task 11),
3. takes `m_annot = dataset.reference_snp_count_totals[count_key]` and `m_tot = overlap.total_common_reference_snps if use_common else overlap.total_all_reference_snps`,
4. builds the full table once via `overlap_aware_category_table(hsq, O_R, m_annot, m_tot, dataset.retained_ld_columns)`,
5. selects the requested `annotation_columns` rows in the single `PARTITIONED_H2_COLUMNS` schema.

```python
def summarize_partitioned_h2(hsq, dataset, annotation_columns):
    from .overlap_matrix import assemble_model_overlap, overlap_aware_category_table
    use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY
    overlap = dataset.ldscore_overlap
    O_R = assemble_model_overlap(overlap, list(dataset.retained_ld_columns), use_common)
    m_annot = np.asarray(dataset.reference_snp_count_totals[dataset.count_key_used_for_regression], float)
    m_tot = overlap.total_common_reference_snps if use_common else overlap.total_all_reference_snps
    table = overlap_aware_category_table(hsq, O_R, m_annot, float(m_tot), dataset.retained_ld_columns)
    table = table.set_index("Category")
    rows = table.loc[list(annotation_columns)].reset_index()
    return rows.loc[:, PARTITIONED_H2_COLUMNS].reset_index(drop=True)
```

`estimate_partitioned_h2_batch`'s cell-type loop now calls
`summarize_partitioned_h2(hsq, dataset, [query_column])` for the per-query primary
row and `summarize_partitioned_h2(hsq, dataset, dataset.retained_ld_columns)` for
the optional per-query model table (both yield the single schema; drop the old
`include_full_columns=` argument at every call site).

Delete the now-dead `_enrichment_p_value`, `_safe_vector_divide`, `_retained_snp_counts`, `_coefficient_covariance_matrix`, and the old disjoint loop if no longer referenced (verify with grep before deleting). Keep `_coefficient_p_value` only if referenced elsewhere; the one-sided p now lives in `overlap_aware_category_table`.

- [ ] **Step 5: Add `ldscore_overlap` to the dataset (see Task 11 for the field) — stub here**

If Task 11 is executed after this, temporarily thread `dataset.ldscore_overlap` by reading `ldscore_result.overlap` in `build_dataset`. (Task 11 finalizes this.)

- [ ] **Step 6: Run to verify pass**

Run: `pytest tests/test_regression_workflow.py -k overlap_aware -v`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(partitioned-h2): emit overlap-aware category summaries"
```

### Task 11: Carry overlap into `RegressionDataset`; add the functional regime; fail fast without overlap

**Files:**
- Modify: `src/ldsc/regression_runner.py` (`RegressionDataset` ~154; `build_dataset` ~265; `estimate_partitioned_h2_batch` ~693; `_validate_partitioned_query_columns` ~1128; `run_partitioned_h2_from_args` ~1788)
- Test: `tests/test_regression_workflow.py` (replace `test_estimate_partitioned_h2_requires_query_annotations` at line 1443)

- [ ] **Step 1: Add the dataset field + backcompat error**

Add to `RegressionDataset`: `ldscore_overlap: "LDScoreOverlap | None" = None`. In `build_dataset`, set it from `ldscore_result.overlap`. Add a module constant:

```python
PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE = (
    "partitioned-h2 needs the annotation overlap matrix, but the LD-score directory has no "
    "ldscore.overlap.parquet. Most likely it was produced by an older `ldsc ldscore`. "
    "Regenerate the LD-score directory with the current version. "
    "Other causes & fixes: docs/troubleshooting.md#partitioned-h2-missing-overlap-matrix"
)
```

- [ ] **Step 2: Write the failing tests**

Replace `test_estimate_partitioned_h2_requires_query_annotations` with:

```python
def test_partitioned_h2_functional_regime_runs_joint_baseline_fit(self):
    # baseline-only ldscore result (query_columns == []) now produces one row per baseline category
    result = self.make_functional_ldscore_result()  # baseline=["base","catA"], query=[], with overlap
    summary = RegressionRunner().estimate_partitioned_h2_batch(
        self.make_sumstats_table(), result, SimpleNamespace(query_columns=[]),
    )
    assert set(summary["Category"]) == {"base", "catA"}

def test_partitioned_h2_without_overlap_sidecar_fails_fast(self):
    result = self.make_ldscore_result()  # overlap=None
    object.__setattr__(result, "overlap", None)
    with self.assertRaisesRegex(LDSCInputError, "overlap matrix"):
        RegressionRunner().estimate_partitioned_h2_batch(
            self.make_sumstats_table(), result, SimpleNamespace(query_columns=["query1"]),
        )
```

- [ ] **Step 3: Run to verify failure**

Run: `pytest tests/test_regression_workflow.py -k "functional_regime or without_overlap_sidecar" -v`
Expected: FAIL.

- [ ] **Step 4: Implement the regimes in `estimate_partitioned_h2_batch`**

```python
    if ldscore_result.overlap is None:
        raise LDSCInputError(PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE)
    query_columns = list(getattr(annotation_bundle, "query_columns", []) or [])
    if not query_columns:
        # Functional regime: one joint fit of all baseline annotations.
        dataset = self.build_dataset(sumstats_table, ldscore_result, config=config, query_columns=[])
        hsq = self.estimate_h2(dataset, config=config)
        summary = summarize_partitioned_h2(hsq, dataset, dataset.retained_ld_columns)
        if include_full_partitioned_h2:
            return PartitionedH2BatchResult(summary=summary, per_query_category_tables={}, per_query_metadata={})
        return summary
    # Cell-type regime: existing baseline + one query loop (unchanged below).
    query_columns = _validate_partitioned_query_columns(ldscore_result, query_columns)
    ...
```

Update `_validate_partitioned_query_columns` so it is only called in the cell-type branch (do not reject empty query columns globally). In `run_partitioned_h2_from_args`, build the `query_bundle` from `ldscore_result.query_columns` directly (which may be empty → functional regime) and do not pre-reject.

- [ ] **Step 5: Run to verify pass**

Run: `pytest tests/test_regression_workflow.py -k "functional_regime or without_overlap_sidecar or estimate_partitioned" -v`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(partitioned-h2): add functional regime and overlap fail-fast"
```

### Task 12: Regime-aware `auto` default for `--summary-sort-by`

**Files:**
- Modify: `src/ldsc/regression_runner.py` (`add_partitioned_h2_arguments` ~1699; `_sort_partitioned_h2_summary` ~1658; `run_partitioned_h2_from_args` ~1839)
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
def test_summary_sort_auto_is_regime_aware(self):
    from ldsc.regression_runner import _resolve_summary_sort
    assert _resolve_summary_sort("auto", has_queries=True) == "coefficient-p"
    assert _resolve_summary_sort("auto", has_queries=False) == "category"
    assert _resolve_summary_sort("enrichment", has_queries=True) == "enrichment"
```

- [ ] **Step 2: Run red, implement**

Add to `add_partitioned_h2_arguments`: change `--summary-sort-by` `default="category"` → `default="auto"` and add `"auto"` to `choices`. Add:

```python
def _resolve_summary_sort(sort_by: str, *, has_queries: bool) -> str:
    if sort_by != "auto":
        return sort_by
    return "coefficient-p" if has_queries else "category"
```

In `run_partitioned_h2_from_args`, replace `getattr(args, "summary_sort_by", "category")` with `_resolve_summary_sort(getattr(args, "summary_sort_by", "auto"), has_queries=bool(ldscore_result.query_columns))`.

- [ ] **Step 3: Run green and commit**

Run: `pytest tests/test_regression_workflow.py -k summary_sort_auto -v`
Expected: PASS

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(partitioned-h2): regime-aware auto summary sort default"
```

---

## Milestone 5 — Output presentation

### Task 13: Confirm the single schema flows through the writer in both regimes

With one schema (Task 10), `PartitionedH2DirectoryWriter` already writes `PARTITIONED_H2_COLUMNS` to `partitioned_h2.tsv` for both regimes — no `primary_schema` switch is needed. This task verifies the functional regime (a multi-row summary of baseline categories, no queries) flows through the writer's non-per-query path correctly and creates no per-query tree.

**Files:**
- Test: `tests/test_output.py`
- Modify (only if a test fails): `src/ldsc/outputs.py`

- [ ] **Step 1: Write the test**

```python
def test_partitioned_writer_functional_regime_writes_all_baseline_rows(tmp_path):
    summary = make_schema_summary(categories=["base", "catA", "catB"])  # PARTITIONED_H2_COLUMNS
    PartitionedH2DirectoryWriter().write(summary, PartitionedH2OutputConfig(output_dir=str(tmp_path)))
    written = pd.read_csv(tmp_path / "partitioned_h2.tsv", sep="\t")
    assert list(written.columns) == PARTITIONED_H2_COLUMNS
    assert set(written["Category"]) == {"base", "catA", "catB"}
    assert not (tmp_path / "diagnostics" / "query_annotations").exists()
```

- [ ] **Step 2: Run and fix if needed**

Run: `pytest tests/test_output.py -k functional_regime_writes_all_baseline_rows -v`
Expected: PASS (the existing non-per-query path writes the full summary as-is). If it fails because the writer still references an old column constant, finish the reference swap from Task 10 Step 1.

- [ ] **Step 3: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "test(partitioned-h2): verify single-schema functional output"
```

### Task 14: Log banner + self-describing `diagnostics/metadata.json`

**Files:**
- Modify: `src/ldsc/regression_runner.py` (`run_partitioned_h2_from_args` ~1811-1862)
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing test**

```python
def test_partitioned_h2_metadata_records_regime_and_headline(self, tmp_path):
    # run cell-type regime with --output-dir; read diagnostics/metadata.json
    run_partitioned_h2_from_args(args_with_queries(tmp_path))
    meta = json.loads((tmp_path / "diagnostics" / "metadata.json").read_text())
    assert meta["analysis_type"] == "cell_type_specific"
    assert meta["headline_metric"] == "coefficient"
    assert meta["coefficient_p_test"] == "one_sided_greater"
    assert meta["enrichment_p_test"] == "two_sided_t"
```

- [ ] **Step 2: Run red, implement**

In `run_partitioned_h2_from_args`, after determining the regime, compute:

```python
        has_queries = bool(ldscore_result.query_columns)
        analysis_type = "cell_type_specific" if has_queries else "functional_category"
        headline_metric = "coefficient" if has_queries else "enrichment"
        if has_queries:
            LOGGER.info(
                "Cell-type-specific regime: baseline + one query per model (%d queries). "
                "Focus on `Coefficient` together with the one-sided `Coefficient_p` "
                "(a one-sided test of whether `Coefficient` > 0): a positive `Coefficient` "
                "with a small `Coefficient_p` means the query annotation contributes heritability "
                "beyond the baseline annotations. The `Enrichment` column is confounded by the "
                "query annotation's overlap with the baseline annotations, so use `Coefficient_p` "
                "to judge whether the additional contribution is significant.",
                len(ldscore_result.query_columns),
            )
        else:
            LOGGER.info(
                "Functional-category regime: joint fit of %d baseline categories. "
                "Focus on `Enrichment` and `Enrichment_p`: `Enrichment > 1` means the category's "
                "SNPs explain a larger share of heritability than their share of SNPs (`< 1` means "
                "a smaller share); a small `Enrichment_p` indicates the enrichment is significantly "
                "different from 1, i.e. significantly larger or smaller.",
                len(ldscore_result.baseline_columns),
            )
```

Extend the writer `metadata=` dict with `"analysis_type": analysis_type, "headline_metric": headline_metric, "enrichment_p_test": "two_sided_t", "coefficient_p_test": "one_sided_greater"`.

- [ ] **Step 3: Run green and commit**

Run: `pytest tests/test_regression_workflow.py -k metadata_records_regime -v`
Expected: PASS

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(partitioned-h2): regime banner and self-describing metadata"
```

---

## Milestone 6 — Collinearity warning and docs

### Task 15: Collinearity warning at regression time

**Files:**
- Modify: `src/ldsc/overlap_matrix.py` (`model_collinearity_warning`)
- Modify: `src/ldsc/regression_runner.py` (`estimate_h2` ~614 or `build_dataset`, call site)
- Test: `tests/test_overlap_matrix.py`

- [ ] **Step 1: Write the failing test**

```python
def test_model_collinearity_warning_names_worst_pair():
    from ldsc.overlap_matrix import model_collinearity_warning
    # near-duplicate LD-score columns -> high condition number
    X = np.array([[1.0, 1.0001, 2.0], [2.0, 2.0001, 1.0], [3.0, 3.0001, 0.5], [4.0, 4.0002, 0.1]])
    columns = ["catA", "catA_dup", "other"]
    O = np.array([[10.0, 9.99, 1.0], [9.99, 10.0, 1.0], [1.0, 1.0, 10.0]])
    msg = model_collinearity_warning(X, columns, O, threshold=1e5)
    assert msg is not None and "catA" in msg and "catA_dup" in msg

def test_model_collinearity_warning_none_when_well_conditioned():
    from ldsc.overlap_matrix import model_collinearity_warning
    X = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    assert model_collinearity_warning(X, ["a", "b"], np.eye(2), threshold=1e5) is None
```

- [ ] **Step 2: Run red, implement (append to `overlap_matrix.py`)**

```python
def model_collinearity_warning(X, columns, overlap_matrix, threshold=1e5):
    """Return a warning string when the LD-score design matrix is ill-conditioned.

    Uses the legacy condition-number threshold on the design matrix, and names the
    most correlated annotation pair (from the overlap Gram matrix) as the likely
    culprit. Returns ``None`` when well conditioned.
    """
    X = np.asarray(X, dtype=np.float64)
    if X.shape[1] < 2 or np.linalg.cond(X) <= threshold:
        return None
    O = np.asarray(overlap_matrix, dtype=np.float64)
    d = np.sqrt(np.clip(np.diag(O), 1e-300, None))
    corr = O / np.outer(d, d)
    np.fill_diagonal(corr, 0.0)
    i, j = np.unravel_index(np.argmax(np.abs(corr)), corr.shape)
    return (
        f"LD-score design matrix is nearly collinear (condition number "
        f"{int(np.linalg.cond(X))} > {int(threshold)}); annotations '{columns[i]}' and "
        f"'{columns[j]}' overlap almost entirely (r={corr[i, j]:.3f}). Consider dropping "
        f"one of them or coarsening the annotation."
    )
```

- [ ] **Step 3: Call it from the regression path (deduplicated)**

In `estimate_h2`, after `x = np.asarray(merged[dataset.retained_ld_columns])` and before fitting, when `dataset.ldscore_overlap is not None and len(dataset.retained_ld_columns) > 1`:

```python
        from .overlap_matrix import assemble_model_overlap, model_collinearity_warning
        use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY
        try:
            O_R = assemble_model_overlap(dataset.ldscore_overlap, dataset.retained_ld_columns, use_common)
            warning = model_collinearity_warning(x, dataset.retained_ld_columns, O_R)
        except (ValueError, KeyError):
            warning = None
        if warning and warning not in _SEEN_COLLINEARITY_WARNINGS:
            _SEEN_COLLINEARITY_WARNINGS.add(warning)
            LOGGER.warning(warning)
```

Add a module-level `_SEEN_COLLINEARITY_WARNINGS: set[str] = set()` for per-run de-duplication.

- [ ] **Step 4: Run green and commit**

Run: `pytest tests/test_overlap_matrix.py -k collinearity -v`
Expected: PASS

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "feat(partitioned-h2): warn on collinear annotations at regression time"
```

### Task 16: Documentation + design map

**Files:**
- Modify: `docs/current/partitioned-ldsc-workflow.md`, `docs/troubleshooting.md`, `design_map.md`

- [ ] **Step 1: Update the workflow doc**

In `docs/current/partitioned-ldsc-workflow.md`: replace §6's "partitioned-h2 requires `ldscore.query.parquet`" paragraph with the two regimes; add the LD-score directory listing of `ldscore.overlap.parquet`; document `overlap_config`, the column schemas, the regime-aware sort default, and the headline interpretation (lift wording from the spec §7).

- [ ] **Step 2: Add a troubleshooting section**

In `docs/troubleshooting.md`, add `### partitioned-h2: missing overlap matrix` (slug `partitioned-h2-missing-overlap-matrix`) with the regenerate remedy, matching the in-code message link from Task 11.

- [ ] **Step 3: Update `design_map.md`**

Add rows mapping the spec to `src/ldsc/_kernel/overlap.py`, `src/ldsc/overlap_matrix.py`, and the modified functions.

- [ ] **Step 4: Commit**

```bash
git -C "$REPO" add -A && git -C "$REPO" commit -m "docs(partitioned-h2): document overlap-aware regimes and sidecar"
```

---

## Final verification

- [ ] **Run the full suite**

Run: `source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev && pytest -q`
Expected: all green. Investigate and fix any regression (especially `tests/test_ldscore_workflow.py`, `tests/test_output.py`, `tests/test_regression_workflow.py`).

- [ ] **End-to-end smoke (both regimes)**

Using an existing integration fixture or a tutorial dataset, run `ldsc ldscore` (baseline-only, then baseline+query) and `ldsc partitioned-h2 --output-dir ...`; confirm `ldscore.overlap.parquet` exists, both runs write the single `PARTITIONED_H2_COLUMNS` schema (functional rows = baseline categories, cell-type rows = queries sorted by `Coefficient_p`), and `diagnostics/metadata.json` records the regime.

---

## Self-Review Notes

- **Spec coverage:** §3 math → Tasks 8-9; §4 storage/reuse → Tasks 2-7,9; §5 strict MAF → Task 1; §6 continuous (no binarization) → Task 2 (float64 matmul, tested with a continuous column); §7 schema/files/sort → Tasks 10,12,13; §7.1 signaling → Task 14; §8 collinearity → Task 15; §9 backcompat → Task 11; §10 tests → Tasks 1-9,11; §11 touch points → all.
- **Type consistency:** `OverlapContribution` (kernel) vs `LDScoreOverlap` (public) are distinct by design; `LDScoreOverlap.from_contribution` bridges them. `assemble_model_overlap` consumes `LDScoreOverlap`; `overlap_aware_category_table` consumes a fitted `hsq` and a NumPy `O_R`. `dataset.ldscore_overlap` is the carrier from Task 11 used by Tasks 10 and 15.
- **Ordering caveat:** Task 10 references `dataset.ldscore_overlap` finalized in Task 11; if executed strictly in order, complete Task 11's Step 1 (add the field + set it in `build_dataset`) before running Task 10's test. The plan notes this in Task 10 Step 5.
