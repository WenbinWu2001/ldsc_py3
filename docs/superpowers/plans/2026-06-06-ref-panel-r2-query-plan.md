# Reference-Panel R² Pair Query + R²→r Conversion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a vectorized API + `ldsc query-r2` CLI that, given a list of SNP pairs, returns the stored adjusted R² (plus a query-orientation sign and a per-row status) from a `build-ref-panel` index-format parquet panel, and a standalone adjusted-R²→signed-Pearson-r converter.

**Architecture:** A new public module `src/ldsc/r2_query.py` exposes an `R2Panel` handle (opens a panel once, loads each chromosome lazily, validates the sidecar binding, builds the identity-key→row map), a one-shot `query_r2(...)` wrapper, and a pure `unbiased_r2_to_pearson_r(...)`. A private kernel `src/ldsc/_kernel/r2_query.py` does the point-lookup compute (int64-key match over parquet row groups, with auto-selected random-access pruning vs. streaming). Both reuse the existing parquet-adapter and identity helpers; the LD-score streaming reader is untouched.

**Tech Stack:** Python 3, pandas, numpy, pyarrow (parquet). Tests with pytest. Conventional Commits.

**Design source of truth:** `docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md`.

**Environment for every command:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
```
(Prefix each `pytest`/`git` command with the `source && conda activate` line.)

---

## File Structure

- **Create** `src/ldsc/_kernel/r2_query.py` — point-lookup compute: arrow→numpy helper, `IDX_1` column index, row-group pruning, `lookup_pairs_in_parquet(...)` (int64-key match + strategy switch). No user-facing API, no pandas-dataframe assembly.
- **Create** `src/ldsc/r2_query.py` — public surface: `unbiased_r2_to_pearson_r`, `R2Panel` (open / lazy chrom state / resolve / `query_pairs`), `query_r2`, and the argparse `build_parser` / `run_query_r2_from_args` / `main`.
- **Modify** `src/ldsc/cli.py` — register + dispatch the `query-r2` subcommand.
- **Modify** `src/ldsc/__init__.py` — export `R2Panel`, `query_r2`, `unbiased_r2_to_pearson_r`.
- **Modify** `CLAUDE.md` — add `query-r2` to the documented CLI-surface invariant.
- **Create** `tests/test_r2_query.py` — all tests (converter, kernel, handle, query, CLI) + a shared real-panel fixture builder.
- **Create** `docs/current/ref-panel-r2-query.md` — user-facing reference; **Modify** `design_map.md` — map spec/plan ↔ new modules.

---

## Task 1: Pure converter `unbiased_r2_to_pearson_r`

**Files:**
- Create: `src/ldsc/r2_query.py`
- Test: `tests/test_r2_query.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_r2_query.py
import numpy as np
import pandas as pd
import pytest

from ldsc.r2_query import unbiased_r2_to_pearson_r


class TestUnbiasedR2ToPearsonR:
    def test_inverts_unbiased_correction_then_roots(self):
        # Build adjusted R2 from a known raw R2, then recover signed r.
        n = 1000
        r2_raw = 0.49
        r2_adj = r2_raw - (1.0 - r2_raw) / (n - 2)
        r = unbiased_r2_to_pearson_r(r2_adj, n, sign=-1)
        assert r == pytest.approx(-0.7, abs=1e-3)

    def test_negative_adjusted_value_still_yields_real_r(self):
        # Near-zero LD: adjusted R2 can be slightly negative; r must be real (not NaN).
        r = unbiased_r2_to_pearson_r(-0.0005, 1000, sign=1)
        assert np.isfinite(r)
        assert r >= 0.0

    def test_sign_none_returns_magnitude(self):
        r = unbiased_r2_to_pearson_r(0.25, 1000, sign=None)
        assert r >= 0.0

    def test_vectorized_arrays(self):
        r2 = np.array([0.25, 0.04], dtype=np.float64)
        sign = np.array([1, -1], dtype=np.int8)
        out = unbiased_r2_to_pearson_r(r2, 1000, sign=sign)
        assert out.shape == (2,)
        assert out[0] > 0 and out[1] < 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestUnbiasedR2ToPearsonR -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'ldsc.r2_query'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/ldsc/r2_query.py
"""Query adjusted R² for SNP pairs from index-format reference panels.

Public surface for reading a ``ldsc build-ref-panel`` parquet panel by SNP pair:
the :class:`R2Panel` handle, the one-shot :func:`query_r2` wrapper, and the pure
:func:`unbiased_r2_to_pearson_r` converter. See
``docs/current/ref-panel-r2-query.md`` and the design spec
``docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md``.
"""
from __future__ import annotations

import numpy as np


def unbiased_r2_to_pearson_r(r2_adj, n, sign=None):
    """Convert adjusted (unbiased) R² to a signed Pearson correlation r.

    Inverts the unbiased correction ``r2_adj = r2_raw - (1 - r2_raw) / (n - 2)``
    to recover the biased squared correlation, then takes the square root and
    applies ``sign``. Vectorized: ``r2_adj`` and ``sign`` may be scalars or
    array-likes, and ``n`` may be a scalar or array-like reference sample size.

    Parameters
    ----------
    r2_adj : float or array-like
        Adjusted (unbiased) squared correlation as stored in the panel.
    n : int or array-like
        Reference sample size used to estimate R² (``ldsc:n_samples``).
    sign : int, array-like, or None, optional
        Sign of the Pearson correlation (``+1``/``-1``). ``None`` returns the
        non-negative magnitude ``|r|``.

    Returns
    -------
    float or numpy.ndarray
        Signed Pearson r (or ``|r|`` when ``sign`` is ``None``).
    """
    r2_adj = np.asarray(r2_adj, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    r2_raw = (r2_adj * (n - 2.0) + 1.0) / (n - 1.0)
    r2_raw = np.clip(r2_raw, 0.0, 1.0)
    r = np.sqrt(r2_raw)
    if sign is not None:
        r = r * np.asarray(sign, dtype=np.float64)
    return r if r.ndim else r.item()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestUnbiasedR2ToPearsonR -v`
Expected: PASS (4 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/r2_query.py tests/test_r2_query.py
git commit -m "feat(r2-query): add adjusted-R2 to signed Pearson r converter"
```

---

## Task 2: Kernel point-lookup primitive + shared panel fixture

**Files:**
- Create: `src/ldsc/_kernel/r2_query.py`
- Test: `tests/test_r2_query.py`

This task adds the compute primitive that, given a `ParquetFile` and canonical
`(i, j)` index arrays for one chromosome, returns `(r2, sign)` aligned to the
inputs (`r2=NaN`, `sign=0` where the pair is not stored). It also adds the shared
fixture that builds a **real, binding-valid** panel on disk for all later tasks.

- [ ] **Step 1: Write the failing test (fixture + kernel)**

```python
# tests/test_r2_query.py  (append)
import gzip
from pathlib import Path

import ldsc._kernel.ref_panel_builder as kernel_builder
from ldsc._kernel.snp_identity import sidecar_identity_sha256

pytest.importorskip("pyarrow")


def build_test_panel(
    panel_dir: Path,
    *,
    chrom: str = "1",
    snp_identifier: str = "chr_pos_allele_aware",
    genome_build: str = "hg38",
    n_samples: int = 1000,
):
    """Write one binding-valid index panel (chrN_r2.parquet + chrN_meta.tsv.gz).

    Sidecar rows (0-based index → identity):
        0: 1:100 A/G   1: 1:200 C/T   2: 1:300 A/C   3: 1:400 G/T
    Stored pairs (canonical i<j), R2 (adjusted) and SIGN (Pearson r >= 0):
        (0,1) R2=0.64 sign=+   (0,2) R2=0.04 sign=-   (1,2) R2=0.25 sign=+
    Pair (0,3) is intentionally absent (out-of-window placeholder).
    """
    build_dir = panel_dir / genome_build
    build_dir.mkdir(parents=True, exist_ok=True)
    sidecar = pd.DataFrame(
        {
            "CHR": [chrom] * 4,
            "POS": [100, 200, 300, 400],
            "SNP": ["rs1", "rs2", "rs3", "rs4"],
            "A1": ["A", "C", "A", "G"],
            "A2": ["G", "T", "C", "T"],
            "CM": [0.0, 0.1, 0.2, 0.3],
            "MAF": [0.2, 0.3, 0.25, 0.4],
        }
    )
    meta_path = build_dir / f"chr{chrom}_meta.tsv.gz"
    kernel_builder.write_runtime_metadata_sidecar(
        sidecar, meta_path, genome_build=genome_build, snp_identifier=snp_identifier
    )
    identity_hash = sidecar_identity_sha256(sidecar)
    i = np.array([0, 0, 1], dtype=np.int64)
    j = np.array([1, 2, 2], dtype=np.int64)
    r2 = np.array([0.64, 0.04, 0.25], dtype=np.float32)
    sign = np.array([1, -1, 1], dtype=np.int8)
    r2_path = build_dir / f"chr{chrom}_r2.parquet"
    kernel_builder.write_r2_parquet(
        pair_chunks=[(i, j, r2, sign)],
        path=r2_path,
        genome_build=genome_build,
        n_samples=n_samples,
        snp_identifier=snp_identifier,
        n_snps=len(sidecar),
        sidecar_identity_sha256=identity_hash,
        min_r2=0.0,
    )
    return build_dir, meta_path, r2_path


class TestLookupPairsInParquet:
    def _open(self, r2_path):
        import pyarrow.parquet as pq

        return pq.ParquetFile(str(r2_path))

    def test_returns_stored_r2_and_sign(self, tmp_path):
        from ldsc._kernel.r2_query import lookup_pairs_in_parquet

        _, _, r2_path = build_test_panel(tmp_path)
        pf = self._open(r2_path)
        i = np.array([0, 0, 1], dtype=np.int64)
        j = np.array([1, 2, 2], dtype=np.int64)
        r2, sign = lookup_pairs_in_parquet(
            pf, i, j, n_snps=4, r2_scale=32767.0, strategy="stream"
        )
        np.testing.assert_allclose(r2, [0.64, 0.04, 0.25], atol=1.5e-5)
        np.testing.assert_array_equal(sign, [1, -1, 1])

    def test_absent_pair_is_nan_with_zero_sign(self, tmp_path):
        from ldsc._kernel.r2_query import lookup_pairs_in_parquet

        _, _, r2_path = build_test_panel(tmp_path)
        pf = self._open(r2_path)
        i = np.array([0], dtype=np.int64)
        j = np.array([3], dtype=np.int64)  # (0,3) not stored
        r2, sign = lookup_pairs_in_parquet(pf, i, j, n_snps=4, r2_scale=32767.0)
        assert np.isnan(r2[0])
        assert sign[0] == 0

    def test_random_and_stream_strategies_match(self, tmp_path):
        from ldsc._kernel.r2_query import lookup_pairs_in_parquet

        _, _, r2_path = build_test_panel(tmp_path)
        pf = self._open(r2_path)
        i = np.array([0, 0, 1, 0], dtype=np.int64)
        j = np.array([1, 2, 2, 3], dtype=np.int64)
        rg_r2, rg_sign = lookup_pairs_in_parquet(
            pf, i, j, n_snps=4, r2_scale=32767.0, strategy="random"
        )
        st_r2, st_sign = lookup_pairs_in_parquet(
            pf, i, j, n_snps=4, r2_scale=32767.0, strategy="stream"
        )
        np.testing.assert_array_equal(np.nan_to_num(rg_r2, nan=-1), np.nan_to_num(st_r2, nan=-1))
        np.testing.assert_array_equal(rg_sign, st_sign)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestLookupPairsInParquet -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'ldsc._kernel.r2_query'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/ldsc/_kernel/r2_query.py
"""Point-lookup compute for index-format reference-panel R² parquet files.

Given a per-chromosome ``pyarrow.parquet.ParquetFile`` and canonical ``(i, j)``
sidecar-row index arrays (``i < j``), return the stored adjusted R² and the
panel-orientation sign for each pair. Pairs are matched by an int64 key
``i * n_snps + j`` against decoded row groups; absent pairs yield ``NaN`` R² and
sign ``0``. Two strategies share the match primitive: random-access (row-group
pruning on ``IDX_1`` footer statistics) and streaming (scan every row group).

This module is private kernel code. The user-facing surface is ``ldsc.r2_query``.
"""
from __future__ import annotations

import numpy as np


def _arrow_to_numpy(column):
    """Convert an Arrow column to NumPy across PyArrow versions."""
    try:
        return column.to_numpy(zero_copy_only=False)
    except TypeError:
        return column.to_numpy()


def _idx1_column_index(parquet_file) -> int:
    """Return the parquet column index of ``IDX_1`` (for footer statistics)."""
    return parquet_file.schema_arrow.names.index("IDX_1")


def _prune_row_groups(parquet_file, left_sorted: np.ndarray, idx1_col: int):
    """Return row groups whose ``IDX_1`` [min,max] covers a needed left index.

    ``left_sorted`` is the sorted unique set of canonical left indices queried.
    A row group qualifies when at least one queried left index falls within its
    ``IDX_1`` footer range. Groups without usable statistics are kept (safe).
    """
    md = parquet_file.metadata
    keep: list[int] = []
    for rg in range(md.num_row_groups):
        stats = md.row_group(rg).column(idx1_col).statistics
        if stats is None or not stats.has_min_max:
            keep.append(rg)
            continue
        lo, hi = int(stats.min), int(stats.max)
        pos = int(np.searchsorted(left_sorted, lo))
        if pos < left_sorted.shape[0] and int(left_sorted[pos]) <= hi:
            keep.append(rg)
    return keep


def lookup_pairs_in_parquet(
    parquet_file,
    i: np.ndarray,
    j: np.ndarray,
    *,
    n_snps: int,
    r2_scale: float | None,
    strategy: str = "auto",
    strategy_threshold: int = 50_000,
):
    """Return ``(r2, sign)`` for canonical pairs ``(i, j)`` from one parquet.

    ``r2`` is float32 with ``NaN`` where the pair is not stored; ``sign`` is int8
    with ``+1``/``-1`` for stored pairs (panel orientation) and ``0`` where the
    pair is not stored. When ``r2_scale`` is set the int16 R² column is divided
    by it (dequantization); a float R² column passes through unscaled.
    """
    i = np.asarray(i, dtype=np.int64)
    j = np.asarray(j, dtype=np.int64)
    n = i.shape[0]
    r2_out = np.full(n, np.nan, dtype=np.float32)
    sign_out = np.zeros(n, dtype=np.int8)
    if n == 0:
        return r2_out, sign_out

    n_snps = np.int64(n_snps)
    qkey = i * n_snps + j
    uniq_key = np.unique(qkey)
    uniq_r2 = np.full(uniq_key.shape[0], np.nan, dtype=np.float32)
    uniq_sign = np.zeros(uniq_key.shape[0], dtype=np.int8)

    chosen = strategy
    if strategy == "auto":
        chosen = "random" if n <= strategy_threshold else "stream"

    idx1_col = _idx1_column_index(parquet_file)
    if chosen == "random":
        row_groups = _prune_row_groups(parquet_file, np.unique(i), idx1_col)
    else:
        row_groups = range(parquet_file.metadata.num_row_groups)

    for rg in row_groups:
        table = parquet_file.read_row_group(int(rg), columns=["IDX_1", "IDX_2", "R2", "SIGN"])
        gi = _arrow_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        gj = _arrow_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        gkey = gi * n_snps + gj
        pos = np.searchsorted(uniq_key, gkey)
        in_range = pos < uniq_key.shape[0]
        matched = np.zeros(gkey.shape[0], dtype=bool)
        matched[in_range] = uniq_key[pos[in_range]] == gkey[in_range]
        if not matched.any():
            continue
        mpos = pos[matched]
        r2_vals = _arrow_to_numpy(table.column("R2")).astype(np.float32, copy=False)[matched]
        if r2_scale is not None:
            r2_vals = r2_vals / np.float32(r2_scale)
        sign_vals = _arrow_to_numpy(table.column("SIGN")).astype(bool, copy=False)[matched]
        uniq_r2[mpos] = r2_vals
        uniq_sign[mpos] = np.where(sign_vals, np.int8(1), np.int8(-1))

    qpos = np.searchsorted(uniq_key, qkey)
    r2_out[:] = uniq_r2[qpos]
    sign_out[:] = uniq_sign[qpos]
    return r2_out, sign_out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestLookupPairsInParquet -v`
Expected: PASS (3 tests). The fixture writing a real binding-valid panel also exercises `write_r2_parquet` + `write_runtime_metadata_sidecar`.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/r2_query.py tests/test_r2_query.py
git commit -m "feat(r2-query): add parquet point-lookup kernel with strategy switch"
```

---

## Task 3: `R2Panel.open` + lazy chromosome state + binding validation

**Files:**
- Modify: `src/ldsc/r2_query.py`
- Test: `tests/test_r2_query.py`

Adds the handle that resolves panel paths, loads each chromosome on first touch
(sidecar, schema metadata, binding check, identity-key→row map, parquet handle +
dequant scale), and exposes `chromosomes` / `snp_identifier` / `n_samples`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_r2_query.py  (append)
class TestR2PanelOpen:
    def test_open_panel_dir_exposes_metadata(self, tmp_path):
        from ldsc.r2_query import R2Panel

        build_test_panel(tmp_path, snp_identifier="chr_pos_allele_aware", n_samples=1234)
        panel = R2Panel.open(tmp_path, genome_build="hg38")
        assert panel.chromosomes == ["1"]
        assert panel.snp_identifier == "chr_pos_allele_aware"
        assert panel.n_samples == 1234

    def test_open_explicit_meta_and_parquet(self, tmp_path):
        from ldsc.r2_query import R2Panel

        _, meta_path, r2_path = build_test_panel(tmp_path)
        panel = R2Panel.open(meta_path=meta_path, parquet_path=r2_path)
        assert panel.chromosomes == ["1"]

    def test_open_requires_exactly_one_input_mode(self, tmp_path):
        from ldsc.errors import LDSCUsageError
        from ldsc.r2_query import R2Panel

        with pytest.raises(LDSCUsageError):
            R2Panel.open()  # neither
        _, meta_path, r2_path = build_test_panel(tmp_path)
        with pytest.raises(LDSCUsageError):
            R2Panel.open(tmp_path, meta_path=meta_path, parquet_path=r2_path)  # both

    def test_binding_mismatch_is_hard_error(self, tmp_path):
        from ldsc.errors import LDSCInputError
        from ldsc.r2_query import R2Panel

        build_dir, meta_path, _ = build_test_panel(tmp_path)
        # Corrupt the sidecar identity by rewriting one allele (breaks the hash).
        with gzip.open(meta_path, "rt") as handle:
            lines = handle.readlines()
        body = [ln for ln in lines if not ln.startswith("#")]
        header = [ln for ln in lines if ln.startswith("#")]
        body[1] = body[1].replace("\tA\t", "\tT\t", 1)  # first data row allele edit
        with gzip.open(meta_path, "wt") as handle:
            handle.writelines(header + body)
        panel = R2Panel.open(tmp_path, genome_build="hg38")
        with pytest.raises(LDSCInputError):
            panel.query_pairs(pd.DataFrame({
                "CHR_1": [1], "POS_1": [100], "A1_1": ["A"], "A2_1": ["G"],
                "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"],
            }))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestR2PanelOpen -v`
Expected: FAIL with `AttributeError`/`ImportError` for `R2Panel`.

- [ ] **Step 3: Write minimal implementation**

Add these imports at the top of `src/ldsc/r2_query.py` (below `import numpy as np`):

```python
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ._kernel.identifiers import build_snp_id_series
from ._kernel.r2_query import lookup_pairs_in_parquet
from ._kernel.ref_panel import (
    _r2_dir_metadata_paths,
    _r2_dir_r2_paths,
    _read_r2_schema_meta,
)
from ._kernel.ldscore import _load_full_panel_sidecar, _validate_index_binding
from ._kernel.snp_identity import (
    COMPLEMENT,
    identity_mode_family,
    is_allele_aware_mode,
    normalize_snp_identifier_mode,
)
from .chromosome_inference import normalize_chromosome
from .errors import LDSCInputError, LDSCUsageError
from ._logging import get_logger

LOGGER = get_logger(__name__)
```

Then add the chromosome-state container and the handle:

```python
@dataclass
class _ChromState:
    """Cached per-chromosome state for one open panel chromosome."""

    sidecar: pd.DataFrame
    key_index: pd.Index          # identity key -> row position (unique)
    a1: np.ndarray               # panel A1 (object), aligned to sidecar rows
    a2: np.ndarray               # panel A2 (object), aligned to sidecar rows
    parquet_file: object         # pyarrow.parquet.ParquetFile
    n_snps: int
    r2_scale: float | None


class R2Panel:
    """Handle over a build-ref-panel index-format R² panel for pair queries.

    Open once with :meth:`open`, then call :meth:`query_pairs` many times. Each
    chromosome is loaded lazily on first touch (sidecar, schema metadata,
    sidecar-binding check, identity-key map, parquet handle), then cached.
    """

    def __init__(self, paths: dict[str, tuple[str, str]], snp_identifier: str) -> None:
        """Store the resolved ``{chrom: (meta, parquet)}`` map and active mode."""
        self._paths = paths
        self._mode = snp_identifier
        self._state: dict[str, _ChromState] = {}
        self._n_samples: int | None = None
        self._n_samples_loaded = False

    @classmethod
    def open(
        cls,
        panel_dir=None,
        *,
        meta_path=None,
        parquet_path=None,
        snp_identifier: str | None = None,
        genome_build: str | None = None,
        validate_binding: bool = True,
    ) -> "R2Panel":
        """Open a panel from a directory or an explicit meta/parquet pair."""
        import pyarrow.parquet as pq

        dir_given = panel_dir is not None
        explicit_given = meta_path is not None or parquet_path is not None
        if dir_given == explicit_given:
            raise LDSCUsageError(
                "R2Panel.open requires exactly one input mode: either `panel_dir`, or "
                "both `meta_path` and `parquet_path`. Most likely neither or both were "
                "passed. Pass a build-ref-panel directory, or a single chromosome's "
                "meta+parquet pair."
            )

        paths: dict[str, tuple[str, str]] = {}
        if dir_given:
            r2_paths = _r2_dir_r2_paths(panel_dir, genome_build=genome_build, chrom=None)
            for r2 in r2_paths:
                chrom = _chrom_from_r2_path(r2)
                meta = _r2_dir_metadata_paths(panel_dir, genome_build=genome_build, chrom=chrom)
                if not meta:
                    raise LDSCInputError(
                        f"R2Panel could not find the metadata sidecar for chromosome {chrom} "
                        f"next to '{r2}'. Most likely the `chr{chrom}_meta.tsv.gz` sidecar was "
                        "not copied with the parquet. Restore the sidecar or regenerate the panel."
                    )
                paths[chrom] = (meta[0], r2)
        else:
            if meta_path is None or parquet_path is None:
                raise LDSCUsageError(
                    "R2Panel.open with explicit paths requires both `meta_path` and "
                    "`parquet_path`. Most likely only one was given. Pass both."
                )
            chrom = _chrom_from_r2_path(str(parquet_path))
            paths[chrom] = (str(meta_path), str(parquet_path))

        if not paths:
            raise LDSCInputError(
                "R2Panel.open found no `chr*_r2.parquet` files. Most likely `panel_dir` is "
                "not a build-ref-panel output directory or the wrong genome-build sub-dir "
                "was selected. Pass the directory containing canonical R2 parquet files."
            )

        if snp_identifier is None:
            first_r2 = next(iter(paths.values()))[1]
            schema_meta = pq.ParquetFile(first_r2).schema_arrow.metadata or {}
            raw = schema_meta.get(b"ldsc:snp_identifier")
            if raw is None:
                raise LDSCInputError(
                    f"R2Panel could not read `ldsc:snp_identifier` from '{first_r2}'. Most "
                    "likely the panel predates identity provenance. Pass `snp_identifier=` "
                    "explicitly (rsid, rsid_allele_aware, chr_pos, or chr_pos_allele_aware)."
                )
            snp_identifier = raw.decode("utf-8")

        panel = cls(paths, normalize_snp_identifier_mode(snp_identifier))
        panel._validate_binding = validate_binding
        return panel

    @property
    def chromosomes(self) -> list[str]:
        """Chromosomes available in this panel (in resolved order)."""
        return list(self._paths.keys())

    @property
    def snp_identifier(self) -> str:
        """Active SNP identifier mode used to resolve query SNPs."""
        return self._mode

    @property
    def n_samples(self) -> int | None:
        """Reference sample size from parquet metadata, or ``None`` if absent."""
        if not self._n_samples_loaded:
            first_r2 = next(iter(self._paths.values()))[1]
            self._n_samples = _read_r2_schema_meta(first_r2).n_samples
            self._n_samples_loaded = True
        return self._n_samples

    def _chrom_state(self, chrom: str) -> _ChromState:
        """Load and cache per-chromosome state on first touch."""
        chrom = normalize_chromosome(chrom)
        if chrom in self._state:
            return self._state[chrom]
        if chrom not in self._paths:
            raise LDSCInputError(
                f"R2Panel has no chromosome {chrom}. Available: {self.chromosomes}. Most "
                "likely the query references a chromosome not built into this panel."
            )
        meta_path, r2_path = self._paths[chrom]
        import pyarrow as pa
        import pyarrow.parquet as pq

        sidecar = _load_full_panel_sidecar(r2_path)
        pf = pq.ParquetFile(r2_path)
        schema_meta = pf.schema_arrow.metadata or {}
        n_snps = int(schema_meta[b"ldsc:n_snps"].decode("utf-8"))
        if getattr(self, "_validate_binding", True):
            identity_hash = schema_meta[b"ldsc:sidecar_identity_sha256"].decode("utf-8")
            _validate_index_binding(
                sidecar, n_snps=n_snps, identity_hash=identity_hash, context=f"R2Panel[{chrom}] {r2_path}"
            )
        key = build_snp_id_series(sidecar, self._mode)
        key_index = pd.Index(key)
        if key_index.has_duplicates:
            raise LDSCInputError(
                f"R2Panel cannot key chromosome {chrom} in mode '{self._mode}': the sidecar "
                "has duplicate SNP identity keys, so a query SNP could resolve ambiguously. "
                "Most likely the mode is too coarse for the panel (e.g. duplicate positions "
                "in chr_pos). Use an allele-aware mode."
            )
        r2_scale = None
        if pa.types.is_integer(pf.schema_arrow.field("R2").type):
            scale_raw = schema_meta.get(b"ldsc:r2_scale")
            r2_scale = float(scale_raw.decode("utf-8")) if scale_raw is not None else 32767.0
        state = _ChromState(
            sidecar=sidecar,
            key_index=key_index,
            a1=sidecar["A1"].to_numpy(dtype=object),
            a2=sidecar["A2"].to_numpy(dtype=object),
            parquet_file=pf,
            n_snps=n_snps,
            r2_scale=r2_scale,
        )
        self._state[chrom] = state
        return state
```

Add the small path helper at module end:

```python
def _chrom_from_r2_path(r2_path: str) -> str:
    """Extract the chromosome label from a ``chr{N}_r2.parquet`` path."""
    name = Path(r2_path).name
    if not name.startswith("chr") or not name.endswith("_r2.parquet"):
        raise LDSCInputError(
            f"R2Panel expected a canonical `chr{{N}}_r2.parquet` filename but got '{name}'. "
            "Most likely a non-panel parquet was passed. Pass a build-ref-panel R2 file."
        )
    return normalize_chromosome(name[len("chr") : -len("_r2.parquet")])
```

(Defining `query_pairs` is Task 4; the binding-mismatch test reaches it but only
exercises `_chrom_state`'s validation, which raises before any lookup.)

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestR2PanelOpen -v`
Expected: PASS (4 tests). `test_binding_mismatch_is_hard_error` requires
`query_pairs` to exist; if running this task in isolation, implement Task 4 first
or temporarily assert via `panel._chrom_state("1")`. **When executing the full
plan in order, implement Task 4 before running this test.**

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/r2_query.py tests/test_r2_query.py
git commit -m "feat(r2-query): add R2Panel handle with lazy chromosome loading"
```

---

## Task 4: `R2Panel.query_pairs` — resolution, sign, status, output

**Files:**
- Modify: `src/ldsc/r2_query.py`
- Test: `tests/test_r2_query.py`

Implements the end-to-end pair query: build per-endpoint keys, resolve to
`(chrom, row, panel A1/A2)`, classify each pair (diagonal / stored / absent /
cross-chromosome / not-in-panel), run the kernel lookup per chromosome, harmonize
the sign in allele-aware modes (always `NA` in base modes), and assemble the
output DataFrame with the always-on `status` column.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_r2_query.py  (append)
class TestQueryPairs:
    def _panel(self, tmp_path, mode="chr_pos_allele_aware"):
        from ldsc.r2_query import R2Panel

        build_test_panel(tmp_path, snp_identifier=mode)
        return R2Panel.open(tmp_path, genome_build="hg38")

    def test_stored_diagonal_absent_cross_and_missing(self, tmp_path):
        panel = self._panel(tmp_path)
        pairs = pd.DataFrame(
            {
                # stored (0,1); diagonal (0,0); absent (0,3); not-in-panel; cross-chr
                "CHR_1": [1, 1, 1, 1, 1],
                "POS_1": [100, 100, 100, 100, 100],
                "A1_1": ["A", "A", "A", "A", "A"],
                "A2_1": ["G", "G", "G", "G", "G"],
                "CHR_2": [1, 1, 1, 1, 2],
                "POS_2": [200, 100, 400, 999, 200],
                "A1_2": ["C", "A", "G", "X", "C"],
                "A2_2": ["T", "G", "T", "Y", "T"],
            }
        )
        out = panel.query_pairs(pairs)
        assert out["r2"].iloc[0] == pytest.approx(0.64, abs=1.5e-5)
        assert out["status"].iloc[0] == ""
        assert out["r2"].iloc[1] == 1.0           # diagonal
        assert out["status"].iloc[1] == ""
        assert np.isnan(out["r2"].iloc[2]) and out["status"].iloc[2] == "absent"
        assert np.isnan(out["r2"].iloc[3]) and out["status"].iloc[3] == "not_in_panel"
        assert np.isnan(out["r2"].iloc[4]) and out["status"].iloc[4] == "cross_chromosome"

    def test_sign_harmonized_in_allele_aware_mode(self, tmp_path):
        panel = self._panel(tmp_path)
        # Panel pair (0,1)=SNP1 A/G, SNP2 C/T, stored SIGN=+ (r>=0).
        pairs = pd.DataFrame(
            {
                "CHR_1": [1, 1, 1],
                "POS_1": [100, 100, 100],
                "A1_1": ["A", "G", "G"],   # aligned / swapped / swapped
                "A2_1": ["G", "A", "A"],
                "CHR_2": [1, 1, 1],
                "POS_2": [200, 200, 200],
                "A1_2": ["C", "C", "T"],   # aligned / aligned / swapped
                "A2_2": ["T", "T", "C"],
            }
        )
        out = panel.query_pairs(pairs)
        assert out["sign"].tolist() == [1, -1, 1]  # 0 / 1 / 2 swaps -> +,-,+

    def test_base_mode_ignores_alleles_and_sign_is_na(self, tmp_path):
        panel = self._panel(tmp_path, mode="chr_pos")
        # With alleles supplied vs not supplied: identical r2, sign always NA.
        with_alleles = pd.DataFrame(
            {"CHR_1": [1], "POS_1": [100], "A1_1": ["A"], "A2_1": ["G"],
             "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"]}
        )
        without = pd.DataFrame({"CHR_1": [1], "POS_1": [100], "CHR_2": [1], "POS_2": [200]})
        out_a = panel.query_pairs(with_alleles)
        out_b = panel.query_pairs(without)
        assert out_a["r2"].iloc[0] == pytest.approx(out_b["r2"].iloc[0], abs=1.5e-5)
        assert pd.isna(out_a["sign"].iloc[0]) and pd.isna(out_b["sign"].iloc[0])
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestQueryPairs -v`
Expected: FAIL with `AttributeError: 'R2Panel' object has no attribute 'query_pairs'`.

- [ ] **Step 3: Write minimal implementation**

Add the resolution + query methods to `R2Panel` (inside the class):

```python
    def _resolve_endpoint(self, endpoint: pd.DataFrame):
        """Resolve one endpoint table to (chrom, row, panel_a1, panel_a2) arrays.

        ``endpoint`` carries canonical columns ``CHR/POS/SNP/A1/A2`` (subset by
        mode). Returns parallel arrays the length of ``endpoint``: chromosome
        label (object, ``None`` if unresolved), row index (int64, ``-1`` if
        unresolved), and panel A1/A2 (object, ``None`` if unresolved).
        """
        n = len(endpoint)
        chrom_out = np.full(n, None, dtype=object)
        idx_out = np.full(n, -1, dtype=np.int64)
        pa1_out = np.full(n, None, dtype=object)
        pa2_out = np.full(n, None, dtype=object)
        family = identity_mode_family(self._mode)
        keys = build_snp_id_series(endpoint, self._mode).to_numpy(dtype=object)

        if family == "chr_pos" or "CHR" in endpoint.columns:
            chroms = endpoint["CHR"].map(lambda v: normalize_chromosome(v)).to_numpy(dtype=object)
            for chrom in {c for c in chroms if c in self._paths}:
                rows = np.where(chroms == chrom)[0]
                state = self._chrom_state(chrom)
                positions = state.key_index.get_indexer(keys[rows])
                found = positions >= 0
                sel = rows[found]
                idx_out[sel] = positions[found]
                chrom_out[sel] = chrom
                pa1_out[sel] = state.a1[positions[found]]
                pa2_out[sel] = state.a2[positions[found]]
        else:
            # rsid family without CHR column: search every chromosome.
            for chrom in self._paths:
                state = self._chrom_state(chrom)
                positions = state.key_index.get_indexer(keys)
                found = (positions >= 0) & (idx_out < 0)
                idx_out[found] = positions[found]
                chrom_out[found] = chrom
                pa1_out[found] = state.a1[positions[found]]
                pa2_out[found] = state.a2[positions[found]]
        return chrom_out, idx_out, pa1_out, pa2_out

    def query_pairs(self, pairs: pd.DataFrame, *, with_r: bool = False,
                    strategy: str = "auto", strategy_threshold: int = 50_000) -> pd.DataFrame:
        """Return adjusted R², sign, and status for each SNP pair in ``pairs``.

        ``pairs`` has per-endpoint columns suffixed ``_1``/``_2``
        (``CHR/POS/A1/A2`` and/or ``SNP``). See the design spec for the full
        contract. ``with_r=True`` appends a signed Pearson ``r`` column.
        """
        ep1 = _endpoint_frame(pairs, "1", self._mode)
        ep2 = _endpoint_frame(pairs, "2", self._mode)
        c1, i1, qa1_1, qa2_1, pa1_1, pa2_1 = self._resolve_with_query_alleles(ep1)
        c2, i2, qa1_2, qa2_2, pa1_2, pa2_2 = self._resolve_with_query_alleles(ep2)

        n = len(pairs)
        r2 = np.full(n, np.nan, dtype=np.float32)
        sign = np.zeros(n, dtype=np.int8)
        status = np.array([""] * n, dtype=object)

        resolved = (i1 >= 0) & (i2 >= 0)
        status[~resolved] = "not_in_panel"
        cross = resolved & (c1 != c2)
        status[cross] = "cross_chromosome"
        same = resolved & (c1 == c2)
        diagonal = same & (i1 == i2)
        r2[diagonal] = 1.0
        sign[diagonal] = 1  # panel orientation; harmonized below
        lookup = same & (i1 != i2)

        for chrom in {c for c in c1[lookup]}:
            rows = np.where(lookup & (c1 == chrom))[0]
            state = self._chrom_state(chrom)
            lo = np.minimum(i1[rows], i2[rows])
            hi = np.maximum(i1[rows], i2[rows])
            chrom_r2, chrom_sign = lookup_pairs_in_parquet(
                state.parquet_file, lo, hi, n_snps=state.n_snps, r2_scale=state.r2_scale,
                strategy=strategy, strategy_threshold=strategy_threshold,
            )
            r2[rows] = chrom_r2
            sign[rows] = chrom_sign

        found_numeric = ~np.isnan(r2)
        absent = lookup & ~found_numeric
        status[absent] = "absent"

        if is_allele_aware_mode(self._mode):
            mult1 = _orientation_multiplier(qa1_1, pa1_1, pa2_1)
            mult2 = _orientation_multiplier(qa1_2, pa1_2, pa2_2)
            harmonized = (sign.astype(np.int16) * mult1 * mult2)
            sign = np.where(found_numeric, harmonized, 0).astype(np.int8)
        else:
            sign = np.zeros(n, dtype=np.int8)  # base mode: always NA

        out = pairs.copy()
        out["r2"] = r2
        sign_col = pd.array(sign, dtype="Int8")
        sign_col[sign == 0] = pd.NA
        out["sign"] = sign_col
        out["status"] = pd.array(status, dtype="string")
        if with_r:
            out["r"] = self._signed_r(r2, sign_col)
        return out

    def _resolve_with_query_alleles(self, endpoint: pd.DataFrame):
        """Resolve an endpoint and also return its query A1/A2 (or None columns)."""
        chrom, idx, pa1, pa2 = self._resolve_endpoint(endpoint)
        qa1 = endpoint["A1"].to_numpy(dtype=object) if "A1" in endpoint.columns else np.full(len(endpoint), None, dtype=object)
        qa2 = endpoint["A2"].to_numpy(dtype=object) if "A2" in endpoint.columns else np.full(len(endpoint), None, dtype=object)
        return chrom, idx, qa1, qa2, pa1, pa2

    def _signed_r(self, r2: np.ndarray, sign_col) -> np.ndarray:
        """Compute the signed r column from adjusted r2 and the sign column."""
        n_samples = self.n_samples
        if n_samples is None:
            raise LDSCInputError(
                "R2Panel cannot compute signed r without `ldsc:n_samples` in the panel "
                "metadata. Most likely a legacy panel was used. Omit `with_r`."
            )
        sign_float = np.where(sign_col.isna(), np.nan, sign_col.to_numpy(dtype="float64", na_value=np.nan))
        if not is_allele_aware_mode(self._mode):
            LOGGER.warning(
                "query-r2 `with_r` in base mode produces an all-NaN r column because base "
                "modes carry no sign. Use an allele-aware panel/mode for signed r."
            )
        return unbiased_r2_to_pearson_r(r2.astype(np.float64), n_samples, sign=sign_float)
```

Add the module-level helpers (`_endpoint_frame`, `_orientation_multiplier`) near
`_chrom_from_r2_path`:

```python
_ENDPOINT_COLUMNS = {"CHR": "CHR", "POS": "POS", "SNP": "SNP", "A1": "A1", "A2": "A2"}


def _endpoint_frame(pairs: pd.DataFrame, suffix: str, mode: str) -> pd.DataFrame:
    """Extract one endpoint's canonical columns from suffixed pair columns.

    Reads ``<NAME>_<suffix>`` columns and renames to canonical ``CHR/POS/SNP/A1/A2``.
    In base modes (``rsid``/``chr_pos``) allele columns are dropped before use, so
    a base-mode query behaves identically with or without alleles supplied.
    """
    allele_aware = is_allele_aware_mode(mode)
    out = {}
    for canonical in _ENDPOINT_COLUMNS:
        if canonical in ("A1", "A2") and not allele_aware:
            continue
        col = f"{canonical}_{suffix}"
        if col in pairs.columns:
            out[canonical] = pairs[col].to_numpy()
    if not out:
        raise LDSCInputError(
            f"query-r2 found no endpoint-{suffix} columns. Most likely the pair table is "
            f"missing `CHR_{suffix}/POS_{suffix}` or `SNP_{suffix}`. Provide the columns "
            "required by the panel's SNP identifier mode."
        )
    return pd.DataFrame(out)


def _orientation_multiplier(query_a1: np.ndarray, panel_a1: np.ndarray, panel_a2: np.ndarray) -> np.ndarray:
    """Return per-endpoint sign multiplier (+1 aligned, -1 swapped, 0 unknown).

    Aligned when the query A1 equals panel A1 (or its strand complement); swapped
    when it equals panel A2 (or its complement). Strand-ambiguous panel SNPs are
    excluded from allele-aware matching, so a matched endpoint always classifies
    cleanly; ``0`` only occurs for unresolved endpoints (panel allele ``None``).
    """
    n = query_a1.shape[0]
    out = np.zeros(n, dtype=np.int16)
    for k in range(n):
        q, p1, p2 = query_a1[k], panel_a1[k], panel_a2[k]
        if q is None or p1 is None or p2 is None or (isinstance(q, float) and np.isnan(q)):
            continue
        q = str(q).upper()
        p1u, p2u = str(p1).upper(), str(p2).upper()
        if q == p1u or q == p1u.translate(COMPLEMENT):
            out[k] = 1
        elif q == p2u or q == p2u.translate(COMPLEMENT):
            out[k] = -1
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestQueryPairs tests/test_r2_query.py::TestR2PanelOpen -v`
Expected: PASS (all). Re-run the full file: `pytest tests/test_r2_query.py -v`.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/r2_query.py tests/test_r2_query.py
git commit -m "feat(r2-query): implement R2Panel.query_pairs with status and sign"
```

---

## Task 5: `query_r2` one-shot wrapper + `with_r`

**Files:**
- Modify: `src/ldsc/r2_query.py`
- Test: `tests/test_r2_query.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_r2_query.py  (append)
class TestQueryR2Wrapper:
    def test_query_r2_one_shot_matches_handle(self, tmp_path):
        from ldsc.r2_query import query_r2

        build_test_panel(tmp_path, snp_identifier="chr_pos_allele_aware")
        pairs = pd.DataFrame(
            {"CHR_1": [1], "POS_1": [100], "A1_1": ["A"], "A2_1": ["G"],
             "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"]}
        )
        out = query_r2(pairs, panel_dir=tmp_path, genome_build="hg38")
        assert out["r2"].iloc[0] == pytest.approx(0.64, abs=1.5e-5)
        assert "status" in out.columns

    def test_with_r_adds_signed_correlation(self, tmp_path):
        from ldsc.r2_query import query_r2

        build_test_panel(tmp_path, snp_identifier="chr_pos_allele_aware", n_samples=1000)
        pairs = pd.DataFrame(
            {"CHR_1": [1], "POS_1": [100], "A1_1": ["G"], "A2_1": ["A"],  # swapped -> negative r
             "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"]}
        )
        out = query_r2(pairs, panel_dir=tmp_path, genome_build="hg38", with_r=True)
        assert "r" in out.columns
        assert out["r"].iloc[0] < 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestQueryR2Wrapper -v`
Expected: FAIL with `ImportError: cannot import name 'query_r2'`.

- [ ] **Step 3: Write minimal implementation**

```python
# src/ldsc/r2_query.py  (append, module level)
def query_r2(
    pairs: pd.DataFrame,
    *,
    panel_dir=None,
    meta_path=None,
    parquet_path=None,
    snp_identifier: str | None = None,
    genome_build: str | None = None,
    with_r: bool = False,
    strategy: str = "auto",
    strategy_threshold: int = 50_000,
) -> pd.DataFrame:
    """Open a panel, query ``pairs`` once, and return the annotated DataFrame.

    Thin one-shot wrapper over :meth:`R2Panel.open` + :meth:`R2Panel.query_pairs`.
    """
    panel = R2Panel.open(
        panel_dir,
        meta_path=meta_path,
        parquet_path=parquet_path,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
    )
    return panel.query_pairs(
        pairs, with_r=with_r, strategy=strategy, strategy_threshold=strategy_threshold
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestQueryR2Wrapper -v`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/r2_query.py tests/test_r2_query.py
git commit -m "feat(r2-query): add query_r2 one-shot wrapper"
```

---

## Task 6: CLI `ldsc query-r2`

**Files:**
- Modify: `src/ldsc/r2_query.py` (argparse + run-from-args + main)
- Modify: `src/ldsc/cli.py` (register + dispatch)
- Modify: `CLAUDE.md` (CLI-surface invariant)
- Test: `tests/test_r2_query.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_r2_query.py  (append)
class TestQueryR2CLI:
    def test_cli_reads_tsv_and_writes_tsv(self, tmp_path):
        from ldsc.cli import main as cli_main

        build_test_panel(tmp_path, snp_identifier="chr_pos_allele_aware")
        pairs_path = tmp_path / "pairs.tsv"
        pd.DataFrame(
            {"CHR_1": [1], "POS_1": [100], "A1_1": ["A"], "A2_1": ["G"],
             "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"]}
        ).to_csv(pairs_path, sep="\t", index=False)
        out_path = tmp_path / "out.tsv"
        cli_main([
            "query-r2", "--panel-dir", str(tmp_path), "--genome-build", "hg38",
            "--pairs", str(pairs_path), "--out", str(out_path), "--with-r",
        ])
        result = pd.read_csv(out_path, sep="\t")
        assert result["r2"].iloc[0] == pytest.approx(0.64, abs=1e-4)
        assert "status" in result.columns and "r" in result.columns
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestQueryR2CLI -v`
Expected: FAIL — `query-r2` is not a recognized command (SystemExit from argparse).

- [ ] **Step 3a: Add the parser/runner/main to `src/ldsc/r2_query.py`**

```python
# src/ldsc/r2_query.py  (append, module level)
import argparse
import sys
from typing import Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the ``ldsc query-r2`` argument parser."""
    parser = argparse.ArgumentParser(
        prog="ldsc query-r2",
        description="Query adjusted R2 (and signed r) for SNP pairs from a reference panel.",
        allow_abbrev=False,
    )
    parser.add_argument("--panel-dir", default=None, help="build-ref-panel output directory.")
    parser.add_argument("--meta", default=None, help="Explicit chrN_meta.tsv.gz path (with --parquet).")
    parser.add_argument("--parquet", default=None, help="Explicit chrN_r2.parquet path (with --meta).")
    parser.add_argument("--pairs", required=True, help="TSV/CSV of pairs with _1/_2 endpoint columns ('-' = stdin).")
    parser.add_argument("--out", default=None, help="Output TSV path (default: stdout).")
    parser.add_argument("--snp-identifier", default=None, help="Override panel SNP identifier mode.")
    parser.add_argument("--genome-build", choices=["hg19", "hg38"], default=None, help="Genome build for sub-dir resolution.")
    parser.add_argument("--with-r", action="store_true", help="Add a signed Pearson r column.")
    parser.add_argument("--strategy", choices=["auto", "random", "stream"], default="auto", help="Lookup strategy.")
    parser.add_argument("--strategy-threshold", type=int, default=50_000, help="Auto-select pair-count threshold.")
    return parser


def run_query_r2_from_args(args: argparse.Namespace) -> pd.DataFrame:
    """Run ``query-r2`` from parsed CLI arguments and write the output table."""
    pairs_handle = sys.stdin if args.pairs == "-" else args.pairs
    sep = "," if str(args.pairs).endswith(".csv") else "\t"
    pairs = pd.read_csv(pairs_handle, sep=sep)
    result = query_r2(
        pairs,
        panel_dir=args.panel_dir,
        meta_path=args.meta,
        parquet_path=args.parquet,
        snp_identifier=args.snp_identifier,
        genome_build=args.genome_build,
        with_r=args.with_r,
        strategy=args.strategy,
        strategy_threshold=args.strategy_threshold,
    )
    if args.out is None:
        result.to_csv(sys.stdout, sep="\t", index=False)
    else:
        result.to_csv(args.out, sep="\t", index=False)
    return result


def main(argv: Sequence[str] | None = None) -> pd.DataFrame:
    """Command-line entry point for the R2 pair-query workflow."""
    args = build_parser().parse_args(argv)
    return run_query_r2_from_args(args)
```

- [ ] **Step 3b: Register and dispatch in `src/ldsc/cli.py`**

In `build_parser()` (after the `rg_parser` block, before `return parser`), add:

```python
    query_r2_parser = subparsers.add_parser("query-r2", help="Query R2 for SNP pairs from a reference panel.")
    _copy_actions(query_r2_parser, r2_query.build_parser())
```

Add the import near the other workflow imports at the top of `cli.py`:

```python
from . import r2_query
```

In `main()`, add an early-dispatch branch alongside the others (after the
`build-ref-panel` branch in the `if argv:` block):

```python
        if command == "query-r2":
            return r2_query.main(subargv)
```

And in the parsed-args dispatch (after the `build-ref-panel` branch):

```python
    if args.command == "query-r2":
        return r2_query.run_query_r2_from_args(args)
```

- [ ] **Step 3c: Update the CLI-surface invariant in `CLAUDE.md`**

Change the invariant line:

```
- Keep one CLI surface: `ldsc` with subcommands `annotate`, `build-ref-panel`,
  `ldscore`, `munge-sumstats`, `h2`, `partitioned-h2`, and `rg`.
```

to include `query-r2`:

```
- Keep one CLI surface: `ldsc` with subcommands `annotate`, `build-ref-panel`,
  `ldscore`, `munge-sumstats`, `h2`, `partitioned-h2`, `rg`, and `query-r2`.
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestQueryR2CLI -v`
Then sanity-check help: `ldsc query-r2 --help`
Expected: test PASS; help prints the options.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/r2_query.py src/ldsc/cli.py CLAUDE.md tests/test_r2_query.py
git commit -m "feat(r2-query): add ldsc query-r2 CLI subcommand"
```

---

## Task 7: Public package exports

**Files:**
- Modify: `src/ldsc/__init__.py`
- Test: `tests/test_r2_query.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_r2_query.py  (append)
class TestPackageExports:
    def test_public_symbols_importable_from_ldsc(self):
        import ldsc

        assert hasattr(ldsc, "R2Panel")
        assert hasattr(ldsc, "query_r2")
        assert hasattr(ldsc, "unbiased_r2_to_pearson_r")
        for name in ("R2Panel", "query_r2", "unbiased_r2_to_pearson_r"):
            assert name in ldsc.__all__
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_r2_query.py::TestPackageExports -v`
Expected: FAIL with `AssertionError` (symbols absent).

- [ ] **Step 3: Write minimal implementation**

In `src/ldsc/__init__.py`, add the import alongside the other concrete imports
(e.g. after the `ref_panel` import line):

```python
from .r2_query import R2Panel, query_r2, unbiased_r2_to_pearson_r
```

Add the three names to the `__all__` list (alphabetical-ish, near the others):

```python
    "R2Panel",
    "query_r2",
    "unbiased_r2_to_pearson_r",
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_r2_query.py::TestPackageExports -v`
Then check no import cycle / surface regression:
`pytest tests/test_package_layout.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/__init__.py tests/test_r2_query.py
git commit -m "feat(r2-query): export R2Panel, query_r2, converter from package"
```

---

## Task 8: Documentation

**Files:**
- Create: `docs/current/ref-panel-r2-query.md`
- Modify: `design_map.md`

- [ ] **Step 1: Write the user-facing doc**

Create `docs/current/ref-panel-r2-query.md` with these sections (prose, no
placeholders): purpose; Python API (`R2Panel.open`, `query_pairs`, `query_r2`,
`unbiased_r2_to_pearson_r`) with a runnable snippet; the pair-input column
contract (`_1`/`_2`, modes, base-mode-ignores-alleles); the output schema
(`r2`, `sign`, `status`, optional `r`) with the status vocabulary
(`not_in_panel`/`cross_chromosome`/`absent`); the sign rule (allele-aware only;
base modes `NA`); the lookup strategies and `--strategy`/`--strategy-threshold`;
and the CLI usage with a worked `ldsc query-r2` example. Mirror the spec
`docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md` for exact
semantics.

Python snippet to include verbatim:

```python
import pandas as pd
from ldsc import query_r2, unbiased_r2_to_pearson_r

pairs = pd.DataFrame({
    "CHR_1": [1, 1], "POS_1": [752721, 776546], "A1_1": ["A", "G"], "A2_1": ["G", "A"],
    "CHR_2": [1, 1], "POS_2": [777122, 798959], "A1_2": ["A", "T"], "A2_2": ["C", "G"],
})
result = query_r2(pairs, panel_dir="ref_panel", genome_build="hg38", with_r=True)
print(result[["r2", "sign", "status", "r"]])
```

- [ ] **Step 2: Update `design_map.md`**

Add a row/section mapping:
- spec `docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md` and plan
  `docs/superpowers/plans/2026-06-06-ref-panel-r2-query-plan.md`
- → `src/ldsc/r2_query.py` (`R2Panel`, `query_r2`, `unbiased_r2_to_pearson_r`,
  CLI) and `src/ldsc/_kernel/r2_query.py` (`lookup_pairs_in_parquet`).

Match the existing `design_map.md` formatting.

- [ ] **Step 3: Verify docs reference real symbols**

Run: `grep -n "R2Panel\|query_r2\|unbiased_r2_to_pearson_r\|lookup_pairs_in_parquet" docs/current/ref-panel-r2-query.md design_map.md`
Expected: every referenced symbol exists in the source modules.

- [ ] **Step 4: Run the full suite**

Run: `pytest tests/test_r2_query.py -v && pytest -q`
Expected: new file fully green; full suite no new failures.

- [ ] **Step 5: Commit**

```bash
git add docs/current/ref-panel-r2-query.md design_map.md
git commit -m "docs(r2-query): document R2 pair query API and CLI"
```

---

## Self-Review

**1. Spec coverage:**
- §4 module layout → Tasks 1–8 create exactly the listed files. ✓
- §5 `R2Panel.open`, two input modes, lazy load, binding, attributes → Task 3. ✓
- §6 `_1`/`_2` input, base-mode drops alleles → Task 4 (`_endpoint_frame`). ✓
- §7 routing + status vocabulary (diagonal/stored/absent/cross/not-in-panel) → Task 4 test `test_stored_diagonal_absent_cross_and_missing`. ✓
- §8 int64-key match, random/stream/auto, row-group pruning → Task 2. ✓
- §9 sign harmonization allele-aware only, base → `NA`, parity rule → Task 4 tests. ✓
- §10 output schema (`r2`/`sign`/`status`/optional `r`) → Task 4/5. ✓
- §11 converter + `with_r` + base-mode warning + missing-`n_samples` error → Tasks 1, 5, and `_signed_r`. ✓
- §12 CLI → Task 6. ✓
- §13 errors (both/neither mode, binding, key collapse, missing alleles, `with_r` no `n_samples`) → Tasks 3/4 (`_chrom_state`, `_endpoint_frame`, `_signed_r`). ✓
- §14 tests (converter, key resolution, sign, integration, parity, CLI) → Tasks 1–6. ✓
- §16 reused helpers → imported in Task 3. ✓

**2. Placeholder scan:** No TBD/TODO; every code step contains complete code; the
doc task (Task 8) lists exact sections + a verbatim snippet rather than "write
docs". ✓

**3. Type consistency:** `lookup_pairs_in_parquet(pf, i, j, *, n_snps, r2_scale, strategy, strategy_threshold)` is called with the same keywords in Tasks 2 and 4. `R2Panel.open` / `query_pairs` / `query_r2` signatures match across Tasks 3–6. `_orientation_multiplier(query_a1, panel_a1, panel_a2)` returns int16 ±1/0, consumed as `mult1`/`mult2` in Task 4. `sign` sentinel `0`→`pd.NA` is consistent (kernel returns 0 for absent; handle maps 0→NA). ✓

**Note on task ordering:** Task 3's `test_binding_mismatch_is_hard_error` calls
`query_pairs`, implemented in Task 4. When executing in order, implement Task 4
before running Task 3's full test set (both modify the same file and are committed
sequentially; the binding test passes once Task 4 lands).
