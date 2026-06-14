# Observed-to-Liability Scale Conversion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Report binary-trait heritability and genetic covariance on both the observed and liability scales across the `h2`, `partitioned-h2`, and `rg` modules, driven by `--samp-prev` / `--pop-prev` (and an `rg` `--prevalence-manifest`), matching legacy ldsc2 numerics.

**Architecture:** A single vectorizable kernel primitive `liability_conversion_factor(P, K)` computes the conversion factor `c`; `h2_obs_to_liab` / `gencov_obs_to_liab` delegate to it. A new `src/ldsc/prevalence.py` workflow module parses/validates the CLI inputs (scalar for h2/partitioned; positional comma-list or manifest TSV for rg) into a normalized per-trait structure. The workflow summary functions (`summarize_total_h2`, `summarize_partitioned_h2`, `_summarize_rg_pair`) gain `_obs`/`_liab` column pairs plus prevalence columns; conversion is a deterministic scalar multiply applied there.

**Tech Stack:** Python 3, numpy, scipy.stats, pandas, argparse, pytest (+ stdlib unittest compatibility).

**Environment for every command below:**
```bash
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
```

**Reference spec:** `docs/superpowers/specs/2026-06-13-liability-scale-conversion-design.md`

**Caveat:** `src/ldsc/regression_runner.py` and `tests/test_regression_workflow.py` carry unrelated uncommitted edits from a parallel session (post-chisq-filter `n_snps` reporting). Build on top of them; stage only the files each task changes. Do not revert them.

---

## File Structure

| File | Responsibility | Action |
| --- | --- | --- |
| `src/ldsc/_kernel/regression.py` | `liability_conversion_factor` primitive; `h2_obs_to_liab`/`gencov_obs_to_liab` delegate to it | Modify |
| `src/ldsc/prevalence.py` | Parse/validate scalar + positional + manifest prevalence into normalized per-trait `(P, K)` | Create |
| `src/ldsc/outputs.py` | `PARTITIONED_H2_COLUMNS`, `RG_FULL_COLUMNS` schema constants | Modify |
| `src/ldsc/regression_runner.py` | Summary functions, conversion helper, CLI flags, runner wiring, metadata, logging | Modify |
| `src/ldsc/config.py` | `RegressionConfig.samp_prev`/`pop_prev` docstring (valid range) | Modify |
| `tests/test_kernel_regression.py` | Primitive + delegate tests | Modify |
| `tests/test_prevalence.py` | Parser/validation tests | Create |
| `tests/test_output.py` | Column-schema tests for renamed/added columns | Modify |
| `tests/test_regression_workflow.py` | Summary + CLI + rg threading tests | Modify |
| `docs/current/regression-configuration.md` | Replace 1.4 "reserved" row with real behavior | Modify |
| `docs/current/io-argument-inventory.md` | Add `--samp-prev`/`--pop-prev`/`--prevalence-manifest` rows | Modify |
| `tutorials/liability-scale-conversion.md` | User-facing how-to: worked examples for every prevalence input form | Create |

---

## Task 1: Vectorizable `liability_conversion_factor` primitive

**Files:**
- Modify: `src/ldsc/_kernel/regression.py` (replace `h2_obs_to_liab` body at lines 138-176; add new primitive above `gencov_obs_to_liab` at line 106)
- Test: `tests/test_kernel_regression.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_kernel_regression.py` (top-level imports: `import math`, `import numpy as np`, `import pytest`, `from ldsc._kernel import regression as reg`):

```python
def test_liability_conversion_factor_known_value():
    # c(P=0.5, K=0.01) hand-computed from K^2(1-K)^2 / (P(1-P) phi(isf(K))^2).
    c = reg.liability_conversion_factor(0.5, 0.01)
    assert math.isclose(c, 0.5519, rel_tol=0, abs_tol=1e-3)


def test_liability_conversion_factor_both_nan_is_one():
    assert reg.liability_conversion_factor(float("nan"), float("nan")) == 1.0


def test_liability_conversion_factor_vectorizes_over_K():
    K = np.array([0.02, 0.05, 0.10, 0.20])
    c = reg.liability_conversion_factor(0.5, K)
    assert c.shape == K.shape
    # Liability h2 = c * h2_obs is monotone increasing in K over this range.
    assert np.all(np.diff(c) > 0)


def test_liability_conversion_factor_rejects_out_of_range():
    with pytest.raises(Exception):
        reg.liability_conversion_factor(0.5, 1.0)
    with pytest.raises(Exception):
        reg.liability_conversion_factor(0.0, 0.01)


def test_liability_conversion_factor_rejects_half_specified():
    # One of (P, K) NaN while the other is finite is underspecified.
    with pytest.raises(Exception):
        reg.liability_conversion_factor(float("nan"), 0.01)
    with pytest.raises(Exception):
        reg.liability_conversion_factor(0.5, float("nan"))


def test_h2_obs_to_liab_matches_primitive():
    assert math.isclose(
        reg.h2_obs_to_liab(0.3, 0.5, 0.01),
        0.3 * reg.liability_conversion_factor(0.5, 0.01),
        rel_tol=1e-12,
    )


def test_gencov_obs_to_liab_uses_sqrt_product():
    g = reg.gencov_obs_to_liab(0.1, 0.5, 0.5, 0.01, 0.02)
    expected = 0.1 * math.sqrt(reg.liability_conversion_factor(0.5, 0.01)) * math.sqrt(
        reg.liability_conversion_factor(0.5, 0.02)
    )
    assert math.isclose(g, expected, rel_tol=1e-12)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `PYTHONPATH=src pytest tests/test_kernel_regression.py -k "liability or obs_to_liab or sqrt_product" -v`
Expected: FAIL -- `AttributeError: module ... has no attribute 'liability_conversion_factor'`.

- [ ] **Step 3: Implement the primitive and delegate**

In `src/ldsc/_kernel/regression.py`, add this function immediately above `def gencov_obs_to_liab` (line 106). Use `my-skills:fun-doc` numpydoc style:

```python
def liability_conversion_factor(samp_prev, pop_prev):
    r"""Observed-to-liability heritability conversion factor :math:`c(P, K)`.

    Computes the multiplicative factor that rescales an observed-scale
    heritability (estimated in an ascertained case-control sample) to the
    liability scale in the population, following Lee et al. 2011 as used by
    legacy LDSC:

    .. math::
        c(P, K) = \frac{K^2 (1 - K)^2}{P (1 - P)\, \phi\!\big(\Phi^{-1}(1 - K)\big)^2}

    where :math:`P` is the sample (case) prevalence, :math:`K` the population
    prevalence, :math:`\phi` the standard normal pdf, and :math:`\Phi^{-1}` the
    inverse CDF. Vectorized: ``samp_prev`` and ``pop_prev`` broadcast, so a fixed
    ``P`` against an array of ``K`` yields the liability-h2 curve used by
    sensitivity analyses.

    Parameters
    ----------
    samp_prev, pop_prev : float or array_like
        Sample and population prevalences, each a probability in the open
        interval ``(0, 1)`` for a binary trait, or NaN for a quantitative trait
        (which yields a factor of ``1.0``). Per element, both must be finite or
        both NaN.

    Returns
    -------
    float or numpy.ndarray
        The factor ``c``; ``1.0`` where both inputs are NaN. Scalar inputs
        return a Python float.

    Raises
    ------
    LDSCConfigError
        If any element of ``samp_prev`` or ``pop_prev`` is outside ``(0, 1)``
        and not NaN, or if a single element has exactly one of ``P``/``K`` NaN.
    """
    P = np.asarray(samp_prev, dtype=np.float64)
    K = np.asarray(pop_prev, dtype=np.float64)
    P_b, K_b = np.broadcast_arrays(P, K)
    nan_P, nan_K = np.isnan(P_b), np.isnan(K_b)
    if np.any(nan_P != nan_K):
        raise LDSCConfigError(
            "Liability-scale conversion received a trait with exactly one of sample/population "
            f"prevalence missing (samp_prev={samp_prev!r}, pop_prev={pop_prev!r}). "
            "Most likely a `nan` was placed in only one of `--samp-prev`/`--pop-prev`. "
            "Set both to numbers in (0, 1) for a binary trait, or both to `nan` for a quantitative trait."
        )
    finite = ~nan_P
    bad = finite & (((P_b <= 0) | (P_b >= 1)) | ((K_b <= 0) | (K_b >= 1)))
    if np.any(bad):
        raise LDSCConfigError(
            "Liability-scale conversion received a prevalence outside the open interval (0, 1) "
            f"(samp_prev={samp_prev!r}, pop_prev={pop_prev!r}). "
            "Most likely `--samp-prev`/`--pop-prev` was given a value <= 0 or >= 1. "
            "Pass a probability strictly between 0 and 1, or `nan` for a quantitative trait."
        )
    out = np.ones(K_b.shape, dtype=np.float64)
    if np.any(finite):
        Pf, Kf = P_b[finite], K_b[finite]
        thresh = norm.isf(Kf)
        out[finite] = Kf**2 * (1 - Kf) ** 2 / (Pf * (1 - Pf) * norm.pdf(thresh) ** 2)
    if out.ndim == 0 or out.size == 1 and np.ndim(samp_prev) == 0 and np.ndim(pop_prev) == 0:
        return float(out.reshape(()))
    return out
```

Then replace the body of `h2_obs_to_liab` (lines 138-176) so it delegates, keeping its public signature and legacy NaN semantics:

```python
def h2_obs_to_liab(h2_obs, P, K):
    """Convert observed-scale heritability to the liability scale.

    Thin wrapper over :func:`liability_conversion_factor` preserving the legacy
    LDSC signature: ``h2_liab = h2_obs * c(P, K)``. When both ``P`` and ``K`` are
    NaN the observed value is returned unchanged.

    Parameters
    ----------
    h2_obs : float
        Observed-scale heritability in an ascertained sample.
    P : float in (0, 1)
        Sample (case) prevalence, or NaN for a quantitative trait.
    K : float in (0, 1)
        Population prevalence, or NaN for a quantitative trait.

    Returns
    -------
    float
        Liability-scale heritability.
    """
    return h2_obs * liability_conversion_factor(P, K)
```

Leave `gencov_obs_to_liab` (lines 106-135) as-is structurally -- it already calls `h2_obs_to_liab` and therefore now delegates transitively. (Confirm its `import` of `norm` is unchanged; `norm` is imported at module top, line 13.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src pytest tests/test_kernel_regression.py -v`
Expected: PASS (all, including pre-existing).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/regression.py tests/test_kernel_regression.py
git commit -m "feat(regression): add vectorizable liability_conversion_factor primitive"
```

---

## Task 2: Prevalence parsing/validation module

**Files:**
- Create: `src/ldsc/prevalence.py`
- Test: `tests/test_prevalence.py`

The normalized type is `TraitPrevalence = tuple[float | None, float | None]` (samp_prev, pop_prev); `None` marks a quantitative trait. Range/consistency validation reuses `liability_conversion_factor` (calling it once per resolved trait fails fast on bad values).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_prevalence.py`:

```python
import math

import pytest

from ldsc import prevalence
from ldsc.errors import LDSCConfigError, LDSCUsageError


def test_scalar_both_none_is_quantitative():
    assert prevalence.parse_scalar_prevalence(None, None) == (None, None)


def test_scalar_pair():
    assert prevalence.parse_scalar_prevalence(0.5, 0.01) == (0.5, 0.01)


def test_scalar_half_specified_raises():
    with pytest.raises(LDSCUsageError):
        prevalence.parse_scalar_prevalence(0.5, None)


def test_scalar_out_of_range_raises():
    with pytest.raises(LDSCConfigError):
        prevalence.parse_scalar_prevalence(0.5, 1.0)


def test_positional_list_aligns_and_marks_qt():
    out = prevalence.parse_positional_prevalences("0.5,0.5,nan", "0.01,0.02,nan", n_traits=3)
    assert out == [(0.5, 0.01), (0.5, 0.02), (None, None)]


def test_positional_empty_token_is_quantitative():
    out = prevalence.parse_positional_prevalences("0.5,,0.3", "0.01,,0.02", n_traits=3)
    assert out == [(0.5, 0.01), (None, None), (0.3, 0.02)]


def test_positional_length_mismatch_raises():
    with pytest.raises(LDSCUsageError):
        prevalence.parse_positional_prevalences("0.5,0.5", "0.01,0.01", n_traits=3)


def test_positional_half_specified_trait_raises():
    with pytest.raises((LDSCUsageError, LDSCConfigError)):
        prevalence.parse_positional_prevalences("0.5,nan", "0.01,0.01", n_traits=2)


def test_manifest_superset_lookup_by_name(tmp_path):
    manifest = tmp_path / "prev.tsv"
    manifest.write_text(
        "# standing repository\n"
        "trait_name\tsamp_prev\tpop_prev\n"
        "scz2022\t0.43\t0.01\n"
        "scz2014\t0.40\t0.01\n"
        "height\tnan\tnan\n",
        encoding="utf-8",
    )
    out = prevalence.parse_manifest_prevalences(
        str(manifest), trait_names=["scz2022", "height"], trait_paths=["a", "b"]
    )
    assert out == [(0.43, 0.01), (None, None)]


def test_manifest_whitespace_delimiter(tmp_path):
    manifest = tmp_path / "prev.txt"
    manifest.write_text("trait_name samp_prev pop_prev\nscz2022   0.43   0.01\n", encoding="utf-8")
    out = prevalence.parse_manifest_prevalences(
        str(manifest), trait_names=["scz2022"], trait_paths=["a"]
    )
    assert out == [(0.43, 0.01)]


def test_manifest_missing_trait_raises(tmp_path):
    manifest = tmp_path / "prev.tsv"
    manifest.write_text("trait_name\tsamp_prev\tpop_prev\nscz2022\t0.43\t0.01\n", encoding="utf-8")
    with pytest.raises(LDSCUsageError):
        prevalence.parse_manifest_prevalences(
            str(manifest), trait_names=["scz2022", "bip"], trait_paths=["a", "b"]
        )


def test_resolve_rg_rejects_both_schemes():
    with pytest.raises(LDSCUsageError):
        prevalence.resolve_rg_prevalences(
            samp_prev="0.5,0.5",
            pop_prev="0.01,0.01",
            manifest_path="x.tsv",
            trait_names=["a", "b"],
            trait_paths=["a", "b"],
        )


def test_resolve_rg_none_returns_none():
    assert (
        prevalence.resolve_rg_prevalences(
            samp_prev=None, pop_prev=None, manifest_path=None,
            trait_names=["a", "b"], trait_paths=["a", "b"],
        )
        is None
    )


def test_resolve_rg_manifest_duplicate_name_raises(tmp_path):
    manifest = tmp_path / "prev.tsv"
    manifest.write_text("trait_name\tsamp_prev\tpop_prev\nscz\t0.43\t0.01\n", encoding="utf-8")
    with pytest.raises(LDSCUsageError):
        prevalence.resolve_rg_prevalences(
            samp_prev=None, pop_prev=None, manifest_path=str(manifest),
            trait_names=["scz", "scz"], trait_paths=["/data/a.sumstats", "/data/b.sumstats"],
        )
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `PYTHONPATH=src pytest tests/test_prevalence.py -v`
Expected: FAIL -- `ModuleNotFoundError: No module named 'ldsc.prevalence'`.

- [ ] **Step 3: Implement `src/ldsc/prevalence.py`**

Create the module (module header via `my-skills:fun-doc`):

```python
"""
prevalence.py

Overview
--------
Parse and validate the binary-trait prevalence inputs (`--samp-prev`,
`--pop-prev`, and the rg `--prevalence-manifest`) into a normalized per-trait
structure consumed by the regression summary layer for observed-to-liability
scale conversion. Validation is fail-fast: range and consistency are checked
here, before LD scores are loaded.

Key Functions
-------------
parse_scalar_prevalence : single-trait (h2, partitioned-h2) scalar pair.
resolve_rg_prevalences : dispatch rg inputs (positional comma-list xor manifest)
    to a per-trait list aligned to the resolved trait order.

Design Notes
------------
- Normalized unit is ``TraitPrevalence = (samp_prev | None, pop_prev | None)``;
  ``None`` marks a quantitative trait (no liability scale).
- Range/consistency validation delegates to
  ``_kernel.regression.liability_conversion_factor`` so the (0, 1)-or-NaN and
  both-or-neither rules live in one place.
- Manifest lookup is by exact munged trait name; supersets are allowed, missing
  resolved traits and (under the manifest scheme) duplicate resolved names are
  hard errors.

Dependencies
------------
- numpy
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .errors import LDSCConfigError, LDSCUsageError
from ._kernel.regression import liability_conversion_factor

TraitPrevalence = tuple[float | None, float | None]

_MISSING_TOKENS = {"", "nan", "na", "none"}


def _token_to_value(token: str) -> float | None:
    """Map one comma-list/manifest token to a float or None (quantitative)."""
    text = token.strip()
    if text.lower() in _MISSING_TOKENS:
        return None
    try:
        return float(text)
    except ValueError as exc:
        raise LDSCUsageError(
            f"Could not parse prevalence value {token!r} as a number. "
            "Most likely a malformed `--samp-prev`/`--pop-prev` entry or manifest cell. "
            "Use a probability in (0, 1), or `nan` for a quantitative trait."
        ) from exc


def _validate_pair(samp_prev: float | None, pop_prev: float | None) -> TraitPrevalence:
    """Enforce both-or-neither and (0, 1) range for one trait; return normalized pair."""
    both_present = samp_prev is not None and pop_prev is not None
    both_absent = samp_prev is None and pop_prev is None
    if not (both_present or both_absent):
        raise LDSCUsageError(
            f"A trait has exactly one of sample/population prevalence specified "
            f"(samp_prev={samp_prev!r}, pop_prev={pop_prev!r}). "
            "Set both for a binary trait, or leave both unset (`nan`) for a quantitative trait."
        )
    if both_absent:
        return (None, None)
    # Reuse the kernel's (0, 1) validation; raises LDSCConfigError on bad input.
    liability_conversion_factor(samp_prev, pop_prev)
    return (float(samp_prev), float(pop_prev))


def parse_scalar_prevalence(samp_prev, pop_prev) -> TraitPrevalence:
    """Validate the scalar `--samp-prev`/`--pop-prev` pair for h2/partitioned-h2.

    ``samp_prev``/``pop_prev`` are floats (possibly NaN) or None. NaN is treated
    as None (quantitative). Returns the normalized ``(P, K)`` pair.
    """
    sp = None if samp_prev is None or (isinstance(samp_prev, float) and math.isnan(samp_prev)) else samp_prev
    pp = None if pop_prev is None or (isinstance(pop_prev, float) and math.isnan(pop_prev)) else pop_prev
    return _validate_pair(sp, pp)


def parse_positional_prevalences(samp_prev: str, pop_prev: str, n_traits: int) -> list[TraitPrevalence]:
    """Parse aligned comma-separated `--samp-prev`/`--pop-prev` lists for rg.

    Both strings must split into exactly ``n_traits`` tokens. Returns a list of
    normalized ``(P, K)`` pairs aligned to the resolved trait order.
    """
    if samp_prev is None or pop_prev is None:
        raise LDSCUsageError(
            "rg liability conversion needs both `--samp-prev` and `--pop-prev` (or neither). "
            "Most likely only one positional list was supplied."
        )
    sp = [_token_to_value(t) for t in samp_prev.split(",")]
    pp = [_token_to_value(t) for t in pop_prev.split(",")]
    for name, values in (("--samp-prev", sp), ("--pop-prev", pp)):
        if len(values) != n_traits:
            raise LDSCUsageError(
                f"rg {name} has {len(values)} values but {n_traits} traits resolved from "
                "`--sumstats-sources`. The comma-separated list must have exactly one entry per "
                "resolved trait, in resolved order. Re-check the order printed in the run log, or "
                "use `--prevalence-manifest` to match by trait name."
            )
    return [_validate_pair(s, p) for s, p in zip(sp, pp)]


def parse_manifest_prevalences(
    manifest_path: str, trait_names: list[str], trait_paths: list[str]
) -> list[TraitPrevalence]:
    """Look up `(P, K)` per resolved trait from a prevalence manifest TSV.

    The manifest is whitespace/tab-delimited with a ``trait_name samp_prev
    pop_prev`` header; ``#`` comment lines are ignored. Extra rows (traits not in
    this run) are allowed; every resolved trait must have a matching row, matched
    by exact name. Returns pairs aligned to ``trait_names``.
    """
    table: dict[str, TraitPrevalence] = {}
    path = Path(manifest_path)
    if not path.is_file():
        raise LDSCUsageError(
            f"rg could not read `--prevalence-manifest` at '{manifest_path}': not a file. "
            "Pass a whitespace- or tab-delimited table with columns trait_name, samp_prev, pop_prev."
        )
    header: list[str] | None = None
    idx: dict[str, int] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.lstrip().startswith("#") or not raw.strip():
            continue
        fields = raw.split()
        if header is None:
            header = [f.strip().lower() for f in fields]
            for required in ("trait_name", "samp_prev", "pop_prev"):
                if required not in header:
                    raise LDSCUsageError(
                        f"Prevalence manifest '{manifest_path}' is missing a '{required}' column. "
                        "Required header columns are trait_name, samp_prev, pop_prev."
                    )
            idx = {col: header.index(col) for col in ("trait_name", "samp_prev", "pop_prev")}
            continue
        name = fields[idx["trait_name"]].strip()
        pair = _validate_pair(
            _token_to_value(fields[idx["samp_prev"]]),
            _token_to_value(fields[idx["pop_prev"]]),
        )
        table[name] = pair
    missing = [t for t in trait_names if t not in table]
    if missing:
        annotated = ", ".join(f"{t} ({p})" for t, p in zip(trait_names, trait_paths) if t in missing)
        raise LDSCUsageError(
            f"Prevalence manifest '{manifest_path}' is missing rows for resolved trait(s): {annotated}. "
            "Add a row per trait (extra rows for other traits are allowed), or use the positional "
            "`--samp-prev`/`--pop-prev` lists instead."
        )
    return [table[t] for t in trait_names]


def _reject_duplicate_names(trait_names: list[str], trait_paths: list[str]) -> None:
    """Raise if resolved munged trait names are not unique (ambiguous by-name lookup)."""
    seen: dict[str, list[str]] = {}
    for name, path in zip(trait_names, trait_paths):
        seen.setdefault(name, []).append(path)
    dupes = {name: paths for name, paths in seen.items() if len(paths) > 1}
    if dupes:
        detail = "; ".join(f"{name}: {', '.join(paths)}" for name, paths in dupes.items())
        raise LDSCUsageError(
            f"rg `--prevalence-manifest` cannot match by name because resolved traits have duplicate "
            f"munged names ({detail}). Either use the positional `--samp-prev`/`--pop-prev` lists "
            "(matched by position), or re-munge so each input has a distinct `trait_name` metadata value."
        )


def resolve_rg_prevalences(
    samp_prev, pop_prev, manifest_path, trait_names: list[str], trait_paths: list[str]
) -> list[TraitPrevalence] | None:
    """Resolve rg prevalence from exactly one supplied scheme, or None.

    Schemes are mutually exclusive: positional comma-lists (``samp_prev`` and
    ``pop_prev``) xor ``manifest_path``. Returns a per-trait list aligned to
    ``trait_names``, or ``None`` when no scheme is supplied (observed scale).
    """
    positional = samp_prev is not None or pop_prev is not None
    if positional and manifest_path is not None:
        raise LDSCUsageError(
            "rg received both `--prevalence-manifest` and `--samp-prev`/`--pop-prev`. "
            "These schemes are mutually exclusive; choose one."
        )
    if manifest_path is not None:
        _reject_duplicate_names(trait_names, trait_paths)
        return parse_manifest_prevalences(manifest_path, trait_names, trait_paths)
    if positional:
        return parse_positional_prevalences(samp_prev, pop_prev, len(trait_names))
    return None
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `PYTHONPATH=src pytest tests/test_prevalence.py -v`
Expected: PASS (all).

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/prevalence.py tests/test_prevalence.py
git commit -m "feat(prevalence): add scalar/positional/manifest prevalence parsing"
```

---

## Task 3: Both-scale columns in `summarize_total_h2` (h2 + h2_per_trait)

**Files:**
- Modify: `src/ldsc/regression_runner.py` (`summarize_total_h2` ~lines 1253-1282; add `_liability_pair` helper nearby; add import)
- Test: `tests/test_regression_workflow.py`, `tests/test_output.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_regression_workflow.py` (uses the existing fake-`hsq`/`dataset` fixtures in that file; mirror `test_summarize_total_h2_reports_post_filter_count` at line 1120 for fixture construction):

```python
def test_summarize_total_h2_observed_when_no_prevalence(self):
    hsq = _make_fake_hsq(tot=0.2, tot_se=0.03)        # reuse the file's fixture builder
    dataset = _make_fake_h2_dataset(n=10)
    summary = regression_runner.summarize_total_h2(hsq, dataset, trait_name="t", n_snps_used=10)
    row = summary.iloc[0]
    assert row["total_h2_obs"] == 0.2
    assert row["total_h2_obs_se"] == 0.03
    assert math.isnan(row["total_h2_liab"])
    assert math.isnan(row["total_h2_liab_se"])
    assert math.isnan(row["samp_prev"]) and math.isnan(row["pop_prev"])


def test_summarize_total_h2_liability_scaled(self):
    from ldsc._kernel.regression import liability_conversion_factor
    hsq = _make_fake_hsq(tot=0.2, tot_se=0.03)
    dataset = _make_fake_h2_dataset(n=10)
    summary = regression_runner.summarize_total_h2(
        hsq, dataset, trait_name="t", n_snps_used=10, samp_prev=0.5, pop_prev=0.01
    )
    c = liability_conversion_factor(0.5, 0.01)
    row = summary.iloc[0]
    assert math.isclose(row["total_h2_liab"], 0.2 * c, rel_tol=1e-12)
    assert math.isclose(row["total_h2_liab_se"], 0.03 * c, rel_tol=1e-12)
    assert row["total_h2_obs"] == 0.2
    assert row["samp_prev"] == 0.5 and row["pop_prev"] == 0.01
```

> If `tests/test_regression_workflow.py` lacks reusable `_make_fake_hsq`/`_make_fake_h2_dataset` builders, instead extend the existing test at line 1120 (`test_summarize_total_h2_reports_post_filter_count`) which already constructs a usable `hsq`/`dataset`; copy its setup into the two new tests.

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "summarize_total_h2" -v`
Expected: FAIL -- `KeyError: 'total_h2_obs'` (current column is `total_h2`).

- [ ] **Step 3: Implement**

In `src/ldsc/regression_runner.py`, add near the other small helpers (e.g. just below `_scalar`), importing the primitive at the top of the file (`from ._kernel.regression import liability_conversion_factor, gencov_obs_to_liab`):

```python
def _liability_pair(obs, obs_se, samp_prev, pop_prev):
    """Return (liab, liab_se) for an absolute h2 estimate, or (nan, nan).

    Liability is defined only when both prevalences are finite; a quantitative
    trait (None/NaN) or an unset prevalence yields NaN. ``obs``/``obs_se`` are
    observed-scale point estimate and SE.
    """
    if samp_prev is None or pop_prev is None:
        return float("nan"), float("nan")
    if isinstance(samp_prev, float) and math.isnan(samp_prev):
        return float("nan"), float("nan")
    if isinstance(pop_prev, float) and math.isnan(pop_prev):
        return float("nan"), float("nan")
    c = float(liability_conversion_factor(samp_prev, pop_prev))
    return obs * c, obs_se * c


def _prev_cell(value):
    """Normalize a prevalence for a table cell: None -> NaN, else float."""
    return float("nan") if value is None else float(value)
```

Update `summarize_total_h2` signature and body (replace lines 1253-1282):

```python
def summarize_total_h2(
    hsq,
    dataset: RegressionDataset,
    trait_name: str | None = None,
    n_snps_used: int | None = None,
    samp_prev: float | None = None,
    pop_prev: float | None = None,
) -> pd.DataFrame:
    """Build the one-row total-heritability summary from a fitted ``Hsq`` result.

    Reports total heritability on both the observed (``total_h2_obs``) and
    liability (``total_h2_liab``) scales. Liability columns are NaN unless both
    ``samp_prev`` and ``pop_prev`` are finite probabilities in (0, 1). ``n_snps_used``
    is the post-chi-square-filter SNP count; when omitted it falls back to the
    full merged count.
    """
    obs = _scalar(hsq.tot)
    obs_se = _scalar(hsq.tot_se)
    liab, liab_se = _liability_pair(obs, obs_se, samp_prev, pop_prev)
    return pd.DataFrame(
        [
            {
                "trait_name": trait_name,
                "n_snps": int(n_snps_used) if n_snps_used is not None else len(dataset.merged),
                "total_h2_obs": obs,
                "total_h2_obs_se": obs_se,
                "total_h2_liab": liab,
                "total_h2_liab_se": liab_se,
                "samp_prev": _prev_cell(samp_prev),
                "pop_prev": _prev_cell(pop_prev),
                "intercept": _scalar_or_value(hsq.intercept),
                "intercept_se": getattr(hsq, "intercept_se", None),
                "mean_chisq": _scalar_or_value(hsq.mean_chisq),
                "lambda_gc": _scalar_or_value(hsq.lambda_gc),
                "ratio": getattr(hsq, "ratio", None),
                "ratio_se": getattr(hsq, "ratio_se", None),
            }
        ]
    )
```

Update `tests/test_output.py` h2 fixtures/assertions that use the old names (lines 447-448, 530, 864-866): change `"total_h2"`->`"total_h2_obs"`, `"total_h2_se"`->`"total_h2_obs_se"`, and the `assertIn("total_h2", ...)` to `assertIn("total_h2_obs", ...)`.

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "summarize_total_h2" tests/test_output.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/regression_runner.py tests/test_regression_workflow.py tests/test_output.py
git commit -m "feat(h2): report observed and liability total h2 with prevalence columns"
```

---

## Task 4: Both-scale columns in `summarize_partitioned_h2`

**Files:**
- Modify: `src/ldsc/outputs.py` (`PARTITIONED_H2_COLUMNS`, lines 102-117)
- Modify: `src/ldsc/regression_runner.py` (`summarize_partitioned_h2`, lines 1285-1314; callers pass `samp_prev`/`pop_prev`)
- Test: `tests/test_regression_workflow.py`, `tests/test_output.py`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_regression_workflow.py` (reuse the partitioned fixture used by the test at line 1797 that asserts `PARTITIONED_H2_COLUMNS`):

```python
def test_summarize_partitioned_h2_both_scales(self):
    from ldsc._kernel.regression import liability_conversion_factor
    hsq, dataset, cols = _make_partitioned_fixture()      # reuse existing fixture path
    obs = regression_runner.summarize_partitioned_h2(hsq, dataset, cols)
    liab = regression_runner.summarize_partitioned_h2(hsq, dataset, cols, samp_prev=0.5, pop_prev=0.01)
    c = liability_conversion_factor(0.5, 0.01)
    # Observed column preserved; liability is c * observed; proportions invariant.
    assert (liab["Category_h2_obs"] == obs["Category_h2_obs"]).all()
    assert math.isclose(liab["Category_h2_liab"].iloc[0], obs["Category_h2_obs"].iloc[0] * c, rel_tol=1e-12)
    assert (liab["Prop._h2"] == obs["Prop._h2"]).all()
    assert math.isnan(obs["Category_h2_liab"].iloc[0])
    assert (liab["samp_prev"] == 0.5).all()
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "partitioned_h2_both_scales" -v`
Expected: FAIL -- `TypeError: summarize_partitioned_h2() got an unexpected keyword argument 'samp_prev'`.

- [ ] **Step 3: Implement**

Replace `PARTITIONED_H2_COLUMNS` in `src/ldsc/outputs.py` (lines 102-117):

```python
PARTITIONED_H2_COLUMNS = [
    "Category",
    "Prop._SNPs",
    "Category_h2_obs",
    "Category_h2_obs_std_error",
    "Category_h2_liab",
    "Category_h2_liab_std_error",
    "samp_prev",
    "pop_prev",
    "Prop._h2",
    "Prop._h2_std_error",
    "Enrichment",
    "Enrichment_std_error",
    "Enrichment_p",
    "Coefficient",
    "Coefficient_std_error",
    "Coefficient_z",
    "Coefficient_p",
    "overlap_aware",
]
```

Update `summarize_partitioned_h2` in `src/ldsc/regression_runner.py` (signature + body, lines 1285-1314):

```python
def summarize_partitioned_h2(
    hsq,
    dataset: RegressionDataset,
    annotation_columns: Sequence[str],
    samp_prev: float | None = None,
    pop_prev: float | None = None,
) -> pd.DataFrame:
    """Build overlap-aware partitioned-h2 rows from a fitted ``Hsq`` result.

    Per-category absolute heritability is reported on both the observed
    (``Category_h2_obs``) and liability (``Category_h2_liab``) scales; liability
    columns are NaN unless both prevalences are finite probabilities in (0, 1).
    Proportions, enrichment, and coefficients are scale-invariant and unchanged.
    ``annotation_columns`` must be a subset of ``dataset.retained_ld_columns``.
    """
    from .overlap_matrix import assemble_model_overlap, overlap_aware_category_table

    overlap = dataset.ldscore_overlap
    if overlap is None:
        raise LDSCInputError(PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE)
    use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY
    model_overlap = assemble_model_overlap(overlap, list(dataset.retained_ld_columns), use_common)
    m_annot = np.asarray(
        dataset.reference_snp_count_totals[dataset.count_key_used_for_regression], dtype=np.float64
    )
    m_tot = overlap.total_common_reference_snps if use_common else overlap.total_all_reference_snps
    table = overlap_aware_category_table(
        hsq, model_overlap, m_annot, float(m_tot), dataset.retained_ld_columns
    )
    table = table.rename(
        columns={
            "Category_h2": "Category_h2_obs",
            "Category_h2_std_error": "Category_h2_obs_std_error",
        }
    )
    c = (
        float(liability_conversion_factor(samp_prev, pop_prev))
        if samp_prev is not None and pop_prev is not None
        and not (isinstance(samp_prev, float) and math.isnan(samp_prev))
        and not (isinstance(pop_prev, float) and math.isnan(pop_prev))
        else float("nan")
    )
    table["Category_h2_liab"] = table["Category_h2_obs"] * c
    table["Category_h2_liab_std_error"] = table["Category_h2_obs_std_error"] * c
    table["samp_prev"] = _prev_cell(samp_prev)
    table["pop_prev"] = _prev_cell(pop_prev)
    rows = table.set_index("Category").loc[list(annotation_columns)].reset_index()
    return rows.loc[:, PARTITIONED_H2_COLUMNS].reset_index(drop=True)
```

Update `tests/test_output.py` partitioned fixtures (lines 559-560, 580-581, 598-599) to use `Category_h2_obs`/`Category_h2_obs_std_error` and add the new columns where full `PARTITIONED_H2_COLUMNS` frames are constructed (so `_select_columns` finds them). For frames that the writer selects via `PARTITIONED_H2_COLUMNS`, add `Category_h2_liab`, `Category_h2_liab_std_error`, `samp_prev`, `pop_prev` keys (values may be NaN) to the fixture dicts.

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "partitioned" tests/test_output.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/outputs.py src/ldsc/regression_runner.py tests/test_regression_workflow.py tests/test_output.py
git commit -m "feat(partitioned-h2): report observed and liability category h2"
```

---

## Task 5: Both-scale columns in rg (`_summarize_rg_pair`, `RG_FULL_COLUMNS`, threading)

**Files:**
- Modify: `src/ldsc/outputs.py` (`RG_FULL_COLUMNS`, lines 127-158)
- Modify: `src/ldsc/regression_runner.py` (`_summarize_rg_pair` ~1320-1374; `estimate_rg_pairs` ~878-968; `_concise_rg_row` unaffected)
- Test: `tests/test_regression_workflow.py`, `tests/test_output.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_regression_workflow.py` (reuse the rg fixture that builds a fake `rg_result` with `hsq1`/`hsq2`/`gencov`; mirror existing rg-pair tests in the file):

```python
def test_summarize_rg_pair_liability_and_invariant_ratio(self):
    from ldsc._kernel.regression import liability_conversion_factor, gencov_obs_to_liab
    rg_result, dataset = _make_fake_rg_result(h2_1=0.2, h2_2=0.3, gencov=0.1, rg=0.4, rg_se=0.05)
    obs_row = regression_runner._summarize_rg_pair(
        rg_result, dataset, trait_1="A", trait_2="B", pair_kind="all_pairs"
    )
    liab_row = regression_runner._summarize_rg_pair(
        rg_result, dataset, trait_1="A", trait_2="B", pair_kind="all_pairs",
        samp_prev_1=0.5, pop_prev_1=0.01, samp_prev_2=0.5, pop_prev_2=0.02,
    )
    c1 = liability_conversion_factor(0.5, 0.01)
    c2 = liability_conversion_factor(0.5, 0.02)
    # rg ratio is scale-invariant.
    assert liab_row["rg"] == obs_row["rg"]
    assert liab_row["rg_se"] == obs_row["rg_se"]
    # Per-trait h2 and gencov go liability.
    assert math.isclose(liab_row["h2_1_liab"], 0.2 * c1, rel_tol=1e-12)
    assert math.isclose(liab_row["h2_2_liab"], 0.3 * c2, rel_tol=1e-12)
    assert math.isclose(liab_row["gencov_liab"], gencov_obs_to_liab(0.1, 0.5, 0.5, 0.01, 0.02), rel_tol=1e-12)
    assert liab_row["h2_1_obs"] == 0.2
    assert math.isnan(obs_row["h2_1_liab"])
    assert liab_row["samp_prev_1"] == 0.5 and liab_row["pop_prev_2"] == 0.02
```

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "summarize_rg_pair_liability" -v`
Expected: FAIL -- `TypeError: _summarize_rg_pair() got an unexpected keyword argument 'samp_prev_1'`.

- [ ] **Step 3: Implement**

Replace `RG_FULL_COLUMNS` in `src/ldsc/outputs.py` (lines 127-158) with (rename `h2_1`/`h2_2`/`gencov` and their `_se` to `_obs`, insert `_liab` + 4 prevalence columns after the gencov block):

```python
RG_FULL_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "z",
    "p",
    "h2_1_obs",
    "h2_1_obs_se",
    "h2_1_liab",
    "h2_1_liab_se",
    "h2_2_obs",
    "h2_2_obs_se",
    "h2_2_liab",
    "h2_2_liab_se",
    "gencov_obs",
    "gencov_obs_se",
    "gencov_liab",
    "gencov_liab_se",
    "samp_prev_1",
    "pop_prev_1",
    "samp_prev_2",
    "pop_prev_2",
    "intercept_h2_1",
    "intercept_h2_1_se",
    "intercept_h2_2",
    "intercept_h2_2_se",
    "intercept_gencov",
    "intercept_gencov_se",
    "ratio_1",
    "ratio_1_se",
    "ratio_2",
    "ratio_2_se",
    "lambda_gc_1",
    "lambda_gc_2",
    "mean_chisq_1",
    "mean_chisq_2",
    "pair_kind",
    "status",
    "error",
]
```

Update `_summarize_rg_pair` in `src/ldsc/regression_runner.py` (signature adds four prevalence params; rename the three obs keys; add liab + prevalence keys). Change the `h2_1`/`h2_1_se`/`h2_2`/`h2_2_se`/`gencov`/`gencov_se` dict entries to `_obs` names and insert after them:

```python
def _summarize_rg_pair(
    rg_result,
    dataset: RGRegressionDataset,
    *,
    trait_1: str,
    trait_2: str,
    pair_kind: str,
    n_snps_used: int | None = None,
    samp_prev_1: float | None = None,
    pop_prev_1: float | None = None,
    samp_prev_2: float | None = None,
    pop_prev_2: float | None = None,
) -> dict[str, object]:
    """Build the full public rg row from one fitted kernel result.

    Per-trait heritability and the pair's genetic covariance are reported on both
    the observed (``*_obs``) and liability (``*_liab``) scales. The rg ratio
    (``rg``/``rg_se``/``z``/``p``) is scale-invariant and unconverted. Liability
    columns are NaN for traits without finite prevalence; ``gencov_liab`` is NaN
    only when no prevalence was supplied for the run.
    """
    h2_1_obs = _numeric_attr(getattr(rg_result, "hsq1", None), "tot", "h2_1_obs")
    h2_1_obs_se = _numeric_attr(getattr(rg_result, "hsq1", None), "tot_se", "h2_1_obs_se")
    h2_2_obs = _numeric_attr(getattr(rg_result, "hsq2", None), "tot", "h2_2_obs")
    h2_2_obs_se = _numeric_attr(getattr(rg_result, "hsq2", None), "tot_se", "h2_2_obs_se")
    gencov_obs = _numeric_attr(getattr(rg_result, "gencov", None), "tot", "gencov_obs")
    gencov_obs_se = _numeric_attr(getattr(rg_result, "gencov", None), "tot_se", "gencov_obs_se")
    h2_1_liab, h2_1_liab_se = _liability_pair(h2_1_obs, h2_1_obs_se, samp_prev_1, pop_prev_1)
    h2_2_liab, h2_2_liab_se = _liability_pair(h2_2_obs, h2_2_obs_se, samp_prev_2, pop_prev_2)
    any_prev = any(v is not None for v in (samp_prev_1, pop_prev_1, samp_prev_2, pop_prev_2))
    if any_prev:
        gfac = gencov_obs_to_liab(1.0, samp_prev_1, samp_prev_2, pop_prev_1, pop_prev_2)
        gencov_liab, gencov_liab_se = gencov_obs * gfac, gencov_obs_se * gfac
    else:
        gencov_liab, gencov_liab_se = float("nan"), float("nan")
    return {
        "trait_1": trait_1,
        "trait_2": trait_2,
        "n_snps_used": int(n_snps_used) if n_snps_used is not None else int(len(dataset.merged)),
        "rg": _required_numeric_scalar(getattr(rg_result, "rg_ratio", None), "rg"),
        "rg_se": _required_numeric_scalar(getattr(rg_result, "rg_se", None), "rg_se"),
        "z": _required_numeric_scalar(getattr(rg_result, "z", None), "z"),
        "p": _required_numeric_scalar(getattr(rg_result, "p", None), "p"),
        "h2_1_obs": h2_1_obs,
        "h2_1_obs_se": h2_1_obs_se,
        "h2_1_liab": h2_1_liab,
        "h2_1_liab_se": h2_1_liab_se,
        "h2_2_obs": h2_2_obs,
        "h2_2_obs_se": h2_2_obs_se,
        "h2_2_liab": h2_2_liab,
        "h2_2_liab_se": h2_2_liab_se,
        "gencov_obs": gencov_obs,
        "gencov_obs_se": gencov_obs_se,
        "gencov_liab": gencov_liab,
        "gencov_liab_se": gencov_liab_se,
        "samp_prev_1": _prev_cell(samp_prev_1),
        "pop_prev_1": _prev_cell(pop_prev_1),
        "samp_prev_2": _prev_cell(samp_prev_2),
        "pop_prev_2": _prev_cell(pop_prev_2),
        "intercept_h2_1": _numeric_attr(getattr(rg_result, "hsq1", None), "intercept", "intercept_h2_1"),
        "intercept_h2_1_se": _numeric_attr(getattr(rg_result, "hsq1", None), "intercept_se", "intercept_h2_1_se"),
        "intercept_h2_2": _numeric_attr(getattr(rg_result, "hsq2", None), "intercept", "intercept_h2_2"),
        "intercept_h2_2_se": _numeric_attr(getattr(rg_result, "hsq2", None), "intercept_se", "intercept_h2_2_se"),
        "intercept_gencov": _numeric_attr(getattr(rg_result, "gencov", None), "intercept", "intercept_gencov"),
        "intercept_gencov_se": _numeric_attr(getattr(rg_result, "gencov", None), "intercept_se", "intercept_gencov_se"),
        "ratio_1": _numeric_attr(getattr(rg_result, "hsq1", None), "ratio", "ratio_1"),
        "ratio_1_se": _numeric_attr(getattr(rg_result, "hsq1", None), "ratio_se", "ratio_1_se"),
        "ratio_2": _numeric_attr(getattr(rg_result, "hsq2", None), "ratio", "ratio_2"),
        "ratio_2_se": _numeric_attr(getattr(rg_result, "hsq2", None), "ratio_se", "ratio_2_se"),
        "lambda_gc_1": _numeric_attr(getattr(rg_result, "hsq1", None), "lambda_gc", "lambda_gc_1"),
        "lambda_gc_2": _numeric_attr(getattr(rg_result, "hsq2", None), "lambda_gc", "lambda_gc_2"),
        "mean_chisq_1": _numeric_attr(getattr(rg_result, "hsq1", None), "mean_chisq", "mean_chisq_1"),
        "mean_chisq_2": _numeric_attr(getattr(rg_result, "hsq2", None), "mean_chisq", "mean_chisq_2"),
        "pair_kind": pair_kind,
        "status": "ok",
        "error": "",
    }
```

Thread per-trait prevalence through `estimate_rg_pairs` (lines 878-968). Add `prevalences: list[tuple[float | None, float | None]] | None = None` to the signature; in the per-trait h2 loop pass that trait's pair into `summarize_total_h2`; in the pair loop pass both traits' pairs into `_summarize_rg_pair`:

```python
    def estimate_rg_pairs(
        self,
        sumstats_tables,
        ldscore_result,
        anchor_index: int | None = None,
        config: RegressionConfig | None = None,
        prevalences: list[tuple[float | None, float | None]] | None = None,
    ) -> RgResultFamily:
        ...
        config = config or self.regression_config
        tables = list(sumstats_tables)
        ...
        prev = prevalences if prevalences is not None else [(None, None)] * len(tables)
        h2_rows = []
        for table, (sp, pp) in zip(tables, prev):
            h2_dataset = self.build_dataset(table, ldscore_result, config=config)
            hsq = self.estimate_h2(h2_dataset, config=config)
            _, n_snps_used = _effective_regression_filter(h2_dataset, config)
            h2_rows.append(
                summarize_total_h2(
                    hsq, h2_dataset, trait_name=_trait_label(table),
                    n_snps_used=n_snps_used, samp_prev=sp, pop_prev=pp,
                )
            )
        ...
        for i, j in _iter_rg_pairs(len(tables), anchor_index):
            ...
                full_row = _summarize_rg_pair(
                    fitted, dataset, trait_1=trait_1, trait_2=trait_2, pair_kind=pair_kind,
                    n_snps_used=n_snps_used,
                    samp_prev_1=prev[i][0], pop_prev_1=prev[i][1],
                    samp_prev_2=prev[j][0], pop_prev_2=prev[j][1],
                )
```

(Keep the existing `_effective_regression_filter` / `_effective_rg_filter` calls intact -- only add the prevalence arguments.)

Update `tests/test_output.py` rg fixtures (lines 804-845): rename `h2_1`/`h2_1_se`/`h2_2`/`h2_2_se`/`gencov`/`gencov_se` to the `_obs` names and add the `_liab` + four `*_prev_*` keys (NaN allowed) so frames match `RG_FULL_COLUMNS`. The `_failed_rg_full_row` helper already fills every `RG_FULL_COLUMNS` key with NaN, so failed rows need no change.

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "rg" tests/test_output.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/outputs.py src/ldsc/regression_runner.py tests/test_regression_workflow.py tests/test_output.py
git commit -m "feat(rg): report observed and liability per-trait h2 and gencov"
```

---

## Task 6: CLI flags + runner wiring + metadata + logging

**Files:**
- Modify: `src/ldsc/regression_runner.py` (`add_h2_arguments` 1746, `add_partitioned_h2_arguments` 1753, `add_rg_arguments` 1776; `_runner_from_args` 2195; `run_h2_from_args` 1803; `run_partitioned_h2_from_args` 1852; `run_rg_from_args` 1941; `_h2_metadata` ~1404; `_rg_pair_metadata` ~1436)
- Modify: `src/ldsc/config.py` (RegressionConfig docstring lines 994-995, valid range)
- Test: `tests/test_regression_workflow.py`

- [ ] **Step 1: Write the failing tests**

Add to `tests/test_regression_workflow.py` (these exercise argument plumbing end-to-end against the existing small fixture LD-score + sumstats used by other workflow tests in the file; follow the existing `run_h2_from_args`/`run_rg_from_args` test patterns):

```python
def test_h2_cli_liability_scale(self):
    args = _make_h2_args(samp_prev=0.5, pop_prev=0.01)   # mirror existing _make_h2_args helper
    summary = regression_runner.run_h2_from_args(args)
    row = summary.iloc[0]
    assert not math.isnan(row["total_h2_liab"])
    assert row["samp_prev"] == 0.5 and row["pop_prev"] == 0.01


def test_h2_cli_half_specified_raises(self):
    args = _make_h2_args(samp_prev=0.5, pop_prev=None)
    with self.assertRaises((LDSCUsageError, LDSCConfigError)):
        regression_runner.run_h2_from_args(args)


def test_rg_cli_positional_prevalences(self):
    args = _make_rg_args(samp_prev="0.5,0.5", pop_prev="0.01,0.02")  # 2 resolved traits
    result = regression_runner.run_rg_from_args(args)
    full = result.rg_full.iloc[0]
    assert not math.isnan(full["h2_1_liab"])
    assert not math.isnan(full["gencov_liab"])
    # rg ratio unchanged vs an observed-only run.
    obs = regression_runner.run_rg_from_args(_make_rg_args(samp_prev=None, pop_prev=None))
    assert math.isclose(result.rg.iloc[0]["rg"], obs.rg.iloc[0]["rg"], rel_tol=1e-12)
```

> If `_make_h2_args` / `_make_rg_args` helpers do not already exist in the file, build the `argparse.Namespace` inline copying the fields the existing workflow tests set, adding `samp_prev`, `pop_prev`, and (for rg) `prevalence_manifest`.

- [ ] **Step 2: Run to verify failure**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -k "cli_liability or cli_positional or cli_half_specified" -v`
Expected: FAIL -- `AttributeError`/`KeyError` (flags/columns not yet wired).

- [ ] **Step 3: Implement**

Add the scalar flags in `add_h2_arguments` and `add_partitioned_h2_arguments` (after their `--trait-name` lines):

```python
    parser.add_argument(
        "--samp-prev",
        type=float,
        default=None,
        help="Sample (case) prevalence for liability-scale conversion: a probability in (0, 1). "
        "Requires --pop-prev. Omit both for observed scale.",
    )
    parser.add_argument(
        "--pop-prev",
        type=float,
        default=None,
        help="Population prevalence for liability-scale conversion: a probability in (0, 1). "
        "Requires --samp-prev. Omit both for observed scale.",
    )
```

Add the rg flags in `add_rg_arguments` (after `--anchor-trait`):

```python
    parser.add_argument(
        "--samp-prev",
        default=None,
        help="Comma-separated sample prevalences aligned to resolved --sumstats-sources order, one per "
        "trait; each a probability in (0, 1) or `nan` for a quantitative trait. Mutually exclusive with "
        "--prevalence-manifest.",
    )
    parser.add_argument(
        "--pop-prev",
        default=None,
        help="Comma-separated population prevalences aligned to resolved --sumstats-sources order, one per "
        "trait; each a probability in (0, 1) or `nan`. Mutually exclusive with --prevalence-manifest.",
    )
    parser.add_argument(
        "--prevalence-manifest",
        default=None,
        help="Whitespace/tab-delimited TSV with columns trait_name, samp_prev, pop_prev (`#` comments "
        "ignored). Looked up by exact munged trait name; may contain extra traits. Mutually exclusive "
        "with --samp-prev/--pop-prev.",
    )
```

In `_runner_from_args` (lines 2195-2210), set scalar prevalence on the config **only** for the scalar (float) case (rg passes strings, handled separately). Insert before `config = RegressionConfig(...)`:

```python
    from .prevalence import parse_scalar_prevalence

    samp_prev = getattr(args, "samp_prev", None)
    pop_prev = getattr(args, "pop_prev", None)
    if isinstance(samp_prev, str) or isinstance(pop_prev, str):
        scalar_samp, scalar_pop = None, None  # rg comma-lists: resolved in run_rg_from_args
    else:
        scalar_samp, scalar_pop = parse_scalar_prevalence(samp_prev, pop_prev)
```

and pass `samp_prev=scalar_samp, pop_prev=scalar_pop` into the `RegressionConfig(...)` constructor.

In `run_h2_from_args` (line 1828), pass the config prevalence into the summary:

```python
        summary = summarize_total_h2(
            hsq, dataset, trait_name=sumstats_table.trait_name, n_snps_used=n_snps_used,
            samp_prev=config.samp_prev, pop_prev=config.pop_prev,
        )
```

In `run_partitioned_h2_from_args`, locate each `summarize_partitioned_h2(hsq, dataset, ...)` call (in `RegressionRunner.estimate_partitioned_h2_batch`, lines ~764 and ~783) and pass `samp_prev=config.samp_prev, pop_prev=config.pop_prev`.

In `run_rg_from_args` (after `sumstats_tables = _disambiguate_trait_names(...)`, line ~2013), resolve and thread prevalences, and log them:

```python
        from .prevalence import resolve_rg_prevalences

        trait_names = [t.trait_name for t in sumstats_tables]
        trait_paths = [str(p) for p in sumstats_paths]
        prevalences = resolve_rg_prevalences(
            samp_prev=getattr(args, "samp_prev", None),
            pop_prev=getattr(args, "pop_prev", None),
            manifest_path=getattr(args, "prevalence_manifest", None),
            trait_names=trait_names,
            trait_paths=trait_paths,
        )
        if prevalences is not None:
            LOGGER.info(
                "Liability-scale prevalences applied: "
                + ", ".join(
                    f"{name}=(P={sp}, K={pp})" for name, (sp, pp) in zip(trait_names, prevalences)
                )
            )
        ...
            result = runner.estimate_rg_pairs(
                sumstats_tables, ldscore_result, anchor_index=anchor_index, config=config,
                prevalences=prevalences,
            )
```

Record prevalence in metadata. In `_h2_metadata` (lines ~1404-1435) add to the returned dict:

```python
        "samp_prev": config_snapshot.samp_prev if hasattr(config_snapshot, "samp_prev") else None,
        "pop_prev": config_snapshot.pop_prev if hasattr(config_snapshot, "pop_prev") else None,
        "scale": "liability" if getattr(config_snapshot, "pop_prev", None) is not None else "observed",
```

(If `config_snapshot` does not expose these, read them from `args` instead: `getattr(args, "samp_prev", None)` / `getattr(args, "pop_prev", None)`.) For `_rg_pair_metadata`, add the pair's `samp_prev_1/pop_prev_1/samp_prev_2/pop_prev_2` from the resolved `prevalences` (thread the resolved list into the metadata builder analogously to the summary).

Update `src/ldsc/config.py` RegressionConfig docstring (lines 994-995):

```python
    samp_prev, pop_prev : float, list of float, or None, optional
        Liability-scale prevalence inputs for binary traits: each a probability
        in the open interval ``(0, 1)``, or NaN/None for a quantitative trait.
        For ``h2``/``partitioned-h2`` these are scalars; for ``rg`` the resolved
        per-trait list is passed to ``estimate_rg_pairs`` rather than stored here.
        Both must be supplied together. Defaults are ``None`` (observed scale).
```

- [ ] **Step 4: Run to verify pass**

Run: `PYTHONPATH=src pytest tests/test_regression_workflow.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/regression_runner.py src/ldsc/config.py tests/test_regression_workflow.py
git commit -m "feat(regression): wire --samp-prev/--pop-prev/--prevalence-manifest CLI flags"
```

---

## Task 7: Documentation

**Files:**
- Modify: `docs/current/regression-configuration.md` (replace the 1.4 "reserved" row at line 111)
- Modify: `docs/current/io-argument-inventory.md` (add flag rows for h2, partitioned-h2, rg)

- [ ] **Step 1: Update `regression-configuration.md`**

Replace the reserved row at line 111 with real behavior. New text (plain-text markers, no emojis per repo rule):

```markdown
| `samp_prev`, `pop_prev` | Liability-scale conversion inputs for binary (case-control) traits. `samp_prev` = case fraction in the GWAS sample; `pop_prev` = disease prevalence in the population. Each a probability in the open interval (0, 1), or `nan` for a quantitative trait. Required only for binary traits; both must be supplied together. | Wired via `--samp-prev` / `--pop-prev` (scalars for `h2` / `partitioned-h2`; comma-separated per-trait lists for `rg`, or an `rg` `--prevalence-manifest` TSV). Output tables report both `*_obs` and `*_liab` columns plus the applied prevalences; the rg ratio is scale-invariant and unchanged. Default (no prevalence) is observed scale, with `*_liab` columns NaN. |
```

Add a short subsection below the table documenting the rg manifest format (TSV columns `trait_name`, `samp_prev`, `pop_prev`; flexible whitespace/tab delimiter; `#` comments ignored; matched by exact trait name; superset allowed; missing resolved trait or duplicate resolved names are errors).

- [ ] **Step 2: Update `io-argument-inventory.md`**

Add rows for `--samp-prev` and `--pop-prev` under the `h2`, `partitioned-h2`, and `rg` subcommands, and `--prevalence-manifest` under `rg`, matching the inventory's existing table format (flag, type, default, description, valid range "(0, 1) or nan").

- [ ] **Step 3: Verify docs render / no broken references**

Run: `PYTHONPATH=src pytest -q` (full suite, confirms nothing references stale docs)
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add docs/current/regression-configuration.md docs/current/io-argument-inventory.md
git commit -m "docs(regression): document liability-scale prevalence flags and manifest"
```

---

## Task 8: User-facing prevalence how-to (tutorial)

**Files:**
- Create: `tutorials/liability-scale-conversion.md`

- [ ] **Step 1: Write the tutorial**

Create `tutorials/liability-scale-conversion.md`. Keep it concise and example-led (plain-text markers, no emojis per repo rule). Use this content:

````markdown
# Liability-Scale Heritability and Genetic Covariance

LDSC estimates heritability and genetic covariance on the **observed scale** by
default. For binary (case-control) traits, the observed-scale value depends on
the case fraction in your sample, so it is not comparable across studies. The
**liability scale** removes that dependence and is the standard scale to report.

You enable the conversion by giving LDSC two prevalences per binary trait:

- `--samp-prev` (P): the **case fraction in the GWAS sample** (e.g. 0.5 for a
  balanced case-control study).
- `--pop-prev` (K): the **disease prevalence in the population** (e.g. 0.01 for
  schizophrenia).

Each value is a probability in the open interval **(0, 1)**, or `nan` for a
**quantitative trait** (which has no liability scale). Supply both or neither.

The genetic correlation (`rg`) is scale-invariant -- it is identical on either
scale -- so conversion changes only heritability and genetic covariance. Output
tables always report both scales side by side (`*_obs` and `*_liab` columns) plus
the prevalences applied; `*_liab` columns are `NaN` when no prevalence is given.

## h2: one binary trait

```bash
ldsc h2 \
  --sumstats-file scz.sumstats.gz \
  --ldscore-dir eur_ldscores/ \
  --output-dir scz_h2/ \
  --samp-prev 0.5 \
  --pop-prev 0.01
```

`h2.tsv` reports `total_h2_obs` and `total_h2_liab` (= observed x conversion
factor), with `samp_prev` / `pop_prev` columns. Omit both flags for observed
scale only (`total_h2_liab` is then `NaN`).

## partitioned-h2: one binary trait

```bash
ldsc partitioned-h2 \
  --sumstats-file scz.sumstats.gz \
  --ldscore-dir baseline_ldscores/ \
  --output-dir scz_part/ \
  --samp-prev 0.5 \
  --pop-prev 0.01
```

`partitioned_h2.tsv` reports `Category_h2_obs` and `Category_h2_liab` per
category. Proportions, enrichment, and coefficients are scale-invariant.

## rg: two binary traits (ordered list)

Give a comma-separated value per trait, in the same order as
`--sumstats-sources`. Population and sample prevalence of schizophrenia and
bipolar disorder are about 1-2% and 50% here:

```bash
ldsc rg \
  --sumstats-sources scz.sumstats.gz bip.sumstats.gz \
  --ldscore-dir eur_ldscores/ \
  --output-dir scz_bip_rg/ \
  --samp-prev 0.5,0.5 \
  --pop-prev 0.01,0.02
```

`rg_full.tsv` reports `h2_1_liab`, `h2_2_liab`, and `gencov_liab`; `rg` is
unchanged. The run log echoes the prevalence applied to each trait.

## rg: one binary trait and one quantitative trait

Mark the quantitative trait with `nan` in both lists (its heritability has no
liability scale, but the genetic covariance is still partially rescaled):

```bash
ldsc rg \
  --sumstats-sources scz.sumstats.gz height.sumstats.gz \
  --ldscore-dir eur_ldscores/ \
  --samp-prev 0.5,nan \
  --pop-prev 0.01,nan
```

## rg: many traits via a prevalence manifest

For multi-trait runs (especially with glob inputs where the resolved order is
not obvious), provide a manifest looked up by **trait name** instead of position:

`prevalences.tsv`:

```text
# standing prevalence repository; nan = quantitative trait
trait_name   samp_prev   pop_prev
scz2022      0.43        0.01
bip2021      0.41        0.02
height       nan         nan
```

```bash
ldsc rg \
  --sumstats-sources 'sumstats/*.sumstats.gz' \
  --ldscore-dir eur_ldscores/ \
  --output-dir panel_rg/ \
  --prevalence-manifest prevalences.tsv
```

Manifest notes:

- Whitespace- or tab-delimited; header columns `trait_name`, `samp_prev`,
  `pop_prev` in any order; lines starting with `#` are ignored.
- Matched by **exact** munged trait name, so `scz2022` and `scz2014` never
  collide. The munged trait name is the `trait_name` recorded when you ran
  `ldsc munge-sumstats`.
- The manifest may list **more** traits than this run uses; extras are ignored.
  Every resolved trait must have a row, or the run aborts listing the missing
  name(s).
- `--prevalence-manifest` and `--samp-prev`/`--pop-prev` are mutually exclusive.

If two inputs share the same munged trait name, the manifest cannot tell them
apart; LDSC aborts and asks you to use the ordered list, or re-munge with
distinct trait names.
````

- [ ] **Step 2: Verify the example commands match the implemented flags**

Run: `PYTHONPATH=src ldsc rg --help` and `PYTHONPATH=src ldsc h2 --help`
Expected: `--samp-prev`, `--pop-prev`, and (rg) `--prevalence-manifest` appear with the (0, 1)-or-`nan` help text used in the tutorial.

- [ ] **Step 3: Commit**

```bash
git add tutorials/liability-scale-conversion.md
git commit -m "docs(tutorials): add liability-scale prevalence how-to with examples"
```

---

## Task 9: Full-suite verification

- [ ] **Step 1: Run the entire suite**

Run: `PYTHONPATH=src pytest -q`
Expected: PASS (all green, pristine output, no warnings introduced by this work).

- [ ] **Step 2: Spot-check a real liability run end-to-end (optional manual)**

Run an `h2` and a 2-trait `rg` against the repo's fixture LD-score directory with `--samp-prev 0.5 --pop-prev 0.01` (h2) and `--samp-prev 0.5,0.5 --pop-prev 0.01,0.02` (rg); confirm `total_h2_liab` / `h2_1_liab` / `gencov_liab` are populated, `rg` matches the observed-only run, and the log echoes the applied prevalences.

- [ ] **Step 3: Final commit if any fixups were needed**

```bash
git add -p
git commit -m "test(regression): finalize liability-scale conversion suite"
```

---

## Self-Review Notes (coverage map)

- Vectorizable primitive + (0,1)/NaN validation + both-NaN no-op -> Task 1.
- Scalar + positional + manifest parsing, mutual exclusivity, superset, duplicate-name guard, fail-fast -> Task 2.
- `h2` both-scale columns -> Task 3; `partitioned-h2` -> Task 4; `rg` (per-trait h2 + gencov, invariant ratio) -> Task 5.
- CLI flags for all three modules, runner/config wiring, metadata + log provenance, valid-range docstrings -> Task 6.
- Reference docs (`regression-configuration.md`, `io-argument-inventory.md`) -> Task 7.
- Separate user-facing how-to with worked examples for every prevalence input form -> Task 8.
- TDD numeric cross-check (`c(0.5,0.01)~=0.5519`), observed-unchanged-when-unset, rg ratio invariant, full suite green -> Tasks 1, 3, 5, 9.
- Known follow-up not in scope: K-sweep sensitivity output/plot module (primitive is vectorizable and ready).
