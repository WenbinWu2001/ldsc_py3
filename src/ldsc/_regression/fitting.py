"""Fit regression models and report their results.

This module holds the fit policy (``_FitPolicy``: intercepts, two-step
cutoff, and the default chi-square cap), the SNP heritability fit
(``_fit_h2_dataset``), the genetic-correlation fit (``_fit_rg_dataset``),
and the complete per-query partitioned fit that every query worker runs
(``fit_partitioned_query``). Both standalone fits are module-level
``serial_calculation`` callables, which the controlled child re-imports by
module and qualified name. Before an estimator runs, the collinearity guard
rejects a rank-deficient LD-score design; after a singular jackknife
deletion, the error gains the model's support counts and block spans.

It also turns fitted estimators into result rows: the h2 summary and its
regression-bin diagnostic, the partitioned-h2 summary, the rg pair rows and
the rows of failed fits, and the dataset, fit and run metadata recorded with
each result. Dataset assembly belongs to ``alignment``; batch and pair
orchestration to ``partitioned`` and ``rg``.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
import json
import logging
import math
import time
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .._runtime.paths import normalize_path_token
from .._runtime.serial_calculation import serial_calculation
from ..config import RegressionConfig
from ..errors import LDSCInputError, LDSCInternalError
from .._ldscore.results import LDScoreResult
from .._ldscore.source import LDScoreSource
from .._sumstats.table import SumstatsTable
from . import estimators as reg
from .jackknife import JackknifeIdentifiabilityError
from .overlap_model import assemble_model_overlap, overlap_aware_category_table, model_collinearity_error
from .outputs import H2_REGRESSION_BIN_COLUMNS, PARTITIONED_H2_COLUMNS, RG_FULL_COLUMNS, stage_partitioned_h2_fit
from .results import RegressionDataset, RegressionFit
from .alignment import COMMON_COUNT_KEY, _dataset_from_query_values, _trait_label


LOGGER = logging.getLogger("LDSC.regression_runner")
PARTITIONED_H2_REQUIRES_OVERLAP_MESSAGE = (
    "partitioned-h2 needs the annotation overlap matrix, but the LD-score directory has no "
    "ldscore.overlap.parquet. Most likely it is an unpartitioned LD-score directory (a single "
    "annotation such as the synthetic `base`), which has no overlap matrix and cannot be "
    "partitioned -- regenerate it with two or more annotation columns. It may also predate the "
    "overlap sidecar; if so, regenerate with the current `ldsc ldscore`. "
    "Other causes & fixes: docs/user-guide/troubleshooting.md#partitioned-h2-missing-overlap-matrix"
)
_RG_H2_OUT_OF_BOUNDS_DOC = "docs/user-guide/troubleshooting.md#rg-a-snp-heritability-estimate-is-not-positive"


@contextmanager
def _log_phase_timing(phase: str, *, level=logging.INFO):
    """Log elapsed wall time for one observable workflow phase."""
    started = time.perf_counter()
    LOGGER.log(level, "Phase start: %s.", phase)
    try:
        yield
    except Exception:
        LOGGER.log(level, "Phase timing: %s failed after %.3fs.", phase, time.perf_counter() - started)
        raise
    LOGGER.log(level, "Phase timing: %s completed in %.3fs.", phase, time.perf_counter() - started)


@serial_calculation
def _fit_h2_dataset(
    dataset: RegressionDataset,
    config: RegressionConfig,
    *, log_level=logging.INFO, retain_rows=True,
) -> RegressionFit:
    """Fit SNP heritability once and retain its exact SNP population and diagnostic bins.

    The aligned frame can borrow read-only worker maps. A multi-annotation
    design is filled column by column from the rows the chi-square filter
    keeps, so the frame itself is not copied before the fit; the filtered
    frame is built afterwards, for the returned dataset. With
    ``retain_rows=False`` (per-query partitioned fits, which read only the
    row count) the returned dataset's ``merged`` is an empty frame with one
    row per fitted SNP. The estimator never modifies the borrowed input columns.
    """
    merged = dataset.merged
    n_snp = len(merged)
    n_blocks = min(n_snp, config.num_blocks)
    chisq = np.asarray(merged["Z"] ** 2).reshape((n_snp, 1))
    # N.max() is read pre-filter so the partitioned default cap matches legacy.
    chisq_max = _resolve_default_chisq_max(
        config.chisq_max, len(dataset.retained_ld_columns), merged["N"].max()
    )
    keep = None
    if chisq_max is not None:
        mask = np.ravel(chisq <= chisq_max)
        n_removed = n_snp - int(mask.sum())
        if n_removed:
            keep = mask
            n_snp = int(mask.sum())
            n_blocks = min(n_snp, config.num_blocks)
            chisq = chisq[keep]
            LOGGER.log(log_level,
                "Removed %d SNPs with chi^2 > %g (%d SNPs remain).",
                n_removed,
                chisq_max,
                n_snp,
            )
    columns = dataset.retained_ld_columns
    if keep is not None and len(columns) == 1:
        # Single-annotation fits read the filtered frame for their diagnostic bins.
        merged, keep = _kept_rows(merged, keep), None
    x = _ld_design(merged, columns, keep)
    policy = _fit_policy(config)
    intercept, two_step = policy.intercept_h2, policy.two_step_cutoff(len(columns))
    old_weights = len(columns) > 1
    _raise_on_model_collinearity(dataset, x)
    reference_snp_counts = np.asarray(
        dataset.reference_snp_count_totals[dataset.count_key_used_for_regression]
    ).reshape((1, -1))
    try:
        hsq = reg.Hsq(
            chisq,
            x,
            np.asarray(_kept_rows(merged[[dataset.weight_column]], keep)),
            np.asarray(_kept_rows(merged[["N"]], keep)),
            reference_snp_counts,
            n_blocks=n_blocks,
            intercept=intercept,
            twostep=two_step,
            old_weights=old_weights,
        )
    except JackknifeIdentifiabilityError as error:
        # Two-step jackknives may use a different row population. The
        # partitioned model uses these exact post-filter rows throughout.
        if two_step is None:
            _add_jackknife_model_context(error, _kept_rows(merged, keep), columns)
        raise
    diagnostic_bins = None
    if len(columns) == 1:
        diagnostic_bins = summarize_ld_score_regression_bins(
            hsq,
            ld_score=x[:, 0],
            chi_square=chisq,
            sample_size=np.asarray(merged[["N"]]),
            regression_ld_score=np.asarray(merged[[dataset.weight_column]]),
            reference_snp_count=float(reference_snp_counts[0, 0]),
        )
    del x
    merged = _kept_rows(merged, keep) if retain_rows else pd.DataFrame(index=pd.RangeIndex(n_snp))
    return RegressionFit(
        estimator=hsq, dataset=replace(dataset, merged=merged), n_blocks=n_blocks,
        effective_chisq_max=chisq_max, diagnostic_bins=diagnostic_bins, effective_two_step_cutoff=two_step,
    )


@serial_calculation
def _fit_rg_dataset(dataset: RegressionDataset, config: RegressionConfig) -> RegressionFit:
    """Fit genetic correlation once on a prepared pair dataset; see ``RegressionRunner.estimate_rg``.

    Any error raised by the fit carries ``fit_summary``: the pair's SNP and
    block counts, the cap and the dataset settings.
    """
    merged = dataset.merged
    n_snp = len(merged)
    if config.chisq_max is not None:
        # Legacy rg filter (sumstats.py `_rg`): drop SNPs with
        # Z1^2 * Z2^2 > chisq_max^2, kept inclusive (<=) per the -max convention.
        keep = np.ravel((merged["Z1"] ** 2) * (merged["Z2"] ** 2) <= config.chisq_max ** 2)
        if not keep.all():
            merged = merged.take(np.flatnonzero(keep))
            merged.index = pd.RangeIndex(len(merged))
        n_snp = len(merged)
    n_blocks = min(n_snp, config.num_blocks)
    policy = _fit_policy(config)
    intercept_hsq, intercept_gencov = policy.intercept_h2, policy.intercept_gencov
    two_step = policy.two_step_cutoff(len(dataset.retained_ld_columns))
    x = np.asarray(merged[dataset.retained_ld_columns])
    try:
        _raise_on_model_collinearity(dataset, x)
        fitted = reg.RG(
            np.asarray(merged[["Z1"]]),
            np.asarray(merged[["Z2"]]),
            x,
            np.asarray(merged[[dataset.weight_column]]),
            np.asarray(merged[["N1"]]),
            np.asarray(merged[["N2"]]),
            np.asarray(dataset.reference_snp_count_totals[dataset.count_key_used_for_regression]).reshape((1, -1)),
            intercept_hsq1=intercept_hsq,
            intercept_hsq2=intercept_hsq,
            intercept_gencov=intercept_gencov,
            n_blocks=n_blocks,
            twostep=two_step,
        )
    except Exception as error:
        # As for h2, a two-step first stage fits a subset of these rows, so
        # spans and support counts are added only without it.
        if isinstance(error, JackknifeIdentifiabilityError) and two_step is None:
            _add_jackknife_model_context(error, merged, dataset.retained_ld_columns)
        # The pair's rows and settings, for failed-pair reporting.
        error.fit_summary = {"n_snps": n_snp, "n_blocks": n_blocks, "effective_chisq_max": config.chisq_max,
                             "effective_two_step_cutoff": two_step, **_dataset_metadata(dataset)}
        raise

    return RegressionFit(
        estimator=fitted, dataset=replace(dataset, merged=merged), n_blocks=n_blocks,
        effective_chisq_max=config.chisq_max, effective_two_step_cutoff=two_step,
    )


def fit_partitioned_query(prepared, query, values, directory, config):
    """Fit, summarize and privately stage one complete genome-wide query.

    Every query worker calls this module-level function, which spawn pickles
    by name; only model errors become query outcomes, while staging errors
    and interrupts still abort runs. Returned objects contain no estimator,
    SNP matrix or traceback frames.
    """
    stage = "model_preparation"
    try:
        with _log_phase_timing("query model preparation", level=logging.DEBUG):
            dataset = _dataset_from_query_values(prepared, [query], values, config)
        stage = "estimator"
        with _log_phase_timing("estimator execution", level=logging.DEBUG):
            outcome = _fit_h2_dataset(dataset, config, log_level=logging.DEBUG, retain_rows=False)
        hsq, dataset = outcome.estimator, outcome.dataset
        stage = "summary"
        with _log_phase_timing("query summary", level=logging.DEBUG):
            samp_prev, pop_prev = _config_prevalence(config)
            full = summarize_partitioned_h2(hsq, dataset, dataset.retained_ld_columns,
                                           samp_prev=samp_prev, pop_prev=pop_prev)
            # The summary computes every retained row; the focal row is one of them.
            focal = (full.set_index("category").loc[[query]].reset_index()
                     .loc[:, PARTITIONED_H2_COLUMNS].reset_index(drop=True))
            delete_values = _coefficient_delete_frame(hsq, dataset.retained_ld_columns)
    except Exception as error:
        status = "unestimable" if isinstance(error, JackknifeIdentifiabilityError) else "failed"
        record = {"query_annotation": query, "status": status, "stage": stage,
                  "error_type": type(error).__name__, "error_message": str(error)}
        LOGGER.exception("Query %r failed at %s; marked %s: %s", query, stage, status, error)
        if isinstance(error, JackknifeIdentifiabilityError):
            LOGGER.error("Query %r jackknife diagnostics: %s", query, json.dumps(error.failures))
        error.add_note(f"partitioned-h2 query {query!r} failed at {stage}.")
        error.__traceback__ = error.__cause__ = error.__context__ = None
        return None, None, record, error
    with _log_phase_timing("query output staging", level=logging.DEBUG):
        artifacts = stage_partitioned_h2_fit(
            directory, full, delete_values,
            {"dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
             **_fit_metadata(outcome),
             "effective_snp_identifier": dataset.effective_snp_identifier,
             "identity_downgrade_applied": dataset.identity_downgrade_applied},
        )
    return focal, artifacts, {"query_annotation": query, "status": "success", "stage": "complete",
                              "error_type": "", "error_message": ""}, None


def _kept_rows(frame: pd.DataFrame, keep: np.ndarray | None) -> pd.DataFrame:
    """Return the rows ``keep`` marks with a fresh RangeIndex, or ``frame`` itself when ``keep`` is None."""
    if keep is None:
        return frame
    kept = frame.loc[keep]
    kept.index = pd.RangeIndex(len(kept))
    return kept


def _ld_design(merged: pd.DataFrame, columns: Sequence[str], keep: np.ndarray | None) -> np.ndarray:
    """Return the LD-score design matrix of ``columns`` over the kept rows.

    A multi-column float design is filled column by column into the
    column-major array ``DataFrame.to_numpy()`` returns, without first copying
    the frame; other designs use ``to_numpy()`` on the kept rows.
    """
    columns = list(columns)
    dtypes = [merged[column].dtype for column in columns]
    if len(columns) < 2 or not all(isinstance(dtype, np.dtype) and dtype.kind == "f" for dtype in dtypes):
        return np.asarray(_kept_rows(merged, keep)[columns])
    n_rows = len(merged) if keep is None else int(keep.sum())
    design = np.empty((n_rows, len(columns)), dtype=np.result_type(*dtypes), order="F")
    for index, column in enumerate(columns):
        values = merged[column].to_numpy()
        # Assignment widens float32 columns exactly, as to_numpy() would.
        design[:, index] = values if keep is None else values[keep]
    return design


def _liability_pair(obs, obs_se, samp_prev, pop_prev):
    """Return (liab, liab_se) for an absolute SNP heritability estimate, or (nan, nan).

    The prevalences come checked from ``RegressionConfig`` or
    ``estimate_rg_pairs``, where a quantitative trait is ``(None, None)``;
    it yields NaN. ``obs`` and ``obs_se`` are the observed-scale point
    estimate and standard error.
    """
    if samp_prev is None or pop_prev is None:
        return float("nan"), float("nan")
    c = float(reg.liability_conversion_factor(samp_prev, pop_prev))
    return obs * c, obs_se * c


def _prev_cell(value):
    """Normalize a prevalence for a table cell: None -> NaN, else float."""
    return float("nan") if value is None else float(value)


def _config_prevalence(config: "RegressionConfig | None") -> tuple[float | None, float | None]:
    """Scalar (samp_prev, pop_prev) from a regression config, or (None, None)."""
    if config is None:
        return None, None
    return getattr(config, "samp_prev", None), getattr(config, "pop_prev", None)


def _scale_label(*prevalences) -> str:
    """Return 'liability' if any supplied prevalence is finite, else 'observed'."""
    return "liability" if any(p is not None for p in prevalences) else "observed"


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
    liability (``total_h2_liab``) scales, plus the prevalences applied. Liability
    columns are NaN unless both ``samp_prev`` and ``pop_prev`` are finite
    probabilities in (0, 1). ``intercept_se`` is NaN for a fixed intercept, and
    ``ratio``/``ratio_se`` are NaN when the mean chi-square is at most 1.
    ``n_snps_used`` is the post-chi-square-filter SNP
    count from the regression fit; when omitted it falls back
    to the supplied dataset's row count. Workflows supply the fitted dataset.
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
                "intercept_se": _numeric_scalar(getattr(hsq, "intercept_se", None), "intercept_se"),
                "mean_chisq": _scalar_or_value(hsq.mean_chisq),
                "lambda_gc": _scalar_or_value(hsq.lambda_gc),
                "ratio": _numeric_scalar(getattr(hsq, "ratio", None), "ratio"),
                "ratio_se": _numeric_scalar(getattr(hsq, "ratio_se", None), "ratio_se"),
            }
        ]
    )


def summarize_ld_score_regression_bins(
    hsq,
    *,
    ld_score,
    chi_square,
    sample_size,
    regression_ld_score,
    reference_snp_count: float,
    max_bins: int = 50,
) -> pd.DataFrame:
    """Summarize the exact fitted SNP population into LD-score rank bins.

    Parameters
    ----------
    hsq : ldsc._regression.estimators.Hsq-like object
        Final unpartitioned fit exposing ``coef``, ``tot``, and ``intercept``.
    ld_score, chi_square, sample_size, regression_ld_score : array_like
        Equal-length per-SNP values after the h2 chi-square filter. The order is
        preserved among tied LD scores.
    reference_snp_count : float
        Reference-SNP count used by the fitted h2 model.
    max_bins : int, optional
        Maximum number of nonempty equal-count rank bins. Default is 50.

    Returns
    -------
    pandas.DataFrame
        One row per ascending LD-score bin in
        :data:`ldsc._regression.outputs.H2_REGRESSION_BIN_COLUMNS` order. Chi-square SD uses
        ``ddof=1`` and is NaN for a singleton bin.

    Raises
    ------
    LDSCInternalError
        If the fitted arrays are empty or differ in length, or if ``max_bins``
        is less than one.

    Notes
    -----
    Fitted expectations use the final LDSC slope and intercept with each SNP's
    own sample size. Weights are the existing heritability regression weights
    evaluated at the final model; this function does not fit another model.
    """
    arrays = [
        np.ravel(np.asarray(values, dtype=np.float64))
        for values in (ld_score, chi_square, sample_size, regression_ld_score)
    ]
    lengths = {len(values) for values in arrays}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) == 0:
        raise LDSCInternalError(
            "h2 regression diagnostics require nonempty, equal-length LD-score, chi-square, "
            "sample-size, and regression-weight LD-score arrays."
        )
    if max_bins < 1:
        raise LDSCInternalError("h2 regression diagnostics require max_bins >= 1.")
    ld_values, chisq_values, n_values, w_ld_values = arrays
    n_snps = len(ld_values)
    n_bins = min(int(max_bins), n_snps)
    slope = _scalar(hsq.coef)
    total_h2 = _scalar(hsq.tot)
    intercept = _scalar(hsq.intercept)
    fitted = intercept + np.multiply(np.multiply(n_values, slope), ld_values)
    weights = np.ravel(
        reg.Hsq.weights(
            ld_values.reshape((-1, 1)),
            w_ld_values.reshape((-1, 1)),
            n_values.reshape((-1, 1)),
            float(reference_snp_count),
            total_h2,
            intercept,
        )
    )
    order = np.argsort(ld_values, kind="stable")
    rows: list[dict[str, float | int]] = []
    for bin_index, indices in enumerate(np.array_split(order, n_bins), start=1):
        bin_chisq = chisq_values[indices]
        rows.append(
            {
                "bin": bin_index,
                "n_snps": len(indices),
                "ld_score_min": float(np.min(ld_values[indices])),
                "ld_score_max": float(np.max(ld_values[indices])),
                "mean_ld_score": float(np.mean(ld_values[indices])),
                "mean_chi_square": float(np.mean(bin_chisq)),
                "sd_chi_square": float(np.std(bin_chisq, ddof=1)) if len(indices) > 1 else float("nan"),
                "mean_sample_size": float(np.mean(n_values[indices])),
                "mean_fitted_chi_square": float(np.mean(fitted[indices])),
                "mean_regression_weight": float(np.mean(weights[indices])),
            }
        )
    return pd.DataFrame(rows, columns=H2_REGRESSION_BIN_COLUMNS)


@serial_calculation
def summarize_partitioned_h2(
    hsq,
    dataset: RegressionDataset,
    annotation_columns: Sequence[str],
    samp_prev: float | None = None,
    pop_prev: float | None = None,
) -> pd.DataFrame:
    """Build overlap-aware partitioned-h2 rows from a fitted ``Hsq`` result.

    Assembles the model overlap matrix for the retained LD-score columns, runs
    the ported legacy ``_overlap_output`` (augmented with the one-sided
    coefficient p-value, conditional annotation heritability, and the overlap
    flag), and returns the requested ``annotation_columns`` rows in the single
    lowercase ``PARTITIONED_H2_COLUMNS`` schema. The whole-model total heritability
    is carried on every row (``total_h2_obs``, constant across one fitted model)
    so each row's ``category_h2_obs`` can be read against it. Per-annotation and
    total heritability are reported on both the observed (``*_obs``) and liability
    (``*_liab``) scales; liability columns are NaN unless both prevalences are
    finite probabilities in (0, 1). Proportions, enrichment, and coefficients are
    scale-invariant. ``annotation_columns`` must be a subset of
    ``dataset.retained_ld_columns``.
    """
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
            "Category": "category",
            "Prop._SNPs": "prop_snps",
            "Prop._h2": "prop_h2",
            "Prop._h2_std_error": "prop_h2_se",
            "Enrichment": "enrichment",
            "Enrichment_std_error": "enrichment_se",
            "Enrichment_p": "enrichment_p",
            "Coefficient": "coefficient",
            "Coefficient_std_error": "coefficient_se",
            "Coefficient_z": "coefficient_z",
            "Coefficient_p": "coefficient_p",
            "overlap_aware": "overlap_annot",
            "Category_h2": "category_h2_obs",
            "Category_h2_std_error": "category_h2_obs_se",
        }
    )
    c = (
        float(reg.liability_conversion_factor(samp_prev, pop_prev))
        if samp_prev is not None and pop_prev is not None
        else float("nan")
    )
    table["category_h2_liab"] = table["category_h2_obs"] * c
    table["category_h2_liab_se"] = table["category_h2_obs_se"] * c
    # Whole-model total heritability: a single scalar broadcast to every row of one
    # fitted model (constant across the functional regime; each query's own model
    # total in the cell-type regime).
    total_obs = _scalar(hsq.tot)
    total_obs_se = _scalar(hsq.tot_se)
    table["total_h2_obs"] = total_obs
    table["total_h2_obs_se"] = total_obs_se
    table["total_h2_liab"] = total_obs * c
    table["total_h2_liab_se"] = total_obs_se * c
    table["samp_prev"] = _prev_cell(samp_prev)
    table["pop_prev"] = _prev_cell(pop_prev)
    rows = table.set_index("category").loc[list(annotation_columns)].reset_index()
    return rows.loc[:, PARTITIONED_H2_COLUMNS].reset_index(drop=True)


def _summarize_rg_pair(
    rg_result,
    dataset: RegressionDataset,
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
    the observed (``*_obs``) and liability (``*_liab``) scales. The genetic correlation
    (``rg``/``rg_se``/``z``/``p``) is scale-invariant and unconverted. A trait's
    ``h2_*_liab`` is NaN unless its prevalence is finite; ``gencov_liab`` is NaN
    only when no prevalence was supplied for the run. ``n_snps_used`` is the
    post-product-filter SNP count from the regression fit; it falls
    back to the supplied dataset's row count. Workflows supply the fitted dataset.

    ``status`` is ``'ok'``, or ``'unestimable'`` with an ``error`` naming the
    estimate that is not positive when either heritability estimate is;
    that row keeps its count and its heritability, genetic covariance and
    intercept estimates, and its ``rg``/``rg_se``/``z``/``p`` are NaN.
    """
    n_snps = int(n_snps_used) if n_snps_used is not None else int(len(dataset.merged))
    estimate, h2_out_of_bounds = _kernel_rg_estimate(rg_result)
    h2_1_obs = _numeric_attr(getattr(rg_result, "hsq1", None), "tot", "h2_1_obs")
    h2_1_obs_se = _numeric_attr(getattr(rg_result, "hsq1", None), "tot_se", "h2_1_obs_se")
    h2_2_obs = _numeric_attr(getattr(rg_result, "hsq2", None), "tot", "h2_2_obs")
    h2_2_obs_se = _numeric_attr(getattr(rg_result, "hsq2", None), "tot_se", "h2_2_obs_se")
    gencov_obs = _numeric_attr(getattr(rg_result, "gencov", None), "tot", "gencov_obs")
    gencov_obs_se = _numeric_attr(getattr(rg_result, "gencov", None), "tot_se", "gencov_obs_se")
    h2_1_liab, h2_1_liab_se = _liability_pair(h2_1_obs, h2_1_obs_se, samp_prev_1, pop_prev_1)
    h2_2_liab, h2_2_liab_se = _liability_pair(h2_2_obs, h2_2_obs_se, samp_prev_2, pop_prev_2)
    if any(v is not None for v in (samp_prev_1, pop_prev_1, samp_prev_2, pop_prev_2)):
        gfac = reg.gencov_obs_to_liab(1.0, samp_prev_1, samp_prev_2, pop_prev_1, pop_prev_2)
        gencov_liab, gencov_liab_se = gencov_obs * gfac, gencov_obs_se * gfac
    else:
        gencov_liab, gencov_liab_se = float("nan"), float("nan")
    return {
        "trait_1": trait_1,
        "trait_2": trait_2,
        "n_snps_used": n_snps,
        **estimate,
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
        "intercept_gencov_se": _numeric_attr(
            getattr(rg_result, "gencov", None),
            "intercept_se",
            "intercept_gencov_se",
        ),
        "ratio_1": _numeric_attr(getattr(rg_result, "hsq1", None), "ratio", "ratio_1"),
        "ratio_1_se": _numeric_attr(getattr(rg_result, "hsq1", None), "ratio_se", "ratio_1_se"),
        "ratio_2": _numeric_attr(getattr(rg_result, "hsq2", None), "ratio", "ratio_2"),
        "ratio_2_se": _numeric_attr(getattr(rg_result, "hsq2", None), "ratio_se", "ratio_2_se"),
        "lambda_gc_1": _numeric_attr(getattr(rg_result, "hsq1", None), "lambda_gc", "lambda_gc_1"),
        "lambda_gc_2": _numeric_attr(getattr(rg_result, "hsq2", None), "lambda_gc", "lambda_gc_2"),
        "mean_chisq_1": _numeric_attr(getattr(rg_result, "hsq1", None), "mean_chisq", "mean_chisq_1"),
        "mean_chisq_2": _numeric_attr(getattr(rg_result, "hsq2", None), "mean_chisq", "mean_chisq_2"),
        "pair_kind": pair_kind,
        "status": "unestimable" if h2_out_of_bounds else "ok",
        "error": _h2_out_of_bounds_error(trait_1, trait_2, h2_1_obs, h2_2_obs, n_snps) if h2_out_of_bounds else "",
    }


def _kernel_rg_estimate(rg_result) -> tuple[dict[str, float], bool]:
    """Return a fitted pair's ``rg``, ``rg_se``, ``z`` and ``p``, and whether a heritability estimate is not positive.

    When either trait's heritability estimate is at most 0, the ``RG`` kernel
    ported from legacy LDSC sets ``_negative_hsq`` and stores the string
    ``'NA'`` in these four fields (legacy ``RG.summary()`` reports the pair as
    "h2 out of bounds"); they are returned as NaN with the flag set.
    """
    fields = {"rg": "rg_ratio", "rg_se": "rg_se", "z": "z", "p": "p"}
    if rg_result._negative_hsq:
        return dict.fromkeys(fields, math.nan), True
    return {column: _required_numeric_scalar(getattr(rg_result, attr), column) for column, attr in fields.items()}, False


def _h2_out_of_bounds_error(trait_1: str, trait_2: str, h2_1_obs: float, h2_2_obs: float, n_snps: int) -> str:
    """Explain an rg pair whose genetic correlation is undefined because a heritability estimate is not positive."""
    values = " and ".join(
        f"h2_{index}_obs (trait '{trait}') is {value:.3g}"
        for index, trait, value in ((1, trait_1, h2_1_obs), (2, trait_2, h2_2_obs)) if value <= 0
    )
    return (
        f"SNP heritability estimate not positive: on the pair's {n_snps:,} SNPs, {values}. "
        "The genetic correlation divides by the square root of the product of the two heritability estimates, "
        "so rg, rg_se, z and p are NaN; the heritability and genetic covariance estimates are kept. "
        "Most likely the trait has too little heritability signal on these SNPs for its sample size; "
        "a munging error, such as a wrong N or signed-statistic column, can also cause it. "
        f"Other causes & fixes: {_RG_H2_OUT_OF_BOUNDS_DOC}"
    )


def _failed_h2_row(
    *, trait_name: str, error: str, status: str, n_snps: int, samp_prev=None, pop_prev=None,
) -> pd.DataFrame:
    """Build the ``h2_per_trait`` row of a trait whose own ``rg`` heritability fit failed.

    The columns are those of :func:`summarize_total_h2` plus ``status`` and
    ``error``. ``n_snps`` is the trait's SNP count when its dataset was built,
    else 0; the estimates are NaN and the prevalences are the trait's inputs.
    """
    estimates = ("total_h2_obs", "total_h2_obs_se", "total_h2_liab", "total_h2_liab_se")
    diagnostics = ("intercept", "intercept_se", "mean_chisq", "lambda_gc", "ratio", "ratio_se")
    return pd.DataFrame([{
        "trait_name": trait_name, "n_snps": n_snps, **dict.fromkeys(estimates, math.nan),
        "samp_prev": _prev_cell(samp_prev), "pop_prev": _prev_cell(pop_prev),
        **dict.fromkeys(diagnostics, math.nan), "status": status, "error": error,
    }])


def _failed_rg_full_row(
    *, trait_1: str, trait_2: str, pair_kind: str, error: str, status: str = "failed", n_snps_used: int = 0,
    samp_prev_1=None, pop_prev_1=None, samp_prev_2=None, pop_prev_2=None,
) -> dict[str, object]:
    """Build a full rg row for a failed or unestimable pair.

    ``n_snps_used`` is the pair's SNP count when its dataset was built, else 0.
    The prevalences are the traits' inputs, which a failed pair keeps although
    it has no estimates.
    """
    row = {column: math.nan for column in RG_FULL_COLUMNS}
    row.update(
        {
            "trait_1": trait_1,
            "trait_2": trait_2,
            "n_snps_used": n_snps_used,
            "samp_prev_1": _prev_cell(samp_prev_1),
            "pop_prev_1": _prev_cell(pop_prev_1),
            "samp_prev_2": _prev_cell(samp_prev_2),
            "pop_prev_2": _prev_cell(pop_prev_2),
            "pair_kind": pair_kind,
            "status": status,
            "error": error,
        }
    )
    return row


def _concise_rg_row(full_row: dict[str, object], note: str) -> dict[str, object]:
    """Build the concise public rg row from a full row and its note: empty for an ``ok`` pair."""
    return {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "n_snps_used": full_row["n_snps_used"],
        "rg": full_row["rg"],
        "rg_se": full_row["rg_se"],
        "p": full_row["p"],
        "note": note,
    }


def _absolute_ldscore_dir(value) -> str | None:
    """Return the absolute LD-score directory to record, resolving a relative value against the working directory.

    Saved results link to their LD-score directory, which ``quantile-h2``
    reopens from any working directory, so a relative path would be ambiguous.
    """
    return None if value is None else str(Path(normalize_path_token(str(value))).resolve())


def _ldscore_dir_of(ldscore_result: LDScoreSource | LDScoreResult) -> str | None:
    """Return the absolute LD-score directory a result was read from, or ``None`` for one built in memory."""
    source_metadata = ldscore_result.output_paths.get("metadata")
    return str(Path(source_metadata).parent.resolve()) if source_metadata else None


def _genome_build_of(dataset: RegressionDataset) -> str | None:
    """Return the genome build of a dataset's LD scores, or ``None`` when they record none."""
    return None if dataset.config_snapshot is None else getattr(dataset.config_snapshot, "genome_build", None)


def _dataset_metadata(dataset: RegressionDataset) -> dict[str, object]:
    """Record a fitted dataset's SNP identity, count vector and LD-score columns in result metadata."""
    return {
        "effective_snp_identifier": dataset.effective_snp_identifier,
        "identity_downgrade_applied": dataset.identity_downgrade_applied,
        "count_key_used_for_regression": dataset.count_key_used_for_regression,
        "retained_ld_columns": list(dataset.retained_ld_columns),
        "dropped_zero_variance_ld_columns": list(dataset.dropped_zero_variance_ld_columns),
    }


def _fit_metadata(fit: RegressionFit, snps_key: str = "n_snps") -> dict[str, object]:
    """Record a fit's SNP and block counts, chi-square cap and two-step cutoff.

    ``rg`` pairs name the SNP count ``n_snps_used``.
    """
    return {snps_key: fit.n_snps, "n_blocks_used": fit.n_blocks, "effective_chisq_max": fit.effective_chisq_max,
            "effective_two_step_cutoff": fit.effective_two_step_cutoff}


def _h2_metadata(
    sumstats_table: SumstatsTable,
    ldscore_dir: str | None,
    outcome: RegressionFit,
    config: RegressionConfig,
    dropped_nonfinite_statistic_rows: int = 0,
) -> dict[str, object]:
    """Build h2 metadata from the realized fit and the config's checked prevalences.

    ``sumstats_file`` is the resolved path the table was read from, as in
    ``partitioned-h2`` and ``rg`` metadata.
    """
    return {
        "artifact_type": "h2_result",
        "trait_name": _trait_label(sumstats_table),
        "sumstats_file": sumstats_table.source_path,
        "ldscore_dir": ldscore_dir,
        "genome_build": _genome_build_of(outcome.dataset),
        **_dataset_metadata(outcome.dataset),
        "dropped_nonfinite_statistic_rows": dropped_nonfinite_statistic_rows,
        **_fit_metadata(outcome),
        "intercept_h2_policy": _intercept_policies(config)["intercept_h2_policy"],
        "samp_prev": config.samp_prev,
        "pop_prev": config.pop_prev,
        "scale": _scale_label(config.samp_prev, config.pop_prev),
    }


def _rg_prevalence_metadata(prevalences) -> dict[str, object]:
    """Per-pair prevalences and scale from the pair's two checked ``(samp_prev, pop_prev)`` inputs."""
    (samp_prev_1, pop_prev_1), (samp_prev_2, pop_prev_2) = prevalences
    return {"samp_prev_1": samp_prev_1, "pop_prev_1": pop_prev_1, "samp_prev_2": samp_prev_2,
            "pop_prev_2": pop_prev_2, "scale": _scale_label(samp_prev_1, pop_prev_1, samp_prev_2, pop_prev_2)}


def _rg_pair_metadata(
    table_1: SumstatsTable,
    table_2: SumstatsTable,
    full_row: dict[str, object],
    config: RegressionConfig,
    pair_kind: str,
    prevalences,
    fit: dict[str, object] | None = None,
) -> dict[str, object]:
    """Build the metadata of one rg pair; ``prevalences`` holds the two traits' checked pairs.

    ``fit`` holds the pair's SNP and block counts, cap and dataset settings:
    from its ``RegressionFit``, or for a failed pair from the ``fit_summary``
    its error carries. Without it the pair failed before its dataset was
    built, and the counts are 0 with the requested cap.
    """
    metadata = {
        "trait_1": full_row["trait_1"],
        "trait_2": full_row["trait_2"],
        "source_1": table_1.source_path,
        "source_2": table_2.source_path,
        "pair_kind": pair_kind,
        "status": full_row["status"],
        "error": full_row["error"],
        "n_snps_used": 0,
        "n_blocks_used": 0,
        "effective_chisq_max": config.chisq_max,
        "effective_two_step_cutoff": None,
        **_rg_prevalence_metadata(prevalences),
        **_intercept_policies(config),
    }
    if fit is not None:
        metadata.update(fit)
    return metadata


def _intercept_policies(config: RegressionConfig) -> dict[str, str]:
    """Describe the heritability and genetic-covariance intercepts: ``free`` or ``fixed:<value>``."""
    policy = _fit_policy(config)
    return {
        f"intercept_{name}_policy": "free" if value is None else f"fixed:{float(value):g}"
        for name, value in (("h2", policy.intercept_h2), ("gencov", policy.intercept_gencov))
    }


def _numeric_attr(obj: object | None, attr: str, field: str) -> float:
    """Return a numeric scalar attribute or NaN when unavailable."""
    if obj is None:
        return math.nan
    return _numeric_scalar(getattr(obj, attr, None), field)


def _required_numeric_scalar(value: object, field: str) -> float:
    """Return a numeric scalar or raise for kernel string/invalid headline values."""
    return _numeric_scalar(value, field, fail_on_non_numeric=True)


def _numeric_scalar(value: object, field: str, *, fail_on_non_numeric: bool = False) -> float:
    """Return the first numeric scalar from common kernel return shapes."""
    if value is None:
        return math.nan
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.upper() == "NA":
            if fail_on_non_numeric:
                raise LDSCInternalError(
                    f"Regression summary field {field} is the non-numeric kernel value 'NA'. "
                    "Most likely the estimator returned an unavailable headline statistic where a number was required. "
                    "Re-run with `--log-level DEBUG` and report the traceback."
                )
            return math.nan
        try:
            return float(stripped)
        except ValueError as exc:
            if fail_on_non_numeric:
                raise LDSCInternalError(
                    f"Regression summary field {field} is non-numeric kernel value {value!r}. "
                    "Most likely the estimator returned an unavailable headline statistic where a number was required. "
                    "Re-run with `--log-level DEBUG` and report the traceback."
                ) from exc
            return math.nan
    if hasattr(value, "__array__") or isinstance(value, (list, tuple)):
        array = np.ravel(value)
        if array.size == 0:
            return math.nan
        return _numeric_scalar(array[0].item() if hasattr(array[0], "item") else array[0], field, fail_on_non_numeric=fail_on_non_numeric)
    try:
        if pd.isna(value):
            return math.nan
    except TypeError:
        pass
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        if fail_on_non_numeric:
            raise LDSCInternalError(
                f"Regression summary field {field} is non-numeric kernel value {value!r}. "
                "Most likely the estimator returned an unavailable headline statistic where a number was required. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            ) from exc
        return math.nan


def _format_exception(exc: Exception) -> str:
    """Return compact user-facing error text for a pair failure."""
    message = str(exc).strip()
    return f"{exc.__class__.__name__}: {message}" if message else exc.__class__.__name__


def _scalar(value) -> float:
    """Return the first scalar numeric value from a possibly array-like object."""
    return float(np.ravel(value)[0])


def _coefficient_delete_frame(hsq, annotation_columns: Sequence[str]) -> pd.DataFrame:
    """Return float64 delete-one-block coefficient values in fitted order.

    ``Hsq.part_delete_values`` already contains annotation coefficients on the
    per-SNP heritability scale. The intercept delete column is intentionally
    excluded because post-fit annotation projection consumes only ``tau``.
    """
    raw_values = getattr(hsq, "part_delete_values", None)
    if not isinstance(raw_values, (np.ndarray, list, tuple, pd.DataFrame)):
        raw_values = np.empty((0, len(annotation_columns)))
    values = np.asarray(raw_values, dtype=np.float64)
    if values.ndim != 2 or values.shape[1] != len(annotation_columns):
        raise LDSCInternalError(
            "partitioned-h2 coefficient delete values do not match the retained annotation order. "
            f"Expected {len(annotation_columns)} columns, got shape {values.shape}. "
            "Re-run with `--log-level DEBUG` and report the traceback."
        )
    frame = pd.DataFrame(values, columns=list(annotation_columns), dtype=np.float64)
    frame.insert(0, "delete_block", np.arange(len(frame), dtype=np.int64))
    return frame


def _scalar_or_value(value):
    """Return a scalar float for array-like values, otherwise preserve the original value."""
    if hasattr(value, "__array__") or isinstance(value, (list, tuple)):
        return _scalar(value)
    return value


def _add_jackknife_model_context(error, merged, columns):
    """Attach filtered-row support and genomic spans to failed deletions.

    A failure with one more parameter than ``columns`` had a free intercept,
    whose column follows the LD-score columns. An rg pair's three fits can
    differ here: the genetic covariance intercept may be fixed while the
    heritability intercepts are free.
    """
    support = np.count_nonzero(merged[columns].to_numpy(), axis=0)
    for failure in error.failures:
        names = [*columns, *(["intercept"] if failure["n_parameters"] == len(columns) + 1 else [])]
        start, end = failure["row_start"], failure["row_end"]
        deleted = merged.iloc[start:end]
        remaining = support - np.count_nonzero(deleted[columns].to_numpy(), axis=0)
        failure["columns"] = names
        failure["zero_columns"] = [names[i] for i in failure["zero_column_indices"]]
        failure["support_counts"] = {
            name: {"full": int(total), "remaining": int(left)}
            for name, total, left in zip(columns, support, remaining)
        }
        chr_column = "_ldsc_panel_CHR" if "_ldsc_panel_CHR" in deleted else "CHR"
        pos_column = "_ldsc_panel_POS" if "_ldsc_panel_POS" in deleted else "POS"
        if {chr_column, pos_column}.issubset(deleted.columns):
            failure["genomic_spans"] = [
                {"CHR": str(chrom), "POS_start": int(rows[pos_column].min()), "POS_end": int(rows[pos_column].max())}
                for chrom, rows in deleted.dropna(subset=[pos_column]).groupby(chr_column, sort=False)
            ]


def _raise_on_model_collinearity(dataset: RegressionDataset, x: np.ndarray) -> None:
    """Abort the fit when the LD-score design matrix is near-collinear.

    Matches legacy LDSC's hard stop on an ill-conditioned design matrix: rather
    than warn and emit a partition the user must second-guess, collinear
    annotations raise ``LDSCInputError`` so the only results produced are
    trustworthy ones. h2 fits and each rg trait pair apply it on their own
    regression rows; it is skipped when the LD-score directory stores no
    overlap or the model has fewer than two LD-score columns.
    """
    overlap = getattr(dataset, "ldscore_overlap", None)
    if overlap is None or len(dataset.retained_ld_columns) < 2:
        return

    use_common = dataset.count_key_used_for_regression == COMMON_COUNT_KEY
    try:
        model_overlap = assemble_model_overlap(overlap, list(dataset.retained_ld_columns), use_common)
    except (ValueError, KeyError):
        return
    message = model_collinearity_error(x, list(dataset.retained_ld_columns), model_overlap)
    if message:
        raise LDSCInputError(message)


@dataclass(frozen=True)
class _FitPolicy:
    """Intercepts and two-step setting of a regression run, resolved once from its config.

    ``None`` marks an intercept that the fit estimates.
    """

    intercept_h2: float | None
    intercept_gencov: float | None
    requested_two_step_cutoff: float | None

    def two_step_cutoff(self, n_ld_columns: int) -> float | None:
        """Two-step cutoff of a model with ``n_ld_columns`` LD-score columns, or ``None`` for no two-step fit.

        A requested cutoff is used as given. Otherwise a single-annotation
        model with an estimated heritability intercept uses 30, as legacy
        LDSC's ``estimate_h2`` and ``estimate_rg`` do (``sumstats.py``).
        """
        if self.requested_two_step_cutoff is not None:
            return self.requested_two_step_cutoff
        if self.intercept_h2 is None and n_ld_columns == 1:
            return 30
        return None


def _fit_policy(config: RegressionConfig) -> _FitPolicy:
    """Resolve a config's intercepts and two-step request for ``h2``, ``partitioned-h2`` and ``rg`` fits.

    ``use_intercept=False`` fixes the heritability intercept at 1 and the
    genetic-covariance intercept at 0; otherwise ``intercept_h2`` and
    ``intercept_gencov`` fix them when set.
    """
    if not config.use_intercept:
        return _FitPolicy(1, 0, config.two_step_cutoff)
    return _FitPolicy(
        None if config.intercept_h2 is None else float(config.intercept_h2),
        None if config.intercept_gencov is None else float(config.intercept_gencov),
        config.two_step_cutoff,
    )


def _resolve_default_chisq_max(chisq_max: float | None, n_annot: int, n_max: float) -> float | None:
    """Effective chi-square ceiling for an h2 fit, matching legacy ``estimate_h2``.

    Legacy LDSC (``ldscore/sumstats.py`` ``estimate_h2``) caps multi-annotation
    (partitioned) models at ``max(0.001 * N.max(), 80)`` when ``--chisq-max`` is
    unset, to down-weight outlier SNPs whose extreme chi-square would otherwise
    dominate the regression. Single-annotation models stay uncapped and rely on
    the two-step estimator instead. An explicit ``chisq_max`` always takes
    precedence over the default.
    """
    if chisq_max is not None:
        return float(chisq_max)
    if n_annot > 1:
        return max(0.001 * float(n_max), 80.0)
    return None
