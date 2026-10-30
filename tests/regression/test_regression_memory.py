"""Borrowed model inputs and independent weighted least-squares expectations."""

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from ldsc._regression import runner as workflow
from ldsc._regression import estimators as regression
from ldsc._regression.partitioned_parallel import _stage_frame
from ldsc._runtime.serial_calculation import serial_calculation
from tests.regression.test_regression_streaming import batch_inputs


@serial_calculation
def _check_shared_models(tmp_path, chisq_max):
    table, source = batch_inputs(tmp_path / "ld")
    runner = workflow.RegressionRunner(source.config_snapshot, workflow.RegressionConfig(num_blocks=6, chisq_max=chisq_max))
    prepared = runner.prepare_inputs(table, source)
    snapshot = _stage_frame(prepared.dataset.merged, tmp_path / "shared")
    with snapshot.open() as shared:
        prepared = replace(prepared, dataset=replace(prepared.dataset, merged=shared))
        before = shared.copy(deep=True)
        first_fit = None
        for query in ["first", "second", "third", "first"]:
            values = runner._read_query_batch(prepared, [query])
            if query == "second":
                values[query] = 1.0
            query_before = values.copy(deep=True)
            dataset = runner._dataset_from_query_values(prepared, [query], values)
            retained = ["base"] if query == "second" else ["base", query]
            assert dataset.retained_ld_columns == retained
            for column in snapshot.numeric:
                actual = dataset.merged[column].to_numpy(copy=False)
                assert np.shares_memory(actual, shared[column].to_numpy(copy=False))
                assert not actual.flags.writeable
            outcome = runner.estimate_h2(dataset)
            keep = before.Z.to_numpy() ** 2 <= chisq_max
            assert outcome.dataset.merged.SNP.tolist() == before.loc[keep, "SNP"].tolist()
            assert outcome.n_blocks == 6
            if keep.all():
                assert outcome.dataset.merged is dataset.merged
            pd.testing.assert_frame_equal(shared, before, check_exact=True)
            pd.testing.assert_frame_equal(values, query_before, check_exact=True)
            if query == "first":
                if first_fit is None:
                    first_fit = outcome.estimator
                else:
                    for field in ("coef", "coef_cov", "part_delete_values", "tot_delete_values"):
                        np.testing.assert_array_equal(getattr(outcome.estimator, field), getattr(first_fit, field))


@pytest.mark.parametrize("chisq_max", [80, 1e6])
def test_repeated_models_borrow_readonly_inputs_without_mutating_them(tmp_path, chisq_max):
    _check_shared_models(tmp_path, chisq_max)


@pytest.mark.parametrize("input_dtype", [np.float64, np.float32, np.int64])
@pytest.mark.parametrize("intercept", [1, None])
def test_partitioned_kernel_matches_independent_weighted_lstsq(input_dtype, intercept):
    rng = np.random.default_rng(7)
    x = rng.uniform(1, 8, (300, 3)).astype(input_dtype)
    n = np.full((300, 1), 10000, dtype=input_dtype)
    w = rng.uniform(1, 3, (300, 1))
    y = rng.uniform(1, 3, (300, 1))
    m = np.full((1, 3), 1000.)
    before = [array.copy() for array in (x, n, w, y, m)]
    for array in (x, n, w, y, m):
        array.flags.writeable = False
    fitted = regression.Hsq(y, x, w, n, m, n_blocks=10, old_weights=True, intercept=intercept)

    # The established initial LDSC weights, followed by independent QR/SVD
    # least-squares fits on the full population and explicit leave-block-out rows.
    total_ld = np.sum(x, axis=1).reshape((-1, 1))
    aggregate = regression.Hsq.aggregate(y, total_ld, n, float(m.sum()), intercept)
    weights = np.sqrt(regression.Hsq.weights(total_ld, w, n, float(m.sum()), aggregate, intercept))
    weights = weights / float(weights.sum())
    design = np.multiply(n, x) / np.mean(n)
    response = y - intercept if intercept is not None else y
    if intercept is None:
        design = np.column_stack([design, np.ones(len(design))])
    design, response = design * weights, response * weights
    expected = np.linalg.lstsq(design, response, rcond=None)[0].T
    deletes = []
    for block in range(10):
        keep = np.ones(300, dtype=bool)
        keep[block * 30:(block + 1) * 30] = False
        deletes.append(np.linalg.lstsq(design[keep], response[keep], rcond=None)[0][:, 0])
    # Well-conditioned small systems: float64 solvers agree well below 1e-10.
    np.testing.assert_allclose(fitted.jknife.est, expected, rtol=1e-10, atol=1e-14)
    np.testing.assert_allclose(fitted.jknife.delete_values, deletes, rtol=1e-10, atol=1e-14)
    assert fitted.jknife.est.dtype == np.float64
    for actual, original in zip((x, n, w, y, m), before):
        np.testing.assert_array_equal(actual, original)
