"""Fit population and diagnostics stay attached to the actual estimator inputs."""

from types import SimpleNamespace
import json

import numpy as np
import pandas as pd
import pytest

from ldsc import GlobalConfig, RegressionConfig
from ldsc import regression_runner as workflow
from ldsc.outputs import H2DirectoryWriter, H2OutputConfig
from ldsc.sumstats_munger import SumstatsTable


def h2_dataset():
    ld = np.arange(1., 61.)
    n = 1000 + 10 * ld
    chisq = 1 + .0002 * n * ld + .02 * np.sin(ld)
    chisq[0], chisq[-1] = 200., 100.
    return workflow.RegressionDataset(
        merged=pd.DataFrame({'SNP': [f'rs{i}' for i in range(60)], 'Z': np.sqrt(chisq), 'N': n,
                             'base': ld, 'weight': np.ones(60)}),
        ref_ld_columns=['base'], retained_ld_columns=['base'], dropped_zero_variance_ld_columns=[],
        weight_column='weight', reference_snp_count_totals={'common_reference_snp_counts': np.array([1000.])},
        count_key_used_for_regression='common_reference_snp_counts', trait_names=['trait'], chromosomes_aggregated=['1'],
        config_snapshot=GlobalConfig(snp_identifier='rsid'), effective_snp_identifier='rsid',
    )


@pytest.mark.parametrize('intercept', [None, 1.0])
def test_h2_fit_outcome_and_written_bins_use_exact_final_population(tmp_path, intercept):
    dataset = h2_dataset()
    config = RegressionConfig(chisq_max=100, n_blocks=6, intercept_h2=intercept)
    runner = workflow.RegressionRunner(GlobalConfig(snp_identifier='rsid'), config)
    outcome = runner._fit_h2_dataset(dataset)
    assert isinstance(outcome.estimator, workflow.reg.Hsq)
    assert outcome.dataset.merged.SNP.tolist() == [f'rs{i}' for i in range(1, 60)]
    assert len(dataset.merged) == 60
    assert outcome.n_snps == 59
    assert outcome.n_blocks == 6
    assert outcome.effective_chisq_max == 100
    assert outcome.dataset.retained_ld_columns == ['base']
    bins = outcome.diagnostic_bins
    assert bins.n_snps.sum() == 59
    assert bins.ld_score_min.min() == 2
    assert bins.ld_score_max.max() == 60
    np.testing.assert_allclose(np.average(bins.mean_chi_square, weights=bins.n_snps),
                               np.mean(dataset.merged.Z.to_numpy()[1:] ** 2))
    np.testing.assert_allclose(np.average(bins.mean_sample_size, weights=bins.n_snps), 1310.)
    table = SumstatsTable(data=dataset.merged, trait_name='trait', has_alleles=False, source_path='trait')
    metadata = workflow._h2_metadata(SimpleNamespace(), table, outcome)
    summary = workflow.summarize_total_h2(outcome.estimator, outcome.dataset)
    H2DirectoryWriter().write(summary, H2OutputConfig(output_dir=tmp_path), metadata=metadata, diagnostic_bins=bins)
    saved = json.loads((tmp_path / 'diagnostics/metadata.json').read_text())
    assert saved['n_snps'] == 59
    assert saved['effective_chisq_max'] == 100
    assert pd.read_csv(tmp_path / 'h2.tsv', sep='\t').n_snps.tolist() == [59]
    assert pd.read_csv(tmp_path / 'diagnostics/ld_score_regression_bins.tsv', sep='\t').n_snps.sum() == 59
    public = runner.estimate_h2(dataset)
    assert isinstance(public, workflow.reg.Hsq)
    np.testing.assert_allclose(public.coef, outcome.estimator.coef)


def test_rg_outcome_uses_explicit_product_filter_not_h2_cap():
    h2 = h2_dataset()
    merged = h2.merged.rename(columns={'Z': 'Z1', 'N': 'N1'}).copy()
    merged['Z1'] = np.sqrt(1.2 + .01 * merged.base)
    merged['Z2'] = np.sqrt(1.1 + .015 * merged.base)
    merged['N2'] = 1500 + 5 * merged.base
    merged.loc[0, ['Z1', 'Z2']] = [20., .1]  # one large chi-square, product only 4
    merged.loc[1, ['Z1', 'Z2']] = [4., 4.]  # product 256: excluded
    merged.loc[2, ['Z1', 'Z2']] = [2., 5.]  # product 100: inclusive boundary
    fields = {name: getattr(h2, name) for name in (
        'ref_ld_columns', 'retained_ld_columns', 'dropped_zero_variance_ld_columns', 'weight_column',
        'reference_snp_count_totals', 'count_key_used_for_regression', 'chromosomes_aggregated', 'config_snapshot',
    )}
    dataset = workflow.RGRegressionDataset(merged=merged, trait_names=['one', 'two'], **fields)
    runner = workflow.RegressionRunner(GlobalConfig(snp_identifier='rsid'), RegressionConfig(chisq_max=10, n_blocks=5, use_intercept=False))
    outcome = runner._fit_rg_dataset(dataset)
    assert isinstance(outcome.estimator, workflow.reg.RG)
    assert outcome.dataset.merged.SNP.tolist() == ['rs0'] + [f'rs{i}' for i in range(2, 60)]
    assert outcome.n_snps == 59
    assert outcome.n_blocks == 5
    assert outcome.effective_chisq_max == 10
    assert len(dataset.merged) == 60
    uncapped = runner._fit_rg_dataset(dataset, RegressionConfig(n_blocks=5, use_intercept=False))
    assert uncapped.n_snps == 60
    assert uncapped.effective_chisq_max is None


def test_partitioned_default_cap_uses_prefilter_sample_size_and_retains_annotation_order():
    from dataclasses import replace

    dataset = h2_dataset()
    merged = dataset.merged.copy()
    merged['query'] = 2 + np.sin(merged.base)
    merged.loc[0, ['N', 'Z']] = [400000., 21.]
    merged.loc[59, 'Z'] = 20.
    dataset = replace(
        dataset, merged=merged, ref_ld_columns=['base', 'query'], retained_ld_columns=['base', 'query'],
        reference_snp_count_totals={'common_reference_snp_counts': np.array([1000., 100.])},
    )
    runner = workflow.RegressionRunner(GlobalConfig(snp_identifier='rsid'), RegressionConfig(n_blocks=6))
    outcome = runner._fit_h2_dataset(dataset)
    assert outcome.effective_chisq_max == 400
    assert outcome.dataset.merged.SNP.tolist() == [f'rs{i}' for i in range(1, 60)]
    assert outcome.dataset.retained_ld_columns == ['base', 'query']
    assert outcome.diagnostic_bins is None
    assert outcome.n_blocks == 6
    # Independent kernel invocation on the known retained rows verifies numerical parity.
    rows = merged.iloc[1:]
    expected = workflow.reg.Hsq(
        rows[['Z']].to_numpy() ** 2, rows[['base', 'query']].to_numpy(),
        rows[['weight']].to_numpy(), rows[['N']].to_numpy(), np.array([[1000., 100.]]),
        n_blocks=6, intercept=None, twostep=None, old_weights=True,
    )
    np.testing.assert_allclose(outcome.estimator.coef, expected.coef)
    np.testing.assert_allclose(outcome.estimator.coef_cov, expected.coef_cov)
    np.testing.assert_allclose(outcome.estimator.part_delete_values, expected.part_delete_values)
