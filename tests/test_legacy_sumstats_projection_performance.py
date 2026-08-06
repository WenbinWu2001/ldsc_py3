"""Focused performance regression coverage for legacy sumstats projection."""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
import pytest

from ldsc.config import GlobalConfig
from ldsc.ldscore_calculator import LDScoreResult
from ldsc.regression_runner import _project_legacy_sumstats_to_panel
from ldsc.sumstats_munger import SumstatsTable


@pytest.mark.slow
def test_allele_unaware_legacy_projection_scales_without_per_row_pandas_work() -> None:
    """A 20k-row projection should stay comfortably below the scalar baseline."""
    n_source = 20_000
    n_missing_allele = 2_500
    n_missing_panel = 1_000
    n_retained = n_source - n_missing_allele - n_missing_panel
    source_rsids = np.asarray([f"rs{index}" for index in range(n_source)], dtype=object)
    source = pd.DataFrame(
        {
            "SNP": source_rsids,
            "A1": np.asarray(["A"] * n_source, dtype=object),
            "A2": np.asarray(["C"] * n_source, dtype=object),
            "Z": np.ones(n_source),
            "N": np.full(n_source, 1000.0),
            "FRQ": np.full(n_source, 0.25),
        }
    )
    source.loc[n_retained : n_retained + n_missing_allele - 1, "A1"] = None
    panel_rsids = source_rsids[: n_source - n_missing_panel]
    panel_frame = pd.DataFrame(
        {
            "CHR": np.asarray(["1"] * len(panel_rsids), dtype=object),
            "POS": np.arange(1, len(panel_rsids) + 1),
            "SNP": panel_rsids,
            "regression_ld_scores": np.ones(len(panel_rsids)),
            "base": np.ones(len(panel_rsids)),
        }
    )
    panel = LDScoreResult(
        baseline_table=panel_frame,
        query_table=None,
        count_records=[],
        baseline_columns=["base"],
        query_columns=[],
        ld_reference_snps=frozenset(),
        ld_regression_snps=frozenset(panel_rsids),
        chromosome_results=[],
        config_snapshot=GlobalConfig(snp_identifier="rsid", genome_build=None),
    )
    legacy = SumstatsTable(
        data=source,
        has_alleles=True,
        source_path="benchmark.sumstats.gz",
        trait_name="benchmark",
        provenance={"source_format": "ldsc2_sumstats", "legacy_ldsc2": True},
        config_snapshot=None,
    )

    started = time.perf_counter()
    projected, drops = _project_legacy_sumstats_to_panel(legacy, panel)
    elapsed = time.perf_counter() - started

    assert len(projected.data) == n_retained
    assert len(drops) == n_missing_allele + n_missing_panel
    assert elapsed < 3.0, f"20k-row projection took {elapsed:.3f}s; expected vectorized execution"
