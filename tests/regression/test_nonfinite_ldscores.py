"""A missing or infinite LD score or regression weight stops regression before any fit."""

import gzip

import numpy as np
import pandas as pd
import pytest

from ldsc import load_ldscore_from_dir
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc._regression import runner as workflow
from ldsc.config import RegressionConfig
from ldsc.errors import LDSCInputError
from tests.regression.test_regression_streaming import batch_inputs
from tests.regression.test_trait_preparation import N_SNPS, munged, panel, signed_statistics
from tests.worker_fixtures import calculation_body

ANCHOR = "troubleshooting.md#regression-ld-scores-or-regression-weights-are-missing-or-infinite"


def damaged_panel(**values):
    """Return the shared panel with ``column=(row, value)`` replacements."""
    source = panel()
    for column, (row, value) in values.items():
        source.baseline_table.loc[row, column] = value
    return source


def test_h2_command_stops_on_an_infinite_ld_score_even_for_a_snp_the_trait_lacks(tmp_path):
    LDScoreDirectoryWriter().write(damaged_panel(base=(N_SNPS - 1, np.inf)),
                                   LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    z = signed_statistics(np.random.default_rng(1), panel().baseline_table["base"])
    rows = [f"rs{i}\tA\tC\t{z[i]:.5f}\t10000" for i in range(N_SNPS - 1)]  # no row for the damaged SNP
    sumstats = tmp_path / "trait.sumstats.gz"
    with gzip.open(sumstats, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tA1\tA2\tZ\tN\n" + "\n".join(rows) + "\n")

    with pytest.raises(LDSCInputError) as caught:
        workflow.run_h2(sumstats_file=str(sumstats), ldscore_dir=str(tmp_path / "ld"),
                        output_dir=str(tmp_path / "out"), num_blocks=10)

    message = str(caught.value)
    assert "'base' (1 row; for example rs159)" in message
    assert ANCHOR in message
    assert not (tmp_path / "out/h2.tsv").exists()


def test_one_error_names_every_damaged_column():
    source = damaged_panel(base=(3, -np.inf), regression_ld_scores=(7, np.nan))
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    with pytest.raises(LDSCInputError) as caught:
        runner.build_dataset(munged("A", "chr_pos_allele_aware", 2), source)

    message = str(caught.value)
    assert "'base' (1 row; for example rs3)" in message
    assert "'regression_ld_scores' (1 row; for example rs7)" in message


def test_rg_stops_before_any_pair_instead_of_dropping_rows():
    source = damaged_panel(base=(5, np.nan))
    tables = [munged("A", "chr_pos_allele_aware", 3), munged("B", "chr_pos_allele_aware", 4)]
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    with calculation_body(), pytest.raises(LDSCInputError, match="'base' \\(1 row; for example rs5\\)"):
        runner.estimate_rg_pairs(tables, source)


def damage_parquet(path, column, row, value):
    frame = pd.read_parquet(path)
    frame[column] = frame[column].astype("float64")
    frame.loc[row, column] = value
    frame.to_parquet(path, index=False, row_group_size=30)


def test_partitioned_h2_fails_only_the_query_with_an_infinite_ld_score(tmp_path):
    table, _ = batch_inputs(tmp_path / "ld")
    damage_parquet(tmp_path / "ld/ldscore.query.parquet", "second", 12, np.inf)
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=6))

    result = runner.estimate_partitioned_h2_batch(table, load_ldscore_from_dir(str(tmp_path / "ld")),
                                                  output_dir=tmp_path / "out", continue_on_query_error=True)

    status = result.query_status.set_index("query_annotation")
    assert status.loc[["first", "third"], "status"].tolist() == ["success", "success"]
    assert status.loc["second", "status"] == "failed"
    assert status.loc["second", "stage"] == "model_preparation"
    assert "query annotation 'second'" in status.loc["second", "error_message"]
    assert "1 row" in status.loc["second", "error_message"]


def test_partitioned_h2_stops_on_a_missing_baseline_ld_score(tmp_path):
    table, _ = batch_inputs(tmp_path / "ld")
    damage_parquet(tmp_path / "ld/ldscore.baseline.parquet", "base", 4, np.nan)
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=6))

    with pytest.raises(LDSCInputError, match="'base' \\(1 row"):
        runner.estimate_partitioned_h2_batch(table, load_ldscore_from_dir(str(tmp_path / "ld")),
                                             output_dir=tmp_path / "out", continue_on_query_error=True)
    assert not (tmp_path / "out/partitioned_h2.tsv").exists()
