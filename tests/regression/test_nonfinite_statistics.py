"""Rows whose Z or N is missing or infinite are dropped once per trait, with diagnostics."""

import gzip
import json
import logging

import numpy as np
import pandas as pd

from ldsc._annotation.overlap import OverlapContribution
from ldsc._ldscore.calculator import LDScoreResult
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc._ldscore.overlap_table import LDScoreOverlap
from ldsc._regression import runner as workflow
from ldsc._regression.outputs import RgDirectoryWriter, RgOutputConfig
from ldsc.config import RegressionConfig
from tests.regression.test_trait_preparation import N_SNPS, munged, panel, signed_statistics
from tests.worker_fixtures import calculation_body


def with_nonfinite(table, z_rows=(), n_rows=()):
    data = table.data.copy()
    data.loc[list(z_rows), "Z"] = np.nan
    data.loc[list(n_rows), "N"] = np.inf
    return workflow.replace(table, data=data)


def test_h2_command_drops_nonfinite_legacy_rows_and_records_them(tmp_path):
    LDScoreDirectoryWriter().write(panel(), LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    rng = np.random.default_rng(11)
    z = signed_statistics(rng, panel().baseline_table["base"])
    rows = [f"rs{i}\tA\tC\t{z[i]:.5f}\t10000" for i in range(N_SNPS)]
    rows[3] = "rs3\tA\tC\tnan\t10000"
    rows[7] = f"rs7\tA\tC\t{z[7]:.5f}\tinf"
    sumstats = tmp_path / "trait.sumstats.gz"
    with gzip.open(sumstats, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tA1\tA2\tZ\tN\n" + "\n".join(rows) + "\n")

    workflow.run_h2(sumstats_file=str(sumstats), ldscore_dir=str(tmp_path / "ld"),
                    output_dir=str(tmp_path / "out"), num_blocks=10)

    audit = pd.read_csv(tmp_path / "out/diagnostics/dropped_snps/legacy_sumstats.tsv.gz", sep="\t")
    assert audit[["SNP", "reason"]].values.tolist() == [["rs3", "nonfinite_statistic"],
                                                         ["rs7", "nonfinite_statistic"]]
    metadata = json.loads((tmp_path / "out/diagnostics/metadata.json").read_text())
    assert metadata["dropped_nonfinite_statistic_rows"] == 2
    assert metadata["n_snps"] == N_SNPS - 2
    assert "whose Z or N is missing or infinite" in (tmp_path / "out/diagnostics/h2.log").read_text()


def test_h2_drops_nonfinite_rows_of_a_python_table():
    table = with_nonfinite(munged("A", "chr_pos_allele_aware", 12), z_rows=[0], n_rows=[5, 9])
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    dataset = runner.build_dataset(table, panel())
    hsq = runner.estimate_h2(dataset)

    assert len(dataset.merged) == N_SNPS - 3
    assert not {"rs0", "rs5", "rs9"} & set(dataset.merged["SNP"])
    assert np.isfinite(hsq.tot).all()


def functional_panel():
    source = panel()
    baseline = source.baseline_table.assign(b2=lambda frame: frame["base"] * (np.arange(N_SNPS) % 3 == 0))
    block = np.array([[1000.0, 300.0], [300.0, 300.0]])
    return LDScoreResult(
        baseline, None,
        [{"group": "baseline", "column": column, "all_reference_snp_count": count,
          "common_reference_snp_count": count} for column, count in (("base", 1000.0), ("b2", 300.0))],
        ["base", "b2"], [], frozenset(), frozenset(baseline.SNP), [], config_snapshot=source.config_snapshot,
        overlap=LDScoreOverlap.from_contribution(
            OverlapContribution(block, block, np.array([]), np.array([]), 1000, 1000),
            baseline_columns=["base", "b2"], query_columns=[]),
    )


def test_capped_partitioned_h2_counts_a_missing_z_as_nonfinite_not_as_an_outlier(tmp_path, caplog):
    complete = munged("A", "chr_pos_allele_aware", 13)
    without_row = workflow.replace(complete, data=complete.data.drop(index=4).reset_index(drop=True))
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    runs = {}
    for name, table in (("missing", with_nonfinite(complete, z_rows=[4])), ("absent", without_row)):
        caplog.clear()
        with calculation_body(), caplog.at_level(logging.INFO, logger="LDSC"):
            runner.estimate_partitioned_h2_batch(table, functional_panel(), output_dir=tmp_path / name)
        metadata = json.loads((tmp_path / name / "diagnostics/metadata.json").read_text())
        removed = [record.getMessage() for record in caplog.records if "SNPs with chi^2" in record.getMessage()]
        runs[name] = (metadata["dropped_nonfinite_statistic_rows"], metadata["n_snps"], removed)

    # The capped model sees the same rows as if the row were absent; the cap's
    # outlier count no longer includes the missing Z.
    assert runs["missing"] == (1, runs["absent"][1], runs["absent"][2])
    assert runs["absent"][0] == 0


def test_rg_finishes_with_a_nonfinite_trait_and_records_counts_per_trait(tmp_path):
    tables = [with_nonfinite(munged("A", "chr_pos_allele_aware", 14), z_rows=[2]),
              munged("B", "chr_pos_allele_aware", 15)]
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    result = runner.estimate_rg_pairs(tables, panel())
    RgDirectoryWriter().write(result, RgOutputConfig(output_dir=tmp_path / "out"))

    assert result.h2_per_trait["n_snps"].tolist() == [N_SNPS - 1, N_SNPS]
    assert result.rg_full["status"].tolist() == ["ok"]
    assert result.dropped_nonfinite_statistic_rows == [1, 0]
    metadata = json.loads((tmp_path / "out/diagnostics/metadata.json").read_text())
    assert metadata["dropped_nonfinite_statistic_rows"] == [1, 0]
