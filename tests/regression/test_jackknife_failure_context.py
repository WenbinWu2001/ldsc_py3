"""Singular jackknife deletions name their component, block, genomic span and unsupported columns."""

import json

import numpy as np
import pandas as pd
import pytest

from ldsc import cli
from ldsc._ldscore.calculator import LDScoreResult
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc._regression import runner as workflow
from ldsc._regression.jackknife import JackknifeIdentifiabilityError
from ldsc._sumstats.munger import SumstatsMunger, SumstatsTable
from ldsc.config import GlobalConfig, RegressionConfig

N_SNPS = 200
# The "confined" annotation is supported in two regions. Each trait below keeps
# support in both, but traits A and B share only the first, which lies inside
# one jackknife block of their pair.
FIRST_REGION, SECOND_A, SECOND_B = range(10, 16), range(150, 156), range(156, 162)


def panel():
    rng = np.random.default_rng(5)
    frame = pd.DataFrame({
        "CHR": "1", "POS": np.arange(1, N_SNPS + 1) * 1000, "SNP": [f"rs{i}" for i in range(N_SNPS)],
        "A1": "A", "A2": "C", "regression_ld_scores": rng.uniform(1.0, 3.0, N_SNPS),
    })
    frame["base"] = frame["regression_ld_scores"] * rng.uniform(0.9, 1.1, N_SNPS)
    frame["confined"] = 0.0
    frame.loc[[*FIRST_REGION, *SECOND_A, *SECOND_B], "confined"] = rng.uniform(1.0, 2.0, 18)
    counts = (("base", 1000.0), ("confined", 100.0))
    return LDScoreResult(
        frame, None, [{"group": "baseline", "column": column, "all_reference_snp_count": count,
                       "common_reference_snp_count": count} for column, count in counts],
        ["base", "confined"], [], frozenset(), frozenset(frame.SNP), [],
        config_snapshot=GlobalConfig(snp_identifier="rsid"),
    )


def trait(name, rows, seed, *, coordinates=True):
    rng = np.random.default_rng(seed)
    source = panel().baseline_table.iloc[list(rows)]
    data = source[["SNP", "A1", "A2", *(["CHR", "POS"] if coordinates else [])]].reset_index(drop=True)
    data["Z"] = rng.standard_normal(len(data)) * np.sqrt(1.0 + 20_000 * 0.2 * source["base"].to_numpy() / 1000)
    data["N"] = 20_000.0
    return SumstatsTable(data=data, has_alleles=True, source_path=None, trait_name=name, provenance={},
                         config_snapshot=GlobalConfig(snp_identifier="rsid"))


def traits():
    every = range(N_SNPS)
    a = [row for row in every if row not in SECOND_B]
    b = [row for row in every if row not in SECOND_A]
    # Trait A has no coordinates, so a genomic span can only come from the panel.
    return trait("A", a, 1, coordinates=False), trait("B", b, 2), trait("C", every, 3)


def only_first_region(seed=4):
    return trait("D", [row for row in range(N_SNPS) if row not in (*SECOND_A, *SECOND_B)], seed)


def test_rg_marks_a_singular_pair_unestimable_with_its_context(tmp_path):
    LDScoreDirectoryWriter().write(panel(), LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    sources = [SumstatsMunger().write_output(table, tmp_path / "sumstats") for table in traits()]

    result = workflow.run_rg(sumstats_sources=sources, ldscore_dir=str(tmp_path / "ld"),
                             output_dir=str(tmp_path / "out"), num_blocks=10, write_per_pair_detail=True)

    out = tmp_path / "out"
    full = pd.read_csv(out / "rg_full.tsv", sep="\t", keep_default_na=False)
    assert full["status"].tolist() == ["unestimable", "ok", "ok"]
    error = full.loc[0, "error"]
    assert "trait 1 heritability" in error
    assert "no remaining support in confined" in error
    assert "1:1000-" in error  # panel coordinates of the deleted block, although trait A has none
    shared_snps = N_SNPS - len(SECOND_A) - len(SECOND_B)
    assert int(full.loc[0, "n_snps_used"]) == shared_snps
    concise = pd.read_csv(out / "rg.tsv", sep="\t", keep_default_na=False)
    assert concise.loc[0, "note"] == workflow.UNESTIMABLE_RG_NOTE
    assert not [column for column in [*full.columns, *concise.columns] if column.startswith("_ldsc")]
    log = (out / "diagnostics/rg.log").read_text()
    assert "'A' vs 'B' jackknife diagnostics" in log
    assert '"component": "trait_1_h2"' in log
    assert '"confined": {"full": 6, "remaining": 0}' in log
    manifest = pd.read_csv(out / "diagnostics/pairs/manifest.tsv", sep="\t")
    pair = json.loads((out / manifest.loc[0, "metadata_path"]).read_text())
    assert pair["status"] == "unestimable"
    assert pair["n_snps_used"] == shared_snps
    assert pair["n_blocks_used"] == 10
    assert pair["retained_ld_columns"] == ["base", "confined"]


def test_h2_command_prints_the_jackknife_context(tmp_path, caplog):
    LDScoreDirectoryWriter().write(panel(), LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    source = SumstatsMunger().write_output(only_first_region(), tmp_path / "sumstats")

    code = cli.run_cli(["h2", "--sumstats-file", source, "--ldscore-dir", str(tmp_path / "ld"),
                        "--output-dir", str(tmp_path / "out"), "--num-blocks", "10"])

    assert code == 1
    printed = [record.getMessage() for record in caplog.records if record.getMessage().startswith("Error:")]
    assert len(printed) == 1
    assert "Jackknife regression is unestimable" in printed[0]
    assert "no remaining support in confined" in printed[0]
    assert "1:1000-" in printed[0]


def test_context_survives_the_calculation_child():
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))
    source = panel()

    with pytest.raises(JackknifeIdentifiabilityError) as h2_error:
        runner.estimate_h2(runner.build_dataset(only_first_region(), source))
    a, b, _ = traits()
    with pytest.raises(JackknifeIdentifiabilityError) as rg_error:
        runner.estimate_rg(runner.build_rg_dataset(a, b, source))

    assert "no remaining support in confined" in str(h2_error.value)
    assert "1:1000-" in str(h2_error.value)
    assert rg_error.value.failures[0]["component"] == "trait_1_h2"
    assert "trait 1 heritability" in str(rg_error.value)
    assert rg_error.value.fit_summary["n_snps"] == N_SNPS - len(SECOND_A) - len(SECOND_B)
    assert rg_error.value.fit_summary["n_blocks"] == 10


def test_message_lists_three_failed_blocks_then_counts_the_rest():
    failures = [{"block_index": index, "rank": 1, "n_parameters": 2, "component": "genetic_covariance",
                 "zero_columns": ["q"], "genomic_spans": [{"CHR": "2", "POS_start": index, "POS_end": index + 5}]}
                for index in range(5)]
    error = JackknifeIdentifiabilityError(failures, 10, separators=np.arange(0, 110, 10))

    message = str(error)

    assert message.count("genetic covariance") == 3
    assert "block 3/10 (zero-based index 2, rows [20, 30)) of the genetic covariance fit" in message
    assert "no remaining support in q" in message
    assert "genomic span 2:2-7" in message
    assert "block 4/10" not in message
    assert "and 2 more failed blocks" in message
