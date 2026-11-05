"""An rg pair with a nonpositive SNP heritability estimate is unestimable and keeps its fitted estimates."""

import json
import logging

import numpy as np
import pandas as pd

from ldsc import cli
from ldsc._ldscore.calculator import LDScoreResult
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc._regression import estimators as reg
from ldsc._regression import runner as workflow
from ldsc._sumstats.munger import SumstatsMunger, SumstatsTable
from ldsc.config import GlobalConfig, RegressionConfig

N_SNPS = 3000
N_SAMPLE = 20_000.0
N_BLOCKS = 50
RNG = np.random.default_rng(0)
LD = RNG.uniform(1.0, 50.0, N_SNPS)
# Expected chi-square: trait T1 has SNP heritability 0.3. T2's chi-square
# falls with LD score (slope -0.01, heritability about -0.0015, several
# standard errors below 0). T3 has SNP heritability 0.2.
Z = {
    "T1": np.sqrt(1.0 + N_SAMPLE * 0.3 * LD / N_SNPS) * RNG.standard_normal(N_SNPS),
    "T2": np.sqrt(1.5 - 0.01 * LD) * RNG.standard_normal(N_SNPS),
    "T3": np.sqrt(1.0 + N_SAMPLE * 0.2 * LD / N_SNPS) * RNG.standard_normal(N_SNPS),
}
SNPS = [f"rs{i}" for i in range(N_SNPS)]


def panel():
    frame = pd.DataFrame({
        "CHR": "1", "POS": np.arange(1, N_SNPS + 1) * 100, "SNP": SNPS, "A1": "A", "A2": "C",
        "regression_ld_scores": LD, "base": LD,
    })
    return LDScoreResult(
        frame, None, [{"group": "baseline", "column": "base", "all_reference_snp_count": float(N_SNPS),
                       "common_reference_snp_count": float(N_SNPS)}],
        ["base"], [], config_snapshot=GlobalConfig(snp_identifier="rsid"),
    )


def trait(name):
    data = pd.DataFrame({"SNP": SNPS, "A1": "A", "A2": "C", "Z": Z[name], "N": N_SAMPLE})
    return SumstatsTable(data=data, has_alleles=True, source_path=None, trait_name=name, provenance={},
                         config_snapshot=GlobalConfig(snp_identifier="rsid"))


def kernel_fit(name_1, name_2):
    """Fit the kernel directly on the pair's arrays: one LD-score column, which is also the weight."""
    column = lambda values: np.asarray(values, dtype=float).reshape((-1, 1))
    n = column(np.full(N_SNPS, N_SAMPLE))
    return reg.RG(column(Z[name_1]), column(Z[name_2]), column(LD), column(LD), n, n,
                  np.array([[float(N_SNPS)]]), n_blocks=N_BLOCKS, twostep=30)


def test_a_pair_with_a_negative_heritability_is_unestimable_and_keeps_its_fit(caplog):
    runner = workflow.RegressionRunner(GlobalConfig(snp_identifier="rsid"), RegressionConfig(num_blocks=N_BLOCKS))

    with caplog.at_level(logging.WARNING, logger="LDSC.regression_runner"):
        result = runner.estimate_rg_pairs([trait("T1"), trait("T2")], panel())

    h2 = result.h2_per_trait.set_index("trait_name")
    assert h2.loc["T2", "total_h2_obs"] < 0  # the simulation's premise
    kernel = kernel_fit("T1", "T2")
    assert kernel.hsq2.tot < 0 < kernel.hsq1.tot
    full = result.rg_full.iloc[0]
    assert full["status"] == "unestimable"
    assert full["error"].startswith("SNP heritability estimate not positive:")
    assert "h2_2_obs" in full["error"] and "'T2'" in full["error"] and "h2_1_obs" not in full["error"]
    assert "LDSCInternalError" not in full["error"]
    assert full["n_snps_used"] == N_SNPS
    # The fit's heritability, genetic covariance and intercept estimates are kept.
    expected = {
        "h2_1_obs": kernel.hsq1.tot, "h2_1_obs_se": kernel.hsq1.tot_se,
        "h2_2_obs": kernel.hsq2.tot, "h2_2_obs_se": kernel.hsq2.tot_se,
        "gencov_obs": kernel.gencov.tot, "gencov_obs_se": kernel.gencov.tot_se,
        "intercept_h2_2": kernel.hsq2.intercept, "intercept_gencov": kernel.gencov.intercept,
    }
    for column, value in expected.items():
        np.testing.assert_allclose(full[column], float(np.ravel(value)[0]), rtol=1e-12, err_msg=column)
    for column in ("rg", "rg_se", "z", "p"):
        assert np.isnan(full[column]), column
    concise = result.rg.iloc[0]
    assert np.isnan(concise["rg"]) and np.isnan(concise["rg_se"]) and np.isnan(concise["p"])
    assert concise["n_snps_used"] == N_SNPS
    assert concise["note"] == workflow.H2_OUT_OF_BOUNDS_RG_NOTE
    metadata = result.per_pair_metadata[0]
    assert metadata["status"] == "unestimable"
    assert (metadata["n_snps_used"], metadata["n_blocks_used"]) == (N_SNPS, N_BLOCKS)
    assert metadata["effective_snp_identifier"] == "rsid"
    assert metadata["error"] == full["error"]
    assert any("rg for pair 'T1' vs 'T2' marked unestimable: SNP heritability estimate not positive" in record.getMessage()
               for record in caplog.records)


def test_the_rg_command_counts_the_pair_as_unestimable_and_fits_the_others(tmp_path):
    LDScoreDirectoryWriter().write_batches([panel()], LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    sources = [SumstatsMunger().write_output(trait(name), tmp_path / "sumstats") for name in ("T1", "T2", "T3")]

    code = cli.run_cli(["rg", "--sumstats-sources", *sources, "--ldscore-dir", str(tmp_path / "ld"),
                        "--output-dir", str(tmp_path / "out"), "--num-blocks", str(N_BLOCKS),
                        "--write-per-pair-detail"])

    assert code == 0
    out = tmp_path / "out"
    full = pd.read_csv(out / "rg_full.tsv", sep="\t", keep_default_na=False)
    assert full[["trait_1", "trait_2"]].values.tolist() == [["T1", "T2"], ["T1", "T3"], ["T2", "T3"]]
    assert full["status"].tolist() == ["unestimable", "ok", "unestimable"]
    assert full.loc[[0, 2], "rg"].tolist() == ["NaN", "NaN"]
    assert float(full.loc[1, "rg"]) == float(full.loc[1, "rg"])  # a number
    assert full.loc[0, "error"].startswith("SNP heritability estimate not positive:")
    assert float(full.loc[0, "h2_2_obs"]) < 0 < float(full.loc[0, "h2_1_obs"])
    assert float(full.loc[2, "h2_1_obs"]) < 0 < float(full.loc[2, "h2_2_obs"])
    concise = pd.read_csv(out / "rg.tsv", sep="\t", keep_default_na=False)
    assert concise["note"].tolist() == [workflow.H2_OUT_OF_BOUNDS_RG_NOTE, "", workflow.H2_OUT_OF_BOUNDS_RG_NOTE]
    assert concise["n_snps_used"].tolist() == [N_SNPS] * 3
    manifest = pd.read_csv(out / "diagnostics/pairs/manifest.tsv", sep="\t")
    pair = json.loads((out / manifest.loc[0, "metadata_path"]).read_text())
    assert (pair["status"], pair["n_snps_used"], pair["n_blocks_used"]) == ("unestimable", N_SNPS, N_BLOCKS)
    log = (out / "diagnostics/rg.log").read_text()
    assert "ok=1; unestimable=2" in log
    assert "LDSCInternalError" not in log
