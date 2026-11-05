"""Each trait's summary statistics are prepared once per regression run."""

import numpy as np
import pandas as pd

from ldsc import cli
from ldsc._ldscore.results import LDScoreResult
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc._regression import alignment, rg, runner as workflow
from ldsc._sumstats.munger import SumstatsMunger
from ldsc._sumstats.table import SumstatsTable
from ldsc.config import GlobalConfig, RegressionConfig
from tests.worker_fixtures import calculation_body

N_SNPS = 160


def panel(mode="chr_pos_allele_aware"):
    rng = np.random.default_rng(3)
    frame = pd.DataFrame({
        "CHR": "1", "POS": np.arange(1, N_SNPS + 1) * 100, "SNP": [f"rs{i}" for i in range(N_SNPS)],
        "A1": "A", "A2": "C", "regression_ld_scores": 1.0 + rng.gamma(2.0, 5.0, N_SNPS),
    })
    frame["base"] = frame["regression_ld_scores"] * (0.9 + 0.2 * rng.random(N_SNPS))
    return LDScoreResult(
        frame, None,
        [{"group": "baseline", "column": "base", "all_reference_snp_count": 1000.0,
          "common_reference_snp_count": 900.0}],
        ["base"], [],
        config_snapshot=GlobalConfig(snp_identifier=mode, genome_build="hg38"),
    )


def munged(name, mode, seed):
    rng = np.random.default_rng(seed)
    ld = panel().baseline_table
    data = ld[["CHR", "POS", "SNP", "A1", "A2"]].copy()
    data["Z"] = signed_statistics(rng, ld["base"])
    data["N"] = 10_000.0
    return SumstatsTable(data=data, has_alleles=True, source_path=f"{name}.parquet", trait_name=name,
                         provenance={}, config_snapshot=GlobalConfig(snp_identifier=mode, genome_build="hg38"))


def legacy(name, seed, *, extra_row):
    rng = np.random.default_rng(seed)
    data = pd.DataFrame({"SNP": [f"rs{i}" for i in range(N_SNPS)], "A1": "A", "A2": "C",
                         "Z": signed_statistics(rng, panel().baseline_table["base"]), "N": 10_000.0})
    data = pd.concat([data, pd.DataFrame([extra_row])], ignore_index=True)
    return SumstatsTable(data=data, has_alleles=True, source_path=f"{name}.sumstats.gz", trait_name=name,
                         provenance={"legacy_ldsc": True}, config_snapshot=None)


def signed_statistics(rng, ld_scores):
    """Draw Z with chi-square mean 1 + N h2 l / M (N = 10,000, h2 = 0.5, M = 900)."""
    return rng.standard_normal(len(ld_scores)) * np.sqrt(1.0 + 10_000 * 0.5 * ld_scores.to_numpy() / 900)


def counting(monkeypatch, name):
    """Record the context (or trait name) of every call to a module-level alignment function."""
    original, calls = getattr(alignment, name), []

    def wrapper(*args, **kwargs):
        calls.append(kwargs["context"] if "context" in kwargs else args[0].trait_name)
        return original(*args, **kwargs)

    monkeypatch.setattr(alignment, name, wrapper)
    return calls


def test_rg_run_projects_each_legacy_trait_once_and_returns_its_drops(monkeypatch):
    projections = counting(monkeypatch, "_project_legacy_sumstats_to_panel")
    tables = [
        legacy("L1", 1, extra_row={"SNP": "rs1", "A1": "A", "A2": "T", "Z": 0.1, "N": 10_000.0}),
        munged("A", "chr_pos_allele_aware", 2),
        legacy("L2", 3, extra_row={"SNP": "rs_absent", "A1": "A", "A2": "C", "Z": 0.1, "N": 10_000.0}),
    ]
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    with calculation_body():
        result = runner.estimate_rg_pairs(tables, panel())

    assert projections == ["L1", "L2"]  # three pairs, two legacy traits, one projection each
    assert result.rg_full["status"].tolist() == ["ok", "ok", "ok"]
    drops = result.legacy_sumstats_drops
    assert drops["trait_name"].tolist() == ["L1", "L1", "L2"]
    assert drops["reason"].tolist() == ["duplicate_source_rsid", "duplicate_source_rsid", "missing_panel_rsid"]


def test_mixed_mode_rg_keys_each_trait_at_most_twice_and_matches_separate_pairs(monkeypatch):
    tables = [munged("A", "chr_pos_allele_aware", 4), munged("B", "chr_pos", 5),
              munged("C", "chr_pos_allele_aware", 6)]
    source = panel()
    config = RegressionConfig(num_blocks=10, allow_identity_downgrade=True)
    runner = workflow.RegressionRunner(regression_config=config)
    fitted = []
    original_fit = rg._fit_rg_dataset

    def record(dataset, config):
        fitted.append(dataset)
        return original_fit(dataset, config)

    monkeypatch.setattr(rg, "_fit_rg_dataset", record)
    keyings = counting(monkeypatch, "_key_regression_rows")

    with calculation_body():
        result = runner.estimate_rg_pairs(tables, source)

    # A and C run allele-aware with each other and in base mode with B; B always runs in base mode.
    assert sorted(keyings) == ["A.parquet", "A.parquet", "B.parquet", "C.parquet", "C.parquet",
                               "LD-score table", "LD-score table"]
    assert [dataset.identity_downgrade_applied for dataset in fitted] == [True, False, True]
    assert result.rg_full["status"].tolist() == ["ok", "ok", "ok"]
    monkeypatch.undo()
    for dataset, (i, j) in zip(fitted, [(0, 1), (0, 2), (1, 2)]):
        with calculation_body():
            separate = runner.build_rg_dataset(tables[i], tables[j], source, config=config)
        pd.testing.assert_frame_equal(dataset.merged, separate.merged)
        assert dataset.effective_snp_identifier == separate.effective_snp_identifier
        assert dataset.retained_ld_columns == separate.retained_ld_columns


def test_base_mode_keying_error_fails_only_its_downgraded_pair(monkeypatch):
    tables = [munged("A", "chr_pos_allele_aware", 7), munged("B", "chr_pos_allele_aware", 8),
              munged("C", "chr_pos", 9)]
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10,
                                                                          allow_identity_downgrade=True))
    original = alignment._PreparedTraitSumstats.keyed_rows

    def keyed_rows(self, mode, *, cleaned):
        # Trait A's base-mode rows are first needed by its downgraded pair with C.
        if self.table.trait_name == "A" and cleaned:
            raise ValueError("injected base-mode keying failure")
        return original(self, mode, cleaned=cleaned)

    monkeypatch.setattr(alignment._PreparedTraitSumstats, "keyed_rows", keyed_rows)

    with calculation_body():
        result = runner.estimate_rg_pairs(tables, panel())

    assert result.rg_full[["trait_1", "trait_2", "status"]].values.tolist() == [
        ["A", "B", "ok"], ["A", "C", "failed"], ["B", "C", "ok"],
    ]
    assert "injected base-mode keying failure" in result.rg_full.loc[1, "error"]


def test_build_dataset_projects_a_legacy_table_once(monkeypatch):
    projections = counting(monkeypatch, "_project_legacy_sumstats_to_panel")
    table = legacy("L1", 10, extra_row={"SNP": "rs_absent", "A1": "A", "A2": "C", "Z": 0.1, "N": 10_000.0})
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    with calculation_body():
        dataset = runner.build_dataset(table, panel())

    assert projections == ["L1"]
    assert len(dataset.merged) == N_SNPS


def identity_lines(log_path):
    return [line for line in log_path.read_text(encoding="utf-8").splitlines() if "SNP identifier" in line
            and "used for this regression" in line]


def run_rg_command(tmp_path, tables, *flags):
    LDScoreDirectoryWriter().write_batches([panel()], LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    sources = [SumstatsMunger().write_output(table, tmp_path / table.trait_name) for table in tables]
    code = cli.run_cli(["rg", "--sumstats-sources", *sources, "--ldscore-dir", str(tmp_path / "ld"),
                        "--output-dir", str(tmp_path / "out"), "--num-blocks", "10", *flags])
    assert code == 0
    return identity_lines(tmp_path / "out/diagnostics/rg.log")


def test_rg_command_logs_the_snp_identifier_each_pair_used(tmp_path):
    # A and B match the allele-aware LD scores; C is allele-unaware, so only the pairs with C are downgraded.
    tables = [munged("A", "chr_pos_allele_aware", 11), munged("B", "chr_pos_allele_aware", 12),
              munged("C", "chr_pos", 13)]

    lines = run_rg_command(tmp_path, tables, "--allow-identity-downgrade")

    # Pairs A-B, A-C and B-C run in chr_pos_allele_aware, chr_pos (downgraded) and chr_pos (downgraded).
    assert len(lines) == 1
    assert ("SNP identifiers used for this regression: chr_pos_allele_aware for 1 pair; "
            "chr_pos, after an allele-awareness downgrade, for 2 pairs.") in lines[0]


def test_rg_command_without_a_downgrade_reports_none(tmp_path):
    tables = [munged(name, "chr_pos_allele_aware", seed) for name, seed in (("A", 14), ("B", 15))]

    lines = run_rg_command(tmp_path, tables, "--allow-identity-downgrade")

    assert len(lines) == 1
    assert "SNP identifier used for this regression: chr_pos_allele_aware. The GlobalConfig" in lines[0]
    assert "downgrade" not in lines[0]
