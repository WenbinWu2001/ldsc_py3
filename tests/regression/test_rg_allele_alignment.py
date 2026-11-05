"""rg drops the SNPs it cannot align between two traits and records why.

Two kinds of SNP cannot be aligned. An A/T or C/G SNP is one SNP in both
traits, but its alleles are each other's strand complements, so the sign of
its statistic cannot be aligned: an expected loss of the method. A SNP whose
two traits list different allele pairs is usually an input error; it can only
reach alignment in a base mode, where SNP identity ignores alleles.
"""

import gzip
import json

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, load_sumstats
from ldsc.cli import run_cli
from ldsc._ldscore.outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ldsc._regression.runner import RegressionRunner
from ldsc.config import RegressionConfig
from tests.regression.test_trait_preparation import N_SNPS, panel, signed_statistics

STRAND_UNKNOWN_ROWS = (0, 1, 2)
INCOMPATIBLE_ROWS = (5, 6)


def write_legacy_trait(path, seed, *, second_trait):
    rng = np.random.default_rng(seed)
    z = signed_statistics(rng, panel("chr_pos").baseline_table["base"])
    rows = []
    for i in range(N_SNPS):
        a1, a2 = ("A", "T") if i in STRAND_UNKNOWN_ROWS else ("A", "C")
        if second_trait and i in INCOMPATIBLE_ROWS:
            a1, a2 = "A", "G"
        rows.append(f"rs{i}\t{a1}\t{a2}\t{z[i]:.5f}\t10000")
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("SNP\tA1\tA2\tZ\tN\n" + "\n".join(rows) + "\n")


def test_rg_records_and_reports_snps_it_cannot_align(tmp_path, capsys):
    LDScoreDirectoryWriter().write_batches([panel("chr_pos")], LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    write_legacy_trait(tmp_path / "a.sumstats.gz", 21, second_trait=False)
    write_legacy_trait(tmp_path / "b.sumstats.gz", 22, second_trait=True)

    status = run_cli(["rg", "--sumstats-sources", str(tmp_path / "a.sumstats.gz"), str(tmp_path / "b.sumstats.gz"),
                      "--ldscore-dir", str(tmp_path / "ld"), "--output-dir", str(tmp_path / "out"),
                      "--num-blocks", "10"])

    assert status == 0
    metadata = json.loads((tmp_path / "out/diagnostics/metadata.json").read_text())
    assert metadata["dropped_orientation_unknown_rows"] == [3]
    assert metadata["dropped_incompatible_allele_rows"] == [2]
    rg_full = (tmp_path / "out/rg_full.tsv").read_text().splitlines()
    assert rg_full[1].split("\t")[2] == str(N_SNPS - 5)
    stderr = capsys.readouterr().err
    assert "WARNING" in stderr
    assert "3 A/T or C/G SNPs" in stderr and "2 SNPs" in stderr
    assert "diagnostics/metadata.json" in stderr


def test_rg_records_zero_when_every_snp_aligns(tmp_path, capsys):
    LDScoreDirectoryWriter().write_batches([panel("chr_pos")], LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    for name, seed in (("a", 31), ("b", 32)):
        rng = np.random.default_rng(seed)
        z = signed_statistics(rng, panel("chr_pos").baseline_table["base"])
        with gzip.open(tmp_path / f"{name}.sumstats.gz", "wt", encoding="utf-8") as handle:
            handle.write("SNP\tA1\tA2\tZ\tN\n" + "\n".join(f"rs{i}\tA\tC\t{z[i]:.5f}\t10000" for i in range(N_SNPS)))

    status = run_cli(["rg", "--sumstats-sources", str(tmp_path / "a.sumstats.gz"), str(tmp_path / "b.sumstats.gz"),
                      "--ldscore-dir", str(tmp_path / "ld"), "--output-dir", str(tmp_path / "out"),
                      "--num-blocks", "10"])

    assert status == 0
    metadata = json.loads((tmp_path / "out/diagnostics/metadata.json").read_text())
    assert metadata["dropped_orientation_unknown_rows"] == [0]
    assert metadata["dropped_incompatible_allele_rows"] == [0]
    assert "could not align" not in capsys.readouterr().err


def test_python_rg_records_the_counts_without_printing(tmp_path, capsys):
    from ldsc._regression.runner import run_rg

    LDScoreDirectoryWriter().write_batches([panel("chr_pos")], LDScoreOutputConfig(output_dir=tmp_path / "ld"))
    write_legacy_trait(tmp_path / "a.sumstats.gz", 21, second_trait=False)
    write_legacy_trait(tmp_path / "b.sumstats.gz", 22, second_trait=True)

    result = run_rg(sumstats_sources=[str(tmp_path / "a.sumstats.gz"), str(tmp_path / "b.sumstats.gz")],
                    ldscore_dir=str(tmp_path / "ld"), output_dir=str(tmp_path / "out"), num_blocks=10)

    assert result.dropped_orientation_unknown_rows == [3]
    assert result.dropped_incompatible_allele_rows == [2]
    assert capsys.readouterr().err == ""


def test_rg_keeps_every_snp_from_both_outputs_of_a_lowercase_base_mode_munge(tmp_path):
    # Both traits list the same pair a/c for every SNP, so rg can align all of them.
    ld_scores = panel("rsid")
    loaded = {"parquet": [], "sumstats.gz": []}
    for name, seed in (("t1", 41), ("t2", 42)):
        # Each magnitude appears with both signs, so the median Z, which munging checks, is 0.
        magnitude = np.abs(signed_statistics(np.random.default_rng(seed), ld_scores.baseline_table["base"]))
        z = np.repeat(magnitude[: N_SNPS // 2], 2) * np.tile([1, -1], N_SNPS // 2)
        raw = pd.DataFrame({"SNP": [f"rs{i}" for i in range(N_SNPS)], "A1": "a", "A2": "c",
                            "P": 2 * norm.sf(np.abs(z)), "Z": z, "N": 10_000})
        raw.to_csv(tmp_path / f"{name}.txt", sep="\t", index=False)
        SumstatsMunger().run(
            MungeConfig(raw_sumstats_file=tmp_path / f"{name}.txt", output_dir=tmp_path / name, trait_name=name,
                        no_snp_restriction=True, output_format="both"),
            global_config=GlobalConfig(snp_identifier="rsid"),
        )
        for suffix in loaded:
            loaded[suffix].append(load_sumstats(tmp_path / name / f"{name}.{suffix}"))
    runner = RegressionRunner(global_config=GlobalConfig(snp_identifier="rsid"),
                              regression_config=RegressionConfig(num_blocks=10))

    datasets = {suffix: runner.build_rg_dataset(*tables, ld_scores) for suffix, tables in loaded.items()}

    for dataset in datasets.values():
        assert dataset.dropped_orientation_unknown_rows == 0
        assert dataset.dropped_incompatible_allele_rows == 0
        assert sorted(dataset.merged["SNP"]) == sorted(f"rs{i}" for i in range(N_SNPS))
    assert datasets["parquet"].merged["SNP"].tolist() == datasets["sumstats.gz"].merged["SNP"].tolist()


@pytest.mark.parametrize("mode", ["rsid", "chr_pos"])
def test_rg_drops_missing_alleles_of_base_mode_munge_tables_passed_in_memory(tmp_path, mode):
    # Every seventh SNP has no A1 in either trait, so rg cannot align it; every other SNP lists A/C in both.
    ld_scores = panel(mode)
    missing = np.arange(N_SNPS) % 7 == 0
    coordinates = {} if mode == "rsid" else {"CHR": "1", "POS": ld_scores.baseline_table["POS"]}
    builds = {} if mode == "rsid" else {"source_genome_build": "hg38", "output_genome_build": "hg38"}
    global_config = GlobalConfig(snp_identifier=mode, **({} if mode == "rsid" else {"genome_build": "hg38"}))
    tables = []
    for name, seed in (("t1", 51), ("t2", 52)):
        # Each magnitude appears with both signs, so the median Z, which munging checks, is 0.
        magnitude = np.abs(signed_statistics(np.random.default_rng(seed), ld_scores.baseline_table["base"]))
        z = np.repeat(magnitude[: N_SNPS // 2], 2) * np.tile([1, -1], N_SNPS // 2)
        raw = pd.DataFrame({**coordinates, "SNP": [f"rs{i}" for i in range(N_SNPS)],
                            "A1": np.where(missing, "NA", "A"), "A2": "C",
                            "P": 2 * norm.sf(np.abs(z)), "Z": z, "N": 10_000})
        raw.to_csv(tmp_path / f"{name}.txt", sep="\t", index=False)
        tables.append(SumstatsMunger().run(
            MungeConfig(raw_sumstats_file=tmp_path / f"{name}.txt", output_dir=tmp_path / name, trait_name=name,
                        no_snp_restriction=True, **builds),
            global_config=global_config,
        ))
    runner = RegressionRunner(global_config=global_config, regression_config=RegressionConfig(num_blocks=10))

    result = runner.estimate_rg_pairs(tables, ld_scores)

    assert result.rg_full["status"].tolist() == ["ok"]
    assert result.rg_full["n_snps_used"].tolist() == [int((~missing).sum())]
