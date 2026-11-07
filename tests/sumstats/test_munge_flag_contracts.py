"""Public munging contracts for format selection and retained input values."""

import pandas as pd
import pytest

from ldsc import LDSCInputError, load_sumstats
from ldsc._sumstats.munger import main


@pytest.mark.parametrize("frequency_column", ["FRQ", "MAF", None])
def test_frequency_is_preserved_automatically_after_maf_filtering(tmp_path, frequency_column):
    raw = tmp_path / "raw.tsv"
    frame = pd.DataFrame({
        "SNP": ["rs1", "rs2", "rs3"], "P": [.05, .05, .05],
        "BETA": [-.05, 0, .05], "N": [1000, 1000, 1000],
    })
    if frequency_column:
        frame[frequency_column] = [.8123456789, .005, .01]
    frame.to_csv(raw, sep="\t", index=False)
    output = tmp_path / "out"
    table = main([
        "--raw-sumstats-file", str(raw), "--output-dir", str(output),
        "--snp-identifier", "rsid", "--no-snp-restriction", "--output-format", "both",
    ])
    for actual in (table, load_sumstats(output / "sumstats.parquet"),
                   load_sumstats(output / "sumstats.gz")):
        if frequency_column:
            assert actual.data.SNP.tolist() == ["rs1", "rs3"]
            assert actual.data.FRQ.tolist() == pytest.approx([.8123456789, .01], abs=1e-12, rel=0)
        else:
            assert actual.data.SNP.tolist() == ["rs1", "rs2", "rs3"]
            assert "FRQ" not in actual.data


@pytest.mark.parametrize("count_columns", [("Nca", "Nco"), ("NCAS", "NCON"), ("nca", "nco")])
@pytest.mark.parametrize("frequency_column", [None, "FRQ_U_60"])
def test_daner_auto_and_explicit_share_aliases_and_optional_frequency(tmp_path, count_columns, frequency_column):
    raw = tmp_path / "daner.tsv"
    cases, controls = count_columns
    frame = pd.DataFrame({
        "SNP": ["rs1", "rs2"], "P": [.05, .05], "BETA": [-.05, .05],
        cases: [40, 60], controls: [60, 40],
    })
    if frequency_column:
        frame[frequency_column] = [.8, .2]
    frame.to_csv(raw, sep="\t", index=False)
    outputs = []
    for profile in ("auto", "daner-new"):
        table = main([
            "--raw-sumstats-file", str(raw), "--output-dir", str(tmp_path / profile),
            "--snp-identifier", "rsid", "--no-snp-restriction", "--input-format", profile,
        ])
        assert table.data.N.tolist() == [80, 120]
        if frequency_column:
            assert table.data.FRQ.tolist() == pytest.approx([.8, .2])
        else:
            assert "FRQ" not in table.data
        outputs.append(table.data)
    pd.testing.assert_frame_equal(*outputs)


def test_daner_auto_and_explicit_share_missing_sample_size_validation(tmp_path):
    raw = tmp_path / "incomplete.tsv"
    raw.write_text("SNP P BETA Nca\nrs1 .05 0 40\n")
    for profile in ("auto", "daner-new"):
        with pytest.raises(LDSCInputError, match="could not determine sample size"):
            main([
                "--raw-sumstats-file", str(raw), "--output-dir", str(tmp_path / profile),
                "--snp-identifier", "rsid", "--no-snp-restriction", "--input-format", profile,
            ])


@pytest.mark.parametrize("options", [[], ["--N-min", "0"]])
def test_legacy_n_threshold_and_constant_precedence_are_preserved(tmp_path, options):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP P BETA N\nrs1 .05 -.05 550\nrs2 .05 .05 1000\n")
    table = main([
        "--raw-sumstats-file", str(raw), "--output-dir", str(tmp_path / "out"),
        "--snp-identifier", "rsid", "--no-snp-restriction", "--N", "9999", *options,
    ])
    # N90=955: dividing by 1.5 removes rs1; dividing by 2 would retain it.
    assert table.data.SNP.tolist() == ["rs2"]
    assert table.data.N.tolist() == [1000]


@pytest.mark.parametrize("sample_size", [["--N", "1000"], ["--N-cas", "400", "--N-con", "600"]])
def test_constant_n_bypasses_n_min_as_in_legacy(tmp_path, sample_size):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP P BETA\nrs1 .05 0\n")
    table = main([
        "--raw-sumstats-file", str(raw), "--output-dir", str(tmp_path / "out"),
        "--snp-identifier", "rsid", "--no-snp-restriction", "--N-min", "2000", "--chunksize", "1", *sample_size,
    ])
    assert table.data.N.tolist() == [1000]


_BUILD_CONTRACT_CASES = {
    'rsID with a source build': (['--snp-identifier', 'rsid', '--source-genome-build', 'hg19'],
                                 {'source_genome_build': 'hg19'}, 'rsid',
                                 '`--source-genome-build` (`MungeConfig.source_genome_build`)'),
    'rsID with an output build': (['--snp-identifier', 'rsid', '--output-genome-build', 'hg19'],
                                  {'output_genome_build': 'hg19'}, 'rsid',
                                  '`--output-genome-build` (`MungeConfig.output_genome_build`)'),
    'rsID with a chain file': (['--snp-identifier', 'rsid', '--liftover-chain-file', 'x.chain'],
                               {'liftover_chain_file': 'x.chain'}, 'rsid',
                               '`--liftover-chain-file` (`MungeConfig.liftover_chain_file`)'),
    'chr_pos without an output build': (['--snp-identifier', 'chr_pos'], {}, 'chr_pos',
                                        '`--output-genome-build` (`MungeConfig.output_genome_build`)'),
}


@pytest.mark.parametrize('case', list(_BUILD_CONTRACT_CASES))
def test_build_options_are_checked_once_with_flag_and_field_names(tmp_path, case):
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger
    from ldsc.errors import LDSCUsageError

    flags, fields, mode, names = _BUILD_CONTRACT_CASES[case]
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP CHR POS A1 A2 P BETA N\nrs1 1 100 A G 0.05 0.1 1000\n')
    base = ['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'), '--no-snp-restriction', *flags]
    messages = []
    for argv in (base, [*base, '--infer-only']):
        with pytest.raises(LDSCUsageError) as caught:
            main(argv)
        messages.append(str(caught.value))
    with pytest.raises(LDSCUsageError) as caught:
        SumstatsMunger().run(MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / 'py', no_snp_restriction=True, **fields),
                             global_config=GlobalConfig(snp_identifier=mode))
    messages.append(str(caught.value))
    assert names in messages[0]
    assert messages[0] == messages[1] == messages[2]


def test_munge_config_rejects_an_unknown_column_hint_key():
    from ldsc.config import MungeConfig
    from ldsc.errors import LDSCConfigError

    with pytest.raises(LDSCConfigError) as failure:
        MungeConfig(column_hints={"N_column": "NEFF", "snp": "rsid"})
    message = str(failure.value)
    assert "'N_column'" in message and "'snp'" not in message.split("Accepted keys")[0]
    assert "N_col" in message and "signed_sumstats" in message and "info_list" in message


def test_munge_config_accepts_every_documented_column_hint_key():
    from ldsc.config import MungeConfig

    keys = ["snp", "chr", "pos", "a1", "a2", "p", "frq", "info", "nstudy", "N_col", "N_cas_col", "N_con_col",
            "N", "N_cas", "N_con", "signed_sumstats", "info_list"]
    assert set(MungeConfig(column_hints={key: "X" for key in keys}).column_hints) == set(keys)
