"""Public munging contracts for format selection and retained input values."""

import pandas as pd
import pytest

from ldsc import LDSCInputError, load_sumstats
from ldsc.sumstats_munger import build_parser, main


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


@pytest.mark.parametrize("removed", [["--format", "auto"], ["--keep-maf"], ["--keep-frequency"]])
def test_removed_munge_flags_are_rejected(removed):
    with pytest.raises(SystemExit) as caught:
        build_parser().parse_args(["--raw-sumstats-file", "unused.tsv", "--output-dir", "unused", *removed])
    assert caught.value.code == 2


@pytest.mark.parametrize("options", [[], ["--n-min", "0"]])
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
        "--snp-identifier", "rsid", "--no-snp-restriction", "--n-min", "2000", "--chunksize", "1", *sample_size,
    ])
    assert table.data.N.tolist() == [1000]
