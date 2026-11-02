"""munge-sumstats allele policy: A1 is the effect allele, and its orientation is never guessed."""
import pandas as pd
import pytest

from ldsc import GlobalConfig, LDSCInputError, MungeConfig, SumstatsMunger, load_sumstats
from ldsc._sumstats import munger as workflow

ROWS = "rs1 A G 0.05 0.1 1000\nrs2 C T 0.5 -0.2 1000\n"
EFFECT_REMEDY = "--a1 <effect-allele column> --a2 <other column>"


def _raw(tmp_path, header, rows=ROWS):
    path = tmp_path / "raw.tsv"
    path.write_text(header + "\n" + rows)
    return path


def _run(path, tmp_path, mode, **options):
    return SumstatsMunger().run(
        MungeConfig(raw_sumstats_file=path, output_dir=tmp_path / "out", no_snp_restriction=True, **options),
        global_config=GlobalConfig(snp_identifier=mode),
    )


def _cli(path, tmp_path, mode, *extra, infer=False):
    argv = ["--raw-sumstats-file", str(path), "--snp-identifier", mode, "--no-snp-restriction",
            "--output-dir", str(tmp_path / ("infer" if infer else "cli")), *extra]
    return workflow.main([*argv, "--infer-only"] if infer else argv)


@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware"])
@pytest.mark.parametrize("header,ambiguous", [
    ("SNP REF ALT P BETA N", ["REF", "ALT"]),
    ("SNP ALLELE1 ALLELE2 P BETA N", ["ALLELE1", "ALLELE2"]),
    ("SNP ALLELE1 ALLELE0 P BETA N", ["ALLELE1", "ALLELE0"]),
    ("SNP Reference_Allele allele_2 P BETA N", ["Reference_Allele", "allele_2"]),
])
def test_orientation_ambiguous_allele_columns_are_never_mapped(tmp_path, mode, header, ambiguous):
    with pytest.raises(LDSCInputError) as caught:
        _run(_raw(tmp_path, header), tmp_path, mode)
    message = str(caught.value)
    assert all(column in message for column in ambiguous)
    assert "A1 (the effect allele, the allele the signed statistic refers to)" in message
    assert EFFECT_REMEDY in message
    assert (f"--ignore {','.join(ambiguous)}" in message) is (mode == "rsid")


def test_explicit_orientation_munges_ambiguous_columns(tmp_path):
    table = _run(_raw(tmp_path, "SNP REF ALT P BETA N"), tmp_path, "rsid_allele_aware",
                 column_hints={"a1": "ALT", "a2": "REF"})
    assert table.data.A1.tolist() == ["G", "T"]
    assert table.data.A2.tolist() == ["A", "C"]


def test_base_mode_munges_without_ignored_allele_columns(tmp_path, capsys):
    table = _cli(_raw(tmp_path, "SNP REF ALT P BETA N"), tmp_path, "rsid", "--ignore", "REF,ALT", "--output-format", "both")
    assert not {"A1", "A2"} & set(table.data.columns)
    assert not {"A1", "A2"} & set(load_sumstats(tmp_path / "cli" / "sumstats.parquet").data.columns)
    assert not {"A1", "A2"} & set(pd.read_csv(tmp_path / "cli" / "sumstats.gz", sep="\t").columns)
    assert "  No A1/A2 columns written: rg cannot use this output." in capsys.readouterr().out


def test_summary_line_appears_only_without_alleles(tmp_path, capsys):
    _cli(_raw(tmp_path, "SNP A1 A2 P BETA N"), tmp_path, "rsid")
    assert "No A1/A2 columns written" not in capsys.readouterr().out


@pytest.mark.parametrize("mode,label", [("rsid", "A1/A2 orientation"), ("rsid_allele_aware", "A1/A2 orientation")])
def test_infer_only_reports_missing_orientation_without_guessing(tmp_path, capsys, mode, label):
    report = _cli(_raw(tmp_path, "SNP REF ALT P BETA N"), tmp_path, mode, infer=True)
    output = capsys.readouterr().out
    assert not report.runnable
    assert label in report.missing_fields
    assert report.column_hints == {}
    assert not {"--a1", "--a2", "--ignore"} & set(report.suggested_args)
    assert "--a1 " not in output.split("Suggested command:")[1]
    assert any(EFFECT_REMEDY in note for note in report.notes)
    assert any("--ignore REF,ALT" in note for note in report.notes) is (mode == "rsid")


@pytest.mark.parametrize("mode", ["rsid", "rsid_allele_aware"])
def test_a_single_allele_column_is_an_error_in_every_mode(tmp_path, mode):
    path = _raw(tmp_path, "SNP EA P BETA N", rows="rs1 A 0.05 0.1 1000\n")
    with pytest.raises(LDSCInputError) as caught:
        _run(path, tmp_path, mode)
    message = str(caught.value)
    assert "EA is read as A1" in message and "no column is read as A2" in message
    assert ("--ignore EA" in message) is (mode == "rsid")
    report = _cli(path, tmp_path, mode, infer=True)
    assert "A1/A2 pair" in report.missing_fields


def test_allele_aware_mode_without_allele_columns_names_the_effect_allele(tmp_path):
    with pytest.raises(LDSCInputError) as caught:
        _run(_raw(tmp_path, "SNP P BETA N", rows="rs1 0.05 0.1 1000\n"), tmp_path, "rsid_allele_aware")
    message = str(caught.value)
    assert EFFECT_REMEDY in message
    assert "REF" not in message


# Effects stay small: the kernel's median check rejects a signed statistic whose median is far from null.
@pytest.mark.parametrize("header,rows,hints", [
    ("SNP INC_ALLELE DEC_ALLELE P BETA N", "rs1 A G 0.05 0.01 1000\nrs2 C T 0.5 0 1000\n", {}),
    ("SNP inc_allele dec_allele P OR N", "rs1 A G 0.05 1.01 1000\nrs2 C T 0.5 1 1000\n", {}),
    ("SNP INC_ALLELE OTHER P Z N", "rs1 A G 0.05 0.05 1000\nrs2 C T 0.5 0.03 1000\n", {"a2": "OTHER"}),
])
def test_increasing_allele_columns_that_agree_with_the_signed_statistic_are_used(tmp_path, header, rows, hints):
    table = _run(_raw(tmp_path, header, rows), tmp_path, "rsid_allele_aware", column_hints=hints)
    assert table.data.A1.tolist() == ["A", "C"]


def test_increasing_allele_columns_that_contradict_the_signed_statistic_stop_the_run(tmp_path):
    rows = "rs1 A G 0.05 0.01 1000\nrs2 C T 0.5 -0.02 1000\nrs3 G T 0.5 -0.04 1000\nrs4 A C 0.5 NA 1000\n"
    with pytest.raises(LDSCInputError) as caught:
        _run(_raw(tmp_path, "SNP INC_ALLELE DEC_ALLELE P BETA N", rows), tmp_path, "rsid_allele_aware")
    message = str(caught.value)
    assert "2 rows" in message
    assert "rs2" in message and "rs3" in message and "rs4" not in message
    assert EFFECT_REMEDY in message


def test_a1_inc_runs_with_increasing_allele_columns(tmp_path):
    path = _raw(tmp_path, "SNP INC_ALLELE DEC_ALLELE P N", rows="rs1 A G 0.05 1000\n")
    table = _run(path, tmp_path, "rsid_allele_aware", a1_inc=True)
    assert table.data.A1.tolist() == ["A"]
    assert (table.data.Z > 0).all()
