"""Restriction and automatic liftover choices through the public munger."""

import pytest

from ldsc import GlobalConfig, LDSCConfigError, LDSCInputError, LDSCUsageError, MungeConfig, SumstatsMunger
from ldsc import sumstats_munger as workflow


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    hm3 = tmp_path / "hm3.tsv"
    hm3.write_text("CHR\thg19_POS\thg38_POS\tSNP\n1\t100\t1000\tkeep\n")
    monkeypatch.setattr(workflow, "packaged_hm3_curated_map_path", lambda: str(hm3))
    raw = tmp_path / "raw.tsv"
    raw.write_text("CHR POS SNP A1 A2 P BETA N\n1 100 keep A G .05 .1 1000\n1 500 outside A G .05 -.1 1000\n")
    return raw


def test_packaged_hm3_is_default_for_cli_and_python(inputs, tmp_path):
    table = workflow.main([
        "--raw-sumstats-file", str(inputs), "--output-dir", str(tmp_path / "cli"),
        "--snp-identifier", "rsid",
    ])
    assert table.data.SNP.tolist() == ["keep"]
    munger = SumstatsMunger()
    table = munger.run(
        MungeConfig(raw_sumstats_file=inputs, output_dir=tmp_path / "api"),
        global_config=GlobalConfig(snp_identifier="rsid"),
    )
    assert table.data.SNP.tolist() == ["keep"]
    assert munger.build_run_summary().drop_counts["sumstats_snps"] == 1


@pytest.mark.parametrize("restriction,expected", [("none", ["keep", "outside"]), ("custom", ["outside"])])
def test_restriction_overrides_replace_hm3(inputs, tmp_path, restriction, expected):
    keep = tmp_path / "custom.tsv"
    keep.write_text("SNP\noutside\n")
    flags = ["--no-snp-restriction"] if restriction == "none" else ["--sumstats-snps-file", str(keep)]
    table = workflow.main([
        "--raw-sumstats-file", str(inputs), "--output-dir", str(tmp_path / "out"),
        "--snp-identifier", "rsid", *flags,
    ])
    assert table.data.SNP.tolist() == expected


@pytest.mark.parametrize("target,chain,position,method", [
    ("hg19", False, 100, None),
    ("hg19", True, 100, None),
    ("hg38", False, 1000, "hm3_curated"),
    ("hg38", True, 200, "chain_file"),
])
def test_automatic_liftover_and_chain_override(inputs, tmp_path, target, chain, position, method):
    if chain and target != "hg19":
        pytest.importorskip("pyliftover")
    chain_path = tmp_path / "map.chain"
    chain_path.write_text("chain 1 chr1 10000 + 0 10000 chr1 20000 + 100 10100 1\n10000\n")
    options = dict(source_genome_build="hg19", output_genome_build=target)
    if chain:
        options["liftover_chain_file"] = chain_path
    table = SumstatsMunger().run(
        MungeConfig(raw_sumstats_file=inputs, output_dir=tmp_path / "out", **options),
        global_config=GlobalConfig(snp_identifier="chr_pos"),
    )
    assert table.data.SNP.tolist() == ["keep"]
    assert table.data.POS.tolist() == [position]
    assert table.config_snapshot.genome_build == target
    assert table.provenance["coordinate_provenance"]["liftover"]["method"] == method


@pytest.mark.parametrize("restriction", ["hm3", "none", "custom"])
def test_infer_only_reports_resolved_policy_without_writes(inputs, tmp_path, restriction, capsys):
    flags = []
    if restriction == "none":
        flags = ["--no-snp-restriction"]
    elif restriction == "custom":
        keep = tmp_path / "custom.tsv"
        keep.write_text("CHR POS\n1 500\n")
        flags = ["--sumstats-snps-file", str(keep)]
    output = tmp_path / "unused"
    args = [
        "--raw-sumstats-file", str(inputs), "--output-dir", str(output),
        "--snp-identifier", "chr_pos", "--source-genome-build", "hg19",
        "--output-genome-build", "hg38", *flags,
    ]
    report = workflow.main([*args, "--infer-only"])
    assert report.runnable == (restriction == "hm3")
    assert report.liftover_method == ("hm3 quick" if restriction == "hm3" else "missing; chain file required")
    assert all(flag in report.suggested_args for flag in flags)
    assert "SNP restriction:" in capsys.readouterr().out
    assert not output.exists()
    if restriction != "hm3":
        with pytest.raises(LDSCUsageError, match="chain"):
            workflow.main(args)


def test_restriction_flags_conflict_and_hm3_switches_are_removed():
    parser = workflow.build_parser()
    base = ["--raw-sumstats-file", "raw", "--output-dir", "out"]
    for flags in (["--no-snp-restriction", "--sumstats-snps-file", "keep"],
                  ["--use-hm3-snps"], ["--use-hm3-quick-liftover"]):
        with pytest.raises(SystemExit):
            parser.parse_args([*base, *flags])
    with pytest.raises(LDSCConfigError, match="mutually exclusive"):
        MungeConfig(no_snp_restriction=True, sumstats_snps_file="keep")


def test_unresolved_source_and_missing_output_never_choose_a_build(inputs, tmp_path):
    for options, error, message in [
        ({}, LDSCUsageError, "output_genome_build"),
        ({"output_genome_build": "hg38"}, LDSCInputError, "source-genome-build"),
    ]:
        with pytest.raises(error, match=message):
            SumstatsMunger().run(
                MungeConfig(raw_sumstats_file=inputs, output_dir=tmp_path / "out", **options),
                global_config=GlobalConfig(snp_identifier="chr_pos"),
            )


def test_real_packaged_hm3_infers_source_before_automatic_mapping(tmp_path):
    from ldsc.genome_build_inference import load_packaged_reference_table

    reference = load_packaged_reference_table().head(300)
    reference = reference.assign(SNP=[f"label_{i}" for i in range(len(reference))])
    raw_frame = reference[["CHR", "hg19_POS", "SNP"]].rename(columns={"hg19_POS": "POS"})
    raw_frame = raw_frame.assign(P=.05, BETA=.1, N=1000)
    raw = tmp_path / "raw.tsv"
    raw_frame.to_csv(raw, sep="\t", index=False)
    table = SumstatsMunger().run(
        MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / "out", output_genome_build="hg38"),
        global_config=GlobalConfig(snp_identifier="chr_pos"),
    )
    assert table.data.POS.tolist() == reference.hg38_POS.tolist()
    assert table.data.SNP.tolist() == reference.SNP.tolist()
    assert table.provenance["coordinate_provenance"]["liftover"]["method"] == "hm3_curated"


@pytest.mark.parametrize("level", ["INFO", "ERROR"])
@pytest.mark.parametrize("mode,target,chain,method", [
    ("chr_pos", "hg38", False, "HM3 quick liftover (automatic; package-bundled HM3 metadata)"),
    ("chr_pos", "hg38", True, "chain file (explicit; HM3 quick liftover disabled)"),
    ("chr_pos", "hg19", True, "none (source and output builds match; supplied chain ignored)"),
    ("rsid", None, False, "none (rsID identity; coordinate conversion not applicable)"),
])
def test_cli_reports_selected_method_and_counts_in_stdout_and_log(
    inputs, tmp_path, capsys, level, mode, target, chain, method,
):
    from ldsc import cli

    flags = []
    if target:
        flags = ["--source-genome-build", "hg19", "--output-genome-build", target]
    if chain:
        pytest.importorskip("pyliftover")
        chain_path = tmp_path / "map.chain"
        chain_path.write_text("chain 1 chr1 10000 + 0 10000 chr1 20000 + 100 10100 1\n10000\n")
        flags += ["--liftover-chain-file", str(chain_path)]
    output = tmp_path / "out"
    assert cli.run_cli([
        "munge-sumstats", "--raw-sumstats-file", str(inputs), "--output-dir", str(output),
        "--snp-identifier", mode, "--log-level", level, *flags,
    ]) == 0
    captured = capsys.readouterr()
    report = captured.out.strip()
    assert "Munge-sumstats summary:" in report
    assert "SNP restriction: packaged HM3 (default)" in report
    assert f"Liftover: {method}" in report
    assert "Rows: 2 input; 1 retained; 1 dropped" in report
    assert "sumstats_snps=1" in report
    assert "liftover=0" in report
    if target == "hg38":
        assert "hg19 -> hg38" in report
        assert "Mapping: 1 input; 1 mapped and retained; 0 dropped" in report
    else:
        assert "Mapping: not performed; 0 liftover drops" in report
    assert report in (output / "diagnostics" / "sumstats.log").read_text()
    assert captured.err == ""


def test_mapping_drop_report_separates_keep_list_and_liftover_removals(inputs, tmp_path, capsys):
    (tmp_path / "hm3.tsv").write_text("CHR\thg19_POS\thg38_POS\tSNP\n1\t100\t1000\tkeep\n1\t200\t2000\tdup\n")
    with inputs.open("a") as stream:
        stream.write("1 200 dup A G .05 .1 1000\n1 200 dup A G .05 .1 1000\n")
    workflow.main([
        "--raw-sumstats-file", str(inputs), "--output-dir", str(tmp_path / "out"),
        "--snp-identifier", "chr_pos", "--source-genome-build", "hg19",
        "--output-genome-build", "hg38",
    ])
    report = capsys.readouterr().out
    assert "Mapping: 3 input; 1 mapped and retained; 2 dropped" in report
    assert "duplicate source=2" in report
    assert "unmapped=0" in report
    assert "Rows: 4 input; 1 retained; 3 dropped" in report
    assert "sumstats_snps=1" in report
    assert "liftover=2" in report


def test_python_api_records_summary_without_console_output(inputs, tmp_path, capsys):
    SumstatsMunger().run(
        MungeConfig(raw_sumstats_file=inputs, output_dir=tmp_path / "api"),
        global_config=GlobalConfig(snp_identifier="rsid", log_level="ERROR"),
    )
    assert capsys.readouterr().out == ""
    log = (tmp_path / "api" / "diagnostics" / "sumstats.log").read_text()
    assert "Munge-sumstats summary:" in log
    assert "Rows: 2 input; 1 retained; 1 dropped" in log
