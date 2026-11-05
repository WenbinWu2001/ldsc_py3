"""The munging output family and its publication through SumstatsDirectoryWriter."""

import gzip

import pandas as pd
import pytest

from ldsc import GlobalConfig, LDSCInputError, SumstatsMunger, SumstatsTable
from ldsc._sumstats.outputs import SumstatsDirectoryWriter


def _table(chrom="1"):
    return SumstatsTable(
        data=pd.DataFrame({"SNP": ["rs1", "rs2"], "CHR": [chrom, "1"], "POS": [10, 20],
                           "A1": ["A", "C"], "A2": ["G", "T"], "Z": [1.0, 2.0], "N": [100.0, 100.0]}),
        has_alleles=True, source_path=None, trait_name="trait", config_snapshot=GlobalConfig(snp_identifier="rsid"),
    )


def test_a_run_owns_both_formats_and_its_diagnostics(tmp_path):
    family = SumstatsDirectoryWriter.artifact_family(tmp_path, trait_name="my trait", output_format="tsv.gz", run=True)

    assert {name: path.relative_to(tmp_path).as_posix() for name, path in family.paths.items()} == {
        "tsv.gz": "my_trait.sumstats.gz",
        "log": "diagnostics/sumstats.log",
        "dropped_snps": "diagnostics/dropped_snps/dropped.tsv.gz",
    }
    assert [path.relative_to(tmp_path).as_posix() for path in family.owned_paths] == [
        "my_trait.parquet", "my_trait.sumstats.gz", "diagnostics/sumstats.log",
        "diagnostics/dropped_snps/dropped.tsv.gz", "diagnostics/input_issues.tsv",
    ]


def test_write_output_owns_only_the_data_files(tmp_path):
    family = SumstatsDirectoryWriter.artifact_family(tmp_path, output_format="parquet")

    assert [path.name for path in family.owned_paths] == ["sumstats.parquet", "sumstats.gz"]
    dropped = tmp_path / "diagnostics" / "dropped_snps" / "dropped.tsv.gz"
    dropped.parent.mkdir(parents=True)
    dropped.write_bytes(gzip.compress(b"CHR\n"))
    (tmp_path / "trait.sumstats.gz").write_text("stale")

    primary = SumstatsMunger().write_output(_table(), tmp_path, overwrite=True)

    assert primary == str(tmp_path / "trait.parquet")
    assert not (tmp_path / "trait.sumstats.gz").exists()
    assert dropped.exists()


def test_an_unsupported_contig_is_rejected_before_any_file(tmp_path):
    # The curated frame normalizes CHR before the writer stages anything.
    with pytest.raises(LDSCInputError, match="Un_gl000220"):
        SumstatsMunger().write_output(_table(chrom="chrUn_gl000220"), tmp_path, output_format="both")

    assert sorted(path.name for path in tmp_path.iterdir()) == []


def test_a_failed_write_publishes_no_file(tmp_path, monkeypatch):
    from ldsc._sumstats import outputs

    def failing_parquet(*args, **kwargs):
        raise OSError("disk full")

    # The text file is complete when the Parquet write fails.
    monkeypatch.setattr(outputs, "write_sumstats_parquet", failing_parquet)
    with pytest.raises(OSError, match="disk full"):
        SumstatsMunger().write_output(_table(), tmp_path, output_format="both")

    assert sorted(path.name for path in tmp_path.iterdir()) == []
