"""Public output naming and collision contracts for munged sumstats."""

from dataclasses import replace

import pytest

from ldsc import SumstatsMunger, load_sumstats
from ldsc.sumstats_munger import main


@pytest.mark.parametrize(
    "trait_name,parquet_name,gzip_name",
    [
        ("BMI", "BMI.parquet", "BMI.sumstats.gz"),
        (None, "sumstats.parquet", "sumstats.gz"),
        ("../BMI / adult:*?", "BMI_adult.parquet", "BMI_adult.sumstats.gz"),
        ("...", "trait.parquet", "trait.sumstats.gz"),
    ],
)
def test_cli_and_python_writer_use_safe_trait_filenames(
    tmp_path, trait_name, parquet_name, gzip_name
):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP A1 A2 P BETA N\nrs1 A G .05 .1 1000\n")
    output = tmp_path / "munged"
    args = [
        "--raw-sumstats-file", str(raw), "--output-dir", str(output),
        "--snp-identifier", "rsid", "--no-snp-restriction", "--output-format", "both",
    ]
    if trait_name is not None:
        args += ["--trait-name", trait_name]
    table = main(args)
    for destination in (output, tmp_path / "written"):
        if destination != output:
            path = SumstatsMunger().write_output(table, destination, output_format="both")
            assert path == str(destination / parquet_name)
        assert {p.name for p in destination.iterdir() if p.is_file()} == {
            parquet_name, gzip_name,
        }
        loaded = load_sumstats(destination / parquet_name)
        assert loaded.trait_name == (trait_name or parquet_name)
        assert loaded.data.SNP.tolist() == ["rs1"]
        assert load_sumstats(destination / gzip_name).data.SNP.tolist() == ["rs1"]


def test_sanitized_name_collision_requires_overwrite_and_cleans_only_sibling(tmp_path):
    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP A1 A2 P BETA N\nrs1 A G .05 .1 1000\n")
    table = main([
        "--raw-sumstats-file", str(raw), "--output-dir", str(tmp_path / "input"),
        "--snp-identifier", "rsid", "--no-snp-restriction", "--trait-name", "BMI/adult",
    ])
    output = tmp_path / "written"
    munger = SumstatsMunger()
    munger.write_output(table, output, output_format="both")
    unrelated = output / "other.parquet"
    unrelated.write_text("unrelated")
    table = replace(table, trait_name="BMI:adult")
    with pytest.raises(FileExistsError, match="overwrite"):
        munger.write_output(table, output)
    path = munger.write_output(table, output, overwrite=True)
    assert path == str(output / "BMI_adult.parquet")
    assert not (output / "BMI_adult.sumstats.gz").exists()
    assert unrelated.read_text() == "unrelated"
    assert load_sumstats(path).trait_name == "BMI:adult"
