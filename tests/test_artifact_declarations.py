"""Artifact-family contracts exercised through public workflows and writers."""
from pathlib import Path
import json

import pandas as pd
import pytest

from ldsc.ldscore_calculator import build_parser, run_ldscore_from_args


def plink_arguments(tmp_path):
    pytest.importorskip("bitarray")
    prefix = Path(__file__).parent / "fixtures" / "plink" / "plink"
    bim = pd.read_csv(prefix.with_suffix(".bim"), sep=r"\s+", header=None)
    snps = tmp_path / "snps.tsv"
    pd.DataFrame({"SNP": bim[1]}).to_csv(snps, sep="\t", index=False)
    return [
        "--output-dir", str(tmp_path / "result"), "--plink-prefix", str(prefix),
        "--snp-identifier", "rsid", "--regr-snps-file", str(snps),
        "--regr-snps-exclude-regions", "none", "--ld-wind-snps", "10", "--yes-really",
    ]


def test_existing_chromosome_drop_report_blocks_ldscore_before_log(tmp_path):
    argv = plink_arguments(tmp_path)
    output = tmp_path / "result"
    report = output / "diagnostics" / "dropped_snps" / "chr1_dropped.tsv.gz"
    report.parent.mkdir(parents=True)
    report.write_bytes(b"existing report")

    with pytest.raises(FileExistsError, match="chr1_dropped"):
        run_ldscore_from_args(build_parser().parse_args(argv))

    assert report.read_bytes() == b"existing report"
    assert not (output / "diagnostics" / "ldscore.log").exists()
    assert not (output / "metadata.json").exists()


def test_ldscore_overwrite_keeps_produced_reports_and_reconciles_optional_outputs(tmp_path):
    from ldsc.regression_runner import load_ldscore_from_dir

    argv = plink_arguments(tmp_path)
    source = Path(argv[argv.index("--plink-prefix") + 1])
    prefix = tmp_path / "panel"
    for suffix in (".bed", ".fam"):
        prefix.with_suffix(suffix).write_bytes(source.with_suffix(suffix).read_bytes())
    bim = pd.read_csv(source.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    # Duplicate two monomorphic source identities to produce a real drop report.
    bim.loc[0, "SNP"] = bim.loc[1, "SNP"]
    bim.to_csv(prefix.with_suffix(".bim"), sep="\t", index=False, header=False)
    argv[argv.index("--plink-prefix") + 1] = str(prefix)
    annotation = tmp_path / "baseline.annot"
    bim.assign(base=1.0, feature=[0, 1] * 4).to_csv(annotation, sep="\t", index=False)
    partitioned = [*argv, "--baseline-annot-sources", str(annotation)]

    run_ldscore_from_args(build_parser().parse_args(partitioned))
    output = tmp_path / "result"
    unrelated = output / "notes.txt"
    unrelated.write_text("keep my notes")
    stale_report = output / "diagnostics" / "dropped_snps" / "chr2_dropped.tsv.gz"
    stale_report.write_bytes(b"previous chromosome")
    for arguments, has_overlap in ((partitioned, True), (argv, False)):
        result = run_ldscore_from_args(build_parser().parse_args([*arguments, "--overwrite"]))
        metadata = json.loads((output / "metadata.json").read_text())
        for relative_path in metadata["files"].values():
            assert (output / relative_path).is_file(), relative_path
        assert "dropped_snps_chr1" in metadata["files"]
        assert not stale_report.exists()
        assert (output / "ldscore.overlap.parquet").exists() == has_overlap
        assert (load_ldscore_from_dir(output).overlap is not None) == has_overlap
        assert output / "diagnostics" / "ldscore.log" not in {Path(path) for path in result.output_paths.values()}
        assert unrelated.read_text() == "keep my notes"
