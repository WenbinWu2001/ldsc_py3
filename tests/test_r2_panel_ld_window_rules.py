"""LD-window rules for ldscore runs on parquet R2 panels.

Fixture: the deterministic 32-SNP chromosome-22 1KG-30x subset (3202 samples,
hg38 positions; see ``tests/fixtures/minimal_external_resources/README.md``).
A panel stores only the pairs inside its build window, so a run may use only
windows the panel is known to cover from its recorded settings.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import cli

pytest.importorskip("bitarray")
pytest.importorskip("pyarrow")

PLINK = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
IDENTITY = ["--snp-identifier", "chr_pos"]


@pytest.fixture(scope="module")
def panels(tmp_path_factory):
    root = tmp_path_factory.mktemp("panels")
    built = {}
    windows = {"kb200": ["--ld-wind-kb", "200"], "snps3": ["--ld-wind-snps", "3"], "snps6": ["--ld-wind-snps", "6"]}
    for name, window in windows.items():
        built[name] = root / name
        assert cli.run_cli(["build-r2-panel", "--plink-prefix", str(PLINK), "--output-dir", str(built[name]),
                            *IDENTITY, "--source-genome-build", "hg38", *window]) == 0
    return built


def _ldscore(output_dir, *args):
    return cli.run_cli(["ldscore", "--output-dir", str(output_dir), *IDENTITY, "--genome-build", "hg38", *args])


def test_window_in_a_different_unit_is_rejected_before_any_chromosome_is_computed(tmp_path, panels):
    output = tmp_path / "out"
    assert _ldscore(output, "--r2-dir", str(panels["kb200"]), "--ld-wind-snps", "3") != 0

    issues = pd.read_csv(output / "diagnostics" / "input_issues.tsv", sep="\t")
    assert issues["chrom"].astype(str).tolist() == ["22"]
    assert "windows in different units cannot be checked for coverage" in issues["details"].iloc[0]
    assert "the panel was built with `--ld-wind-kb 200.0`" in issues["details"].iloc[0]
    assert not (output / "diagnostics" / "chromosome-logs").exists()
    assert not (output / "ldscore.baseline.parquet").exists()


def test_wider_window_in_the_same_unit_is_rejected_before_any_chromosome_is_computed(tmp_path, panels):
    output = tmp_path / "out"
    assert _ldscore(output, "--r2-dir", str(panels["snps3"]), "--ld-wind-snps", "4") != 0

    issues = pd.read_csv(output / "diagnostics" / "input_issues.tsv", sep="\t")
    assert "is wider than the input R2 parquet panel window" in issues["details"].iloc[0]
    assert not (output / "diagnostics" / "chromosome-logs").exists()


def _dosages(prefix):
    """Decode SNP-major PLINK genotypes as A1 counts, one row per SNP."""
    n = len(Path(f"{prefix}.fam").read_text().splitlines())
    m = len(Path(f"{prefix}.bim").read_text().splitlines())
    raw = np.fromfile(f"{prefix}.bed", dtype=np.uint8)
    assert raw[:3].tolist() == [0x6C, 0x1B, 0x01]
    bits = np.unpackbits(raw[3:].reshape(m, -1), axis=1, bitorder="little").reshape(m, -1, 2)[:, :n]
    code = bits[..., 0] + 2 * bits[..., 1]  # 0 hom A1, 1 missing, 2 het, 3 hom A2
    assert not (code == 1).any(), "the fixture has no missing calls"
    return np.choose(code, [2, -1, 1, 0]).astype(float)


@pytest.mark.parametrize(("panel", "window"), [("snps3", 3), ("snps6", 2)])
def test_snp_count_window_counts_panel_rows_when_the_run_keeps_fewer_snps(tmp_path, panels, panel, window):
    """``--ld-wind-snps K`` spans K panel rows per side, even after ``--maf-min``.

    Counting only the kept SNPs would reach past K panel rows: beyond the build
    window for an equal K (pairs the panel never stored), or further than the
    request for a narrower K. The oracle sums adjusted r2 (Pearson r from
    decoded dosages, n = 3202) over kept partners within K panel rows.
    """
    output = tmp_path / "out"
    assert _ldscore(output, "--r2-dir", str(panels[panel]), "--ld-wind-snps", str(window), "--maf-min", "0.3") == 0

    sidecar = pd.read_csv(panels[panel] / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#")
    bim = pd.read_csv(f"{PLINK}.bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    dosage = _dosages(PLINK)[bim.reset_index().set_index("POS").loc[sidecar["POS"], "index"].to_numpy()]
    n = dosage.shape[1]
    r2 = np.corrcoef(dosage) ** 2
    adjusted = r2 - (1 - r2) / (n - 2)
    rows = np.arange(len(sidecar))
    kept = sidecar["MAF"].to_numpy() >= 0.3
    assert 0 < kept.sum() < len(sidecar)
    within = np.abs(rows[:, None] - rows[None, :]) <= window
    expected = pd.Series((np.where(within, adjusted, 0.0)[:, kept]).sum(axis=1)[kept],
                         index=sidecar.loc[kept, "POS"].to_numpy())

    scores = pd.read_parquet(output / "ldscore.baseline.parquet").set_index("POS")["base"]
    assert len(scores) > 0 and scores.index.isin(expected.index).all()
    # int16 R2 storage rounds each stored pair by at most 0.5/32767; a SNP has at
    # most 2 * window <= 6 off-diagonal partners here, so scores agree within 6 * 1.6e-5.
    np.testing.assert_allclose(scores.to_numpy(), expected.loc[scores.index].to_numpy(), rtol=0, atol=1e-4)


def test_snp_count_window_rejects_kept_snps_out_of_panel_row_order():
    """Two SNPs at one position listed in the opposite order to the panel sidecar."""
    from argparse import Namespace

    from ldsc._kernel.ldscore import build_window_coordinates
    from ldsc.errors import LDSCInputError

    metadata = pd.DataFrame({"CHR": ["22"] * 3, "SNP": ["rs1", "rs3", "rs2"], "POS": [10, 20, 20]})
    args = Namespace(ld_wind_snps=1, ld_wind_kb=None, ld_wind_cm=None)
    with pytest.raises(LDSCInputError, match="first out of order: rs2"):
        build_window_coordinates(metadata, args, snp_positions=np.array([0, 2, 1]))
    coords, max_dist = build_window_coordinates(metadata, args, snp_positions=np.array([0, 1, 4]))
    np.testing.assert_array_equal(coords, [0, 1, 4])
    assert max_dist == 1


WARNING_TEXT = "counts positions in the parquet R2 panel"


def test_snp_count_window_records_panel_and_kept_counts_and_warns_once(tmp_path, panels):
    output = tmp_path / "out"
    assert _ldscore(output, "--r2-dir", str(panels["snps6"]), "--ld-wind-snps", "2", "--maf-min", "0.3") == 0

    sidecar = pd.read_csv(panels["snps6"] / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#")
    kept = int((sidecar["MAF"] >= 0.3).sum())
    counts = pd.read_csv(output / "diagnostics" / "ld_window_snp_counts.tsv", sep="\t", dtype={"CHR": str})
    assert counts.to_dict("records") == [{"CHR": "22", "panel_snps": len(sidecar), "kept_snps": kept,
                                         "ld_wind_snps": 2, "counting_basis": "panel_rows"}]
    files = json.loads((output / "metadata.json").read_text())["files"]
    assert files["ld_window_snp_counts"] == "diagnostics/ld_window_snp_counts.tsv"
    warnings = [line for line in (output / "diagnostics" / "ldscore.log").read_text().splitlines()
                if "[WARNING]" in line and WARNING_TEXT in line]
    assert len(warnings) == 1
    assert f"chr22: {kept} of {len(sidecar)}" in warnings[0]
    assert "diagnostics/ld_window_snp_counts.tsv" in warnings[0]


@pytest.mark.parametrize(("source", "window", "table"), [
    ("snps3", ["--ld-wind-snps", "3"], True),
    ("kb200", ["--ld-wind-kb", "200"], False),
    ("plink", ["--ld-wind-snps", "3"], False),
])
def test_no_ld_window_warning_without_a_snp_subset(tmp_path, panels, source, window, table):
    output = tmp_path / "out"
    reference = ["--plink-prefix", str(PLINK)] if source == "plink" else ["--r2-dir", str(panels[source])]
    assert _ldscore(output, *reference, *window) == 0

    assert WARNING_TEXT not in (output / "diagnostics" / "ldscore.log").read_text()
    path = output / "diagnostics" / "ld_window_snp_counts.tsv"
    assert path.exists() == table
    if table:
        counts = pd.read_csv(path, sep="\t")
        assert (counts["kept_snps"] == counts["panel_snps"]).all()
