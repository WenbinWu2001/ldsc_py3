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
    assert "is wider than the input R2 panel window" in issues["details"].iloc[0]
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


def _tied_snp_panel(build_dir):
    """A 4-SNP rsID panel whose rows 1 and 2 (rs3, then rs2) share position 200."""
    from ldsc._refpanel import build_kernel as kb
    from ldsc._core.snp_identity import sidecar_identity_sha256
    from ldsc._refpanel.r2_format import (
        build_runtime_metadata_table,
        write_r2_parquet,
        write_runtime_metadata_sidecar,
    )

    panel = pd.DataFrame({"CHR": ["1"] * 4, "POS": [100, 200, 200, 300], "SNP": ["rs1", "rs3", "rs2", "rs4"],
                          "A1": ["A", "A", "A", "A"], "A2": ["C", "C", "G", "C"], "MAF": [0.3, 0.25, 0.4, 0.35]})
    build_dir.mkdir(parents=True)
    sidecar = build_runtime_metadata_table(metadata=panel, positions=panel["POS"].to_numpy(), cm_values=None)
    pairs = [(0, 1, 0.5), (0, 2, 0.1), (1, 2, 0.4), (1, 3, 0.2), (2, 3, 0.6)]
    write_r2_parquet(
        pair_chunks=[(np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]),
                      np.array([p[2] for p in pairs], dtype=np.float32), np.ones(len(pairs), dtype=np.int8))],
        path=build_dir / "chr1_r2.parquet", genome_build="hg19", n_samples=200, snp_identifier="rsid",
        n_snps=len(panel), sidecar_identity_sha256=sidecar_identity_sha256(sidecar),
        ld_window_mode="snps", ld_window_value=1.0,
    )
    write_runtime_metadata_sidecar(sidecar, build_dir / "chr1_meta.tsv.gz", genome_build="hg19",
                                      snp_identifier="rsid")
    return sidecar


def test_snp_count_window_computes_same_position_snps_in_panel_order(tmp_path):
    """An annotation may list same-position SNPs in another order than the panel sidecar.

    ``--ld-wind-snps`` counts panel rows, so the SNPs are computed in panel-row
    order. The written LD scores, and their row order (CHR, POS, then SNP), equal
    those of the same run with the annotation listed in sidecar order.
    """
    from ldsc import AnnotationBundle, GlobalConfig, LDScoreCalculator, LDScoreConfig, RefPanelConfig
    from ldsc._refpanel.adapters import ParquetR2RefPanel
    from ldsc.outputs import LDScoreOutputConfig

    sidecar = _tied_snp_panel(tmp_path / "panel")
    config = GlobalConfig(snp_identifier="rsid")
    panel = ParquetR2RefPanel(config, RefPanelConfig(backend="parquet_r2", r2_dir=str(tmp_path / "panel")))
    query = pd.DataFrame({"q": [0.1, 0.2, 0.3, 0.4]})

    def written_scores(order, name):
        metadata = sidecar[["CHR", "POS", "SNP", "A1", "A2"]].iloc[order].reset_index(drop=True)
        with AnnotationBundle.from_frames(metadata, pd.DataFrame({"base": np.ones(4)}),
                                          query.iloc[order].reset_index(drop=True), config_snapshot=config) as bundle:
            LDScoreCalculator().run(bundle, panel, LDScoreConfig(ld_wind_snps=1), config,
                                    LDScoreOutputConfig(tmp_path / name))
        return [pd.read_parquet(tmp_path / name / f"ldscore.{kind}.parquet") for kind in ("baseline", "query")]

    in_panel_order = written_scores([0, 1, 2, 3], "panel-order")
    swapped = written_scores([0, 2, 1, 3], "swapped")

    for expected, observed in zip(in_panel_order, swapped):
        pd.testing.assert_frame_equal(observed, expected)
    assert swapped[0]["SNP"].tolist() == ["rs1", "rs2", "rs3", "rs4"]


WARNING_TEXT = "counts positions in the R2 panel"


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


def test_cm_panel_stores_every_pair_inside_the_window_its_sidecar_records(tmp_path):
    """Build-time windowing uses the CM values the sidecar stores, so the two agree exactly.

    ``ldscore`` rebuilds a parquet cM window from the sidecar CM. The map places two
    SNP pairs on the 1 cM edge: 20.0000004 to 21.0000008 is outside it, but reads
    as 20 to 21 at 6 significant digits; 30.0000000004 to 31.0000000009 is outside
    it at full precision but reads as 30 to 31 at 9 digits, so it is inside the
    window the sidecar records.
    """
    bim = pd.read_csv(f"{PLINK}.bim", sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    cm = 5.0 + 0.7 * np.arange(len(bim))
    cm[10:20] = 20.0000004 + 0.7 * np.arange(10)
    cm[11] = 21.0000008
    cm[20:] = 30.0000000004 + 0.7 * np.arange(len(bim) - 20)
    cm[21] = 31.0000000009
    genetic_map = tmp_path / "hg38.map"
    pd.DataFrame({"CHR": bim["CHR"], "POS": bim["POS"], "CM": cm}).to_csv(
        genetic_map, sep="\t", index=False, float_format="%.13f")
    out = tmp_path / "panel"
    assert cli.run_cli(["build-r2-panel", "--plink-prefix", str(PLINK), "--output-dir", str(out), *IDENTITY,
                        "--source-genome-build", "hg38", "--genetic-map-hg38-sources", str(genetic_map),
                        "--ld-wind-cm", "1"]) == 0

    stored_cm = pd.read_csv(out / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#")["CM"].to_numpy()
    rows = np.arange(len(stored_cm))
    upper = rows[:, None] < rows[None, :]
    within = upper & ~(np.abs(stored_cm[:, None] - stored_cm[None, :]) > 1.0)
    expected = set(zip(*np.nonzero(within)))
    pairs = pd.read_parquet(out / "hg38" / "chr22_r2.parquet", columns=["IDX_1", "IDX_2"])
    stored = {(min(i, j), max(i, j)) for i, j in zip(pairs["IDX_1"], pairs["IDX_2"]) if i != j}
    assert stored == expected
    assert (10, 11) not in expected and (20, 21) in expected


def _genotype_maf(prefix):
    """Minor-allele frequency of each SNP, from decoded dosages (no missing calls)."""
    freq = _dosages(prefix).mean(axis=1) / 2
    return np.minimum(freq, 1 - freq)


def test_both_sidecar_writers_keep_nine_significant_digits_of_maf(tmp_path, panels):
    """MAF = k / 6404 here, so 6 significant digits would be off by up to 5e-6 relative."""
    maf = pd.Series(_genotype_maf(PLINK), index=pd.read_csv(
        f"{PLINK}.bim", sep=r"\s+", header=None, usecols=[3]).iloc[:, 0].to_numpy())
    output = tmp_path / "out"
    assert _ldscore(output, "--plink-prefix", str(PLINK), "--ld-wind-snps", "3", "--export-ref-metadata") == 0
    for path in (panels["snps3"] / "hg38" / "chr22_meta.tsv.gz", output / "ref_metadata" / "chr22_meta.tsv.gz"):
        sidecar = pd.read_csv(path, sep="\t", comment="#")
        np.testing.assert_allclose(sidecar["MAF"].to_numpy(), maf.loc[sidecar["POS"]].to_numpy(), rtol=5e-9, atol=0,
                                   err_msg=str(path))


def test_build_maf_filter_holds_for_the_stored_maf(tmp_path):
    """Rounding is monotone, so every stored MAF still meets the build --maf-min."""
    out = tmp_path / "panel"
    assert cli.run_cli(["build-r2-panel", "--plink-prefix", str(PLINK), "--output-dir", str(out), *IDENTITY,
                        "--source-genome-build", "hg38", "--ld-wind-snps", "3", "--maf-min", "0.3"]) == 0

    stored = pd.read_csv(out / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#")["MAF"]
    assert 0 < len(stored) < 32 and (stored >= 0.3).all()
