"""Bit-for-bit guardrails for the PLINK genotype reader refactor.

These tests pin the current reader and ``build-ref-panel`` R2 output so the
selective-read (O1), fused individual filter (O3), and streaming (O2) changes
cannot alter any observable value. If a golden assertion fails, the change is
wrong -- fix the change, not the golden.
"""
import glob
import importlib.util
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

_HAS_BA = importlib.util.find_spec("bitarray") is not None
pytestmark = pytest.mark.skipif(not _HAS_BA, reason="bitarray required")

from ldsc._kernel import formats as legacy_parse
from ldsc._kernel.plink_bed import PlinkBEDFile

FIX = Path("tests/fixtures/plink/plink")
READER_GOLDEN = Path("tests/fixtures/golden/reader_golden.npz")
BUILD_GOLDEN = Path("tests/fixtures/golden/build_ref_panel_r2.npz")


def _build(keep_snps=None, keep_indivs=None, maf_min=None):
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(
        str(FIX) + ".bed", len(fam.IDList), bim,
        keep_snps=keep_snps, keep_indivs=keep_indivs, mafMin=maf_min,
    )


def _decode_all(geno):
    geno._currentSNP = 0
    cols = [geno.nextSNPs(1, dtype=np.float64).ravel() for _ in range(geno.m)]
    return np.column_stack(cols) if cols else np.zeros((geno.n, 0))


def test_in_ram_reader_matches_golden():
    g = np.load(READER_GOLDEN)
    geno = _build()
    assert int(geno.m) == int(g["m"])
    assert int(geno.n) == int(g["n"])
    np.testing.assert_array_equal(np.asarray(geno.kept_snps), g["kept_snps"])
    np.testing.assert_allclose(np.asarray(geno.freq), g["freq"], rtol=0, atol=0)
    np.testing.assert_allclose(np.asarray(geno.maf), g["maf"], rtol=0, atol=0)
    np.testing.assert_allclose(_decode_all(geno), g["decoded"], rtol=0, atol=0)


def _run_build_ref_panel(tmp_path):
    """Run an unrestricted source-only build on the real fixture; return parquet path."""
    from ldsc.ref_panel_builder import ReferencePanelBuilder, ReferencePanelBuildConfig
    from ldsc.config import GlobalConfig

    for ext in (".bed", ".bim", ".fam"):
        shutil.copy(str(FIX) + ext, tmp_path / f"panel.1{ext}")
    cfg = ReferencePanelBuildConfig(
        plink_prefix=tmp_path / "panel.@",
        source_genome_build="hg19",
        output_dir=tmp_path / "out",
        ld_wind_snps=4,
    )
    ReferencePanelBuilder(global_config=GlobalConfig(snp_identifier="rsid")).run(cfg)
    paths = sorted(glob.glob(str(tmp_path / "out" / "**" / "*_r2.parquet"), recursive=True))
    assert paths, "no parquet produced"
    return paths[0]


def _read_r2(parquet_path):
    import pyarrow.parquet as pq

    table = pq.read_table(parquet_path)
    return {name: table.column(name).to_numpy() for name in table.column_names}


def test_build_ref_panel_parquet_golden(tmp_path):
    out = _read_r2(_run_build_ref_panel(tmp_path))
    if not BUILD_GOLDEN.exists():
        np.savez(BUILD_GOLDEN, **out)
        pytest.skip("captured build-ref-panel golden; rerun to assert")
    g = np.load(BUILD_GOLDEN)
    assert set(out) == set(g.files)
    for name in out:
        np.testing.assert_array_equal(out[name], g[name], err_msg=f"column {name} drifted")
