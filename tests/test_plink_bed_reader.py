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


def test_keep_indivs_matches_golden():
    g = np.load(READER_GOLDEN)
    geno = _build(keep_indivs=[0, 1, 2, 3])
    assert int(geno.m) == int(g["ki_m"])
    assert int(geno.n) == int(g["ki_n"])
    np.testing.assert_array_equal(np.asarray(geno.kept_snps), g["ki_kept_snps"])
    np.testing.assert_allclose(np.asarray(geno.freq), g["ki_freq"], rtol=0, atol=0)
    np.testing.assert_allclose(_decode_all(geno), g["ki_decoded"], rtol=0, atol=0)


def test_selective_read_subset_matches_full():
    full = _build()
    full_decoded = _decode_all(full)
    keep = list(full.kept_snps[:3])  # guaranteed-polymorphic source indices
    sub = _build(keep_snps=keep)
    assert list(sub.kept_snps) == keep
    assert int(sub.m) == len(keep)
    np.testing.assert_allclose(np.asarray(sub.freq), np.asarray(full.freq[:3]),
                               rtol=0, atol=0)
    np.testing.assert_allclose(_decode_all(sub), full_decoded[:, :3], rtol=0, atol=0)


def test_no_whole_chromosome_bitarray_after_init():
    full = _build()
    keep = list(full.kept_snps[:2])
    geno = _build(keep_snps=keep)
    # retained bitarray holds only kept SNPs (2 * nru * m bits), not the chromosome
    assert len(geno.geno) == 2 * geno.nru * geno.m


def _build_streaming(keep_indivs=None, maf_min=None):
    bim = legacy_parse.PlinkBIMFile(str(FIX) + ".bim")
    fam = legacy_parse.PlinkFAMFile(str(FIX) + ".fam")
    return PlinkBEDFile(
        str(FIX) + ".bed", len(fam.IDList), bim,
        keep_snps=None, keep_indivs=keep_indivs, mafMin=maf_min, streaming=True,
    )


def test_streaming_matches_in_ram_unrestricted():
    in_ram = _build()
    streamed = _build_streaming()
    assert int(streamed.m) == int(in_ram.m)
    assert int(streamed.n) == int(in_ram.n)
    assert streamed.geno is None  # streaming never materializes the bitarray
    np.testing.assert_array_equal(np.asarray(streamed.kept_snps),
                                  np.asarray(in_ram.kept_snps))
    np.testing.assert_allclose(np.asarray(streamed.freq),
                               np.asarray(in_ram.freq), rtol=0, atol=0)
    np.testing.assert_allclose(_decode_all(streamed), _decode_all(in_ram),
                               rtol=0, atol=0)


def test_streaming_matches_in_ram_with_keep_indivs():
    in_ram = _build(keep_indivs=[0, 1, 2, 3])
    streamed = _build_streaming(keep_indivs=[0, 1, 2, 3])
    np.testing.assert_array_equal(np.asarray(streamed.kept_snps),
                                  np.asarray(in_ram.kept_snps))
    np.testing.assert_allclose(_decode_all(streamed), _decode_all(in_ram),
                               rtol=0, atol=0)


def test_streaming_rewind_replays_identically():
    streamed = _build_streaming()
    first = _decode_all(streamed)   # _decode_all resets _currentSNP to 0
    second = _decode_all(streamed)  # rewind + replay
    np.testing.assert_allclose(first, second, rtol=0, atol=0)


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
