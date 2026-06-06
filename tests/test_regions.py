from importlib import resources

import numpy as np
import pytest

from ldsc._kernel import regions
from ldsc.errors import LDSCConfigError, LDSCInputError, LDSCUsageError


def _read_bed(name: str) -> list[tuple[str, int, int]]:
    path = resources.files("ldsc").joinpath("data", "regions", name)
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        chrom, start, end = line.split("\t")[:3]
        rows.append((chrom, int(start), int(end)))
    return rows


@pytest.mark.parametrize("build", ["hg19", "hg38"])
def test_mhc_bed_is_single_chr6_interval(build):
    rows = _read_bed(f"mhc.{build}.bed")
    assert len(rows) == 1
    chrom, start, end = rows[0]
    assert chrom == "6"
    assert 0 < start < end


@pytest.mark.parametrize("build", ["hg19", "hg38"])
def test_centromere_bed_covers_all_autosomes(build):
    rows = _read_bed(f"centromeres.{build}.bed")
    chroms = {chrom for chrom, _, _ in rows}
    assert ({str(i) for i in range(1, 23)} | {"X"}).issubset(chroms)
    for chrom, start, end in rows:
        assert start < end, f"{chrom}: start {start} not < end {end}"


def test_load_preset_intervals_mhc_hg19():
    result = regions.load_preset_intervals(["mhc"], "hg19")
    assert set(result.intervals) == {"6"}
    bounds = result.intervals["6"]
    assert bounds.shape == (1, 2)
    assert bounds[0, 0] == 25_000_000 and bounds[0, 1] == 35_000_000
    assert result.source_labels == ("preset:mhc[hg19]",)


def test_load_preset_intervals_unions_presets():
    result = regions.load_preset_intervals(["mhc", "centromeres"], "hg19")
    assert "6" in result.intervals  # chr6 has both MHC and a centromere
    assert "1" in result.intervals  # centromere only
    assert set(result.source_labels) == {"preset:mhc[hg19]", "preset:centromeres[hg19]"}


def test_load_preset_unknown_name_raises():
    with pytest.raises(LDSCUsageError, match="unknown region preset"):
        regions.load_preset_intervals(["telomeres"], "hg19")


def test_load_preset_unknown_build_raises():
    with pytest.raises(LDSCConfigError, match="genome build"):
        regions.load_preset_intervals(["mhc"], "hg37")


def test_load_bed_intervals_normalizes_chrom_and_merges(tmp_path):
    bed = tmp_path / "user.bed"
    bed.write_text("chr1\t100\t200\nchr1\t150\t300\n2\t10\t20\n", encoding="utf-8")
    result = regions.load_bed_intervals([str(bed)])
    # chr1 -> "1"; overlapping [100,200)+[150,300) coalesce to [100,300)
    assert result.intervals["1"].tolist() == [[100, 300]]
    assert result.intervals["2"].tolist() == [[10, 20]]
    assert result.source_labels == (f"bed:{bed}",)


def test_load_bed_intervals_merges_adjacent(tmp_path):
    bed = tmp_path / "adjacent.bed"
    bed.write_text("1\t100\t200\n1\t200\t300\n", encoding="utf-8")
    result = regions.load_bed_intervals([str(bed)])
    # Adjacent half-open intervals coalesce into one contiguous block.
    assert result.intervals["1"].tolist() == [[100, 300]]


def test_load_bed_intervals_malformed_row_raises(tmp_path):
    bed = tmp_path / "bad.bed"
    bed.write_text("1\t100\n", encoding="utf-8")
    with pytest.raises(LDSCInputError, match="bad.bed"):
        regions.load_bed_intervals([str(bed)])


def test_merge_intervals_unions_sources():
    a = regions.load_preset_intervals(["mhc"], "hg19")
    b = regions.RegionIntervals(intervals={"6": np.array([[40_000_000, 41_000_000]], dtype=np.int64)}, source_labels=("bed:x",))
    merged = regions.merge_intervals(a, b)
    # chr6 now has the MHC interval plus the extra disjoint interval, sorted
    assert merged.intervals["6"].tolist() == [[25_000_000, 35_000_000], [40_000_000, 41_000_000]]
    assert set(merged.source_labels) == {"preset:mhc[hg19]", "bed:x"}
