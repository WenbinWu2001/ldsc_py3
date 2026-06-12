import gzip
from importlib import resources

import numpy as np
import pandas as pd
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


def test_load_bed_intervals_skips_non_interval_lines_and_standard_header(tmp_path):
    bed = tmp_path / "user.bed"
    bed.write_text(
        "\n# comment line\ntrack name=\"regions\"\nbrowser position chr1:1-1000\n"
        "chrom\tstart\tend\tname\nchr1\t100\t200\tregion_a\n",
        encoding="utf-8",
    )

    result = regions.load_bed_intervals([str(bed)])

    assert result.intervals["1"].tolist() == [[100, 200]]


def test_load_bed_intervals_skips_chrom_start_end_header_variant(tmp_path):
    bed = tmp_path / "user.bed"
    bed.write_text("chrom\tchromStart\tchromEnd\nchr2\t10\t20\n", encoding="utf-8")

    result = regions.load_bed_intervals([str(bed)])

    assert result.intervals["2"].tolist() == [[10, 20]]


def test_load_bed_intervals_accepts_gzip_bed(tmp_path):
    bed = tmp_path / "user.bed.gz"
    with gzip.open(bed, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tstart\tend\nchr1\t100\t200\n")

    result = regions.load_bed_intervals([str(bed)])

    assert result.intervals["1"].tolist() == [[100, 200]]


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


def test_load_bed_intervals_rejects_header_after_data(tmp_path):
    bed = tmp_path / "bad_header.bed"
    bed.write_text("chr1\t100\t200\nchrom\tstart\tend\n", encoding="utf-8")

    with pytest.raises(LDSCInputError, match="header"):
        regions.load_bed_intervals([str(bed)])


def test_load_bed_intervals_rejects_multiple_leading_headers(tmp_path):
    bed = tmp_path / "bad_headers.bed"
    bed.write_text("chrom\tstart\tend\nchrom\tchromStart\tchromEnd\nchr1\t100\t200\n", encoding="utf-8")

    with pytest.raises(LDSCInputError, match="header"):
        regions.load_bed_intervals([str(bed)])


def test_load_bed_intervals_rejects_invalid_standard_bed_bounds(tmp_path):
    cases = {
        "negative_start.bed": "chr1\t-1\t200\n",
        "zero_length.bed": "chr1\t200\t200\n",
        "reversed.bed": "chr1\t200\t100\n",
    }
    for filename, content in cases.items():
        bed = tmp_path / filename
        bed.write_text(content, encoding="utf-8")

        with pytest.raises(LDSCInputError, match="start"):
            regions.load_bed_intervals([str(bed)])


def test_merge_intervals_unions_sources():
    a = regions.load_preset_intervals(["mhc"], "hg19")
    b = regions.RegionIntervals(intervals={"6": np.array([[40_000_000, 41_000_000]], dtype=np.int64)}, source_labels=("bed:x",))
    merged = regions.merge_intervals(a, b)
    # chr6 now has the MHC interval plus the extra disjoint interval, sorted
    assert merged.intervals["6"].tolist() == [[25_000_000, 35_000_000], [40_000_000, 41_000_000]]
    assert set(merged.source_labels) == {"preset:mhc[hg19]", "bed:x"}


def _meta(chrom_pos: list[tuple[str, int]]) -> pd.DataFrame:
    return pd.DataFrame({"CHR": [c for c, _ in chrom_pos], "POS": [p for _, p in chrom_pos]})


def test_keep_mask_half_open_boundaries():
    # BED [100, 200) 0-based excludes 1-based POS 101..200 inclusive.
    intervals = regions.RegionIntervals(
        intervals={"1": np.array([[100, 200]], dtype=np.int64)}, source_labels=("bed:x",)
    )
    meta = _meta([("1", 100), ("1", 101), ("1", 200), ("1", 201)])
    keep = regions.region_exclusion_keep_mask(meta, intervals)
    assert keep.tolist() == [True, False, False, True]


def test_keep_mask_other_chromosome_untouched():
    intervals = regions.RegionIntervals(
        intervals={"1": np.array([[100, 200]], dtype=np.int64)}, source_labels=("bed:x",)
    )
    meta = _meta([("2", 150), ("1", 150)])
    keep = regions.region_exclusion_keep_mask(meta, intervals)
    assert keep.tolist() == [True, False]


def test_keep_mask_empty_intervals_keeps_all():
    intervals = regions.RegionIntervals(intervals={}, source_labels=())
    meta = _meta([("1", 10), ("2", 20)])
    keep = regions.region_exclusion_keep_mask(meta, intervals)
    assert keep.tolist() == [True, True]


def test_keep_mask_multiple_disjoint_intervals():
    intervals = regions.RegionIntervals(
        intervals={"1": np.array([[100, 200], [300, 400]], dtype=np.int64)}, source_labels=("bed:x",)
    )
    meta = _meta([("1", 150), ("1", 250), ("1", 350)])
    keep = regions.region_exclusion_keep_mask(meta, intervals)
    assert keep.tolist() == [False, True, False]


def test_keep_mask_uses_bp_column_when_requested():
    intervals = regions.RegionIntervals(
        intervals={"1": np.array([[100, 200]], dtype=np.int64)}, source_labels=("bed:x",)
    )
    meta = pd.DataFrame({"CHR": ["1"], "BP": [150]})
    keep = regions.region_exclusion_keep_mask(meta, intervals, pos_col="BP")
    assert keep.tolist() == [False]


def test_keep_mask_empty_metadata_with_intervals():
    intervals = regions.RegionIntervals(
        intervals={"1": np.array([[100, 200]], dtype=np.int64)}, source_labels=("bed:x",)
    )
    keep = regions.region_exclusion_keep_mask(_meta([]), intervals)
    assert keep.tolist() == []


def test_keep_mask_interval_chrom_absent_from_metadata():
    intervals = regions.RegionIntervals(
        intervals={"9": np.array([[100, 200]], dtype=np.int64)}, source_labels=("bed:x",)
    )
    keep = regions.region_exclusion_keep_mask(_meta([("1", 150), ("2", 150)]), intervals)
    assert keep.tolist() == [True, True]


def test_exclude_regions_choice_to_presets_maps_all_choices():
    assert regions.exclude_regions_choice_to_presets("none") == ()
    assert regions.exclude_regions_choice_to_presets("mhc") == ("mhc",)
    assert regions.exclude_regions_choice_to_presets("centromeres") == ("centromeres",)
    assert regions.exclude_regions_choice_to_presets("mhc-and-centromeres") == ("mhc", "centromeres")
    # Mapped presets are valid REGION_PRESETS members.
    assert set(regions.exclude_regions_choice_to_presets("mhc-and-centromeres")) <= regions.REGION_PRESETS


def test_exclude_regions_choice_rejects_unknown_token():
    with pytest.raises(LDSCConfigError, match="Invalid --exclude-regions choice"):
        regions.exclude_regions_choice_to_presets("mhc,centromeres")


@pytest.mark.parametrize("build", ["hg19", "hg38"])
def test_centromeres_core_bed_covers_all_autosomes(build):
    rows = _read_bed(f"centromeres_core.{build}.bed")
    chroms = {chrom for chrom, _, _ in rows}
    assert {str(i) for i in range(1, 23)}.issubset(chroms)
    for chrom, start, end in rows:
        assert start < end, f"{chrom}: start {start} not < end {end}"


def test_centromeres_core_in_preset_menu():
    assert "centromeres_core" in regions.REGION_PRESETS
    result = regions.load_preset_intervals(["centromeres_core"], "hg38")
    assert result.source_labels == ("preset:centromeres_core[hg38]",)
    assert "1" in result.intervals


@pytest.mark.parametrize("build", ["hg19", "hg38"])
def test_centromeres_preset_is_pericentromeric_superset_of_core(build):
    # The active `centromeres` preset ships the pericentromeric (+/-3 cM) region,
    # which must always contain the raw `centromeres_core` span -- including on
    # acrocentric chromosomes (13/14/15/21/22) whose p-arm is absent from the map.
    core = regions.load_preset_intervals(["centromeres_core"], build)
    active = regions.load_preset_intervals(["centromeres"], build)
    for chrom, core_bounds in core.intervals.items():
        assert chrom in active.intervals, chrom
        active_bounds = active.intervals[chrom].tolist()
        for cstart, cend in core_bounds.tolist():
            covered = any(astart <= cstart and aend >= cend for astart, aend in active_bounds)
            assert covered, f"chr{chrom} core centromere {cstart}-{cend} not within the padded 'centromeres' region"
    assert active.source_labels == (f"preset:centromeres[{build}]",)
