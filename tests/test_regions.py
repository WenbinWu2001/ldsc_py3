from importlib import resources

import pytest


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
