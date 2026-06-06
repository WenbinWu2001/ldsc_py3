"""One-time curation tool: regenerate packaged region BED files from UCSC.

Run manually from the repo root:
    python tools/regions/build_region_beds.py

Outputs standard 0-based half-open BED (chrom start end) to
src/ldsc/data/regions/. Not shipped with the package.

Sources:
  - Centromeres: UCSC REST API. hg19 uses the `gap` track (type == centromere);
    hg38 uses the `centromeres` track (multiple models per chromosome, merged to
    one [min(start), max(end)) span per chromosome).
  - MHC: pinned broad-exclusion windows (hg19 chr6:25-35Mb; hg38 lifted core).
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parents[2] / "src" / "ldsc" / "data" / "regions"
AUTOSOMES_PLUS_X = [str(i) for i in range(1, 23)] + ["X"]

MHC = {
    "hg19": [("6", 25_000_000, 35_000_000)],
    "hg38": [("6", 28_477_797, 33_448_354)],
}


def _fetch(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return json.load(resp)


def _norm_chrom(token: str) -> str | None:
    token = token[3:] if token.startswith("chr") else token
    return token if token in AUTOSOMES_PLUS_X else None


def centromeres_hg19() -> list[tuple[str, int, int]]:
    data = _fetch("https://api.genome.ucsc.edu/getData/track?genome=hg19;track=gap")
    rows: list[tuple[str, int, int]] = []
    for chrom_key, items in data["gap"].items():
        chrom = _norm_chrom(chrom_key)
        if chrom is None:
            continue
        for item in items:
            if item.get("type") == "centromere":
                rows.append((chrom, int(item["chromStart"]), int(item["chromEnd"])))
    return rows


def centromeres_hg38() -> list[tuple[str, int, int]]:
    data = _fetch("https://api.genome.ucsc.edu/getData/track?genome=hg38;track=centromeres")
    spans: dict[str, tuple[int, int]] = {}
    for chrom_key, items in data["centromeres"].items():
        chrom = _norm_chrom(chrom_key)
        if chrom is None:
            continue
        starts = [int(it["chromStart"]) for it in items]
        ends = [int(it["chromEnd"]) for it in items]
        spans[chrom] = (min(starts), max(ends))
    return [(chrom, *spans[chrom]) for chrom in AUTOSOMES_PLUS_X if chrom in spans]


def _chrom_key(chrom: str) -> tuple[int, str]:
    return (int(chrom), "") if chrom.isdigit() else (99, chrom)


def write_bed(path: Path, rows: list[tuple[str, int, int]], provenance: str) -> None:
    rows = sorted(rows, key=lambda r: (_chrom_key(r[0]), r[1]))
    lines = [f"# {provenance}"]
    lines += [f"{chrom}\t{start}\t{end}" for chrom, start, end in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {len(rows)} rows -> {path}")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_bed(OUT_DIR / "mhc.hg19.bed", MHC["hg19"], "MHC broad exclusion window, hg19 (chr6:25-35Mb), pinned constant; 0-based half-open")
    write_bed(OUT_DIR / "mhc.hg38.bed", MHC["hg38"], "MHC core exclusion window, hg38 (chr6:28477797-33448354), pinned constant; 0-based half-open")
    write_bed(OUT_DIR / "centromeres.hg19.bed", centromeres_hg19(), "UCSC hg19 gap track, type=centromere (0-based half-open)")
    write_bed(OUT_DIR / "centromeres.hg38.bed", centromeres_hg38(), "UCSC hg38 centromeres track, merged per-chromosome span (0-based half-open)")


if __name__ == "__main__":
    main()
