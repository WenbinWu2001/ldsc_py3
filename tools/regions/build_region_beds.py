"""One-time curation tool: regenerate packaged region BED files.

Run manually from the repo root:
    python tools/regions/build_region_beds.py

Outputs standard 0-based half-open BED (chrom start end) to
src/ldsc/data/regions/. Not shipped with the package.

Sources:
  - `centromeres_core.{build}.bed`: raw centromere span from UCSC REST API. hg19
    uses the `gap` track (type == centromere); hg38 uses the `centromeres` track
    (multiple models per chromosome, merged to one [min(start), max(end)) span).
    Kept for reference; NOT wired to any CLI choice.
  - `centromeres.{build}.bed` (the active `centromeres` preset): each core span
    padded by +/-3 cM, matching the LD Score regression pericentromeric exclusion
    of Bulik-Sullivan et al. 2015 Nat Genet (Online Methods). cM is read from the
    Alkes-group recombination map and inverted back to base pairs; the padded
    interval is unioned with the core span so it always contains it (see
    `pericentromere_intervals`).
  - MHC: pinned broad-exclusion windows (hg19 chr6:25-35Mb; hg38 lifted core).

The genetic maps are a CURATION-TIME input only (not bundled in the package);
they live in the workspace `resources/` tree by default, overridable with the
`LDSC_GENETIC_MAP_DIR` environment variable.
"""
from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

OUT_DIR = Path(__file__).resolve().parents[2] / "src" / "ldsc" / "data" / "regions"
AUTOSOMES_PLUS_X = [str(i) for i in range(1, 23)] + ["X"]

# Genetic maps used for pericentromere padding. Curation-time input only; the
# default points at the Alkes-group maps in the workspace resources tree.
GENETIC_MAP_DIR = Path(
    os.environ.get(
        "LDSC_GENETIC_MAP_DIR",
        str(Path(__file__).resolve().parents[3] / "resources" / "genetic_maps" / "genetic_map_alkesgroup"),
    )
)
_GENETIC_MAP_FILE = {"hg19": "genetic_map_hg19_withX.txt", "hg38": "genetic_map_hg38_withX.txt"}
PERICENTROMERE_PAD_CM = 3.0

# The active ``centromeres`` preset ships the pericentromeric (+/-3 cM) regions,
# matching the LD Score regression exclusion. The raw centromere gap is preserved
# as ``centromeres_core`` for reference but is not wired to any CLI choice.
_CENTROMERE_PROVENANCE = {
    ("core", "hg19"): "UCSC hg19 gap track, type=centromere; RAW core gap, kept for reference, NOT the active 'centromeres' exclusion; 0-based half-open",
    ("core", "hg38"): "UCSC hg38 centromeres track, merged per-chromosome span; RAW core, kept for reference, NOT the active 'centromeres' exclusion; 0-based half-open",
    ("padded", "hg19"): "Pericentromeric +/-3 cM (LDSC, Bulik-Sullivan 2015 Nat Genet): hg19 centromere padded via Alkes-group genetic_map_hg19_withX; active 'centromeres' exclusion; 0-based half-open",
    ("padded", "hg38"): "Pericentromeric +/-3 cM (LDSC, Bulik-Sullivan 2015 Nat Genet): hg38 centromere padded via Alkes-group genetic_map_hg38_withX; active 'centromeres' exclusion; 0-based half-open",
}

# Active `mhc` preset: broad GWAS exclusion window, build-consistent (chr6:25-35 Mb
# in both builds; hg19 25-35 Mb lifts to hg38 ~25.00-35.03 Mb, so the same numeric
# window covers the same biology). `mhc_core` is the narrow classical-HLA core,
# kept for reference only and not wired to any CLI choice.
MHC_BROAD = {
    "hg19": [("6", 25_000_000, 35_000_000)],
    "hg38": [("6", 25_000_000, 35_000_000)],
}
MHC_CORE = {
    "hg38": [("6", 28_477_797, 33_448_354)],   # pinned hg38 classical-HLA core
    "hg19": [("6", 28_445_574, 33_416_131)],   # liftOver of the hg38 core
}
_MHC_PROVENANCE = {
    ("broad", "hg19"): "MHC broad GWAS exclusion window, hg19 chr6:25-35Mb, pinned; active 'mhc' exclusion; 0-based half-open",
    ("broad", "hg38"): "MHC broad GWAS exclusion window, hg38 chr6:25-35Mb, pinned (hg19 25-35Mb lifts here); active 'mhc' exclusion; 0-based half-open",
    ("core", "hg38"): "MHC core (classical HLA), hg38 chr6:28477797-33448354, pinned; reference only, NOT the active 'mhc' exclusion; 0-based half-open",
    ("core", "hg19"): "MHC core (classical HLA), hg19 chr6:28445574-33416131, liftOver of the hg38 core; reference only, NOT the active 'mhc' exclusion; 0-based half-open",
}


def _fetch(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as resp:
        return json.load(resp)


def _norm_chrom(token: str) -> str | None:
    token = token[3:] if token.startswith("chr") else token
    return token if token in AUTOSOMES_PLUS_X else None


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    """Read a standard BED file into ``(chrom, start, end)`` rows."""
    rows: list[tuple[str, int, int]] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        chrom, start, end = line.split("\t")[:3]
        rows.append((chrom, int(start), int(end)))
    return rows


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


def _load_genetic_map(build: str):
    """Return ``{chrom: (positions, cM)}`` from the Alkes-group genetic map.

    The map encodes chrX as ``23``; this normalizes it to ``X`` to match BED
    chromosome tokens. Per chromosome, ``positions`` and ``cM`` are
    sorted-ascending float arrays.
    """
    import numpy as np
    import pandas as pd

    path = GENETIC_MAP_DIR / _GENETIC_MAP_FILE[build]
    frame = pd.read_csv(path, sep=r"\s+")
    frame.columns = ["chr", "pos", "rate", "cM"]
    frame["chr"] = frame["chr"].astype(str).replace({"23": "X"})
    out: dict[str, tuple] = {}
    for chrom, sub in frame.groupby("chr"):
        sub = sub.sort_values("pos")
        out[chrom] = (sub["pos"].to_numpy(dtype=float), sub["cM"].to_numpy(dtype=float))
    return out


def pericentromere_intervals(
    centromere_rows: list[tuple[str, int, int]],
    build: str,
    pad_cm: float = PERICENTROMERE_PAD_CM,
) -> list[tuple[str, int, int]]:
    """Pad each centromere span by ``pad_cm`` cM per side (LDSC pericentromeric).

    For each centromere ``[c0, c1]``, interpolate cM at the endpoints from the
    genetic map, extend by ``+/- pad_cm`` in genetic distance, and invert back to
    base pairs. The padded interval is unioned with ``[c0, c1]`` so the result
    always contains the centromere. The union is what gives the correct
    conservative behavior on the acrocentric chromosomes (13/14/15/21/22): their
    p-arm is absent from the genetic map, so the negative-side padding underflows
    and clamps, and ``min(c0, .)`` keeps the centromere start.
    """
    import numpy as np

    gmap = _load_genetic_map(build)
    out: list[tuple[str, int, int]] = []
    for chrom, start, end in centromere_rows:
        pos, cm = gmap[chrom]
        cm0 = float(np.interp(start, pos, cm))
        cm1 = float(np.interp(end, pos, cm))
        lo_pad = int(round(float(np.interp(cm0 - pad_cm, cm, pos))))
        hi_pad = int(round(float(np.interp(cm1 + pad_cm, cm, pos))))
        lo = max(min(start, lo_pad), 0)
        hi = max(end, hi_pad)
        out.append((chrom, lo, hi))
    return out


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
    for build in ("hg19", "hg38"):
        write_bed(OUT_DIR / f"mhc.{build}.bed", MHC_BROAD[build], _MHC_PROVENANCE[("broad", build)])
        write_bed(OUT_DIR / f"mhc_core.{build}.bed", MHC_CORE[build], _MHC_PROVENANCE[("core", build)])
    for build, fetch in (("hg19", centromeres_hg19), ("hg38", centromeres_hg38)):
        core = fetch()
        # Raw gap kept for reference; the active `centromeres` preset is the +/-3 cM padded region.
        write_bed(OUT_DIR / f"centromeres_core.{build}.bed", core, _CENTROMERE_PROVENANCE[("core", build)])
        write_bed(OUT_DIR / f"centromeres.{build}.bed", pericentromere_intervals(core, build), _CENTROMERE_PROVENANCE[("padded", build)])


if __name__ == "__main__":
    main()
