# Region-Exclusion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add named region presets (`mhc`, `centromeres`) and user-BED region exclusion that drop all SNPs inside the regions before R2/LD computation in `build-ref-panel` and `ldscore`.

**Architecture:** A new private kernel module `src/ldsc/_kernel/regions.py` loads interval sets (curated per-build BED presets + user BEDs) and produces a boolean keep-mask over CHR/POS metadata. Both workflows apply that mask at their existing single SNP-restriction chokepoint, so no numerical kernel changes are needed. Presets are build-resolved (`ldscore` requires an explicit `--exclude-regions-build`; `build-ref-panel` reuses its already-resolved `source_genome_build`); user BEDs are build-agnostic.

**Tech Stack:** Python 3, NumPy, pandas, argparse, `importlib.resources`, pytest.

**Design doc:** `docs/superpowers/specs/2026-06-06-region-exclusion-design.md` (read it first).

**Environment for every command:**
```bash
source /Users/wenbinwu/miniforge3/etc/profile.d/conda.sh && conda activate ldsc3-dev
cd /Users/wenbinwu/Documents_local/Research/SullivanLab/LDSC/repos/ldsc_py3_Jerry_workspace/ldsc_py3_restructured
```
All `pytest` / `git` commands below assume this prefix and CWD. Confirm the branch is `restructure` before starting: `git branch --show-current`.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/ldsc/data/regions/{mhc,centromeres}.{hg19,hg38}.bed` | Curated static region data (4 standard BED files) |
| `tools/regions/build_region_beds.py` | One-time dev curation script that regenerates the BEDs from UCSC (not shipped) |
| `setup.py` | Package the BEDs (`data/regions/*.bed`) |
| `src/ldsc/_kernel/regions.py` | Interval loading (presets + user BED + merge) and the keep-mask engine |
| `src/ldsc/config.py` | New config fields + validation on `ReferencePanelBuildConfig` and `RefPanelConfig` |
| `src/ldsc/_kernel/ref_panel.py` | Apply region exclusion inside `_apply_snp_restriction` (covers both `ldscore` backends) |
| `src/ldsc/ldscore_calculator.py` | `ldscore` CLI flags + arg→`RefPanelConfig` wiring |
| `src/ldsc/ref_panel_builder.py` | `build-ref-panel` CLI flags + `_BuildState` + `_prepare_build_state` + `_build_chromosome` wiring |
| `tests/test_regions.py` | Unit tests for the kernel module |
| `tests/test_config_identifiers.py` | Extend with new config-validation tests |
| `tests/test_ref_panel.py` | `ldscore` region-exclusion behavior |
| `tests/test_ref_panel_builder.py` | `build-ref-panel` two-build region drop |

---

## Task 1: Curate and package the region BED files

**Files:**
- Create: `tools/regions/build_region_beds.py`
- Create: `src/ldsc/data/regions/mhc.hg19.bed`
- Create: `src/ldsc/data/regions/mhc.hg38.bed`
- Create: `src/ldsc/data/regions/centromeres.hg19.bed`
- Create: `src/ldsc/data/regions/centromeres.hg38.bed`
- Modify: `setup.py` (the `package_data` line)
- Test: `tests/test_regions.py` (structural data test)

Centromere coordinates must come from the authoritative UCSC source, not be
transcribed by hand. The curation script below fetches them deterministically;
its committed output is the shipped data. MHC bounds are pinned directly (stable,
widely-used broad windows).

- [ ] **Step 1: Write the curation script**

Create `tools/regions/build_region_beds.py`:

```python
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
    write_bed(OUT_DIR / "mhc.hg19.bed", MHC["hg19"], "MHC broad exclusion window, hg19 (chr6:25-35Mb, 0-based half-open)")
    write_bed(OUT_DIR / "mhc.hg38.bed", MHC["hg38"], "MHC core exclusion window, hg38 (0-based half-open)")
    write_bed(OUT_DIR / "centromeres.hg19.bed", centromeres_hg19(), "UCSC hg19 gap track, type=centromere (0-based half-open)")
    write_bed(OUT_DIR / "centromeres.hg38.bed", centromeres_hg38(), "UCSC hg38 centromeres track, merged per-chromosome span (0-based half-open)")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the curation script to generate the BED files**

Run: `python tools/regions/build_region_beds.py`
Expected output (row counts may vary slightly if UCSC adds models; X present):
```
wrote 1 rows -> .../mhc.hg19.bed
wrote 1 rows -> .../mhc.hg38.bed
wrote 23 rows -> .../centromeres.hg19.bed
wrote 23 rows -> .../centromeres.hg38.bed
```

If the network is unavailable in the build environment, generate the files on a
networked machine and copy them in — they are static committed data.

- [ ] **Step 3: Sanity-inspect the generated files**

Run: `head -3 src/ldsc/data/regions/mhc.hg19.bed src/ldsc/data/regions/centromeres.hg19.bed`
Expected: `mhc.hg19.bed` shows a `# ...` comment then `6<TAB>25000000<TAB>35000000`; centromeres show a comment then `1<TAB><start><TAB><end>` with `start < end`.

- [ ] **Step 4: Write the structural data test**

Create `tests/test_regions.py`:

```python
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
    assert {str(i) for i in range(1, 23)}.issubset(chroms)
    for chrom, start, end in rows:
        assert start < end, f"{chrom}: start {start} not < end {end}"
```

- [ ] **Step 5: Run the data test**

Run: `pytest tests/test_regions.py -v`
Expected: PASS (4 parametrized cases).

> Note: this requires the package to see the new data dir. Because the project is
> installed editable (`pip install -e`), `importlib.resources` reads directly
> from `src/ldsc/data/regions/`; no reinstall is needed. If resolution fails, run
> `pip install -e ".[dev]"` once.

- [ ] **Step 6: Package the BED files in setup.py**

In `setup.py`, change:
```python
package_data={"ldsc": ["data/*.tsv.gz"]},
```
to:
```python
package_data={"ldsc": ["data/*.tsv.gz", "data/regions/*.bed"]},
```

- [ ] **Step 7: Commit**

```bash
git add tools/regions/build_region_beds.py src/ldsc/data/regions/ setup.py tests/test_regions.py
git commit -m "feat(regions): add curated MHC/centromere BED presets and packaging"
```

---

## Task 2: Interval loading in the kernel module

**Files:**
- Create: `src/ldsc/_kernel/regions.py`
- Test: `tests/test_regions.py` (extend)

- [ ] **Step 1: Write failing tests for interval loading**

Append to `tests/test_regions.py`:

```python
import numpy as np

from ldsc._kernel import regions
from ldsc.errors import LDSCConfigError, LDSCInputError, LDSCUsageError


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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/test_regions.py -k "load_preset or load_bed or merge_intervals" -v`
Expected: FAIL with `ModuleNotFoundError` / `AttributeError: module 'ldsc._kernel.regions' has no attribute ...`.

- [ ] **Step 3: Implement interval loading**

Create `src/ldsc/_kernel/regions.py`:

```python
"""Region-exclusion interval loading and SNP keep-masking (private kernel).

Intervals use the standard BED convention: 0-based, half-open ``[start, end)``.
SNP positions in panel metadata are 1-based (PLINK convention). A 1-based
position ``p`` is excluded iff some interval contains it, i.e.
``start < p <= end``. See ``region_exclusion_keep_mask``.
"""
from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
from os import PathLike
from typing import Sequence

import numpy as np
import pandas as pd

from ..chromosome_inference import normalize_chromosome
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError

REGION_PRESETS: frozenset[str] = frozenset({"mhc", "centromeres"})
_PRESET_BUILDS: frozenset[str] = frozenset({"hg19", "hg38"})


@dataclass(frozen=True)
class RegionIntervals:
    """Per-chromosome exclusion intervals in 0-based half-open BED coordinates.

    Parameters
    ----------
    intervals : dict of str to numpy.ndarray
        Maps normalized chromosome token to a sorted, disjoint ``(k, 2)`` int64
        array of ``[start, end)`` rows. An empty mapping excludes nothing.
    source_labels : tuple of str
        Human-readable provenance for logging, e.g. ``("preset:mhc[hg19]",)``.
    """

    intervals: dict[str, np.ndarray]
    source_labels: tuple[str, ...]


def _coalesce(rows: list[tuple[int, int]]) -> np.ndarray:
    """Sort and union possibly-overlapping intervals into disjoint sorted rows."""
    if not rows:
        return np.empty((0, 2), dtype=np.int64)
    merged: list[list[int]] = []
    for start, end in sorted(rows):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return np.asarray(merged, dtype=np.int64).reshape(-1, 2)


def _intervals_from_rows(rows: Sequence[tuple[str, int, int]]) -> dict[str, np.ndarray]:
    """Group ``(chrom, start, end)`` rows by normalized chromosome and coalesce."""
    grouped: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in rows:
        grouped.setdefault(normalize_chromosome(chrom), []).append((int(start), int(end)))
    return {chrom: _coalesce(pairs) for chrom, pairs in grouped.items()}


def _parse_bed_text(text: str, label: str) -> list[tuple[str, int, int]]:
    """Parse standard BED text into ``(chrom, start, end)`` rows."""
    rows: list[tuple[str, int, int]] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) < 3:
            raise LDSCInputError(
                f"Region BED parsing failed in '{label}' line {lineno}: expected at least 3 "
                f"columns (chrom start end), found {len(fields)}. Most likely the file is not "
                "tab/space-delimited BED. Provide a standard BED file."
            )
        chrom, start_token, end_token = fields[0], fields[1], fields[2]
        try:
            start, end = int(start_token), int(end_token)
        except ValueError as exc:
            raise LDSCInputError(
                f"Region BED parsing failed in '{label}' line {lineno}: start/end "
                f"('{start_token}', '{end_token}') are not integers. Provide 0-based "
                "half-open integer BED coordinates."
            ) from exc
        if start >= end:
            raise LDSCInputError(
                f"Region BED parsing failed in '{label}' line {lineno}: start {start} is not "
                f"less than end {end}. BED intervals are 0-based half-open with start < end."
            )
        rows.append((chrom, start, end))
    return rows


def load_preset_intervals(names: Sequence[str], build: str) -> RegionIntervals:
    """Read packaged preset BEDs for ``names`` in genome ``build``.

    Parameters
    ----------
    names : sequence of str
        Preset names; each must be in :data:`REGION_PRESETS`.
    build : str
        Genome build, either ``"hg19"`` or ``"hg38"``.

    Returns
    -------
    RegionIntervals
        Union of the requested presets' intervals.
    """
    if build not in _PRESET_BUILDS:
        raise LDSCConfigError(
            f"Region preset loading received unsupported genome build '{build}'. Most likely "
            "an invalid --exclude-regions-build was passed. Use 'hg19' or 'hg38'."
        )
    unknown = sorted(set(names) - REGION_PRESETS)
    if unknown:
        raise LDSCUsageError(
            f"Region preset loading received unknown region preset(s): {unknown}. Most likely "
            f"a name was misspelled. Valid presets are {sorted(REGION_PRESETS)}."
        )
    rows: list[tuple[str, int, int]] = []
    labels: list[str] = []
    for name in names:
        path = resources.files("ldsc").joinpath("data", "regions", f"{name}.{build}.bed")
        try:
            text = path.read_text(encoding="utf-8")
        except (FileNotFoundError, OSError) as exc:
            raise LDSCInternalError(
                f"Region preset loading could not read packaged BED for '{name}.{build}'. "
                "Most likely the package data was not installed. Reinstall the package or "
                "report the traceback."
            ) from exc
        rows.extend(_parse_bed_text(text, label=f"{name}.{build}.bed"))
        labels.append(f"preset:{name}[{build}]")
    return RegionIntervals(intervals=_intervals_from_rows(rows), source_labels=tuple(labels))


def load_bed_intervals(paths: Sequence[str | PathLike[str]]) -> RegionIntervals:
    """Read one or more user BED files, applied as-is on panel CHR/POS."""
    rows: list[tuple[str, int, int]] = []
    labels: list[str] = []
    for path in paths:
        path_str = str(path)
        try:
            text = open(path_str, encoding="utf-8").read()
        except (FileNotFoundError, OSError) as exc:
            raise LDSCInputError(
                f"Region BED loading could not read user BED file '{path_str}'. Most likely "
                "the path is wrong or unreadable. Provide an existing standard BED file."
            ) from exc
        rows.extend(_parse_bed_text(text, label=path_str))
        labels.append(f"bed:{path_str}")
    return RegionIntervals(intervals=_intervals_from_rows(rows), source_labels=tuple(labels))


def merge_intervals(*groups: RegionIntervals) -> RegionIntervals:
    """Union several RegionIntervals into one (coalescing overlaps per chromosome)."""
    rows: list[tuple[str, int, int]] = []
    labels: list[str] = []
    for group in groups:
        for chrom, bounds in group.intervals.items():
            for start, end in bounds.tolist():
                rows.append((chrom, int(start), int(end)))
        labels.extend(group.source_labels)
    return RegionIntervals(intervals=_intervals_from_rows(rows), source_labels=tuple(labels))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/test_regions.py -k "load_preset or load_bed or merge_intervals" -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ldsc/_kernel/regions.py tests/test_regions.py
git commit -m "feat(regions): interval loading for presets and user BEDs"
```

---

## Task 3: The keep-mask engine

**Files:**
- Modify: `src/ldsc/_kernel/regions.py`
- Test: `tests/test_regions.py` (extend)

- [ ] **Step 1: Write failing tests for the keep-mask (boundary-focused)**

Append to `tests/test_regions.py`:

```python
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/test_regions.py -k keep_mask -v`
Expected: FAIL with `AttributeError: module 'ldsc._kernel.regions' has no attribute 'region_exclusion_keep_mask'`.

- [ ] **Step 3: Implement the keep-mask**

Append to `src/ldsc/_kernel/regions.py`:

```python
def region_exclusion_keep_mask(
    metadata: pd.DataFrame,
    intervals: RegionIntervals,
    *,
    chr_col: str = "CHR",
    pos_col: str = "POS",
) -> np.ndarray:
    """Return a boolean KEEP mask over ``metadata`` rows.

    A 1-based SNP position ``p`` is excluded iff some interval ``[start, end)``
    (0-based half-open) satisfies ``start < p <= end``. Intervals are assumed
    sorted and disjoint per chromosome (guaranteed by the loaders). Rows whose
    chromosome has no intervals are kept. Returns all-True when ``intervals`` is
    empty.

    Parameters
    ----------
    metadata : pandas.DataFrame
        Must contain ``chr_col`` and ``pos_col``. Positions are 1-based.
    intervals : RegionIntervals
        Exclusion intervals to apply.
    chr_col, pos_col : str, optional
        Column names for chromosome and 1-based position. Defaults ``"CHR"`` /
        ``"POS"``.

    Returns
    -------
    numpy.ndarray
        Boolean array, ``True`` to keep the row.
    """
    n = len(metadata)
    keep = np.ones(n, dtype=bool)
    if not intervals.intervals or n == 0:
        return keep
    chrom_tokens = metadata[chr_col].astype(str).map(normalize_chromosome).to_numpy()
    pos = pd.to_numeric(metadata[pos_col], errors="raise").to_numpy()
    for chrom, bounds in intervals.intervals.items():
        sel = chrom_tokens == chrom
        if not sel.any():
            continue
        p = pos[sel]
        starts = bounds[:, 0]
        ends = bounds[:, 1]
        # Candidate interval for each p: the last one with start < p.
        idx = np.searchsorted(starts, p, side="left") - 1
        excluded = np.zeros(p.shape, dtype=bool)
        valid = idx >= 0
        excluded[valid] = ends[idx[valid]] >= p[valid]
        local = keep[sel]
        local[excluded] = False
        keep[sel] = local
    return keep
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest tests/test_regions.py -k keep_mask -v`
Expected: PASS.

- [ ] **Step 5: Run the full module test file**

Run: `pytest tests/test_regions.py -v`
Expected: PASS (all data, loading, and mask tests).

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/_kernel/regions.py tests/test_regions.py
git commit -m "feat(regions): half-open keep-mask engine over CHR/POS"
```

---

## Task 4: Config fields and validation

**Files:**
- Modify: `src/ldsc/config.py` (`RefPanelConfig` at line ~360; `ReferencePanelBuildConfig` at line ~513)
- Test: `tests/test_config_identifiers.py`

- [ ] **Step 1: Write failing config-validation tests**

Append to `tests/test_config_identifiers.py`:

```python
import pytest

from ldsc.config import RefPanelConfig, ReferencePanelBuildConfig
from ldsc.errors import LDSCConfigError


def test_refpanelconfig_preset_requires_build():
    with pytest.raises(LDSCConfigError, match="exclude-regions-build"):
        RefPanelConfig(backend="plink", plink_prefix="x", exclude_regions=("mhc",))


def test_refpanelconfig_preset_with_build_ok():
    cfg = RefPanelConfig(
        backend="plink", plink_prefix="x", exclude_regions=("mhc",), exclude_regions_build="hg19"
    )
    assert cfg.exclude_regions == ("mhc",)
    assert cfg.exclude_regions_build == "hg19"


def test_refpanelconfig_unknown_preset_raises():
    with pytest.raises(LDSCConfigError, match="unknown region preset"):
        RefPanelConfig(backend="plink", plink_prefix="x", exclude_regions=("telomeres",), exclude_regions_build="hg19")


def test_refpanelconfig_user_bed_without_build_ok():
    cfg = RefPanelConfig(backend="plink", plink_prefix="x", exclude_regions_bed=("/tmp/a.bed",))
    assert cfg.exclude_regions_bed == ("/tmp/a.bed",)


def test_build_config_unknown_preset_raises():
    with pytest.raises(LDSCConfigError, match="unknown region preset"):
        ReferencePanelBuildConfig(plink_prefix="x", output_dir="o", ld_wind_snps=1, exclude_regions=("foo",))


def test_build_config_preset_no_build_field_required():
    cfg = ReferencePanelBuildConfig(
        plink_prefix="x", output_dir="o", ld_wind_snps=1, exclude_regions=("mhc", "centromeres")
    )
    assert cfg.exclude_regions == ("mhc", "centromeres")
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_config_identifiers.py -k "exclude or preset or region or user_bed" -v`
Expected: FAIL with `TypeError: __init__() got an unexpected keyword argument 'exclude_regions'`.

- [ ] **Step 3: Add a shared validation helper**

In `src/ldsc/config.py`, near the other `_*_message` helpers at the top of the
module, add:

```python
def _validate_region_presets(class_name: str, names: tuple[str, ...]) -> None:
    """Raise if any preset name is outside the supported region menu."""
    from ._kernel.regions import REGION_PRESETS

    unknown = sorted(set(names) - REGION_PRESETS)
    if unknown:
        raise LDSCConfigError(
            f"Could not construct {class_name}: unknown region preset(s) {unknown}. "
            f"Most likely a name was misspelled. Valid presets are {sorted(REGION_PRESETS)}."
        )
```

> The import is function-local to avoid any import-time cycle between `config`
> and `_kernel.regions` (which imports `errors`, not `config`).

- [ ] **Step 4: Add fields and validation to `RefPanelConfig`**

In `RefPanelConfig` (after `use_hm3_ref_panel_snps: bool = False`), add:
```python
    exclude_regions: tuple[str, ...] = ()
    exclude_regions_bed: tuple[str, ...] = ()
    exclude_regions_build: Literal["hg19", "hg38"] | None = None
```

> `Literal` is already imported in `config.py` (used by `GenomeBuild`). If not,
> add `from typing import Literal`.

At the end of `RefPanelConfig.__post_init__`, add:
```python
        object.__setattr__(self, "exclude_regions", tuple(self.exclude_regions))
        object.__setattr__(
            self,
            "exclude_regions_bed",
            tuple(_normalize_optional_path(token) for token in self.exclude_regions_bed if token),
        )
        _validate_region_presets("RefPanelConfig", self.exclude_regions)
        if self.exclude_regions_build not in {None, "hg19", "hg38"}:
            raise LDSCConfigError(
                _invalid_choice_message(
                    "RefPanelConfig", "exclude_regions_build", self.exclude_regions_build, "None, 'hg19', or 'hg38'"
                )
            )
        if self.exclude_regions and self.exclude_regions_build is None:
            raise LDSCConfigError(
                "RefPanelConfig received exclude_regions presets without exclude_regions_build. "
                "Most likely --exclude-regions was passed without --exclude-regions-build. "
                "Region presets are build-specific; pass --exclude-regions-build hg19 or hg38."
            )
```

- [ ] **Step 5: Add fields and validation to `ReferencePanelBuildConfig`**

In `ReferencePanelBuildConfig` (after `overwrite: bool = False`), add:
```python
    exclude_regions: tuple[str, ...] = ()
    exclude_regions_bed: tuple[str, ...] = ()
```

At the end of `ReferencePanelBuildConfig.__post_init__`, add:
```python
        object.__setattr__(self, "exclude_regions", tuple(self.exclude_regions))
        object.__setattr__(
            self,
            "exclude_regions_bed",
            tuple(_normalize_optional_path(token) for token in self.exclude_regions_bed if token),
        )
        _validate_region_presets("ReferencePanelBuildConfig", self.exclude_regions)
```

> No build field here: `build-ref-panel` reuses the resolved `source_genome_build`.

- [ ] **Step 6: Run the config tests to verify they pass**

Run: `pytest tests/test_config_identifiers.py -k "exclude or preset or region or user_bed" -v`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/config.py tests/test_config_identifiers.py
git commit -m "feat(regions): config fields and validation for region exclusion"
```

---

## Task 5: Apply region exclusion in the ldscore reference panel

**Files:**
- Modify: `src/ldsc/_kernel/ref_panel.py` (`RefPanel.__init__` ~line 209; `_apply_snp_restriction` ~line 267)
- Test: `tests/test_ref_panel.py`

- [ ] **Step 1: Write a failing behavior test**

This uses the real on-disk chr22 PLINK fixture (`hm3_chr22_subset`, hg38) and the
`PlinkRefPanel(global_config, spec).load_metadata("22")` idiom already used by
`test_available_chromosomes_and_metadata` in this file. It discovers a real SNP
position from a no-exclusion baseline, then excludes exactly that position via a
user BED — no hard-coded coordinates. Append to `tests/test_ref_panel.py`:

```python
from pathlib import Path

import pytest

from ldsc.config import GlobalConfig, RefPanelConfig
from ldsc._kernel.ref_panel import PlinkRefPanel

_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


def _chr22_available() -> bool:
    return all(Path(str(_CHR22) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def test_ldscore_panel_excludes_user_bed_region(tmp_path):
    if not _chr22_available():
        pytest.skip("chr22 PLINK fixture unavailable; run tests/fixtures/generate_minimal_external_resources.py")
    gc = GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38")

    base_panel = PlinkRefPanel(gc, RefPanelConfig(backend="plink", plink_prefix=str(_CHR22)))
    base_meta = base_panel.load_metadata("22")
    target_pos = int(base_meta["POS"].iloc[0])

    # BED [target_pos-1, target_pos) (0-based half-open) excludes 1-based POS == target_pos.
    bed = tmp_path / "exclude.bed"
    bed.write_text(f"22\t{target_pos - 1}\t{target_pos}\n", encoding="utf-8")

    panel = PlinkRefPanel(
        gc, RefPanelConfig(backend="plink", plink_prefix=str(_CHR22), exclude_regions_bed=(str(bed),))
    )
    meta = panel.load_metadata("22")
    assert target_pos not in set(meta["POS"])
    assert len(meta) == len(base_meta) - int((base_meta["POS"] == target_pos).sum())
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ref_panel.py -k excludes_user_bed_region -v`
Expected: FAIL — `target_pos` still present (region exclusion not wired yet), or
`TypeError` if `exclude_regions_bed` is not yet a config field (it is, from Task 4).

- [ ] **Step 3: Add an interval cache to `RefPanel.__init__`**

In `RefPanel.__init__` (where `self._metadata_cache` is set), add:
```python
        self._region_intervals_cache: "regions.RegionIntervals | None" = None
```
And add the import near the top of `ref_panel.py` with the other `from .` imports:
```python
from . import regions
from ..path_resolution import resolve_plink_prefix, resolve_plink_prefix_group, resolve_scalar_path
```

> `resolve_plink_prefix`/`resolve_plink_prefix_group` are already imported on the
> existing line; just add `resolve_scalar_path` to that existing import (do not
> duplicate the line).

- [ ] **Step 4: Refactor `_apply_snp_restriction` to also apply region exclusion**

Replace the body of `_apply_snp_restriction` so the keep-restriction and region
exclusion are separate steps and the region step always runs:

```python
    def _apply_snp_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Filter metadata to the keep-restriction, then drop excluded regions."""
        metadata = self._apply_keep_restriction(metadata)
        metadata = self._apply_region_exclusion(metadata)
        return metadata

    def _apply_keep_restriction(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Filter metadata rows to ``RefPanelConfig.ref_panel_snps_file`` when set."""
        restrict_path = packaged_hm3_curated_map_path() if self.spec.use_hm3_ref_panel_snps else self.spec.ref_panel_snps_file
        if restrict_path is None or len(metadata) == 0:
            return metadata
        if self.spec.use_hm3_ref_panel_snps:
            LOGGER.info(f"Using packaged curated HM3 map for reference-panel SNP restriction: {restrict_path}.")
        restriction = read_snp_restriction_keys(
            restrict_path,
            self.global_config.snp_identifier,
            genome_build=self.global_config.genome_build,
            logger=LOGGER,
        )
        keep = restriction_membership_mask(
            metadata,
            restriction,
            self.global_config.snp_identifier,
            context=f"reference-panel restriction matching for {restrict_path}",
        )
        return metadata.loc[keep].reset_index(drop=True)

    def _region_intervals(self) -> regions.RegionIntervals:
        """Resolve and cache region-exclusion intervals for this panel."""
        if self._region_intervals_cache is not None:
            return self._region_intervals_cache
        groups: list[regions.RegionIntervals] = []
        if self.spec.exclude_regions:
            groups.append(
                regions.load_preset_intervals(list(self.spec.exclude_regions), self.spec.exclude_regions_build)
            )
        if self.spec.exclude_regions_bed:
            paths = [resolve_scalar_path(token, suffixes=("", ".bed"), label="region exclusion BED") for token in self.spec.exclude_regions_bed]
            groups.append(regions.load_bed_intervals(paths))
        merged = regions.merge_intervals(*groups) if groups else regions.RegionIntervals(intervals={}, source_labels=())
        self._region_intervals_cache = merged
        return merged

    def _apply_region_exclusion(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """Drop SNPs inside any configured exclusion region (CHR/POS based)."""
        intervals = self._region_intervals()
        if not intervals.intervals or len(metadata) == 0:
            return metadata
        pos_col = "POS" if "POS" in metadata.columns else "BP"
        keep = regions.region_exclusion_keep_mask(metadata, intervals, pos_col=pos_col)
        dropped = int((~keep).sum())
        if dropped:
            LOGGER.info(
                f"Region exclusion dropped {dropped} reference-panel SNPs via "
                f"{', '.join(intervals.source_labels)}."
            )
        return metadata.loc[keep].reset_index(drop=True)
```

- [ ] **Step 5: Run the behavior test to verify it passes**

Run: `pytest tests/test_ref_panel.py -k excludes_user_bed_region -v`
Expected: PASS.

- [ ] **Step 6: Run the full ref_panel test file (regression check)**

Run: `pytest tests/test_ref_panel.py -v`
Expected: PASS (no existing behavior broken; the refactor preserves the
keep-restriction path).

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/_kernel/ref_panel.py tests/test_ref_panel.py
git commit -m "feat(regions): apply region exclusion in ldscore reference panel"
```

---

## Task 6: ldscore CLI flags and wiring

**Files:**
- Modify: `src/ldsc/ldscore_calculator.py` (`build_parser` ~line 801; `_normalize_run_args` ~line 1239; `_ref_panel_from_args` ~line 1432)
- Test: `tests/test_ldscore_workflow.py`

- [ ] **Step 1: Write a failing CLI/workflow test**

Append to `tests/test_ldscore_workflow.py`:

```python
import pytest

from ldsc.ldscore_calculator import build_parser, _ref_panel_from_args, _normalize_run_args
from ldsc.errors import LDSCConfigError


def test_ldscore_parser_accepts_region_flags():
    parser = build_parser()
    args = parser.parse_args(
        ["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
         "--exclude-regions", "mhc,centromeres", "--exclude-regions-build", "hg19",
         "--exclude-regions-bed", "/tmp/a.bed,/tmp/b.bed"]
    )
    assert args.exclude_regions == "mhc,centromeres"
    assert args.exclude_regions_build == "hg19"
    assert args.exclude_regions_bed == "/tmp/a.bed,/tmp/b.bed"


def test_ldscore_preset_without_build_rejected():
    parser = build_parser()
    args = parser.parse_args(
        ["--output-dir", "o", "--plink-prefix", "p", "--ld-wind-cm", "1",
         "--snp-identifier", "chr_pos_allele_aware", "--genome-build", "hg19",
         "--exclude-regions", "mhc"]
    )
    normalized, global_config = _normalize_run_args(args)
    with pytest.raises(LDSCConfigError, match="exclude-regions-build"):
        _ref_panel_from_args(normalized, global_config)
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ldscore_workflow.py -k "region_flags or preset_without_build" -v`
Expected: FAIL — argparse rejects unknown `--exclude-regions`.

- [ ] **Step 3: Add the parser flags**

In `ldscore_calculator.py` `build_parser`, after the `--use-hm3-ref-panel-snps`
argument block, add:
```python
    parser.add_argument(
        "--exclude-regions",
        default=None,
        help="Comma-separated curated region presets to exclude before LD computation (choices: mhc, centromeres). Requires --exclude-regions-build.",
    )
    parser.add_argument(
        "--exclude-regions-build",
        choices=("hg19", "hg38"),
        default=None,
        help="Genome build of the panel coordinates, used to select preset BEDs. Required whenever --exclude-regions is given.",
    )
    parser.add_argument(
        "--exclude-regions-bed",
        default=None,
        help="Comma-separated user BED file path tokens of regions to exclude, applied as-is on panel CHR/POS (0-based half-open).",
    )
```

> Preset-name validation happens in `RefPanelConfig.__post_init__`, so the menu
> is not duplicated as argparse `choices` here (a comma-joined value can't use
> `choices`). The help text documents the menu.

- [ ] **Step 4: Normalize the new tokens in `_normalize_run_args`**

In `_normalize_run_args`, alongside the existing `ref_panel_snps_file`
normalization, add:
```python
    normalized_args.exclude_regions = getattr(args, "exclude_regions", None)
    normalized_args.exclude_regions_build = getattr(args, "exclude_regions_build", None)
    normalized_args.exclude_regions_bed = getattr(args, "exclude_regions_bed", None)
```

> If `_normalize_run_args` constructs a fresh namespace, copy these through the
> same way the existing fields are copied. If it mutates `args` in place, the
> `getattr` defaults already cover programmatic namespaces.

- [ ] **Step 5: Pass the fields into `RefPanelConfig` in `_ref_panel_from_args`**

In `_ref_panel_from_args`, both `RefPanelConfig(...)` construction sites (the
parquet branch ~line 1439 and the PLINK branch ~line 1450) gain:
```python
            exclude_regions=tuple(split_cli_path_tokens(getattr(args, "exclude_regions", None))),
            exclude_regions_bed=tuple(split_cli_path_tokens(getattr(args, "exclude_regions_bed", None))),
            exclude_regions_build=getattr(args, "exclude_regions_build", None),
```

> `split_cli_path_tokens` is already imported in this module and turns a
> comma-separated string (or `None`) into a tuple of tokens.

- [ ] **Step 6: Run the workflow tests to verify they pass**

Run: `pytest tests/test_ldscore_workflow.py -k "region_flags or preset_without_build" -v`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/ldsc/ldscore_calculator.py tests/test_ldscore_workflow.py
git commit -m "feat(regions): ldscore CLI flags and RefPanelConfig wiring"
```

---

## Task 7: build-ref-panel wiring (source-build exclusion, both emitted builds)

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py` (`_BuildState` ~line 129; `_prepare_build_state` ~line 555; `_build_chromosome` ~line 734, near the restriction mask at line 789)
- Test: `tests/test_ref_panel_builder.py`

- [ ] **Step 1: Write a failing behavior test**

The real chr22 fixture (`hm3_chr22_subset`, hg38) emits a single build, which is
enough to prove `_build_chromosome` applies the exclusion. The "removed from
**both** emitted builds" property is guaranteed by construction (the mask runs on
`chrom_metadata` *before* liftover, so one filtered `keep_snps` feeds every
emitted build) and is additionally covered by the kernel mask test in Task 3. This
test models on the real-fixture builder smoke test
`test_hm3_chr22_subset_builds_source_only_without_liftover_chain_or_map` in this
file. Append to `tests/test_ref_panel_builder.py`:

```python
import gzip
from pathlib import Path

import pandas as pd
import pytest

from ldsc import ref_panel_builder
from ldsc.config import set_global_config, GlobalConfig

_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


def _chr22_available() -> bool:
    return all(Path(str(_CHR22) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def _read_meta(path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return pd.read_csv(handle, sep="\t", comment="#")


def test_build_ref_panel_excludes_user_bed_region(tmp_path):
    if not _chr22_available():
        pytest.skip("chr22 PLINK fixture unavailable; run tests/fixtures/generate_minimal_external_resources.py")
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg38"))

    base = ref_panel_builder.run_build_ref_panel(
        plink_prefix=str(_CHR22), source_genome_build="hg38",
        output_dir=str(tmp_path / "base"), ld_wind_snps=10, snp_batch_size=64,
    )
    base_meta = _read_meta(base.output_paths["meta_hg38"][0])
    target_pos = int(base_meta["POS"].iloc[0])

    bed = tmp_path / "exclude.bed"
    bed.write_text(f"22\t{target_pos - 1}\t{target_pos}\n", encoding="utf-8")

    excluded = ref_panel_builder.run_build_ref_panel(
        plink_prefix=str(_CHR22), source_genome_build="hg38",
        output_dir=str(tmp_path / "excl"), ld_wind_snps=10, snp_batch_size=64,
        exclude_regions_bed=(str(bed),),
    )
    excl_meta = _read_meta(excluded.output_paths["meta_hg38"][0])
    assert target_pos not in set(excl_meta["POS"])
    assert len(excl_meta) == len(base_meta) - int((base_meta["POS"] == target_pos).sum())
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ref_panel_builder.py -k excludes_user_bed_region -v`
Expected: FAIL — `target_pos` still present in the excluded build's sidecar.

- [ ] **Step 3: Add the interval field to `_BuildState`**

In the `_BuildState` dataclass, add a field:
```python
    region_intervals: "kernel_regions.RegionIntervals | None" = None
```
Add the import near the other `from ._kernel import ... as ...` imports at the top
of `ref_panel_builder.py`:
```python
from ._kernel import regions as kernel_regions
```

- [ ] **Step 4: Load intervals in `_prepare_build_state`**

`source_genome_build` is concrete here (resolved at line 422 before this call).
At the end of `_prepare_build_state`, build the intervals and pass them into the
returned `_BuildState(...)` by adding this before the `return` and a
`region_intervals=` argument to the constructor:

```python
        region_groups: list[kernel_regions.RegionIntervals] = []
        if config.exclude_regions:
            region_groups.append(kernel_regions.load_preset_intervals(list(config.exclude_regions), source_build))
        if config.exclude_regions_bed:
            bed_paths = [
                resolve_scalar_path(token, suffixes=("", ".bed"), label="region exclusion BED")
                for token in config.exclude_regions_bed
            ]
            region_groups.append(kernel_regions.load_bed_intervals(bed_paths))
        region_intervals = (
            kernel_regions.merge_intervals(*region_groups)
            if region_groups
            else kernel_regions.RegionIntervals(intervals={}, source_labels=())
        )
        if region_intervals.intervals:
            LOGGER.info(f"Reference-panel region exclusion active: {', '.join(region_intervals.source_labels)}.")
```

Then add `region_intervals=region_intervals,` to the `_BuildState(...)` constructor
call.

- [ ] **Step 5: Apply the mask in `_build_chromosome`**

In `_build_chromosome`, immediately after `chrom_metadata` is constructed
(around line 782, after `chrom_metadata["_plink_row_index"] = ...`) and before the
restriction-mask block at line 784, fold in the region keep-mask:

```python
        if build_state.region_intervals is not None and build_state.region_intervals.intervals:
            region_keep = kernel_regions.region_exclusion_keep_mask(
                chrom_metadata, build_state.region_intervals, pos_col="POS"
            )
            dropped = int((~region_keep).sum())
            if dropped:
                LOGGER.info(f"Region exclusion dropped {dropped} SNPs on chromosome {chrom} (source build).")
            chrom_metadata = chrom_metadata.loc[region_keep].reset_index(drop=True)
            keep_snps = chrom_metadata["_plink_row_index"].to_numpy(dtype=int)
            if len(keep_snps) == 0:
                _write_dropped_sidecar(_empty_unified_drop_frame(), sidecar_path, chrom)
                LOGGER.info(f"Skipping chromosome {chrom}: no SNPs remain after region exclusion.")
                return None
```

> This runs on source-build `POS`, before identity cleanup and liftover, so the
> dropped SNPs are excluded from every emitted build. `keep_snps` is recomputed
> from the filtered `chrom_metadata` so the existing restriction-mask block (which
> indexes `chrom_metadata.reset_index(drop=True)`) stays consistent.

- [ ] **Step 6: Run the behavior test to verify it passes**

Run: `pytest tests/test_ref_panel_builder.py -k excludes_user_bed_region -v`
Expected: PASS.

- [ ] **Step 7: Run the full builder test file (regression check)**

Run: `pytest tests/test_ref_panel_builder.py -v`
Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(regions): build-ref-panel source-build exclusion across emitted builds"
```

---

## Task 8: build-ref-panel CLI flags and wiring

**Files:**
- Modify: `src/ldsc/ref_panel_builder.py` (`build_parser` ~line 1682; `config_from_args` ~line 1798)
- Test: `tests/test_ref_panel_builder.py`

- [ ] **Step 1: Write a failing parser test**

Append to `tests/test_ref_panel_builder.py`:

```python
from ldsc.ref_panel_builder import build_parser, config_from_args


def test_build_ref_panel_parser_region_flags():
    parser = build_parser()
    args = parser.parse_args(
        ["--plink-prefix", "p", "--output-dir", "o", "--ld-wind-kb", "1000",
         "--exclude-regions", "mhc,centromeres", "--exclude-regions-bed", "/tmp/x.bed"]
    )
    build_config, _global = config_from_args(args)
    assert build_config.exclude_regions == ("mhc", "centromeres")
    assert build_config.exclude_regions_bed == ("/tmp/x.bed",)
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_ref_panel_builder.py -k parser_region_flags -v`
Expected: FAIL — argparse rejects unknown `--exclude-regions`.

- [ ] **Step 3: Add the parser flags**

In `ref_panel_builder.py` `build_parser`, after the `--maf-min` argument, add:
```python
    parser.add_argument(
        "--exclude-regions",
        default=None,
        help="Comma-separated curated region presets to exclude before R2 computation (choices: mhc, centromeres). Resolved in the panel's source build and removed from every emitted build.",
    )
    parser.add_argument(
        "--exclude-regions-bed",
        default=None,
        help="Comma-separated user BED file path tokens of regions to exclude, applied as-is on source-build CHR/POS (0-based half-open).",
    )
```

> No `--exclude-regions-build`: `build-ref-panel` reuses `source_genome_build`.

- [ ] **Step 4: Map the tokens in `config_from_args`**

In `config_from_args`, add to the `ReferencePanelBuildConfig(...)` constructor:
```python
        exclude_regions=tuple(split_cli_path_tokens(getattr(args, "exclude_regions", None))),
        exclude_regions_bed=tuple(split_cli_path_tokens(getattr(args, "exclude_regions_bed", None))),
```

> Ensure `split_cli_path_tokens` is imported in `ref_panel_builder.py`. If it is
> not already imported from `.path_resolution`, add it to that import line.

- [ ] **Step 5: Run the parser test to verify it passes**

Run: `pytest tests/test_ref_panel_builder.py -k parser_region_flags -v`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/ldsc/ref_panel_builder.py tests/test_ref_panel_builder.py
git commit -m "feat(regions): build-ref-panel CLI flags for region exclusion"
```

---

## Task 9: End-to-end ldscore exclusion + full-suite verification

**Files:**
- Test: `tests/test_ldscore_workflow.py` (one end-to-end test via the public
  `run_ldscore` entry on the real chr22 fixture)
- No production changes expected (this task is verification; if it surfaces a
  bug, fix under TDD back in the relevant task).

`run_ldscore` accepts CLI-equivalent keyword names (see its docstring). This test
drives the full public path (`run_ldscore` → `RefPanelLoader` →
`_apply_region_exclusion` → aggregated result) on the real chr22 fixture, using
the discover-then-exclude pattern so no coordinates are hard-coded.

- [ ] **Step 1: Write an end-to-end ldscore exclusion test**

Append to `tests/test_ldscore_workflow.py`:

```python
import gzip
from pathlib import Path

import pytest

from ldsc.config import set_global_config, GlobalConfig
from ldsc.ldscore_calculator import run_ldscore
from ldsc._kernel.ref_panel import PlinkRefPanel
from ldsc.config import RefPanelConfig

_CHR22 = Path(__file__).resolve().parent / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


def _chr22_available() -> bool:
    return all(Path(str(_CHR22) + ext).exists() for ext in (".bed", ".bim", ".fam"))


def test_ldscore_end_to_end_drops_user_bed_snp(tmp_path):
    if not _chr22_available():
        pytest.skip("chr22 PLINK fixture unavailable; run tests/fixtures/generate_minimal_external_resources.py")
    set_global_config(GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"))

    # Discover a real emitted POS from a no-exclusion panel load.
    base_meta = PlinkRefPanel(
        GlobalConfig(snp_identifier="chr_pos_allele_aware", genome_build="hg38"),
        RefPanelConfig(backend="plink", plink_prefix=str(_CHR22)),
    ).load_metadata("22")
    target_pos = int(base_meta["POS"].iloc[0])

    bed = tmp_path / "exclude.bed"
    bed.write_text(f"22\t{target_pos - 1}\t{target_pos}\n", encoding="utf-8")

    result = run_ldscore(
        plink_prefix=str(_CHR22),
        output_dir=str(tmp_path / "ld"),
        ld_wind_snps=10,
        exclude_regions_bed=(str(bed),),
    )
    assert target_pos not in set(result.baseline_table["POS"])
```

- [ ] **Step 2: Run the integration test**

Run: `pytest tests/test_ldscore_workflow.py -k end_to_end_drops_user_bed_snp -v`
Expected: PASS. If it fails because `run_ldscore` does not thread the keyword,
trace it back to Task 6's wiring and fix under TDD.

- [ ] **Step 3: Run the entire test suite**

Run: `pytest -q`
Expected: PASS (all green). Investigate and fix any regression before continuing.

- [ ] **Step 4: Verify CLI help surfaces the flags**

Run: `ldsc build-ref-panel --help | grep exclude-regions`
Expected: shows `--exclude-regions` and `--exclude-regions-bed`, and **no**
`--exclude-regions-build`.

Run: `ldsc ldscore --help | grep exclude-regions`
Expected: shows `--exclude-regions`, `--exclude-regions-build`, and
`--exclude-regions-bed`.

- [ ] **Step 5: Commit**

```bash
git add tests/test_ldscore_workflow.py
git commit -m "test(regions): end-to-end ldscore user-BED exclusion"
```

---

## Task 10: Documentation

**Files:**
- Modify: `docs/current/data-flow.md` and/or `docs/current/class-and-features.md`
  (wherever the SNP-restriction chokepoints are described)
- Modify: `design_map.md` (if present at repo root)
- Modify: `tutorials/` (only if a build-ref-panel/ldscore tutorial enumerates flags)

- [ ] **Step 1: Document the feature in the current-design docs**

Add a short subsection to the most relevant `docs/current/*.md` describing:
- the two new inputs (curated presets vs user BED) and their build semantics;
- the two chokepoints (`RefPanel._apply_region_exclusion`,
  `ReferencePanelBuilder._build_chromosome`);
- the asymmetry (`ldscore` requires `--exclude-regions-build`; `build-ref-panel`
  reuses `source_genome_build`).
Cite the design doc `docs/superpowers/specs/2026-06-06-region-exclusion-design.md`.

- [ ] **Step 2: Update `design_map.md` if it exists**

Run: `test -f design_map.md && echo present || echo absent`
If present, add rows mapping `src/ldsc/_kernel/regions.py` and the new config
fields/CLI flags to this design doc.

- [ ] **Step 3: Verify docstrings are complete**

Confirm every new public-ish function (`load_preset_intervals`,
`load_bed_intervals`, `merge_intervals`, `region_exclusion_keep_mask`,
`RefPanel._apply_region_exclusion`) has a NumPy-style docstring (they do per the
code above). Use `my-skills:fun-doc` for any that are thin.

- [ ] **Step 4: Commit**

```bash
git add docs/ design_map.md
git commit -m "docs(regions): document region-exclusion feature and chokepoints"
```

---

## Final verification checklist

- [ ] `pytest -q` is fully green.
- [ ] `ldsc build-ref-panel --help` shows `--exclude-regions` and
      `--exclude-regions-bed` (and **not** `--exclude-regions-build`).
- [ ] `ldsc ldscore --help` shows all three flags.
- [ ] A preset on `ldscore` without `--exclude-regions-build` errors clearly.
- [ ] An in-MHC SNP is absent from both emitted builds of a dual-build
      `build-ref-panel` run.
- [ ] All commits use Conventional Commits and contain no AI-tool authorship.

---

*End of plan.*
