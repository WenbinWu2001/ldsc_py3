"""Render the partitioned-LDSC workflow figure with code-mapped labels.

This script keeps the figure layout deterministic so text-fit fixes can be
made in one place instead of hand-editing the generated SVG.
"""

from __future__ import annotations

import os
from pathlib import Path


# Keep matplotlib and fontconfig from trying to write into the home directory.
_cache_root = Path("/tmp/codex-mpl-cache")
_cache_root.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "mpl"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root / "xdg"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parent
SVG_PATH = ROOT / "partitioned-ldsc-workflow-code-mapped.svg"
PNG_PATH = ROOT / "partitioned-ldsc-workflow-code-mapped.png"

WIDTH = 1600
HEIGHT = 1900

BG = "#f7f4eb"
INK = "#16202a"
MUTED = "#5b6470"
BLUE = "#2563eb"
ORANGE = "#ea580c"
RED = "#dc2626"
MAGENTA = "#d946ef"
GREEN = "#16a34a"
SLATE = "#475569"
PAPER = "#fffdf8"


def lane(ax, x: int, y: int, w: int, h: int, color: str) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=28",
        linewidth=1.4,
        edgecolor=color,
        facecolor=color,
        alpha=0.06,
        zorder=0,
    )
    ax.add_patch(patch)


def box(
    ax,
    x: int,
    y: int,
    w: int,
    h: int,
    title: str,
    color: str,
    body: str,
    *,
    dashed: bool = False,
    fill: str = PAPER,
    title_fs: float = 15.0,
    body_fs: float = 11.5,
) -> None:
    outer = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=18",
        linewidth=2.0,
        edgecolor=INK,
        facecolor=fill,
        linestyle="--" if dashed else "-",
        zorder=2,
    )
    ax.add_patch(outer)

    band_h = 38
    band = Rectangle((x, y + h - band_h), w, band_h, linewidth=0, facecolor=color, alpha=0.18, zorder=2.1)
    ax.add_patch(band)

    ax.text(
        x + 14,
        y + h - 10,
        title,
        va="top",
        ha="left",
        fontsize=title_fs,
        fontweight="bold",
        color=color,
        family="DejaVu Sans",
        linespacing=1.05,
        zorder=3,
    )
    ax.text(
        x + 14,
        y + h - band_h - 12,
        body,
        va="top",
        ha="left",
        fontsize=body_fs,
        color=INK,
        family="DejaVu Sans",
        linespacing=1.22,
        zorder=3,
    )


def tag(ax, x: int, y: int, w: int, h: int, text: str, color: str) -> None:
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=14",
        linewidth=1.2,
        edgecolor=color,
        facecolor=color,
        alpha=0.12,
        zorder=1,
    )
    ax.add_patch(patch)
    ax.text(
        x + 10,
        y + h / 2,
        text,
        va="center",
        ha="left",
        fontsize=11.0,
        color=INK,
        family="DejaVu Sans",
        zorder=2,
    )


def arrow(
    ax,
    x1: int,
    y1: int,
    x2: int,
    y2: int,
    label: str | None = None,
    *,
    color: str = INK,
    rad: float = 0.0,
    label_dx: int = 0,
    label_dy: int = 0,
    label_fs: float = 9.8,
) -> None:
    arr = FancyArrowPatch(
        (x1, y1),
        (x2, y2),
        arrowstyle="-|>",
        mutation_scale=18,
        linewidth=2.0,
        color=color,
        connectionstyle=f"arc3,rad={rad}",
        zorder=1.5,
    )
    ax.add_patch(arr)
    if label:
        mx = (x1 + x2) / 2 + label_dx
        my = (y1 + y2) / 2 + label_dy
        ax.text(
            mx,
            my,
            label,
            ha="center",
            va="center",
            fontsize=label_fs,
            color=color,
            family="DejaVu Sans",
            bbox={"boxstyle": "round,pad=0.16", "fc": BG, "ec": "none", "alpha": 0.96},
            zorder=4,
        )


def render() -> None:
    fig = plt.figure(figsize=(16, 19), facecolor=BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, WIDTH)
    ax.set_ylim(0, HEIGHT)
    ax.axis("off")

    ax.text(
        60,
        1840,
        "Partitioned LDSC Workflow: code-mapped flow chart",
        fontsize=28,
        fontweight="bold",
        color=INK,
        family="DejaVu Sans",
    )
    ax.text(
        60,
        1800,
        "Design doc + code cross-reference. Nodes are code objects; arrow labels use owner functions replacing sketch prose.",
        fontsize=14,
        color=MUTED,
        family="DejaVu Sans",
    )

    tag(ax, 60, 1710, 250, 42, "A = raw reference-panel SNPs", BLUE)
    tag(ax, 330, 1710, 330, 42, "A' = panel after `ref_panel_snps_file`", BLUE)
    tag(ax, 680, 1710, 290, 42, "B = `AnnotationBundle` SNP rows", ORANGE)
    tag(ax, 990, 1710, 220, 42, "C = regression SNP mask", GREEN)
    tag(ax, 60, 1658, 320, 42, "`ld_reference_snps = B ∩ A'`", MAGENTA)
    tag(ax, 400, 1658, 380, 42, "`ld_regression_snps = B ∩ A' ∩ C`", MAGENTA)

    lane(ax, 50, 1040, 700, 560, ORANGE)
    lane(ax, 850, 1040, 700, 560, BLUE)
    lane(ax, 420, 690, 760, 280, MAGENTA)
    lane(ax, 420, 400, 760, 230, SLATE)
    lane(ax, 420, 60, 760, 280, GREEN)

    box(
        ax,
        90,
        1310,
        250,
        170,
        "Baseline annotation\ninput",
        RED,
        "`baseline_annot_sources`\nLDSC `.annot[.gz]` shards",
        title_fs=16,
        body_fs=11.8,
    )
    box(
        ax,
        90,
        1080,
        250,
        185,
        "Query annotation\ninput",
        ORANGE,
        "`query_annot_sources`\nor `bed_paths`\n\nBED paths stay in-memory\nin the main `ldscore`\nworkflow.",
        title_fs=16,
        body_fs=11.2,
    )
    box(
        ax,
        390,
        1080,
        300,
        130,
        "Optional explicit\nmaterialization",
        ORANGE,
        "`run_bed_to_annot(output_dir=...)`\n-> `query.<chrom>.annot.gz`",
        dashed=True,
        title_fs=15,
        body_fs=11.0,
    )
    box(
        ax,
        390,
        1270,
        300,
        250,
        "Annotation bundle\n(`B`)",
        ORANGE,
        "`AnnotationBuilder.run()` parses `.annot`\nwith `parse_annotation_file()` and\nprojects BED with\n`_compute_bed_query_columns()`.\n\nBuilds `AnnotationBundle` with\n`metadata`, `baseline_annotations`,\n`query_annotations`.",
        title_fs=15,
        body_fs=11.2,
    )

    box(
        ax,
        900,
        1310,
        260,
        190,
        "Reference-panel\nsource",
        BLUE,
        "`--bfile` -> PLINK `.bed/.bim/.fam`\n\nor\n\n`--r2-table` -> parquet R2\n+ metadata sidecars",
        title_fs=16,
        body_fs=11.2,
    )
    box(
        ax,
        1200,
        1450,
        280,
        100,
        "Reference-panel\nspec",
        BLUE,
        "`RefPanelConfig`\n`ref_panel_snps_file`",
        title_fs=15,
        body_fs=11.8,
    )
    box(
        ax,
        1200,
        1265,
        280,
        150,
        "Reference-panel\nadapter",
        BLUE,
        "`_ref_panel_from_args()`\n-> `RefPanelLoader.load()`\n-> `PlinkRefPanel` or\n   `ParquetR2RefPanel`",
        title_fs=15,
        body_fs=11.2,
    )
    box(
        ax,
        900,
        1070,
        340,
        200,
        "Restricted metadata\n(`A'`)",
        BLUE,
        "`load_metadata(chrom)`\n+ `RefPanel._apply_snp_restriction()`\n\nReturned metadata already reflects\n`ref_panel_snps_file`.",
        title_fs=15,
        body_fs=11.4,
    )
    box(
        ax,
        1280,
        1070,
        250,
        160,
        "Regression SNP mask\n(`C`)",
        GREEN,
        "`LDScoreConfig.regression_snps_file`\nloaded by `_load_regression_snps()`",
        title_fs=15,
        body_fs=11.2,
    )

    box(
        ax,
        470,
        730,
        660,
        210,
        "Per-chromosome LD computation",
        MAGENTA,
        "`LDScoreCalculator.run()` iterates chromosomes.\n`compute_chromosome()` dispatches to\n`kernel_ldscore.compute_chrom_from_plink()` or\n`kernel_ldscore.compute_chrom_from_parquet()`.\n\nInternal compute-time universe:\n`ld_reference_snps = B ∩ A'`",
        title_fs=15,
        body_fs=11.4,
    )

    box(
        ax,
        470,
        430,
        320,
        175,
        "Aggregated normalized\nresult",
        MAGENTA,
        "`LDScoreCalculator._aggregate_chromosome_results()`\n-> `LDScoreResult`\n\nRows = `ld_regression_snps = B ∩ A' ∩ C`\nColumns = `CHR`, `SNP`, `BP`, baseline L2,\nquery L2, `regr_weight`.\n\nPublic normalized results keep\n`ld_reference_snps = frozenset()`.",
        title_fs=15,
        body_fs=10.8,
    )
    box(
        ax,
        820,
        430,
        310,
        175,
        "Written artifacts",
        SLATE,
        "`LDScoreDirectoryWriter.write()`\n\n`baseline.parquet`, optional `query.parquet`,\n`manifest.json`\n\nManifest counts:\n`all_reference_snp_count`\noptional `common_reference_snp_count`\nplus top-level `count_config`",
        title_fs=15,
        body_fs=10.8,
    )

    box(
        ax,
        470,
        110,
        210,
        170,
        "Summary-stat\ninput",
        GREEN,
        "Curated `sumstats.parquet`\n`load_sumstats()` -> `SumstatsTable`",
        title_fs=15,
        body_fs=11.6,
    )
    box(
        ax,
        720,
        110,
        210,
        170,
        "Reload path",
        SLATE,
        "`load_ldscore_from_files()`\nrebuilds normalized\n`LDScoreResult` from disk",
        title_fs=15,
        body_fs=11.6,
    )
    box(
        ax,
        970,
        110,
        250,
        170,
        "Merged regression\ndataset",
        GREEN,
        "`RegressionRunner.build_dataset()`\nmerges sumstats + `ldscore_table`\n+ embedded `regr_weight`",
        title_fs=15,
        body_fs=11.2,
    )
    box(
        ax,
        1260,
        110,
        240,
        190,
        "Partitioned h²\noutput",
        GREEN,
        "`RegressionRunner.estimate_partitioned_h2_batch()`\nloops over `query_columns` and is used by\n`run_partitioned_h2_from_args()`.\n\nWrites `<out>.partitioned_h2.tsv`.",
        title_fs=15,
        body_fs=11.0,
    )

    arrow(ax, 340, 1395, 390, 1395, "input rows", color=RED, label_dy=14)
    arrow(ax, 340, 1170, 390, 1170, "`AnnotationBuilder.run()`", color=ORANGE, label_dy=14)
    arrow(ax, 215, 1265, 215, 1310, color=ORANGE)
    arrow(ax, 340, 1160, 390, 1125, "`run_bed_to_annot()`", color=ORANGE, rad=-0.05, label_dy=-12)

    arrow(ax, 1160, 1405, 1200, 1485, "`_ref_panel_from_args()`", color=BLUE, label_dy=16)
    arrow(ax, 1340, 1450, 1340, 1415, color=BLUE)
    arrow(ax, 1200, 1335, 1240, 1335, "`load()`", color=BLUE, label_dy=16)
    arrow(ax, 1160, 1170, 1240, 1170, "`load_metadata()`", color=BLUE, label_dy=16)
    arrow(ax, 1280, 1450, 1070, 1270, "`RefPanel._apply_snp_restriction()`", color=BLUE, rad=0.05, label_dx=-20, label_dy=18)

    arrow(ax, 690, 1395, 790, 835, "`AnnotationBundle` = `B`", color=ORANGE, label_dx=-80, label_dy=26)
    arrow(ax, 1070, 1070, 930, 940, "A'", color=BLUE, label_dx=-20, label_dy=10)
    arrow(ax, 1405, 1070, 1040, 940, "`_load_regression_snps()`", color=GREEN, label_dx=30, label_dy=24)

    arrow(ax, 800, 730, 770, 605, "`_wrap_legacy_chrom_result()`", color=MAGENTA, label_dx=120, label_dy=10)
    arrow(ax, 790, 520, 820, 520, "`write_outputs()`", color=SLATE, label_dy=14)
    arrow(ax, 975, 430, 825, 280, "`load_ldscore_from_files()`", color=SLATE, label_dx=40, label_dy=18)
    arrow(ax, 575, 280, 1095, 280, "`load_sumstats()`", color=GREEN, label_dy=18)
    arrow(ax, 825, 110, 1095, 280, color=SLATE)
    arrow(ax, 1220, 195, 1260, 195, "`estimate_partitioned_h2_batch()`", color=GREEN, label_dy=18)

    note = FancyBboxPatch(
        (60, 1518),
        1490,
        108,
        boxstyle="round,pad=0.02,rounding_size=18",
        linewidth=1.2,
        edgecolor=INK,
        facecolor="#fff8e8",
        zorder=1,
    )
    ax.add_patch(note)
    ax.text(
        80,
        1594,
        "Notes:\n"
        "1. Main `ldscore` uses `AnnotationBuilder.run()` for BED input; `run_bed_to_annot()` is the optional explicit materialization path.\n"
        "2. `ref_panel_snps_file` defines `A'`; `regression_snps_file` defines `C`.\n"
        "3. Manifest counts use `ld_reference_snps`; parquet LD-score rows use `ld_regression_snps`.",
        fontsize=10.2,
        color=INK,
        family="DejaVu Sans",
        linespacing=1.22,
    )

    fig.savefig(SVG_PATH, facecolor=fig.get_facecolor())
    fig.savefig(PNG_PATH, dpi=220, facecolor=fig.get_facecolor())


if __name__ == "__main__":
    render()
