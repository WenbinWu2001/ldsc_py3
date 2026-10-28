"""Render the current workflow from direct LD scores to partitioned regression.

The SVG and PNG are generated together from the same code-mapped labels.
Indexed scoring shares the publication and reader boundary shown here; its
separate preparation and atom assembly are documented in
docs/reference/gene-ldscore-index.md. Both files are written to docs/assets/.
"""

from __future__ import annotations

import os
from pathlib import Path
import tempfile


_cache_root = Path(tempfile.gettempdir()) / "ldsc-workflow-figure-cache"
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "mpl"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root / "xdg"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[2] / "docs" / "assets"
INK = "#172c3c"
MUTED = "#496273"
BLUE = "#2d658f"
GREEN = "#2b7563"
ORANGE = "#9a602b"


def box(ax, x, y, width, height, title, lines, color):
    """Draw one responsibility with fixed, readable line spacing."""
    ax.add_patch(FancyBboxPatch(
        (x, y), width, height, boxstyle="round,pad=0,rounding_size=10",
        linewidth=1.4, edgecolor=color, facecolor="white", zorder=2,
    ))
    ax.text(x + 20, y + height - 20, title, va="top", color=color,
            fontsize=15, fontweight="bold", zorder=3)
    ax.text(x + 20, y + height - 55, "\n".join(lines), va="top", color=INK,
            fontsize=11.5, linespacing=1.65, zorder=3)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=14, linewidth=1.5,
        color=MUTED, shrinkA=4, shrinkB=4, zorder=1,
    ))


def render():
    """Write the checked-in SVG and PNG with the same deterministic layout."""
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "svg.fonttype": "none",
        "svg.hashsalt": "ldsc-partitioned-workflow",
    })
    fig, ax = plt.subplots(figsize=(12, 16.08), dpi=160)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.set_facecolor("#f4f7f8")
    ax.set(xlim=(0, 1000), ylim=(0, 1340))
    ax.axis("off")
    ax.text(50, 1300, "Direct LD scoring to partitioned regression", color=INK,
            fontsize=22, fontweight="bold", va="top")
    ax.text(50, 1260, "Current module boundaries and canonical artifacts  |  2026-10-15",
            fontsize=12, color=MUTED, va="top")

    box(ax, 50, 1060, 430, 160, "Annotation inputs", [
        "Baseline .annot(.gz) + optional query tables",
        "or BED / resolved gene intervals",
        "AnnotationBuilder -> chromosome shards",
        "iter_bed_rows -> ChromosomeProjector",
    ], BLUE)
    box(ax, 520, 1060, 430, 160, "Reference inputs", [
        "PLINK: --plink-prefix",
        "or R$^2$ matrix: --r2-matrix-dir",
        "Resolved paths and chromosome scope",
        "Shared identity and genome-build policy",
    ], BLUE)
    box(ax, 50, 870, 900, 145, "1. Prepare one chromosome", [
        "RefPanel.prepare_chromosome() -> PreparedChromosome",
        "Align annotation/reference rows; apply reference QC, restrictions and LD windows.",
        "Own the PLINK/R$^2$ reader and preserve physical PLINK row mapping.",
    ], GREEN)
    box(ax, 50, 680, 900, 145, "2. Compute direct LD scores", [
        "LDScoreCalculator.compute_chromosome() -> _ldscore.compute.compute_chromosome()",
        "Accumulate LD scores, annotation counts and overlaps in float64.",
        "Workers stage float32 query fragments; the parent combines chromosome results.",
    ], GREEN)
    box(ax, 50, 470, 900, 165, "3. Publish the LD-score directory", [
        "LDScoreDirectoryWriter.write_batches() (also used by indexed scoring)",
        "metadata.json + ldscore.baseline.parquet + ldscore.overlap.parquet when required",
        "Optional query files: ldscore.query.parquet or ldscore.query.batchNNNNN.parquet",
        "Ordered query manifest in metadata.json; chromosome row groups in baseline/query tables.",
    ], ORANGE)
    box(ax, 50, 280, 430, 145, "4. Open saved LD scores", [
        "load_ldscore_from_dir() -> LDScoreSource",
        "Owner: _ldscore/source.py",
        "Shared baseline, identities, counts, overlap",
        "read_queries(names): uncached selections",
    ], BLUE)
    box(ax, 520, 280, 430, 145, "Munge or load GWAS sumstats", [
        "SumstatsMunger / load_sumstats()",
        "Munged Parquet or legacy .sumstats(.gz)",
        "Canonical identities, Z, N and allele metadata",
        "Validated SumstatsTable for trait alignment",
    ], BLUE)
    box(ax, 50, 85, 900, 150, "5. Fit and write partitioned heritability", [
        "RegressionRunner -> _regression.estimators.Hsq",
        "Fit the complete baseline model, or one baseline-plus-query model per query.",
        "PartitionedH2DirectoryWriter -> partitioned_h2.tsv + diagnostics/",
    ], GREEN)
    for start, end in [
        ((265, 1060), (265, 1015)), ((735, 1060), (735, 1015)),
        ((500, 870), (500, 825)), ((500, 680), (500, 635)),
        ((265, 470), (265, 425)),
        ((265, 280), (265, 235)), ((735, 280), (735, 235)),
    ]:
        arrow(ax, start, end)

    ax.text(50, 48, "Counts and overlap use retained reference SNPs; regression SNPs and weights use a separate filtered set.",
            fontsize=10.2, color=MUTED, va="top")
    ax.text(50, 25, "Public annotation reads and saved LD scores are float32. Numerical accumulation remains float64.",
            fontsize=10.2, color=MUTED, va="top")
    svg_path = ROOT / "partitioned-ldsc-workflow-code-mapped.svg"
    fig.savefig(svg_path, metadata={"Date": None})
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
    fig.savefig(ROOT / "partitioned-ldsc-workflow-code-mapped.png", metadata={"Software": "LDSC"})
    plt.close(fig)


if __name__ == "__main__":
    render()
