"""Count the SNPs each reference-preparation step removes, per chromosome.

A reference-SNP funnel follows one chromosome through ``RefPanel.prepare_chromosome``
in the backend's own execution order. The ``reference_panel`` funnel starts from
the panel's rows on the chromosome; the ``annotation_snp_grid`` funnel starts
from the annotation rows. Both end at the reference SNPs that enter LD scores,
SNP counts and overlap. The annotation files define the annotation SNP grid, so
panel SNPs that are not on the grid are excluded at the ``annotation_snp_grid``
step; that exclusion is recorded here and reported, never silent
(``docs/reference/ldscore-snp-universes.md``).

Every backend writes the same rows. A step that did not run, because its option
was not set, the identifier mode skips it, or preparation stopped earlier, has
``removed`` missing and leaves ``remaining`` unchanged.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .._core.chromosomes import chrom_sort_key

REFERENCE_PANEL = "reference_panel"
ANNOTATION_SNP_GRID = "annotation_snp_grid"
FUNNEL_COLUMNS = ["funnel", "chrom", "backend", "step_index", "step", "removed", "remaining"]

_STEPS = {
    ("plink", REFERENCE_PANEL): (
        "panel_rows", "identity_cleanup", "ref_panel_snps_file", "allele_aware_identity_cleanup",
        ANNOTATION_SNP_GRID, "genotype_qc", "maf_floor",
    ),
    ("parquet_r2", REFERENCE_PANEL): (
        "panel_rows", "ref_panel_snps_file", "maf_floor", "allele_aware_identity_cleanup", ANNOTATION_SNP_GRID,
    ),
    ("plink", ANNOTATION_SNP_GRID): ("grid_rows", ANNOTATION_SNP_GRID, "genotype_qc", "maf_floor"),
    ("parquet_r2", ANNOTATION_SNP_GRID): ("grid_rows", ANNOTATION_SNP_GRID),
}


class ReferenceSNPFunnel:
    """Record one chromosome's funnel while its reference preparation runs.

    ``record`` takes the starting row count for a funnel's first step and the
    number of SNPs removed for every later step. Steps never recorded did not run.
    """

    def __init__(self, chrom: str, backend: str) -> None:
        self.chrom = str(chrom)
        self.backend = backend
        self._counts: dict[tuple[str, str], int] = {}

    def record(self, funnel: str, step: str, count: int) -> None:
        """Record the start count (first step) or the removed count (later steps)."""
        self._counts[(funnel, step)] = int(count)

    def frame(self) -> pd.DataFrame:
        """Return the fixed rows of both funnels with running ``remaining`` counts."""
        rows = []
        for funnel in (REFERENCE_PANEL, ANNOTATION_SNP_GRID):
            steps = _STEPS[(self.backend, funnel)]
            remaining = self._counts.get((funnel, steps[0]))
            rows.append((funnel, self.chrom, self.backend, 0, steps[0], pd.NA,
                         pd.NA if remaining is None else remaining))
            for index, step in enumerate(steps[1:], 1):
                removed = self._counts.get((funnel, step))
                if removed is not None and remaining is not None:
                    remaining -= removed
                rows.append((funnel, self.chrom, self.backend, index, step,
                             pd.NA if removed is None else removed, pd.NA if remaining is None else remaining))
        return _typed(pd.DataFrame(rows, columns=FUNNEL_COLUMNS))


def _typed(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.astype({"funnel": str, "chrom": str, "backend": str, "step_index": "int64",
                         "step": str, "removed": "Int64", "remaining": "Int64"})


def concat_funnels(frames) -> pd.DataFrame | None:
    """Combine chromosome funnels in genomic, funnel and step order; ``None`` if there are none."""
    frames = [frame for frame in frames if frame is not None]
    if not frames:
        return None
    combined = pd.concat(frames, ignore_index=True)
    order = {REFERENCE_PANEL: 0, ANNOTATION_SNP_GRID: 1}
    combined = combined.assign(_chrom=combined["chrom"].map(chrom_sort_key), _funnel=combined["funnel"].map(order))
    combined = combined.sort_values(["_chrom", "_funnel", "step_index"], kind="stable")
    return _typed(combined.drop(columns=["_chrom", "_funnel"]).reset_index(drop=True))


def _grid_step(frame: pd.DataFrame, funnel: str) -> pd.DataFrame:
    return frame.loc[frame["funnel"].eq(funnel) & frame["step"].eq(ANNOTATION_SNP_GRID)]


def annotation_snp_grid_emptied(frame: pd.DataFrame) -> bool:
    """Return whether the grid step removed every panel SNP still left on the chromosome.

    A chromosome emptied by an earlier step, such as ``--ref-panel-snps-file``,
    reaches the grid step with no SNPs and is not emptied by it.
    """
    step = _grid_step(frame, REFERENCE_PANEL)
    removed, remaining = step["removed"].fillna(0), step["remaining"].fillna(0)
    return bool(((removed > 0) & (remaining == 0)).any())


def excluded_reference_snp_counts(frame: pd.DataFrame) -> pd.Series:
    """Return panel SNPs removed at the ``annotation_snp_grid`` step, by chromosome, when nonzero."""
    step = _grid_step(frame, REFERENCE_PANEL)
    counts = step.set_index("chrom")["removed"].fillna(0).astype("int64")
    return counts.loc[counts > 0]


def with_grid_rows(frame: pd.DataFrame, grid_rows: int) -> pd.DataFrame:
    """Restart the ``annotation_snp_grid`` funnel at ``grid_rows`` rows, before an earlier intersection.

    Rows dropped between ``grid_rows`` and the rows that reached reference
    preparation are counted at the grid step; later ``remaining`` values are
    unchanged.
    """
    frame = frame.copy()
    grid = frame["funnel"].eq(ANNOTATION_SNP_GRID)
    start = grid & frame["step_index"].eq(0)
    extra = int(grid_rows) - int(frame.loc[start, "remaining"].iloc[0])
    frame.loc[start, "remaining"] = int(grid_rows)
    step = grid & frame["step"].eq(ANNOTATION_SNP_GRID)
    frame.loc[step, "removed"] = frame.loc[step, "removed"].fillna(0) + extra
    return _typed(frame)


def reference_snp_funnel_warning(frame: pd.DataFrame | None, *, command: str, diagnostics: str) -> str | None:
    """Return the run-level warning when the grid step excluded any panel SNP, else ``None``."""
    if frame is None:
        return None
    counts = excluded_reference_snp_counts(frame)
    if counts.empty:
        return None
    total, chromosomes = int(counts.sum()), len(counts)
    return (
        f"{command} excluded {total:,} reference-panel SNPs on {chromosomes} "
        f"chromosome{'s' if chromosomes != 1 else ''} that are not on the annotation SNP grid; "
        "LD scores, SNP counts and output rows use only SNPs on the grid. This is expected for "
        "legacy baseline-LD annotations on a newer reference panel; otherwise check for a genome-build, "
        f"allele or SNP-identifier mismatch. See {diagnostics}."
    )


def write_reference_snp_funnel(frame: pd.DataFrame, path: Path) -> None:
    """Write the funnel table as TSV, with missing counts as ``NA``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def read_reference_snp_funnel(path: str | Path) -> pd.DataFrame:
    """Read a funnel table written by ``write_reference_snp_funnel``."""
    frame = pd.read_csv(path, sep="\t", dtype={"chrom": str}, keep_default_na=False, na_values=["NA"])
    return _typed(frame)


def funnel_summary(frame: pd.DataFrame) -> str:
    """Render one chromosome's reference-panel funnel as ``step=removed`` pairs for logs."""
    rows = frame.loc[frame["funnel"].eq(REFERENCE_PANEL)].sort_values("step_index")
    parts = [f"{row.step}={'NA' if pd.isna(row.removed) else int(row.removed)}" for row in rows.itertuples()]
    parts[0] = f"panel_rows={rows['remaining'].iloc[0]}"
    return ", ".join(parts) + f"; reference SNPs={rows['remaining'].iloc[-1]}"


def funnel_records(frame: pd.DataFrame) -> list[dict]:
    """Return JSON-ready funnel rows, with missing counts as ``None``."""
    integers = {"step_index", "removed", "remaining"}
    return [{column: (None if pd.isna(value) else int(value) if column in integers else str(value))
             for column, value in row.items()} for row in frame.to_dict("records")]


def funnel_from_records(records) -> pd.DataFrame:
    """Rebuild a funnel table from ``funnel_records`` output."""
    return _typed(pd.DataFrame(list(records), columns=FUNNEL_COLUMNS))
