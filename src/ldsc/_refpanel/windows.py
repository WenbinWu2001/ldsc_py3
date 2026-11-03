"""LD-window geometry for LD-score calculation and R2 matrix building.

Window coordinates in SNPs, kilobases or centimorgans, the block-left index
of each SNP's window, and the checks that a window fits the reference data.
``window_coordinates`` and ``get_block_lefts`` serve both ``ldscore`` and
``build-r2-matrix``; each command checks its own CM values first.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import LDWindow
from ..errors import LDSCInputError, LDSCInternalError, LDSCUsageError


def get_block_lefts(coords: np.ndarray, max_dist: float) -> np.ndarray:
    """
    Compute the standard LDSC ``block_left`` window start array.

    ``block_left[i]`` stores the first index retained in the LD window centered
    on SNP ``i`` for the chosen distance metric.
    """
    M = len(coords)
    j = 0
    block_left = np.zeros(M, dtype=np.int64)
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1
        block_left[i] = j
    return block_left


def validate_window_positions_sorted(metadata: pd.DataFrame, chrom: str) -> None:
    """Tripwire: base-pair positions must be non-decreasing within a chromosome.

    :func:`get_block_lefts` is a forward-only two-pointer scan whose ``coords``
    must be non-decreasing, or the LD-window blocks are silently wrong. Every
    caller already supplies coordinate-sorted rows by construction: annotation
    metadata is sorted in :func:`parse_annotation_file`, and the PLINK path
    inherits that order through ``keep_snps`` (the ``.bim`` file order is
    discarded). This guard makes that implicit invariant explicit so a future
    change that breaks the upstream sort fails loudly here instead of producing
    wrong LD scores. Position drops across chromosome boundaries are expected
    (the window is per chromosome) and are not flagged.
    """
    pos = pd.to_numeric(metadata["POS"], errors="raise").to_numpy()
    chrom_vals = metadata["CHR"].to_numpy()
    same_chrom = chrom_vals[1:] == chrom_vals[:-1]
    out_of_order = (np.diff(pos) < 0) & same_chrom
    if out_of_order.any():
        first = int(np.flatnonzero(out_of_order)[0])
        raise LDSCInternalError(
            f"ldscore LD-window construction failed on chromosome {chrom}: reference "
            f"SNP positions are not in non-decreasing order (POS {int(pos[first + 1])} "
            f"follows POS {int(pos[first])}). Callers supply coordinate-sorted rows by "
            "construction, so most likely an upstream sorting step was removed or "
            "bypassed. Re-run with DEBUG logging and report the traceback."
        )


def assert_cm_usable(cm: pd.Series, chrom: str) -> None:
    """Reject a CM column that cannot order SNPs on a chromosome.

    Unusable means fewer than two distinct finite values (all zero, all identical,
    or all missing): every SNP collapses to one coordinate, so any positive cM
    window spans the whole chromosome -- a meaningless result, not an intentional
    choice. This check is unconditional; ``--yes-really`` does not bypass it.
    """
    finite = pd.to_numeric(cm, errors="coerce").dropna().unique()
    if len(finite) < 2:
        raise LDSCInputError(
            f"ldscore cannot use `--ld-wind-cm` on chromosome {chrom}: the reference panel "
            "CM column is all zero or otherwise uninformative (fewer than two distinct "
            "values), so it cannot define a genetic-distance window. Provide real "
            "genetic-map positions: pass a PLINK `.bim` with informative CM, add "
            "`--genetic-map-hg38-sources` / `--genetic-map-hg19-sources` for the panel's "
            "build, or use `--ld-wind-kb` / `--ld-wind-snps`. "
            "See docs/user-guide/troubleshooting.md#ldscore-unusable-cm-for-ld-wind-cm."
        )


def window_coordinates(
    positions: np.ndarray,
    cm: np.ndarray | None,
    window: LDWindow,
    *,
    snp_positions: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    """Build the coordinate array and maximum distance for one LD window.

    A SNP-count window counts retained rows unless ``snp_positions`` gives each
    retained SNP's row in a larger SNP set (the R2 matrix sidecar), in which case
    it counts rows of that set. ``snp_positions`` must be increasing; parquet
    preparation orders retained SNPs by panel row. ``cm`` must be complete for a
    cM window; callers raise their own command-specific error first.
    """
    if window.unit == "snps":
        rows = np.arange(len(positions)) if snp_positions is None else snp_positions
        return np.asarray(rows, dtype=float), float(window.value)
    if window.unit == "kb":
        return np.asarray(positions, dtype=float), float(window.value) * 1000.0
    return np.asarray(cm, dtype=float), float(window.value)


class WholeChromosomeWindowError(LDSCUsageError):
    """An LD window spans a whole chromosome and no override allows it."""


def check_whole_chromosome_window(block_left: np.ndarray, chrom: str, *, allowed: bool) -> None:
    """Reject a window that spans the whole chromosome unless ``allowed``.

    The window spans the chromosome when the last SNP's window still starts at
    the first SNP. The error is worded for ``ldscore`` and its ``--yes-really``
    override; the gene-index build, which has no override, catches
    ``WholeChromosomeWindowError`` and raises its own error.
    """
    if len(block_left) == 0 or allowed:
        return
    if block_left[-1] == 0:
        raise WholeChromosomeWindowError(
            f"ldscore would use a whole-chromosome LD window on chromosome {chrom}. "
            "Most likely the LD-window setting is too large for this chromosome. "
            "Use a smaller LD window, or rerun with `--yes-really` if a whole-chromosome "
            "window is intentional."
        )
