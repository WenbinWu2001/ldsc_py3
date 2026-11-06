"""Coordinate liftover mechanics shared by summary-statistics munging and reference-panel building.

The helpers in this module are deliberately package-private mechanics rather
than a public liftover API. This module owns chain-file translation of 1-based
positions, duplicate-coordinate detection, readable drop reports and their
DEBUG-only examples. Each workflow owns its own liftover policy, drop-row
shaping and messages: summary-statistics liftover (HM3 quick liftover and the
munging drop audit) lives in ``ldsc._sumstats.liftover``, and
reference-panel liftover in ``ldsc._refpanel.builder``.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from os import PathLike
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._core.chromosomes import normalize_chromosome
from ..errors import LDSCDependencyError, LDSCInputError, LDSCInternalError, LDSCUsageError


@dataclass(frozen=True)
class LiftOverMappingResult:
    """Partial 1-based liftover result for one chromosome.

    ``translated_positions`` contains only retained same-chromosome mappings,
    while ``keep_mask`` and the two reason masks are aligned to the input
    position sequence. Chain hits keep the first same-chromosome candidate and
    treat all-cross-chromosome hits as drops.
    """

    translated_positions: np.ndarray
    keep_mask: np.ndarray
    unmapped_count: int
    cross_chrom_count: int
    unmapped_mask: np.ndarray
    cross_chrom_mask: np.ndarray


@dataclass(frozen=True)
class LiftoverDropReport:
    """Readable drop report shared by liftover-aware workflows."""

    reason: str
    n_dropped: int
    examples: list[dict[str, Any]]


@dataclass(frozen=True)
class DuplicateCoordinateDropResult:
    """Result of dropping all rows in duplicate coordinate groups."""

    keep_mask: np.ndarray
    report: LiftoverDropReport


class LiftOverTranslator:
    """Translate one chromosome's 1-based positions between explicit genome-build chain files."""

    def __init__(
        self,
        *,
        source_build: str,
        target_build: str,
        chain_path: str | PathLike[str] | None,
        chain_flag_hint: str,
        workflow_label: str,
    ) -> None:
        """Initialize a build-to-build position translator for one run."""
        self.source_build = source_build
        self.target_build = target_build
        self.chain_path = None if chain_path is None else Path(chain_path)
        if source_build == target_build:
            self._identity = True
            self._liftover = None
            return
        self._identity = False
        if self.chain_path is None:
            raise LDSCUsageError(
                f"Cannot run {workflow_label}: liftover from {source_build} to {target_build} requires a chain file. "
                "Most likely a build-changing liftover was requested without the matching chain path. "
                f"Provide the chain file with {chain_flag_hint}."
            )
        try:
            from pyliftover import LiftOver
        except ImportError as exc:
            raise LDSCDependencyError(
                f"{workflow_label} from {source_build} to {target_build} could not run chain-file liftover "
                "because pyliftover is not installed. Most likely a build-changing liftover was requested "
                "in an environment missing the optional dependency. Install pyliftover or use a workflow "
                "that does not request chain-file liftover."
            ) from exc
        self._liftover = LiftOver(str(self.chain_path))

    def map_positions(self, chrom: str, positions: Sequence[int] | np.ndarray) -> LiftOverMappingResult:
        """Translate 1-based positions and retain only same-chromosome mappings."""
        array = np.asarray(positions, dtype=np.int64)
        if self._identity:
            return LiftOverMappingResult(
                translated_positions=array.copy(),
                keep_mask=np.ones(len(array), dtype=bool),
                unmapped_count=0,
                cross_chrom_count=0,
                unmapped_mask=np.zeros(len(array), dtype=bool),
                cross_chrom_mask=np.zeros(len(array), dtype=bool),
            )
        chrom = _normalize_liftover_chromosome(chrom)
        query_chrom = f"chr{chrom}"
        translated = np.zeros(len(array), dtype=np.int64)
        keep_mask = np.zeros(len(array), dtype=bool)
        unmapped_mask = np.zeros(len(array), dtype=bool)
        cross_chrom_mask = np.zeros(len(array), dtype=bool)
        unmapped_count = 0
        cross_chrom_count = 0
        for idx, pos in enumerate(array):
            hits = self._liftover.convert_coordinate(query_chrom, int(pos) - 1)
            if not hits:
                unmapped_count += 1
                unmapped_mask[idx] = True
                continue
            hit = next(
                (candidate for candidate in hits if _try_normalize_liftover_chromosome(candidate[0]) == chrom),
                None,
            )
            if hit is None:
                cross_chrom_count += 1
                cross_chrom_mask[idx] = True
                continue
            translated[idx] = int(hit[1]) + 1
            keep_mask[idx] = True
        return LiftOverMappingResult(
            translated_positions=translated[keep_mask],
            keep_mask=keep_mask,
            unmapped_count=unmapped_count,
            cross_chrom_count=cross_chrom_count,
            unmapped_mask=unmapped_mask,
            cross_chrom_mask=cross_chrom_mask,
        )


def duplicate_coordinate_drop_result(
    frame: pd.DataFrame,
    *,
    chrom_col: str = "CHR",
    pos_col: str = "POS",
    source_pos_col: str | None = None,
    target_pos_col: str | None = None,
    snp_col: str = "SNP",
    reason: str,
    max_examples: int = 5,
) -> DuplicateCoordinateDropResult:
    """Return rows kept after dropping every row in duplicate coordinate groups.

    The policy is intentionally "drop all duplicates" only; callers decide
    whether to use this result, raise instead, or ignore coordinate duplicates
    for workflows such as rsID source-only reference-panel builds.
    """
    if chrom_col not in frame.columns or pos_col not in frame.columns:
        missing = [column for column in (chrom_col, pos_col) if column not in frame.columns]
        raise LDSCInternalError(
            "liftover duplicate-coordinate detection could not run because the prepared table is missing "
            f"required columns: {missing}. Most likely an upstream liftover caller passed a frame before "
            "coordinate normalization. Re-run with DEBUG logging and report the traceback."
        )
    duplicate_mask = frame.duplicated(subset=[chrom_col, pos_col], keep=False).to_numpy(dtype=bool)
    report = liftover_drop_report(
        frame,
        duplicate_mask,
        reason=reason,
        chrom_col=chrom_col,
        source_pos_col=source_pos_col or pos_col,
        target_pos_col=target_pos_col,
        snp_col=snp_col,
        max_examples=max_examples,
    )
    return DuplicateCoordinateDropResult(keep_mask=~duplicate_mask, report=report)


def liftover_drop_report(
    frame: pd.DataFrame,
    drop_mask: Sequence[bool] | np.ndarray,
    *,
    reason: str,
    chrom_col: str = "CHR",
    source_pos_col: str | None = "POS",
    target_pos_col: str | None = None,
    snp_col: str = "SNP",
    max_examples: int = 5,
) -> LiftoverDropReport:
    """Build a compact readable drop report with shared core fields."""
    mask = np.asarray(drop_mask, dtype=bool)
    if len(mask) != len(frame):
        raise LDSCInternalError(
            "liftover drop-report construction failed because drop_mask length does not match the input frame. "
            "Most likely an upstream liftover filter produced a mask from a different table version. "
            "Re-run with DEBUG logging and report the traceback."
        )
    dropped = frame.loc[mask]
    examples = _drop_examples(
        dropped,
        reason=reason,
        chrom_col=chrom_col,
        source_pos_col=source_pos_col,
        target_pos_col=target_pos_col,
        snp_col=snp_col,
        max_examples=max_examples,
    )
    return LiftoverDropReport(reason=reason, n_dropped=int(mask.sum()), examples=examples)


def log_liftover_drop_report(
    logger: logging.Logger,
    report: LiftoverDropReport,
    *,
    workflow_label: str,
    sidecar_path: str | PathLike[str] | None = None,
) -> None:
    """Log one shared liftover drop report with count-only normal verbosity."""
    if report.n_dropped == 0:
        return
    sidecar_suffix = (
        f" See '{sidecar_path}' for row-level provenance."
        if sidecar_path is not None
        else ""
    )
    logger.info(
        f"{workflow_label} dropped {report.n_dropped} SNPs for "
        f"{report.reason}.{sidecar_suffix}"
    )
    if report.examples:
        logger.debug(f"{workflow_label} example dropped ({report.reason}): {report.examples}")


def _drop_examples(
    dropped: pd.DataFrame,
    *,
    reason: str,
    chrom_col: str,
    source_pos_col: str | None,
    target_pos_col: str | None,
    snp_col: str,
    max_examples: int,
) -> list[dict[str, Any]]:
    """Format up to ``max_examples`` dropped rows using the shared core fields."""
    examples: list[dict[str, Any]] = []
    for _, row in dropped.head(max_examples).iterrows():
        example: dict[str, Any] = {}
        if snp_col in dropped.columns:
            example["SNP"] = _example_value(row[snp_col])
        if chrom_col in dropped.columns:
            example["CHR"] = _example_value(row[chrom_col])
        if source_pos_col is not None and source_pos_col in dropped.columns:
            example["source_POS"] = _example_value(row[source_pos_col])
        if target_pos_col is not None and target_pos_col in dropped.columns:
            example["target_POS"] = _example_value(row[target_pos_col])
        example["reason"] = reason
        examples.append(example)
    return examples


def _example_value(value: object) -> object:
    """Return a JSON-like scalar for human-readable example logging."""
    if pd.isna(value):
        return None
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _normalize_liftover_chromosome(value: object) -> str:
    """Normalize one liftover chromosome label to the package canon."""
    return normalize_chromosome(value)


def _try_normalize_liftover_chromosome(value: object) -> str | None:
    """Normalize a liftover hit chromosome, ignoring unsupported auxiliary contigs."""
    try:
        return _normalize_liftover_chromosome(value)
    except (ValueError, LDSCInputError):
        return None


