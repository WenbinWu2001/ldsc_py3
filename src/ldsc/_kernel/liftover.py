"""Internal coordinate liftover helpers shared by LDSC workflows.

The helpers in this module are deliberately package-private mechanics rather
than a public liftover API. Workflow modules own user-facing validation and
output contracts; this module owns reusable coordinate normalization, chain
translation, duplicate-coordinate detection, drop accounting, sidecar row
shaping, and DEBUG-only readable drop examples for sumstats munging and
reference-panel building.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import logging
from os import PathLike
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._coordinates import complete_coordinate_mask, normalize_chr_pos_frame, positive_int_position_series
from ..chromosome_inference import normalize_chromosome
from ..column_inference import (
    CHR_COLUMN_SPEC,
    SNP_COLUMN_SPEC,
    ColumnSpec,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_required_column,
)
from ..errors import LDSCDependencyError

LOGGER = logging.getLogger("LDSC.liftover")

SUPPORTED_LIFTOVER_BUILDS = {"hg19", "hg38"}
LIFTOVER_DROP_COLUMNS = ["CHR", "SNP", "source_pos", "target_pos", "reason"]

HM3_HG19_POS_SPEC = ColumnSpec(
    "hg19_POS",
    (
        "hg19_POS",
        "hg19_BP",
        "hg19_POSITION",
        "hg37_POS",
        "hg37_BP",
        "GRCh37_POS",
        "GRCh37_BP",
        "GRCh37_POSITION",
    ),
    "hg19 position",
    allow_suffix_match=False,
)
HM3_HG38_POS_SPEC = ColumnSpec(
    "hg38_POS",
    (
        "hg38_POS",
        "hg38_BP",
        "hg38_POSITION",
        "GRCh38_POS",
        "GRCh38_BP",
        "GRCh38_POSITION",
    ),
    "hg38 position",
    allow_suffix_match=False,
)


@dataclass(frozen=True)
class LiftOverMappingResult:
    """Partial 1-based liftover result for one chromosome.

    ``translated_positions`` contains only retained same-chromosome mappings,
    while ``keep_mask`` and the optional reason masks are aligned to the input
    position sequence. Chain hits keep the first same-chromosome candidate and
    treat all-cross-chromosome hits as drops.
    """

    translated_positions: np.ndarray
    keep_mask: np.ndarray
    unmapped_count: int
    cross_chrom_count: int
    unmapped_mask: np.ndarray | None = None
    cross_chrom_mask: np.ndarray | None = None


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
        chain_path: str | PathLike[str] | None = None,
        chain_flag_hint: str | None = None,
        workflow_label: str = "liftover",
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
            flag = chain_flag_hint or _default_chain_flag_hint(source_build, target_build)
            raise ValueError(
                f"{workflow_label} requires an explicit liftover chain path for {source_build} -> {target_build}. "
                f"Provide it via {flag}."
            )
        try:
            from pyliftover import LiftOver
        except ImportError as exc:
            raise LDSCDependencyError(
                f"{workflow_label} from {source_build} to {target_build} requires the optional dependency "
                "'pyliftover'. Install pyliftover or use a workflow that does not request chain-file liftover."
            ) from exc
        self.chain_path = Path(self.chain_path)
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

    def translate_positions(self, chrom: str, positions: Sequence[int] | np.ndarray) -> np.ndarray:
        """Translate 1-based positions and require complete same-chromosome mapping."""
        result = self.map_positions(chrom, positions)
        if not result.keep_mask.all():
            failed = np.asarray(positions, dtype=np.int64)[~result.keep_mask]
            preview = ", ".join(str(int(value)) for value in failed[:5])
            raise ValueError(
                f"Failed to liftover {len(failed)} positions on chromosome {chrom} "
                f"from {self.source_build} to {self.target_build}: {preview}"
            )
        return result.translated_positions


@dataclass(frozen=True)
class SumstatsLiftoverRequest:
    """Resolved liftover request passed from the munger workflow into the kernel."""

    target_build: str | None = None
    liftover_chain_file: str | PathLike[str] | None = None
    use_hm3_quick_liftover: bool = False
    hm3_map_file: str | PathLike[str] | None = None

    def __post_init__(self) -> None:
        """Normalize build and path fields, then validate method-level shape."""
        target_build = normalize_genome_build(self.target_build)
        if target_build == "auto":
            raise ValueError("target_genome_build must be hg19 or hg38; 'auto' is not a valid liftover target.")
        if target_build is not None and target_build not in SUPPORTED_LIFTOVER_BUILDS:
            raise ValueError("target_genome_build must be hg19 or hg38.")
        chain_file = None if self.liftover_chain_file is None else str(self.liftover_chain_file)
        hm3_map_file = None if self.hm3_map_file is None else str(self.hm3_map_file)
        object.__setattr__(self, "target_build", target_build)
        object.__setattr__(self, "liftover_chain_file", chain_file)
        object.__setattr__(self, "hm3_map_file", hm3_map_file)
        if chain_file is not None and self.use_hm3_quick_liftover:
            raise ValueError("liftover_chain_file and use_hm3_quick_liftover are mutually exclusive.")
        if self.method is not None and target_build is None:
            raise ValueError("target_genome_build is required when a liftover method is specified.")

    @property
    def method(self) -> str | None:
        """Return the selected liftover method token, if any."""
        if self.liftover_chain_file is not None:
            return "chain_file"
        if self.use_hm3_quick_liftover:
            return "hm3_curated"
        return None

    @property
    def requested(self) -> bool:
        """Return whether any target or method flag was supplied."""
        return self.target_build is not None or self.method is not None


@dataclass(frozen=True)
class Hm3DualBuildLifter:
    """Coordinate-only lifter backed by the curated dual-build HM3 map."""

    source_build: str
    target_build: str
    map_path: str | PathLike[str]

    def __post_init__(self) -> None:
        """Normalize build labels and validate the supported HM3 pair."""
        source = normalize_genome_build(self.source_build)
        target = normalize_genome_build(self.target_build)
        if source not in SUPPORTED_LIFTOVER_BUILDS or target not in SUPPORTED_LIFTOVER_BUILDS:
            raise ValueError("HM3 quick liftover supports only hg19 and hg38.")
        object.__setattr__(self, "source_build", source)
        object.__setattr__(self, "target_build", target)
        object.__setattr__(self, "map_path", str(self.map_path))

    def lift(self, frame: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
        """Return a lifted frame and the original indices that failed to map."""
        if self.source_build == self.target_build:
            return frame.copy(), np.asarray([], dtype=np.int64)
        work, _report = _normalized_coordinate_frame(frame)
        request = SumstatsLiftoverRequest(
            target_build=self.target_build,
            use_hm3_quick_liftover=True,
            hm3_map_file=self.map_path,
        )
        lifted, _unmapped_count, _drop_frame = _apply_hm3_liftover(work, request, source_build=self.source_build)
        unmapped_indices = np.asarray(frame.index[~frame.index.isin(lifted.index)], dtype=np.int64)
        return lifted, unmapped_indices


def default_liftover_metadata(*, source_build: str | None, snp_identifier: str) -> dict[str, Any]:
    """Return the no-op liftover report used when liftover is not applicable."""
    normalize_snp_identifier_mode(snp_identifier)
    normalize_genome_build(source_build)
    return {
        "applied": False,
        "source_build": None,
        "target_build": None,
        "method": None,
        "chain_file": None,
        "hm3_map_file": None,
        "n_input": None,
        "n_lifted": None,
        "n_dropped": None,
        "n_missing_chr_pos_dropped": None,
        "n_unmapped": None,
        "n_cross_chrom": None,
        "n_duplicate_source_dropped": None,
        "n_duplicate_target_dropped": None,
    }


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
        raise ValueError(f"Duplicate coordinate detection requires columns: {missing}.")
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
        raise ValueError("drop_mask length must match frame length.")
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


def mapping_reason_masks(result: LiftOverMappingResult) -> tuple[np.ndarray, np.ndarray]:
    """Return unmapped and cross-chromosome masks aligned to a mapping input."""
    dropped_mask = ~np.asarray(result.keep_mask, dtype=bool)
    if result.unmapped_mask is not None:
        unmapped_mask = np.asarray(result.unmapped_mask, dtype=bool)
    elif result.unmapped_count == int(dropped_mask.sum()) and result.cross_chrom_count == 0:
        unmapped_mask = dropped_mask
    else:
        unmapped_mask = np.zeros_like(dropped_mask, dtype=bool)
    if result.cross_chrom_mask is not None:
        cross_chrom_mask = np.asarray(result.cross_chrom_mask, dtype=bool)
    elif result.cross_chrom_count == int(dropped_mask.sum()) and result.unmapped_count == 0:
        cross_chrom_mask = dropped_mask
    else:
        cross_chrom_mask = np.zeros_like(dropped_mask, dtype=bool)
    return unmapped_mask, cross_chrom_mask


def _empty_liftover_drop_frame() -> pd.DataFrame:
    """Return the canonical empty dropped-SNP sidecar frame."""
    return pd.DataFrame(
        {
            "CHR": pd.Series(dtype="string"),
            "SNP": pd.Series(dtype="string"),
            "source_pos": pd.Series(dtype="Int64"),
            "target_pos": pd.Series(dtype="Int64"),
            "reason": pd.Series(dtype="string"),
        },
        columns=LIFTOVER_DROP_COLUMNS,
    )


def _concat_liftover_drop_frames(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate drop frames while preserving the unified nullable schema."""
    non_empty = [frame for frame in frames if frame is not None and not frame.empty]
    if not non_empty:
        return _empty_liftover_drop_frame()
    return _coerce_liftover_drop_frame(pd.concat(non_empty, ignore_index=True))


def _liftover_drop_frame(
    frame: pd.DataFrame,
    drop_mask: Sequence[bool] | np.ndarray,
    *,
    reason: str,
    chrom_col: str = "CHR",
    source_pos_col: str | None = "POS",
    target_pos_col: str | None = None,
    snp_col: str = "SNP",
) -> pd.DataFrame:
    """Build unified dropped-SNP sidecar rows from one filtered frame slice."""
    mask = np.asarray(drop_mask, dtype=bool)
    if len(mask) != len(frame):
        raise ValueError("drop_mask length must match frame length.")
    if not bool(mask.any()):
        return _empty_liftover_drop_frame()
    dropped = frame.loc[mask]
    output = pd.DataFrame(index=dropped.index)
    output["CHR"] = dropped[chrom_col] if chrom_col in dropped.columns else pd.NA
    output["SNP"] = dropped[snp_col] if snp_col in dropped.columns else pd.NA
    if source_pos_col is not None and source_pos_col in dropped.columns:
        output["source_pos"] = dropped[source_pos_col]
    else:
        output["source_pos"] = pd.NA
    if target_pos_col is not None and target_pos_col in dropped.columns:
        output["target_pos"] = dropped[target_pos_col]
    else:
        output["target_pos"] = pd.NA
    output["reason"] = reason
    return _coerce_liftover_drop_frame(output)


def _coerce_liftover_drop_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``frame`` with canonical dropped-SNP sidecar columns and dtypes."""
    output = frame.copy()
    for column in LIFTOVER_DROP_COLUMNS:
        if column not in output.columns:
            output[column] = pd.NA
    output = output.loc[:, LIFTOVER_DROP_COLUMNS]
    output["CHR"] = output["CHR"].astype("string")
    output["SNP"] = output["SNP"].astype("string")
    output["source_pos"] = pd.to_numeric(output["source_pos"], errors="coerce").astype("Int64")
    output["target_pos"] = pd.to_numeric(output["target_pos"], errors="coerce").astype("Int64")
    output["reason"] = output["reason"].astype("string")
    return output.reset_index(drop=True)


def load_hm3_curated_map(path: str | PathLike[str]) -> pd.DataFrame:
    """Load the curated dual-build HM3 map into canonical columns."""
    return _load_hm3_curated_map_cached(str(path)).copy()


@lru_cache(maxsize=8)
def _load_hm3_curated_map_cached(path: str) -> pd.DataFrame:
    """Cached implementation for :func:`load_hm3_curated_map`."""
    df = pd.read_csv(path, sep=None, engine="python")
    context = str(path)
    chr_col = resolve_required_column(df.columns, CHR_COLUMN_SPEC, context=context)
    hg19_col = resolve_required_column(df.columns, HM3_HG19_POS_SPEC, context=context)
    hg38_col = resolve_required_column(df.columns, HM3_HG38_POS_SPEC, context=context)
    snp_col = resolve_required_column(df.columns, SNP_COLUMN_SPEC, context=context)
    out = pd.DataFrame(
        {
            "CHR": [normalize_chromosome(value, context=context) for value in df[chr_col]],
            "hg19_POS": _positive_int_position(df[hg19_col], label="hg19_POS", context=context),
            "hg38_POS": _positive_int_position(df[hg38_col], label="hg38_POS", context=context),
            "SNP": df[snp_col].astype(str),
        }
    )
    _reject_duplicate_hm3_coordinates(out, build="hg19")
    _reject_duplicate_hm3_coordinates(out, build="hg38")
    return out


def apply_sumstats_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest | None,
    *,
    source_build: str | None,
    snp_identifier: str,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Apply a summary-statistics liftover request and return metadata plus drops.

    Liftover is meaningful only in ``chr_pos`` mode. The munger drops missing
    coordinates, source duplicate ``CHR/POS`` groups, unmapped/cross-chromosome
    chain hits, and target duplicate ``CHR/POS`` groups only when a liftover
    request reaches the mapping stage. It updates ``CHR``/``POS`` but preserves
    ``SNP`` as a label field. The third return value is a unified nullable
    dropped-SNP frame with columns ``CHR``, ``SNP``, ``source_pos``,
    ``target_pos``, and ``reason``; clean/no-op runs return the same schema
    with zero rows.
    """
    logger = LOGGER if logger is None else logger
    request = request or SumstatsLiftoverRequest()
    mode = normalize_snp_identifier_mode(snp_identifier)
    source = normalize_genome_build(source_build)
    if request.requested and mode != "chr_pos":
        raise ValueError("Summary-statistics liftover is only valid when snp_identifier='chr_pos'.")
    if not request.requested:
        return frame, default_liftover_metadata(source_build=source, snp_identifier=mode), _empty_liftover_drop_frame()
    if source not in SUPPORTED_LIFTOVER_BUILDS:
        raise ValueError(
            "Cannot apply summary-statistics liftover because the source genome build is unresolved. "
            "Use --genome-build hg19/hg38 or --genome-build auto with inferable coordinates."
        )
    if request.target_build == source:
        if request.method is not None:
            raise ValueError("A liftover method was specified, but target_genome_build equals the source genome build.")
        return frame, default_liftover_metadata(source_build=source, snp_identifier=mode), _empty_liftover_drop_frame()
    if request.method is None:
        raise ValueError("target_genome_build differs from the source genome build, but no liftover method was specified.")

    _require_supported_pair(source, request.target_build)
    missing_mask = ~complete_coordinate_mask(frame)
    missing_drop_frame = _liftover_drop_frame(
        frame,
        missing_mask.to_numpy(dtype=bool),
        reason="missing_coordinate",
        source_pos_col="POS",
    )
    work, missing_report = _normalized_coordinate_frame(frame, drop_missing=True, logger=logger)
    work = work.copy()
    work["_ldsc_source_POS"] = work["POS"].astype(np.int64)
    source_duplicate_result = duplicate_coordinate_drop_result(
        work,
        source_pos_col="_ldsc_source_POS",
        reason="source_duplicate",
    )
    source_duplicate_count = source_duplicate_result.report.n_dropped
    source_duplicate_frame = _liftover_drop_frame(
        work,
        ~source_duplicate_result.keep_mask,
        reason="source_duplicate",
        source_pos_col="_ldsc_source_POS",
    )
    log_liftover_drop_report(
        logger,
        source_duplicate_result.report,
        workflow_label="Summary-statistics liftover",
    )
    work = work.loc[source_duplicate_result.keep_mask].copy()
    if len(work) == 0:
        _raise_sumstats_all_dropped(
            source=source,
            target=str(request.target_build),
            method=str(request.method),
            n_input=len(frame),
            n_missing=missing_report.n_dropped,
            n_source_duplicate=source_duplicate_count,
            n_unmapped=0,
            n_cross_chrom=0,
            n_target_duplicate=0,
        )
    if request.method == "hm3_curated":
        lifted, unmapped_count, method_drop_frame = _apply_hm3_liftover(
            work,
            request,
            source_build=source,
            logger=logger,
        )
        cross_chrom_count = 0
    else:
        lifted, unmapped_count, cross_chrom_count, method_drop_frame = _apply_chain_liftover(
            work,
            request,
            source_build=source,
            logger=logger,
        )

    before_duplicate_drop = len(lifted)
    target_duplicate_result = duplicate_coordinate_drop_result(
        lifted,
        source_pos_col="_ldsc_source_POS",
        target_pos_col="POS",
        reason="target_collision",
    )
    target_duplicate_count = target_duplicate_result.report.n_dropped
    target_duplicate_frame = _liftover_drop_frame(
        lifted,
        ~target_duplicate_result.keep_mask,
        reason="target_collision",
        source_pos_col="_ldsc_source_POS",
        target_pos_col="POS",
    )
    log_liftover_drop_report(
        logger,
        target_duplicate_result.report,
        workflow_label="Summary-statistics liftover",
    )
    if target_duplicate_count:
        lifted = lifted.loc[target_duplicate_result.keep_mask].copy()
    lifted = lifted.drop(columns=["_ldsc_source_POS"], errors="ignore").reset_index(drop=True)
    if len(lifted) == 0:
        _raise_sumstats_all_dropped(
            source=source,
            target=str(request.target_build),
            method=str(request.method),
            n_input=len(frame),
            n_missing=missing_report.n_dropped,
            n_source_duplicate=source_duplicate_count,
            n_unmapped=unmapped_count,
            n_cross_chrom=cross_chrom_count,
            n_target_duplicate=target_duplicate_count,
        )
    report = {
        "applied": True,
        "source_build": source,
        "target_build": request.target_build,
        "method": request.method,
        "chain_file": request.liftover_chain_file if request.method == "chain_file" else None,
        "hm3_map_file": request.hm3_map_file if request.method == "hm3_curated" else None,
        "n_input": int(len(frame)),
        "n_lifted": int(len(lifted)),
        "n_dropped": int(len(frame) - len(lifted)),
        "n_missing_chr_pos_dropped": int(missing_report.n_dropped),
        "n_unmapped": int(unmapped_count),
        "n_cross_chrom": int(cross_chrom_count),
        "n_duplicate_source_dropped": int(source_duplicate_count),
        "n_duplicate_target_dropped": int(target_duplicate_count),
    }
    if target_duplicate_count:
        logger.info(
            f"Removed {target_duplicate_count} summary-statistics rows in duplicate target CHR/POS groups "
            f"after liftover ({before_duplicate_drop} rows before duplicate filtering; {len(lifted)} rows retained)."
        )
    logger.info(
        f"Applied summary-statistics liftover {source} -> {request.target_build} using {request.method}: "
        f"{len(frame)} input rows, {len(lifted)} retained, {report['n_dropped']} dropped."
    )
    drop_frame = _concat_liftover_drop_frames(
        [missing_drop_frame, source_duplicate_frame, method_drop_frame, target_duplicate_frame]
    )
    return lifted, report, drop_frame


def _positive_int_position(values: pd.Series, *, label: str, context: str) -> pd.Series:
    """Return positive integer coordinates with useful validation errors."""
    return positive_int_position_series(values, label=label, context=context).astype(np.int64)


def _default_chain_flag_hint(source_build: str, target_build: str) -> str:
    """Return the reference-panel chain flag for a build pair."""
    return f"--liftover-chain-{source_build}-to-{target_build}-file"


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
    except ValueError:
        return None


def _reject_duplicate_hm3_coordinates(frame: pd.DataFrame, *, build: str) -> None:
    """Raise when a curated HM3 map contains duplicate coordinates for one build."""
    pos_column = f"{build}_POS"
    duplicated = frame.duplicated(subset=["CHR", pos_column], keep=False)
    if duplicated.any():
        examples = frame.loc[duplicated, ["CHR", pos_column, "SNP"]].head(5).to_dict("records")
        raise ValueError(f"HM3 curated map contains duplicate {build} coordinates: {examples}")


def _require_supported_pair(source_build: str, target_build: str | None) -> None:
    """Raise for unsupported build pairs before applying a method."""
    if source_build not in SUPPORTED_LIFTOVER_BUILDS or target_build not in SUPPORTED_LIFTOVER_BUILDS:
        raise ValueError("Summary-statistics liftover supports only hg19 and hg38.")
    if source_build == target_build:
        raise ValueError("source and target genome builds are identical.")


def _normalized_coordinate_frame(
    frame: pd.DataFrame,
    *,
    drop_missing: bool = False,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, Any]:
    """Return a copy with complete normalized CHR/POS columns for liftover."""
    missing_columns = [column for column in ("CHR", "POS") if column not in frame.columns]
    if missing_columns:
        raise ValueError(f"Summary-statistics liftover requires CHR/POS columns; missing {missing_columns}.")
    return normalize_chr_pos_frame(
        frame,
        context="sumstats liftover",
        drop_missing=drop_missing,
        logger=logger,
    )


def _raise_sumstats_all_dropped(
    *,
    source: str,
    target: str,
    method: str,
    n_input: int,
    n_missing: int,
    n_source_duplicate: int,
    n_unmapped: int,
    n_cross_chrom: int,
    n_target_duplicate: int,
) -> None:
    """Raise the standardized hard error when liftover removes every row."""
    raise ValueError(
        "Summary-statistics liftover dropped all rows. "
        f"source_build={source}; target_build={target}; method={method}; "
        f"n_input={n_input}; n_missing_chr_pos_dropped={n_missing}; "
        f"n_duplicate_source_dropped={n_source_duplicate}; "
        f"n_unmapped={n_unmapped}; n_cross_chrom={n_cross_chrom}; "
        f"n_duplicate_target_dropped={n_target_duplicate}."
    )


def _apply_hm3_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    *,
    source_build: str,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, int, pd.DataFrame]:
    """Apply coordinate-only HM3 map liftover and return unmapped sidecar rows."""
    if request.hm3_map_file is None:
        raise ValueError("hm3_map_file is required for HM3 quick liftover.")
    mapping = load_hm3_curated_map(request.hm3_map_file)
    source_column = f"{source_build}_POS"
    target_column = f"{request.target_build}_POS"
    lookup = mapping.loc[:, ["CHR", source_column, target_column]].rename(columns={source_column: "POS"})
    joined = frame.merge(lookup, on=["CHR", "POS"], how="left", sort=False)
    unmapped = joined[target_column].isna()
    unmapped_mask = unmapped.to_numpy(dtype=bool)
    drop_frame = _liftover_drop_frame(
        frame,
        unmapped_mask,
        reason="unmapped_liftover",
        source_pos_col="_ldsc_source_POS" if "_ldsc_source_POS" in frame.columns else "POS",
    )
    if logger is not None:
        log_liftover_drop_report(
            logger,
            liftover_drop_report(
                frame,
                unmapped_mask,
                reason="unmapped_liftover",
                source_pos_col="_ldsc_source_POS" if "_ldsc_source_POS" in frame.columns else "POS",
            ),
            workflow_label="Summary-statistics liftover",
        )
    lifted = frame.loc[~unmapped_mask].copy()
    lifted["POS"] = joined.loc[~unmapped, target_column].to_numpy(dtype=np.int64)
    return lifted, int(unmapped.sum()), drop_frame


def _apply_chain_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    *,
    source_build: str,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, int, int, pd.DataFrame]:
    """Apply chain-file liftover and return unmapped/cross-chromosome sidecar rows."""
    translator = LiftOverTranslator(
        source_build=source_build,
        target_build=str(request.target_build),
        chain_path=request.liftover_chain_file,
        chain_flag_hint="--liftover-chain-file",
        workflow_label="summary-statistics liftover",
    )
    lifted_parts: list[pd.DataFrame] = []
    drop_frames: list[pd.DataFrame] = []
    unmapped_count = 0
    cross_chrom_count = 0
    for chrom, group in frame.groupby("CHR", sort=False):
        result = translator.map_positions(str(chrom), group["POS"].to_numpy(dtype=np.int64))
        unmapped_count += int(result.unmapped_count)
        cross_chrom_count += int(result.cross_chrom_count)
        unmapped_mask, cross_chrom_mask = mapping_reason_masks(result)
        drop_frames.append(
            _liftover_drop_frame(
                group,
                unmapped_mask,
                reason="unmapped_liftover",
                source_pos_col="_ldsc_source_POS" if "_ldsc_source_POS" in group.columns else "POS",
            )
        )
        drop_frames.append(
            _liftover_drop_frame(
                group,
                cross_chrom_mask,
                reason="cross_chromosome_liftover",
                source_pos_col="_ldsc_source_POS" if "_ldsc_source_POS" in group.columns else "POS",
            )
        )
        if logger is not None:
            log_liftover_drop_report(
                logger,
                liftover_drop_report(
                    group,
                    unmapped_mask,
                    reason="unmapped_liftover",
                    source_pos_col="_ldsc_source_POS" if "_ldsc_source_POS" in group.columns else "POS",
                ),
                workflow_label="Summary-statistics liftover",
            )
            log_liftover_drop_report(
                logger,
                liftover_drop_report(
                    group,
                    cross_chrom_mask,
                    reason="cross_chromosome_liftover",
                    source_pos_col="_ldsc_source_POS" if "_ldsc_source_POS" in group.columns else "POS",
                ),
                workflow_label="Summary-statistics liftover",
            )
        kept = group.loc[result.keep_mask].copy()
        kept["POS"] = result.translated_positions.astype(np.int64)
        lifted_parts.append(kept)
    if not lifted_parts:
        return frame.iloc[0:0].copy(), unmapped_count, cross_chrom_count, _concat_liftover_drop_frames(drop_frames)
    return pd.concat(lifted_parts, axis=0), unmapped_count, cross_chrom_count, _concat_liftover_drop_frames(drop_frames)


def _log_drop_examples(
    original: pd.DataFrame,
    retained_index: pd.Index,
    category: str,
    logger: logging.Logger,
) -> None:
    """Log up to five dropped rows for one category at DEBUG verbosity."""
    dropped = original.loc[~original.index.isin(retained_index)]
    if dropped.empty:
        return
    columns = [column for column in ("SNP", "CHR", "POS") if column in dropped.columns]
    examples = dropped.loc[:, columns].head(5).to_dict("records")
    logger.debug(
        f"Dropped {len(dropped)} rows for {category} during summary-statistics liftover. "
        f"Examples: {examples}"
    )
