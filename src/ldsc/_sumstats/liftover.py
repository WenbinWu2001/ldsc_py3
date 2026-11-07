"""Summary-statistics liftover for ``munge-sumstats``.

Core functionality:
    Apply a resolved :class:`SumstatsLiftoverRequest` to munged rows: HM3
    quick liftover through the dual-build HM3 map, or chain-file liftover
    through ``ldsc._genome.liftover.LiftOverTranslator``.

Overview
--------
Liftover is meaningful only in coordinate-family SNP identifier modes. Rows
are dropped in a fixed order (missing coordinates, source duplicate
``CHR``/``POS`` groups, unmapped and cross-chromosome positions, target
duplicate groups); every dropped row enters the munging drop audit with its
reason, and a liftover that drops every row stops the run. The munging
plan chooses the method; the chain translation, duplicate-coordinate
detection and readable drop reports are shared with reference-panel building
in ``ldsc._genome.liftover``.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Any, Sequence

import numpy as np
import pandas as pd

from .._core.coordinates import complete_coordinate_mask, normalize_chr_pos_frame
from .._core.snp_identity import (
    coerce_identity_drop_frame,
    concat_identity_drop_frames,
    empty_identity_drop_frame,
)
from .._genome.liftover import (
    LiftOverTranslator,
    duplicate_coordinate_drop_result,
    liftover_drop_report,
    log_liftover_drop_report,
)
from ..errors import LDSCInputError, LDSCInternalError

LOGGER = logging.getLogger("LDSC.sumstats_liftover")

SUPPORTED_LIFTOVER_BUILDS = {"hg19", "hg38"}
_SUMSTATS_LIFTOVER_TROUBLESHOOTING = (
    "https://wenbinwu2001.github.io/ldsc_py3/user-guide/troubleshooting.html#munge-sumstats-liftover-dropped-all-rows"
)
_WORKFLOW_LABEL = "Summary-statistics liftover"


@dataclass(frozen=True)
class SumstatsLiftoverRequest:
    """Liftover request resolved by the munging plan and executed by the kernel.

    ``ldsc._sumstats.munge_input.liftover_request`` builds it from a
    ``MungeConfig``, which has already normalized and checked the output
    build, so the request holds no option checks of its own. ``target_build``
    is the output build (``None`` in rsID-family modes). A method is set only
    when the resolved source build and the target differ: HM3 quick liftover
    under the packaged HM3 restriction, otherwise the chain file, whose path
    the munging plan resolves.

    ``hm3_map_file`` names the HM3 map file that HM3 quick liftover reads; it
    is provenance for the liftover report, because the positions come from the
    loaded map that :func:`apply_sumstats_liftover` receives as ``hm3_map``.
    ``ldsc._sumstats.munge_input.prepare_munge_input`` sets it to the packaged
    keep-list path, the file it loads that map from.
    """

    target_build: str | None = None
    liftover_chain_file: str | None = None
    use_hm3_quick_liftover: bool = False
    hm3_map_file: str | None = None

    @property
    def method(self) -> str | None:
        """Return the selected liftover method token, if any."""
        if self.liftover_chain_file is not None:
            return "chain_file"
        if self.use_hm3_quick_liftover:
            return "hm3_curated"
        return None


def default_liftover_metadata() -> dict[str, Any]:
    """Return the no-op liftover report used when liftover is not applicable."""
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


def _liftover_drop_frame(
    frame: pd.DataFrame,
    drop_mask: Sequence[bool] | np.ndarray,
    *,
    reason: str,
    source_pos_col: str,
    target_pos_col: str | None = None,
) -> pd.DataFrame:
    """Build dropped-SNP table rows, stage ``liftover``, from one filtered frame slice."""
    mask = np.asarray(drop_mask, dtype=bool)
    if len(mask) != len(frame):
        raise LDSCInternalError(
            "liftover dropped-SNP sidecar construction failed because drop_mask length does not match the input frame. "
            "Most likely an upstream liftover filter produced a mask from a different table version. "
            "Re-run with DEBUG logging and report the traceback."
        )
    if not bool(mask.any()):
        return empty_identity_drop_frame()
    dropped = frame.loc[mask]
    output = pd.DataFrame(index=dropped.index)
    for column, source in (("CHR", "CHR"), ("SNP", "SNP"), ("source_pos", source_pos_col),
                           ("target_pos", target_pos_col)):
        output[column] = dropped[source] if source is not None and source in dropped.columns else pd.NA
    output["reason"] = reason
    return coerce_identity_drop_frame(output, default_stage="liftover")


def apply_sumstats_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    *,
    source_build: str | None,
    hm3_map: pd.DataFrame | None = None,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Apply a summary-statistics liftover request and return metadata plus drops.

    ``request`` comes from the munging plan, which has already stopped on
    every option problem: a request without a method is a no-op (rsID-family
    modes, or matching source and output builds), and a request with a method
    lifts between hg19 and hg38. ``source_build`` is the resolved source
    build. The munger drops
    missing coordinates, source duplicate ``CHR/POS`` groups, unmapped/cross-
    chromosome chain hits, and target duplicate ``CHR/POS`` groups only when a
    liftover request reaches the mapping stage. It updates ``CHR``/``POS`` but
    preserves ``SNP`` as a label field. The third return value lists the
    dropped rows in the dropped-SNP table schema of
    ``ldsc._core.snp_identity.coerce_identity_drop_frame``, with stage
    ``liftover``; clean/no-op runs return the same schema with zero rows.
    HM3 quick liftover maps positions through ``hm3_map``, the dual-build map
    with canonical ``CHR``, ``hg19_POS``, ``hg38_POS`` and ``SNP`` columns
    that ``ldsc._genome.hm3._load_hm3_curated_map_from_path`` loads; the
    munging workflow loads it once and derives the packaged keep-list from
    the same frame.
    """
    logger = LOGGER if logger is None else logger
    source = source_build
    _require_resolved_request(request, source)
    if request.method is None:
        return frame, default_liftover_metadata(), empty_identity_drop_frame()

    missing_mask = ~complete_coordinate_mask(frame)
    missing_drop_frame = _liftover_drop_frame(
        frame,
        missing_mask.to_numpy(dtype=bool),
        reason="missing_coordinate",
        source_pos_col="POS",
    )
    work, missing_report = normalize_chr_pos_frame(frame, context="sumstats liftover", logger=logger)
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
        workflow_label=_WORKFLOW_LABEL,
    )
    work = work.take(np.flatnonzero(source_duplicate_result.keep_mask))
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
            hm3_map,
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
        workflow_label=_WORKFLOW_LABEL,
    )
    if target_duplicate_count:
        lifted = lifted.take(np.flatnonzero(target_duplicate_result.keep_mask))
    lifted = lifted.drop(columns=["_ldsc_source_POS"], errors="ignore")
    lifted.index = pd.RangeIndex(len(lifted))
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
    drop_frame = concat_identity_drop_frames(
        [missing_drop_frame, source_duplicate_frame, method_drop_frame, target_duplicate_frame]
    )
    return lifted, report, drop_frame


def _require_resolved_request(request: SumstatsLiftoverRequest, source_build: str | None) -> None:
    """Assert the shape the munging plan gives a request: a method exactly when hg19 and hg38 differ."""
    lifts = {source_build, request.target_build} == SUPPORTED_LIFTOVER_BUILDS
    if lifts != (request.method is not None):
        raise LDSCInternalError(
            f"Cannot apply summary-statistics liftover from {source_build!r} to {request.target_build!r} "
            f"with method {request.method!r}: a method must be set exactly when the builds are hg19 and hg38 and differ. "
            "Most likely an internal caller built the request without the munging plan. "
            "Build the request with ldsc._sumstats.munge_input.liftover_request() after resolving the source build."
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
    raise LDSCInputError(
        f"Summary-statistics liftover from {source} to {target} using {method} dropped all rows. "
        "Most likely the declared source build does not match the input CHR/POS coordinates, or the liftover method cannot map this SNP set. "
        "Fix the source CHR/POS coordinates, pass the correct --source-genome-build, or use chain-file liftover for non-HM3 SNP sets. "
        f"n_input={n_input}; n_missing_chr_pos_dropped={n_missing}; "
        f"n_duplicate_source_dropped={n_source_duplicate}; "
        f"n_unmapped={n_unmapped}; n_cross_chrom={n_cross_chrom}; "
        f"n_duplicate_target_dropped={n_target_duplicate}. "
        f"Other causes & fixes: {_SUMSTATS_LIFTOVER_TROUBLESHOOTING}"
    )


def _mapping_drops(frame: pd.DataFrame, mask: np.ndarray, reason: str, logger: logging.Logger) -> pd.DataFrame:
    """Return the dropped-SNP rows of positions ``mask`` failed to map, logging their count-only report."""
    drop_frame = _liftover_drop_frame(frame, mask, reason=reason, source_pos_col="_ldsc_source_POS")
    log_liftover_drop_report(
        logger,
        liftover_drop_report(frame, mask, reason=reason, source_pos_col="_ldsc_source_POS"),
        workflow_label=_WORKFLOW_LABEL,
    )
    return drop_frame


def _apply_hm3_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    mapping: pd.DataFrame | None,
    *,
    source_build: str,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, int, pd.DataFrame]:
    """Apply coordinate-only HM3 map liftover and return unmapped sidecar rows."""
    if mapping is None:
        raise LDSCInternalError(
            "munge-sumstats reached HM3 quick liftover without the loaded HM3 map; this is an unexpected internal state. "
            "Report the command, LDSC3 version, and saved traceback. "
            "--log-level DEBUG adds troubleshooting detail."
        )
    source_column = f"{source_build}_POS"
    target_column = f"{request.target_build}_POS"
    lookup = mapping.loc[:, ["CHR", source_column, target_column]].rename(columns={source_column: "POS"})
    # Only the target positions are needed from the map; merge the keys alone.
    joined = frame.loc[:, ["CHR", "POS"]].merge(lookup, on=["CHR", "POS"], how="left", sort=False)
    unmapped = joined[target_column].isna()
    unmapped_mask = unmapped.to_numpy(dtype=bool)
    drop_frame = _mapping_drops(frame, unmapped_mask, "unmapped_liftover", logger)
    lifted = frame.take(np.flatnonzero(~unmapped_mask))
    lifted["POS"] = joined.loc[~unmapped, target_column].to_numpy(dtype=np.int64)
    return lifted, int(unmapped.sum()), drop_frame


def _apply_chain_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    *,
    source_build: str,
    logger: logging.Logger,
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
        unmapped_mask, cross_chrom_mask = result.unmapped_mask, result.cross_chrom_mask
        drop_frames.append(
            _mapping_drops(group, unmapped_mask, "unmapped_liftover", logger)
        )
        drop_frames.append(_mapping_drops(group, cross_chrom_mask, "cross_chromosome_liftover", logger))
        kept = group.loc[result.keep_mask].copy()
        kept["POS"] = result.translated_positions.astype(np.int64)
        lifted_parts.append(kept)
    return pd.concat(lifted_parts, axis=0), unmapped_count, cross_chrom_count, concat_identity_drop_frames(drop_frames)
