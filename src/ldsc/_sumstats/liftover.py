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
workflow chooses the method; the chain translation, duplicate-coordinate
detection and readable drop reports are shared with reference-panel building
in ``ldsc._genome.liftover``.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from os import PathLike
from typing import Any, Sequence
import warnings

import numpy as np
import pandas as pd

from .._core.columns import normalize_genome_build
from .._core.coordinates import complete_coordinate_mask, normalize_chr_pos_frame
from .._core.snp_identity import (
    coerce_identity_drop_frame,
    concat_identity_drop_frames,
    empty_identity_drop_frame,
    identity_mode_family,
    normalize_snp_identifier_mode,
)
from .._genome.hm3 import _load_hm3_curated_map_from_path
from .._genome.liftover import (
    LiftOverTranslator,
    duplicate_coordinate_drop_result,
    liftover_drop_report,
    log_liftover_drop_report,
    mapping_reason_masks,
)
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError

LOGGER = logging.getLogger("LDSC.sumstats_liftover")

SUPPORTED_LIFTOVER_BUILDS = {"hg19", "hg38"}
_SUMSTATS_LIFTOVER_TROUBLESHOOTING = (
    "docs/user-guide/troubleshooting.md#munge-sumstats-liftover-dropped-all-rows"
)


@dataclass(frozen=True)
class SumstatsLiftoverRequest:
    """Resolved liftover request passed from the munger workflow into the kernel.

    ``hm3_map_file`` names the HM3 map file that HM3 quick liftover reads; it
    is provenance for the liftover report, because the positions come from the
    loaded map that :func:`apply_sumstats_liftover` receives as ``hm3_map``.
    ``ldsc._sumstats.munge_input.prepare_munge_input`` sets it to the packaged
    keep-list path, the file it loads that map from.
    """

    target_build: str | None = None
    liftover_chain_file: str | PathLike[str] | None = None
    use_hm3_quick_liftover: bool = False
    hm3_map_file: str | PathLike[str] | None = None

    def __post_init__(self) -> None:
        """Normalize build and path fields, then validate method-level shape."""
        try:
            target_build = normalize_genome_build(self.target_build)
        except LDSCConfigError as exc:
            raise LDSCConfigError(
                f"Could not configure summary-statistics liftover: output_genome_build={self.target_build!r} is invalid. "
                "Most likely the target build is misspelled or unsupported. "
                "Use output_genome_build='hg19' or output_genome_build='hg38'."
            ) from exc
        if target_build == "auto":
            raise LDSCUsageError(
                "Could not configure summary-statistics liftover: output_genome_build='auto' is not a valid liftover target. "
                "Most likely genome-build inference was requested where an explicit target build is required. "
                "Use output_genome_build='hg19' or output_genome_build='hg38'."
            )
        if target_build is not None and target_build not in SUPPORTED_LIFTOVER_BUILDS:
            raise LDSCConfigError(
                f"Could not configure summary-statistics liftover: output_genome_build={target_build!r} is unsupported. "
                "Most likely the requested target build is outside LDSC3 liftover support. "
                "Use output_genome_build='hg19' or output_genome_build='hg38'."
            )
        chain_file = None if self.liftover_chain_file is None else str(self.liftover_chain_file)
        hm3_map_file = None if self.hm3_map_file is None else str(self.hm3_map_file)
        object.__setattr__(self, "target_build", target_build)
        object.__setattr__(self, "liftover_chain_file", chain_file)
        object.__setattr__(self, "hm3_map_file", hm3_map_file)
        if chain_file is not None and self.use_hm3_quick_liftover:
            raise LDSCUsageError(
                "Could not configure summary-statistics liftover: liftover_chain_file and use_hm3_quick_liftover are mutually exclusive. "
                "Most likely two liftover methods were selected for one run. "
                "Choose either chain-file liftover or HM3 quick liftover, not both."
            )
        if self.method is not None and target_build is None:
            raise LDSCUsageError(
                "Could not configure summary-statistics liftover: output_genome_build is required when a liftover method is specified. "
                "Most likely a method flag was set without the target build. "
                "Set output_genome_build='hg19' or output_genome_build='hg38'."
            )

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
    chrom_col: str = "CHR",
    source_pos_col: str | None = "POS",
    target_pos_col: str | None = None,
    snp_col: str = "SNP",
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
    return coerce_identity_drop_frame(output, default_stage="liftover")


def load_hm3_map_file(path: str | PathLike[str], *, alleles: bool = False) -> pd.DataFrame:
    """Load a curated dual-build HM3 map file into canonical columns.

    The frame has ``CHR``, ``hg19_POS``, ``hg38_POS`` and ``SNP``.
    ``alleles=True`` also requires ``A1``/``A2``, returns them uppercase and
    keeps the file's other columns. Unlike ``ldsc.load_hm3_curated_map``,
    which loads the packaged map, this reads any map file.
    """
    return _load_hm3_curated_map_from_path(path, preserve_extra_columns=alleles)


def apply_sumstats_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest | None,
    *,
    source_build: str | None,
    snp_identifier: str,
    hm3_map: pd.DataFrame | None = None,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    """Apply a summary-statistics liftover request and return metadata plus drops.

    Liftover is meaningful only in coordinate-family modes. The munger drops
    missing coordinates, source duplicate ``CHR/POS`` groups, unmapped/cross-
    chromosome chain hits, and target duplicate ``CHR/POS`` groups only when a
    liftover request reaches the mapping stage. It updates ``CHR``/``POS`` but
    preserves ``SNP`` as a label field. The third return value lists the
    dropped rows in the dropped-SNP table schema of
    ``ldsc._core.snp_identity.coerce_identity_drop_frame``, with stage
    ``liftover``; clean/no-op runs return the same schema with zero rows.
    HM3 quick liftover maps positions through ``hm3_map``, the dual-build map
    from :func:`load_hm3_map_file`; the munging workflow loads it once and
    derives the packaged keep-list from the same frame.
    """
    logger = LOGGER if logger is None else logger
    request = request or SumstatsLiftoverRequest()
    mode = normalize_snp_identifier_mode(snp_identifier)
    source = normalize_genome_build(source_build)
    if request.requested and identity_mode_family(mode) != "chr_pos":
        raise LDSCUsageError(
            f"Cannot apply summary-statistics liftover with snp_identifier={mode!r}: liftover is only valid for chr_pos-family modes. "
            "Most likely liftover was requested for rsID-only summary statistics. "
            "Use chr_pos or chr_pos_allele_aware identifiers, or omit liftover for rsID-family inputs."
        )
    if not request.requested:
        return frame, default_liftover_metadata(), empty_identity_drop_frame()
    if source not in SUPPORTED_LIFTOVER_BUILDS:
        raise LDSCUsageError(
            "Cannot apply summary-statistics liftover: the source genome build is unresolved. "
            "Most likely source genome-build inference was skipped or had insufficient coordinate evidence. "
            "Use --source-genome-build hg19/hg38 or --source-genome-build auto with inferable coordinates."
        )
    if request.target_build == source:
        if request.method is not None:
            message = (
                "A summary-statistics liftover method was specified, but output_genome_build "
                "equals the resolved source genome build; the liftover method will be ignored."
            )
            warnings.warn(message, UserWarning, stacklevel=2)
            logger.warning(message)
        return frame, default_liftover_metadata(), empty_identity_drop_frame()
    if request.method is None:
        raise LDSCUsageError(
            f"Cannot apply summary-statistics liftover from {source} to {request.target_build}: no liftover method was specified. "
            "Most likely the target build differs from the source build but neither a chain file nor HM3 quick liftover was selected. "
            "Provide --liftover-chain-file for custom-list or unrestricted SNPs. "
            "The public munge-sumstats workflow selects quick liftover automatically under packaged HM3 restriction."
        )

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
        workflow_label="Summary-statistics liftover",
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


def _require_supported_pair(source_build: str, target_build: str | None) -> None:
    """Raise for unsupported build pairs before applying a method."""
    if source_build not in SUPPORTED_LIFTOVER_BUILDS or target_build not in SUPPORTED_LIFTOVER_BUILDS:
        raise LDSCUsageError(
            f"Cannot apply summary-statistics liftover from {source_build!r} to {target_build!r}: only hg19 and hg38 are supported. "
            "Most likely an unsupported source or target build was requested. "
            "Use hg19/hg38 or a workflow-specific chain liftover path that supports the requested pair."
        )
    if source_build == target_build:
        raise LDSCInternalError(
            f"Cannot apply summary-statistics liftover from {source_build} to {target_build}: source and target builds are identical. "
            "Most likely an internal caller failed to short-circuit a no-op liftover. "
            "Return no-op liftover metadata before calling the mapping stage."
        )


def _normalized_coordinate_frame(
    frame: pd.DataFrame,
    *,
    drop_missing: bool = False,
    logger: logging.Logger | None = None,
) -> tuple[pd.DataFrame, Any]:
    """Return a copy with complete normalized CHR/POS columns for liftover."""
    missing_columns = [column for column in ("CHR", "POS") if column not in frame.columns]
    if missing_columns:
        raise LDSCInputError(
            f"Cannot apply summary-statistics liftover: missing CHR/POS column(s) {missing_columns}. "
            "Most likely the raw summary statistics were parsed without coordinate columns. "
            "Provide CHR and POS columns or use a chr_pos-compatible input format."
        )
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


def _apply_hm3_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    mapping: pd.DataFrame | None,
    *,
    source_build: str,
    logger: logging.Logger | None = None,
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
    lifted = frame.take(np.flatnonzero(~unmapped_mask))
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
        return frame.iloc[0:0].copy(), unmapped_count, cross_chrom_count, concat_identity_drop_frames(drop_frames)
    return pd.concat(lifted_parts, axis=0), unmapped_count, cross_chrom_count, concat_identity_drop_frames(drop_frames)
