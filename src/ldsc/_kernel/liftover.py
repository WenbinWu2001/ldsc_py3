"""Shared coordinate liftover helpers for LDSC workflows."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import logging
from os import PathLike
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd

from ..chromosome_inference import normalize_chromosome, normalize_chromosome_series
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
    """Partial liftover result for one chromosome."""

    translated_positions: np.ndarray
    keep_mask: np.ndarray
    unmapped_count: int
    cross_chrom_count: int


class LiftOverTranslator:
    """Translate one chromosome's 1-based positions between explicit genome-build chain files."""

    def __init__(
        self,
        *,
        source_build: str,
        target_build: str,
        chain_path: str | PathLike[str] | None = None,
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
            flag = (
                "--liftover-chain-hg19-to-hg38"
                if (source_build, target_build) == ("hg19", "hg38")
                else "--liftover-chain-hg38-to-hg19"
            )
            raise ValueError(
                f"An explicit liftover chain path is required for {source_build} -> {target_build}. "
                f"Provide it via {flag}."
            )
        try:
            from pyliftover import LiftOver
        except ImportError as exc:
            raise LDSCDependencyError(
                "Building reference panels or munging summary statistics across hg19/hg38 requires "
                "the optional dependency 'pyliftover'."
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
            )
        chrom = _normalize_liftover_chromosome(chrom)
        query_chrom = f"chr{chrom}"
        translated = np.zeros(len(array), dtype=np.int64)
        keep_mask = np.zeros(len(array), dtype=bool)
        unmapped_count = 0
        cross_chrom_count = 0
        for idx, pos in enumerate(array):
            hits = self._liftover.convert_coordinate(query_chrom, int(pos) - 1)
            if not hits:
                unmapped_count += 1
                continue
            hit = next(
                (candidate for candidate in hits if _try_normalize_liftover_chromosome(candidate[0]) == chrom),
                None,
            )
            if hit is None:
                cross_chrom_count += 1
                continue
            translated[idx] = int(hit[1]) + 1
            keep_mask[idx] = True
        return LiftOverMappingResult(
            translated_positions=translated[keep_mask],
            keep_mask=keep_mask,
            unmapped_count=unmapped_count,
            cross_chrom_count=cross_chrom_count,
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
        work = _normalized_coordinate_frame(frame)
        request = SumstatsLiftoverRequest(
            target_build=self.target_build,
            use_hm3_quick_liftover=True,
            hm3_map_file=self.map_path,
        )
        lifted, _unmapped_count = _apply_hm3_liftover(work, request, source_build=self.source_build)
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
        "n_unmapped": None,
        "n_cross_chrom": None,
        "n_duplicate_target_dropped": None,
    }


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
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Apply a summary-statistics liftover request and return updated metadata."""
    logger = LOGGER if logger is None else logger
    request = request or SumstatsLiftoverRequest()
    mode = normalize_snp_identifier_mode(snp_identifier)
    source = normalize_genome_build(source_build)
    if request.requested and mode != "chr_pos":
        raise ValueError("Summary-statistics liftover is only valid when snp_identifier='chr_pos'.")
    if not request.requested:
        return frame, default_liftover_metadata(source_build=source, snp_identifier=mode)
    if source not in SUPPORTED_LIFTOVER_BUILDS:
        raise ValueError(
            "Cannot apply summary-statistics liftover because the source genome build is unresolved. "
            "Use --genome-build hg19/hg38 or --genome-build auto with inferable coordinates."
        )
    if request.target_build == source:
        if request.method is not None:
            raise ValueError("A liftover method was specified, but target_genome_build equals the source genome build.")
        return frame, default_liftover_metadata(source_build=source, snp_identifier=mode)
    if request.method is None:
        raise ValueError("target_genome_build differs from the source genome build, but no liftover method was specified.")

    _require_supported_pair(source, request.target_build)
    work = _normalized_coordinate_frame(frame)
    if request.method == "hm3_curated":
        lifted, unmapped_count = _apply_hm3_liftover(work, request, source_build=source)
        cross_chrom_count = 0
    else:
        lifted, unmapped_count, cross_chrom_count = _apply_chain_liftover(work, request, source_build=source)

    _log_drop_examples(frame, lifted.index, "unmapped/cross-chromosome liftover", logger)
    before_duplicate_drop = len(lifted)
    duplicate_mask = lifted.duplicated(subset=["CHR", "POS"], keep=False)
    duplicate_count = int(duplicate_mask.sum())
    if duplicate_count:
        _log_drop_examples(lifted, lifted.index[~duplicate_mask], "duplicate target coordinate", logger)
        lifted = lifted.loc[~duplicate_mask].copy()
    lifted = lifted.reset_index(drop=True)
    if len(lifted) == 0:
        raise ValueError(
            "Summary-statistics liftover dropped all rows. "
            f"source_build={source}; target_build={request.target_build}; method={request.method}; "
            f"n_input={len(frame)}; n_unmapped={unmapped_count}; n_cross_chrom={cross_chrom_count}; "
            f"n_duplicate_target_dropped={duplicate_count}."
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
        "n_unmapped": int(unmapped_count),
        "n_cross_chrom": int(cross_chrom_count),
        "n_duplicate_target_dropped": int(duplicate_count),
    }
    if duplicate_count:
        logger.info(
            "Removed %d summary-statistics rows in duplicate target CHR/POS groups after liftover "
            "(%d rows before duplicate filtering; %d rows retained).",
            duplicate_count,
            before_duplicate_drop,
            len(lifted),
        )
    logger.info(
        "Applied summary-statistics liftover %s -> %s using %s: %d input rows, %d retained, %d dropped.",
        source,
        request.target_build,
        request.method,
        len(frame),
        len(lifted),
        report["n_dropped"],
    )
    return lifted, report


def _positive_int_position(values: pd.Series, *, label: str, context: str) -> pd.Series:
    """Return positive integer coordinates with useful validation errors."""
    numeric = pd.to_numeric(values, errors="raise")
    if numeric.isna().any():
        raise ValueError(f"{context} contains missing {label} values.")
    non_integral = (numeric % 1) != 0
    if bool(non_integral.any()):
        raise ValueError(f"{context} contains non-integer {label} values.")
    if bool((numeric <= 0).any()):
        raise ValueError(f"{context} contains non-positive {label} values.")
    return numeric.astype(np.int64)


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


def _normalized_coordinate_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with complete normalized CHR/POS columns for liftover."""
    missing_columns = [column for column in ("CHR", "POS") if column not in frame.columns]
    if missing_columns:
        raise ValueError(f"Summary-statistics liftover requires CHR/POS columns; missing {missing_columns}.")
    missing_chr = _missing_coordinate_mask(frame["CHR"])
    missing_pos = _missing_coordinate_mask(frame["POS"])
    if bool((missing_chr | missing_pos).any()):
        raise ValueError("Summary-statistics liftover requires complete CHR/POS coordinates for every retained row.")
    work = frame.copy()
    work["CHR"] = normalize_chromosome_series(work["CHR"], context="sumstats liftover").astype(object)
    work["POS"] = _positive_int_position(work["POS"], label="POS", context="sumstats liftover")
    return work


def _missing_coordinate_mask(values: pd.Series) -> pd.Series:
    """Return missing or NA-like coordinate tokens."""
    tokens = values.astype("string")
    return tokens.isna() | tokens.str.strip().str.lower().isin({"", "na", "nan", "none"})


def _apply_hm3_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    *,
    source_build: str,
) -> tuple[pd.DataFrame, int]:
    """Apply coordinate-only HM3 map liftover."""
    if request.hm3_map_file is None:
        raise ValueError("hm3_map_file is required for HM3 quick liftover.")
    mapping = load_hm3_curated_map(request.hm3_map_file)
    source_column = f"{source_build}_POS"
    target_column = f"{request.target_build}_POS"
    lookup = mapping.loc[:, ["CHR", source_column, target_column]].rename(columns={source_column: "POS"})
    joined = frame.merge(lookup, on=["CHR", "POS"], how="left", sort=False)
    unmapped = joined[target_column].isna()
    lifted = frame.loc[~unmapped.to_numpy()].copy()
    lifted["POS"] = joined.loc[~unmapped, target_column].to_numpy(dtype=np.int64)
    return lifted, int(unmapped.sum())


def _apply_chain_liftover(
    frame: pd.DataFrame,
    request: SumstatsLiftoverRequest,
    *,
    source_build: str,
) -> tuple[pd.DataFrame, int, int]:
    """Apply chain-file liftover and drop unmapped or cross-chromosome hits."""
    translator = LiftOverTranslator(
        source_build=source_build,
        target_build=str(request.target_build),
        chain_path=request.liftover_chain_file,
    )
    lifted_parts: list[pd.DataFrame] = []
    unmapped_count = 0
    cross_chrom_count = 0
    for chrom, group in frame.groupby("CHR", sort=False):
        result = translator.map_positions(str(chrom), group["POS"].to_numpy(dtype=np.int64))
        unmapped_count += int(result.unmapped_count)
        cross_chrom_count += int(result.cross_chrom_count)
        kept = group.loc[result.keep_mask].copy()
        kept["POS"] = result.translated_positions.astype(np.int64)
        lifted_parts.append(kept)
    if not lifted_parts:
        return frame.iloc[0:0].copy(), unmapped_count, cross_chrom_count
    return pd.concat(lifted_parts, axis=0), unmapped_count, cross_chrom_count


def _log_drop_examples(
    original: pd.DataFrame,
    retained_index: pd.Index,
    category: str,
    logger: logging.Logger,
) -> None:
    """Log up to five dropped rows for one category."""
    dropped = original.loc[~original.index.isin(retained_index)]
    if dropped.empty:
        return
    columns = [column for column in ("SNP", "CHR", "POS") if column in dropped.columns]
    examples = dropped.loc[:, columns].head(5).to_dict("records")
    logger.info("Dropped %d rows for %s during summary-statistics liftover. Examples: %s", len(dropped), category, examples)
