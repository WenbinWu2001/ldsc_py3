# Includes primitives derived from legacy LDSC ldscore/ldscore.py.
# GPLv3 (LICENSE); upstream authors and Python 3 port credit are in NOTICE.
# Modified 2026 for prepared reference panels, streaming R2, and projection.
# See NOTICE for attribution and dated modifications.

"""LD-score projection primitives and streaming R2 panel reader.

``RefPanel.prepare_chromosome`` owns source resolution, reference filtering,
annotation alignment, LD-window validation, and reader lifetime. This module's
``compute_chromosome`` consumes that ``PreparedChromosome`` state and projects
baseline/query annotations and regression weights in one traversal. It also
provides annotation parsing, window geometry, and count primitives used by the
workflow adapters. BED-to-SNP projection belongs to ``ldsc.annotation_builder``.

Scientific contracts
--------------------
Prepared metadata, selected annotation rows, and window bounds describe the
same retained reference SNPs in genomic order. Binary annotations can remain
packed until selected reads; continuous values retain float32 normalization.
Projection accumulates in float64. Reference-panel CM and MAF are
authoritative. PLINK sample selection precedes genotype QC and inclusive MAF
filtering; chromosome LD scores retain the existing float32 result conversion.
Parquet projection accumulates in float64. Common-SNP counts and overlap use
``MAF >= common_maf_min`` (0.05 by default), independently of panel filtering.

R2 panel Parquet files contain sidecar-row indices ``IDX_1``/``IDX_2``,
``R2``, and ``SIGN_R``. The required sidecar is bound by its identity digest and
row count. Each chromosome builds one identity remap; pairs stream through
``iter_all_pairs`` without a dense chromosome matrix or row-group cache. Every
file must match the ``build-r2-panel`` metadata contract checked by
``read_r2_parquet_contract``: int16 ``R2`` quantized with scale 32767 and already
bias-corrected, so decoding divides by the scale and applies no further
correction. Each unordered off-diagonal pair contributes in both directions,
and the diagonal is added internally as one.

This module returns in-memory results and writes no artifacts. Public result
normalization and aggregation live in ``ldsc.ldscore_calculator``; canonical
output layout belongs to ``ldsc.outputs``. See
``docs/reference/formats/r2-panels.md`` for the file contract.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .._coordinates import CHR_POS_KEY_COLUMN, build_chr_pos_key_frame
from ..column_inference import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    CHR_COLUMN_ALIASES,
    CHR_COLUMN_SPEC,
    CM_COLUMN_ALIASES,
    CM_COLUMN_SPEC,
    MAF_COLUMN_SPEC,
    POS_COLUMN_ALIASES,
    POS_COLUMN_SPEC,
    REFERENCE_METADATA_SPEC_MAP,
    SNP_COLUMN_ALIASES,
    SNP_COLUMN_SPEC,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_optional_column,
    resolve_required_column,
)
from ..chromosome_inference import chrom_sort_key, normalize_chromosome
from ..errors import LDSCConfigError, LDSCDependencyError, LDSCInputError, LDSCInternalError, LDSCUsageError
from . import formats as legacy_parse
from .identifiers import build_snp_id_series
from .overlap import OverlapContribution, compute_overlap, annotation_statistics
from .plink_bed import PlinkBEDFile  # Re-exported for the reference-panel adapters.
from .ldscore_projection import ArrayAnnotations, ProjectionAccumulator
from .snp_identity import (
    RestrictionIdentityKeys,
    empty_identity_drop_frame,
    effective_merge_key_series,
    identity_mode_family,
    is_allele_aware_mode,
    restriction_membership_mask,
    sidecar_identity_sha256,
)
from .regions import RegionIntervals, region_exclusion_keep_mask


LOGGER = logging.getLogger("LDSC.ldscore")
_LDSCORE_INTERSECTION_DOC = (
    "docs/user-guide/troubleshooting.md#ldscore-no-annotation-snps-remain-after-reference-panel-intersection"
)
_LDSCORE_PARQUET_DOC = "docs/user-guide/troubleshooting.md#ldscore-r2-panel-is-incompatible"
# K bounds buffered pairs; CSR construction and projection add working arrays.
_CSR_CHUNK_PAIRS = 16_000_000
CHROM_ALIASES = CHR_COLUMN_ALIASES
POS_ALIASES = POS_COLUMN_ALIASES
SNP_ALIASES = SNP_COLUMN_ALIASES
CM_ALIASES = CM_COLUMN_ALIASES
MAF_ALIASES = ("MAF", "FRQ", "FREQ", "FREQUENCY")


@dataclass
class AnnotationBundle:
    """Chromosome metadata and a selected-read annotation source.

    Explicitly supplied DataFrames are numerical inputs owned by the caller;
    workflow preparation supplies a shard descriptor instead.
    """
    metadata: pd.DataFrame
    annotations: object
    baseline_columns: list[str]
    query_columns: list[str]

    def __post_init__(self):
        if isinstance(self.annotations, pd.DataFrame):
            self.annotations = ArrayAnnotations(
                self.annotations.to_numpy(dtype=np.float32, copy=False),
                tuple(self.baseline_columns + self.query_columns))


@dataclass
class ChromComputationResult:
    """Chromosome-level LD-score outputs before cross-chromosome aggregation."""
    chrom: str
    metadata: pd.DataFrame
    ld_scores: np.ndarray
    w_ld: np.ndarray
    M: np.ndarray
    M_5_50: np.ndarray
    ldscore_columns: list[str]
    baseline_columns: list[str]
    query_columns: list[str]
    reference_snp_count: int = 0
    regression_selected_snp_count: int = 0
    regression_region_removed_snp_count: int = 0
    annotation_types: dict = field(default_factory=dict)
    overlap: OverlapContribution | None = None
    identity_drops: pd.DataFrame = field(default_factory=empty_identity_drop_frame)


@dataclass
class PreparedChromosome:
    """Aligned reference state whose reader belongs to one chromosome calculation.

    ``annotations`` maps selected reads to source rows in ``metadata`` order.
    ``block_left`` describes all of those reference SNPs. Use as a context
    manager; closing releases the reader and chromosome working state.
    """

    backend: str
    reader: object
    metadata: pd.DataFrame
    annotations: object
    block_left: np.ndarray
    baseline_columns: list[str]
    query_columns: list[str]
    reference_rows_before_genotype_qc: int = 0
    genotype_qc_removed: int = 0
    maf_removed: int = 0
    selected_individual_count: int = 0
    cm_source: str = "bim_cm"
    identity_drops: pd.DataFrame = field(default_factory=empty_identity_drop_frame)

    def close(self):
        """Release chromosome state once, including when reader closure fails."""
        reader = self.reader
        self.reader = None
        self.annotations = None
        self.metadata = None
        self.block_left = None
        self.identity_drops = None
        if reader is not None:
            reader.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


# Shared computational helpers.


def find_column(columns: Iterable[str], aliases: Sequence[str]) -> str | None:
    """Return the first column whose normalized name matches one of ``aliases``."""
    spec = {
        tuple(CHROM_ALIASES): CHR_COLUMN_SPEC,
        tuple(POS_ALIASES): POS_COLUMN_SPEC,
        tuple(SNP_ALIASES): SNP_COLUMN_SPEC,
        tuple(CM_ALIASES): CM_COLUMN_SPEC,
        tuple(MAF_ALIASES): MAF_COLUMN_SPEC,
    }.get(tuple(aliases))
    if spec is not None:
        return resolve_optional_column(columns, spec)
    return None


def get_block_lefts(coords: np.ndarray, max_dist: float) -> np.ndarray:
    """
    Compute the standard LDSC ``block_left`` window start array.

    ``block_left[i]`` stores the first index retained in the LD window centered
    on SNP ``i`` for the chosen distance metric.
    """
    M = len(coords)
    j = 0
    block_left = np.zeros(M)
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


def identifier_keys(df: pd.DataFrame, mode: str) -> pd.Series:
    """Build the canonical SNP identifier series used for matching within the kernel."""
    mode = normalize_snp_identifier_mode(mode)
    if identity_mode_family(mode) == "rsid" or is_allele_aware_mode(mode):
        return build_snp_id_series(df, mode)
    keyed, _report = build_chr_pos_key_frame(
        df,
        context="LD-score SNP matching",
        drop_missing=True,
        logger=LOGGER,
    )
    keys = pd.Series(pd.NA, index=df.index, dtype="object")
    keys.loc[keyed.index] = keyed[CHR_POS_KEY_COLUMN].astype(str)
    return keys


def build_index_remap(
    full_sidecar: pd.DataFrame,
    retained_metadata: pd.DataFrame,
    identifier_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Map panel (build) indices to retained matrix indices for an R2 panel.

    ``full_sidecar`` is the complete panel in build order (the parquet IDX
    space). ``retained_metadata`` lists the retained reference SNPs in
    matrix order. Returns ``(remap, retained_build_idx)`` where
    ``remap[build_idx]`` is the retained matrix index or ``-1``, and
    ``retained_build_idx[matrix_idx]`` is the originating build index (ascending,
    used for IDX_1 row-group pruning). Matching uses the same mode-dependent
    identity keys as the legacy per-pair decode, so the result is bit-identical.

    The remap must be injective over retained indices: two distinct panel rows
    sharing an identity key under ``identifier_mode`` (a collapse) is a hard
    error, which also keeps ``retained_build_idx`` a well-defined permutation.
    """
    mode = normalize_snp_identifier_mode(identifier_mode)
    full_keys = effective_merge_key_series(full_sidecar, mode, context="panel sidecar index remap").to_numpy()
    retained_keys = effective_merge_key_series(retained_metadata, mode, context="retained metadata index remap").to_numpy()
    retained_index = pd.Index(retained_keys)
    if retained_index.has_duplicates:
        raise LDSCInputError(
            "ldscore could not align retained reference metadata to the R2 panel sidecar. "
            f"The retained metadata has duplicate SNP identity keys under mode '{mode}'. "
            "Most likely the SNP identifier mode is too coarse for this reference panel. "
            "Use an allele-aware `--snp-identifier` mode or rebuild the reference panel "
            "after removing duplicate SNP identities."
        )
    remap = retained_index.get_indexer(full_keys).astype(np.int32, copy=False)

    valid = remap >= 0
    valid_targets = remap[valid]
    if np.unique(valid_targets).size != valid_targets.size:
        raise LDSCInputError(
            "ldscore could not align the reference panel to its sidecar: distinct "
            f"panel SNPs collapse onto the same retained SNP under mode '{mode}'. "
            "Most likely the panel sidecar has duplicate SNP identities at this "
            "identifier resolution. Use an allele-aware `--snp-identifier` mode, or "
            "rebuild the reference panel after removing duplicate SNP identities. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )

    m = len(retained_metadata)
    retained_build_idx = np.empty(m, dtype=np.int64)
    retained_build_idx[valid_targets] = np.nonzero(valid)[0]
    return remap, retained_build_idx


def sort_frame_by_genomic_position(df: pd.DataFrame) -> pd.DataFrame:
    """Sort by chromosome/position and use a present SNP label as a stable tie-breaker."""
    pos_col = find_column(df.columns, POS_ALIASES)
    if pos_col is None:
        raise KeyError("No POS-like column available for genomic sorting.")
    sort_df = df.copy()
    sort_df["_chrom_key"] = sort_df["CHR"].map(chrom_sort_key)
    sort_columns = ["_chrom_key", pos_col]
    if "SNP" in sort_df.columns:
        sort_columns.append("SNP")
    sort_df = sort_df.sort_values(by=sort_columns, kind="mergesort")
    return sort_df.drop(columns="_chrom_key").reset_index(drop=True)


def _panel_sidecar_path_for_r2(r2_path: str) -> Path:
    """Return the ``chrN_meta.tsv.gz`` sidecar path paired with an R2 Parquet file."""
    p = Path(r2_path)
    if not p.name.endswith("_r2.parquet"):
        raise LDSCInputError(
            f"ldscore could not locate the sidecar for R2 Parquet file '{r2_path}'. Most likely "
            "the parquet file does not use the canonical `chrN_r2.parquet` filename written "
            "by `ldsc build-r2-panel`. Regenerate the reference panel or pass the "
            "R2 panel directory. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    return p.with_name(p.name[: -len("_r2.parquet")] + "_meta.tsv.gz")


def _load_full_panel_sidecar(r2_path: str) -> pd.DataFrame:
    """Load the complete (unrestricted) panel sidecar that defines the index space."""
    sidecar_path = _panel_sidecar_path_for_r2(r2_path)
    if not sidecar_path.exists():
        raise LDSCInputError(
            f"ldscore could not load R2 Parquet file '{r2_path}' because the "
            f"required sidecar '{sidecar_path}' is missing. Most likely the parquet file "
            "was copied without its matching `chrN_meta.tsv.gz` sidecar. Keep the R2 "
            "Parquet file and sidecar together or regenerate with `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    df = pd.read_csv(sidecar_path, sep="\t", comment="#")
    context = f"panel sidecar {sidecar_path}"
    renamed = {
        resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=context): "CHR",
        resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=context): "POS",
        resolve_required_column(df.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=context): "SNP",
        resolve_required_column(df.columns, A1_COLUMN_SPEC, context=context): "A1",
        resolve_required_column(df.columns, A2_COLUMN_SPEC, context=context): "A2",
    }
    cm_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=context)
    if cm_col is not None:
        renamed[cm_col] = "CM"
    return df.rename(columns=renamed)


def _validate_index_binding(full_sidecar: pd.DataFrame, *, n_snps: int, identity_hash: str, context: str) -> None:
    """Hard-fail if the sidecar does not match the parquet's recorded binding."""
    if len(full_sidecar) != int(n_snps):
        raise LDSCInputError(
            f"ldscore could not use R2 panel input at {context}: the sidecar has "
            f"{len(full_sidecar)} rows but the parquet records n_snps={n_snps}. Most "
            "likely the parquet and sidecar are not a matched pair. Restore the matching "
            "sidecar or regenerate the reference panel. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    actual = sidecar_identity_sha256(full_sidecar, context=context)
    if actual != identity_hash:
        raise LDSCInputError(
            f"ldscore could not use R2 panel input at {context}: the sidecar identity "
            f"hash {actual} does not match parquet ldsc:sidecar_identity_sha256 "
            f"{identity_hash}. Most likely the sidecar is wrong, reordered, or edited. "
            "Restore the matching sidecar or regenerate the reference panel. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )


def _arrow_column_to_numpy(column):
    """Convert an Arrow column to NumPy across PyArrow array API versions."""
    try:
        return column.to_numpy(zero_copy_only=False)
    except TypeError:
        return column.to_numpy()


def validate_retained_identifier_uniqueness(metadata: pd.DataFrame, identifier_mode: str, chrom: str) -> None:
    """Reject ambiguous retained SNP identifiers before parquet matching begins."""
    if is_allele_aware_mode(identifier_mode):
        keys = effective_merge_key_series(
            metadata,
            identifier_mode,
            context=f"retained SNP metadata for chromosome {chrom}",
        )
        duplicated = keys.duplicated(keep=False)
        if duplicated.any():
            raise LDSCInputError(
                f"Chromosome {chrom} has duplicate retained SNP identities. "
                f"This is ambiguous in {identifier_mode} mode. Most likely the selected "
                "SNP identifier mode is too coarse for the retained reference panel. "
                "Use an allele-aware mode or remove duplicate retained SNP identities."
            )
        return

    if identity_mode_family(identifier_mode) == "chr_pos":
        duplicated = metadata.duplicated(subset=["CHR", "POS"], keep=False)
        if duplicated.any():
            raise LDSCInputError(
                f"Chromosome {chrom} has duplicate retained SNP positions. "
                "This is ambiguous in base chr_pos mode. Most likely multiple alleles "
                "share the same CHR/POS coordinate. Use an allele-aware chr_pos mode or "
                "remove duplicate retained positions."
            )
        return

    duplicated = metadata["SNP"].duplicated(keep=False)
    if duplicated.any():
        raise LDSCInputError(
            f"Chromosome {chrom} has duplicate retained SNP identities. "
            "This is ambiguous in `rsid` mode. Most likely the retained reference "
            "panel contains duplicate rsIDs. Use `rsid_allele_aware` or remove "
            "duplicate retained SNP identities."
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


def build_window_coordinates(
    metadata: pd.DataFrame,
    args: argparse.Namespace,
    *,
    snp_positions: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    """Build the coordinate array and maximum distance for the active LD window mode.

    ``--ld-wind-snps`` counts retained rows unless ``snp_positions`` gives each
    retained SNP's row in a larger SNP set (the R2 panel sidecar),
    in which case it counts rows of that set. ``snp_positions`` must be
    increasing; parquet preparation orders retained SNPs by panel row.
    """
    selectors = np.array([args.ld_wind_snps is not None, args.ld_wind_kb is not None, args.ld_wind_cm is not None], dtype=bool)
    if selectors.sum() != 1:
        raise LDSCUsageError(
            "ldscore could not choose an LD-window mode. Most likely zero or multiple "
            "LD-window options were supplied. Specify exactly one of `--ld-wind-snps`, "
            "`--ld-wind-kb`, or `--ld-wind-cm`."
        )

    if args.ld_wind_snps is not None:
        if snp_positions is None:
            return np.arange(len(metadata), dtype=float), float(args.ld_wind_snps)
        return np.asarray(snp_positions, dtype=float), float(args.ld_wind_snps)
    if args.ld_wind_kb is not None:
        return metadata["POS"].to_numpy(dtype=float), float(args.ld_wind_kb) * 1000.0
    if metadata["CM"].isna().any():
        raise LDSCInputError(
            "ldscore cannot use `--ld-wind-cm` because retained SNP metadata contains "
            "missing CM values. Most likely the annotation or frequency sidecar lacks "
            "genetic-map positions for at least one retained SNP. Provide complete CM "
            "metadata or use `--ld-wind-kb` / `--ld-wind-snps`."
        )
    return metadata["CM"].to_numpy(dtype=float), float(args.ld_wind_cm)


def check_whole_chromosome_window(block_left: np.ndarray, args: argparse.Namespace, chrom: str) -> None:
    """Guard against accidental whole-chromosome windows unless explicitly allowed."""
    if len(block_left) == 0:
        return
    if block_left[-1] == 0 and not args.yes_really:
        raise LDSCUsageError(
            f"ldscore would use a whole-chromosome LD window on chromosome {chrom}. "
            "Most likely the LD-window setting is too large for this chromosome. "
            "Use a smaller LD window, or rerun with `--yes-really` if a whole-chromosome "
            "window is intentional."
        )


@dataclass(frozen=True)
class _LDWindowSpec:
    mode: str
    value: float

    @property
    def flag(self) -> str:
        return {
            "snps": "--ld-wind-snps",
            "kb": "--ld-wind-kb",
            "cm": "--ld-wind-cm",
        }[self.mode]


def _format_ld_window_spec(spec: _LDWindowSpec) -> str:
    """Return the CLI spelling of one LD-window specification."""
    value = int(spec.value) if spec.mode == "snps" else spec.value
    return f"{spec.flag} {value}"


def _requested_ld_window_spec(args: argparse.Namespace) -> _LDWindowSpec:
    """Return the user-requested LD-window mode and value from parsed args."""
    if args.ld_wind_snps is not None:
        return _LDWindowSpec("snps", float(args.ld_wind_snps))
    if args.ld_wind_kb is not None:
        return _LDWindowSpec("kb", float(args.ld_wind_kb))
    if args.ld_wind_cm is not None:
        return _LDWindowSpec("cm", float(args.ld_wind_cm))
    raise LDSCUsageError(
        "ldscore could not validate the R2 panel LD window because no "
        "user LD-window option was set. Most likely argument validation was "
        "bypassed. Specify exactly one of `--ld-wind-snps`, `--ld-wind-kb`, or "
        "`--ld-wind-cm`."
    )


@dataclass(frozen=True)
class R2ParquetContract:
    """Decode, binding, and identity metadata of one validated R2 Parquet file."""

    path: str
    snp_identifier: str
    genome_build: str
    sorted_by_build: str
    n_snps: int
    sidecar_identity_sha256: str
    n_samples: int
    r2_scale: float
    ld_window: _LDWindowSpec
    # Build-time adjusted-R2 floor; pairs below a positive value were not stored.
    min_r2: float


def read_r2_parquet_contract(path) -> R2ParquetContract:
    """Validate one R2 Parquet file footer against the current ``build-r2-panel`` format.

    Readers accept exactly what the writer produces: ``IDX_1``/``IDX_2``/
    ``R2``/``SIGN_R`` columns with an int16 ``R2``, the writer's
    ``ldsc:r2_encoding`` and ``ldsc:r2_scale``, already-unbiased values
    (``ldsc:r2_bias="unbiased"``) with an ``ldsc:n_samples`` of at least 3, and
    complete binding, identity, LD-window, and ``ldsc:min_r2`` metadata.
    Decoding is ``R2 = stored / r2_scale`` with no further correction. No
    value is defaulted; see :func:`read_r2_parquet_contracts` for the error.
    """
    return read_r2_parquet_contracts([path])[0]


def read_r2_parquet_contracts(paths: Sequence) -> list[R2ParquetContract]:
    """Validate several R2 Parquet file footers, reporting every defect in one error.

    Returns one :class:`R2ParquetContract` per path, in order. Every defect of
    every file is collected before a single ``LDSCInputError`` is raised.
    """
    contracts, defects = [], []
    for path in paths:
        contract, problems = _check_r2_parquet_contract(path)
        contracts.append(contract)
        if problems:
            defects.append((path, problems))
    if defects:
        if len(defects) == 1:
            path, problems = defects[0]
            summary = f"R2 Parquet file '{path}' does not match"
            details = "; ".join(problems)
        else:
            summary = f"{len(defects)} R2 Parquet files do not match"
            details = " | ".join(f"'{path}': {'; '.join(problems)}" for path, problems in defects)
        raise LDSCInputError(
            f"{summary} the format written by the current `ldsc build-r2-panel`: {details}. "
            "Most likely the panel was built by an older LDSC3 version or its metadata was "
            "altered. Regenerate the panel with `ldsc build-r2-panel`. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    return contracts


def _check_r2_parquet_contract(path) -> tuple[R2ParquetContract | None, list[str]]:
    """Return the parsed contract of one R2 Parquet file, or ``None`` and its defects."""
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "Reading R2 panels requires pyarrow. Most likely R2 "
            "panel input was requested in an environment missing pyarrow. "
            "Install pyarrow or use PLINK reference-panel input instead."
        ) from exc
    from .ref_panel_builder import R2_ENCODING, R2_QUANT_SCALE

    try:
        schema = pq.read_schema(str(path))
    except (pa.ArrowException, OSError) as exc:
        return None, [f"file cannot be read as parquet ({exc})"]
    raw = {key.decode("utf-8"): value.decode("utf-8", errors="replace")
           for key, value in (schema.metadata or {}).items()}
    problems: list[str] = []

    def text(key, expected=None):
        found = raw.get(f"ldsc:{key}")
        if found is None:
            problems.append(f"ldsc:{key} is missing")
        elif expected is not None and found != expected:
            problems.append(f"ldsc:{key}={found!r} (expected {expected!r})")
        return found

    def parsed(key, parse, valid, expectation):
        found = text(key)
        if found is None:
            return None
        try:
            result = parse(found)
        except (ValueError, LDSCConfigError):
            result = None
        if result is None or not valid(result):
            problems.append(f"ldsc:{key}={found!r} (expected {expectation})")
            return None
        return result

    missing_columns = [name for name in ("IDX_1", "IDX_2", "R2", "SIGN_R") if name not in schema.names]
    if missing_columns:
        problems.append("missing column(s) " + ", ".join(missing_columns))
    elif schema.field("R2").type != pa.int16():
        problems.append(f"R2 column has type {schema.field('R2').type} (expected int16)")
    text("r2_encoding", R2_ENCODING)
    r2_scale = parsed("r2_scale", float, lambda scale: scale == R2_QUANT_SCALE, repr(str(R2_QUANT_SCALE)))
    text("r2_bias", "unbiased")
    # The writer's unbiased correction divides by n - 2, so it needs at least 3 samples.
    n_samples = parsed("n_samples", int, lambda n: n >= 3, "an integer of at least 3")
    n_snps = parsed("n_snps", int, lambda n: n >= 0, "a nonnegative integer")
    identity_hash = text("sidecar_identity_sha256")
    sorted_by_build = parsed("sorted_by_build", normalize_genome_build, lambda build: build in {"hg19", "hg38"},
                             "hg19 or hg38")
    text("artifact_type", "ref_panel_r2")
    snp_identifier = parsed("snp_identifier", normalize_snp_identifier_mode, lambda mode: True,
                            "an LDSC3 SNP identifier mode")
    genome_build = text("genome_build")
    window_mode = parsed("ld_window_mode", str, lambda mode: mode in {"snps", "kb", "cm"}, "snps, kb, or cm")
    window_value = parsed("ld_window_value", float, lambda value: bool(np.isfinite(value)) and value > 0,
                          "a positive number")
    min_r2 = parsed("min_r2", float, lambda value: bool(np.isfinite(value)) and value <= 1,
                    "a number no greater than 1")
    if problems:
        return None, problems
    return R2ParquetContract(
        path=str(path), snp_identifier=snp_identifier, genome_build=genome_build,
        sorted_by_build=sorted_by_build, n_snps=n_snps, sidecar_identity_sha256=identity_hash,
        n_samples=n_samples, r2_scale=r2_scale, ld_window=_LDWindowSpec(window_mode, window_value),
        min_r2=min_r2,
    ), []


def _r2_panel_ld_window_spec(parquet_paths: Sequence[str]) -> _LDWindowSpec | None:
    """Return the common recorded LD-window spec for an R2 panel input set."""
    panel_spec: _LDWindowSpec | None = None
    for path in parquet_paths:
        current = read_r2_parquet_contract(path).ld_window
        if panel_spec is None:
            panel_spec = current
            continue
        if current != panel_spec:
            raise LDSCInputError(
                "ldscore could not validate the R2 panel LD window because "
                "the resolved parquet files record conflicting build windows: "
                f"`{_format_ld_window_spec(panel_spec)}` and "
                f"`{_format_ld_window_spec(current)}`. Most likely files from "
                "different reference-panel builds were mixed in one input. Use a "
                "single coherent `--r2-dir` or regenerate the panel."
            )
    return panel_spec


def _raise_ldscore_window_exceeds_panel_window(
    *,
    requested: _LDWindowSpec,
    panel: _LDWindowSpec,
    chrom: str,
) -> None:
    """Raise the user-facing parquet-window mismatch error."""
    raise LDSCUsageError(
        f"ldscore cannot compute chromosome {chrom} from the R2 "
        f"panel: the user-requested LD window `{_format_ld_window_spec(requested)}` "
        f"is wider than the input R2 panel window "
        f"`{_format_ld_window_spec(panel)}`. Most likely the R2 "
        "panel was built with a smaller LD window, so it does not contain all "
        "SNP pairs required by this ldscore run. Rebuild the R2 panel "
        f"with `{_format_ld_window_spec(requested)}` or a wider build window, or "
        f"rerun ldscore with an LD window no wider than the recorded panel window "
        f"`{_format_ld_window_spec(panel)}`."
    )


def validate_ldscore_window_within_r2_panel_window(
    args: argparse.Namespace,
    *,
    parquet_paths: Sequence[str],
    chrom: str,
) -> None:
    """Reject ldscore windows the R2 panel is not known to cover.

    The rule uses only the recorded parquet LD-window mode and value, so it is
    decided before any computation: a window in another unit is rejected, and
    a same-unit window must be no wider than the panel's. Coverage is not
    inferred from per-SNP coordinates. For ``--ld-wind-snps`` this relies on
    counting panel rows (see :func:`build_window_coordinates`).
    """
    panel_spec = _r2_panel_ld_window_spec(parquet_paths)
    if panel_spec is None:
        return
    requested_spec = _requested_ld_window_spec(args)

    if requested_spec.mode != panel_spec.mode:
        requested, panel = _format_ld_window_spec(requested_spec), _format_ld_window_spec(panel_spec)
        raise LDSCUsageError(
            f"ldscore cannot use the R2 panel for chromosome {chrom} with "
            f"`{requested}`: the panel was built with `{panel}`, and windows in different units "
            "cannot be checked for coverage from their settings, so the run could silently miss "
            f"SNP pairs the panel never stored. Rerun with the panel's window `{panel}` (or a "
            f"narrower `{panel_spec.flag}` value), or rebuild the panel with "
            f"`ldsc build-r2-panel {requested}`."
        )

    if requested_spec.value > panel_spec.value:
        _raise_ldscore_window_exceeds_panel_window(
            requested=requested_spec,
            panel=panel_spec,
            chrom=chrom,
        )


def reject_truncated_r2_panel(parquet_paths: Sequence[str], chrom: str) -> None:
    """Reject an R2 panel that omitted pairs below a positive ``ldsc:min_r2``.

    LD scores sum every SNP pair in the window, so pairs the panel did not
    store would count as zero. Such panels stay valid for ``query-r2``, which
    reports unstored pairs as absent.
    """
    for path in parquet_paths:
        min_r2 = read_r2_parquet_contract(path).min_r2
        if min_r2 > 0:
            raise LDSCInputError(
                f"ldscore cannot use R2 Parquet file '{path}' (chromosome {chrom}): the panel omits SNP pairs with "
                f"unbiased R2 below {min_r2:g} (ldsc:min_r2, set by `build-r2-panel --min-r2`). LD scores sum "
                "every SNP pair in the window, so the omitted pairs would count as zero and change the scores. "
                "Rebuild the panel without `--min-r2` (the default keeps every pair). Panels built with a "
                "positive `--min-r2` remain usable for `ldsc query-r2`."
            )


# Parquet R2 adapter.
@dataclass(frozen=True)
class _DecodedR2RowGroup:
    """Decoded index-format parquet row group as retained-matrix index arrays."""

    row_group_index: int
    i: np.ndarray
    j: np.ndarray
    r2: np.ndarray


class SortedR2BlockReader:
    """
    Stream retained R2 pairs from one chromosome's R2 Parquet file.

    R2 Parquet files use logical fields ``IDX_1``, ``IDX_2``, ``R2``
    (int16 on-disk, dequantized to float32 by dividing by ``ldsc:r2_scale``),
    and ``SIGN_R``. The file must satisfy :func:`read_r2_parquet_contract`; its
    values are already unbiased and are not corrected again. Each row group is
    decoded once and remapped to retained SNP indices.
    """

    def __init__(
        self,
        paths: Sequence[str],
        chrom: str,
        metadata: pd.DataFrame,
        identifier_mode: str,
        genome_build: str | None = None,
    ) -> None:
        """Open one chromosome's sorted R2 Parquet file and build index maps."""
        if not paths:
            raise LDSCInputError(
                f"ldscore could not find a sorted R2 Parquet file for chromosome {chrom}. "
                "Most likely the R2 panel directory is missing this chromosome or the path token "
                "does not resolve to `chrN_r2.parquet`. Pass the correct `--r2-dir` or "
                "regenerate the reference panel. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        self.chrom = normalize_chromosome(chrom)
        self.identifier_mode = normalize_snp_identifier_mode(identifier_mode)
        self.genome_build = genome_build
        assert self.genome_build in {"hg19", "hg38", None}, (
            f"genome_build={self.genome_build!r} must be concrete by this point."
        )
        self._pf = None
        metadata = metadata.copy()
        metadata_context = f"SortedR2BlockReader[{self.chrom}] metadata"
        renamed = {
            resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=metadata_context): "CHR",
            resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=metadata_context): "POS",
            resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=metadata_context): "SNP",
        }
        optional_cm = resolve_optional_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=metadata_context)
        if optional_cm is not None:
            renamed[optional_cm] = "CM"
        metadata = metadata.rename(columns=renamed)

        try:
            import pyarrow.parquet as pq
        except ImportError as exc:
            raise LDSCDependencyError(
                "ldscore could not read R2 panel input because pyarrow is not installed. "
                "Most likely R2 panel input was requested in an environment missing pyarrow. "
                "Install pyarrow or use PLINK reference-panel input instead."
            ) from exc

        self.contract = read_r2_parquet_contract(paths[0])

        if is_allele_aware_mode(self.identifier_mode):
            a1_col = resolve_required_column(metadata.columns, A1_COLUMN_SPEC, context=metadata_context)
            a2_col = resolve_required_column(metadata.columns, A2_COLUMN_SPEC, context=metadata_context)
            metadata = metadata.rename(columns={a1_col: "A1", a2_col: "A2"})
        validate_retained_identifier_uniqueness(metadata, self.identifier_mode, chrom)

        if len(paths) != 1:
            raise LDSCInputError(
                f"ldscore resolved {len(paths)} R2 Parquet files for chromosome {self.chrom}, "
                "but an R2 panel requires exactly one file per chromosome. Most likely "
                "the R2 panel path/glob is too broad or the directory contains duplicate chromosome "
                "artifacts. Narrow `--r2-dir` or remove duplicate parquet files."
            )
        self._pf = pq.ParquetFile(paths[0])
        try:
            self._init_index_path(paths[0], metadata)
        except BaseException:
            self.close()
            raise

    def close(self):
        """Release the parquet handle owned by this chromosome reader."""
        if self._pf is not None:
            self._pf.close()

    @property
    def panel_rows(self) -> np.ndarray:
        """Panel sidecar row of each retained SNP, in retained order."""
        return self._retained_build_idx

    def _init_index_path(self, path: str, retained_metadata: pd.DataFrame) -> None:
        """Validate index-format metadata, the sidecar binding, and build the remap."""
        if self._pf is None:
            raise LDSCInternalError(
                "LD-score parquet reader failed in SortedR2BlockReader._init_index_path(): "
                "the parquet file handle is not initialized. Most likely reader construction "
                "was bypassed. Re-run with DEBUG logging and report the traceback."
            )
        parquet_build = self.contract.sorted_by_build
        if self.genome_build not in {None, parquet_build}:
            raise LDSCInputError(
                f"ldscore could not use R2 Parquet file '{path}': it is sorted for "
                f"{parquet_build}, but the analysis uses {self.genome_build}. Most likely "
                "the reference panel and annotation inputs use different genome builds. "
                "Use the matching reference file or regenerate with "
                f"`--genome-build {self.genome_build}`. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        self.genome_build = parquet_build

        full_sidecar = _load_full_panel_sidecar(path)
        _validate_index_binding(
            full_sidecar, n_snps=self.contract.n_snps, identity_hash=self.contract.sidecar_identity_sha256,
            context=f"SortedR2BlockReader[{self.chrom}] {path}",
        )

        self._remap, self._retained_build_idx = build_index_remap(
            full_sidecar, retained_metadata, self.identifier_mode
        )
        self.m = len(retained_metadata)

        meta = self._pf.metadata
        if meta.num_row_groups > 0:
            avg = meta.num_rows / meta.num_row_groups
            if avg > 500_000:
                LOGGER.warning(
                    f"R2 panel '{path}' stores {meta.num_rows} SNP-pair rows in "
                    f"{meta.num_row_groups} large file blocks (average {avg:.0f} rows per block). "
                    "Reading continues, but queries may be slower and require more memory. "
                    "For repeated use, rebuild the panel with the current `ldsc build-r2-panel`, "
                    "which writes smaller file blocks by default."
                )

    def _decode_index_batch(self, table, row_group_index: int) -> _DecodedR2RowGroup:
        """Dequantize and remap one bounded Arrow batch of stored pairs."""
        idx1 = _arrow_column_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        idx2 = _arrow_column_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        stored = _arrow_column_to_numpy(table.column("R2")).astype(np.float32, copy=False)
        r2 = stored / np.float32(self.contract.r2_scale)
        i = self._remap[idx1]
        j = self._remap[idx2]
        keep = (i >= 0) & (j >= 0)
        return _DecodedR2RowGroup(
            row_group_index=int(row_group_index),
            i=i[keep].astype(np.int32, copy=False),
            j=j[keep].astype(np.int32, copy=False),
            r2=r2[keep].astype(np.float32, copy=False),
        )

    def iter_all_pairs(self, *, batch_rows=65_536):
        """Yield each pair once in bounded batches, including oversized row groups.

        Panel indices are remapped to the retained reference SNPs before
        yielding. SIGN_R is unused. No decoded row group is cached.
        """
        if self._pf is None:
            raise LDSCInternalError(
                "LD-score parquet reader failed in SortedR2BlockReader.iter_all_pairs(): "
                "the parquet file handle is not initialized. Re-run with DEBUG logging "
                "and report the traceback."
            )
        for rg_index in range(self._pf.metadata.num_row_groups):
            for table in self._pf.iter_batches(batch_size=batch_rows,
                    row_groups=[rg_index], columns=["IDX_1", "IDX_2", "R2"]):
                group = self._decode_index_batch(table, rg_index)
                if group.i.size:
                    yield group.i, group.j, group.r2


def ld_score_streaming_from_r2_reader(
    block_left: np.ndarray,
    annot,
    block_reader: SortedR2BlockReader,
    chunk_pairs: int = _CSR_CHUNK_PAIRS,
    *, output_rows=None, n_baseline=None, query_batch_size=100, weight_mask=None,
) -> np.ndarray:
    """Project selected annotation batches through each stored R2 chunk once.

    Only output SNPs own score buffers, but both pair directions retain every
    eligible reference SNP. Diagonals and the separate regression
    weight mask use the same float64 accumulator; direct results are float32.
    """
    block_left = np.asarray(block_left, dtype=np.int64)
    if isinstance(annot, np.ndarray):
        annot = ArrayAnnotations(annot, tuple(str(i) for i in range(annot.shape[1])))
    accumulator = ProjectionAccumulator(annot,
        n_baseline=annot.shape[1] if n_baseline is None else n_baseline,
        output_rows=output_rows, query_batch_size=query_batch_size, weight_mask=weight_mask)
    accumulator.add_diagonal()
    buf_i, buf_j, buf_r2 = [], [], []
    buffered = 0

    def flush():
        nonlocal buffered
        if buffered:
            accumulator.add_pairs(np.concatenate(buf_i), np.concatenate(buf_j),
                                  np.concatenate(buf_r2), block_left)
            buf_i.clear()
            buf_j.clear()
            buf_r2.clear()
            buffered = 0

    for i, j, r2 in block_reader.iter_all_pairs():
        start = 0
        while start < len(i):
            stop = min(len(i), start + chunk_pairs-buffered)
            buf_i.append(i[start:stop])
            buf_j.append(j[start:stop])
            buf_r2.append(r2[start:stop])
            buffered += stop-start
            start = stop
            if buffered == chunk_pairs:
                flush()
    flush()
    return np.asarray(accumulator.values, dtype=np.float32)


def compute_counts(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    common_maf_min: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute LDSC-style annotation count vectors ``M`` and common-SNP counts.

    ``M`` is the column-wise sum over the retained reference SNPs.
    The common-SNP count vector is the same sum restricted to rows with
    ``MAF >= common_maf_min`` (inclusive; deviates from legacy LDSC's strict
    ``0.05 < FRQ < 0.95`` on canonical folded MAF). ``metadata`` must carry
    MAF for every row; both reference-panel adapters guarantee it.
    """
    if "MAF" not in metadata.columns or metadata["MAF"].isna().any():
        raise LDSCInternalError(
            "LD-score counts received reference metadata without MAF for every SNP, so the "
            "common-SNP counts cannot be computed. Most likely a custom reference-panel adapter "
            "dropped the MAF column or returned missing values. Supply MAF for every retained "
            "SNP; for the built-in panels, re-run with DEBUG logging and report the traceback."
        )
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=False)
    M = np.asarray(annot_matrix.sum(axis=0), dtype=np.float64)
    common = metadata["MAF"] >= common_maf_min
    M_5_50 = np.asarray(annot_matrix[common.to_numpy(), :].sum(axis=0), dtype=np.float64)
    return M, M_5_50


def regression_mask_from_keys(
    metadata: pd.DataFrame,
    regression_keys: set[str] | RestrictionIdentityKeys | None,
    identifier_mode: str,
    region_intervals: RegionIntervals | None = None,
) -> np.ndarray:
    """Build the selected-regression-set mask after named region subtraction.

    This is deliberately separate from the reference SNPs: callers use the
    returned mask only for written regression SNPs and the SNPs that contribute to ``w_ld``.
    """
    if regression_keys is None:
        keep = np.ones(len(metadata), dtype=bool)
    elif isinstance(regression_keys, RestrictionIdentityKeys):
        keep = restriction_membership_mask(
            metadata,
            regression_keys,
            identifier_mode,
            context="regression-SNP restriction matching for LD scores",
        ).to_numpy(dtype=bool)
    else:
        keys = identifier_keys(metadata, identifier_mode)
        keep = keys.isin(regression_keys).to_numpy(dtype=bool)
    if region_intervals is not None and region_intervals.intervals:
        keep &= np.asarray(region_exclusion_keep_mask(metadata, region_intervals), dtype=bool)
    return keep.astype(np.float32)


def compute_chromosome(
    chrom: str, prepared: PreparedChromosome, *, snp_identifier: str,
    snp_batch_size: int, common_maf_min: float = 0.05, query_batch_size: int = 100,
    regression_keys: set[str] | RestrictionIdentityKeys | None = None,
    regression_regions: RegionIntervals | None = None,
) -> ChromComputationResult:
    """Project annotations and regression weights through one prepared chromosome.

    Preparation owns identity, sample/MAF selection, reader policy, and window
    validation. This function only projects the aligned matrices and derives
    counts and overlap; it never resolves paths or reopens a reference panel.
    """
    metadata = prepared.metadata
    regression_selected = regression_mask_from_keys(metadata, regression_keys, snp_identifier)
    regression_mask = regression_mask_from_keys(
        metadata, regression_keys, snp_identifier, region_intervals=regression_regions,
    )
    output_rows = np.flatnonzero(regression_mask)
    options = dict(output_rows=output_rows, n_baseline=len(prepared.baseline_columns),
                   query_batch_size=query_batch_size, weight_mask=regression_mask)
    if prepared.backend == "plink":
        prepared.reader._currentSNP = 0
        combined_scores = np.asarray(prepared.reader.ldScoreVarBlocks(
            prepared.block_left, snp_batch_size, annot=prepared.annotations, **options,
        ), dtype=np.float32)
    else:
        combined_scores = ld_score_streaming_from_r2_reader(
            block_left=prepared.block_left, annot=prepared.annotations, block_reader=prepared.reader, **options,
        )
    n_columns = prepared.annotations.shape[1]
    M, M_5_50, overlap, kinds = annotation_statistics(metadata, prepared.annotations,
        len(prepared.baseline_columns), common_maf_min=common_maf_min, query_batch_size=query_batch_size)
    return ChromComputationResult(
        chrom=chrom, metadata=metadata.iloc[output_rows].reset_index(drop=True),
        ld_scores=combined_scores[:, :n_columns], w_ld=combined_scores[:, n_columns:],
        M=M, M_5_50=M_5_50,
        ldscore_columns=prepared.baseline_columns + prepared.query_columns,
        baseline_columns=prepared.baseline_columns, query_columns=prepared.query_columns,
        reference_snp_count=len(metadata), regression_selected_snp_count=int(regression_selected.sum()),
        regression_region_removed_snp_count=int(regression_selected.sum()-len(output_rows)),
        annotation_types=kinds, overlap=overlap, identity_drops=prepared.identity_drops,
    )


def resolve_keep_individuals(keep_path: str | None, fam) -> list[int] | None:
    """Return FAM row indices retained by ``keep_path`` using IID matching."""
    if not keep_path:
        return None
    keep_indivs = fam.loj(legacy_parse.FilterFile(keep_path).IDList)
    if len(keep_indivs) == 0:
        raise LDSCInputError(
            f"ldscore retained no PLINK individuals after applying keep file '{keep_path}'. "
            "Most likely the keep-list IDs do not match the `.fam` file. Check IID/FID "
            "columns in the keep file and pass a keep list from the same PLINK sample set."
        )
    return keep_indivs.tolist()
