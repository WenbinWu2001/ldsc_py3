# Includes primitives derived from upstream LDSC ldscore/ldscore.py.
# GPLv3 (LICENSE); upstream authors and Python 3 port credit are in NOTICE.
# Modified 2026 for prepared reference panels, streaming R2, and projection.
# See NOTICE for attribution and dated modifications.

"""LD-score projection primitives and streaming index-format R2 reader.

``RefPanel.prepare_chromosome`` owns source resolution, reference filtering,
annotation alignment, LD-window validation, and reader lifetime. This module's
``compute_chromosome`` consumes that ``PreparedChromosome`` state and projects
baseline/query annotations and regression weights in one traversal. It also
provides annotation parsing, window geometry, and count primitives used by the
workflow adapters. BED-to-SNP projection belongs to ``ldsc.annotation_builder``.

Scientific contracts
--------------------
Prepared metadata, float32 annotation rows, and window bounds describe the
same retained SNP universe in genomic order. Reference-panel CM and MAF are
authoritative. PLINK sample selection precedes genotype QC and inclusive MAF
filtering; chromosome LD scores retain the existing float32 result conversion.
Parquet projection accumulates in float64. Common-SNP counts and overlap use
``MAF >= common_maf_min`` (0.05 by default), independently of panel filtering.

Canonical parquet R2 inputs contain sidecar-row indices ``IDX_1``/``IDX_2``,
``R2``, and ``SIGN_R``. The required sidecar is bound by its identity digest and
row count. Each chromosome builds one identity remap; pairs stream through
``iter_all_pairs`` without a dense chromosome matrix or row-group cache. Int16
R2 values use scale 32767; floating-point inputs are also accepted. Each
unordered off-diagonal pair contributes in both directions, and the diagonal
is added internally as one. The panel adapter resolves raw/unbiased metadata
and sample size before reader construction. Raw values receive the correction
``r2_raw - (1 - r2_raw) / (n - 2)``; unbiased values receive no second correction.

This module returns in-memory results and writes no artifacts. Public result
normalization and aggregation live in ``ldsc.ldscore_calculator``; canonical
output layout belongs to ``ldsc.outputs``. See
``docs/current/parquet-r2-format-and-read-pipeline.md`` for the file contract.
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
    ANNOTATION_METADATA_SPEC_MAP,
    CHR_COLUMN_ALIASES,
    CHR_COLUMN_SPEC,
    CM_COLUMN_ALIASES,
    CM_COLUMN_SPEC,
    ColumnSpec,
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
from ..path_resolution import (
    resolve_plink_prefix,
)
from .._row_alignment import assert_same_snp_rows
from . import formats as legacy_parse
from .annotation import _annotation_parse_error_message, _validate_annotation_values
from .identifiers import build_snp_id_series, read_snp_restriction_keys
from .overlap import OverlapContribution, compute_overlap, annotation_statistics
from .plink_bed import PlinkBEDFile  # Re-exported for the reference-panel adapters.
from .ldscore_projection import ArrayAnnotations, ProjectionAccumulator
from .snp_identity import (
    RestrictionIdentityKeys,
    empty_identity_drop_frame,
    effective_merge_key_series,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    restriction_membership_mask,
    sidecar_identity_sha256,
)
from .regions import RegionIntervals, region_exclusion_keep_mask


LOGGER = logging.getLogger("LDSC.ldscore")
_LDSCORE_INTERSECTION_DOC = (
    "docs/troubleshooting.md#ldscore-no-annotation-snps-remain-after-reference-panel-intersection"
)
_LDSCORE_PARQUET_DOC = "docs/troubleshooting.md#ldscore-parquet-r2-input-is-incompatible"
# K: stored R2 pairs buffered before one float64 CSR SpMM step (~0.45 GiB/chunk).
_CSR_CHUNK_PAIRS = 16_000_000
REQUIRED_ANNOT_COLUMNS = ("CHR", "POS", "SNP")
ANNOT_META_COLUMNS = ("CHR", "POS", "SNP", "A1", "A2", "CM", "MAF")
ANNOTATION_A1_COLUMN_SPEC = ColumnSpec(
    A1_COLUMN_SPEC.canonical,
    A1_COLUMN_SPEC.aliases,
    A1_COLUMN_SPEC.label,
    allow_suffix_match=False,
)
ANNOTATION_A2_COLUMN_SPEC = ColumnSpec(
    A2_COLUMN_SPEC.canonical,
    A2_COLUMN_SPEC.aliases,
    A2_COLUMN_SPEC.label,
    allow_suffix_match=False,
)
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
    M_5_50: np.ndarray | None
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
    ``block_left`` describes that full contributor universe. Use as a context
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


def block_left_to_right(block_left):
    """Convert block-left coordinates to block-right coordinates."""
    M = len(block_left)
    j = 0
    block_right = np.zeros(M)
    for i in range(M):
        while j < M and block_left[j] <= i:
            j += 1
        block_right[i] = j
    return block_right


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
    """Map panel (build) indices to retained matrix indices for the index format.

    ``full_sidecar`` is the complete panel in build order (the parquet IDX
    space). ``retained_metadata`` is the analysis-restricted matrix universe in
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
            "ldscore could not align retained reference metadata to the parquet R2 sidecar. "
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




def resolve_bfile_prefix(args: argparse.Namespace, chrom: str | None = None) -> str | None:
    """Resolve the PLINK prefix for the requested chromosome, if any."""
    if args.bfile is None:
        return None
    return resolve_plink_prefix(args.bfile, chrom=chrom)


def read_text_table(path: str) -> pd.DataFrame:
    """Read a whitespace-delimited kernel input table with optional gzip compression."""
    compression = "gzip" if path.endswith(".gz") else None
    return pd.read_csv(path, sep=r"\s+", compression=compression, comment="#")


def _panel_sidecar_path_for_r2(r2_path: str) -> Path:
    """Return the ``chrN_meta.tsv.gz`` sidecar path paired with an R2 parquet."""
    p = Path(r2_path)
    if not p.name.endswith("_r2.parquet"):
        raise LDSCInputError(
            f"ldscore could not locate the sidecar for R2 parquet '{r2_path}'. Most likely "
            "the parquet file does not use the canonical `chrN_r2.parquet` filename written "
            "by `ldsc build-r2-panel`. Regenerate the reference panel or pass the canonical "
            "R2 directory. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    return p.with_name(p.name[: -len("_r2.parquet")] + "_meta.tsv.gz")


def _load_full_panel_sidecar(r2_path: str) -> pd.DataFrame:
    """Load the complete (unrestricted) panel sidecar that defines the index space."""
    sidecar_path = _panel_sidecar_path_for_r2(r2_path)
    if not sidecar_path.exists():
        raise LDSCInputError(
            f"ldscore could not load index-format R2 parquet '{r2_path}' because the "
            f"required sidecar '{sidecar_path}' is missing. Most likely the parquet file "
            "was copied without its matching `chrN_meta.tsv.gz` sidecar. Keep the R2 "
            "parquet and sidecar together or regenerate with `ldsc build-r2-panel`. "
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
            f"ldscore could not use parquet R2 input at {context}: the sidecar has "
            f"{len(full_sidecar)} rows but the parquet records n_snps={n_snps}. Most "
            "likely the parquet and sidecar are not a matched pair. Restore the matching "
            "sidecar or regenerate the reference panel. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    actual = sidecar_identity_sha256(full_sidecar, context=context)
    if actual != identity_hash:
        raise LDSCInputError(
            f"ldscore could not use parquet R2 input at {context}: the sidecar identity "
            f"hash {actual} does not match parquet ldsc:sidecar_identity_sha256 "
            f"{identity_hash}. Most likely the sidecar is wrong, reordered, or edited. "
            "Restore the matching sidecar or regenerate the reference panel. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )


def _parquet_schema_layout(schema_names: Sequence[str]) -> str:
    """Classify a runtime parquet schema as index format or unsupported."""
    if {"IDX_1", "IDX_2", "R2", "SIGN_R"}.issubset(set(schema_names)):
        return "index"
    return "unsupported"


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
            f"Chromosome {chrom} has duplicate retained SNP IDs. "
            "This is ambiguous in base rsID mode. Most likely the retained reference "
            "panel contains duplicate rsIDs. Use an allele-aware rsID mode or remove "
            "duplicate retained SNP IDs."
        )


# Annotation loading and normalization.
def parse_annotation_file(
    path: str,
    chrom: str | None = None,
    identifier_mode: str = "rsid",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse one SNP-level annotation table into normalized metadata and values.

    When ``chrom`` is provided, the returned tables are restricted to rows whose
    normalized ``CHR`` matches that chromosome. This is the final safeguard for
    group-style path tokens that resolved to shared multi-chromosome files or
    to globs whose filenames did not encode chromosome labels clearly enough for
    earlier path-level filtering. For LD-score calculation, annotation ``CM``
    is discarded and the reference panel supplies genetic-map metadata.
    """
    try:
        df = read_text_table(path)
    except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
        raise LDSCInputError(_annotation_parse_error_message(path, details=str(exc))) from exc
    context = path
    chr_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CHR"], context=context)
    pos_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["POS"], context=context)
    if identity_mode_family(identifier_mode) == "rsid":
        snp_col = resolve_required_column(
            df.columns, ANNOTATION_METADATA_SPEC_MAP["SNP"], context=context
        )
    else:
        snp_col = resolve_optional_column(
            df.columns, ANNOTATION_METADATA_SPEC_MAP["SNP"], context=context
        )
    # CM/MAF are population-specific; the reference panel is authoritative. They are
    # resolved only to exclude them from annotation value columns, never used as values.
    cm_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CM"], context=context)
    maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)
    a1_col = resolve_optional_column(df.columns, ANNOTATION_A1_COLUMN_SPEC, context=context)
    a2_col = resolve_optional_column(df.columns, ANNOTATION_A2_COLUMN_SPEC, context=context)
    if (a1_col is None) ^ (a2_col is None):
        raise LDSCInputError(
            f"ldscore could not parse annotation file '{path}': it has only one allele "
            "column. Most likely the file contains A1 without A2, or A2 without A1. "
            "Provide both allele columns, or remove both to match annotations by base SNP identity. "
            "Allele-free annotations are also supported in allele-aware modes."
        )

    metadata_columns = {"CHR": df[chr_col], "POS": df[pos_col]}
    if snp_col is not None:
        metadata_columns["SNP"] = df[snp_col]
    meta = pd.DataFrame(metadata_columns)
    meta["CHR"] = meta["CHR"].map(lambda value: normalize_chromosome(value, context=path))
    meta["POS"] = pd.to_numeric(meta["POS"], errors="raise").astype(np.int64)
    if "SNP" in meta.columns:
        meta["SNP"] = meta["SNP"].astype(str)
    # CM is a population-agnostic placeholder (NaN); any input value is discarded.
    meta["CM"] = np.nan
    if a1_col is not None and a2_col is not None:
        meta["A1"] = df[a1_col]
        meta["A2"] = df[a2_col]

    if chrom is not None:
        keep = meta["CHR"] == normalize_chromosome(chrom, context=path)
        meta = meta.loc[keep].reset_index(drop=True)
        df = df.loc[keep].reset_index(drop=True)
    else:
        meta = meta.reset_index(drop=True)

    if len(meta) == 0:
        return meta, pd.DataFrame(index=meta.index)

    metadata_source_columns = {chr_col, pos_col, snp_col, cm_col, maf_col, a1_col, a2_col}
    annotation_columns = [col for col in df.columns if col not in metadata_source_columns]
    if not annotation_columns:
        raise LDSCInputError(
            f"ldscore could not parse annotation file '{path}': no annotation value "
            "columns remain after metadata columns were removed. Most likely the file "
            "contains only CHR/POS/SNP/CM metadata, or the annotation columns were named "
            "as metadata aliases. Add at least one annotation column with numeric values."
        )
    annotation_values = _validate_annotation_values(df, annotation_columns, path=path)

    # Sort metadata and annotation values together through the one shared genomic
    # sort, then split the result. A single permutation keeps the two tables
    # aligned by construction -- there is no second, hand-rolled sort that could
    # silently diverge from the metadata order.
    combined = pd.concat([meta, annotation_values.reset_index(drop=True)], axis=1)
    combined = sort_frame_by_genomic_position(combined)
    annotations = combined.loc[:, annotation_columns].reset_index(drop=True)
    meta = combined.drop(columns=annotation_columns).reset_index(drop=True)
    return meta, annotations


def combine_annotation_groups(
    baseline_files: Sequence[str],
    query_files: Sequence[str],
    chrom: str,
    identifier_mode: str,
) -> AnnotationBundle | None:
    """
    Read and column-bind baseline/query annotation files for one chromosome.

    Main steps:
    1. Parse each SNP-level annotation file and restrict it to ``chrom``.
    2. Normalize the chosen SNP identifier key.
    3. Require identical SNP rows across all files for that chromosome.
    4. Merge annotation columns into one dense matrix while preserving the
       baseline/query grouping.

    "Identical SNP rows" means that every contributing file for ``chrom`` must
    produce the same ordered SNP universe after chromosome filtering and
    sorting. In practice, the ``CHR``, ``POS``, and ``SNP`` columns must match
    row-for-row, and the identifier key derived from ``identifier_mode`` must
    also match row-for-row. If two files contribute columns for the same
    chromosome but disagree on their retained SNP rows, the function raises a
    ``ValueError`` instead of attempting a partial merge.
    """
    frames: list[pd.DataFrame] = []
    blocks: list[pd.DataFrame] = []
    baseline_columns: list[str] = []
    query_columns: list[str] = []
    seen_columns: set[str] = set()

    for group_name, files in (("baseline", baseline_files), ("query", query_files)):
        for path in files:
            meta, annotations = parse_annotation_file(
                path, chrom=chrom, identifier_mode=identifier_mode
            )
            if len(meta) == 0:
                continue
            meta = meta.copy()
            meta["_key"] = identifier_keys(meta, _annotation_available_precision_mode(meta, identifier_mode))
            if not frames:
                frames.append(meta)
            else:
                reference = frames[0]
                alignment_mode = _annotation_pair_alignment_mode(reference, meta, identifier_mode)
                assert_same_snp_rows(
                    reference,
                    meta,
                    context=f"Annotation SNP rows do not match across files for chromosome {chrom}: {path}",
                    snp_identifier=alignment_mode,
                )
                if {"A1", "A2"}.issubset(meta.columns) and not {"A1", "A2"}.issubset(frames[0].columns):
                    frames[0]["A1"] = meta["A1"].to_numpy()
                    frames[0]["A2"] = meta["A2"].to_numpy()
                missing_cm = reference["CM"].isna() & meta["CM"].notna()
                if missing_cm.any():
                    frames[0].loc[missing_cm, "CM"] = meta.loc[missing_cm, "CM"].to_numpy()
                if "MAF" in meta.columns:
                    if "MAF" not in frames[0].columns:
                        frames[0]["MAF"] = np.nan
                    missing_maf = frames[0]["MAF"].isna() & meta["MAF"].notna()
                    if missing_maf.any():
                        frames[0].loc[missing_maf, "MAF"] = meta.loc[missing_maf, "MAF"].to_numpy()

            for column in annotations.columns:
                if column in seen_columns:
                    raise LDSCInputError(
                        f"ldscore could not combine annotation files for chromosome {chrom}: "
                        f"duplicate annotation column '{column}' was found. Most likely two "
                        "baseline/query files define the same annotation name. Rename one "
                        "column or remove the duplicate annotation input."
                    )
                seen_columns.add(column)

            blocks.append(annotations.reset_index(drop=True))
            if group_name == "baseline":
                baseline_columns.extend(annotations.columns.tolist())
            else:
                query_columns.extend(annotations.columns.tolist())

    if not frames:
        return None

    metadata = frames[0].drop(columns="_key").reset_index(drop=True)
    annotations = pd.concat(blocks, axis=1)
    ordered_columns = baseline_columns + query_columns
    annotations = annotations.loc[:, ordered_columns].reset_index(drop=True)
    return AnnotationBundle(
        metadata=metadata,
        annotations=annotations,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
    )


def _annotation_available_precision_mode(metadata: pd.DataFrame, identifier_mode: str) -> str:
    """Return allele-aware annotation identity only when metadata has A1/A2."""
    mode = normalize_snp_identifier_mode(identifier_mode)
    if is_allele_aware_mode(mode) and not {"A1", "A2"}.issubset(metadata.columns):
        return identity_base_mode(mode)
    return mode


def _annotation_pair_alignment_mode(left: pd.DataFrame, right: pd.DataFrame, identifier_mode: str) -> str:
    """Return the row-alignment mode supported by both annotation tables."""
    mode = normalize_snp_identifier_mode(identifier_mode)
    if is_allele_aware_mode(mode) and not (
        {"A1", "A2"}.issubset(left.columns) and {"A1", "A2"}.issubset(right.columns)
    ):
        return identity_base_mode(mode)
    return mode


def read_identifier_list(path: str, mode: str) -> RestrictionIdentityKeys:
    """Read a SNP list file into base or allele-aware restriction keys for ``mode``."""
    return read_snp_restriction_keys(path, normalize_snp_identifier_mode(mode))


def load_regression_keys(args: argparse.Namespace) -> RestrictionIdentityKeys | None:
    """Load the optional regression SNP universe from CLI arguments."""
    if not getattr(args, "regr_snps_file", None):
        return None
    return read_identifier_list(args.regr_snps_file, args.snp_identifier)




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
            "See docs/troubleshooting.md#ldscore-unusable-cm-for-ld-wind-cm."
        )


def require_reference_maf(metadata: pd.DataFrame, chrom: str) -> None:
    """Raise when the reference panel supplies no usable MAF for a chromosome.

    MAF is mandatory: M_5_50 common-SNP counts, the partitioned-h2 common-overlap
    correction, and any ``--maf-min`` filter are meaningless without it.
    """
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        raise LDSCInputError(
            f"ldscore requires MAF for the reference panel but none is available on "
            f"chromosome {chrom}. Most likely a parquet panel was built without allele "
            "frequencies (sidecar MAF=NA). Rebuild the reference panel with MAF, or use a "
            "PLINK panel (MAF is computed from genotypes)."
        )


def build_window_coordinates(metadata: pd.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray, float]:
    """Build the coordinate array and maximum distance for the active LD window mode."""
    selectors = np.array([args.ld_wind_snps is not None, args.ld_wind_kb is not None, args.ld_wind_cm is not None], dtype=bool)
    if selectors.sum() != 1:
        raise LDSCUsageError(
            "ldscore could not choose an LD-window mode. Most likely zero or multiple "
            "LD-window options were supplied. Specify exactly one of `--ld-wind-snps`, "
            "`--ld-wind-kb`, or `--ld-wind-cm`."
        )

    if args.ld_wind_snps is not None:
        return np.arange(len(metadata), dtype=float), float(args.ld_wind_snps)
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
    return f"{spec.flag} {spec.value}"


def _requested_ld_window_spec(args: argparse.Namespace) -> _LDWindowSpec:
    """Return the user-requested LD-window mode and value from parsed args."""
    if args.ld_wind_snps is not None:
        return _LDWindowSpec("snps", float(args.ld_wind_snps))
    if args.ld_wind_kb is not None:
        return _LDWindowSpec("kb", float(args.ld_wind_kb))
    if args.ld_wind_cm is not None:
        return _LDWindowSpec("cm", float(args.ld_wind_cm))
    raise LDSCUsageError(
        "ldscore could not validate the parquet R2 panel LD window because no "
        "user LD-window option was set. Most likely argument validation was "
        "bypassed. Specify exactly one of `--ld-wind-snps`, `--ld-wind-kb`, or "
        "`--ld-wind-cm`."
    )


def _read_r2_panel_ld_window_spec(path: str) -> _LDWindowSpec | None:
    """Read build-time LD-window metadata from one canonical R2 parquet."""
    try:
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "ldscore could not validate the parquet R2 panel LD window because "
            "pyarrow is not installed. Most likely parquet reference-panel mode "
            "was requested in an environment missing pyarrow. Install pyarrow or "
            "use PLINK reference-panel input instead."
        ) from exc

    schema_meta = pq.read_schema(str(path)).metadata or {}
    mode_raw = schema_meta.get(b"ldsc:ld_window_mode")
    value_raw = schema_meta.get(b"ldsc:ld_window_value")
    if mode_raw is None and value_raw is None:
        return None
    if mode_raw is None or value_raw is None:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            "the parquet has incomplete ldsc:ld_window_* schema metadata. Most "
            "likely the artifact was edited or written by an incompatible "
            "development version. Regenerate the reference panel with "
            "`ldsc build-r2-panel`."
        )
    mode = mode_raw.decode("utf-8")
    if mode not in {"snps", "kb", "cm"}:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            f"unsupported parquet ldsc:ld_window_mode={mode!r}. Most likely the "
            "artifact was edited or written by an incompatible development "
            "version. Regenerate the reference panel with `ldsc build-r2-panel`."
        )
    try:
        value = float(value_raw.decode("utf-8"))
    except ValueError as exc:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            f"parquet ldsc:ld_window_value={value_raw.decode('utf-8', errors='replace')!r} "
            "is not numeric. Most likely the artifact schema metadata was edited "
            "or corrupted. Regenerate the reference panel with `ldsc build-r2-panel`."
        ) from exc
    if value <= 0:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            f"parquet ldsc:ld_window_value={value} is not positive. Most likely "
            "the artifact schema metadata was edited or corrupted. Regenerate the "
            "reference panel with `ldsc build-r2-panel`."
        )
    return _LDWindowSpec(mode, value)


def _r2_panel_ld_window_spec(parquet_paths: Sequence[str]) -> _LDWindowSpec | None:
    """Return the common recorded LD-window spec for a parquet R2 input set."""
    panel_spec: _LDWindowSpec | None = None
    for path in parquet_paths:
        current = _read_r2_panel_ld_window_spec(path)
        if current is None:
            continue
        if panel_spec is None:
            panel_spec = current
            continue
        if current != panel_spec:
            raise LDSCInputError(
                "ldscore could not validate the parquet R2 panel LD window because "
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
        f"ldscore cannot compute chromosome {chrom} from the R2 parquet reference "
        f"panel: the user-requested LD window `{_format_ld_window_spec(requested)}` "
        f"is wider than the input R2 parquet panel window "
        f"`{_format_ld_window_spec(panel)}`. Most likely the parquet reference "
        "panel was built with a smaller LD window, so it does not contain all "
        "SNP pairs required by this ldscore run. Rebuild the R2 parquet panel "
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
    """Reject ldscore windows wider than the recorded R2 parquet build window.

    Old parquet artifacts without ``ldsc:ld_window_*`` metadata cannot be checked
    reliably, so they retain their previous behavior. The check is deliberately
    limited to the recorded parquet LD-window mode and value; it does not infer
    pair coverage from per-SNP coordinates.
    """
    panel_spec = _r2_panel_ld_window_spec(parquet_paths)
    if panel_spec is None:
        return
    requested_spec = _requested_ld_window_spec(args)

    if requested_spec.mode != panel_spec.mode:
        LOGGER.warning(
            f"ldscore cannot directly compare the requested LD window to the R2 parquet "
            f"panel window because their modes differ: the R2 parquet panel was built "
            f"with `{_format_ld_window_spec(panel_spec)}`, but this ldscore run requested "
            f"`{_format_ld_window_spec(requested_spec)}`. Most likely the run is mixing "
            "a reference panel built under one LD-window coordinate system with an "
            "ldscore request in another. Continuing, but the requested window may require "
            "SNP pairs not stored in the parquet, so LD scores may miss some SNP-pair "
            "contributions. To avoid this risk, rebuild the R2 parquet panel with the "
            "same LD-window mode as this ldscore run, or rerun ldscore with the panel's "
            f"recorded mode `{panel_spec.flag}`."
        )
        return

    if requested_spec.value > panel_spec.value:
        _raise_ldscore_window_exceeds_panel_window(
            requested=requested_spec,
            panel=panel_spec,
            chrom=chrom,
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
    Stream retained R2 pairs from a chromosome-scoped index parquet table.

    Index-format parquet files use logical fields ``IDX_1``, ``IDX_2``, ``R2``
    (int16 on-disk, dequantized to float32 by dividing by ``ldsc:r2_scale``),
    and ``SIGN_R``. Each row group is decoded once and remapped to retained SNP
    indices. Float32 ``R2`` columns (absent ``ldsc:r2_encoding`` metadata)
    are read unscaled for backward compatibility.
    """

    def __init__(
        self,
        paths: Sequence[str],
        chrom: str,
        metadata: pd.DataFrame,
        identifier_mode: str,
        r2_bias_mode: str,
        r2_sample_size: float | None,
        genome_build: str | None = None,
    ) -> None:
        """Open one chromosome's sorted parquet R2 tables and build index maps."""
        if not paths:
            raise LDSCInputError(
                f"ldscore could not find sorted parquet R2 input for chromosome {chrom}. "
                "Most likely the R2 directory is missing this chromosome or the path token "
                "does not resolve to `chrN_r2.parquet`. Pass the correct `--r2-dir` or "
                "regenerate the reference panel. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        self.chrom = normalize_chromosome(chrom)
        self.identifier_mode = normalize_snp_identifier_mode(identifier_mode)
        self.r2_bias_mode = r2_bias_mode
        self.r2_sample_size = r2_sample_size
        self.genome_build = genome_build
        assert self.genome_build in {"hg19", "hg38", None}, (
            f"genome_build={self.genome_build!r} must be concrete by this point."
        )
        self._pf = None
        self._r2_scale: float | None = None
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
                "ldscore could not read index parquet R2 input because pyarrow is not installed. "
                "Most likely parquet reference-panel mode was requested in an environment missing pyarrow. "
                "Install pyarrow or use PLINK reference-panel input instead."
            ) from exc

        layout = _parquet_schema_layout(pq.read_schema(paths[0]).names)
        if layout != "index":
            raise LDSCInputError(
                f"ldscore could not use R2 parquet '{paths[0]}': it is not an index-format "
                "R2 parquet with columns IDX_1/IDX_2/R2/SIGN_R. Most likely this file was "
                "written by an old LDSC version or is not an LDSC R2 artifact. Regenerate "
                "the reference panel with `ldsc build-r2-panel`. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )

        if is_allele_aware_mode(self.identifier_mode):
            a1_col = resolve_required_column(metadata.columns, A1_COLUMN_SPEC, context=metadata_context)
            a2_col = resolve_required_column(metadata.columns, A2_COLUMN_SPEC, context=metadata_context)
            metadata = metadata.rename(columns={a1_col: "A1", a2_col: "A2"})
        validate_retained_identifier_uniqueness(metadata, self.identifier_mode, chrom)

        if len(paths) != 1:
            raise LDSCInputError(
                f"ldscore resolved {len(paths)} parquet R2 files for chromosome {self.chrom}, "
                "but index parquet mode requires exactly one file per chromosome. Most likely "
                "the R2 path/glob is too broad or the directory contains duplicate chromosome "
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

    def _init_index_path(self, path: str, retained_metadata: pd.DataFrame) -> None:
        """Validate index-format metadata, the sidecar binding, and build the remap."""
        if self._pf is None:
            raise LDSCInternalError(
                "LD-score parquet reader failed in SortedR2BlockReader._init_index_path(): "
                "the parquet file handle is not initialized. Most likely reader construction "
                "was bypassed. Re-run with DEBUG logging and report the traceback."
            )
        schema_meta = self._pf.schema_arrow.metadata or {}
        self._r2_scale = self._resolve_r2_scale(schema_meta, path)

        build_raw = schema_meta.get(b"ldsc:sorted_by_build")
        if build_raw is None:
            raise LDSCInputError(
                f"ldscore could not use index-format R2 parquet '{path}': it has no "
                "ldsc:sorted_by_build metadata. Most likely the artifact was written by "
                "an old LDSC version or had schema metadata stripped. Regenerate with "
                "`ldsc build-r2-panel`. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        parquet_build = normalize_genome_build(build_raw.decode("utf-8"))
        if self.genome_build not in {None, parquet_build}:
            raise LDSCInputError(
                f"ldscore could not use R2 parquet '{path}': it is sorted for "
                f"{parquet_build}, but the analysis uses {self.genome_build}. Most likely "
                "the reference panel and annotation inputs use different genome builds. "
                "Use the matching reference file or regenerate with "
                f"`--genome-build {self.genome_build}`. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        self.genome_build = parquet_build

        n_snps_raw = schema_meta.get(b"ldsc:n_snps")
        hash_raw = schema_meta.get(b"ldsc:sidecar_identity_sha256")
        if n_snps_raw is None or hash_raw is None:
            raise LDSCInputError(
                f"ldscore could not use R2 parquet '{path}': it is missing ldsc:n_snps "
                "or ldsc:sidecar_identity_sha256 binding metadata. Most likely the file "
                "was written by an old LDSC version or metadata was stripped. Regenerate "
                "with `ldsc build-r2-panel`. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        n_snps = int(n_snps_raw.decode("utf-8"))

        full_sidecar = _load_full_panel_sidecar(path)
        _validate_index_binding(
            full_sidecar, n_snps=n_snps, identity_hash=hash_raw.decode("utf-8"),
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
                    f"'{path}' has {meta.num_row_groups} row group(s) (avg {avg:.0f} rows/group); "
                    "query performance will be degraded. Regenerate with row_group_size=50000."
                )

    def _resolve_r2_scale(self, schema_meta: dict, path: str) -> float | None:
        """Return the dequant scale when R2 is stored as quantized integers.

        Detection is by on-disk column dtype: an integer ``R2`` column is
        quantized and the scale comes from ``ldsc:r2_scale`` (defaulting to
        32767 with a warning if the key is absent). A float ``R2`` column is the
        legacy/unquantized path and returns ``None``.
        """
        import pyarrow as pa

        if self._pf is None or not pa.types.is_integer(self._pf.schema_arrow.field("R2").type):
            return None
        scale_raw = schema_meta.get(b"ldsc:r2_scale")
        if scale_raw is None:
            LOGGER.warning(
                f"'{path}' has an integer R2 column but no ldsc:r2_scale; defaulting to 32767."
            )
            return 32767.0
        return float(scale_raw.decode("utf-8"))

    def _transform_r2(self, values: np.ndarray) -> np.ndarray:
        """Apply the configured raw-to-unbiased R2 correction when required."""
        values = values.astype(np.float32, copy=False)
        if self.r2_bias_mode == "raw":
            if self.r2_sample_size is None:
                raise LDSCUsageError(
                    "ldscore cannot apply raw R2 bias correction without a sample size. "
                    "Most likely the parquet panel declares `ldsc:r2_bias=raw` but omits "
                    "`ldsc:n_samples`. Rebuild the panel with the current build-r2-panel, "
                    "or record `ldsc:n_samples` (and `ldsc:r2_bias`) in the parquet metadata."
                )
            denom = self.r2_sample_size - 2
            if denom <= 0:
                raise LDSCConfigError(
                    f"ldscore received invalid R2 sample size ({self.r2_sample_size}) "
                    "for raw R2 correction. Most likely the panel's `ldsc:n_samples` is too "
                    "small for the unbiased correction denominator. Use a panel whose sample "
                    "size exceeds 2, or one that stores pre-corrected (`ldsc:r2_bias=unbiased`) R2."
                )
            values = values - (1.0 - values) / denom
            # Share the writer's R2<=1 invariant: roundoff/raw inputs can exceed 1.
            values = np.minimum(values, np.float32(1.0))
        return values

    def _decode_index_batch(self, table, row_group_index: int) -> _DecodedR2RowGroup:
        """Dequantize and remap one bounded Arrow batch of stored pairs."""
        idx1 = _arrow_column_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        idx2 = _arrow_column_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        r2_raw = _arrow_column_to_numpy(table.column("R2")).astype(np.float32, copy=False)
        if self._r2_scale is not None:
            r2_raw = r2_raw / np.float32(self._r2_scale)
        r2 = self._transform_r2(r2_raw)
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

        Panel indices are remapped to the retained reference universe before
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


def _accumulate_pair_contributions(
    cor_sum: np.ndarray,
    i: np.ndarray,
    j: np.ndarray,
    r2: np.ndarray,
    annot: np.ndarray,
    block_left: np.ndarray,
) -> None:
    """Add one chunk of stored pairs ``(i<j, r2)`` to ``cor_sum`` via a float64 CSR SpMM.

    Builds the strict-upper-triangular sparse matrix ``U`` of the window-filtered
    chunk and accumulates ``cor_sum += U @ annot + U.T @ annot`` (R is symmetric).
    The SpMM is float64 (``U.data`` is float64 and the driver passes a float64
    ``annot``): scipy accumulates each row sum in the operand dtype, so a float32
    SpMM would lose ~3e-3 (see ``lessons.md``). Pairs outside the ldscore window are
    dropped via ``i >= block_left[j]``; the CSR sums repeated row indices natively.
    """
    from scipy import sparse

    keep = i >= block_left[j]
    if not keep.all():
        i, j, r2 = i[keep], j[keep], r2[keep]
    if i.size == 0:
        return
    m = cor_sum.shape[0]
    u = sparse.csr_matrix((r2.astype(np.float64, copy=False), (i, j)), shape=(m, m))
    cor_sum += u @ annot
    cor_sum += u.T @ annot


def ld_score_streaming_from_r2_reader(
    block_left: np.ndarray,
    annot,
    block_reader: SortedR2BlockReader,
    chunk_pairs: int = _CSR_CHUNK_PAIRS,
    *, output_rows=None, n_baseline=None, query_batch_size=1000, weight_mask=None,
) -> np.ndarray:
    """Project selected annotation batches through each stored R² chunk once.

    Only output SNPs own score buffers, but both pair directions retain every
    eligible reference contributor. Diagonals and the separate regression
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
) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Compute LDSC-style annotation count vectors ``M`` and optional common counts.

    ``M`` is the column-wise sum over the retained reference SNP universe.
    The common-count vector is the same sum restricted to rows with
    ``MAF >= common_maf_min`` (inclusive; deviates from legacy LDSC's strict
    ``0.05 < FRQ < 0.95`` on canonical folded MAF) when MAF metadata is
    available.
    """
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=False)
    M = np.asarray(annot_matrix.sum(axis=0), dtype=np.float64)
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        return M, None
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

    This is deliberately separate from the LD-reference universe: callers use
    the returned mask only for written regression rows and ``w_ld`` contributors.
    """
    if regression_keys is None:
        keep = np.ones(len(metadata), dtype=bool)
    elif isinstance(regression_keys, RestrictionIdentityKeys):
        keep = restriction_membership_mask(
            metadata,
            regression_keys,
            identifier_mode,
            context="LD-score regression SNP restriction matching",
        ).to_numpy(dtype=bool)
    else:
        keys = identifier_keys(metadata, identifier_mode)
        keep = keys.isin(regression_keys).to_numpy(dtype=bool)
    if region_intervals is not None and region_intervals.intervals:
        keep &= np.asarray(region_exclusion_keep_mask(metadata, region_intervals), dtype=bool)
    return keep.astype(np.float32)


def compute_chromosome(
    chrom: str, prepared: PreparedChromosome, *, snp_identifier: str,
    snp_batch_size: int, common_maf_min: float = 0.05, query_batch_size: int = 1000,
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
