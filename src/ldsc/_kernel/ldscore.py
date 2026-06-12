"""
Description
-----------
New LD-score estimation core for SNP-level annotation files with support for
either PLINK genotypes or sorted per-chromosome parquet `R2` tables as the
reference panel.

This module consumes SNP-level `.annot` / `.annot.gz` inputs only. It does not
perform BED-to-SNP projection; that remains the responsibility of
`utils/BED2annot.py`.

Supported Inputs
----------------
- Query and baseline annotation inputs are provided separately. Each file must
  contain SNP-level rows with required metadata columns `CHR POS SNP CM`, plus
  one or more annotation columns.
- Reference-panel input may be either:
  - PLINK genotype files (`.bed/.bim/.fam`), or
  - one per-chromosome sorted parquet `R2` file.
- Optional frequency metadata may be provided to populate `MAF` in outputs and
  to enable common-count generation when parquet input is used. For LD-score
  calculation, the annotation file's `CM` is the first source. The sidecar
  metadata only fills missing `CM` values; duplicate sidecar identity clusters
  are dropped entirely before filling.
- Optional regression SNP lists may be provided to define the SNP set used for
  `w_ld`. If omitted, the retained reference SNP set is used.

Index Parquet `R2` Format
--------------------------
- Package-written parquet R2 files contain exactly four columns: `IDX_1`,
  `IDX_2` (int32 sidecar-row indices), `R2` (int16, symmetric quantization scale
  32767, dequantized to float32 on read), and `SIGN` (bit-packed bool, ``True``
  for Pearson r >= 0). They carry no SNP identity; identity lives once-per-SNP in
  the paired `chrN_meta.tsv.gz` sidecar. See
  ``docs/current/parquet-r2-format-and-read-pipeline.md``.
- The sort build is recorded under `ldsc:sorted_by_build`; the runtime query
  build must match it. The parquet is bound to its sidecar by `ldsc:n_snps` and
  `ldsc:sidecar_identity_sha256`; the sidecar is mandatory.
- Rows are sorted by non-decreasing `IDX_1`. Footer statistics for `IDX_1` form
  a row-group index so each window reads only overlapping row groups; decoded
  row groups are cached across overlapping sliding-window queries.
- The same index parquet serves all four identifier modes: the reader resolves
  each panel SNP to a retained matrix index once per chromosome (the remap), and
  per-pair decode is a pure integer gather. Legacy 10-column and external raw
  schemas are not supported; such files raise a regenerate error.
- In `chr_pos`-family modes, retained reference SNP rows must have unique
  chromosome positions after active identity cleanup. If two retained SNPs
  share the same `CHR` and `POS` in base `chr_pos`, matching to the `R2` table
  is ambiguous and the code fails fast; allele-aware multi-allelic base-key
  clusters are removed before matching.

Identifier and Genome-Build Rules
---------------------------------
- Supported SNP identifier modes are `rsid`, `rsid_allele_aware`, `chr_pos`,
  and `chr_pos_allele_aware`.
- In chr-pos-family modes, retained reference SNP rows are matched by
  chromosome and position, plus A1/A2 identity when the active mode and input
  schema require it.
- Canonical parquet files are build-specific. If `ldsc:sorted_by_build`
  conflicts with the selected `genome_build`, reader initialization raises a
  `ValueError`. If the metadata key is absent, the reader infers the build from
  the first row group and logs a warning.

Outputs
-------
This module is the internal compute kernel. Its historical writer can still
materialize legacy prefix-based files for compatibility tests, but the public
``ldsc ldscore`` workflow wraps kernel results in ``LDScoreResult`` and writes a
canonical directory:

- ``metadata.json``
- ``ldscore.baseline.parquet``, containing ``CHR``, ``POS``, ``SNP``, ``regression_ld_scores``,
  and baseline LD-score columns
- optional ``ldscore.query.parquet``, containing ``CHR``, ``POS``, ``SNP``, and query
  LD-score columns

Count records are stored in root metadata rather than as public ``.M`` sidecar
files. Outputs retain ``CM`` and ``MAF`` internally when available, with missing
values represented as ``NA`` by the legacy serializers.

Example Usage
-------------
Python
    from ldsc import run_ldscore

    run_ldscore(
        output_dir="results/example_ldscores",
        baseline_annot_sources="data/baseline.@.annot.gz",
        query_annot_sources="data/query.annot.gz",
        r2_dir="data/ref_panel/hg38",
        r2_bias_mode="unbiased",
        ld_wind_cm=1,
    )

CLI
    python -m ldsc ldscore \
        --output-dir results/example_ldscores \
        --baseline-annot-sources data/baseline.@.annot.gz \
        --query-annot-sources data/query.annot.gz \
        --r2-dir data/ref_panel/hg38 \
        --snp-identifier rsid \
        --r2-bias-mode unbiased \
        --ld-wind-cm 1

Computation Overview
--------------------
1. Load query and baseline SNP-level annotation files chromosome by chromosome.
2. Normalize SNP identity using the selected identifier mode and require all
   annotation files for a chromosome to have identical retained SNP rows.
3. Build the canonical reference SNP table from annotation inputs.
4. Intersect annotation SNPs with the reference-panel SNP universe. Metadata
   sidecars define the full universe when present; otherwise the workflow falls
   back to the partial SNP universe recoverable from parquet pair endpoints.
5. Compute LD scores per chromosome:
   - PLINK mode reuses the original genotype-block implementation.
   - parquet mode follows the old LDSC sliding-block accumulation workflow and
     fills block-local dense `R2` matrices from the sorted parquet file.
6. Compute all partitioned LD-score columns and `w_ld` from the same block
   traversal.
7. Aggregate chromosome results; the public wrapper writes the canonical
   LD-score result directory.

Raw-vs-Unbiased `R2` Handling
-----------------------------
- parquet input requires the user to declare whether `R2` is `raw` sample
  `r^2` or already `unbiased`.
- If `R2` is raw, off-diagonal entries are transformed using the original LDSC
  correction:

      r2_unbiased = r2_raw - (1 - r2_raw) / (n - 2)

  where `n` is the LD reference sample size.
- If `R2` is already unbiased, it is used as-is.
- Each unordered SNP pair contributes to two symmetric matrix entries in the
  block-local dense matrices.
- The diagonal is always added internally as `r_jj^2 = 1`.

Window Rules
------------
- `--ld-wind-snps`, `--ld-wind-kb`, and `--ld-wind-cm` follow the same semantics
  as the original codebase.
- In parquet mode, `block_left` and the old sliding-block iteration determine
  the query bounds and accumulation order.
- `--ld-wind-cm` requires non-missing CM coordinates for retained SNPs. Those
  coordinates come first from annotation metadata, with sidecar metadata filling
  only missing `CM` values.

MAF and Common-Count Rules
--------------------------
- In PLINK mode, MAF is taken from the genotype data after SNP filtering.
- In parquet mode, MAF is only available when supplied via an optional
  frequency/metadata file.
- If MAF is unavailable, `MAF` is written as `NA` in LD-score outputs and
  common counts are not emitted.

- Runtime parquet input must be the canonical 4-column index format
  (``IDX_1``, ``IDX_2``, ``R2``, ``SIGN``), paired with its metadata sidecar.
  The same parquet serves all four identifier modes; mode-specific matching
  happens once per chromosome when the reader builds its remap.
- parquet input requires `pyarrow` at runtime.
- Per-query partitioned LDSC wrappers are intentionally out of scope for this
  module; this module computes LD scores only.

Dependencies
------------
- `numpy`
- `pandas`
- `pyarrow` for sorted parquet conversion and runtime queries
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
import gzip
import logging
import os
import sys
from dataclasses import dataclass
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
    MAF_COLUMN_ALIASES,
    MAF_COLUMN_SPEC,
    POS_COLUMN_ALIASES,
    POS_COLUMN_SPEC,
    REFERENCE_METADATA_SPEC_MAP,
    RESTRICTION_CHRPOS_SPEC_MAP,
    RESTRICTION_RSID_SPEC_MAP,
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
    ANNOTATION_SUFFIXES,
    FREQUENCY_SUFFIXES,
    PARQUET_SUFFIXES,
    resolve_chromosome_group,
    resolve_file_group,
    resolve_plink_prefix,
)
from .._row_alignment import assert_same_snp_rows
from . import formats as legacy_parse
from .identifiers import build_snp_id_series, read_snp_restriction_keys
from .overlap import OverlapContribution, compute_overlap
from .plink_bed import __GenotypeArrayInMemory__, PlinkBEDFile  # noqa: F401
from .snp_identity import (
    RestrictionIdentityKeys,
    effective_merge_key_series,
    identity_base_mode,
    identity_mode_family,
    is_allele_aware_mode,
    restriction_membership_mask,
    sidecar_identity_sha256,
)

try:  # pragma: no cover - optional dependency
    import bitarray as ba
except ImportError:  # pragma: no cover - optional dependency
    ba = None


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
# Legacy .l2.ldscore parsing stays permissive: a folded MAF column may appear
# under FRQ/FREQ/FREQUENCY in older files. The registry's MAF_COLUMN_ALIASES is
# narrowed to ("MAF",) for the oriented-vs-folded sidecar distinction, so keep the
# legacy breadth here explicitly.
MAF_ALIASES = ("MAF", "FRQ", "FREQ", "FREQUENCY")


@dataclass
class AnnotationBundle:
    """Per-chromosome annotation payload consumed by the LD-score backends."""
    metadata: pd.DataFrame
    annotations: pd.DataFrame
    baseline_columns: list[str]
    query_columns: list[str]


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
    overlap: OverlapContribution | None = None


# Basic configuration and shared helpers.
def configure_logging(level: str) -> None:
    """Set the ``LDSC`` logger threshold for kernel-side execution.

    Console and file routing is owned by the public workflow boundary
    (``ldsc._logging`` / ``ldsc.cli``); the kernel only adjusts the level so
    direct kernel invocation does not install a root console handler.
    """
    logging.getLogger("LDSC").setLevel(getattr(logging, level.upper()))


def get_legacy_ld_module():
    """Return the bitarray-backed LD kernel or raise a dependency error."""
    if ba is None:
        raise LDSCDependencyError(
            "ldscore could not use PLINK reference-panel mode because the bitarray-backed LD kernel is unavailable. "
            "Most likely the optional dependency 'bitarray' is not installed. Install bitarray and retry."
        )
    return sys.modules[__name__]


def split_arg_list(value: str | None) -> list[str]:
    """Split a comma-delimited CLI argument into normalized path tokens."""
    if not value:
        return []
    out = []
    for token in value.split(","):
        token = os.path.expanduser(os.path.expandvars(token.strip()))
        if token:
            out.append(token)
    return out


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


def getBlockLefts(coords, max_dist):
    """Backward-compatible alias for :func:`get_block_lefts`."""
    return get_block_lefts(coords, max_dist)


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
    """Sort a metadata-like frame by chromosome, position, and SNP name."""
    pos_col = find_column(df.columns, POS_ALIASES)
    if pos_col is None:
        raise KeyError("No POS-like column available for genomic sorting.")
    sort_df = df.copy()
    sort_df["_chrom_key"] = sort_df["CHR"].map(chrom_sort_key)
    sort_df = sort_df.sort_values(by=["_chrom_key", pos_col, "SNP"], kind="mergesort")
    return sort_df.drop(columns="_chrom_key").reset_index(drop=True)


def resolve_annotation_files(spec: str | None) -> list[str]:
    """Resolve comma-delimited annotation tokens into concrete file paths."""
    return resolve_file_group(
        split_arg_list(spec),
        label="annotation input",
        suffixes=ANNOTATION_SUFFIXES,
        allow_chromosome_suite=True,
    )


def resolve_parquet_files(args: argparse.Namespace, chrom: str | None = None) -> list[str]:
    """Resolve the sorted parquet R2 files participating in one LD-score run."""
    r2_table = getattr(args, "r2_table", None)
    if r2_table:
        tokens = split_arg_list(r2_table)
        if chrom is not None:
            return resolve_chromosome_group(
                tokens,
                chrom=chrom,
                label="r2_table",
                suffixes=PARQUET_SUFFIXES,
            )
        return resolve_file_group(
            tokens,
            label="r2_table",
            suffixes=PARQUET_SUFFIXES,
            allow_chromosome_suite=True,
        )
    return []


def resolve_bfile_prefix(args: argparse.Namespace, chrom: str | None = None) -> str | None:
    """Resolve the PLINK prefix for the requested chromosome, if any."""
    if args.bfile is None:
        return None
    return resolve_plink_prefix(args.bfile, chrom=chrom)


def resolve_frequency_files(args: argparse.Namespace, chrom: str | None = None) -> list[str]:
    """Resolve optional frequency or metadata files for one chromosome."""
    tokens = split_arg_list(args.frqfile)
    if not tokens:
        return []
    if chrom is not None:
        return resolve_chromosome_group(
            tokens,
            chrom=chrom,
            label="frqfile",
            suffixes=FREQUENCY_SUFFIXES,
        )
    return resolve_file_group(
        tokens,
        label="frqfile",
        suffixes=FREQUENCY_SUFFIXES,
        allow_chromosome_suite=True,
    )


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
            "by `ldsc build-ref-panel`. Regenerate the reference panel or pass the canonical "
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
            "parquet and sidecar together or regenerate with `ldsc build-ref-panel`. "
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
    if {"IDX_1", "IDX_2", "R2"}.issubset(set(schema_names)):
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
def parse_annotation_file(path: str, chrom: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse one SNP-level annotation table into normalized metadata and values.

    When ``chrom`` is provided, the returned tables are restricted to rows whose
    normalized ``CHR`` matches that chromosome. This is the final safeguard for
    group-style path tokens that resolved to shared multi-chromosome files or
    to globs whose filenames did not encode chromosome labels clearly enough for
    earlier path-level filtering. For LD-score calculation, the annotation
    file's ``CM`` is the first source. The sidecar metadata only fills missing
    ``CM`` values.
    """
    df = read_text_table(path)
    context = path
    chr_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CHR"], context=context)
    pos_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["POS"], context=context)
    snp_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["SNP"], context=context)
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
            "Provide both allele columns, or remove both and run with a base SNP identifier mode."
        )

    meta = pd.DataFrame(
        {
            "CHR": df[chr_col],
            "POS": df[pos_col],
            "SNP": df[snp_col],
        }
    )
    meta["CHR"] = meta["CHR"].map(lambda value: normalize_chromosome(value, context=path))
    meta["POS"] = pd.to_numeric(meta["POS"], errors="raise").astype(np.int64)
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

    meta = sort_frame_by_genomic_position(meta)
    metadata_source_columns = {chr_col, pos_col, snp_col, cm_col, maf_col, a1_col, a2_col}
    annotation_columns = [col for col in df.columns if col not in metadata_source_columns]
    if not annotation_columns:
        raise LDSCInputError(
            f"ldscore could not parse annotation file '{path}': no annotation value "
            "columns remain after metadata columns were removed. Most likely the file "
            "contains only CHR/POS/SNP/CM metadata, or the annotation columns were named "
            "as metadata aliases. Add at least one annotation column with numeric values."
        )

    sorted_df = df.copy()
    sorted_df["_CHR"] = sorted_df[chr_col].map(lambda value: normalize_chromosome(value, context=path))
    sorted_df["_POS"] = pd.to_numeric(sorted_df[pos_col], errors="raise").astype(np.int64)
    sorted_df["_SNP"] = sorted_df[snp_col].astype(str)
    sorted_df["_chrom_key"] = sorted_df["_CHR"].map(chrom_sort_key)
    sorted_df = sorted_df.sort_values(by=["_chrom_key", "_POS", "_SNP"], kind="mergesort").drop(columns="_chrom_key")
    annotations = sorted_df.loc[:, annotation_columns].astype(np.float32).reset_index(drop=True)
    meta = meta.reset_index(drop=True)
    if len(meta) != len(annotations):
        raise LDSCInternalError(
            f"LD-score annotation parsing failed for '{path}': metadata and annotation "
            "lengths diverged after sorting. Most likely an internal sorting step was "
            "applied inconsistently. Re-run with DEBUG logging and report the traceback."
        )

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
            meta, annotations = parse_annotation_file(path, chrom=chrom)
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
    if not getattr(args, "regression_snps_file", None):
        return None
    return read_identifier_list(args.regression_snps_file, args.snp_identifier)


def parse_frequency_metadata(path: str, chrom: str | None, identifier_mode: str) -> pd.DataFrame:
    """Parse an optional sidecar metadata file into ``_key``/``CM``/``MAF`` columns."""
    df = read_text_table(path)
    identifier_mode = normalize_snp_identifier_mode(identifier_mode)
    context = path
    chr_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=context)
    pos_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=context)
    snp_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=context)
    cm_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=context)
    maf_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["MAF"], context=context)
    a1_col = resolve_optional_column(df.columns, A1_COLUMN_SPEC, context=context)
    a2_col = resolve_optional_column(df.columns, A2_COLUMN_SPEC, context=context)
    if (a1_col is None) != (a2_col is None):
        raise LDSCInputError(
            f"ldscore could not parse frequency metadata '{path}': it has only one allele "
            "column. Most likely the file contains A1 without A2, or A2 without A1. "
            "Provide both allele columns, or remove both and use a base SNP identifier mode."
        )

    if chrom is not None and chr_col is not None and identity_mode_family(identifier_mode) == "rsid":
        keep = df[chr_col].map(lambda value: normalize_chromosome(value, context=path)) == normalize_chromosome(chrom, context=path)
        df = df.loc[keep].reset_index(drop=True)

    if len(df) == 0:
        return pd.DataFrame(columns=["_key", "_key_mode", "CM", "MAF"])

    out = pd.DataFrame(index=df.index)
    key_frame = pd.DataFrame(index=df.index)
    if identity_mode_family(identifier_mode) == "rsid":
        if snp_col is None:
            raise LDSCInputError(
                f"ldscore could not parse frequency metadata '{path}' in rsID-family mode: "
                "no SNP column was found. Most likely the sidecar uses CHR/POS-only "
                "identifiers or an unrecognized SNP column name. Add a SNP column or run "
                "with a chr_pos-family SNP identifier mode."
            )
        key_frame["SNP"] = df[snp_col].astype(str)
    else:
        if chr_col is None or pos_col is None:
            raise LDSCInputError(
                f"ldscore could not parse frequency metadata '{path}' in chr_pos-family mode: "
                "CHR and POS columns are required. Most likely the sidecar uses rsID-only "
                "identifiers or unrecognized CHR/POS column names. Add CHR and POS columns "
                "or run with an rsID-family SNP identifier mode."
            )
        key_frame["CHR"] = df[chr_col]
        key_frame["POS"] = df[pos_col]
    has_alleles = a1_col is not None and a2_col is not None
    if has_alleles:
        key_frame["A1"] = df[a1_col].astype(str)
        key_frame["A2"] = df[a2_col].astype(str)
    key_mode = identifier_mode if is_allele_aware_mode(identifier_mode) and has_alleles else identity_base_mode(identifier_mode)
    out["_key"] = build_snp_id_series(key_frame, key_mode)
    out["_key_mode"] = key_mode

    if cm_col is not None:
        out["CM"] = pd.to_numeric(df[cm_col], errors="coerce")
    if maf_col is not None:
        maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
        out["MAF"] = np.minimum(maf, 1.0 - maf)

    out = out.loc[out["_key"].notna()].copy()
    if chrom is not None and chr_col is not None:
        keep_chrom = df.loc[out.index, chr_col].map(lambda value: normalize_chromosome(value, context=path)) == normalize_chromosome(chrom, context=path)
        out = out.loc[keep_chrom.to_numpy()].copy()
    return out.reset_index(drop=True)


def merge_frequency_metadata(
    metadata: pd.DataFrame,
    args: argparse.Namespace,
    chrom: str,
    identifier_mode: str,
) -> pd.DataFrame:
    """
    Merge reference-panel sidecar ``CM`` and ``MAF`` onto retained SNP metadata.

    The reference-panel sidecar is **authoritative** for ``CM`` and ``MAF``: for
    every SNP the sidecar covers, its ``CM``/``MAF`` overwrite any
    annotation-provided values. ``CM``/``MAF`` are population-specific and belong
    to the reference panel, not the annotation. Frequency sidecars are metadata
    providers, not SNP filters: when multiple sidecar rows share the active
    effective identity key, the whole duplicate-key cluster is ignored and a
    warning is logged, so those SNPs keep whatever ``CM``/``MAF`` they already had
    (typically missing, which `--ld-wind-cm` then rejects).
    """
    files = resolve_frequency_files(args, chrom=chrom)
    if not files:
        return metadata

    frames = [parse_frequency_metadata(path, chrom=chrom, identifier_mode=identifier_mode) for path in files]
    freq_df = pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame(columns=["_key", "_key_mode", "CM", "MAF"])
    duplicate_mask = freq_df[["_key_mode", "_key"]].duplicated(keep=False)
    if bool(duplicate_mask.any()):
        message = (
            f"Dropping {int(duplicate_mask.sum())} frequency metadata rows in duplicate SNP identity clusters "
            f"for chromosome {chrom}; CM/MAF will remain missing for those keys unless already present in "
            "annotation metadata."
        )
        LOGGER.warning(message)
        freq_df = freq_df.loc[~duplicate_mask].reset_index(drop=True)

    merged = metadata.copy()
    for key_mode, freq_part in freq_df.groupby("_key_mode", sort=False):
        merged_keys = build_snp_id_series(merged, str(key_mode))
        freq_part = freq_part.set_index("_key")
        # The reference-panel sidecar is authoritative for CM/MAF: overwrite any
        # annotation-provided values for every SNP the sidecar covers, so the
        # reference panel (not the annotation) determines CM and MAF.
        if "CM" in freq_part.columns:
            if "CM" not in merged.columns:
                merged["CM"] = np.nan
            present_cm = merged_keys.isin(freq_part.index)
            merged.loc[present_cm, "CM"] = freq_part.loc[merged_keys.loc[present_cm], "CM"].to_numpy()
        if "MAF" in freq_part.columns:
            if "MAF" not in merged.columns:
                merged["MAF"] = np.nan
            present_maf = merged_keys.isin(freq_part.index)
            merged.loc[present_maf, "MAF"] = freq_part.loc[merged_keys.loc[present_maf], "MAF"].to_numpy()
    return merged


def apply_maf_filter(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    maf_min: float | None,
    context: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply the optional MAF threshold to metadata and annotation rows together."""
    if maf_min is None:
        return metadata, annotations
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        LOGGER.warning(f"Cannot apply --maf-min in {context} because MAF metadata is unavailable.")
        return metadata, annotations
    keep = metadata["MAF"] >= maf_min
    removed = int((~keep).sum())
    if removed:
        LOGGER.info(f"Removed {removed} SNPs with MAF < {maf_min} in {context}.")
    metadata = metadata.loc[keep].reset_index(drop=True)
    annotations = annotations.loc[keep].reset_index(drop=True)
    return metadata, annotations


def chromosome_set_from_annotation_inputs(args: argparse.Namespace) -> list[str]:
    """Discover the chromosome set implied by the supplied annotation inputs."""
    chromosomes: set[str] = set()
    all_files = resolve_annotation_files(args.query_annot) + resolve_annotation_files(args.baseline_annot)
    for path in all_files:
        df = read_text_table(path)
        if "CHR" not in df.columns:
            raise LDSCInputError(
                f"ldscore could not discover chromosomes from annotation file '{path}': "
                "the file is missing a CHR column. Most likely this is not an LDSC "
                "annotation table or the header uses an unrecognized chromosome name. "
                "Provide annotation files with CHR/POS/SNP/CM metadata columns."
            )
        chromosomes.update(df["CHR"].map(lambda value: normalize_chromosome(value, context=path)).unique().tolist())

    if not chromosomes:
        raise LDSCInputError(
            "ldscore could not resolve any chromosomes from the supplied annotation inputs. "
            "Most likely the annotation paths are empty after chromosome filtering or contain "
            "no CHR values. Check the annotation headers and pass files that contain retained "
            "chromosome rows."
        )
    return sorted(chromosomes, key=chrom_sort_key)


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
            "build, or use `--ld-wind-kb` / `--ld-wind-snps`."
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
            "`ldsc build-ref-panel`."
        )
    mode = mode_raw.decode("utf-8")
    if mode not in {"snps", "kb", "cm"}:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            f"unsupported parquet ldsc:ld_window_mode={mode!r}. Most likely the "
            "artifact was edited or written by an incompatible development "
            "version. Regenerate the reference panel with `ldsc build-ref-panel`."
        )
    try:
        value = float(value_raw.decode("utf-8"))
    except ValueError as exc:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            f"parquet ldsc:ld_window_value={value_raw.decode('utf-8', errors='replace')!r} "
            "is not numeric. Most likely the artifact schema metadata was edited "
            "or corrupted. Regenerate the reference panel with `ldsc build-ref-panel`."
        ) from exc
    if value <= 0:
        raise LDSCInputError(
            f"ldscore could not validate the LD window for R2 parquet '{path}': "
            f"parquet ldsc:ld_window_value={value} is not positive. Most likely "
            "the artifact schema metadata was edited or corrupted. Regenerate the "
            "reference panel with `ldsc build-ref-panel`."
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
    Query block-local dense R2 matrices from a per-chromosome parquet table.

    Index-format parquet files use logical fields ``IDX_1``, ``IDX_2``, ``R2``
    (int16 on-disk, dequantized to float32 by dividing by ``ldsc:r2_scale``),
    and ``SIGN``. Files are queried via row-group pruning on sorted ``IDX_1``
    bounds. Legacy float32 ``R2`` columns (absent ``ldsc:r2_encoding`` metadata)
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

        layout = _parquet_schema_layout(pq.ParquetFile(paths[0]).schema_arrow.names)
        if layout != "index":
            raise LDSCInputError(
                f"ldscore could not use R2 parquet '{paths[0]}': it is not an index-format "
                "R2 parquet with columns IDX_1/IDX_2/R2/SIGN. Most likely this file was "
                "written by an old LDSC version or is not an LDSC R2 artifact. Regenerate "
                "the reference panel with `ldsc build-ref-panel`. "
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
        self._runtime_layout = "index"
        self._pf = pq.ParquetFile(paths[0])
        self._init_index_path(paths[0], metadata)

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
                "`ldsc build-ref-panel`. "
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
                "with `ldsc build-ref-panel`. "
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
                    "`ldsc:n_samples`. Rebuild the panel with the current build-ref-panel, "
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

    def _decode_index_row_group(self, row_group_index: int) -> _DecodedR2RowGroup:
        """Decode one index row group: dequantize R2, gather endpoints through the remap.

        ``IDX_1``/``IDX_2`` are panel (build) indices; ``self._remap`` maps each
        to its retained matrix index (or ``-1`` when the endpoint SNP is not in
        the analysis universe). Pairs with either endpoint dropped are removed.
        When ``self._r2_scale`` is set (int16 panels), the raw int16 column is
        divided by the scale to produce float32 before the raw→unbiased transform.
        ``SIGN`` is not read: it is unused by LD-score computation.
        """
        if self._pf is None:
            raise LDSCInternalError(
                "LD-score parquet reader failed in SortedR2BlockReader._decode_index_row_group(): "
                "the parquet file handle is not initialized. Most likely row-group decoding "
                "was called before reader setup completed. Re-run with DEBUG logging and "
                "report the traceback."
            )
        table = self._pf.read_row_group(int(row_group_index), columns=["IDX_1", "IDX_2", "R2"])
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

    def iter_all_pairs(self):
        """Yield ``(i, j, r2)`` arrays for every row group, decoding each once.

        Streams the whole chromosome's stored pairs in IDX_1 order with no
        window pruning and no row-group cache: each row group is decoded a
        single time, remapped to retained indices, and dropped endpoints
        removed. Empty groups (all endpoints outside the analysis universe) are
        skipped.
        """
        if self._pf is None:
            raise LDSCInternalError(
                "LD-score parquet reader failed in SortedR2BlockReader.iter_all_pairs(): "
                "the parquet file handle is not initialized. Re-run with DEBUG logging "
                "and report the traceback."
            )
        for rg_index in range(self._pf.metadata.num_row_groups):
            group = self._decode_index_row_group(rg_index)
            if group.i.size == 0:
                continue
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
    annot: np.ndarray,
    block_reader: SortedR2BlockReader,
    chunk_pairs: int = _CSR_CHUNK_PAIRS,
) -> np.ndarray:
    """Compute ``cor_sum = R @ annot`` by streaming the parquet's stored R2 pairs.

    The diagonal (R2=1) seeds ``cor_sum`` from ``annot``; stored within-window pairs
    are buffered into chunks of ~``chunk_pairs`` and accumulated with a float64 CSR
    SpMM (:func:`_accumulate_pair_contributions`), restricted to the ldscore window
    by ``block_left``. Each stored pair is touched exactly once (O(nnz * n_a)).
    """
    block_left = np.asarray(block_left, dtype=np.int64)
    annot64 = annot.astype(np.float64, copy=False)
    cor_sum = annot64.copy()  # diagonal R2 = 1 for every SNP
    buf_i: list[np.ndarray] = []
    buf_j: list[np.ndarray] = []
    buf_r2: list[np.ndarray] = []
    buffered = 0

    def flush() -> None:
        nonlocal buffered
        if buffered == 0:
            return
        _accumulate_pair_contributions(
            cor_sum,
            np.concatenate(buf_i),
            np.concatenate(buf_j),
            np.concatenate(buf_r2),
            annot64,
            block_left,
        )
        buf_i.clear()
        buf_j.clear()
        buf_r2.clear()
        buffered = 0

    for i, j, r2 in block_reader.iter_all_pairs():
        buf_i.append(i)
        buf_j.append(j)
        buf_r2.append(r2)
        buffered += int(i.size)
        if buffered >= chunk_pairs:
            flush()
    flush()
    return np.asarray(cor_sum, dtype=np.float32)


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
) -> np.ndarray:
    """Build the binary mask column used to compute regression-weight LD scores."""
    if regression_keys is None:
        return np.ones(len(metadata), dtype=np.float32)
    if isinstance(regression_keys, RestrictionIdentityKeys):
        return restriction_membership_mask(
            metadata,
            regression_keys,
            identifier_mode,
            context="LD-score regression SNP restriction matching",
        ).to_numpy(dtype=np.float32)
    keys = identifier_keys(metadata, identifier_mode)
    return keys.isin(regression_keys).to_numpy(dtype=np.float32)


# Per-chromosome compute backends.
def compute_chrom_from_parquet(
    chrom: str,
    bundle: AnnotationBundle,
    args: argparse.Namespace,
    regression_keys: set[str] | RestrictionIdentityKeys | None,
) -> ChromComputationResult:
    """
    Compute all LD-score outputs for one chromosome from sorted parquet R2 input.

    Main steps:
    1. Merge optional CM/MAF metadata and apply any MAF filter. For LD-score
       calculation, the annotation file's CM is the first source. The sidecar
       metadata only fills missing CM values.
    2. Treat the sidecar-aligned metadata rows as the authoritative retained SNP universe.
    3. Validate the LD window against the retained SNP coordinates.
    4. Query block-local dense R2 matrices from the sorted parquet file while
       following the old LDSC sliding-block traversal.
    5. Return chromosome-level LD scores plus all-SNP and common-SNP counts.
    """
    metadata = merge_frequency_metadata(bundle.metadata.copy(), args, chrom=chrom, identifier_mode=args.snp_identifier)
    metadata, annotations = apply_maf_filter(
        metadata,
        bundle.annotations.copy(),
        getattr(args, "maf_min", getattr(args, "maf", None)),
        context="parquet mode",
    )
    if len(metadata) == 0:
        raise LDSCInputError(
            f"ldscore retained no annotation SNPs on chromosome {chrom} after parquet "
            "sidecar alignment. Most likely the annotation SNP identifiers, genome build, "
            "or allele-aware identifier mode do not match the parquet R2 sidecar. Use "
            "annotation and R2 artifacts built with the same SNP identifier mode and "
            "genome build. "
            f"Other causes & fixes: {_LDSCORE_INTERSECTION_DOC}"
        )

    validate_retained_identifier_uniqueness(metadata, args.snp_identifier, chrom)
    require_reference_maf(metadata, chrom)
    coords, max_dist = build_window_coordinates(metadata, args)
    block_left = get_block_lefts(
        coords,
        max_dist,
    )
    regression_mask = regression_mask_from_keys(metadata, regression_keys, args.snp_identifier).reshape(-1, 1)
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=True)
    combined_annot = np.c_[annot_matrix, regression_mask]
    parquet_paths = resolve_parquet_files(args, chrom=chrom)
    validate_ldscore_window_within_r2_panel_window(
        args,
        parquet_paths=parquet_paths,
        chrom=chrom,
    )
    check_whole_chromosome_window(block_left, args, chrom)
    block_reader = SortedR2BlockReader(
        paths=parquet_paths,
        chrom=chrom,
        metadata=metadata,
        identifier_mode=args.snp_identifier,
        r2_bias_mode=getattr(args, "r2_bias_mode", "unbiased"),
        r2_sample_size=getattr(args, "r2_sample_size", None),
        genome_build=args.genome_build,
    )
    combined_scores = ld_score_streaming_from_r2_reader(
        block_left=block_left,
        annot=combined_annot,
        block_reader=block_reader,
    )
    ld_scores = combined_scores[:, :annot_matrix.shape[1]]
    w_ld = combined_scores[:, annot_matrix.shape[1]:]
    out_metadata = metadata.reset_index(drop=True)
    if "MAF" not in out_metadata.columns:
        out_metadata["MAF"] = np.nan
    M, M_5_50 = compute_counts(
        out_metadata,
        annotations,
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
    overlap = compute_overlap(
        out_metadata,
        annotations,
        n_baseline=len(bundle.baseline_columns),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
    return ChromComputationResult(
        chrom=chrom,
        metadata=out_metadata,
        ld_scores=ld_scores,
        w_ld=w_ld,
        M=M,
        M_5_50=M_5_50,
        ldscore_columns=bundle.baseline_columns + bundle.query_columns,
        baseline_columns=bundle.baseline_columns,
        query_columns=bundle.query_columns,
        overlap=overlap,
    )


def compute_chrom_from_plink(
    chrom: str,
    bundle: AnnotationBundle,
    args: argparse.Namespace,
    regression_keys: set[str] | RestrictionIdentityKeys | None,
) -> ChromComputationResult:
    """
    Compute all LD-score outputs for one chromosome from a PLINK reference panel.

    Main steps:
    1. Align annotation SNPs to the PLINK BIM table.
    2. Reuse the legacy PLINK genotype reader and LD-score kernel.
    3. Compute partitioned reference LD scores and one-column regression-universe LD scores.
    4. Return chromosome-level LD scores plus all-SNP and common-SNP counts.
    """
    legacy_ld = get_legacy_ld_module()
    prefix = resolve_bfile_prefix(args, chrom=chrom)
    if prefix is None:
        raise LDSCUsageError(
            "ldscore cannot run PLINK mode without a PLINK prefix. Most likely PLINK "
            "mode was selected but `--plink-prefix`/`--bfile` was omitted. Pass the "
            "prefix shared by the `.bed`, `.bim`, and `.fam` files, or use parquet "
            "mode with `--r2-dir`."
        )

    bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
    fam = legacy_parse.PlinkFAMFile(prefix + ".fam")
    panel_df = bim.df.loc[:, ["CHR", "SNP", "CM", "BP", "A1", "A2"]].copy().rename(columns={"BP": "POS"})
    panel_df["CHR"] = panel_df["CHR"].map(lambda value: normalize_chromosome(value, context=prefix + ".bim"))
    panel_df["SNP"] = panel_df["SNP"].astype(str)
    panel_df["POS"] = pd.to_numeric(panel_df["POS"], errors="raise").astype(np.int64)
    panel_df["CM"] = pd.to_numeric(panel_df["CM"], errors="coerce")
    panel_df["A1"] = panel_df["A1"].astype(str)
    panel_df["A2"] = panel_df["A2"].astype(str)
    panel_df = panel_df.loc[panel_df["CHR"] == normalize_chromosome(chrom, context=prefix + ".bim")].copy()
    if len(panel_df) == 0:
        raise LDSCInputError(
            f"ldscore found no PLINK SNPs for chromosome {chrom} in prefix '{prefix}'. "
            "Most likely the `.bim` file does not contain that chromosome or the chromosome "
            "labels do not match the annotation inputs. Pass the correct chromosome-specific "
            "PLINK prefix or rebuild the reference panel with matching chromosome labels."
        )
    panel_df["_key"] = identifier_keys(panel_df, args.snp_identifier)

    metadata = bundle.metadata.copy()
    annotations = bundle.annotations.copy()
    metadata["_key"] = identifier_keys(metadata, args.snp_identifier)

    key_to_panel_index = {key: idx for key, idx in zip(panel_df["_key"], panel_df.index)}
    keep = metadata["_key"].isin(key_to_panel_index)
    removed = int((~keep).sum())
    if removed:
        LOGGER.warning(
            f"Dropping {removed} annotated SNPs on chromosome {chrom} "
            "because they are absent from the PLINK reference panel."
        )
    metadata = metadata.loc[keep].reset_index(drop=True)
    annotations = annotations.loc[keep].reset_index(drop=True)
    if regression_keys is not None and not isinstance(regression_keys, RestrictionIdentityKeys):
        regression_keys = regression_keys.intersection(set(metadata["_key"]))
    if len(metadata) == 0:
        raise LDSCInputError(
            f"ldscore retained no annotation SNPs on chromosome {chrom} after PLINK "
            "intersection. Most likely the annotation SNP identifiers, genome build, "
            "or allele-aware identifier mode do not match the PLINK `.bim` file. Use "
            "annotation and PLINK reference files built with the same SNP identifier mode "
            "and genome build. "
            f"Other causes & fixes: {_LDSCORE_INTERSECTION_DOC}"
        )

    keep_indivs = resolve_keep_individuals(getattr(args, "keep", None), fam)
    keep_snps = [key_to_panel_index[key] for key in metadata["_key"]]
    geno = legacy_ld.PlinkBEDFile(
        prefix + ".bed",
        len(fam.IDList),
        bim,
        keep_snps=keep_snps,
        keep_indivs=keep_indivs,
        mafMin=getattr(args, "maf_min", getattr(args, "maf", None)),
    )

    geno_meta = pd.DataFrame(geno.df, columns=geno.colnames)
    if "BP" in geno_meta.columns:
        geno_meta = geno_meta.rename(columns={"BP": "POS"})
    geno_meta["CHR"] = geno_meta["CHR"].map(lambda value: normalize_chromosome(value, context=prefix + ".bim"))
    geno_meta["SNP"] = geno_meta["SNP"].astype(str)
    geno_meta["POS"] = pd.to_numeric(geno_meta["POS"], errors="raise").astype(np.int64)
    geno_meta["CM"] = pd.to_numeric(geno_meta["CM"], errors="coerce")
    geno_meta["MAF"] = pd.to_numeric(geno_meta["MAF"], errors="coerce")
    geno_meta = geno_meta.merge(
        panel_df.loc[:, ["CHR", "SNP", "POS", "A1", "A2"]],
        how="left",
        on=["CHR", "SNP", "POS"],
        sort=False,
    )
    geno_meta["_key"] = identifier_keys(geno_meta, args.snp_identifier)

    annotation_matrix = annotations.set_index(metadata["_key"]).loc[geno_meta["_key"]]

    if len(geno_meta) == 0:
        raise LDSCInputError(
            f"ldscore retained no PLINK reference SNPs on chromosome {chrom} after genotype "
            "filtering. Most likely every candidate SNP is monomorphic, --maf-min is too high, "
            "or no annotation SNPs overlap the panel. Lower --maf-min, or check the SNP "
            "overlap and genome build."
        )
    require_reference_maf(geno_meta, chrom)
    if args.ld_wind_cm is not None:
        assert_cm_usable(geno_meta["CM"], chrom)
    coords, max_dist = build_window_coordinates(geno_meta.drop(columns="_key"), args)
    block_left = legacy_ld.getBlockLefts(coords, max_dist)
    check_whole_chromosome_window(block_left, args, chrom)

    ld_scores = geno.ldScoreVarBlocks(block_left, args.snp_batch_size, annot=annotation_matrix.to_numpy(dtype=np.float32))
    regression_mask = regression_mask_from_keys(geno_meta.drop(columns="_key"), regression_keys, args.snp_identifier)
    geno._currentSNP = 0
    w_ld = geno.ldScoreVarBlocks(block_left, args.snp_batch_size, annot=regression_mask.reshape(-1, 1))
    out_metadata = geno_meta.drop(columns="_key").reset_index(drop=True)
    M, M_5_50 = compute_counts(
        out_metadata,
        annotation_matrix.reset_index(drop=True),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
    overlap = compute_overlap(
        out_metadata,
        annotation_matrix.reset_index(drop=True),
        n_baseline=len(bundle.baseline_columns),
        common_maf_min=getattr(args, "common_maf_min", 0.05),
    )
    return ChromComputationResult(
        chrom=chrom,
        metadata=out_metadata,
        ld_scores=np.asarray(ld_scores, dtype=np.float32),
        w_ld=np.asarray(w_ld, dtype=np.float32),
        M=M,
        M_5_50=M_5_50,
        ldscore_columns=bundle.baseline_columns + bundle.query_columns,
        baseline_columns=bundle.baseline_columns,
        query_columns=bundle.query_columns,
        overlap=overlap,
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


def result_to_dataframe(result: ChromComputationResult) -> pd.DataFrame:
    """Materialize one chromosome result as a standard ``.l2.ldscore`` table."""
    df = result.metadata.copy()
    for idx, column in enumerate(result.ldscore_columns):
        df[column + "L2"] = result.ld_scores[:, idx]
    return df


def weight_result_to_dataframe(result: ChromComputationResult) -> pd.DataFrame:
    """Materialize one chromosome regression-weight result as a one-column table."""
    df = result.metadata.copy()
    df["L2"] = np.ravel(result.w_ld)
    return df


# Output assembly.
def write_ldscore_file(df: pd.DataFrame, path: str) -> None:
    """Write one LDSC-compatible LD-score table, preserving metadata columns first."""
    out = df.copy()
    out = out.loc[:, [col for col in ANNOT_META_COLUMNS if col in out.columns] + [col for col in out.columns if col not in ANNOT_META_COLUMNS]]
    with gzip.open(path, "wt") as handle:
        out.to_csv(handle, sep="\t", index=False, na_rep="NA", float_format="%.6g")


def write_counts(path: str, counts: np.ndarray) -> None:
    """Write one LDSC ``.M``-style count vector to disk."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(str(x) for x in counts))


def write_annotation_groups(path: str, baseline_columns: Sequence[str], query_columns: Sequence[str]) -> None:
    """Write the baseline/query annotation-group manifest used by downstream tools."""
    rows = [{"annotation": col, "group": "baseline"} for col in baseline_columns]
    rows.extend({"annotation": col, "group": "query"} for col in query_columns)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def aggregate_results(results: Sequence[ChromComputationResult]) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray | None]:
    """
    Combine chromosome-level results into one genome-wide LD-score output set.

    Returns the aggregated reference LD-score table, regression-weight table,
    and the chromosome-summed all-SNP / common-SNP count vectors.
    """
    ld_frames = [result_to_dataframe(result) for result in results]
    weight_frames = [weight_result_to_dataframe(result) for result in results]
    ld_df = pd.concat(ld_frames, axis=0, ignore_index=True)
    weight_df = pd.concat(weight_frames, axis=0, ignore_index=True)
    ld_df = sort_frame_by_genomic_position(ld_df)
    weight_df = sort_frame_by_genomic_position(weight_df)
    M = np.sum(np.vstack([result.M for result in results]), axis=0)
    all_have_m_5_50 = all(result.M_5_50 is not None for result in results)
    M_5_50 = np.sum(np.vstack([result.M_5_50 for result in results]), axis=0) if all_have_m_5_50 else None
    return ld_df, weight_df, M, M_5_50


def emit_outputs(results: Sequence[ChromComputationResult], args: argparse.Namespace) -> None:
    """Write either per-chromosome outputs or one aggregated LDSC-compatible result set."""
    if not results:
        raise LDSCInternalError(
            "LD-score output emission failed in emit_outputs(): no chromosome results were "
            "provided. Most likely all chromosomes were skipped before output writing. "
            "Re-run with DEBUG logging and report the traceback."
        )

    write_annotation_groups(
        args.out + ".annotation_groups.tsv",
        results[0].baseline_columns,
        results[0].query_columns,
    )

    if args.per_chr_output:
        for result in results:
            prefix = f"{args.out}.{result.chrom}"
            write_ldscore_file(result_to_dataframe(result), prefix + ".l2.ldscore.gz")
            write_ldscore_file(weight_result_to_dataframe(result), prefix + ".w.l2.ldscore.gz")
            write_counts(prefix + ".l2.M", result.M)
            if result.M_5_50 is not None:
                write_counts(prefix + ".l2.M_5_50", result.M_5_50)
            else:
                LOGGER.warning(f"Skipping {prefix}.l2.M_5_50 because MAF is unavailable.")
        return

    ld_df, weight_df, M, M_5_50 = aggregate_results(results)
    write_ldscore_file(ld_df, args.out + ".l2.ldscore.gz")
    write_ldscore_file(weight_df, args.out + ".w.l2.ldscore.gz")
    write_counts(args.out + ".l2.M", M)
    if M_5_50 is not None:
        write_counts(args.out + ".l2.M_5_50", M_5_50)
    else:
        LOGGER.warning(f"Skipping {args.out}.l2.M_5_50 because MAF is unavailable.")


def _first_resolved_r2_parquet(args: argparse.Namespace) -> str | None:
    """
    Return the first existing R2 parquet path available during validation.

    Validation remains permissive for legacy tests and callers that pass paths
    not yet present on disk. When a direct parquet file or package-built R2
    directory is resolvable, the caller can inspect its Arrow schema metadata
    before applying default R2 bias settings.
    """
    try:
        paths = resolve_parquet_files(args)
    except (FileNotFoundError, LDSCInputError):
        return None
    for path in paths:
        candidate = Path(path)
        if candidate.is_dir():
            matches = sorted(candidate.glob("chr*_r2.parquet"))
            if matches:
                return str(matches[0])
            continue
        if candidate.is_file():
            return str(candidate)
    return None


def validate_args(args: argparse.Namespace) -> None:
    """Validate top-level LD-score CLI arguments before any heavy work starts."""
    for attr in ("query_annot_chr", "baseline_annot_chr", "bfile_chr", "r2_table_chr", "frqfile_chr"):
        if not hasattr(args, attr):
            setattr(args, attr, None)
    args.snp_identifier = normalize_snp_identifier_mode(args.snp_identifier)
    args.genome_build = normalize_genome_build(args.genome_build)
    assert args.genome_build in {"hg19", "hg38", None}, (
        f"genome_build={args.genome_build!r} must be concrete by this point."
    )
    keep = getattr(args, "keep", None)
    if not (args.query_annot or args.baseline_annot):
        raise LDSCUsageError(
            "ldscore cannot run without annotation inputs. Most likely neither "
            "`--baseline-annot` nor `--query-annot` was supplied to the standalone "
            "kernel CLI. Pass at least one annotation input path."
        )
    if bool(args.r2_table) == bool(args.bfile):
        raise LDSCUsageError(
            "ldscore could not choose a reference-panel backend. Most likely both parquet "
            "R2 input and PLINK input were supplied, or neither was supplied. Pass exactly "
            "one of `--r2-table` for parquet mode or `--bfile` for PLINK mode."
        )
    if args.r2_table:
        if keep:
            raise LDSCUsageError(
                "ldscore cannot apply `--keep-indivs-file` in parquet R2 mode. Most likely "
                "`--keep-indivs-file` was combined with parquet R2 input, but individual-level "
                "filtering only exists before PLINK genotype LD calculation. Remove "
                "`--keep-indivs-file`, or rerun in PLINK mode with `--bfile`."
            )
        first_r2 = _first_resolved_r2_parquet(args)
        if first_r2 is not None:
            from . import ref_panel as ref_panel_mod

            stored = ref_panel_mod._read_r2_schema_meta(first_r2)
            # R2 bias mode and sample size come from parquet schema metadata
            # (ldsc:r2_bias / ldsc:n_samples); there is no CLI override.
            args.r2_bias_mode, args.r2_sample_size = ref_panel_mod._resolve_r2_bias_from_meta(
                getattr(args, "r2_bias_mode", None),
                getattr(args, "r2_sample_size", None),
                stored,
            )
        if getattr(args, "r2_bias_mode", None) is None:
            args.r2_bias_mode = "unbiased"
            args.r2_sample_size = getattr(args, "r2_sample_size", None)
        if identity_mode_family(args.snp_identifier) == "chr_pos" and args.genome_build is None:
            raise LDSCUsageError(
                "ldscore cannot run parquet R2 mode with chr_pos-family SNP identifiers "
                "without a genome build. Most likely `--snp-identifier chr_pos` or an "
                "allele-aware chr_pos mode was supplied without `--genome-build`. Pass "
                "`--genome-build hg19`, `--genome-build hg38`, or `--genome-build auto`."
            )
    if args.ld_wind_cm is not None and args.ld_wind_cm <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--ld-wind-cm={args.ld_wind_cm}`. Most likely "
            "the LD window was set to zero or a negative value. Pass a positive "
            "centimorgan window."
        )
    if args.ld_wind_kb is not None and args.ld_wind_kb <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--ld-wind-kb={args.ld_wind_kb}`. Most likely "
            "the LD window was set to zero or a negative value. Pass a positive "
            "kilobase window."
        )
    if args.ld_wind_snps is not None and args.ld_wind_snps <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--ld-wind-snps={args.ld_wind_snps}`. Most likely "
            "the LD window was set to zero or a negative SNP count. Pass a positive "
            "SNP-count window."
        )
    if getattr(args, "maf_min", None) is not None and not 0 <= args.maf_min <= 0.5:
        raise LDSCConfigError(
            f"ldscore received invalid `--maf-min={args.maf_min}`. Most likely the "
            "minor-allele frequency threshold was entered outside the valid [0, 0.5] "
            "range. Pass a value between 0 and 0.5, or omit the option."
        )
    if not 0 <= getattr(args, "common_maf_min", 0.05) <= 0.5:
        raise LDSCConfigError(
            f"ldscore received invalid `--common-maf-min={args.common_maf_min}`. Most "
            "likely the common-SNP MAF threshold was entered outside the valid [0, 0.5] "
            "range. Pass a value between 0 and 0.5."
        )
    if args.snp_batch_size <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid `--snp-batch-size={args.snp_batch_size}`. Most "
            "likely the parquet query batch size was set to zero or a negative value. "
            "Pass a positive integer batch size."
        )


def build_parser() -> argparse.ArgumentParser:
    """Build the standalone LD-score kernel CLI parser."""
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
        allow_abbrev=False,
    )
    parser.add_argument("--out", required=True, help="Output prefix.")
    parser.add_argument("--query-annot", default=None, help="Comma-separated SNP-level query annotation inputs. Each token may be an exact path, glob, or explicit @ chromosome-suite token.")
    parser.add_argument("--baseline-annot", default=None, help="Comma-separated SNP-level baseline annotation inputs. Each token may be an exact path, glob, or explicit @ chromosome-suite token.")
    parser.add_argument("--bfile", default=None, help="PLINK reference-panel prefix or explicit @ chromosome-suite token.")
    parser.add_argument("--r2-table", default=None, help="Comma-separated sorted parquet R2 inputs. Each token may be an exact path, glob, or explicit @ chromosome-suite token.")
    parser.add_argument(
        "--snp-identifier",
        default="chr_pos_allele_aware",
        choices=("rsid", "rsid_allele_aware", "chr_pos", "chr_pos_allele_aware"),
        help="Identifier mode used to match annotations to the reference panel.",
    )
    parser.add_argument(
        "--genome-build",
        default=None,
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build assumed for the sorted parquet R2 file and chr_pos-family matching. Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates.",
    )
    parser.add_argument(
        "--regression-snps-file",
        default=None,
        help=(
            "Optional identity-only SNP list defining the regression SNP set for weight LD computation and written LD-score rows. "
            "Duplicate restriction keys collapse to one retained key; non-identity columns such as CM or MAF are ignored."
        ),
    )
    parser.add_argument(
        "--frqfile",
        default=None,
        help=(
            "Optional frequency/metadata inputs for MAF and CM. Each token may be an exact path, glob, or explicit @ "
            "chromosome-suite token. Duplicate effective SNP identity clusters are dropped entirely before metadata fill."
        ),
    )
    parser.add_argument("--keep", default=None, help="File with individuals to include in LD Score estimation. The file should contain one IID per row.")
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf-min", default=None, type=float, help="Optional MAF filter for retained SNPs when MAF is available.")
    parser.add_argument("--common-maf-min", default=0.05, type=float, help="MAF threshold used only for common-SNP annotation count vectors.")
    parser.add_argument("--snp-batch-size", default=128, type=int, help="Genotype batch size for the PLINK reference-panel backend; ignored by the parquet-R2 backend, which streams stored pairs. Defaults to 128.")
    parser.add_argument("--per-chr-output", default=False, action="store_true", help="Emit per-chromosome outputs instead of an aggregated output.")
    parser.add_argument("--yes-really", default=False, action="store_true", help="Allow whole-chromosome LD windows.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser


# Top-level workflow.
def run_ldscore_from_args(args: argparse.Namespace) -> list[ChromComputationResult]:
    """
    Execute the end-to-end LD-score workflow from parsed CLI-style arguments.

    Main steps:
    1. Validate CLI arguments and discover the chromosomes to process.
    2. Load and align baseline/query annotations for each chromosome.
    3. Dispatch to either the parquet or PLINK chromosome backend.
    4. Aggregate and write LDSC-compatible outputs.
    """
    validate_args(args)
    chromosomes = chromosome_set_from_annotation_inputs(args)
    regression_keys = load_regression_keys(args)
    results: list[ChromComputationResult] = []

    for chrom in chromosomes:
        # Resolve all annotation files contributing columns for this chromosome.
        baseline_files = resolve_chromosome_group(
            split_arg_list(args.baseline_annot),
            chrom=chrom,
            label="baseline_annot",
            suffixes=ANNOTATION_SUFFIXES,
        )
        query_files = resolve_chromosome_group(
            split_arg_list(args.query_annot),
            chrom=chrom,
            label="query_annot",
            suffixes=ANNOTATION_SUFFIXES,
        )

        bundle = combine_annotation_groups(
            baseline_files=baseline_files,
            query_files=query_files,
            chrom=chrom,
            identifier_mode=args.snp_identifier,
        )
        if bundle is None:
            continue

        LOGGER.info(
            f"Processing chromosome {chrom} with {len(bundle.baseline_columns)} baseline "
            f"and {len(bundle.query_columns)} query annotation columns."
        )

        # Compute LD scores from the selected reference-panel backend.
        if args.r2_table:
            result = compute_chrom_from_parquet(chrom, bundle, args, regression_keys)
        else:
            result = compute_chrom_from_plink(chrom, bundle, args, regression_keys)
        results.append(result)

    emit_outputs(results, args)
    return results


def run_ldscore(
    *,
    out: str,
    query_annot: str | None = None,
    baseline_annot: str | None = None,
    bfile: str | None = None,
    r2_table: str | None = None,
    snp_identifier: str = "chr_pos_allele_aware",
    genome_build: str | None = None,
    r2_bias_mode: str | None = None,
    r2_sample_size: float | None = None,
    regression_snps_file: str | None = None,
    frqfile: str | None = None,
    keep: str | None = None,
    ld_wind_snps: int | None = None,
    ld_wind_kb: float | None = None,
    ld_wind_cm: float | None = None,
    maf_min: float | None = None,
    common_maf_min: float = 0.05,
    snp_batch_size: int = 128,
    per_chr_output: bool = False,
    yes_really: bool = False,
    log_level: str = "INFO",
) -> list[ChromComputationResult]:
    """
    Execute the end-to-end LD-score workflow through a normal Python API.

    This is the package-level function entrypoint. It mirrors the CLI options
    with keyword arguments and writes the same outputs as the command-line path.
    """
    args = argparse.Namespace(
        out=out,
        query_annot=query_annot,
        baseline_annot=baseline_annot,
        bfile=bfile,
        r2_table=r2_table,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
        r2_bias_mode=r2_bias_mode,
        r2_sample_size=r2_sample_size,
        regression_snps_file=regression_snps_file,
        frqfile=frqfile,
        query_annot_chr=None,
        baseline_annot_chr=None,
        bfile_chr=None,
        r2_table_chr=None,
        frqfile_chr=None,
        keep=keep,
        ld_wind_snps=ld_wind_snps,
        ld_wind_kb=ld_wind_kb,
        ld_wind_cm=ld_wind_cm,
        maf_min=maf_min,
        common_maf_min=common_maf_min,
        snp_batch_size=snp_batch_size,
        per_chr_output=per_chr_output,
        yes_really=yes_really,
        log_level=log_level,
    )
    configure_logging(log_level)
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> list[ChromComputationResult]:
    """CLI entrypoint for the standalone LD-score kernel module."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return run_ldscore_from_args(args)


if __name__ == "__main__":
    main()
