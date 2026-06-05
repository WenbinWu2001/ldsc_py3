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
REQUIRED_ANNOT_COLUMNS = ("CHR", "POS", "SNP", "CM")
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
MAF_ALIASES = MAF_COLUMN_ALIASES


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
    cm_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CM"], context=context)
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
            "CM": df[cm_col],
        }
    )
    meta["CHR"] = meta["CHR"].map(lambda value: normalize_chromosome(value, context=path))
    meta["POS"] = pd.to_numeric(meta["POS"], errors="raise").astype(np.int64)
    meta["SNP"] = meta["SNP"].astype(str)
    meta["CM"] = pd.to_numeric(meta["CM"], errors="coerce")
    if a1_col is not None and a2_col is not None:
        meta["A1"] = df[a1_col]
        meta["A2"] = df[a2_col]
    if maf_col is not None:
        meta["MAF"] = pd.to_numeric(df[maf_col], errors="coerce")

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
    Merge optional sidecar ``CM`` and ``MAF`` onto retained SNP metadata.

    For LD-score calculation, the annotation file's ``CM`` is the first source.
    The sidecar metadata only fills missing ``CM`` values. ``MAF`` is filled by
    the same missing-value rule so annotation-provided metadata remains
    authoritative when present. Frequency sidecars are metadata providers, not
    SNP filters: when multiple sidecar rows share the active effective identity
    key, the whole duplicate-key cluster is ignored and a warning is logged, so
    ``CM``/``MAF`` stay missing unless already supplied by annotation metadata.
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
        if "CM" in freq_part.columns:
            if "CM" not in merged.columns:
                merged["CM"] = np.nan
            missing_cm = merged["CM"].isna() & merged_keys.isin(freq_part.index)
            merged.loc[missing_cm, "CM"] = freq_part.loc[merged_keys.loc[missing_cm], "CM"].to_numpy()
        if "MAF" in freq_part.columns:
            if "MAF" not in merged.columns:
                merged["MAF"] = np.nan
            missing_maf = merged["MAF"].isna() & merged_keys.isin(freq_part.index)
            merged.loc[missing_maf, "MAF"] = freq_part.loc[merged_keys.loc[missing_maf], "MAF"].to_numpy()
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
    keep = metadata["MAF"] > maf_min
    removed = int((~keep).sum())
    if removed:
        LOGGER.info(f"Removed {removed} SNPs with MAF <= {maf_min} in {context}.")
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


# Parquet R2 adapter.
@dataclass(frozen=True)
class _DecodedR2RowGroup:
    """Decoded index-format parquet row group as retained-matrix index arrays."""

    row_group_index: int
    i: np.ndarray
    j: np.ndarray
    r2: np.ndarray


class _RowGroupLRUCache:
    """Small LRU cache for decoded index-format parquet R2 row groups."""

    def __init__(self, capacity: int) -> None:
        if capacity <= 0:
            raise LDSCInternalError(
                f"LD-score row-group cache setup failed: capacity={capacity} is not positive. "
                "Most likely automatic cache sizing produced an invalid value. Re-run with "
                "DEBUG logging and report the traceback."
            )
        self.capacity = int(capacity)
        self._entries: OrderedDict[int, _DecodedR2RowGroup] = OrderedDict()
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self.row_group_reads = 0

    def get(self, row_group_index: int) -> _DecodedR2RowGroup | None:
        """Return a cached decoded row group, updating LRU counters."""
        key = int(row_group_index)
        entry = self._entries.get(key)
        if entry is None:
            self.misses += 1
            return None
        self.hits += 1
        self._entries.move_to_end(key)
        return entry

    def put(self, entry: _DecodedR2RowGroup) -> None:
        """Insert ``entry`` and evict the least recently used group if needed."""
        key = int(entry.row_group_index)
        if key in self._entries:
            self._entries[key] = entry
            self._entries.move_to_end(key)
            return
        self._entries[key] = entry
        while len(self._entries) > self.capacity:
            self._entries.popitem(last=False)
            self.evictions += 1


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
        self._rg_bounds: list[tuple[int, int, int]] = []
        self._row_group_cache: _RowGroupLRUCache | None = None
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

        idx1_col = self._pf.schema_arrow.names.index("IDX_1")
        self._rg_bounds = []
        for rg_index in range(meta.num_row_groups):
            stats = meta.row_group(rg_index).column(idx1_col).statistics
            if stats is None or not getattr(stats, "has_min_max", False):
                continue
            self._rg_bounds.append((int(stats.min), int(stats.max), rg_index))

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
                    "ldscore cannot apply raw R2 bias correction without `--r2-sample-size`. "
                    "Most likely `--r2-bias-mode raw` was selected without the sample size "
                    "used to estimate R2. Pass `--r2-sample-size <N>`, or use "
                    "`--r2-bias-mode unbiased` for pre-corrected R2 values."
                )
            denom = self.r2_sample_size - 2
            if denom <= 0:
                raise LDSCConfigError(
                    f"ldscore received invalid `--r2-sample-size={self.r2_sample_size}` "
                    "for raw R2 correction. Most likely the sample size is too small "
                    "for the unbiased correction denominator. Pass a sample size greater "
                    "than 2, or use `--r2-bias-mode unbiased`."
                )
            values = values - (1.0 - values) / denom
            # Share the writer's R2<=1 invariant: roundoff/raw inputs can exceed 1.
            values = np.minimum(values, np.float32(1.0))
        return values

    def _row_group_indices_for_index_window(self, start: int, stop: int) -> list[int]:
        """Return row groups whose IDX_1 footer bounds overlap a retained-SNP index window."""
        start = max(0, int(start))
        stop = min(int(stop), int(getattr(self, "m", 0)))
        if stop <= start or not self._rg_bounds:
            return []
        lo = int(self._retained_build_idx[start])
        hi = int(self._retained_build_idx[stop - 1])
        return [rg for mn, mx, rg in self._rg_bounds if mn <= hi and mx >= lo]

    @staticmethod
    def _sliding_query_index_windows(block_left: np.ndarray, snp_batch_size: int, m: int) -> list[tuple[int, int]]:
        """Mirror parquet LD-score matrix query windows for cache sizing."""
        if m <= 0:
            return []
        snp_batch_size = int(snp_batch_size)
        if snp_batch_size <= 0:
            raise LDSCConfigError(
                f"ldscore received invalid snp_batch_size={snp_batch_size}. Most likely "
                "the parquet query batch size was set to zero or a negative value. Pass "
                "a positive integer batch size."
            )

        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / snp_batch_size) * snp_batch_size
        windows: list[tuple[int, int]] = []

        positive = np.nonzero(block_left > 0)
        if np.any(positive):
            block_width = int(positive[0][0])
        else:
            block_width = m
        block_width = int(np.ceil(block_width / snp_batch_size) * snp_batch_size)
        if block_width > m:
            snp_batch_size = 1
            block_width = m

        l_A = 0
        for l_B in range(0, block_width, snp_batch_size):
            chunk_width = min(snp_batch_size, m - l_B)
            if chunk_width <= 0:
                continue
            windows.append((min(l_A, l_B), max(l_A + block_width, l_B + chunk_width)))

        b0 = block_width
        md = int(snp_batch_size * np.floor(m / snp_batch_size))
        end = md + 1 if md != m else md
        for l_B in range(b0, end, snp_batch_size):
            old_block_width = block_width
            block_width = int(block_sizes[l_B])
            if l_B > b0 and block_width > 0:
                l_A += old_block_width - block_width + snp_batch_size
            elif l_B == b0 and block_width > 0:
                l_A = b0 - block_width
            elif block_width == 0:
                l_A = l_B

            chunk_width = snp_batch_size
            if l_B == md:
                chunk_width = m - md
            if chunk_width <= 0:
                continue
            if block_width > 0:
                windows.append((min(l_A, l_B), max(l_A + block_width, l_B + chunk_width)))
            windows.append((l_B, l_B + chunk_width))

        return windows

    def configure_auto_row_group_cache(self, block_left: np.ndarray, snp_batch_size: int) -> None:
        """Size the decoded row-group cache from this chromosome's sliding windows."""
        if getattr(self, "_runtime_layout", None) != "index":
            self._row_group_cache = None
            return
        num_row_groups = len(self._rg_bounds)
        if num_row_groups == 0:
            self._row_group_cache = None
            LOGGER.debug(f"Chromosome {self.chrom} row-group cache disabled: no footer bounds available.")
            return

        m = int(getattr(self, "m", 0))
        query_sets = [
            set(self._row_group_indices_for_index_window(start, stop))
            for start, stop in self._sliding_query_index_windows(block_left, snp_batch_size, m)
        ]
        query_sets = [rg_set for rg_set in query_sets if rg_set]
        if not query_sets:
            self._row_group_cache = None
            LOGGER.debug(f"Chromosome {self.chrom} row-group cache disabled: no row groups overlap LD windows.")
            return

        if len(query_sets) == 1:
            max_adjacent_union = len(query_sets[0])
        else:
            max_adjacent_union = max(len(left | right) for left, right in zip(query_sets, query_sets[1:]))
        capacity = min(max_adjacent_union + 1, num_row_groups)
        self._row_group_cache = _RowGroupLRUCache(capacity)
        LOGGER.debug(
            f"Chromosome {self.chrom} row-group cache capacity={capacity}, "
            f"row_groups={num_row_groups}, query_windows={len(query_sets)}."
        )

    def log_row_group_cache_summary(self) -> None:
        """Log DEBUG-only cache diagnostics after a chromosome finishes."""
        cache = self._row_group_cache
        if cache is None:
            return
        LOGGER.debug(
            f"Chromosome {self.chrom} row-group cache summary: capacity={cache.capacity}, "
            f"hits={cache.hits}, misses={cache.misses}, evictions={cache.evictions}, "
            f"row_group_reads={cache.row_group_reads}."
        )

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

    def _get_decoded_row_group(self, row_group_index: int) -> _DecodedR2RowGroup:
        """Return a decoded index row group from cache or parquet."""
        cache = self._row_group_cache
        if cache is not None:
            cached = cache.get(row_group_index)
            if cached is not None:
                return cached

        decoded = self._decode_index_row_group(row_group_index)
        if cache is not None:
            cache.row_group_reads += 1
            cache.put(decoded)
        return decoded

    def _query_union_arrays(self, start: int, stop: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(i, j, r2)`` numpy arrays for retained-index window ``[start, stop)``.

        Each selected decoded row group has ``i`` sorted ascending (the remap
        preserves build order), so the window's rows are located by binary
        search on ``i`` and only the slice with ``j`` in ``[start, stop)`` is
        kept. Canonical artifacts are duplicate-free (parquet pair uniqueness +
        the injective remap from ``build_index_remap``), so no deduplication is
        performed.
        """
        start = max(0, int(start))
        stop = min(int(stop), self.m)
        if stop <= start:
            return (
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float32),
            )

        parts_i: list[np.ndarray] = []
        parts_j: list[np.ndarray] = []
        parts_r2: list[np.ndarray] = []
        for index in self._row_group_indices_for_index_window(start, stop):
            group = self._get_decoded_row_group(index)
            if group.i.size == 0:
                continue
            lo = int(np.searchsorted(group.i, start, side="left"))
            hi = int(np.searchsorted(group.i, stop, side="left"))
            if hi <= lo:
                continue
            gj = group.j[lo:hi]
            keep = (gj >= start) & (gj < stop)
            if not keep.any():
                continue
            parts_i.append(group.i[lo:hi][keep])
            parts_j.append(gj[keep])
            parts_r2.append(group.r2[lo:hi][keep])

        if not parts_i:
            return (
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.int32),
                np.empty(0, dtype=np.float32),
            )
        return (
            np.concatenate(parts_i),
            np.concatenate(parts_j),
            np.concatenate(parts_r2),
        )

    def cross_block_matrix(self, l_A: int, b: int, l_B: int, c: int) -> np.ndarray:
        """Build the dense cross-block R2 matrix used by the parquet backend."""
        if b <= 0 or c <= 0:
            return np.zeros((max(b, 0), max(c, 0)), dtype=np.float32)

        a_start, a_stop = l_A, l_A + b
        b_start, b_stop = l_B, l_B + c
        i, j, values = self._query_union_arrays(min(a_start, b_start), max(a_stop, b_stop))

        matrix = np.zeros((b, c), dtype=np.float32)
        if i.size > 0:
            left_i = (a_start <= i) & (i < a_stop) & (b_start <= j) & (j < b_stop)
            if left_i.any():
                matrix[i[left_i] - a_start, j[left_i] - b_start] = values[left_i]

            left_j = (a_start <= j) & (j < a_stop) & (b_start <= i) & (i < b_stop)
            if left_j.any():
                matrix[j[left_j] - a_start, i[left_j] - b_start] = values[left_j]

        overlap_start = max(a_start, b_start)
        overlap_stop = min(a_stop, b_stop)
        for idx in range(overlap_start, overlap_stop):
            matrix[idx - a_start, idx - b_start] = 1.0

        return matrix

    def within_block_matrix(self, l_B: int, c: int) -> np.ndarray:
        """Build the dense within-block R2 matrix used by the parquet backend."""
        if c <= 0:
            return np.zeros((0, 0), dtype=np.float32)

        b_start, b_stop = l_B, l_B + c
        # _query_union_arrays already restricts both endpoints to [b_start, b_stop).
        i, j, values = self._query_union_arrays(b_start, b_stop)

        matrix = np.zeros((c, c), dtype=np.float32)
        if i.size > 0:
            ii = i - b_start
            jj = j - b_start
            matrix[ii, jj] = values
            matrix[jj, ii] = values

        np.fill_diagonal(matrix, 1.0)
        return matrix


def ld_score_var_blocks_from_r2_reader(
    block_left: np.ndarray,
    snp_batch_size: int,
    annot: np.ndarray,
    block_reader: SortedR2BlockReader,
) -> np.ndarray:
    """
    Mirror the old LDSC sliding-block accumulation while sourcing block-local
    dense R2 matrices from the index-format parquet reader instead of genotype
    blocks. The reader configures a chromosome-local decoded row-group cache
    from ``block_left`` and ``snp_batch_size`` before traversal.
    """
    m = annot.shape[0]
    n_a = annot.shape[1]
    if snp_batch_size <= 0:
        raise LDSCConfigError(
            f"ldscore received invalid snp_batch_size={snp_batch_size}. Most likely "
            "the parquet query batch size was set to zero or a negative value. Pass "
            "a positive integer batch size."
        )
    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / snp_batch_size) * snp_batch_size
    cor_sum = np.zeros((m, n_a), dtype=np.float64)
    block_reader.configure_auto_row_group_cache(block_left, snp_batch_size)

    try:
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b / snp_batch_size) * snp_batch_size)
        if b > m:
            snp_batch_size = 1
            b = m

        l_A = 0
        for l_B in range(0, b, snp_batch_size):
            chunk_width = min(snp_batch_size, m - l_B)
            if chunk_width <= 0:
                continue
            rfuncAB = block_reader.cross_block_matrix(l_A, b, l_B, chunk_width)
            cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + chunk_width, :])

        b0 = b
        md = int(snp_batch_size * np.floor(m / snp_batch_size))
        end = md + 1 if md != m else md
        for l_B in range(b0, end, snp_batch_size):
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                l_A += old_b - b + snp_batch_size
            elif l_B == b0 and b > 0:
                l_A = b0 - b
            elif b == 0:
                l_A = l_B

            chunk_width = snp_batch_size
            if l_B == md:
                chunk_width = m - md
            if chunk_width <= 0:
                continue

            p1 = np.all(annot[l_A:l_A + b, :] == 0)
            p2 = np.all(annot[l_B:l_B + chunk_width, :] == 0)
            if p1 and p2:
                continue

            rfuncAB = block_reader.cross_block_matrix(l_A, b, l_B, chunk_width)
            cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + chunk_width, :])
            cor_sum[l_B:l_B + chunk_width, :] += np.dot(annot[l_A:l_A + b, :].T, rfuncAB).T

            rfuncBB = block_reader.within_block_matrix(l_B, chunk_width)
            cor_sum[l_B:l_B + chunk_width, :] += np.dot(rfuncBB, annot[l_B:l_B + chunk_width, :])
    finally:
        block_reader.log_row_group_cache_summary()

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
    ``MAF >= common_maf_min`` when MAF metadata is available.
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
    coords, max_dist = build_window_coordinates(metadata, args)
    block_left = get_block_lefts(
        coords,
        max_dist,
    )
    check_whole_chromosome_window(block_left, args, chrom)
    regression_mask = regression_mask_from_keys(metadata, regression_keys, args.snp_identifier).reshape(-1, 1)
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=True)
    combined_annot = np.c_[annot_matrix, regression_mask]
    block_reader = SortedR2BlockReader(
        paths=resolve_parquet_files(args, chrom=chrom),
        chrom=chrom,
        metadata=metadata,
        identifier_mode=args.snp_identifier,
        r2_bias_mode=args.r2_bias_mode,
        r2_sample_size=args.r2_sample_size,
        genome_build=args.genome_build,
    )
    combined_scores = ld_score_var_blocks_from_r2_reader(
        block_left=block_left,
        snp_batch_size=args.snp_batch_size,
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
            args.r2_bias_mode, args.r2_sample_size = ref_panel_mod._resolve_r2_bias_from_meta(
                args.r2_bias_mode,
                getattr(args, "r2_sample_size", None),
                stored,
            )
        if args.r2_bias_mode is None:
            args.r2_bias_mode = "unbiased"
        if args.r2_bias_mode == "raw" and args.r2_sample_size is None:
            raise LDSCUsageError(
                "ldscore cannot apply raw R2 bias correction without `--r2-sample-size`. "
                "Most likely `--r2-bias-mode raw` was selected for parquet input without "
                "the sample size used to estimate R2. Pass `--r2-sample-size <N>`, or use "
                "`--r2-bias-mode unbiased` for pre-corrected R2 values."
            )
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
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default="unbiased", help="Whether sorted parquet R2 values are raw sample r^2 or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
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
    parser.add_argument("--snp-batch-size", default=128, type=int, help="Number of SNPs processed per LD-score sliding batch.")
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
