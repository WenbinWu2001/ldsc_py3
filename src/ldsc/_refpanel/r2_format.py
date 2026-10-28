"""The R2 matrix file format: schema, writer, readers and pair lookup.

An R2 matrix stores, for each chromosome, pairwise unbiased R2 values in
``chr<chrom>_r2.parquet`` (int16-quantized, indexed by sidecar row) and
SNP metadata in the ``chr<chrom>_meta.tsv.gz`` sidecar, bound together by
an identity hash. This module writes both files, checks an R2 matrix's metadata
contract, streams its pairs for LD-score calculation and looks up
individual pairs for ``query-r2-matrix``.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import warnings
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .._core.chromosomes import normalize_chromosome
from .._core.columns import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    REFERENCE_METADATA_SPEC_MAP,
    normalize_genome_build,
    resolve_optional_column,
    resolve_required_column,
)
from .._core.identifiers import validate_retained_identifier_uniqueness
from .._core.snp_identity import (
    effective_merge_key_series,
    identity_artifact_metadata,
    is_allele_aware_mode,
    sidecar_identity_sha256,
    normalize_snp_identifier_mode,
)
from .build_kernel import PairColumns
from ..errors import LDSCConfigError, LDSCDependencyError, LDSCInputError, LDSCInternalError, LDSCUsageError
from .windows import _LDWindowSpec, _format_ld_window_spec, _requested_ld_window_spec

LOGGER = logging.getLogger("LDSC.ref_panel")


_LDSCORE_PARQUET_DOC = "docs/user-guide/troubleshooting.md#ldscore-r2-matrix-is-incompatible"


def build_index_remap(
    full_sidecar: pd.DataFrame,
    retained_metadata: pd.DataFrame,
    identifier_mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Map R2 matrix (build) indices to retained SNP indices.

    ``full_sidecar`` is the complete R2 matrix in build order (the parquet IDX
    space). ``retained_metadata`` lists the retained reference SNPs in
    matrix order. Returns ``(remap, retained_build_idx)`` where
    ``remap[build_idx]`` is the retained matrix index or ``-1``, and
    ``retained_build_idx[matrix_idx]`` is the originating build index (ascending,
    used for IDX_1 row-group pruning). Matching uses the same mode-dependent
    identity keys as the legacy per-pair decode, so the result is bit-identical.

    The remap must be injective over retained indices: two distinct R2 matrix rows
    sharing an identity key under ``identifier_mode`` (a collapse) is a hard
    error, which also keeps ``retained_build_idx`` a well-defined permutation.
    """
    mode = normalize_snp_identifier_mode(identifier_mode)
    full_keys = effective_merge_key_series(full_sidecar, mode, context="R2 matrix sidecar index remap").to_numpy()
    retained_keys = effective_merge_key_series(retained_metadata, mode, context="retained metadata index remap").to_numpy()
    retained_index = pd.Index(retained_keys)
    if retained_index.has_duplicates:
        raise LDSCInputError(
            "ldscore could not align retained reference metadata to the R2 matrix sidecar. "
            f"The retained metadata has duplicate SNP identity keys under mode '{mode}'. "
            "Most likely the SNP identifier mode is too coarse for this R2 matrix. "
            "Use an allele-aware `--snp-identifier` mode or rebuild the R2 matrix "
            "after removing duplicate SNP identities."
        )
    remap = retained_index.get_indexer(full_keys).astype(np.int32, copy=False)

    valid = remap >= 0
    valid_targets = remap[valid]
    if np.unique(valid_targets).size != valid_targets.size:
        raise LDSCInputError(
            "ldscore could not align the R2 matrix to its sidecar: distinct "
            f"R2 matrix SNPs collapse onto the same retained SNP under mode '{mode}'. "
            "Most likely the R2 matrix sidecar has duplicate SNP identities at this "
            "identifier resolution. Use an allele-aware `--snp-identifier` mode, or "
            "rebuild the R2 matrix after removing duplicate SNP identities. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )

    m = len(retained_metadata)
    retained_build_idx = np.empty(m, dtype=np.int64)
    retained_build_idx[valid_targets] = np.nonzero(valid)[0]
    return remap, retained_build_idx


def _panel_sidecar_path_for_r2(r2_path: str) -> Path:
    """Return the ``chrN_meta.tsv.gz`` sidecar path paired with an R2 Parquet file."""
    p = Path(r2_path)
    if not p.name.endswith("_r2.parquet"):
        raise LDSCInputError(
            f"ldscore could not locate the sidecar for R2 Parquet file '{r2_path}'. Most likely "
            "the parquet file does not use the canonical `chrN_r2.parquet` filename written "
            "by `ldsc build-r2-matrix`. regenerate the R2 matrix or pass the "
            "R2 matrix directory. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    return p.with_name(p.name[: -len("_r2.parquet")] + "_meta.tsv.gz")


def _load_full_panel_sidecar(r2_path: str) -> pd.DataFrame:
    """Load the complete (unrestricted) R2 matrix sidecar that defines the index space."""
    sidecar_path = _panel_sidecar_path_for_r2(r2_path)
    if not sidecar_path.exists():
        raise LDSCInputError(
            f"ldscore could not load R2 Parquet file '{r2_path}' because the "
            f"required sidecar '{sidecar_path}' is missing. Most likely the parquet file "
            "was copied without its matching `chrN_meta.tsv.gz` sidecar. Keep the R2 "
            "Parquet file and sidecar together or regenerate with `ldsc build-r2-matrix`. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    df = pd.read_csv(sidecar_path, sep="\t", comment="#")
    context = f"R2 matrix sidecar {sidecar_path}"
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
            f"ldscore could not use R2 matrix input at {context}: the sidecar has "
            f"{len(full_sidecar)} rows but the parquet records n_snps={n_snps}. Most "
            "likely the parquet and sidecar are not a matched pair. Restore the matching "
            "sidecar or regenerate the R2 matrix. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )
    actual = sidecar_identity_sha256(full_sidecar, context=context)
    if actual != identity_hash:
        raise LDSCInputError(
            f"ldscore could not use R2 matrix input at {context}: the sidecar identity "
            f"hash {actual} does not match parquet ldsc:sidecar_identity_sha256 "
            f"{identity_hash}. Most likely the sidecar is wrong, reordered, or edited. "
            "Restore the matching sidecar or regenerate the R2 matrix. "
            f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
        )


def _arrow_column_to_numpy(column):
    """Convert an Arrow column to NumPy across PyArrow array API versions."""
    try:
        return column.to_numpy(zero_copy_only=False)
    except TypeError:
        return column.to_numpy()


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
    r2_min: float


def read_r2_parquet_contract(path) -> R2ParquetContract:
    """Validate one R2 Parquet file footer against the current ``build-r2-matrix`` format.

    Readers accept exactly what the writer produces: ``IDX_1``/``IDX_2``/
    ``R2``/``SIGN_R`` columns with an int16 ``R2``, the writer's
    ``ldsc:r2_encoding`` and ``ldsc:r2_scale``, already-unbiased values
    (``ldsc:r2_bias="unbiased"``) with an ``ldsc:n_samples`` of at least 3, and
    complete binding, identity, LD-window, and ``ldsc:r2_min`` metadata.
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
            f"{summary} the format written by the current `ldsc build-r2-matrix`: {details}. "
            "Most likely the R2 matrix was built by an older LDSC3 version or its metadata was "
            "altered. Regenerate the R2 matrix with `ldsc build-r2-matrix`. "
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
            "Reading R2 matrices requires pyarrow. Most likely "
            "R2 matrix input was requested in an environment missing pyarrow. "
            "Install pyarrow or use PLINK reference-panel input instead."
        ) from exc

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
    text("artifact_type", "r2_matrix")
    snp_identifier = parsed("snp_identifier", normalize_snp_identifier_mode, lambda mode: True,
                            "an LDSC3 SNP identifier mode")
    genome_build = text("genome_build")
    window_mode = parsed("ld_window_mode", str, lambda mode: mode in {"snps", "kb", "cm"}, "snps, kb, or cm")
    window_value = parsed("ld_window_value", float, lambda value: bool(np.isfinite(value)) and value > 0,
                          "a positive number")
    r2_min = parsed("r2_min", float, lambda value: bool(np.isfinite(value)) and value <= 1,
                    "a number no greater than 1")
    if problems:
        return None, problems
    return R2ParquetContract(
        path=str(path), snp_identifier=snp_identifier, genome_build=genome_build,
        sorted_by_build=sorted_by_build, n_snps=n_snps, sidecar_identity_sha256=identity_hash,
        n_samples=n_samples, r2_scale=r2_scale, ld_window=_LDWindowSpec(window_mode, window_value),
        r2_min=r2_min,
    ), []


def _r2_matrix_ld_window_spec(parquet_paths: Sequence[str]) -> _LDWindowSpec | None:
    """Return the common recorded LD-window spec for an R2 matrix input set."""
    panel_spec: _LDWindowSpec | None = None
    for path in parquet_paths:
        current = read_r2_parquet_contract(path).ld_window
        if panel_spec is None:
            panel_spec = current
            continue
        if current != panel_spec:
            raise LDSCInputError(
                "ldscore could not validate the R2 matrix LD window because "
                "the resolved parquet files record conflicting build windows: "
                f"`{_format_ld_window_spec(panel_spec)}` and "
                f"`{_format_ld_window_spec(current)}`. Most likely files from "
                "different R2 matrix builds were mixed in one input. Use a "
                "single coherent `--r2-matrix-dir` or regenerate the R2 matrix."
            )
    return panel_spec


def _raise_ldscore_window_exceeds_panel_window(
    *,
    requested: _LDWindowSpec,
    recorded: _LDWindowSpec,
    chrom: str,
) -> None:
    """Raise the user-facing parquet-window mismatch error."""
    raise LDSCUsageError(
        f"ldscore cannot compute chromosome {chrom} from the "
        f"R2 matrix: the user-requested LD window `{_format_ld_window_spec(requested)}` "
        f"is wider than the input R2 matrix window "
        f"`{_format_ld_window_spec(recorded)}`. Most likely the "
        "R2 matrix was built with a smaller LD window, so it does not contain all "
        "SNP pairs required by this ldscore run. Rebuild the R2 matrix "
        f"with `{_format_ld_window_spec(requested)}` or a wider build window, or "
        f"rerun ldscore with an LD window no wider than the recorded R2 matrix window "
        f"`{_format_ld_window_spec(recorded)}`."
    )


def validate_ldscore_window_within_r2_matrix_window(
    args: argparse.Namespace,
    *,
    parquet_paths: Sequence[str],
    chrom: str,
) -> None:
    """Reject ldscore windows the R2 matrix is not known to cover.

    The rule uses only the recorded parquet LD-window mode and value, so it is
    decided before any computation: a window in another unit is rejected, and
    a same-unit window must be no wider than the R2 matrix's. Coverage is not
    inferred from per-SNP coordinates. For ``--ld-wind-snps`` this relies on
    counting R2 matrix rows (see :func:`build_window_coordinates`).
    """
    panel_spec = _r2_matrix_ld_window_spec(parquet_paths)
    if panel_spec is None:
        return
    requested_spec = _requested_ld_window_spec(args)

    if requested_spec.mode != panel_spec.mode:
        requested, recorded = _format_ld_window_spec(requested_spec), _format_ld_window_spec(panel_spec)
        raise LDSCUsageError(
            f"ldscore cannot use the R2 matrix for chromosome {chrom} with "
            f"`{requested}`: the R2 matrix was built with `{recorded}`, and windows in different units "
            "cannot be checked for coverage from their settings, so the run could silently miss "
            f"SNP pairs the R2 matrix never stored. Rerun with the R2 matrix's window `{recorded}` (or a "
            f"narrower `{panel_spec.flag}` value), or rebuild the R2 matrix with "
            f"`ldsc build-r2-matrix {requested}`."
        )

    if requested_spec.value > panel_spec.value:
        _raise_ldscore_window_exceeds_panel_window(
            requested=requested_spec,
            recorded=panel_spec,
            chrom=chrom,
        )


def reject_truncated_r2_matrix(parquet_paths: Sequence[str], chrom: str) -> None:
    """Reject an R2 matrix that omitted pairs below a positive ``ldsc:r2_min``.

    LD scores sum every SNP pair in the window, so pairs the R2 matrix did not
    store would count as zero. Such R2 matrices stay valid for ``query-r2-matrix``, which
    reports unstored pairs as absent.
    """
    for path in parquet_paths:
        r2_min = read_r2_parquet_contract(path).r2_min
        if r2_min > 0:
            raise LDSCInputError(
                f"ldscore cannot use R2 Parquet file '{path}' (chromosome {chrom}): the R2 matrix omits SNP pairs with "
                f"unbiased R2 below {r2_min:g} (ldsc:r2_min, set by `build-r2-matrix --r2-min`). LD scores sum "
                "every SNP pair in the window, so the omitted pairs would count as zero and change the scores. "
                "Rebuild the R2 matrix without `--r2-min` (the default keeps every pair). R2 matrices built with a "
                "positive `--r2-min` remain usable for `ldsc query-r2-matrix`."
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
                "Most likely the R2 matrix directory is missing this chromosome or the path token "
                "does not resolve to `chrN_r2.parquet`. Pass the correct `--r2-matrix-dir` or "
                "regenerate the R2 matrix. "
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
                "ldscore could not read R2 matrix input because pyarrow is not installed. "
                "Most likely R2 matrix input was requested in an environment missing pyarrow. "
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
                "but an R2 matrix requires exactly one file per chromosome. Most likely "
                "the R2 matrix path/glob is too broad or the directory contains duplicate chromosome "
                "artifacts. Narrow `--r2-matrix-dir` or remove duplicate parquet files."
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
        """R2 matrix sidecar row of each retained SNP, in retained order."""
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
                "the R2 matrix and annotation inputs use different genome builds. "
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
                    f"R2 matrix '{path}' stores {meta.num_rows} SNP-pair rows in "
                    f"{meta.num_row_groups} large file blocks (average {avg:.0f} rows per block). "
                    "Reading continues, but queries may be slower and require more memory. "
                    "For repeated use, rebuild the R2 matrix with the current `ldsc build-r2-matrix`, "
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

        R2 matrix indices are remapped to the retained reference SNPs before
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


# Text format of float columns (CM, MAF) in every reference-metadata sidecar.
SIDECAR_FLOAT_FORMAT = "%.9g"


def round_to_sidecar_precision(values) -> np.ndarray:
    """Return ``values`` exactly as a sidecar stores them and readers parse them back.

    Build-time decisions such as the cM window use these values, so the pairs a
    R2 matrix stores match the window ``ldscore`` later rebuilds from the sidecar.
    """
    return np.char.mod(SIDECAR_FLOAT_FORMAT, np.asarray(values, dtype=float)).astype(float)


# int16 symmetric R2 quantization (refactor #11). Scale 32767 = int16 max, so
# the endpoint is exact: round(1.0 * 32767) = 32767, decoded 32767/32767 = 1.0.
R2_QUANT_SCALE = 32767


R2_ENCODING = "int16_symmetric"


def _quantize_r2(r2: np.ndarray) -> np.ndarray:
    """Symmetric int16 quantization of unbiased R2 (scale ``R2_QUANT_SCALE``).

    ``stored = clip(round(R2 * 32767), -32767, 32767)``. The clip is defensive:
    with the upstream ``min(R2, 1.0)`` and ``|negatives| ~ 1/N`` it never fires
    on real data. float64 multiply makes rounding deterministic at half-steps.
    """
    scaled = np.rint(r2.astype(np.float64) * R2_QUANT_SCALE)
    return np.clip(scaled, -R2_QUANT_SCALE, R2_QUANT_SCALE).astype(np.int16)


def build_runtime_metadata_table(
    *,
    metadata: pd.DataFrame,
    positions: np.ndarray,
    cm_values: np.ndarray | None,
) -> pd.DataFrame:
    """Build the LDSC3 runtime metadata sidecar for one emitted build.

    When ``cm_values`` is ``None``, the sidecar keeps the SNP rows and writes a
    nullable ``CM`` column. This is used by SNP- and kb-window builds that omit
    genetic maps.
    """
    cm_column = (
        pd.Series(pd.array([pd.NA] * len(metadata), dtype="Float64"))
        if cm_values is None
        else pd.Series(np.asarray(cm_values, dtype=float))
    )

    values = {
        "CHR": metadata["CHR"].map(normalize_chromosome).astype(str),
        "POS": np.asarray(positions, dtype=np.int64),
        "SNP": metadata["SNP"].astype(str),
        "A1": metadata["A1"].astype(str),
        "A2": metadata["A2"].astype(str),
        "CM": cm_column,
        "MAF": pd.to_numeric(metadata["MAF"], errors="coerce").astype(float),
    }
    return pd.DataFrame(values).reset_index(drop=True)


def _ensure_parent_dir(path: str | PathLike[str]) -> None:
    """Create the parent directory for ``path`` if it does not already exist."""
    parent = Path(path).parent
    parent.mkdir(parents=True, exist_ok=True)


def _standard_r2_index_table(pa, schema, *, i, j, r2, sign):
    """Build one 4-column index ``pyarrow.Table`` from concatenated pair columns.

    ``i``/``j`` are sidecar-row indices (stored int32), ``r2`` is unbiased float32
    stored as symmetric int16 (scale ``R2_QUANT_SCALE``), and ``sign`` is int8
    ``+1``/``-1`` stored as the ``SIGN_R`` bool (``True`` when the correlation r >= 0).
    The indices are stored as-is; no reference-SNP join or identifier expansion.
    """
    if i.size == 0:
        return schema.empty_table()
    return pa.Table.from_arrays(
        [
            pa.array(i.astype(np.int32), type=pa.int32()),
            pa.array(j.astype(np.int32), type=pa.int32()),
            pa.array(_quantize_r2(r2), type=pa.int16()),
            pa.array(sign == 1, type=pa.bool_()),
        ],
        schema=schema,
    )


def write_r2_parquet(
    *,
    pair_chunks: Iterable[PairColumns],
    path: str | PathLike[str],
    genome_build: str,
    n_samples: int,
    snp_identifier: str,
    n_snps: int,
    sidecar_identity_sha256: str,
    batch_size: int = 100_000,
    row_group_size: int = 50_000,
    ld_window_mode: str,
    ld_window_value: float,
    r2_min: float = 0.0,
) -> str:
    """
    Write one canonical 4-column R2 Parquet file with LDSC3 schema metadata.

    Columns are ``IDX_1``/``IDX_2`` (int32 sidecar-row indices), ``R2``
    (int16 symmetric quantization, scale ``R2_QUANT_SCALE``=32767, with
    ``BYTE_STREAM_SPLIT`` encoding; the endpoint ``1.0`` maps to exactly 32767),
    and ``SIGN_R`` (bit-packed bool, ``True`` when the Pearson correlation r >= 0).
    Pairs must arrive in non-decreasing ``IDX_1`` order because row-group pruning
    depends on monotonic footer statistics.

    The parquet is bound to its per-SNP sidecar by ``ldsc:n_snps`` (the index
    space size) and ``ldsc:sidecar_identity_sha256`` (a hash of the sidecar's
    ``CHR:POS:A1:A2`` identity). The stored indices are meaningless without the
    matching sidecar, so the reader hard-fails on a binding mismatch.

    ``r2_min`` records the unbiased-R2 floor applied upstream during pair
    emission; it does not itself filter rows here. The table is
    pairwise-complete only when ``ldsc:r2_min`` is non-positive; ldscore
    rejects any other value, and query-r2-matrix reports absent pairs as absent.

    ``ld_window_mode`` and ``ld_window_value`` record the build-time LD window
    that bounded emitted pairs. LD-score calculation uses these fields to reject
    wider requested windows that would require pairs the parquet does not store.

    Column chunks are compressed with zstd level 9 when the codec is available
    (falling back to snappy otherwise). Compression is recorded per column-chunk
    in the footer and auto-detected on read.
    """

    _ensure_parent_dir(path)
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:
        raise LDSCDependencyError(
            "build-r2-matrix could not write R2 Parquet files because pyarrow is not installed. "
            "Most likely R2 matrix generation was requested in an environment missing pyarrow. "
            "Install pyarrow and rerun build-r2-matrix."
        ) from exc

    if int(n_snps) >= 2**31:
        raise LDSCInputError(
            f"build-r2-matrix cannot write an R2 Parquet file with n_snps={n_snps}: the index "
            "space exceeds the int32 IDX column range. Most likely the retained reference "
            "R2 matrix is too large for the current parquet format. Split the R2 matrix or reduce "
            "the retained reference SNPs."
        )
    if ld_window_mode is None or ld_window_value is None:
        raise LDSCInternalError(
            "build-r2-matrix could not write R2 Parquet file LD-window metadata because "
            "ld_window_mode and ld_window_value were not both provided. Most likely "
            "the R2 matrix builder passed an incomplete build-window record. "
            "Re-run with DEBUG logging and report the traceback."
        )

    pa_meta = {
        **{
            f"ldsc:{key}".encode("utf-8"): str(value).encode("utf-8")
            for key, value in identity_artifact_metadata(
                artifact_type="r2_matrix",
                snp_identifier=snp_identifier,
                genome_build=genome_build,
            ).items()
        },
        b"ldsc:sorted_by_build": genome_build.encode("utf-8"),
        b"ldsc:row_group_size": str(row_group_size).encode("utf-8"),
        b"ldsc:n_samples": str(n_samples).encode("utf-8"),
        b"ldsc:r2_bias": b"unbiased",
        b"ldsc:r2_min": str(r2_min).encode("utf-8"),
        b"ldsc:n_snps": str(int(n_snps)).encode("utf-8"),
        b"ldsc:sidecar_identity_sha256": str(sidecar_identity_sha256).encode("utf-8"),
        b"ldsc:r2_encoding": R2_ENCODING.encode("utf-8"),
        b"ldsc:r2_scale": str(R2_QUANT_SCALE).encode("utf-8"),
        b"ldsc:ld_window_mode": str(ld_window_mode).encode("utf-8"),
        b"ldsc:ld_window_value": str(ld_window_value).encode("utf-8"),
    }
    schema = pa.schema(
        [
            ("IDX_1", pa.int32()),
            ("IDX_2", pa.int32()),
            ("R2", pa.int16()),
            ("SIGN_R", pa.bool_()),
        ]
    ).with_metadata(pa_meta)
    writer = None
    buf_i: list[np.ndarray] = []
    buf_j: list[np.ndarray] = []
    buf_r2: list[np.ndarray] = []
    buf_sign: list[np.ndarray] = []
    buffered = 0
    prev_idx1: int | None = None

    def _flush() -> None:
        nonlocal writer, prev_idx1, buf_i, buf_j, buf_r2, buf_sign, buffered
        i = np.concatenate(buf_i) if buf_i else np.empty(0, dtype=np.int64)
        j = np.concatenate(buf_j) if buf_j else np.empty(0, dtype=np.int64)
        r2 = np.concatenate(buf_r2) if buf_r2 else np.empty(0, dtype=np.float32)
        sign = np.concatenate(buf_sign) if buf_sign else np.empty(0, dtype=np.int8)
        table = _standard_r2_index_table(pa, schema, i=i, j=j, r2=r2, sign=sign)
        if i.size:
            sequence = i if prev_idx1 is None else np.concatenate(([prev_idx1], i))
            bad = np.flatnonzero(np.diff(sequence) < 0)
            if bad.size:
                k = int(bad[0])
                raise LDSCInternalError(
                    "build-r2-matrix R2 Parquet file writing failed: pair rows arrived out of "
                    "non-decreasing IDX_1 order. "
                    f"Received IDX_1={int(sequence[k + 1])} after IDX_1={int(sequence[k])}. "
                    "Most likely retained SNP sorting or pair emission is inconsistent. "
                    "Re-run with DEBUG logging and report the traceback."
                )
            prev_idx1 = int(i[-1])
        if writer is None:
            # IDX_2: DELTA_BINARY_PACKED exploits sorted right-neighbors (median
            # forward-gap = 1), collapsing ~650 MB to ~6 MB vs the default
            # dictionary path. R2: BYTE_STREAM_SPLIT splits the int16 byte planes
            # so zstd compresses the low-entropy high byte (R2 clusters near 0)
            # apart from the noisy low byte. use_dictionary restricts auto-
            # dictionary to IDX_1 only so PyArrow cannot override either encoding.
            enc_kwargs = dict(
                column_encoding={"IDX_2": "DELTA_BINARY_PACKED", "R2": "BYTE_STREAM_SPLIT"},
                use_dictionary=["IDX_1"],
            )
            if pa.Codec.is_available("zstd"):
                writer = pq.ParquetWriter(
                    str(path), schema, compression="zstd", compression_level=9, **enc_kwargs
                )
            else:
                warnings.warn(
                    "zstd codec is not available in this pyarrow build; "
                    "falling back to snappy compression. "
                    "Install pyarrow from conda-forge or PyPI to enable zstd.",
                    UserWarning,
                    stacklevel=2,
                )
                writer = pq.ParquetWriter(str(path), schema, compression="snappy", **enc_kwargs)
        writer.write_table(table, row_group_size=row_group_size)
        buf_i, buf_j, buf_r2, buf_sign = [], [], [], []
        buffered = 0

    try:
        for chunk in pair_chunks:
            ci, cj, cr2, csign = chunk
            if ci.size == 0:
                continue
            buf_i.append(ci)
            buf_j.append(cj)
            buf_r2.append(cr2)
            buf_sign.append(csign)
            buffered += int(ci.size)
            if buffered >= batch_size:
                _flush()
        if buffered or writer is None:
            _flush()
    finally:
        if writer is not None:
            writer.close()
    return str(path)


def write_runtime_metadata_sidecar(
    df: pd.DataFrame,
    path: str | PathLike[str],
    *,
    genome_build: str,
    snp_identifier: str,
) -> str:
    """Write one gzip-compressed LDSC3 runtime metadata sidecar."""

    _ensure_parent_dir(path)
    identity_metadata = identity_artifact_metadata(
        artifact_type="ref_panel_metadata",
        snp_identifier=snp_identifier,
        genome_build=genome_build,
    )
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for key, value in identity_metadata.items():
            handle.write(f"# ldsc:{key}={value}\n")
        df.to_csv(handle, sep="\t", index=False, na_rep="NA", float_format=SIDECAR_FLOAT_FORMAT)
    return str(path)


# Auto-strategy crossover: at or below this many queried pairs, prune+random-access
# row groups; above it, scan every row group sequentially (a large query touches
# most groups anyway, so the linear scan beats per-group seek + pruning overhead).
_RANDOM_STRATEGY_MAX_PAIRS = 50_000


def _arrow_to_numpy(column):
    """Convert an Arrow column to NumPy across PyArrow versions."""
    try:
        return column.to_numpy(zero_copy_only=False)
    except TypeError:
        return column.to_numpy()


def _idx1_column_index(parquet_file) -> int:
    """Return the parquet column index of ``IDX_1`` (for footer statistics)."""
    return parquet_file.schema_arrow.names.index("IDX_1")


def _prune_row_groups(parquet_file, left_sorted: np.ndarray, idx1_col: int):
    """Return row groups whose ``IDX_1`` [min,max] covers a needed left index.

    ``left_sorted`` is the sorted unique set of canonical left indices queried.
    A row group qualifies when at least one queried left index falls within its
    ``IDX_1`` footer range. Groups without usable statistics are kept (safe).
    """
    md = parquet_file.metadata
    keep: list[int] = []
    for rg in range(md.num_row_groups):
        stats = md.row_group(rg).column(idx1_col).statistics
        if stats is None or not stats.has_min_max:
            keep.append(rg)
            continue
        lo, hi = int(stats.min), int(stats.max)
        pos = int(np.searchsorted(left_sorted, lo))
        if pos < left_sorted.shape[0] and int(left_sorted[pos]) <= hi:
            keep.append(rg)
    return keep


def lookup_pairs_in_parquet(
    parquet_file,
    i: np.ndarray,
    j: np.ndarray,
    *,
    n_snps: int,
    r2_scale: float,
):
    """Return ``(r2, sign)`` for canonical pairs ``(i, j)`` from one parquet.

    ``r2`` is float32 with ``NaN`` where the pair is not stored; ``sign`` is int8
    with ``+1``/``-1`` for stored pairs (R2 matrix orientation) and ``0`` where the
    pair is not stored. The int16 R2 column is divided by ``r2_scale``
    (dequantization).
    """
    i = np.asarray(i, dtype=np.int64)
    j = np.asarray(j, dtype=np.int64)
    n = i.shape[0]
    r2_out = np.full(n, np.nan, dtype=np.float32)
    sign_out = np.zeros(n, dtype=np.int8)
    if n == 0:
        return r2_out, sign_out

    n_snps = np.int64(n_snps)
    qkey = i * n_snps + j
    uniq_key = np.unique(qkey)
    uniq_r2 = np.full(uniq_key.shape[0], np.nan, dtype=np.float32)
    uniq_sign = np.zeros(uniq_key.shape[0], dtype=np.int8)

    # Auto strategy: prune+random-access row groups for small queries, full
    # sequential scan once the query is large enough to touch most groups anyway.
    chosen = "random" if n <= _RANDOM_STRATEGY_MAX_PAIRS else "stream"

    idx1_col = _idx1_column_index(parquet_file)
    if chosen == "random":
        row_groups = _prune_row_groups(parquet_file, np.unique(i), idx1_col)
    else:
        row_groups = range(parquet_file.metadata.num_row_groups)

    for rg in row_groups:
        table = parquet_file.read_row_group(int(rg), columns=["IDX_1", "IDX_2", "R2", "SIGN_R"])
        gi = _arrow_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        gj = _arrow_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        gkey = gi * n_snps + gj
        pos = np.searchsorted(uniq_key, gkey)
        in_range = pos < uniq_key.shape[0]
        matched = np.zeros(gkey.shape[0], dtype=bool)
        matched[in_range] = uniq_key[pos[in_range]] == gkey[in_range]
        if not matched.any():
            continue
        mpos = pos[matched]
        r2_vals = _arrow_to_numpy(table.column("R2")).astype(np.float32, copy=False)[matched]
        r2_vals = r2_vals / np.float32(r2_scale)
        sign_vals = _arrow_to_numpy(table.column("SIGN_R")).astype(bool, copy=False)[matched]
        uniq_r2[mpos] = r2_vals
        uniq_sign[mpos] = np.where(sign_vals, np.int8(1), np.int8(-1))

    qpos = np.searchsorted(uniq_key, qkey)
    r2_out[:] = uniq_r2[qpos]
    sign_out[:] = uniq_sign[qpos]
    return r2_out, sign_out
