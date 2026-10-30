"""The R2 matrix file format: names, schema, writer, readers and pair lookup.

An R2 matrix stores, for each chromosome, pairwise unbiased R2 values in
``chr<chrom>_r2.parquet`` (int16-quantized, indexed by sidecar row) and
SNP metadata in the ``chr<chrom>_meta.tsv.gz`` sidecar, bound together by
an identity hash. This module owns the file names, writes both files, opens a
chromosome of an existing R2 matrix through :class:`R2MatrixChromosome` for
``ldscore`` and ``query-r2-matrix``, streams its pairs for LD-score
calculation and looks up individual pairs.
"""

from __future__ import annotations

import gzip
import re
import warnings
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .._core.chromosomes import chrom_sort_key, normalize_chromosome
from .._core.columns import normalize_genome_build
from .._core.identifiers import validate_retained_identifier_uniqueness
from .._core.snp_identity import (
    effective_merge_key_series,
    identity_artifact_metadata,
    sidecar_identity_sha256,
    normalize_snp_identifier_mode,
)
from .build_kernel import PairColumns
from ..errors import LDSCConfigError, LDSCInputError, LDSCInternalError, LDSCUsageError
from .windows import _LDWindowSpec, _format_ld_window_spec, _requested_ld_window_spec

# Troubleshooting section for R2 matrix errors, by the command that reads the R2 matrix.
R2_MATRIX_TROUBLESHOOTING = {
    "ldscore": "docs/user-guide/troubleshooting.md#ldscore-r2-matrix-is-incompatible",
    "query-r2-matrix": "docs/user-guide/troubleshooting.md#query-r2-matrix-r2-matrix-is-incompatible",
}
_LDSCORE_PARQUET_DOC = R2_MATRIX_TROUBLESHOOTING["ldscore"]

# Sidecar columns, in the order the writer emits them and the reader requires.
SIDECAR_COLUMNS = ("CHR", "POS", "SNP", "A1", "A2", "CM", "MAF")

# File names of one chromosome: chr<chrom>_r2.parquet and chr<chrom>_meta.tsv.gz.
R2_MATRIX_FILE_GLOBS = ("chr*_r2.parquet", "chr*_meta.tsv.gz")
_R2_MATRIX_FILE_RE = re.compile(r"^chr(?P<chrom>[^_]+)_(?:r2\.parquet|meta\.tsv\.gz)$")


def r2_matrix_file_paths(build_dir: str | PathLike[str], chrom: str) -> tuple[Path, Path]:
    """Return the R2 Parquet file and sidecar paths of one chromosome in an R2 matrix build directory."""
    chrom = normalize_chromosome(chrom)
    build_dir = Path(build_dir)
    return build_dir / f"chr{chrom}_r2.parquet", build_dir / f"chr{chrom}_meta.tsv.gz"


def r2_matrix_chromosomes(build_dir: str | PathLike[str]) -> list[str]:
    """Return the chromosomes with an R2 Parquet file or a sidecar in one build directory, in genomic order.

    A chromosome with only one of its two files is listed, so that opening it
    reports the missing file instead of silently skipping the chromosome.
    """
    build_dir = Path(build_dir)
    labels = {match.group("chrom") for path in build_dir.iterdir()
              if (match := _R2_MATRIX_FILE_RE.match(path.name)) is not None}
    return sorted({normalize_chromosome(label, context=str(build_dir)) for label in labels}, key=chrom_sort_key)


def resolve_r2_build_dir(r2_matrix_dir: str | PathLike[str], genome_build: str | None) -> Path:
    """Resolve the build directory holding the chromosome files of an R2 matrix.

    ``r2_matrix_dir`` is either a build directory or the ``build-r2-matrix``
    output root with ``hg19/`` and ``hg38/`` children; ``genome_build`` selects
    a child when the root has both.
    """
    root = Path(r2_matrix_dir)
    if not root.exists():
        raise LDSCInputError(
            f"R2 matrix directory does not exist: '{root}'. Most likely the "
            "path is misspelled or relative to a different working directory. Pass the "
            "existing `ldsc build-r2-matrix` output directory."
        )
    if not root.is_dir():
        raise LDSCInputError(
            f"R2 matrix path is not a directory: '{root}'. Most likely a file "
            "path was passed where an R2 matrix directory was expected. Pass the directory "
            "containing `chrN_r2.parquet` files."
        )
    if r2_matrix_chromosomes(root):
        return root

    requested_build = genome_build if genome_build in {"hg19", "hg38"} else None
    build_dirs = {build: root / build for build in ("hg19", "hg38") if (root / build).is_dir()}
    if requested_build is not None and requested_build in build_dirs:
        return build_dirs[requested_build]
    if len(build_dirs) > 1:
        builds = ", ".join(sorted(build_dirs))
        raise LDSCInputError(
            f"Reference-panel directory '{root}' is ambiguous because it contains multiple "
            f"genome-build directories ({builds}). Most likely a parent output directory "
            "was passed while genome_build is unset or auto. Pass the build-specific "
            "directory or set genome_build to hg19 or hg38."
        )
    if len(build_dirs) == 1:
        return next(iter(build_dirs.values()))
    return root


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



def _arrow_to_numpy(column):
    """Convert an Arrow column to a NumPy array, copying when Arrow cannot share memory."""
    return column.to_numpy(zero_copy_only=False)


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



def _check_r2_parquet_contract(path) -> tuple[R2ParquetContract | None, list[str]]:
    """Return the parsed contract of one R2 Parquet file, or ``None`` and its defects."""
    import pyarrow as pa
    import pyarrow.parquet as pq

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



def read_r2_sidecar(path: str | PathLike[str], *, rows: bool = True) -> tuple[dict[str, str], list[str], pd.DataFrame | None]:
    """Parse one sidecar without validating it.

    Returns the leading ``# ldsc:<key>=<value>`` identity lines as a dict, the
    header columns, and the rows. Rows are parsed only when ``rows`` is true and
    the header is :data:`SIDECAR_COLUMNS`; otherwise the frame is ``None``.
    CHR, SNP, A1 and A2 are read as text.
    """
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        identity: dict[str, str] = {}
        line = stream.readline()
        while line.startswith("#"):
            key, sep, value = line[1:].strip().partition("=")
            if sep and key.startswith("ldsc:"):
                identity[key.removeprefix("ldsc:")] = value
            line = stream.readline()
        header = line.rstrip("\r\n").split("\t")
        if not rows or tuple(header) != SIDECAR_COLUMNS:
            return identity, header, None
        try:
            frame = pd.read_csv(stream, sep="\t", header=None, names=list(SIDECAR_COLUMNS),
                                dtype={"CHR": str, "SNP": str, "A1": str, "A2": str})
        except pd.errors.EmptyDataError:
            frame = pd.DataFrame(columns=list(SIDECAR_COLUMNS))
    return identity, header, frame


class R2MatrixChromosome:
    """One chromosome of an R2 matrix: its R2 Parquet file and sidecar.

    ``ldscore`` and ``query-r2-matrix`` read an R2 matrix through this handle.
    Opening checks that both files exist and that the R2 Parquet footer matches
    what the current ``build-r2-matrix`` writes (:attr:`contract`). The sidecar
    is parsed and validated on the first :meth:`sidecar` call and kept, so a
    process parses it once per chromosome however many readers it builds.
    ``command`` names the reading command in every error and selects its
    troubleshooting section.

    Use :meth:`open` or :meth:`open_all`; the constructor does no checks.
    """

    def __init__(self, chrom: str, parquet_path: Path, sidecar_path: Path, contract: R2ParquetContract,
                 command: str) -> None:
        """Bind already checked paths and footer contract; see :meth:`open`."""
        self.chrom = chrom
        self.parquet_path = parquet_path
        self.sidecar_path = sidecar_path
        self.contract = contract
        self.command = command
        self._sidecar: pd.DataFrame | None = None

    @classmethod
    def open(cls, build_dir: str | PathLike[str], chrom: str, *, command: str) -> R2MatrixChromosome:
        """Open one chromosome of the R2 matrix in ``build_dir``; see :meth:`open_all`."""
        return cls.open_all(build_dir, [chrom], command=command)[0]

    @classmethod
    def open_all(cls, build_dir: str | PathLike[str], chromosomes: Sequence[str], *,
                 command: str) -> list[R2MatrixChromosome]:
        """Open several chromosomes, reporting every missing file and footer defect in one error.

        Parameters
        ----------
        build_dir : path-like
            R2 matrix build directory (see :func:`resolve_r2_build_dir`).
        chromosomes : sequence of str
            Chromosome labels to open, in the order to return them.
        command : {"ldscore", "query-r2-matrix"}
            Command reading the R2 matrix, named in errors.

        Raises
        ------
        LDSCInputError
            A chromosome's R2 Parquet file or sidecar is missing, or its footer
            does not match the current ``build-r2-matrix`` format.
        """
        opened, problems = [], []
        for chrom in chromosomes:
            chrom = normalize_chromosome(chrom)
            parquet_path, sidecar_path = r2_matrix_file_paths(build_dir, chrom)
            missing = [f"the {label} '{path}' is missing" for label, path in
                       (("R2 Parquet file", parquet_path), ("sidecar", sidecar_path)) if not path.is_file()]
            if missing:
                problems.append((chrom, "; ".join(missing)))
                continue
            contract, defects = _check_r2_parquet_contract(parquet_path)
            if defects:
                problems.append((chrom, f"R2 Parquet file '{parquet_path}' does not match the format written by "
                                        f"the current `ldsc build-r2-matrix`: {'; '.join(defects)}"))
                continue
            opened.append(cls(chrom, parquet_path, sidecar_path, contract, command))
        if problems:
            if len(problems) == 1:
                subject, details = f"R2 matrix chromosome {problems[0][0]}", problems[0][1]
            else:
                subject = f"{len(problems)} R2 matrix chromosomes"
                details = " | ".join(f"chromosome {chrom}: {text}" for chrom, text in problems)
            raise LDSCInputError(
                f"{command} could not use {subject}: {details}. Most likely the R2 matrix is incomplete, was built "
                "by an older LDSC3 version, or its files were altered or separated. Keep each chromosome's R2 "
                "Parquet file and sidecar together, and rebuild the R2 matrix with the current "
                f"`ldsc build-r2-matrix` if they are damaged. Other causes & fixes: {R2_MATRIX_TROUBLESHOOTING[command]}"
            )
        return opened

    def _error(self, reason: str, likely: str) -> LDSCInputError:
        """Return a run-aborting error about this chromosome, worded for the reading command."""
        return LDSCInputError(
            f"{self.command} could not use R2 matrix chromosome {self.chrom}: {reason}. Most likely {likely}. "
            "Rebuild the R2 matrix with the current `ldsc build-r2-matrix`, and keep each chromosome's R2 "
            f"Parquet file and sidecar together. Other causes & fixes: {R2_MATRIX_TROUBLESHOOTING[self.command]}"
        )

    def check_sidecar_header(self) -> None:
        """Check the sidecar's identity lines and header without parsing its rows."""
        identity, header, _ = read_r2_sidecar(self.sidecar_path, rows=False)
        self._check_preamble(identity, header)

    def _check_preamble(self, identity: dict[str, str], header: list[str]) -> None:
        """Require the writer's header and identity lines that match the R2 Parquet footer."""
        if tuple(header) != SIDECAR_COLUMNS:
            raise self._error(
                f"the sidecar '{self.sidecar_path}' has columns {', '.join(header)}, but the current "
                f"`ldsc build-r2-matrix` writes {', '.join(SIDECAR_COLUMNS)}",
                "the sidecar was edited, written by another tool, or written by an older LDSC3 version",
            )
        expected = {"artifact_type": "ref_panel_metadata", "snp_identifier": self.contract.snp_identifier,
                    "genome_build": self.contract.genome_build}
        found = {key: identity.get(key) for key in expected}
        if found != expected:
            def describe(values):
                return ", ".join(f"ldsc:{key}={value}" for key, value in values.items())

            raise self._error(
                f"the sidecar '{self.sidecar_path}' records {describe(found)}, but its R2 Parquet file "
                f"'{self.parquet_path}' records {describe(expected)}",
                "the two files come from different `build-r2-matrix` runs, or the sidecar was written by an "
                "older LDSC3 version or edited",
            )

    def sidecar(self) -> pd.DataFrame:
        """Return the validated sidecar, parsed on the first call and then kept.

        Returns
        -------
        pandas.DataFrame
            One row per R2 matrix row, in R2 matrix row order, with columns
            CHR (this chromosome), POS (int64), SNP, A1 and A2 (text), CM
            (float64, NaN where the build recorded no genetic map) and MAF
            (float64 in [0, 0.5]). Callers must not modify it.

        Raises
        ------
        LDSCInputError
            The sidecar's header or identity lines differ from what
            ``build-r2-matrix`` writes for this R2 Parquet file, its rows name
            another chromosome or hold invalid POS, CM or MAF values, or its
            row count or CHR:POS:A1:A2 hash differs from the footer's
            ``ldsc:n_snps`` and ``ldsc:sidecar_identity_sha256``.
        """
        if self._sidecar is None:
            identity, header, frame = read_r2_sidecar(self.sidecar_path)
            self._check_preamble(identity, header)
            frame = self._validated_rows(frame)
            self._check_binding(frame)
            self._sidecar = frame
        return self._sidecar

    def _validated_rows(self, frame: pd.DataFrame) -> pd.DataFrame:
        """Check the sidecar rows and return them with normalized dtypes."""
        def chromosome(label):
            try:
                return normalize_chromosome(label)
            except LDSCInputError:
                return str(label)

        labels = sorted({chromosome(label) for label in frame["CHR"].unique()})
        if labels and labels != [self.chrom]:
            raise self._error(
                f"the sidecar '{self.sidecar_path}' lists chromosome(s) {', '.join(labels)}, but it belongs to "
                f"chromosome {self.chrom}",
                "the sidecar was copied from another chromosome or R2 matrix; renaming a sidecar does not change "
                "its SNP identities",
            )
        positions = pd.to_numeric(frame["POS"], errors="coerce")
        cm = pd.to_numeric(frame["CM"], errors="coerce")
        maf = pd.to_numeric(frame["MAF"], errors="coerce").astype(float)
        bad = {
            "POS": positions.isna() | (positions != positions.round()),
            "CM": cm.isna() & frame["CM"].notna(),
            # NaN fails both comparisons, so this also catches missing and non-numeric values.
            "MAF": ~((maf >= 0.0) & (maf <= 0.5)),
        }
        problems = [f"{column} for {int(mask.sum())} of {len(frame)} rows (first SNPs: "
                    f"{', '.join(frame.loc[mask, 'SNP'].astype(str).head(5))})" for column, mask in bad.items() if mask.any()]
        if problems:
            raise self._error(
                f"the sidecar '{self.sidecar_path}' has missing, non-numeric, or out-of-range "
                + "; ".join(problems) + ". `ldsc build-r2-matrix` writes an integer POS, a numeric or NA CM, and "
                "the minor-allele (A1) frequency in [0, 0.5] as MAF for every row",
                "the sidecar was edited or produced by another tool",
            )
        return pd.DataFrame({
            "CHR": self.chrom,
            "POS": positions.astype(np.int64),
            "SNP": frame["SNP"].astype(str),
            "A1": frame["A1"].astype(str),
            "A2": frame["A2"].astype(str),
            "CM": cm.astype(float),
            "MAF": maf,
        }, index=pd.RangeIndex(len(frame)))

    def _check_binding(self, frame: pd.DataFrame) -> None:
        """Require the sidecar rows that the R2 Parquet indices refer to."""
        if len(frame) != self.contract.n_snps:
            raise self._error(
                f"the sidecar '{self.sidecar_path}' has {len(frame)} rows, but its R2 Parquet file "
                f"'{self.parquet_path}' records ldsc:n_snps={self.contract.n_snps}",
                "the R2 Parquet file and sidecar are not a matched pair",
            )
        actual = sidecar_identity_sha256(frame, context=f"R2 matrix sidecar {self.sidecar_path}")
        if actual != self.contract.sidecar_identity_sha256:
            raise self._error(
                f"the CHR:POS:A1:A2 hash {actual} of the sidecar '{self.sidecar_path}' differs from "
                f"ldsc:sidecar_identity_sha256 {self.contract.sidecar_identity_sha256} in its R2 Parquet file "
                f"'{self.parquet_path}'",
                "the sidecar was reordered or edited, or belongs to another R2 matrix",
            )

    def check_ld_window(self, args) -> None:
        """Reject ``ldscore`` LD windows this R2 matrix is not known to cover.

        The rule uses only the recorded footer LD-window mode and value, so it is
        decided before any computation: a window in another unit is rejected,
        and a same-unit window must be no wider than the R2 matrix's. Coverage
        is not inferred from per-SNP coordinates. For ``--ld-wind-snps`` this
        relies on counting R2 matrix rows (see :func:`window_coordinates`).
        """
        recorded = self.contract.ld_window
        requested = _requested_ld_window_spec(args)
        if requested.mode != recorded.mode:
            requested_text, recorded_text = _format_ld_window_spec(requested), _format_ld_window_spec(recorded)
            raise LDSCUsageError(
                f"ldscore cannot use the R2 matrix for chromosome {self.chrom} with "
                f"`{requested_text}`: the R2 matrix was built with `{recorded_text}`, and windows in different units "
                "cannot be checked for coverage from their settings, so the run could silently miss "
                f"SNP pairs the R2 matrix never stored. Rerun with the R2 matrix's window `{recorded_text}` (or a "
                f"narrower `{recorded.flag}` value), or rebuild the R2 matrix with "
                f"`ldsc build-r2-matrix {requested_text}`."
            )
        if requested.value > recorded.value:
            raise LDSCUsageError(
                f"ldscore cannot compute chromosome {self.chrom} from the "
                f"R2 matrix: the user-requested LD window `{_format_ld_window_spec(requested)}` "
                f"is wider than the input R2 matrix window "
                f"`{_format_ld_window_spec(recorded)}`. Most likely the "
                "R2 matrix was built with a smaller LD window, so it does not contain all "
                "SNP pairs required by this ldscore run. Rebuild the R2 matrix "
                f"with `{_format_ld_window_spec(requested)}` or a wider build window, or "
                f"rerun ldscore with an LD window no wider than the recorded R2 matrix window "
                f"`{_format_ld_window_spec(recorded)}`."
            )

    def reject_truncation(self) -> None:
        """Reject, for ``ldscore``, an R2 matrix that omitted pairs below a positive ``ldsc:r2_min``.

        LD scores sum every SNP pair in the window, so pairs the R2 matrix did not
        store would count as zero. Such R2 matrices stay valid for
        ``query-r2-matrix``, which reports unstored pairs as absent.
        """
        r2_min = self.contract.r2_min
        if r2_min > 0:
            raise LDSCInputError(
                f"ldscore cannot use R2 Parquet file '{self.parquet_path}' (chromosome {self.chrom}): the R2 matrix "
                f"omits SNP pairs with unbiased R2 below {r2_min:g} (ldsc:r2_min, set by `build-r2-matrix "
                "--r2-min`). LD scores sum every SNP pair in the window, so the omitted pairs would count as zero "
                "and change the scores. Rebuild the R2 matrix without `--r2-min` (the default keeps every pair). "
                "R2 matrices built with a positive `--r2-min` remain usable for `ldsc query-r2-matrix`."
            )

    def pair_reader(self, metadata: pd.DataFrame, identifier_mode: str,
                    genome_build: str | None = None) -> SortedR2BlockReader:
        """Return a reader that streams this chromosome's pairs among the SNPs in ``metadata``; the caller closes it."""
        return SortedR2BlockReader(self, metadata, identifier_mode, genome_build)


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
    and ``SIGN_R``. Values are already unbiased and are not corrected again.
    Each row group is decoded once and remapped to the retained SNPs. Build it
    with :meth:`R2MatrixChromosome.pair_reader`.
    """

    def __init__(
        self,
        chromosome: R2MatrixChromosome,
        metadata: pd.DataFrame,
        identifier_mode: str,
        genome_build: str | None = None,
    ) -> None:
        """Map the retained SNPs in ``metadata`` (canonical columns) to R2 matrix rows and open the Parquet file."""
        import pyarrow.parquet as pq

        self.chrom = chromosome.chrom
        self.identifier_mode = normalize_snp_identifier_mode(identifier_mode)
        assert genome_build in {"hg19", "hg38", None}, (
            f"genome_build={genome_build!r} must be concrete by this point."
        )
        self.contract = chromosome.contract
        validate_retained_identifier_uniqueness(metadata, self.identifier_mode, self.chrom)
        if genome_build not in {None, self.contract.sorted_by_build}:
            raise LDSCInputError(
                f"ldscore could not use R2 Parquet file '{chromosome.parquet_path}': it is sorted for "
                f"{self.contract.sorted_by_build}, but the analysis uses {genome_build}. Most likely "
                "the R2 matrix and annotation inputs use different genome builds. "
                "Use the matching reference file or regenerate with "
                f"`--genome-build {genome_build}`. "
                f"Other causes & fixes: {_LDSCORE_PARQUET_DOC}"
            )
        self.genome_build = self.contract.sorted_by_build
        self._remap, self._retained_build_idx = build_index_remap(
            chromosome.sidecar(), metadata, self.identifier_mode
        )
        self.m = len(metadata)
        self._pf = pq.ParquetFile(chromosome.parquet_path)

    def close(self):
        """Release the parquet handle owned by this chromosome reader."""
        self._pf.close()

    @property
    def panel_rows(self) -> np.ndarray:
        """R2 matrix sidecar row of each retained SNP, in retained order."""
        return self._retained_build_idx

    def _decode_index_batch(self, table, row_group_index: int) -> _DecodedR2RowGroup:
        """Dequantize and remap one bounded Arrow batch of stored pairs."""
        idx1 = _arrow_to_numpy(table.column("IDX_1")).astype(np.int64, copy=False)
        idx2 = _arrow_to_numpy(table.column("IDX_2")).astype(np.int64, copy=False)
        stored = _arrow_to_numpy(table.column("R2")).astype(np.float32, copy=False)
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
    import pyarrow as pa
    import pyarrow.parquet as pq

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
