"""Writers for the ``ldsc build-r2-matrix`` and ``ldsc query-r2-matrix`` directories.

``R2MatrixBuildDirectoryWriter`` declares the files an R2 matrix build owns,
scoped to one chromosome for a concrete single-chromosome run, and writes the
build metadata and the per-chromosome dropped-SNP tables; the R2 Parquet files
and sidecars are written by ``r2_format``. ``QueryR2MatrixDirectoryWriter``
declares the files a pair query owns and writes the queried pairs with their
metadata.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Sequence

import pandas as pd

from ..config import normalize_required_path
from .._core.snp_identity import coerce_identity_drop_frame
from .._runtime.artifacts import (
    ArtifactFamily,
    atomic_write_dataframe,
    atomic_write_json,
    declare_artifacts,
    result_metadata,
)
from .._runtime.paths import ensure_output_directory
from .._runtime.progress import report_phase
from .._runtime.workflow_logging import FAILURE_MARKER
from .r2_format import R2_MATRIX_FILE_GLOBS, r2_matrix_file_paths

LOGGER = logging.getLogger("LDSC.r2_matrix_builder")
BUILD_R2_MATRIX_LOG = "diagnostics/build-r2-matrix.log"


class R2MatrixBuildDirectoryWriter:
    """Declare and write the ``build-r2-matrix`` output directory.

    The directory holds ``{build}/chr{chrom}_r2.parquet`` and
    ``{build}/chr{chrom}_meta.tsv.gz`` for every emitted build,
    ``diagnostics/dropped_snps/chr{chrom}_dropped.tsv.gz`` for every processed
    chromosome, the build metadata, and the workflow log. A concrete
    single-chromosome run (``scope`` is that chromosome) owns only its own
    chromosome's files and names its metadata, log and failure marker after
    it, so parallel chromosome runs can share one directory; a run over an
    ``@`` chromosome suite (``scope`` is ``None``) owns the whole directory.
    """

    @staticmethod
    def chromosome_scope(plink_prefix, chromosomes: Sequence[str]) -> str | None:
        """Return the one chromosome a concrete single-chromosome run owns, or ``None``."""
        return str(chromosomes[0]) if len(chromosomes) == 1 and "@" not in str(plink_prefix) else None

    @staticmethod
    def metadata_path(output_dir, scope: str | None) -> Path:
        """Return the build metadata path, chromosome-scoped for a single-chromosome run."""
        name = "metadata.json" if scope is None else f"metadata.chr{scope}.json"
        return Path(output_dir) / "diagnostics" / name

    @staticmethod
    def log_path(log_path, scope: str | None) -> Path:
        """Return the workflow log path, chromosome-scoped for a single-chromosome run."""
        log_path = Path(log_path)
        return log_path if scope is None else log_path.with_name(f"{log_path.stem}.chr{scope}{log_path.suffix}")

    @staticmethod
    def failure_marker(scope: str | None) -> str:
        """Return the failed-overwrite marker name, scoped like the metadata."""
        return FAILURE_MARKER if scope is None else f"RUN_FAILED.chr{scope}.txt"

    @staticmethod
    def dropped_snps_path(output_dir, chrom: str) -> Path:
        """Return one chromosome's dropped-SNP table path."""
        return Path(output_dir) / "diagnostics" / "dropped_snps" / f"chr{chrom}_dropped.tsv.gz"

    @classmethod
    def artifact_family(
        cls, output_dir, *, scope: str | None, builds: Sequence[str] | None = None, chromosomes: Sequence[str] = (),
    ) -> ArtifactFamily:
        """Describe the files this run owns and, once the builds are known, the files it writes.

        With ``builds=None`` (before the source build is resolved) nothing is
        selected and only ownership is described. Otherwise the selection is
        the metadata, every emitted build's R2 Parquet file and sidecar, and
        every chromosome's dropped-SNP table, which is written even when
        header-only. Ownership adds the existing files of the same kinds in
        the run's scope, including stale chromosomes and logs; the run's own
        log is passed to ``ArtifactFamily.preflight(additional_paths=...)``.
        """
        output_dir = Path(output_dir)
        paths: dict[str, Path] = {}
        if builds is not None:
            paths["metadata"] = cls.metadata_path(output_dir, scope)
            for chrom in chromosomes:
                for build in builds:
                    r2_path, meta_path = r2_matrix_file_paths(output_dir / build, chrom)
                    paths[f"r2_{build}_chr{chrom}"] = r2_path
                    paths[f"meta_{build}_chr{chrom}"] = meta_path
                paths[f"dropped_snps_chr{chrom}"] = cls.dropped_snps_path(output_dir, chrom)
        return ArtifactFamily(output_dir, paths, (*cls._existing_owned_paths(output_dir, scope), *paths.values()),
                              "R2 matrix output artifact")

    @classmethod
    def _existing_owned_paths(cls, output_dir: Path, scope: str | None) -> list[Path]:
        """Return the owned paths of every chromosome in scope, found on disk for a whole-directory run."""
        diagnostics_dir = output_dir / "diagnostics"
        if scope is None:
            paths: list[Path] = [diagnostics_dir / "metadata.json"]
            paths.extend(sorted(diagnostics_dir.glob("metadata.chr*.json")))
        else:
            paths = [cls.metadata_path(output_dir, scope)]
        for build_subdir_name in ("hg19", "hg38"):
            build_subdir = output_dir / build_subdir_name
            if not build_subdir.is_dir():
                continue
            if scope is None:
                for pattern in R2_MATRIX_FILE_GLOBS:
                    paths.extend(sorted(build_subdir.glob(pattern)))
            else:
                paths.extend(r2_matrix_file_paths(build_subdir, scope))
        dropped_dir = diagnostics_dir / "dropped_snps"
        if dropped_dir.is_dir():
            if scope is None:
                paths.extend(sorted(dropped_dir.glob("chr*_dropped.tsv.gz")))
            else:
                paths.append(cls.dropped_snps_path(output_dir, scope))
        if scope is None:
            paths.extend(sorted(diagnostics_dir.glob("build-r2-matrix*.log")))
        else:
            paths.append(cls.log_path(output_dir / BUILD_R2_MATRIX_LOG, scope))
        return paths

    @staticmethod
    @report_phase(LOGGER, 'publication', 'R2 matrix metadata')
    def write_metadata(path: Path, output_dir, output_paths: dict[str, list[str]], provenance: dict) -> None:
        """Write the build metadata: ``artifact_type``, output-relative ``files``, then ``provenance``."""
        output_dir = Path(output_dir)
        files = {
            key: [_relative_to_output_dir(item, output_dir) for item in values]
            for key, values in output_paths.items()
            if key != "metadata"
        }
        atomic_write_json({"artifact_type": "r2_matrix_build", "files": files, **provenance}, path)

    @staticmethod
    def write_dropped_snps(dropped_df: pd.DataFrame, path: Path, chrom: str) -> None:
        """Write one chromosome's dropped-SNP table, header-only when nothing was dropped.

        Rows use the shared identity-drop columns; rows without a ``stage``
        come from the liftover stage (``ref_panel_liftover``).
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        dropped_df = coerce_identity_drop_frame(dropped_df, default_stage="ref_panel_liftover")
        dropped_df.to_csv(path, sep="\t", index=False, compression="gzip", na_rep="")
        n_dropped = len(dropped_df)
        if n_dropped == 0:
            LOGGER.info(f"No SNPs dropped on chromosome {chrom}; audit sidecar at '{path}'.")
            return
        counts = dropped_df["reason"].value_counts(sort=False)
        count_text = ", ".join(f"{reason}={int(count)}" for reason, count in counts.items())
        LOGGER.info(
            f"R2 matrix dropped SNPs: {n_dropped} SNPs on chromosome {chrom} "
            f"({count_text}); audit sidecar at '{path}'."
        )


def _relative_to_output_dir(path: str, output_dir: Path) -> str:
    """Return output-relative path strings when possible."""
    candidate = Path(path)
    try:
        return str(candidate.relative_to(output_dir))
    except ValueError:
        return str(candidate)


@dataclass(frozen=True)
class QueryR2MatrixOutputConfig:
    """Directory-oriented output config for query-r2-matrix pair results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``query_r2_matrix.tsv`` and ``diagnostics/metadata.json``.
    overwrite : bool, optional
        If ``True``, replace existing query-r2-matrix outputs. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))


class QueryR2MatrixDirectoryWriter:
    """Write query-r2-matrix results and diagnostic provenance.

    The workflow supplies input columns plus ``r2``, nullable ``sign_r``,
    ``r``, and ``status``. The TSV preserves those columns and writes missing
    values as ``NaN``; it does not add compatibility aliases.
    """

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare the pair-query result table and provenance."""
        return declare_artifacts(output_dir, {
            "result": ("query_r2_matrix.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
            "input_issues": ("diagnostics/input_issues.tsv", False),
        }, label="query-r2-matrix output artifact")

    def write(
        self,
        result: pd.DataFrame,
        output_config: QueryR2MatrixOutputConfig,
        *,
        metadata: dict[str, object],
    ) -> dict[str, str]:
        """Write ``query_r2_matrix.tsv`` and ``diagnostics/metadata.json`` to a result directory.

        Existing fixed query-r2-matrix artifacts are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        result_path = family.paths["result"]
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_dataframe(result, result_path, na_rep="NaN")
        atomic_write_json(
            result_metadata(metadata, artifact_type="query_r2_matrix_result", files=family.metadata_files()),
            metadata_path,
        )
        return {name: str(path) for name, path in family.paths.items()}
