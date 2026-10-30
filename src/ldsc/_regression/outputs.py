"""Output writers for the LDSC3 workflows.

The public LD-score workflow writes one LD-score directory containing
``metadata.json``, ``ldscore.baseline.parquet``, optional query files
(``ldscore.query.parquet`` for one batch or numbered
``ldscore.query.batchNNNNN.parquet`` files for multiple batches), and an optional
``ldscore.overlap.parquet`` (the
annotation overlap matrix required by partitioned-h2's overlap-aware annotation
tables and read by the shared h2 collinearity guard when there are two or more
LD-score columns; ``rg`` does not use it). The overlap sidecar is
written only when the run has two or more annotation columns; an unpartitioned
single-annotation run (e.g. the synthetic ``base``) omits it because the matrix
would collapse to a SNP count already in ``metadata.json``. Query runs also
write ``diagnostics/query_annotation_status.tsv`` and gene-list runs write a
row-complete audit plus per-source resolution summary. Run identity comes from the chosen
directory name;
output filenames inside that directory are fixed. The parquet payloads are
written with one row group per chromosome and matching metadata so
downstream readers can load a single chromosome without scanning the whole
table. The writer creates missing directories, reuses existing directories, and
refuses existing canonical family files unless the caller explicitly sets
``overwrite=True``; successful overwrites remove stale owned siblings that the
current result did not produce. Root ``metadata.json`` is part of the
downstream contract for LD-score directories. Regression
metadata emitted by this module is diagnostic provenance and is written below
``diagnostics/`` without legacy top-level ``format`` discriminators.

Each directory writer's ``artifact_family`` declares its owned paths and
conditional outputs. Workflows consult that declaration for early collision
checks, adding their own log/audit paths. Writers use the final result to
select produced paths, metadata file entries, and post-write stale cleanup;
workflows do not retain a predicted stale list across computation.

The ``partitioned-h2`` regression summaries use the same directory-oriented output
policy. ``PartitionedH2DirectoryWriter`` writes the aggregate ``partitioned_h2.tsv``
(the scientific ``PARTITIONED_H2_COLUMNS`` plus ``query_result_dir`` for both
regimes) and can optionally stage a
per-query ``diagnostics/query_annotations`` tree with ``manifest.tsv``, one-row
``partitioned_h2.tsv`` summaries, full ``partitioned_h2_full.tsv`` model
tables, and ``metadata.json`` files before moving it into place. ``RgDirectoryWriter``
uses the same whole-tree staging and replacement policy for optional
``diagnostics/pairs`` details.

Unpartitioned h2 results also carry
``diagnostics/ld_score_regression_bins.tsv``. This table summarizes the exact
SNP population and final fitted model used by h2 and is the authoritative input
for the post-processing binned LD score regression diagnostic plot. The writer
layer persists the table but never imports a plotting library.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import tempfile
import unicodedata
from dataclasses import dataclass
from os import PathLike
from pathlib import Path
from typing import Any

import pandas as pd

from ..config import normalize_required_path
from ..errors import LDSCInternalError
from .._runtime.artifacts import (
    ArtifactFamily,
    atomic_write_dataframe,
    atomic_write_json,
    declare_artifacts,
    result_metadata,
)
from .._runtime.paths import (
    ensure_output_directory,
    remove_output_artifacts,
)


# Fitted-model files of a partitioned-h2 result that other modules locate.
PARTITIONED_H2_SUMMARY_FILE = "partitioned_h2.tsv"
PARTITIONED_H2_FULL_FILE = "partitioned_h2_full.tsv"
COEFFICIENT_DELETE_VALUES_FILE = "coefficient_delete_values.parquet"
RESULT_METADATA_FILE = "diagnostics/metadata.json"
QUERY_RESULT_ROOT = "diagnostics/query_annotations"
QUERY_METADATA_FILE = "metadata.json"
# A query folder's link to the run metadata, relative to that folder: the
# folder sits two levels below RESULT_METADATA_FILE's directory.
QUERY_GLOBAL_METADATA = "../../metadata.json"

H2_REGRESSION_BIN_COLUMNS = [
    "bin",
    "n_snps",
    "ld_score_min",
    "ld_score_max",
    "mean_ld_score",
    "mean_chi_square",
    "sd_chi_square",
    "mean_sample_size",
    "mean_fitted_chi_square",
    "mean_regression_weight",
]
# Single partitioned-h2 schema for both regimes. The functional regime fills it
# with one row per baseline category; the cell-type regime fills it with one row
# per query annotation. Differs only in rows and default sort, never in columns.
PARTITIONED_H2_COLUMNS = [
    "category",
    "prop_snps",
    "prop_h2",
    "prop_h2_se",
    "enrichment",
    "enrichment_se",
    "enrichment_p",
    "coefficient",
    "coefficient_se",
    "coefficient_z",
    "coefficient_p",
    "overlap_annot",
    "total_h2_obs",
    "total_h2_obs_se",
    "total_h2_liab",
    "total_h2_liab_se",
    "category_h2_obs",
    "category_h2_obs_se",
    "category_h2_liab",
    "category_h2_liab_se",
    "samp_prev",
    "pop_prev",
]
QUANTILE_H2_COLUMNS = [
    "quantile", "target_value_lower", "target_value_upper", "n_snps", "prop_snps",
    "h2_obs", "h2_obs_se", "h2_liab", "h2_liab_se", "prop_h2", "prop_h2_se",
    "enrichment", "enrichment_se", "enrichment_p",
]
STANDARDIZED_COEFFICIENT_COLUMNS = [
    "annotation", "annotation_type", "annotation_sd", "tau", "tau_se", "tau_z", "tau_p",
    "tau_star", "tau_star_se", "tau_star_z", "tau_star_p",
]
SNP_ALIGNMENT_ISSUE_COLUMNS = [
    "source_role", "source", "CHR", "POS", "SNP", "A1", "A2", "effective_snp_id",
    "issue", "action", "details",
]
RG_CONCISE_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "p",
    "note",
]
RG_FULL_COLUMNS = [
    "trait_1",
    "trait_2",
    "n_snps_used",
    "rg",
    "rg_se",
    "z",
    "p",
    "h2_1_obs",
    "h2_1_obs_se",
    "h2_1_liab",
    "h2_1_liab_se",
    "h2_2_obs",
    "h2_2_obs_se",
    "h2_2_liab",
    "h2_2_liab_se",
    "gencov_obs",
    "gencov_obs_se",
    "gencov_liab",
    "gencov_liab_se",
    "samp_prev_1",
    "pop_prev_1",
    "samp_prev_2",
    "pop_prev_2",
    "intercept_h2_1",
    "intercept_h2_1_se",
    "intercept_h2_2",
    "intercept_h2_2_se",
    "intercept_gencov",
    "intercept_gencov_se",
    "ratio_1",
    "ratio_1_se",
    "ratio_2",
    "ratio_2_se",
    "lambda_gc_1",
    "lambda_gc_2",
    "mean_chisq_1",
    "mean_chisq_2",
    "pair_kind",
    "status",
    "error",
]


@dataclass(frozen=True)
class H2OutputConfig:
    """Directory-oriented output config for unpartitioned h2 summaries.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``h2.tsv``, the exact LD score regression-bin
        diagnostic, and metadata.
    overwrite : bool, optional
        If ``True``, replace existing fixed h2 outputs. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))


class H2DirectoryWriter:
    """Write the unpartitioned h2 summary, exact fitted-bin diagnostic, and metadata."""

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare h2 outputs and derived roots invalidated by a successful rewrite."""
        return declare_artifacts(output_dir, {
            "summary": ("h2.tsv", True),
            "ld_score_regression_bins": ("diagnostics/ld_score_regression_bins.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
            "plots": ("plots", False),
            "postprocessing": ("postprocessing", False),
        }, label="h2 output artifact")

    def write(
        self,
        summary: pd.DataFrame,
        output_config: H2OutputConfig,
        *,
        metadata: dict[str, object],
        diagnostic_bins: pd.DataFrame,
    ) -> dict[str, str]:
        """Write the h2 summary, regression bins, and metadata.

        Parameters
        ----------
        summary : pandas.DataFrame
            One-row h2 summary in the canonical public schema.
        output_config : H2OutputConfig
            Output directory and overwrite policy.
        metadata : dict
            Scientific and source provenance for the result.
        diagnostic_bins : pandas.DataFrame
            Exact post-filter LD score regression-bin summary in
            :data:`H2_REGRESSION_BIN_COLUMNS` order.

        Returns
        -------
        dict of str to str
            Written summary, bin-diagnostic, and metadata paths.

        Raises
        ------
        ValueError
            If the summary or diagnostic table lacks required columns.
        FileExistsError
            If an owned h2 artifact exists and overwrite is disabled.

        Notes
        -----
        Existing fixed h2 artifacts are checked before any output file is
        written. Replacement requires ``output_config.overwrite=True``. A
        successful replacement also removes default plot and post-processing
        roots derived from the superseded h2 result.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        paths = family.paths
        summary_path = paths["summary"]
        metadata_path = paths["metadata"]
        bins_path = paths["ld_score_regression_bins"]
        diagnostics_dir = metadata_path.parent
        stale = family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_dataframe(summary, summary_path, na_rep="NaN")
        atomic_write_dataframe(
            _select_columns(diagnostic_bins, H2_REGRESSION_BIN_COLUMNS, label="h2 regression-bin diagnostic"),
            bins_path,
            na_rep="NaN",
        )
        atomic_write_json(
            result_metadata(
                metadata,
                artifact_type="h2_result",
                files=family.metadata_files(),
            ),
            metadata_path,
        )
        remove_output_artifacts(stale)
        return {name: str(path) for name, path in family.paths.items()}


@dataclass(frozen=True)
class QuantileH2OutputConfig:
    """Output directory and replacement policy for ``quantile-h2`` results.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory receiving the two scientific TSVs and ``diagnostics/``.
    overwrite : bool, optional
        Replace existing owned quantile-h2 artifacts when ``True``. Default is
        ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False

    def __post_init__(self) -> None:
        """Normalize the required output directory."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))


class QuantileH2DirectoryWriter:
    """Write post-fit quantile summaries, coefficients, and SNP diagnostics.

    A successful overwrite removes the default plot root derived from the
    superseded quantile result.
    """

    @staticmethod
    def artifact_family(output_dir) -> ArtifactFamily:
        """Declare quantile outputs, alignment diagnostics, and the derived plot root."""
        return declare_artifacts(output_dir, {
            "quantile_h2": ("quantile_h2.tsv", True),
            "standardized_coefficients": ("standardized_coefficients.tsv", True),
            "metadata": ("diagnostics/metadata.json", True),
            "snp_alignment_issues": ("diagnostics/snp_alignment_issues.tsv.gz", True),
            "plots": ("plots", False),
        }, label="quantile-h2 output artifact")

    def write(
        self,
        quantile_h2: pd.DataFrame,
        standardized_coefficients: pd.DataFrame,
        alignment_issues_path: str | Path,
        output_config: QuantileH2OutputConfig,
        metadata: dict[str, object],
    ) -> dict[str, str]:
        """Write the complete successful ``quantile-h2`` artifact family.

        Parameters
        ----------
        quantile_h2 : pandas.DataFrame
            Low-to-high quantile rows in the stable public schema.
        standardized_coefficients : pandas.DataFrame
            One row per fitted annotation in coefficient order.
        alignment_issues_path : str or pathlib.Path
            Complete streamed gzip diagnostic table, including its header.
        output_config : QuantileH2OutputConfig
            Output directory and replacement policy.
        metadata : dict
            Scientific and source provenance added to the standard result
            metadata envelope.

        Returns
        -------
        dict of str to str
            Paths for both scientific tables, metadata, and SNP diagnostics.
        """
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir)
        paths = family.paths
        diagnostics_dir = paths["metadata"].parent
        stale = family.preflight(overwrite=output_config.overwrite)
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        atomic_write_dataframe(
            _select_columns(quantile_h2, QUANTILE_H2_COLUMNS, label="quantile-h2 summary"),
            paths["quantile_h2"],
            na_rep="NaN",
        )
        atomic_write_dataframe(
            _select_columns(
                standardized_coefficients,
                STANDARDIZED_COEFFICIENT_COLUMNS,
                label="standardized coefficients",
            ),
            paths["standardized_coefficients"],
            na_rep="NaN",
        )
        if Path(alignment_issues_path).resolve() != paths["snp_alignment_issues"].resolve():
            shutil.copyfile(alignment_issues_path, paths["snp_alignment_issues"])
        payload = result_metadata(
            metadata,
            artifact_type="quantile_h2_result",
            files=family.metadata_files(),
        )
        atomic_write_json(payload, paths["metadata"])
        remove_output_artifacts(stale)
        return {name: str(path) for name, path in paths.items()}


@dataclass(frozen=True)
class PartitionedH2OutputConfig:
    """Directory-oriented output config for partitioned-h2 regression summaries.

    The config keeps the partitioned-h2 output contract explicit and separate
    from LD-score parquet output. The default mode writes only the stable
    aggregate summary. ``write_per_query_results=True`` enables the richer
    per-query tree without changing the aggregate file name.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``partitioned_h2.tsv`` and, when requested,
        the optional ``diagnostics/query_annotations`` per-query result tree.
    overwrite : bool, optional
        If ``True``, replace existing fixed partitioned-h2 outputs and remove
        stale owned siblings after a successful write. Default is ``False``.
    write_per_query_results : bool, optional
        If ``True``, also write one subdirectory per query annotation under
        ``diagnostics/query_annotations``. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    write_per_query_results: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))


@dataclass(frozen=True)
class PartitionedH2FitArtifacts:
    """Paths to one complete fitted model and its publication ownership.

    ``owned=False`` preserves persistent files by copying them at publication.
    ``owned=True`` transfers the full table and delete values from disposable
    staging on the output filesystem; those source paths are consumed even if
    later publication fails. The caller remains responsible for scratch cleanup.
    """

    full: Path
    delete_values: Path
    metadata: Path
    owned: bool = False


def stage_partitioned_h2_fit(directory, category_table, delete_values, metadata):
    """Write one fit privately before releasing its estimator and matrices.

    The caller supplies a unique directory inside its owned output workspace.
    Returned artifacts have ``owned=True`` and remain private until the final
    summary determines folder ordering. The directory writer moves their full
    table and delete values into the complete query tree without rewriting them.
    Concurrent query workers receive distinct ordinal directories from the
    coordinator. Only the coordinator assigns public folder names and publishes
    results, after applying the run's query-error policy.
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    artifacts = _query_fit_files(directory, owned=True)
    atomic_write_dataframe(_select_columns(category_table, PARTITIONED_H2_COLUMNS, label="full partitioned-h2 summary"),
                            artifacts.full, na_rep="NaN")
    delete_values.to_parquet(artifacts.delete_values, index=False)
    atomic_write_json({**metadata, "coefficient_delete_block_count": int(len(delete_values)),
                       "retained_ld_columns": [c for c in delete_values if c != "delete_block"]}, artifacts.metadata)
    return artifacts


def _query_records(summary: pd.DataFrame) -> list[dict[str, object]]:
    """Return ordered query records with deterministic safe folder names."""
    records: list[dict[str, object]] = []
    width = max(4, len(str(len(summary))))
    for idx, query_name in enumerate(summary["category"].astype(str).tolist(), start=1):
        slug = _slugify_query_name(query_name)
        records.append({"ordinal": idx, "query_annotation": query_name, "slug": slug, "folder": f"{idx:0{width}d}_{slug}"})
    return records


def _query_fit_files(folder: Path, *, owned: bool = False) -> PartitionedH2FitArtifacts:
    """Name one fitted model's files inside a staged or published query folder."""
    return PartitionedH2FitArtifacts(folder / PARTITIONED_H2_FULL_FILE, folder / COEFFICIENT_DELETE_VALUES_FILE,
                                     folder / QUERY_METADATA_FILE, owned=owned)


def fitted_model_files(result_dir) -> tuple[str, PartitionedH2FitArtifacts]:
    """Locate one fitted model inside a partitioned-h2 result.

    Parameters
    ----------
    result_dir : str or os.PathLike[str]
        A baseline-only result root or one query folder under
        ``diagnostics/query_annotations``. A directory holding its own
        ``metadata.json`` is read as a query folder.

    Returns
    -------
    tuple of (str, PartitionedH2FitArtifacts)
        The ``artifact_type`` this writer records for that layout, and the
        model's full table, coefficient delete values and metadata. A
        baseline-only root's ``partitioned_h2.tsv`` is its full table. Files
        are named, not checked for existence.
    """
    root = Path(result_dir)
    if (root / QUERY_METADATA_FILE).is_file():
        return "partitioned_h2_query_result", _query_fit_files(root)
    return "partitioned_h2_result", PartitionedH2FitArtifacts(
        root / PARTITIONED_H2_SUMMARY_FILE, root / "diagnostics" / COEFFICIENT_DELETE_VALUES_FILE,
        root / RESULT_METADATA_FILE,
    )


def published_query_artifacts(output_dir, summary: pd.DataFrame) -> dict[str, PartitionedH2FitArtifacts]:
    """Return each query's published fitted-model files, keyed by query name.

    ``summary`` must be the table passed to ``PartitionedH2DirectoryWriter.write``,
    whose row order assigned the query folders. The returned descriptors have
    ``owned=False``, so republishing them copies the files.
    """
    query_root = Path(output_dir) / QUERY_RESULT_ROOT
    return {str(record["query_annotation"]): _query_fit_files(query_root / str(record["folder"]))
            for record in _query_records(summary)}


class PartitionedH2DirectoryWriter:
    """Write aggregate and optional per-query partitioned-h2 outputs.

    The writer owns the fixed regression summary layout:
    ``partitioned_h2.tsv`` at the output root uses the compact scientific schema
    plus ``query_result_dir``, relative to the output root and missing when no
    per-query tree is written. Optional diagnostic query folders contain the
    scientific columns as a one-row summary
    plus a full baseline-plus-query ``partitioned_h2_full.tsv`` table. The
    per-query tree
    is written to a temporary sibling directory before it is moved into the
    final location, so ordinary validation and I/O failures do not expose a
    partially populated final tree. A successful overwrite also removes the
    default plot root derived from the superseded result.
    An optional ``diagnostics/query_status.tsv`` records all attempted queries,
    including failed queries absent from the summary and model tree. The
    regression runner owns the strict/continue decision; this writer publishes
    only the successful model artifacts supplied by that runner.
    """

    @staticmethod
    def artifact_family(output_dir, *, write_per_query_results=False, coefficient_delete_values=False,
                        query_status=False) -> ArtifactFamily:
        """Declare summaries, model artifacts, and an optional attempted-query ledger."""
        return declare_artifacts(output_dir, {
            "summary": (PARTITIONED_H2_SUMMARY_FILE, True),
            "metadata": (RESULT_METADATA_FILE, True),
            "coefficient_delete_values": (f"diagnostics/{COEFFICIENT_DELETE_VALUES_FILE}", coefficient_delete_values),
            "query_annotations": (QUERY_RESULT_ROOT, write_per_query_results),
            "query_status": ("diagnostics/query_status.tsv", query_status),
            "plots": ("plots", False),
        }, label="partitioned-h2 output artifact")

    @staticmethod
    def write_query_status(output_dir, query_status):
        """Write attempted-query statuses after the caller's output preflight.

        This diagnostic can survive a failed run without publishing any
        scientific results. Successful runs write it with their result family.

        Parameters
        ----------
        output_dir : path-like
            Caller-preflighted result directory. This method does not request
            overwrite authorization or write scientific result files.
        query_status : pandas.DataFrame
            Attempted queries in input order, with ``query_annotation``,
            ``status``, ``stage``, ``error_type``, and ``error_message`` columns.

        Returns
        -------
        pathlib.Path
            Path to the written ``diagnostics/query_status.tsv``.
        """
        path = PartitionedH2DirectoryWriter.artifact_family(output_dir, query_status=True).paths["query_status"]
        path.parent.mkdir(parents=True, exist_ok=True)
        atomic_write_dataframe(query_status, path, na_rep="NaN")
        return path

    def write(
        self,
        summary: pd.DataFrame,
        output_config: PartitionedH2OutputConfig,
        per_query_artifacts: dict[str, PartitionedH2FitArtifacts] | None = None,
        metadata: dict[str, object] | None = None,
        coefficient_delete_values: pd.DataFrame | None = None,
        query_status: pd.DataFrame | None = None,
    ) -> dict[str, str]:
        """Write partitioned-h2 summary artifacts.

        Parameters
        ----------
        summary : pandas.DataFrame
            Aggregate partitioned-h2 table with the compact public columns in
            ``PARTITIONED_H2_COLUMNS``. The input is not modified. The saved
            root TSV appends ``query_result_dir`` relative to the output root;
            it is written as ``NaN`` when per-query output is disabled.
            Per-query summaries and full model tables keep scientific columns
            only. Folder ordinals follow the supplied summary row order.
        output_config : PartitionedH2OutputConfig
            Output directory, overwrite policy, and per-query mode.
        per_query_artifacts : dict of str to PartitionedH2FitArtifacts, optional
            Completed private or persistent fits keyed by original query name.
            Owned files are moved into the sorted publication tree; persistent
            files are copied. Neither loads annotation tables or delete-value
            matrices. Owned staging must share the output filesystem and may be
            consumed on failure. Workers may stage fits
            concurrently; the coordinator calls this writer after fitting and
            decides which successful artifacts may be published.
        metadata : dict, optional
            Run-level metadata stored once in ``diagnostics/metadata.json``.
            Each query's metadata references it via ``global_metadata`` relative
            to the query directory; model-specific fields take precedence.
        coefficient_delete_values : pandas.DataFrame, optional
            Baseline-only fitted model's delete-one-block coefficient vectors.
            The first column is ``delete_block`` and remaining columns retain
            fitted annotation order.
        query_status : pandas.DataFrame, optional
            One diagnostic row per attempted query, including failures absent
            from the scientific summary and per-query artifact tree.

        Returns
        -------
        dict of str to str
            Written path map. Always includes ``"summary"`` and additionally
            includes ``"per_query_root"`` and ``"per_query_manifest"`` when
            per-query output is enabled, and ``"query_status"`` when the
            attempted-query table is supplied.

        Raises
        ------
        ValueError
            If ``summary`` lacks the required ``Category`` column.
        FileExistsError
            If a final output path already exists and overwrite is disabled.
        """
        self._validate_summary(summary)
        root_summary = _select_columns(summary, PARTITIONED_H2_COLUMNS, label="partitioned-h2 summary").copy()
        root_summary["query_result_dir"] = pd.NA
        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(
            output_dir, write_per_query_results=output_config.write_per_query_results,
            coefficient_delete_values=coefficient_delete_values is not None,
            query_status=query_status is not None,
        )
        summary_path = family.paths["summary"]
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        query_root = family.paths.get("query_annotations")
        coefficient_delete_path = family.paths.get("coefficient_delete_values")
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        root_files = family.metadata_files()
        paths = {name: str(path) for name, path in family.paths.items() if name != "query_annotations"}
        metadata_payload = dict(metadata or {})
        if coefficient_delete_values is not None:
            metadata_payload["coefficient_delete_block_count"] = int(len(coefficient_delete_values))
            metadata_payload["retained_ld_columns"] = [
                column for column in coefficient_delete_values.columns if column != "delete_block"
            ]
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        if query_status is not None:
            self.write_query_status(output_dir, query_status)
        if not output_config.write_per_query_results:
            atomic_write_dataframe(
                root_summary,
                summary_path,
                na_rep="NaN",
            )
            if coefficient_delete_values is not None:
                coefficient_delete_values.to_parquet(coefficient_delete_path, index=False)
            atomic_write_json(
                result_metadata(metadata_payload, artifact_type="partitioned_h2_result", files=root_files),
                metadata_path,
            )
            remove_output_artifacts(stale_paths)
            return paths

        query_records = _query_records(summary)
        root_summary["query_result_dir"] = [f"{QUERY_RESULT_ROOT}/{record['folder']}" for record in query_records]
        staging_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.tmp.", dir=str(diagnostics_dir)))
        backup_dir: Path | None = None
        try:
            manifest_rows = self._write_staged_query_tree(
                staging_dir,
                query_records,
                summary,
                per_query_artifacts or {},
            )
            atomic_write_dataframe(pd.DataFrame(manifest_rows), staging_dir / "manifest.tsv", na_rep="NaN")
            atomic_write_dataframe(
                root_summary,
                summary_path,
                na_rep="NaN",
            )
            atomic_write_json(
                result_metadata(metadata_payload, artifact_type="partitioned_h2_result", files=root_files),
                metadata_path,
            )
            if output_config.overwrite and query_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".query_annotations.backup.", dir=str(diagnostics_dir)))
                backup_dir.rmdir()
                os.replace(query_root, backup_dir)
            os.replace(staging_dir, query_root)
        except Exception:
            if staging_dir.exists():
                shutil.rmtree(staging_dir, ignore_errors=True)
            if backup_dir is not None and backup_dir.exists() and not query_root.exists():
                os.replace(backup_dir, query_root)
            raise
        else:
            if backup_dir is not None and backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)

        remove_output_artifacts(stale_paths)
        paths["per_query_root"] = str(query_root)
        paths["per_query_manifest"] = str(query_root / "manifest.tsv")
        return paths

    def _validate_summary(self, summary: pd.DataFrame) -> None:
        """Validate the aggregate summary before writing any final files."""
        if "category" not in summary.columns:
            raise LDSCInternalError(
                "partitioned-h2 output writer could not validate the aggregate summary because the required "
                "'category' column is missing. Most likely the regression workflow returned a result table with "
                "the wrong schema. Re-run with DEBUG logging and report the traceback."
            )

    def _write_staged_query_tree(
        self,
        staging_dir: Path,
        query_records: list[dict[str, object]],
        summary: pd.DataFrame,
        per_query_artifacts: dict[str, PartitionedH2FitArtifacts],
    ) -> list[dict[str, object]]:
        """Populate the staged per-query result tree and return manifest rows."""
        manifest_rows: list[dict[str, object]] = []
        for record in query_records:
            query_name = str(record["query_annotation"])
            folder = str(record["folder"])
            query_dir = staging_dir / folder
            query_dir.mkdir(parents=True, exist_ok=False)
            query_summary = summary.loc[summary["category"].astype(str) == query_name].reset_index(drop=True)
            summary_rel = f"{QUERY_RESULT_ROOT}/{folder}/{PARTITIONED_H2_SUMMARY_FILE}"
            full_rel = f"{QUERY_RESULT_ROOT}/{folder}/{PARTITIONED_H2_FULL_FILE}"
            metadata_rel = f"{QUERY_RESULT_ROOT}/{folder}/{QUERY_METADATA_FILE}"
            delete_rel = f"{QUERY_RESULT_ROOT}/{folder}/{COEFFICIENT_DELETE_VALUES_FILE}"
            atomic_write_dataframe(
                _select_columns(query_summary, PARTITIONED_H2_COLUMNS, label="query summary"),
                query_dir / PARTITIONED_H2_SUMMARY_FILE,
                na_rep="NaN",
            )
            artifacts = per_query_artifacts.get(query_name)
            if artifacts is None:
                raise LDSCInternalError(
                    f"partitioned-h2 output writer is missing coefficient delete values for query {query_name!r}. "
                    "Most likely the fitted result family was assembled incompletely."
                )
            transfer = os.replace if artifacts.owned else shutil.copyfile
            transfer(artifacts.delete_values, query_dir / COEFFICIENT_DELETE_VALUES_FILE)
            transfer(artifacts.full, query_dir / PARTITIONED_H2_FULL_FILE)
            fit_metadata = json.loads(artifacts.metadata.read_text(encoding="utf-8"))
            payload = {
                **fit_metadata,
                "global_metadata": QUERY_GLOBAL_METADATA,
                "artifact_type": "partitioned_h2_query_result",
                "files": {"summary": summary_rel, "full": full_rel, "coefficient_delete_values": delete_rel},
                "ordinal": record["ordinal"],
                "query_annotation": query_name,
                "slug": record["slug"],
                "folder": folder,
            }
            atomic_write_json(payload, query_dir / QUERY_METADATA_FILE)
            manifest_rows.append(
                {
                    "ordinal": record["ordinal"],
                    "query_annotation": query_name,
                    "slug": record["slug"],
                    "folder": folder,
                    "summary_path": summary_rel,
                    "partitioned_h2_full_path": full_rel,
                    "metadata_path": metadata_rel,
                    "coefficient_delete_values_path": delete_rel,
                }
            )
        return manifest_rows


@dataclass(frozen=True)
class RgOutputConfig:
    """Directory-oriented output config for genetic-correlation summaries.

    Parameters
    ----------
    output_dir : str or os.PathLike[str]
        Directory that receives ``rg.tsv``, ``rg_full.tsv``, and
        ``h2_per_trait.tsv``.
    overwrite : bool, optional
        If ``True``, replace existing fixed rg outputs and remove stale owned
        siblings after a successful write. Default is ``False``.
    write_per_pair_detail : bool, optional
        If ``True``, also write one subdirectory per tested pair under
        ``diagnostics/pairs/``. Default is ``False``.
    """

    output_dir: str | PathLike[str]
    overwrite: bool = False
    write_per_pair_detail: bool = False

    def __post_init__(self) -> None:
        """Normalize the output directory."""
        object.__setattr__(self, "output_dir", normalize_required_path(self.output_dir))


class RgDirectoryWriter:
    """Write concise, full, per-trait, and optional per-pair rg outputs.

    The workflow-level logger owns ``diagnostics/rg.log``. This writer owns the public data
    artifacts plus diagnostic metadata in the rg output family:
    ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, and optional
    ``diagnostics/pairs/``.
    The optional pair tree is staged before replacement so a failed write does
    not expose a partially populated final tree. A successful overwrite also
    removes the default plot root derived from the superseded result.
    """

    @staticmethod
    def artifact_family(output_dir, *, include_pairs=False) -> ArtifactFamily:
        """Declare rg summaries, optional pair details, and the derived plot root."""
        return declare_artifacts(output_dir, {
            "metadata": ("diagnostics/metadata.json", True),
            "rg": ("rg.tsv", True),
            "rg_full": ("rg_full.tsv", True),
            "h2_per_trait": ("h2_per_trait.tsv", True),
            "pairs": ("diagnostics/pairs", include_pairs),
            "plots": ("plots", False),
        }, label="rg output artifact")

    def write(self, result: Any, output_config: RgOutputConfig) -> dict[str, str]:
        """Write an ``RgResultFamily``-like object to a result directory.

        Parameters
        ----------
        result : object
            Object exposing ``rg``, ``rg_full``, ``h2_per_trait``, and
            ``per_pair_metadata`` attributes. The first two tables must contain
            the canonical rg concise and full schemas.
        output_config : RgOutputConfig
            Output directory, overwrite policy, and optional per-pair detail
            setting.

        Returns
        -------
        dict of str to str
            Written path map. Always includes ``rg``, ``rg_full``, and
            ``h2_per_trait``; also includes ``pairs_root`` and
            ``pairs_manifest`` when per-pair detail is enabled.

        Raises
        ------
        ValueError
            If required result columns are missing or the pair metadata length
            does not match ``rg_full`` when detail output is requested.
        FileExistsError
            If an owned rg artifact already exists and overwrite is disabled.
        """
        rg = _select_columns(getattr(result, "rg"), RG_CONCISE_COLUMNS, label="rg summary")
        rg_full = _select_columns(getattr(result, "rg_full"), RG_FULL_COLUMNS, label="full rg summary")
        h2_per_trait = getattr(result, "h2_per_trait")
        per_pair_metadata = list(getattr(result, "per_pair_metadata", []))
        if output_config.write_per_pair_detail and len(per_pair_metadata) != len(rg_full):
            raise LDSCInternalError(
                "rg output writer could not write per-pair diagnostics because per_pair_metadata does not contain "
                "one record per rg_full row. Most likely the rg workflow returned inconsistent pair-detail "
                "metadata. Re-run with DEBUG logging and report the traceback."
            )

        output_dir = ensure_output_directory(output_config.output_dir, label="output directory")
        family = self.artifact_family(output_dir, include_pairs=output_config.write_per_pair_detail)
        metadata_path = family.paths["metadata"]
        diagnostics_dir = metadata_path.parent
        rg_path = family.paths["rg"]
        full_path = family.paths["rg_full"]
        h2_path = family.paths["h2_per_trait"]
        pairs_root = family.paths.get("pairs")
        stale_paths = family.preflight(overwrite=output_config.overwrite)
        paths = {name: str(path) for name, path in family.paths.items() if name != "pairs"}
        diagnostics_dir.mkdir(parents=True, exist_ok=True)
        if not output_config.write_per_pair_detail:
            atomic_write_dataframe(rg, rg_path, na_rep="NaN")
            atomic_write_dataframe(rg_full, full_path, na_rep="NaN")
            atomic_write_dataframe(h2_per_trait, h2_path, na_rep="NaN")
            atomic_write_json(_rg_root_metadata(result, files=family.metadata_files()), metadata_path)
            remove_output_artifacts(stale_paths)
            return paths

        pair_records = self._pair_records(rg_full, per_pair_metadata)
        staging_dir = Path(tempfile.mkdtemp(prefix=".pairs.tmp.", dir=str(diagnostics_dir)))
        backup_dir: Path | None = None
        try:
            manifest_rows = self._write_staged_pair_tree(staging_dir, pair_records, rg_full)
            atomic_write_dataframe(pd.DataFrame(manifest_rows), staging_dir / "manifest.tsv", na_rep="NaN")
            atomic_write_dataframe(rg, rg_path, na_rep="NaN")
            atomic_write_dataframe(rg_full, full_path, na_rep="NaN")
            atomic_write_dataframe(h2_per_trait, h2_path, na_rep="NaN")
            atomic_write_json(_rg_root_metadata(result, files=family.metadata_files()), metadata_path)
            if output_config.overwrite and pairs_root.exists():
                backup_dir = Path(tempfile.mkdtemp(prefix=".pairs.backup.", dir=str(diagnostics_dir)))
                backup_dir.rmdir()
                os.replace(pairs_root, backup_dir)
            os.replace(staging_dir, pairs_root)
        except Exception:
            if staging_dir.exists():
                shutil.rmtree(staging_dir, ignore_errors=True)
            if backup_dir is not None and backup_dir.exists() and not pairs_root.exists():
                os.replace(backup_dir, pairs_root)
            raise
        else:
            if backup_dir is not None and backup_dir.exists():
                shutil.rmtree(backup_dir, ignore_errors=True)

        remove_output_artifacts(stale_paths)
        paths["pairs_root"] = str(pairs_root)
        paths["pairs_manifest"] = str(pairs_root / "manifest.tsv")
        return paths

    def _pair_records(self, rg_full: pd.DataFrame, per_pair_metadata: list[dict[str, object]]) -> list[dict[str, object]]:
        """Return ordered pair records with deterministic safe folder names."""
        width = max(4, len(str(len(rg_full))))
        records: list[dict[str, object]] = []
        for idx, (_, row) in enumerate(rg_full.iterrows(), start=1):
            metadata = dict(per_pair_metadata[idx - 1])
            trait_1 = str(row["trait_1"])
            trait_2 = str(row["trait_2"])
            slug = f"{_slugify_rg_trait_name(trait_1)}_vs_{_slugify_rg_trait_name(trait_2)}"
            folder = f"{idx:0{width}d}_{slug}"
            records.append(
                {
                    "ordinal": idx,
                    "trait_1": trait_1,
                    "trait_2": trait_2,
                    "slug": slug,
                    "folder": folder,
                    "metadata": metadata,
                }
            )
        return records

    def _write_staged_pair_tree(
        self,
        staging_dir: Path,
        pair_records: list[dict[str, object]],
        rg_full: pd.DataFrame,
    ) -> list[dict[str, object]]:
        """Populate the staged pair result tree and return manifest rows."""
        manifest_rows: list[dict[str, object]] = []
        for record in pair_records:
            folder = str(record["folder"])
            pair_dir = staging_dir / folder
            pair_dir.mkdir(parents=True, exist_ok=False)
            row = rg_full.iloc[[int(record["ordinal"]) - 1]].reset_index(drop=True)
            detail_rel = f"diagnostics/pairs/{folder}/rg_full.tsv"
            metadata_rel = f"diagnostics/pairs/{folder}/metadata.json"
            atomic_write_dataframe(row, pair_dir / "rg_full.tsv", na_rep="NaN")
            payload = {
                **dict(record["metadata"]),
                "artifact_type": "rg_pair_result",
                "files": {"rg_full": detail_rel},
                "ordinal": record["ordinal"],
                "trait_1": record["trait_1"],
                "trait_2": record["trait_2"],
                "slug": record["slug"],
                "folder": folder,
            }
            atomic_write_json(payload, pair_dir / "metadata.json")
            manifest_rows.append(
                {
                    "ordinal": record["ordinal"],
                    "trait_1": record["trait_1"],
                    "trait_2": record["trait_2"],
                    "slug": record["slug"],
                    "folder": folder,
                    "rg_full_path": detail_rel,
                    "metadata_path": metadata_rel,
                }
            )
        return manifest_rows


def _rg_root_metadata(result: Any, *, files: dict[str, str]) -> dict[str, object]:
    """Build diagnostic metadata for an rg result directory from result tables."""
    rg_full = getattr(result, "rg_full", pd.DataFrame())
    h2_per_trait = getattr(result, "h2_per_trait", pd.DataFrame())
    trait_names: list[str] = []
    if "trait_name" in h2_per_trait.columns:
        trait_names = [str(value) for value in h2_per_trait["trait_name"].dropna().tolist()]
    pair_kind = None
    if "pair_kind" in rg_full.columns:
        values = [str(value) for value in rg_full["pair_kind"].dropna().unique().tolist()]
        pair_kind = values[0] if len(values) == 1 else ("mixed" if values else None)
    payload = result_metadata(
        {
            "trait_names": trait_names,
            "pair_kind": pair_kind,
            "n_pairs": int(len(rg_full)),
            "n_traits": int(len(trait_names)),
        },
        artifact_type="rg_result",
        files=files,
    )
    ldscore_dir = getattr(result, "ldscore_dir", None)
    if ldscore_dir is not None:
        payload["ldscore_dir"] = ldscore_dir
    sumstats_files = getattr(result, "sumstats_files", None)
    if sumstats_files is not None:
        payload["sumstats_files"] = list(sumstats_files)
    return payload


def _slugify_query_name(value: str) -> str:
    """Return a filesystem-safe slug for a query annotation name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "annotation"


def _slugify_rg_trait_name(value: str) -> str:
    """Return a filesystem-safe slug for an rg trait name."""
    normalized = unicodedata.normalize("NFKD", str(value)).encode("ascii", "ignore").decode("ascii")
    slug = re.sub(r"[^a-z0-9._-]+", "_", normalized.lower()).strip("._-")
    return slug or "trait"


def _select_columns(df: pd.DataFrame, columns: list[str], *, label: str) -> pd.DataFrame:
    """Return ``df`` with required public output columns in canonical order."""
    if df.empty and len(df.columns) == 0:
        return df
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise LDSCInternalError(
            f"output writer could not select public columns for {label} because required columns are missing: "
            f"{missing}. Most likely the upstream workflow returned a result table with the wrong schema. "
            "Re-run with DEBUG logging and report the traceback."
        )
    return df.loc[:, columns]


