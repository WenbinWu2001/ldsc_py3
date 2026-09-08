"""Optional plotting dispatcher for canonical LDSC result directories.

Core functionality:
    Select and save one approved exploratory figure from the scientific regime
    recorded in a canonical result's metadata.

Overview
--------
``plot_result`` is the only public plotting entry point. It validates the
source artifact and dispatches by metadata rather than by filenames or table
shape. Matplotlib is imported only when this module is asked to build a plot;
the numerical workflows do not depend on the plotting extra.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile
from typing import Any

import pandas as pd

from .._logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from ..errors import LDSCDependencyError, LDSCInputError
from ..path_resolution import ensure_output_directory, preflight_output_artifact_family


@dataclass(frozen=True)
class PlotArtifact:
    """One saved LDSC figure and its live Matplotlib objects.

    Parameters
    ----------
    kind : str
        Stable scientific plot identifier.
    path : pathlib.Path
        Saved 300-dpi PNG path.
    figure : matplotlib.figure.Figure
        Live figure for Python-level customization.
    axes : matplotlib.axes.Axes
        Primary scientific axes. Auxiliary axes, such as a weight scale, remain
        available through ``figure.axes``.
    """

    kind: str
    path: Path
    figure: Any
    axes: Any


@dataclass(frozen=True)
class _PlotContract:
    kind: str
    filename: str
    source_key: str
    uncertainty: str


@materializing_overwrite_guard(
    lambda result_dir, **kwargs: (
        (
            kwargs.get("output_dir") or Path(result_dir).expanduser() / "plots",
            kwargs.get("overwrite", False),
            "RUN_FAILED.txt",
        )
    ),
    command="plot_result(...)",
)
def plot_result(
    result_dir: str | os.PathLike[str],
    *,
    output_dir: str | os.PathLike[str] | None = None,
    overwrite: bool = False,
    log_level: str = "INFO",
) -> PlotArtifact:
    """Create the approved plot for a canonical LDSC result.

    Parameters
    ----------
    result_dir : path-like
        Root directory of a current canonical h2, rg, partitioned-h2, or
        quantile-h2 result.
    output_dir : path-like, optional
        Advanced Python-only destination. By default, write to
        ``<result-dir>/plots``.
    overwrite : bool, optional
        Replace the selected plot family's fixed PNG, metadata, and log.
        Default is ``False``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Workflow log threshold. Default is ``"INFO"``.

    Returns
    -------
    PlotArtifact
        Saved path and live Matplotlib figure and primary axes.

    Raises
    ------
    LDSCInputError
        If source metadata or the selected numerical table violates the
        plotting contract.
    LDSCDependencyError
        If the optional Matplotlib dependency is unavailable.
    FileExistsError
        If a fixed plot artifact exists and ``overwrite`` is ``False``.
    """
    source_dir = _require_result_directory(result_dir)
    metadata = _read_source_metadata(source_dir)
    contract = _select_contract(metadata)
    source_path = _declared_source_file(source_dir, metadata, contract.source_key)
    try:
        table = pd.read_csv(source_path, sep="\t")
    except (OSError, pd.errors.ParserError) as exc:
        raise LDSCInputError(f"Could not read plotting source table '{source_path}': {exc}.") from exc
    builders = _load_builders()
    destination = Path(output_dir).expanduser() if output_dir is not None else source_dir / "plots"
    paths = {
        "plot": destination / contract.filename,
        "metadata": destination / "diagnostics" / "metadata.json",
        "log": destination / "diagnostics" / "plot.log",
    }
    preflight_output_artifact_family(
        paths.values(),
        paths.values(),
        overwrite=overwrite,
        label=f"{contract.kind} plot artifact",
    )
    figure, axes = builders.build_plot(contract.kind, table, metadata)
    plot_metadata = {
        "artifact_type": "plot_result",
        "plot_kind": contract.kind,
        "source_artifact_type": metadata["artifact_type"],
        "source_result_dir": str(source_dir),
        "source_table": str(source_path.relative_to(source_dir)),
        "files": {"plot": contract.filename},
        "uncertainty": contract.uncertainty,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    try:
        ensure_output_directory(destination, label="plot output directory")
        ensure_output_directory(destination / "diagnostics", label="plot diagnostics directory")
        with workflow_logging("plot", paths["log"], log_level=log_level):
            log_inputs(result_dir=source_dir, source_table=source_path, plot_kind=contract.kind)
            figure.savefig(paths["plot"], dpi=300, bbox_inches="tight", facecolor="white")
            _atomic_write_json(plot_metadata, paths["metadata"])
            log_outputs(plot=paths["plot"], metadata=paths["metadata"])
    except Exception:
        builders.close_figure(figure)
        raise
    return PlotArtifact(kind=contract.kind, path=paths["plot"], figure=figure, axes=axes)


def add_plot_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the public ``plot`` arguments to ``parser``."""
    parser.add_argument("--result-dir", required=True, help="Canonical LDSC result root to plot.")
    parser.add_argument("--overwrite", action="store_true", help="Replace the selected fixed plot family.")
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Workflow log threshold (default: INFO).",
    )


def run_plot_from_args(args: argparse.Namespace) -> PlotArtifact:
    """Run plotting from parsed CLI arguments and close the command's figure."""
    artifact = plot_result(
        args.result_dir,
        overwrite=args.overwrite,
        log_level=args.log_level,
    )
    _load_builders().close_figure(artifact.figure)
    return artifact


def _require_result_directory(path: str | os.PathLike[str]) -> Path:
    result = Path(path).expanduser()
    if not result.is_dir():
        raise LDSCInputError(f"Plot result directory does not exist or is not a directory: '{result}'.")
    return result.resolve()


def _read_source_metadata(result_dir: Path) -> dict[str, Any]:
    path = result_dir / "diagnostics" / "metadata.json"
    if not path.is_file():
        query_metadata = result_dir / "metadata.json"
        if query_metadata.is_file():
            try:
                query_payload = json.loads(query_metadata.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                query_payload = {}
            if query_payload.get("artifact_type") == "partitioned_h2_query_result":
                raise LDSCInputError(
                    "Per-query partitioned-h2 result directories are not plotting inputs. "
                    "Pass the root partitioned-h2 result directory that summarizes all query annotations."
                )
        raise LDSCInputError(f"Canonical plotting metadata is missing: '{path}'.")
    try:
        metadata = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LDSCInputError(f"Could not read canonical plotting metadata '{path}': {exc}.") from exc
    if not isinstance(metadata, dict):
        raise LDSCInputError(f"Canonical plotting metadata '{path}' must contain a JSON object.")
    if not isinstance(metadata.get("artifact_type"), str):
        raise LDSCInputError(f"Canonical plotting metadata '{path}' must declare artifact_type.")
    return metadata


def _select_contract(metadata: dict[str, Any]) -> _PlotContract:
    artifact_type = metadata["artifact_type"]
    if artifact_type == "h2_result":
        files = metadata.get("files")
        if not isinstance(files, dict) or not files.get("ld_score_regression_bins"):
            raise LDSCInputError(
                "This h2 result does not contain diagnostics/ld_score_regression_bins.tsv. "
                "The binned regression plot cannot be reconstructed faithfully; rerun the current `ldsc h2` command."
            )
        return _PlotContract(
            "ld_score_regression",
            "ld_score_regression.png",
            "ld_score_regression_bins",
            "saved_bin_dispersion_not_plotted",
        )
    if artifact_type == "rg_result":
        pair_kind = metadata.get("pair_kind")
        if pair_kind == "all_pairs":
            return _PlotContract("rg_heatmap", "rg_heatmap.png", "rg", "block_jackknife_standard_error")
        if pair_kind == "anchor":
            return _PlotContract(
                "rg_anchor_forest", "rg_anchor_forest.png", "rg", "block_jackknife_standard_error"
            )
        raise LDSCInputError(
            "rg plotting requires metadata pair_kind='all_pairs' or pair_kind='anchor'; "
            f"got {pair_kind!r}."
        )
    if artifact_type == "partitioned_h2_result":
        analysis_type = metadata.get("analysis_type")
        headline_metric = metadata.get("headline_metric")
        if analysis_type == "functional_category" and headline_metric == "enrichment":
            return _PlotContract(
                "functional_h2_enrichment",
                "functional_h2_enrichment.png",
                "summary",
                "block_jackknife_standard_error",
            )
        if analysis_type == "cell_type_specific" and headline_metric == "coefficient":
            return _PlotContract(
                "cell_type_query_pvalues",
                "cell_type_query_pvalues.png",
                "summary",
                "nominal_one_sided_conditional_coefficient_p_value",
            )
        raise LDSCInputError(
            "Unsupported partitioned-h2 plotting regime: "
            f"analysis_type={analysis_type!r}, headline_metric={headline_metric!r}."
        )
    if artifact_type == "quantile_h2_result":
        if not isinstance(metadata.get("target_annotation"), str) or not metadata["target_annotation"].strip():
            raise LDSCInputError("quantile-h2 plotting metadata must declare target_annotation.")
        return _PlotContract(
            "continuous_annotation_quantile_enrichment",
            "continuous_annotation_quantile_enrichment.png",
            "quantile_h2",
            "block_jackknife_standard_error",
        )
    raise LDSCInputError(f"Artifact type {artifact_type!r} is not supported by `ldsc plot`.")


def _declared_source_file(result_dir: Path, metadata: dict[str, Any], key: str) -> Path:
    files = metadata.get("files")
    if not isinstance(files, dict) or not isinstance(files.get(key), str) or not files[key].strip():
        raise LDSCInputError(f"Canonical plotting metadata must declare files.{key}.")
    token = Path(files[key])
    if token.is_absolute():
        raise LDSCInputError(f"Canonical plotting metadata files.{key} must be relative to the result directory.")
    path = (result_dir / token).resolve()
    try:
        path.relative_to(result_dir)
    except ValueError as exc:
        raise LDSCInputError(f"Canonical plotting metadata files.{key} escapes the result directory.") from exc
    if not path.is_file():
        raise LDSCInputError(f"Canonical plotting metadata declares files.{key}='{files[key]}', but it is missing.")
    return path


def _load_builders():
    try:
        from . import _builders
    except ImportError as exc:
        if exc.name == "matplotlib" or (exc.name and exc.name.startswith("matplotlib.")):
            raise LDSCDependencyError(
                "Plotting requires Matplotlib. Install the optional dependency with "
                "`pip install 'ldsc[plot]'`, then rerun the command."
            ) from exc
        raise
    return _builders


def _atomic_write_json(payload: dict[str, Any], path: Path) -> None:
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(file_descriptor)
    temporary_path = Path(temporary_name)
    try:
        temporary_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


__all__ = ["PlotArtifact", "plot_result"]
