"""Plotting dispatcher for canonical LDSC result directories.

Core functionality:
    Select and save one approved exploratory figure from the scientific regime
    recorded in a canonical result's metadata.

Overview
--------
``plot_result`` is the only public plotting entry point. It validates the
source artifact and dispatches by metadata rather than by filenames or table
shape. Matplotlib is imported only when this module is asked to build a plot;
numerical workflows do not import the plotting runtime.

Both rg views also read the declared single-trait heritability table, when
available, to annotate observed-scale h2 and jackknife SE. Missing values stay
explicit in the figure; no estimate is refitted or converted by plotting.
"""

from __future__ import annotations

from .._cli_help import CLIHelpFormatter
from .._logging import LOG_LEVEL_HELP
from .._result_files import atomic_write_json, declared_result_file

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
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
        available through ``figure.axes``. In rg plots, the heritability text
        belongs to these primary axes, including the anchor plot's external
        text column; that column has no separate numerical axis.
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
        If the required Matplotlib dependency cannot be imported from the
        active environment.
    FileExistsError
        If a fixed plot artifact exists and ``overwrite`` is ``False``.

    Notes
    -----
    Both rg plots annotate saved single-trait observed-scale heritability and
    block-jackknife SE from the table declared by ``files.h2_per_trait`` in
    source metadata. Rows are matched by ``trait_name`` using ``total_h2_obs``
    and ``total_h2_obs_se``; pair-fit and liability-scale columns are not used. The
    heatmap places these values on the diagonal; the anchor plot uses a text
    column for partners and a subtitle for the anchor. Values use two decimal
    places and finite heritabilities are not clipped to [0, 1]. Missing sources
    or trait rows display ``failed``, as do nonnumeric/nonfinite estimates or
    SEs and negative SEs. Zero SE is valid. Extra traits are omitted from the
    figure; malformed tables, missing required columns, duplicate or empty
    trait names, and unsafe source paths raise ``LDSCInputError``. Plotting
    never refits. Derived metadata records these annotations separately under
    ``heritability_annotations`` while retaining the correlation source table.
    """
    source_dir = _require_result_directory(result_dir)
    metadata = _read_source_metadata(source_dir)
    contract = _select_contract(metadata)
    source_path = declared_result_file(source_dir, metadata, contract.source_key, context="Canonical plotting metadata")
    try:
        table = pd.read_csv(source_path, sep="\t")
    except (OSError, pd.errors.ParserError) as exc:
        raise LDSCInputError(f"Could not read plotting source table '{source_path}': {exc}.") from exc
    h2_per_trait = None
    h2_path = None
    if metadata["artifact_type"] == "rg_result" and "h2_per_trait" in metadata["files"]:
        h2_path = declared_result_file(
            source_dir, metadata, "h2_per_trait", context="Canonical rg plotting metadata", allow_missing=True
        )
        try:
            h2_per_trait = pd.read_csv(h2_path, sep="\t", dtype={"trait_name": str})
        except FileNotFoundError:
            pass
        except (OSError, ValueError) as exc:
            raise LDSCInputError(f"Could not read rg heritability table '{h2_path}': {exc}.") from exc
        if h2_per_trait is not None and not h2_per_trait.empty and not isinstance(h2_per_trait.index, pd.RangeIndex):
            raise LDSCInputError(
                f"Malformed rg heritability table '{h2_path}': rows contain more fields than the header. "
                "Supply a tab-separated table with one field per declared column."
            )
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
    figure, axes = builders.build_plot(contract.kind, table, metadata, h2_per_trait=h2_per_trait)
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
    if metadata["artifact_type"] == "rg_result":
        plot_metadata["heritability_annotations"] = {
            "source_table": str(h2_path.relative_to(source_dir)) if h2_path is not None else None,
            "source_available": h2_per_trait is not None,
            "scale": "observed",
            "estimate_column": "total_h2_obs",
            "se_column": "total_h2_obs_se",
            "uncertainty": "block_jackknife_standard_error",
            "missing_label": "failed",
        }
    try:
        ensure_output_directory(destination, label="plot output directory")
        ensure_output_directory(destination / "diagnostics", label="plot diagnostics directory")
        with workflow_logging("plot", paths["log"], log_level=log_level):
            log_inputs(result_dir=source_dir, source_table=source_path, plot_kind=contract.kind)
            if metadata["artifact_type"] == "rg_result":
                log_inputs(h2_per_trait=h2_path, h2_source_available=h2_per_trait is not None, h2_scale="observed")
            figure.savefig(paths["plot"], dpi=300, bbox_inches="tight", facecolor="white")
            atomic_write_json(plot_metadata, paths["metadata"])
            log_outputs(plot=paths["plot"], metadata=paths["metadata"])
    except Exception:
        builders.close_figure(figure)
        raise
    return PlotArtifact(kind=contract.kind, path=paths["plot"], figure=figure, axes=axes)


def add_plot_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the public ``plot`` arguments to ``parser``."""
    parser.prog = 'ldsc plot'
    parser.formatter_class = CLIHelpFormatter
    parser.description = 'Create a plot from a saved LDSC result directory.'
    inputs = parser.add_argument_group('Input result')
    runtime = parser.add_argument_group('Output and logging')

    inputs.add_argument(
        '--result-dir', required=True, metavar='DIR',
        help=(
            'Required saved LDSC result directory to plot. The result type selects the plot automatically; '
            'files are written below plots/ in this directory.'
        ),
    )

    runtime.add_argument(
        '--overwrite', action='store_true',
        help=(
            'Replace existing plot files for the selected result type. Default: off; stop if these output '
            'files exist.'
        ),
    )
    runtime.add_argument(
        '--log-level', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'), default='INFO',
        help=LOG_LEVEL_HELP,
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




def _load_builders():
    try:
        from . import _builders
    except ImportError as exc:
        if exc.name == "matplotlib" or (exc.name and exc.name.startswith("matplotlib.")):
            raise LDSCDependencyError(
                "Plotting requires Matplotlib, which is a required LDSC dependency. "
                "Repair the active environment with `pip install 'matplotlib>=3.9,<4'`, "
                "then rerun the command."
            ) from exc
        raise
    return _builders




__all__ = ["PlotArtifact", "plot_result"]
