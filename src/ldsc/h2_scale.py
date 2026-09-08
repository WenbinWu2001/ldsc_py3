"""Post-fit observed-to-liability-scale conversion for LDSC h2 results.

Core functionality:
    Convert the observed-scale estimate and block-jackknife standard error in a
    canonical ``h2`` result at either one population prevalence or an inclusive
    prevalence grid.

Overview
--------
The workflow always reads ``total_h2_obs`` and ``total_h2_obs_se`` from the
declared h2 summary. It delegates the numerical conversion to the LDSC kernel,
so post-processing uses the same formula as the regression shortcuts. Exact
conversion has no plotting dependency; range conversion imports Matplotlib
only when it is explicitly requested.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from ._kernel.regression import liability_conversion_factor
from ._logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from .errors import LDSCDependencyError, LDSCInputError, LDSCUsageError
from .path_resolution import (
    ensure_output_directory,
    preflight_output_artifact_family,
    remove_output_artifacts,
)


_TABLE_NAME = "h2_scale_conversion.tsv"
_PLOT_NAME = "h2_prevalence_sensitivity.png"
_METADATA_NAME = "diagnostics/metadata.json"
_LOG_NAME = "diagnostics/convert-h2-scale.log"
_TABLE_COLUMNS = [
    "trait_name",
    "samp_prev",
    "pop_prev",
    "conversion_factor",
    "total_h2_obs",
    "total_h2_obs_se",
    "total_h2_liab",
    "total_h2_liab_se",
]


@dataclass(frozen=True)
class H2ScaleConversionArtifact:
    """Files written by one liability-scale conversion.

    Parameters
    ----------
    table_path : pathlib.Path
        Tidy conversion table.
    metadata_path : pathlib.Path
        Machine-readable provenance for the derived result.
    log_path : pathlib.Path
        Human-readable workflow log.
    plot_path : pathlib.Path or None
        Sensitivity figure for range mode, otherwise ``None``.
    """

    table_path: Path
    metadata_path: Path
    log_path: Path
    plot_path: Path | None


@materializing_overwrite_guard(
    lambda h2_result_dir, **kwargs: (
        (
            kwargs.get("output_dir")
            or Path(h2_result_dir).expanduser() / "postprocessing" / "liability-scale",
            kwargs.get("overwrite", False),
            "RUN_FAILED.txt",
        )
    ),
    command="convert_h2_scale(...)",
)
def convert_h2_scale(
    h2_result_dir: str | os.PathLike[str],
    *,
    samp_prev: float,
    pop_prev: float | None = None,
    pop_prev_range: tuple[float, float] | None = None,
    num_points: int = 201,
    output_dir: str | os.PathLike[str] | None = None,
    overwrite: bool = False,
    log_level: str = "INFO",
) -> H2ScaleConversionArtifact:
    """Convert a saved observed-scale SNP-heritability estimate.

    Parameters
    ----------
    h2_result_dir : path-like
        Canonical unpartitioned h2 result directory.
    samp_prev : float
        Case fraction in the analyzed sample, strictly between zero and one.
    pop_prev : float, optional
        One population prevalence. Mutually exclusive with
        ``pop_prev_range``.
    pop_prev_range : tuple of float, optional
        Inclusive lower and upper population prevalences for a linear grid.
    num_points : int, optional
        Number of grid points in range mode. Default is 201.
    output_dir : path-like, optional
        Advanced Python-only destination. By default, write below
        ``<h2-result>/postprocessing/liability-scale``.
    overwrite : bool, optional
        Replace this conversion's fixed artifact family. Default is ``False``.
    log_level : {"DEBUG", "INFO", "WARNING", "ERROR"}, optional
        Workflow log threshold. Default is ``"INFO"``.

    Returns
    -------
    H2ScaleConversionArtifact
        Paths to the conversion table, metadata, log, and optional figure.

    Raises
    ------
    LDSCInputError
        If the source is not a current canonical h2 result or its observed-scale
        estimate and standard error are unavailable or invalid.
    LDSCUsageError
        If prevalence arguments are missing, mutually incompatible, outside
        ``(0, 1)``, or define an invalid sensitivity grid.
    LDSCDependencyError
        If sensitivity mode is requested without the optional Matplotlib
        dependency. Exact mode does not require Matplotlib.
    FileExistsError
        If a fixed conversion artifact exists and ``overwrite`` is ``False``.

    Notes
    -----
    Both the observed estimate and its block-jackknife standard error are
    multiplied by the same conversion factor. Population prevalence is treated
    as fixed; uncertainty in prevalence is not propagated.
    """
    source_dir = _require_directory(h2_result_dir, label="h2 result directory")
    source_metadata = _read_json(source_dir / "diagnostics" / "metadata.json")
    if source_metadata.get("artifact_type") != "h2_result":
        raise LDSCInputError(
            "convert-h2-scale requires a canonical h2 result with "
            f"artifact_type='h2_result'; got {source_metadata.get('artifact_type')!r}."
        )
    summary_path = _declared_source_file(source_dir, source_metadata, "summary")
    source_row = _read_observed_h2(summary_path)
    prevalences, mode = _population_prevalence_grid(
        pop_prev=pop_prev,
        pop_prev_range=pop_prev_range,
        num_points=num_points,
    )
    factors = np.asarray(liability_conversion_factor(samp_prev, prevalences), dtype=np.float64)

    pyplot = _load_pyplot() if mode == "sensitivity" else None
    destination = (
        Path(output_dir).expanduser()
        if output_dir is not None
        else source_dir / "postprocessing" / "liability-scale"
    )
    paths = {
        "table": destination / _TABLE_NAME,
        "plot": destination / _PLOT_NAME,
        "metadata": destination / _METADATA_NAME,
        "log": destination / _LOG_NAME,
    }
    produced = [paths["table"], paths["metadata"], paths["log"]]
    if mode == "sensitivity":
        produced.append(paths["plot"])
    stale = preflight_output_artifact_family(
        produced,
        paths.values(),
        overwrite=overwrite,
        label="liability-scale conversion artifact",
    )
    ensure_output_directory(destination, label="liability-scale conversion directory")
    ensure_output_directory(destination / "diagnostics", label="conversion diagnostics directory")

    table = pd.DataFrame(
        {
            "trait_name": [source_row["trait_name"]] * len(prevalences),
            "samp_prev": np.full(len(prevalences), float(samp_prev)),
            "pop_prev": prevalences,
            "conversion_factor": factors,
            "total_h2_obs": np.full(len(prevalences), source_row["total_h2_obs"]),
            "total_h2_obs_se": np.full(len(prevalences), source_row["total_h2_obs_se"]),
            "total_h2_liab": source_row["total_h2_obs"] * factors,
            "total_h2_liab_se": source_row["total_h2_obs_se"] * factors,
        },
        columns=_TABLE_COLUMNS,
    )
    files = {"table": _TABLE_NAME}
    if mode == "sensitivity":
        files["plot"] = _PLOT_NAME
    metadata = {
        "artifact_type": "h2_scale_conversion_result",
        "source_artifact_type": "h2_result",
        "source_result_dir": str(source_dir),
        "source_table": str(summary_path.relative_to(source_dir)),
        "mode": mode,
        "samp_prev": float(samp_prev),
        "pop_prev": float(prevalences[0]) if mode == "exact" else None,
        "pop_prev_range": (
            [float(prevalences[0]), float(prevalences[-1])]
            if mode == "sensitivity"
            else None
        ),
        "num_points": int(len(prevalences)),
        "files": files,
        "uncertainty": "block_jackknife_standard_error",
        "population_prevalence_uncertainty_propagated": False,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }

    with workflow_logging("convert-h2-scale", paths["log"], log_level=log_level):
        log_inputs(
            h2_result_dir=source_dir,
            summary=summary_path,
            samp_prev=samp_prev,
            pop_prev=pop_prev,
            pop_prev_range=pop_prev_range,
            num_points=len(prevalences),
        )
        _atomic_write_dataframe(table, paths["table"])
        if mode == "sensitivity":
            _plot_prevalence_sensitivity(table, paths["plot"], pyplot)
        _atomic_write_json(metadata, paths["metadata"])
        log_outputs(
            table=paths["table"],
            metadata=paths["metadata"],
            plot=paths["plot"] if mode == "sensitivity" else None,
        )
    remove_output_artifacts(stale)
    return H2ScaleConversionArtifact(
        table_path=paths["table"],
        metadata_path=paths["metadata"],
        log_path=paths["log"],
        plot_path=paths["plot"] if mode == "sensitivity" else None,
    )


def add_convert_h2_scale_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the public ``convert-h2-scale`` arguments to ``parser``."""
    parser.add_argument("--h2-result-dir", required=True, help="Canonical unpartitioned h2 result directory.")
    parser.add_argument("--samp-prev", required=True, type=float, help="Sample case fraction P in (0, 1).")
    prevalence = parser.add_mutually_exclusive_group(required=True)
    prevalence.add_argument("--pop-prev", type=float, help="One population prevalence K in (0, 1).")
    prevalence.add_argument(
        "--pop-prev-range",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        help="Inclusive population-prevalence range for a sensitivity curve.",
    )
    parser.add_argument("--num-points", type=int, default=201, help="Number of range grid points (default: 201).")
    parser.add_argument("--overwrite", action="store_true", help="Replace this conversion's fixed artifacts.")
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
        help="Workflow log threshold (default: INFO).",
    )


def run_convert_h2_scale_from_args(args: argparse.Namespace) -> H2ScaleConversionArtifact:
    """Run post-fit h2 conversion from parsed CLI arguments."""
    prevalence_range = tuple(args.pop_prev_range) if args.pop_prev_range is not None else None
    return convert_h2_scale(
        args.h2_result_dir,
        samp_prev=args.samp_prev,
        pop_prev=args.pop_prev,
        pop_prev_range=prevalence_range,
        num_points=args.num_points,
        overwrite=args.overwrite,
        log_level=args.log_level,
    )


def _require_directory(path: str | os.PathLike[str], *, label: str) -> Path:
    result = Path(path).expanduser()
    if not result.is_dir():
        raise LDSCInputError(f"{label} does not exist or is not a directory: '{result}'.")
    return result.resolve()


def _read_json(path: Path) -> dict:
    if not path.is_file():
        raise LDSCInputError(f"Required canonical metadata file is missing: '{path}'.")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LDSCInputError(f"Could not read canonical metadata JSON at '{path}': {exc}.") from exc
    if not isinstance(value, dict):
        raise LDSCInputError(f"Canonical metadata at '{path}' must contain a JSON object.")
    return value


def _declared_source_file(result_dir: Path, metadata: dict, key: str) -> Path:
    files = metadata.get("files")
    if not isinstance(files, dict) or not isinstance(files.get(key), str) or not files[key].strip():
        raise LDSCInputError(f"Canonical h2 metadata must declare files.{key}.")
    token = Path(files[key])
    if token.is_absolute():
        raise LDSCInputError(f"Canonical h2 metadata files.{key} must be relative to the result directory.")
    path = (result_dir / token).resolve()
    try:
        path.relative_to(result_dir)
    except ValueError as exc:
        raise LDSCInputError(f"Canonical h2 metadata files.{key} escapes the result directory.") from exc
    if not path.is_file():
        raise LDSCInputError(f"Canonical h2 metadata declares files.{key}='{files[key]}', but that file is missing.")
    return path


def _read_observed_h2(path: Path) -> dict[str, object]:
    try:
        table = pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError) as exc:
        raise LDSCInputError(f"Could not read h2 summary '{path}': {exc}.") from exc
    required = {"trait_name", "total_h2_obs", "total_h2_obs_se"}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise LDSCInputError(f"H2 summary '{path}' is missing required columns: {', '.join(missing)}.")
    if len(table) != 1:
        raise LDSCInputError(f"H2 summary '{path}' must contain exactly one row; found {len(table)}.")
    row = table.iloc[0]
    estimate = _finite_float(row["total_h2_obs"], "total_h2_obs")
    standard_error = _finite_float(row["total_h2_obs_se"], "total_h2_obs_se")
    if standard_error < 0:
        raise LDSCInputError("total_h2_obs_se must be non-negative.")
    return {
        "trait_name": str(row["trait_name"]),
        "total_h2_obs": estimate,
        "total_h2_obs_se": standard_error,
    }


def _finite_float(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise LDSCInputError(f"{label} must be numeric; got {value!r}.") from exc
    if not math.isfinite(result):
        raise LDSCInputError(f"{label} must be finite; got {value!r}.")
    return result


def _population_prevalence_grid(
    *,
    pop_prev: float | None,
    pop_prev_range: tuple[float, float] | None,
    num_points: int,
) -> tuple[np.ndarray, str]:
    if (pop_prev is None) == (pop_prev_range is None):
        raise LDSCUsageError("Specify exactly one of pop_prev or pop_prev_range.")
    if pop_prev is not None:
        value = _finite_float(pop_prev, "pop_prev")
        liability_conversion_factor(0.5, value)
        return np.asarray([value], dtype=np.float64), "exact"
    if len(pop_prev_range) != 2:
        raise LDSCUsageError("pop_prev_range must contain exactly two endpoints.")
    lower = _finite_float(pop_prev_range[0], "pop_prev_range lower endpoint")
    upper = _finite_float(pop_prev_range[1], "pop_prev_range upper endpoint")
    liability_conversion_factor(0.5, np.asarray([lower, upper]))
    if lower >= upper:
        raise LDSCUsageError("pop_prev_range must be strictly increasing.")
    if isinstance(num_points, bool) or not isinstance(num_points, int) or num_points < 2:
        raise LDSCUsageError("num_points must be an integer of at least 2 for a prevalence range.")
    return np.linspace(lower, upper, num_points, dtype=np.float64), "sensitivity"


def _load_pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as pyplot
    except ImportError as exc:
        raise LDSCDependencyError(
            "Prevalence-range conversion requires Matplotlib. Install the optional plotting dependency with "
            "`pip install 'ldsc[plot]'`, then rerun the command."
        ) from exc
    return pyplot


def _plot_prevalence_sensitivity(table: pd.DataFrame, path: Path, pyplot) -> None:
    figure, axis = pyplot.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    x = table["pop_prev"].to_numpy(dtype=float)
    estimate = table["total_h2_liab"].to_numpy(dtype=float)
    standard_error = table["total_h2_liab_se"].to_numpy(dtype=float)
    axis.fill_between(
        x,
        estimate - standard_error,
        estimate + standard_error,
        color="#4C78A8",
        alpha=0.2,
        linewidth=0,
        label="±1 jackknife SE",
    )
    axis.plot(x, estimate, color="#2F5D8A", linewidth=2, label="Estimate")
    axis.set_xlabel(r"Population prevalence ($K$)")
    axis.set_ylabel(r"Liability-scale SNP heritability ($h^2_{\mathrm{SNP}}$)")
    axis.set_title(f"{table['trait_name'].iloc[0]}\nSample prevalence ($P$) = {table['samp_prev'].iloc[0]:g}")
    axis.legend(loc="lower right", frameon=False)
    axis.spines[["top", "right"]].set_visible(False)
    figure.savefig(path, dpi=300, facecolor="white")
    pyplot.close(figure)


def _atomic_write_dataframe(table: pd.DataFrame, path: Path) -> None:
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    os.close(file_descriptor)
    temporary_path = Path(temporary_name)
    try:
        table.to_csv(temporary_path, sep="\t", index=False, na_rep="NaN")
        os.replace(temporary_path, path)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def _atomic_write_json(payload: dict[str, object], path: Path) -> None:
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
