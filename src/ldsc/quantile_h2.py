"""Post-fit heritability projection for continuous annotation quantiles.

The numerical functions operate on fitted partitioned-LDSC coefficients and
annotation sufficient statistics. They do not refit the regression model and
do not assign a coefficient to an external target annotation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import logging
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from ._cli_help import CLIHelpFormatter, CHROMOSOME_PATH_HELP
from ._logging import LOG_LEVEL_HELP
from .errors import LDSCInputError
from ._annotation_storage import AnnotationWorkspace
from ._quantile_inputs import prepare_quantile_statistics
from ._quantile_storage import AlignmentDiagnostics
from .outputs import (
    QuantileH2DirectoryWriter,
    QuantileH2OutputConfig,
)
from .path_resolution import (
    ensure_output_directory,
)
from ._logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging


LOGGER = logging.getLogger("LDSC.quantile_h2")

QUANTILE_H2_COLUMNS = [
    "quantile",
    "target_value_lower",
    "target_value_upper",
    "n_snps",
    "prop_snps",
    "h2_obs",
    "h2_obs_se",
    "h2_liab",
    "h2_liab_se",
    "prop_h2",
    "prop_h2_se",
    "enrichment",
    "enrichment_se",
    "enrichment_p",
]

STANDARDIZED_COEFFICIENT_COLUMNS = [
    "annotation",
    "annotation_type",
    "annotation_sd",
    "tau",
    "tau_se",
    "tau_z",
    "tau_p",
    "tau_star",
    "tau_star_se",
    "tau_star_z",
    "tau_star_p",
]


@dataclass(frozen=True)
class QuantileAssignment:
    """Legacy-compatible quantile membership and realized boundaries.

    Attributes
    ----------
    quantile : (n_snps,) ndarray of int
        One-based quantile labels in the original SNP order.
    lower, upper : (n_quantiles,) ndarray
        Realized target-value bounds ordered from low to high.
    counts : (n_quantiles,) ndarray of int
        Number of eligible SNPs assigned to each quantile.
    """

    quantile: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    counts: np.ndarray


def assign_legacy_quantiles(values: Sequence[float] | np.ndarray, n_quantiles: int = 5) -> QuantileAssignment:
    """Assign finite target values with the LDSC2 continuous-annotation rule.

    Boundary index ``i`` is ``floor(i * (n - 1) / Q + 0.5)``. The first
    interval includes both bounds; subsequent intervals exclude the lower and
    include the upper bound, so ties remain in the lower-valued quantile.

    Parameters
    ----------
    values : (n_snps,) array_like
        Finite numeric target values after explicit missing-token exclusions.
    n_quantiles : int, optional
        Requested number of quantiles. Must be at least two. Default is 5.

    Returns
    -------
    QuantileAssignment
        Membership, bounds, and counts in low-to-high order.

    Raises
    ------
    LDSCInputError
        If values are empty or nonfinite, fewer than two quantiles are
        requested, or the realized boundaries produce an empty quantile.
    """
    if int(n_quantiles) != n_quantiles or n_quantiles < 2:
        raise LDSCInputError("quantile-h2 requires --num-quantiles to be an integer of at least 2.")
    array = np.asarray(values, dtype=np.float64).reshape(-1)
    if not len(array) or not np.isfinite(array).all():
        raise LDSCInputError("quantile-h2 target values must be nonempty, numeric, and finite after exclusions.")
    sorted_values = np.sort(array, kind="stable")
    boundary_indices = np.floor(np.arange(n_quantiles + 1) * (len(array) - 1) / n_quantiles + 0.5).astype(int)
    boundaries = sorted_values[boundary_indices]
    labels = np.searchsorted(boundaries[1:], array, side="left") + 1
    counts = np.bincount(labels, minlength=n_quantiles + 1)[1:]
    if np.any(counts == 0):
        empty = (np.flatnonzero(counts == 0) + 1).tolist()
        raise LDSCInputError(
            f"quantile-h2 produced an empty quantile: {empty}. Target-value ties do not provide enough "
            "distinct realized intervals; reduce --num-quantiles or use a less discrete target annotation."
        )
    return QuantileAssignment(
        quantile=labels.astype(np.int64),
        lower=boundaries[:-1].astype(np.float64),
        upper=boundaries[1:].astype(np.float64),
        counts=counts.astype(np.int64),
    )


def compute_quantile_h2(
    *,
    annotation_sums: np.ndarray,
    tau: np.ndarray,
    tau_delete: np.ndarray,
    snp_counts: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
    liability_factor: float | None = None,
) -> pd.DataFrame:
    """Project one fitted joint model onto target-defined quantiles.

    Parameters
    ----------
    annotation_sums : (n_annotations, n_quantiles) ndarray
        Within-quantile sums for fitted annotations only.
    tau : (n_annotations,) ndarray
        Whole-data fitted coefficients in the same annotation order.
    tau_delete : (n_blocks, n_annotations) ndarray
        Delete-one-block coefficient values in fitted order.
    snp_counts : (n_quantiles,) ndarray
        Fixed eligible-SNP counts.
    lower, upper : (n_quantiles,) ndarray
        Realized target-value bounds.
    liability_factor : float or None, optional
        Observed-to-liability multiplier for absolute heritability. ``None``
        leaves liability columns missing.

    Returns
    -------
    pandas.DataFrame
        Quantile summaries in :data:`QUANTILE_H2_COLUMNS` order. The enrichment
        p-value is NaN when the inside-versus-complement contrast has zero or
        missing jackknife SE.
    """
    sums = np.asarray(annotation_sums, dtype=np.float64)
    point_tau = np.asarray(tau, dtype=np.float64).reshape(-1)
    deletes_tau = np.asarray(tau_delete, dtype=np.float64)
    counts = np.asarray(snp_counts, dtype=np.int64).reshape(-1)
    if sums.ndim != 2 or sums.shape != (len(point_tau), len(counts)):
        raise LDSCInputError(
            "quantile-h2 annotation-sum dimensions do not match fitted coefficients and quantile counts."
        )
    if deletes_tau.ndim != 2 or deletes_tau.shape[1] != len(point_tau) or not len(deletes_tau):
        raise LDSCInputError("quantile-h2 requires a nonempty delete-value matrix in fitted coefficient order.")
    if np.any(counts <= 0) or counts.sum() <= 0:
        raise LDSCInputError("quantile-h2 requires every quantile to contain at least one SNP.")

    h2 = sums.T @ point_tau
    h2_delete = deletes_tau @ sums
    h2_se = _jackknife_se(h2, h2_delete)
    total = float(h2.sum())
    delete_totals = h2_delete.sum(axis=1)
    prop_snps = counts / counts.sum()

    prop_h2 = np.full(len(counts), np.nan)
    prop_h2_se = np.full(len(counts), np.nan)
    enrichment = np.full(len(counts), np.nan)
    enrichment_se = np.full(len(counts), np.nan)
    if total != 0:
        prop_h2 = h2 / total
        enrichment = prop_h2 / prop_snps
        if np.all(delete_totals != 0):
            prop_delete = h2_delete / delete_totals[:, None]
            prop_h2_se = _jackknife_se(prop_h2, prop_delete)
            enrichment_se = prop_h2_se / prop_snps

    complement_counts = counts.sum() - counts
    contrast = h2 / counts - (total - h2) / complement_counts
    contrast_delete = h2_delete / counts - (delete_totals[:, None] - h2_delete) / complement_counts
    contrast_se = _jackknife_se(contrast, contrast_delete)
    contrast_z = np.divide(
        contrast, contrast_se,
        out=np.full_like(contrast, np.nan), where=contrast_se > 0,
    )
    enrichment_p = 2.0 * stats.norm.sf(np.abs(contrast_z))

    factor = float("nan") if liability_factor is None else float(liability_factor)
    result = pd.DataFrame(
        {
            "quantile": np.arange(1, len(counts) + 1, dtype=np.int64),
            "target_value_lower": np.asarray(lower, dtype=np.float64),
            "target_value_upper": np.asarray(upper, dtype=np.float64),
            "n_snps": counts,
            "prop_snps": prop_snps,
            "h2_obs": h2,
            "h2_obs_se": h2_se,
            "h2_liab": h2 * factor,
            "h2_liab_se": h2_se * factor,
            "prop_h2": prop_h2,
            "prop_h2_se": prop_h2_se,
            "enrichment": enrichment,
            "enrichment_se": enrichment_se,
            "enrichment_p": enrichment_p,
        }
    )
    return result.loc[:, QUANTILE_H2_COLUMNS]


def compute_standardized_coefficients(
    *,
    annotation_names: Sequence[str],
    annotation_types: Mapping[str, str],
    annotation_sd: np.ndarray,
    tau: np.ndarray,
    tau_se: np.ndarray,
    total_h2: float,
    n_common_snps: int,
) -> pd.DataFrame:
    """Compute paper-compatible fixed-scale standardized coefficients.

    Population standard deviations (``ddof=0``) and the whole-model
    heritability denominator must be computed over the complete common
    reference-SNP universe. The denominator is not ratio-jackknifed.
    """
    names = [str(name) for name in annotation_names]
    sd = np.asarray(annotation_sd, dtype=np.float64).reshape(-1)
    point_tau = np.asarray(tau, dtype=np.float64).reshape(-1)
    se = np.asarray(tau_se, dtype=np.float64).reshape(-1)
    if not (len(names) == len(sd) == len(point_tau) == len(se)):
        raise LDSCInputError("quantile-h2 standardized-coefficient inputs have inconsistent lengths.")
    tau_z = np.divide(point_tau, se, out=np.full_like(point_tau, np.nan), where=se != 0)
    tau_p = 2.0 * stats.norm.sf(np.abs(tau_z))
    tau_star = np.full(len(names), np.nan)
    tau_star_se = np.full(len(names), np.nan)
    tau_star_z = np.full(len(names), np.nan)
    tau_star_p = np.full(len(names), np.nan)
    if total_h2 > 0 and n_common_snps > 0:
        scale = sd / (float(total_h2) / int(n_common_snps))
        tau_star = point_tau * scale
        tau_star_se = se * scale
        tau_star_z = np.divide(tau_star, tau_star_se, out=np.full_like(tau_star, np.nan), where=tau_star_se != 0)
        tau_star_p = tau_p.copy()
    result = pd.DataFrame(
        {
            "annotation": names,
            "annotation_type": [annotation_types.get(name, "unknown") for name in names],
            "annotation_sd": sd,
            "tau": point_tau,
            "tau_se": se,
            "tau_z": tau_z,
            "tau_p": tau_p,
            "tau_star": tau_star,
            "tau_star_se": tau_star_se,
            "tau_star_z": tau_star_z,
            "tau_star_p": tau_star_p,
        }
    )
    return result.loc[:, STANDARDIZED_COEFFICIENT_COLUMNS]


def _jackknife_se(point: np.ndarray, delete_values: np.ndarray) -> np.ndarray:
    """Return block-jackknife SEs from whole-data and delete-one-block values."""
    point = np.asarray(point, dtype=np.float64).reshape(1, -1)
    deletes = np.asarray(delete_values, dtype=np.float64)
    if not np.isfinite(deletes).all():
        return np.full(point.shape[1], np.nan)
    n_blocks = deletes.shape[0]
    pseudovalues = n_blocks * point - (n_blocks - 1) * deletes
    covariance = np.atleast_2d(np.cov(pseudovalues.T, ddof=1) / n_blocks)
    return np.sqrt(np.maximum(0.0, np.diag(covariance)))


@dataclass(frozen=True)
class FittedPartitionedModel:
    """One complete partitioned-LDSC fitted model loaded from disk.

    Coefficient and delete-value arrays always follow ``annotation_names``.
    ``model_type`` distinguishes a baseline-only joint fit from one
    baseline-plus-query joint fit; it never represents an aggregate of query
    regressions.
    """

    result_dir: Path
    model_type: str
    metadata: dict[str, object]
    annotation_names: list[str]
    summary: pd.DataFrame
    tau: np.ndarray
    tau_se: np.ndarray
    tau_delete: np.ndarray
    ldscore_dir: Path


@dataclass(frozen=True)
class QuantileH2Result:
    """Compact result family returned by :func:`run_quantile_h2_from_args`.

    Attributes
    ----------
    quantile_h2 : pandas.DataFrame
        Low-to-high target quantile summaries.
    standardized_coefficients : pandas.DataFrame
        Raw and standardized coefficients in fitted annotation order.
    alignment_issues_path : str
        Persistent compressed table of every recorded alignment issue.
    metadata : dict
        Model, input, universe, missingness, and statistic provenance.
    output_paths : dict of str to str
        Written result-family paths.
    """

    quantile_h2: pd.DataFrame
    standardized_coefficients: pd.DataFrame
    alignment_issues_path: str
    metadata: dict[str, object]
    output_paths: dict[str, str]


def load_fitted_partitioned_model(result_dir: str | Path) -> FittedPartitionedModel:
    """Load one baseline-only or per-query LDSC3 partitioned-h2 fitted model.

    Parameters
    ----------
    result_dir : str or pathlib.Path
        Baseline-only partitioned-h2 output root or one
        ``diagnostics/query_annotations/<query>`` directory.

    Returns
    -------
    FittedPartitionedModel
        Ordered point coefficients, SEs, delete values, summary rows, and the
        linked LD-score directory.

    Raises
    ------
    LDSCInputError
        If the path is not one complete current LDSC3 model, is an aggregate
        multi-query root, lacks coefficient delete values, or has inconsistent
        fitted annotation order.
    """
    root = Path(result_dir)
    per_query_metadata = root / "metadata.json"
    root_metadata = root / "diagnostics" / "metadata.json"
    if per_query_metadata.is_file():
        metadata_path = per_query_metadata
        model_type = "baseline_plus_query"
        summary_path = root / "partitioned_h2_full.tsv"
        delete_path = root / "coefficient_delete_values.parquet"
    elif root_metadata.is_file():
        metadata_path = root_metadata
        model_type = "baseline_only"
        summary_path = root / "partitioned_h2.tsv"
        delete_path = root / "diagnostics" / "coefficient_delete_values.parquet"
    else:
        raise LDSCInputError(
            f"quantile-h2 could not load fitted model directory '{root}': no LDSC3 partitioned-h2 metadata was found. "
            "Pass a baseline-only partitioned-h2 output root or one diagnostics/query_annotations/<query> directory. "
            "Raw LDSC2 regression prefixes are not accepted."
        )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if model_type == "baseline_only" and metadata.get("analysis_type") == "cell_type_specific":
        raise LDSCInputError(
            "quantile-h2 cannot consume an aggregate multi-query partitioned-h2 root because its rows come from "
            "separate regressions. Select one diagnostics/query_annotations/<query> fitted-model directory."
        )
    if not delete_path.is_file():
        raise LDSCInputError(
            f"quantile-h2 requires '{delete_path.name}', but it is missing from '{root}'. Re-run partitioned-h2 "
            "with the current LDSC3 package; LD scores do not need to be recomputed."
        )
    if not summary_path.is_file():
        raise LDSCInputError(f"quantile-h2 fitted model is missing summary table '{summary_path}'.")
    delete_frame = pd.read_parquet(delete_path)
    if "delete_block" not in delete_frame.columns:
        raise LDSCInputError("quantile-h2 coefficient delete artifact lacks the required delete_block column.")
    annotation_names = [str(column) for column in delete_frame.columns if column != "delete_block"]
    retained = [str(name) for name in metadata.get("retained_ld_columns", annotation_names)]
    if retained != annotation_names:
        raise LDSCInputError(
            "quantile-h2 fitted-model coefficient order disagrees with metadata retained_ld_columns. "
            "Re-run partitioned-h2 to regenerate a consistent result directory."
        )
    summary = pd.read_csv(summary_path, sep="\t")
    indexed = summary.set_index(summary["category"].astype(str))
    missing = [name for name in annotation_names if name not in indexed.index]
    if missing:
        raise LDSCInputError(f"quantile-h2 fitted summary is missing retained annotation row(s): {missing}.")
    ordered = indexed.loc[annotation_names].reset_index(drop=True)
    ldscore_value = metadata.get("ldscore_dir")
    if not ldscore_value:
        raise LDSCInputError("quantile-h2 fitted-model metadata is missing ldscore_dir provenance.")
    ldscore_dir = Path(str(ldscore_value))
    if not ldscore_dir.is_dir():
        candidate = (root / ldscore_dir).resolve()
        if candidate.is_dir():
            ldscore_dir = candidate
    if not ldscore_dir.is_dir():
        raise LDSCInputError(
            f"quantile-h2 cannot resolve linked LD-score directory '{ldscore_value}'. Restore that directory "
            "or update the fitted-model metadata to its current location."
        )
    return FittedPartitionedModel(
        result_dir=root,
        model_type=model_type,
        metadata=metadata,
        annotation_names=annotation_names,
        summary=ordered,
        tau=pd.to_numeric(ordered["coefficient"], errors="raise").to_numpy(dtype=np.float64),
        tau_se=pd.to_numeric(ordered["coefficient_se"], errors="raise").to_numpy(dtype=np.float64),
        tau_delete=delete_frame.loc[:, annotation_names].to_numpy(dtype=np.float64),
        ldscore_dir=ldscore_dir,
    )


def add_quantile_h2_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the public ``quantile-h2`` command-line contract."""
    parser.allow_abbrev = False
    parser.prog = 'ldsc quantile-h2'
    parser.formatter_class = CLIHelpFormatter
    parser.description = 'Summarize an existing fitted model across SNP quantiles of a target annotation, without refitting LDSC.'
    inputs = parser.add_argument_group('Inputs and output')
    original = parser.add_argument_group('Original fitted-model inputs', description='Resupply the same annotations and reference SNP metadata used for the saved model.')
    quantiles = parser.add_argument_group('Target quantiles')
    genes = parser.add_argument_group('Original BED and gene-list settings')
    runtime = parser.add_argument_group('Output and logging')

    inputs.add_argument(
        '--partitioned-h2-result-dir', required=True, metavar='DIR',
        help=(
            'Required saved joint model: a baseline-only partitioned-h2 result or one '
            'diagnostics/query_annotations/<query>/ directory. An aggregate multi-query result root is not '
            'a single model.'
        ),
    )
    inputs.add_argument(
        '--target-annot-sources', nargs='+', required=True, metavar='SOURCES',
        help=(
            'Required annotation files containing --target-annotation, used to divide SNPs into quantiles. '
            'Accepts exact paths. The target may be an external score or a fitted annotation. '
            + CHROMOSOME_PATH_HELP
        ),
    )
    inputs.add_argument(
        '--target-annotation', required=True,
        help=(
            'Required unique column name in --target-annot-sources whose values define the quantiles.'
        ),
    )
    inputs.add_argument(
        '--output-dir', required=True, metavar='DIR',
        help=(
            'Required destination for quantile heritability, standardized coefficients, and diagnostics.'
        ),
    )

    original.add_argument(
        '--baseline-annot-sources', nargs='+', required=True, metavar='SOURCES',
        help=(
            "Required original baseline annotation files used to fit the selected model. Accepts exact paths, '*' "
            "and '@' patterns; see --target-annot-sources for pattern rules. Resupply their original values and "
            'SNP coverage.'
        ),
    )
    original.add_argument(
        '--ref-metadata-sources', nargs='+', required=True, metavar='SOURCES',
        help=(
            'Required reference SNP metadata files with chromosome, position, SNP ID, and MAF; allele-aware '
            "matching also uses A1/A2. Accepts exact paths, '*' and '@' patterns; see --target-annot-sources for "
            'pattern rules. Use original R2 metadata or PLINK metadata exported by ldscore --export-ref-metadata.'
        ),
    )
    query_group = original.add_mutually_exclusive_group()
    query_group.add_argument(
        '--query-annot-sources', nargs='+', metavar='SOURCES',
        help=(
            "Original query annotation tables for the selected per-query fit; accepts paths, '*' and '@' patterns "
            '(see --target-annot-sources). Cannot be combined with --query-annot-bed-sources or '
            '--query-annot-gene-list-sources. Omit for baseline-only fits; per-query fits require their matching '
            'source.'
        ),
    )
    query_group.add_argument(
        '--query-annot-bed-sources', nargs='+', metavar='SOURCES',
        help=(
            "Original BED query files for the selected per-query fit; accepts paths or '*' patterns. Cannot be "
            'combined with --query-annot-sources or --query-annot-gene-list-sources. Omit for other query types; '
            "reproduce the original --padding-bp. '@' is not expanded; see --target-annot-sources for wildcard "
            'and quoting rules.'
        ),
    )
    query_group.add_argument(
        '--query-annot-gene-list-sources', nargs='+', metavar='SOURCES',
        help=(
            'Original one-column gene-list files for the selected per-query fit. Requires --gene-coordinate-file; '
            'cannot be combined with --query-annot-sources or --query-annot-bed-sources. Omit for other query '
            "types; reproduce the original gene settings. Accepts '*' patterns (see --target-annot-sources); '@' "
            'is not expanded.'
        ),
    )

    quantiles.add_argument(
        '--target-missing-value', default=None,
        help=(
            'One target-value token to exclude before forming quantiles, for example NaN or -999. Default: '
            'exclude no value; unselected nonnumeric or nonfinite values cause an error. Zero remains a '
            'valid score unless explicitly selected.'
        ),
    )
    quantiles.add_argument(
        '--num-quantiles', type=int, default=5, metavar='N',
        help=(
            'Number of low-to-high target quantiles; must be at least 2. Default: 5. Ties stay in the lower '
            'quantile; an empty quantile causes an error.'
        ),
    )

    genes.add_argument(
        '--gene-coordinate-file', default=None, metavar='FILE',
        help=(
            'Original one-based inclusive gene-coordinate TSV/TSV.GZ catalog. Requires '
            '--query-annot-gene-list-sources and is required by that flag; no default catalog. Supply an exact '
            'file path; patterns are not expanded.'
        ),
    )
    genes.add_argument(
        '--control-gene-list-file', default=None, metavar='FILE',
        help=(
            'Original fixed-control gene list, when the selected fit used one. Requires '
            '--query-annot-gene-list-sources. If omitted, add no gene control. Supply an exact file path; '
            'patterns are not expanded.'
        ),
    )
    genes.add_argument(
        '--padding-bp', type=int, default=0, metavar='BP',
        help=(
            'Extend BED intervals or gene boundaries by this many base pairs at both ends before finding covered '
            'SNPs. Default: 0. Use the original LD-score setting; not applicable to prebuilt annotation tables. '
            'Use 0 for BED files that already include the padding used in the fitted model; avoid double padding.'
        ),
    )
    genes.add_argument(
        '--gene-list-resolution-policy', choices=('strict', 'resolved-only'), default='strict',
        help=(
            (
            'Handle rejected gene identifiers: strict stops; resolved-only continues with usable genes and '
            'requires --query-annot-gene-list-sources. Default: strict. Reproduce the original LD-score '
            'policy.'
        )
        ),
    )
    genes.add_argument(
        '--gene-exclude-regions', choices=('none', 'mhc'), default='none',
        help=(
            (
            'Exclude MHC-overlapping genes before padding, or none to exclude no genes. Default: none. '
            'Selecting mhc requires --query-annot-gene-list-sources; reproduce the original LD-score '
            'setting.'
        )
        ),
    )

    runtime.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )
    runtime.add_argument(
        '--log-level', choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'), default='INFO',
        help=LOG_LEVEL_HELP,
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the standalone parser used by the unified CLI."""
    parser = argparse.ArgumentParser(
        prog="ldsc quantile-h2",
        description="Project a fitted partitioned-LDSC model onto continuous-annotation quantiles.",
        allow_abbrev=False,
    )
    add_quantile_h2_arguments(parser)
    return parser


@materializing_overwrite_guard(
    lambda args: (
        (getattr(args, "output_dir"), getattr(args, "overwrite", False), "RUN_FAILED.txt")
        if getattr(args, "output_dir", None)
        else None
    ),
    command="run_quantile_h2_from_args(...)",
)
def run_quantile_h2_from_args(args) -> QuantileH2Result:
    """Run the verified post-fit continuous-annotation quantile workflow.

    The workflow reconstructs the fitted model's common reference-SNP
    universe from the original annotation sources and reference metadata,
    assigns eligible SNPs to target-value quantiles, and projects the saved
    whole-data and delete-one-block coefficient vectors onto those quantiles.
    Alignment, counts, and available overlap cross-products are checked; these
    checks do not establish the original annotation value at every SNP.
    It does not refit LDSC. A successful overwrite removes the default plot
    root derived from the superseded quantile result.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command arguments registered by
        :func:`add_quantile_h2_arguments`.

    Returns
    -------
    QuantileH2Result
        Quantile summaries, standardized coefficients, alignment diagnostics,
        provenance metadata, and written output paths.

    Raises
    ------
    LDSCInputError
        If the fitted model is unsupported, required sources cannot be aligned,
        provenance checks fail, target values are invalid, or a requested
        quantile is empty.
    """
    output_dir = ensure_output_directory(args.output_dir, label="output directory")
    diagnostics_dir = output_dir / "diagnostics"
    writer = QuantileH2DirectoryWriter()
    log_path = diagnostics_dir / "quantile-h2.log"
    writer.artifact_family(output_dir).preflight(overwrite=args.overwrite, additional_paths=[log_path])
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    issues_path = diagnostics_dir / "snp_alignment_issues.tsv.gz"
    with workflow_logging("quantile-h2", log_path, log_level=args.log_level):
        from ._input_preflight import inspect_declared_inputs, inspect_artifact_paths
        inspect_declared_inputs(files=[(name, getattr(args, name, None)) for name in (
            'baseline_annot_sources', 'query_annot_sources', 'query_annot_bed_sources', 'ref_metadata_sources')],
            checks=[('fitted model', args.partitioned_h2_result_dir,
                     lambda: inspect_artifact_paths(args.partitioned_h2_result_dir, metadata_file=('metadata.json' if (Path(args.partitioned_h2_result_dir)/'metadata.json').is_file() else 'diagnostics/metadata.json')))],
            issues_path=diagnostics_dir/'input_issues.tsv')
        log_inputs(
            partitioned_h2_result_dir=args.partitioned_h2_result_dir,
            target_annotation=args.target_annotation,
            output_dir=str(output_dir),
        )
        with AnnotationWorkspace(output_dir) as workspace, AlignmentDiagnostics(issues_path) as diagnostics:
            model = load_fitted_partitioned_model(args.partitioned_h2_result_dir)
            inputs = prepare_quantile_statistics(args, model, workspace, diagnostics)
        ldscore_metadata = inputs.ldscore_metadata
        common_maf_min = inputs.common_maf_min
        samp_prev = pd.to_numeric(model.summary.get("samp_prev", pd.Series([np.nan])), errors="coerce").iloc[0]
        pop_prev = pd.to_numeric(model.summary.get("pop_prev", pd.Series([np.nan])), errors="coerce").iloc[0]
        liability_factor = None
        if np.isfinite(samp_prev) and np.isfinite(pop_prev):
            from ._kernel.regression import liability_conversion_factor
            liability_factor = float(liability_conversion_factor(float(samp_prev), float(pop_prev)))
        quantile_result = compute_quantile_h2(
            annotation_sums=inputs.annotation_sums,
            tau=model.tau,
            tau_delete=model.tau_delete,
            snp_counts=inputs.counts,
            lower=inputs.lower,
            upper=inputs.upper,
            liability_factor=liability_factor,
        )
        eligible_total_h2 = float(quantile_result["h2_obs"].sum())
        if eligible_total_h2 < 0:
            LOGGER.warning(
                "Target-eligible fitted-model heritability is negative (%g). Negative quantile estimates and "
                "nonzero-total proportion/enrichment ratios remain visible; interpret them as noisy linear-model "
                "estimates, not literal negative biological variance.",
                eligible_total_h2,
            )
        elif eligible_total_h2 == 0:
            LOGGER.warning(
                "Target-eligible fitted-model heritability is zero. Quantile h2 estimates remain visible, while "
                "prop_h2 and enrichment ratio fields are NaN because their denominator is zero."
            )
        total_h2 = float(model.tau @ inputs.full_sums)
        annotation_types = {
            name: str((ldscore_metadata.get("annotation_types") or {}).get(name, "unknown"))
            for name in model.annotation_names
        }
        standardized = compute_standardized_coefficients(
            annotation_names=model.annotation_names,
            annotation_types=annotation_types,
            annotation_sd=inputs.annotation_sd,
            tau=model.tau,
            tau_se=model.tau_se,
            total_h2=total_h2,
            n_common_snps=inputs.n_common,
        )
        metadata = {
            "selected_model_path": str(model.result_dir),
            "selected_model_type": model.model_type,
            "ldscore_dir": str(model.ldscore_dir),
            "snp_identifier": ldscore_metadata.get("snp_identifier"),
            "genome_build": ldscore_metadata.get("genome_build"),
            "retained_ld_columns": model.annotation_names,
            "target_annotation": args.target_annotation,
            "target_annot_sources": inputs.target_paths,
            "ref_metadata_sources": inputs.reference_paths,
            "common_reference_snp_maf_min": common_maf_min,
            "common_reference_snp_maf_operator": ">=",
            "common_reference_snp_universe_size": int(inputs.n_common),
            "target_eligible_snp_count": int(inputs.counts.sum()),
            "target_missing_exclusion_count": inputs.n_excluded,
            "target_missing_value": args.target_missing_value,
            "num_quantiles": int(args.num_quantiles),
            "quantile_boundary_rule": "round_i_times_n_minus_1_over_q; lower_ties",
            "enrichment_p_test": "two_sided_normal_inside_vs_complement",
            "tau_star_scale": "fixed_full_common_reference_snp_universe_total",
            "samp_prev": None if not np.isfinite(samp_prev) else float(samp_prev),
            "pop_prev": None if not np.isfinite(pop_prev) else float(pop_prev),
        }
        LOGGER.info(
            "Fitted model=%s; target=%s; common MAF rule=MAF >= %g; common reference-SNP universe=%d; "
            "target missing exclusions=%d.",
            model.model_type,
            args.target_annotation,
            common_maf_min,
            inputs.n_common,
            inputs.n_excluded,
        )
        for row in quantile_result.itertuples(index=False):
            LOGGER.info(
                "Quantile %d: target [%g, %g], n_snps=%d.",
                row.quantile,
                row.target_value_lower,
                row.target_value_upper,
                row.n_snps,
            )
        if total_h2 <= 0:
            LOGGER.warning(
                "Full fitted-model heritability over the common reference-SNP universe is nonpositive (%g); "
                "tau_star fields are NaN while raw tau fields remain available.",
                total_h2,
            )
        output_paths = writer.write(
            quantile_result,
            standardized,
            issues_path,
            QuantileH2OutputConfig(output_dir=output_dir, overwrite=True),
            metadata,
        )
        log_outputs(**output_paths)
    return QuantileH2Result(quantile_result, standardized, str(issues_path), metadata, output_paths)


def main(argv: Sequence[str] | None = None) -> QuantileH2Result:
    """Parse standalone arguments and execute ``quantile-h2``."""
    return run_quantile_h2_from_args(build_parser().parse_args(argv))


__all__ = [
    "QUANTILE_H2_COLUMNS",
    "STANDARDIZED_COEFFICIENT_COLUMNS",
    "QuantileAssignment",
    "assign_legacy_quantiles",
    "compute_quantile_h2",
    "compute_standardized_coefficients",
    "add_quantile_h2_arguments",
    "build_parser",
    "load_fitted_partitioned_model",
    "run_quantile_h2_from_args",
]
