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

from .errors import LDSCInputError
from .annotation_builder import AnnotationBuilder
from .annotation_semantics import (
    FINGERPRINT_CANONICALIZATION,
    build_annotation_fingerprint_metadata,
)
from .column_inference import (
    A1_COLUMN_SPEC,
    A2_COLUMN_SPEC,
    CHR_COLUMN_SPEC,
    POS_COLUMN_SPEC,
    SNP_COLUMN_SPEC,
    resolve_optional_column,
    resolve_required_column,
)
from .config import AnnotationBuildConfig, GlobalConfig
from ._kernel import annotation as kernel_annotation
from ._kernel.snp_identity import effective_merge_key_series, identity_base_mode, is_allele_aware_mode
from .outputs import (
    QuantileH2DirectoryWriter,
    QuantileH2OutputConfig,
    SNP_ALIGNMENT_ISSUE_COLUMNS,
)
from .path_resolution import (
    ANNOTATION_SUFFIXES,
    ensure_output_directory,
    preflight_output_artifact_family,
    resolve_file_group,
    split_cli_path_tokens,
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
        Quantile summaries in :data:`QUANTILE_H2_COLUMNS` order.
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
    enrichment_p = np.where(
        contrast_se > 0,
        2.0 * stats.norm.sf(np.abs(contrast / contrast_se)),
        np.nan,
    )

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
    """In-memory result family returned by :func:`run_quantile_h2_from_args`.

    Attributes
    ----------
    quantile_h2 : pandas.DataFrame
        Low-to-high target quantile summaries.
    standardized_coefficients : pandas.DataFrame
        Raw and standardized coefficients in fitted annotation order.
    alignment_issues : pandas.DataFrame
        Nonfatal SNP exclusions recorded during source reconstruction.
    metadata : dict
        Model, input, universe, missingness, and statistic provenance.
    output_paths : dict of str to str
        Written result-family paths.
    """

    quantile_h2: pd.DataFrame
    standardized_coefficients: pd.DataFrame
    alignment_issues: pd.DataFrame
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


def _empty_alignment_issues() -> pd.DataFrame:
    """Return an empty table with the stable SNP-alignment schema."""
    return pd.DataFrame(columns=SNP_ALIGNMENT_ISSUE_COLUMNS)


def _read_reference_metadata(sources: Sequence[str]) -> tuple[pd.DataFrame, list[str]]:
    """Read canonical whole-genome or chromosome-sharded reference metadata."""
    paths = resolve_file_group(
        sources,
        suffixes=("", ".tsv", ".tsv.gz", ".txt", ".txt.gz"),
        label="reference metadata",
        allow_chromosome_suite=True,
    )
    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = kernel_annotation._read_text_table(path)
        aliases = {"BP": "POS", "FRQ": "MAF"}
        frame = frame.rename(columns={old: new for old, new in aliases.items() if old in frame.columns and new not in frame.columns})
        required = [column for column in ("CHR", "POS", "SNP", "MAF") if column not in frame.columns]
        if required:
            raise LDSCInputError(
                f"Reference metadata '{path}' is missing required column(s) {required}. "
                "Expected CHR, POS (or BP), SNP, and MAF."
            )
        keep = ["CHR", "POS", "SNP", "MAF", *[column for column in ("A1", "A2") if column in frame.columns]]
        frames.append(frame.loc[:, keep].copy())
    combined = pd.concat(frames, ignore_index=True)
    combined["CHR"] = combined["CHR"].astype(str).str.replace(r"^chr", "", regex=True)
    combined["POS"] = pd.to_numeric(combined["POS"], errors="raise").astype(np.int64)
    combined["MAF"] = pd.to_numeric(combined["MAF"], errors="coerce")
    return combined, paths


def _read_target_annotation(
    sources: Sequence[str],
    target_annotation: str,
    missing_token: str | None,
) -> tuple[pd.DataFrame, pd.Series, pd.Series, list[str]]:
    """Read target values while preserving raw tokens for explicit exclusion."""
    paths = resolve_file_group(
        sources,
        suffixes=ANNOTATION_SUFFIXES,
        label="target annotation",
        allow_chromosome_suite=True,
    )
    metadata_frames: list[pd.DataFrame] = []
    raw_values: list[pd.Series] = []
    for path in paths:
        try:
            frame = pd.read_csv(path, sep=r"\s+", dtype=str, keep_default_na=False)
        except (pd.errors.ParserError, UnicodeDecodeError, ValueError) as exc:
            raise LDSCInputError(f"quantile-h2 could not parse target annotation '{path}': {exc}") from exc
        chr_column = resolve_required_column(frame.columns, CHR_COLUMN_SPEC, context=path)
        pos_column = resolve_required_column(frame.columns, POS_COLUMN_SPEC, context=path)
        snp_column = resolve_required_column(frame.columns, SNP_COLUMN_SPEC, context=path)
        if target_annotation not in frame.columns:
            raise LDSCInputError(
                f"Target annotation '{target_annotation}' is absent from '{path}'. Choose an existing column "
                "with --target-annotation or supply the correct --target-annot-sources."
            )
        a1_column = resolve_optional_column(frame.columns, A1_COLUMN_SPEC, context=path)
        a2_column = resolve_optional_column(frame.columns, A2_COLUMN_SPEC, context=path)
        metadata = pd.DataFrame(
            {"CHR": frame[chr_column], "POS": frame[pos_column], "SNP": frame[snp_column]}
        )
        if a1_column is not None and a2_column is not None:
            metadata["A1"] = frame[a1_column]
            metadata["A2"] = frame[a2_column]
        metadata_frames.append(metadata)
        raw_values.append(frame[target_annotation].astype("string"))
    metadata = pd.concat(metadata_frames, ignore_index=True)
    raw = pd.concat(raw_values, ignore_index=True)
    excluded = pd.Series(False, index=raw.index)
    if missing_token is not None:
        try:
            numeric_token = float(missing_token)
        except ValueError:
            excluded = raw.str.strip().eq(str(missing_token).strip()).fillna(False)
        else:
            if not np.isfinite(numeric_token):
                excluded = raw.str.strip().eq(str(missing_token).strip()).fillna(False)
            else:
                numeric_raw = pd.to_numeric(raw, errors="coerce")
                excluded = numeric_raw.eq(numeric_token).fillna(False)
    numeric = pd.to_numeric(raw.mask(excluded), errors="coerce")
    invalid = ~excluded & (~np.isfinite(numeric.to_numpy(dtype=float, na_value=np.nan)))
    if bool(invalid.any()):
        examples = sorted(set(raw.loc[invalid].astype(str).head(5)))
        raise LDSCInputError(
            f"Target annotation '{target_annotation}' contains nonnumeric, NaN, or infinite value(s): {examples}. "
            "If one token denotes missingness, pass it with --target-missing-value; otherwise clean the source."
        )
    metadata["CHR"] = metadata["CHR"].astype(str).str.replace(r"^chr", "", regex=True)
    metadata["POS"] = pd.to_numeric(metadata["POS"], errors="raise").astype(np.int64)
    return metadata, numeric.astype(np.float64), excluded.astype(bool), paths


def _effective_keys_with_reference_inference(
    metadata: pd.DataFrame,
    reference: pd.DataFrame,
    mode: str,
    *,
    role: str,
) -> pd.Series:
    """Build effective keys, inferring omitted annotation alleles when unambiguous.

    Rows absent from reference metadata receive deterministic nonmatching keys
    so the caller can report and exclude them. Nonunique reference base
    identities are fatal because allele inference would be ambiguous.
    """
    if not is_allele_aware_mode(mode) or {"A1", "A2"}.issubset(metadata.columns):
        return effective_merge_key_series(metadata, mode, context=role)
    base_mode = identity_base_mode(mode)
    reference_base = effective_merge_key_series(reference, base_mode, context="reference metadata base identity")
    duplicate_base = reference_base.duplicated(keep=False)
    if duplicate_base.any():
        issue_frame = reference.copy()
        issue_frame["effective_snp_id"] = reference_base.astype(str)
        raise AlignmentError(
            f"quantile-h2 cannot infer omitted alleles for {role}: reference metadata has non-unique "
            f"{base_mode} identities. Supply A1/A2 in the annotation source.",
            [
                _alignment_issue_rows(
                    issue_frame,
                    duplicate_base,
                    role="reference_metadata",
                    source="resupplied reference metadata",
                    issue="duplicate_identity",
                    action="fatal",
                    details=f"non-unique {base_mode} identity prevents allele inference for {role}",
                )
            ],
        )
    mapping = pd.Series(
        effective_merge_key_series(reference, mode, context="reference metadata").to_numpy(),
        index=reference_base.astype(str),
    )
    base = effective_merge_key_series(metadata, base_mode, context=role).astype(str)
    keys = base.map(mapping)
    if keys.isna().any():
        missing_positions = np.flatnonzero(keys.isna().to_numpy())
        for position in missing_positions:
            keys.iloc[position] = f"__unresolved_reference_metadata__:{role}:{position}:{base.iloc[position]}"
    return keys.reset_index(drop=True).rename("effective_snp_id")


def _alignment_issue_rows(
    frame: pd.DataFrame,
    mask: pd.Series | np.ndarray,
    *,
    role: str,
    source: str,
    issue: str,
    action: str,
    details: str,
) -> pd.DataFrame:
    """Build diagnostic rows for one alignment condition."""
    subset = frame.loc[np.asarray(mask, dtype=bool)].copy()
    if subset.empty:
        return _empty_alignment_issues()
    rows = pd.DataFrame(
        {
            "source_role": role,
            "source": source,
            "CHR": subset.get("CHR", ""),
            "POS": subset.get("POS", ""),
            "SNP": subset.get("SNP", ""),
            "A1": subset.get("A1", ""),
            "A2": subset.get("A2", ""),
            "effective_snp_id": subset.get("effective_snp_id", ""),
            "issue": issue,
            "action": action,
            "details": details,
        }
    )
    return rows.loc[:, SNP_ALIGNMENT_ISSUE_COLUMNS]


def _prepare_quantile_inputs(args, model: FittedPartitionedModel):
    """Reconstruct and verify fitted/target values over the inherited universe."""
    ldscore_metadata = json.loads((model.ldscore_dir / "metadata.json").read_text(encoding="utf-8"))
    mode = str(ldscore_metadata.get("snp_identifier"))
    genome_build = ldscore_metadata.get("genome_build")
    global_config = GlobalConfig(snp_identifier=mode, genome_build=genome_build)
    annotation_config = AnnotationBuildConfig(
        baseline_annot_sources=tuple(split_cli_path_tokens(args.baseline_annot_sources)),
        query_annot_sources=tuple(split_cli_path_tokens(getattr(args, "query_annot_sources", None))),
        query_annot_bed_sources=tuple(split_cli_path_tokens(getattr(args, "query_annot_bed_sources", None))),
        query_annot_gene_list_sources=tuple(split_cli_path_tokens(getattr(args, "query_annot_gene_list_sources", None))),
        gene_coordinate_file=getattr(args, "gene_coordinate_file", None),
        control_gene_list_file=getattr(args, "control_gene_list_file", None),
        gene_list_resolution_policy=getattr(args, "gene_list_resolution_policy", "strict"),
        gene_exclude_regions=getattr(args, "gene_exclude_regions", "none"),
        padding_bp=getattr(args, "padding_bp", 0),
    )
    bundle = AnnotationBuilder(
        global_config,
        projection_genome_build=genome_build,
    ).run(annotation_config)
    all_annotations = bundle.annotation_matrix(include_query=True)
    missing_annotations = [name for name in model.annotation_names if name not in all_annotations.columns]
    if missing_annotations:
        raise LDSCInputError(
            f"quantile-h2 resupplied annotation sources are missing fitted annotation(s): {missing_annotations}. "
            "Resupply every original baseline/query source for this fitted model."
        )
    fitted_values = all_annotations.loc[:, model.annotation_names].astype(np.float64)
    reference, reference_paths = _read_reference_metadata(split_cli_path_tokens(args.ref_metadata_sources))
    issues: list[pd.DataFrame] = []
    if reference["MAF"].isna().any() or (~np.isfinite(reference["MAF"])).any():
        bad = reference["MAF"].isna() | ~np.isfinite(reference["MAF"])
        issues.append(_alignment_issue_rows(
            reference.assign(effective_snp_id=""), bad, role="reference_metadata", source=",".join(reference_paths),
            issue="missing_maf", action="fatal", details="MAF must be finite for every intersected reference SNP",
        ))
        raise AlignmentError(
            "quantile-h2 reference metadata contains missing or nonfinite MAF values.",
            issues,
        )
    reference = reference.reset_index(drop=True)
    reference["effective_snp_id"] = effective_merge_key_series(reference, mode, context="reference metadata").astype(str)
    duplicate_ref = reference["effective_snp_id"].duplicated(keep=False)
    if duplicate_ref.any():
        issues.append(_alignment_issue_rows(
            reference, duplicate_ref, role="reference_metadata", source=",".join(reference_paths),
            issue="duplicate_identity", action="fatal", details="effective SNP identity occurs more than once",
        ))
        raise AlignmentError("quantile-h2 reference metadata has duplicate effective SNP identities.", issues)

    fitted_metadata = bundle.metadata.reset_index(drop=True).copy()
    fitted_metadata["effective_snp_id"] = _effective_keys_with_reference_inference(
        fitted_metadata, reference, mode, role="fitted annotations"
    ).astype(str)
    fitted = pd.concat([fitted_metadata, fitted_values.reset_index(drop=True)], axis=1)
    reference_ids = set(reference["effective_snp_id"])
    fitted_ids = set(fitted["effective_snp_id"])
    annotation_only = ~fitted["effective_snp_id"].isin(reference_ids)
    reference_only = ~reference["effective_snp_id"].isin(fitted_ids)
    if annotation_only.any():
        issues.append(_alignment_issue_rows(
            fitted, annotation_only, role="fitted_annotation", source="resupplied fitted annotations",
            issue="missing_reference_metadata", action="excluded", details="annotation SNP is absent from reference metadata",
        ))
    if reference_only.any():
        issues.append(_alignment_issue_rows(
            reference, reference_only, role="reference_metadata", source=",".join(reference_paths),
            issue="missing_baseline_annotation", action="excluded", details="reference SNP is outside the fitted annotation grid",
        ))
    common_maf_min = float((ldscore_metadata.get("count_config") or {}).get("common_reference_snp_maf_min", 0.05))
    common_operator = str((ldscore_metadata.get("count_config") or {}).get("common_reference_snp_maf_operator", ">="))
    if common_operator != ">=":
        raise LDSCInputError(f"quantile-h2 does not support inherited common-MAF operator {common_operator!r}.")
    aligned_reference = reference.merge(
        fitted.loc[:, ["effective_snp_id", *model.annotation_names]],
        on="effective_snp_id",
        how="inner",
        sort=False,
        validate="one_to_one",
    )
    expected_all = (ldscore_metadata.get("overlap_config") or {}).get("total_all_reference_snps")
    if expected_all is not None and int(round(float(expected_all))) != len(aligned_reference):
        raise LDSCInputError(
            "quantile-h2 reconstructed a different all-reference SNP universe size than the fitted LD-score artifact: "
            f"{len(aligned_reference)} vs {expected_all}. Check --ref-metadata-sources and annotation sources."
        )
    common = aligned_reference.loc[aligned_reference["MAF"] >= common_maf_min].reset_index(drop=True)

    target_metadata, target_values, target_excluded, target_paths = _read_target_annotation(
        split_cli_path_tokens(args.target_annot_sources),
        args.target_annotation,
        getattr(args, "target_missing_value", None),
    )
    target_metadata["effective_snp_id"] = _effective_keys_with_reference_inference(
        target_metadata, reference, mode, role="target annotation"
    ).astype(str)
    duplicate_target = target_metadata["effective_snp_id"].duplicated(keep=False)
    if duplicate_target.any():
        issues.append(_alignment_issue_rows(
            target_metadata, duplicate_target, role="target_annotation", source=",".join(target_paths),
            issue="duplicate_identity", action="fatal", details="effective SNP identity occurs more than once",
        ))
        raise AlignmentError("quantile-h2 target annotation has duplicate effective SNP identities.", issues)
    target_only = ~target_metadata["effective_snp_id"].isin(reference_ids)
    if target_only.any():
        issues.append(_alignment_issue_rows(
            target_metadata, target_only, role="target_annotation", source=",".join(target_paths),
            issue="missing_reference_metadata", action="excluded", details="target SNP is absent from reference metadata",
        ))
    target = target_metadata.loc[:, ["effective_snp_id"]].copy()
    target["target_value"] = target_values
    target["target_excluded"] = target_excluded.to_numpy(dtype=bool)
    common = common.merge(target, on="effective_snp_id", how="left", sort=False, validate="one_to_one")
    missing_target = common["target_value"].isna() & ~common["target_excluded"].fillna(False)
    absent_target = common["target_excluded"].isna()
    if missing_target.any() or absent_target.any():
        mask = missing_target | absent_target
        issues.append(_alignment_issue_rows(
            common, mask, role="target_annotation", source=",".join(target_paths),
            issue="missing_target_annotation", action="fatal", details="common reference-SNP universe member is absent from target source",
        ))
        raise AlignmentError("quantile-h2 target annotation does not cover the common reference-SNP universe.", issues)

    expected_common = (ldscore_metadata.get("overlap_config") or {}).get("total_common_reference_snps")
    if expected_common is not None and int(round(float(expected_common))) != len(common):
        raise LDSCInputError(
            "quantile-h2 reconstructed a different common reference-SNP universe size than the fitted LD-score artifact: "
            f"{len(common)} vs {expected_common}. Check --ref-metadata-sources and the original annotation sources."
        )
    count_records = {str(record.get("column")): record for record in ldscore_metadata.get("counts", [])}
    for name in model.annotation_names:
        expected_sum = count_records.get(name, {}).get("common_reference_snp_count")
        actual_sum = float(np.sum(common[name].to_numpy(dtype=np.float64), dtype=np.float64))
        if expected_sum is not None and not np.isclose(actual_sum, float(expected_sum), rtol=1e-6, atol=1e-8):
            raise LDSCInputError(
                f"quantile-h2 fitted annotation '{name}' has common reference-SNP universe sum {actual_sum}, expected {expected_sum}. "
                "Resupply the exact annotation sources used for LD-score construction."
            )
    overlap_rel = (ldscore_metadata.get("files") or {}).get("overlap")
    if overlap_rel:
        from .overlap_matrix import assemble_model_overlap, overlap_from_long_frame

        overlap_config = ldscore_metadata.get("overlap_config") or {}
        overlap = overlap_from_long_frame(
            pd.read_parquet(model.ldscore_dir / str(overlap_rel)),
            baseline_columns=[str(name) for name in ldscore_metadata.get("baseline_columns", [])],
            query_columns=[str(name) for name in ldscore_metadata.get("query_columns", [])],
            total_all_reference_snps=overlap_config.get("total_all_reference_snps"),
            total_common_reference_snps=overlap_config.get("total_common_reference_snps"),
        )
        expected_crossproducts = assemble_model_overlap(overlap, model.annotation_names, use_common=True)
        actual_crossproducts = common.loc[:, model.annotation_names].to_numpy(dtype=np.float64).T @ common.loc[:, model.annotation_names].to_numpy(dtype=np.float64)
        if not np.allclose(actual_crossproducts, expected_crossproducts, rtol=1e-6, atol=1e-8):
            raise LDSCInputError(
                "quantile-h2 reconstructed fitted-annotation cross-products that disagree with the linked "
                "LD-score overlap artifact. Resupply the exact original annotations and reference metadata."
            )
    current_fingerprints = build_annotation_fingerprint_metadata(
        common["effective_snp_id"], common.loc[:, model.annotation_names]
    )
    stored_fingerprints = ldscore_metadata.get("annotation_fingerprints")
    verification_level = "aggregate_only"
    if stored_fingerprints:
        verification_level = "exact_common_values"
        fingerprints_match = (
            current_fingerprints["algorithm"] == stored_fingerprints.get("algorithm")
            and current_fingerprints["canonicalization"] == stored_fingerprints.get("canonicalization")
            and current_fingerprints["common_reference_snp_universe"]
            == stored_fingerprints.get("common_reference_snp_universe")
            and all(
                current_fingerprints["annotation_values"][name]
                == (stored_fingerprints.get("annotation_values") or {}).get(name)
                for name in model.annotation_names
            )
        )
        if not fingerprints_match:
            raise LDSCInputError(
                "quantile-h2 SHA256 semantic fingerprint mismatch for the common reference-SNP universe or fitted "
                "annotation values. Resupply the exact inputs used to build the linked LD-score artifact."
            )
    else:
        LOGGER.warning(
            "Linked LD-score artifact has no annotation fingerprints; using aggregate-only verification. "
            "Regenerate LD scores with the current package for exact common-value verification."
        )
    if args.target_annotation in model.annotation_names:
        fitted_target = common[args.target_annotation].to_numpy(dtype=np.float64)
        external_target = common["target_value"].to_numpy(dtype=np.float64)
        comparable = ~common["target_excluded"].to_numpy(dtype=bool)
        if not np.array_equal(fitted_target[comparable].astype(np.float32), external_target[comparable].astype(np.float32)):
            raise LDSCInputError(
                f"Target annotation '{args.target_annotation}' matches a fitted annotation name but its values differ. "
                "Supply the same values or rename the external target annotation."
            )
    return common, issues, ldscore_metadata, reference_paths, target_paths, verification_level, common_maf_min


class AlignmentError(LDSCInputError):
    """Fatal alignment error carrying machine-readable issue tables."""

    def __init__(self, message: str, issue_frames: Sequence[pd.DataFrame]) -> None:
        super().__init__(message)
        self.issues = pd.concat(issue_frames, ignore_index=True) if issue_frames else _empty_alignment_issues()


def add_quantile_h2_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the public ``quantile-h2`` command-line contract."""
    parser.allow_abbrev = False
    parser.add_argument("--partitioned-h2-result-dir", required=True)
    parser.add_argument("--baseline-annot-sources", nargs="+", required=True)
    query_group = parser.add_mutually_exclusive_group()
    query_group.add_argument("--query-annot-sources", nargs="+")
    query_group.add_argument("--query-annot-bed-sources", nargs="+")
    query_group.add_argument("--query-annot-gene-list-sources", nargs="+")
    parser.add_argument("--gene-coordinate-file", default=None)
    parser.add_argument("--control-gene-list-file", default=None)
    parser.add_argument("--gene-list-resolution-policy", choices=("strict", "resolved-only"), default="strict")
    parser.add_argument("--gene-exclude-regions", choices=("none", "mhc"), default="none")
    parser.add_argument("--padding-bp", type=int, default=0)
    parser.add_argument("--target-annot-sources", nargs="+", required=True)
    parser.add_argument("--target-annotation", required=True)
    parser.add_argument("--ref-metadata-sources", nargs="+", required=True)
    parser.add_argument("--target-missing-value", default=None)
    parser.add_argument("--num-quantiles", type=int, default=5)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--overwrite", action="store_true", default=False)
    parser.add_argument("--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"), default="INFO")


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
    paths = [
        output_dir / "quantile_h2.tsv",
        output_dir / "standardized_coefficients.tsv",
        diagnostics_dir / "metadata.json",
        diagnostics_dir / "snp_alignment_issues.tsv.gz",
        diagnostics_dir / "quantile-h2.log",
    ]
    preflight_output_artifact_family(
        paths,
        [*paths, output_dir / "plots"],
        overwrite=args.overwrite,
        label="quantile-h2 output artifact",
    )
    diagnostics_dir.mkdir(parents=True, exist_ok=True)
    writer = QuantileH2DirectoryWriter()
    empty_issues = _empty_alignment_issues()
    writer.write_diagnostics(empty_issues, QuantileH2OutputConfig(output_dir=output_dir, overwrite=True))
    log_path = diagnostics_dir / "quantile-h2.log"
    with workflow_logging("quantile-h2", log_path, log_level=args.log_level):
        log_inputs(
            partitioned_h2_result_dir=args.partitioned_h2_result_dir,
            target_annotation=args.target_annotation,
            output_dir=str(output_dir),
        )
        model = load_fitted_partitioned_model(args.partitioned_h2_result_dir)
        try:
            common, issue_frames, ldscore_metadata, reference_paths, target_paths, verification_level, common_maf_min = (
                _prepare_quantile_inputs(args, model)
            )
        except AlignmentError as exc:
            writer.write_diagnostics(exc.issues, QuantileH2OutputConfig(output_dir=output_dir, overwrite=True))
            raise
        issues = pd.concat(issue_frames, ignore_index=True) if issue_frames else empty_issues
        eligible = ~common["target_excluded"].to_numpy(dtype=bool)
        target_values = common.loc[eligible, "target_value"].to_numpy(dtype=np.float64)
        assignment = assign_legacy_quantiles(target_values, args.num_quantiles)
        fitted_values = common.loc[eligible, model.annotation_names].to_numpy(dtype=np.float64)
        indicator = np.equal.outer(assignment.quantile, np.arange(1, args.num_quantiles + 1)).astype(np.float64)
        annotation_sums = fitted_values.T @ indicator
        samp_prev = pd.to_numeric(model.summary.get("samp_prev", pd.Series([np.nan])), errors="coerce").iloc[0]
        pop_prev = pd.to_numeric(model.summary.get("pop_prev", pd.Series([np.nan])), errors="coerce").iloc[0]
        liability_factor = None
        if np.isfinite(samp_prev) and np.isfinite(pop_prev):
            from ._kernel.regression import liability_conversion_factor
            liability_factor = float(liability_conversion_factor(float(samp_prev), float(pop_prev)))
        quantile_result = compute_quantile_h2(
            annotation_sums=annotation_sums,
            tau=model.tau,
            tau_delete=model.tau_delete,
            snp_counts=assignment.counts,
            lower=assignment.lower,
            upper=assignment.upper,
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
        full_values = common.loc[:, model.annotation_names].to_numpy(dtype=np.float64)
        full_sums = full_values.sum(axis=0, dtype=np.float64)
        total_h2 = float(model.tau @ full_sums)
        annotation_types = {
            name: str((ldscore_metadata.get("annotation_types") or {}).get(name, "unknown"))
            for name in model.annotation_names
        }
        standardized = compute_standardized_coefficients(
            annotation_names=model.annotation_names,
            annotation_types=annotation_types,
            annotation_sd=full_values.std(axis=0, ddof=0),
            tau=model.tau,
            tau_se=model.tau_se,
            total_h2=total_h2,
            n_common_snps=len(common),
        )
        metadata = {
            "selected_model_path": str(model.result_dir),
            "selected_model_type": model.model_type,
            "ldscore_dir": str(model.ldscore_dir),
            "snp_identifier": ldscore_metadata.get("snp_identifier"),
            "genome_build": ldscore_metadata.get("genome_build"),
            "retained_ld_columns": model.annotation_names,
            "target_annotation": args.target_annotation,
            "target_annot_sources": target_paths,
            "ref_metadata_sources": reference_paths,
            "common_reference_snp_maf_min": common_maf_min,
            "common_reference_snp_maf_operator": ">=",
            "common_reference_snp_universe_size": int(len(common)),
            "target_eligible_snp_count": int(eligible.sum()),
            "target_missing_exclusion_count": int((~eligible).sum()),
            "target_missing_value": args.target_missing_value,
            "num_quantiles": int(args.num_quantiles),
            "quantile_boundary_rule": "round_i_times_n_minus_1_over_q; lower_ties",
            "verification_level": verification_level,
            "fingerprint_canonicalization": FINGERPRINT_CANONICALIZATION,
            "enrichment_p_test": "two_sided_normal_inside_vs_complement",
            "tau_star_scale": "fixed_full_common_reference_snp_universe_total",
            "samp_prev": None if not np.isfinite(samp_prev) else float(samp_prev),
            "pop_prev": None if not np.isfinite(pop_prev) else float(pop_prev),
        }
        LOGGER.info(
            "Fitted model=%s; target=%s; common MAF rule=MAF >= %g; common reference-SNP universe=%d; "
            "target missing exclusions=%d; verification=%s.",
            model.model_type,
            args.target_annotation,
            common_maf_min,
            len(common),
            int((~eligible).sum()),
            verification_level,
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
            issues,
            QuantileH2OutputConfig(output_dir=output_dir, overwrite=True),
            metadata,
        )
        log_outputs(**output_paths)
    return QuantileH2Result(quantile_result, standardized, issues, metadata, output_paths)


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
