"""Annotation overlap in partitioned-heritability models.

Assembles the overlap matrix of the annotations in one fitted model,
builds the overlap-aware category table, and detects collinear
annotations before a fit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import estimators as reg
from .kernel_values import numeric_overlap_table
from .._ldscore.overlap_table import LDScoreOverlap


# The one rename of the legacy overlap-aware table to the partitioned-h2 schema.
_LEGACY_OVERLAP_COLUMNS = {
    "Category": "category",
    "Prop._SNPs": "prop_snps",
    "Prop._h2": "prop_h2",
    "Prop._h2_std_error": "prop_h2_se",
    "Enrichment": "enrichment",
    "Enrichment_std_error": "enrichment_se",
    "Enrichment_p": "enrichment_p",
    "Coefficient": "coefficient",
    "Coefficient_std_error": "coefficient_se",
    "Coefficient_z-score": "coefficient_z",
}


def overlap_aware_category_table(hsq, overlap_matrix, m_annot, m_tot, category_names) -> pd.DataFrame:
    """Legacy overlap-aware annotation table plus LDSC3 augmentation columns, in the partitioned-h2 column names.

    The overlap-aware columns are produced verbatim by the ported
    ``Hsq._overlap_output``; ``kernel_values`` turns its ``'NA'``
    ``Enrichment_p`` (a category containing every SNP) into NaN. This function
    then renames them once to the lowercase ``PARTITIONED_H2_COLUMNS`` names
    and appends the one-sided ``coefficient_p`` (``H1: coefficient > 0``), the
    conditional ``category_h2_obs`` / ``category_h2_obs_se`` (from ``hsq.cat`` /
    ``hsq.cat_se``), and a per-model ``overlap_annot`` flag.
    """
    overlap_matrix = np.asarray(overlap_matrix, dtype=np.float64)
    table = numeric_overlap_table(reg.Hsq._overlap_output(
        hsq,
        list(category_names),
        overlap_matrix,
        np.asarray(m_annot, dtype=np.float64).reshape(1, -1),
        float(m_tot),
    ).reset_index(drop=True)).rename(columns=_LEGACY_OVERLAP_COLUMNS)
    coef = np.ravel(hsq.coef).astype(np.float64)
    coef_se = np.ravel(hsq.coef_se).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(coef_se > 0, coef / coef_se, np.nan)
    table["coefficient_p"] = stats.norm.sf(z)
    table["category_h2_obs"] = np.ravel(hsq.cat).astype(np.float64)
    table["category_h2_obs_se"] = np.ravel(hsq.cat_se).astype(np.float64)
    off_diagonal = overlap_matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    table["overlap_annot"] = bool(np.any(off_diagonal > 0))
    return table


def assemble_model_overlap(
    overlap: LDScoreOverlap, retained_columns: list[str], use_common: bool
) -> np.ndarray:
    """Build the ``KxK`` model overlap matrix for ``retained_columns`` by name.

    Baseline-baseline and baseline-query entries come from the baseline-rows
    block (using symmetry for query-baseline); the query-query diagonal comes
    from the stored self-overlap. ``use_common`` selects the common reference SNPs.
    """
    block = overlap.baseline_block_common if use_common else overlap.baseline_block_all
    qdiag = overlap.query_diagonal_common if use_common else overlap.query_diagonal_all
    if block is None:
        raise ValueError("common reference-SNP overlap requested but not stored")
    baseline_index = set(block.index)
    columns = list(retained_columns)
    is_baseline = np.array([column in baseline_index for column in columns], dtype=bool)
    baseline = [column for column, flag in zip(columns, is_baseline) if flag]
    queries = [column for column, flag in zip(columns, is_baseline) if not flag]
    missing = _first_missing_overlap_label(block, qdiag, columns, baseline_index)
    if missing is not None:
        raise KeyError(missing)
    k = len(columns)
    out = np.zeros((k, k), dtype=np.float64)
    base_pos, query_pos = np.flatnonzero(is_baseline), np.flatnonzero(~is_baseline)
    if baseline:
        out[base_pos, :] = block.loc[baseline, columns].to_numpy(dtype=np.float64)
    if queries:
        if baseline:  # query rows read the baseline block by symmetry
            out[np.ix_(query_pos, base_pos)] = block.loc[baseline, queries].to_numpy(dtype=np.float64).T
        # both query -> only the diagonal ca == cb is reachable per model
        out[np.ix_(query_pos, query_pos)] = qdiag.loc[queries].to_numpy(dtype=np.float64)[:, None]
    return out


def _first_missing_overlap_label(block, qdiag, columns, baseline_index):
    """Return the first label an entry-by-entry lookup in row-major order would miss, else None.

    Baseline rows read every column of the block, query rows read their own
    block column and their stored diagonal entry.
    """
    known = set(block.columns)
    queries = [column for column in columns if column not in baseline_index]
    if (len(queries) == len(columns) or set(columns) <= known) and set(queries) <= set(qdiag.index):
        return None
    for ca in columns:
        for cb in columns:
            if ca in baseline_index:
                label = None if cb in known else cb
            elif cb in baseline_index:
                label = None if ca in known else ca
            else:
                label = None if ca in qdiag.index else ca
            if label is not None:
                return label
    return None


def model_collinearity_error(X, columns, overlap_matrix=None, threshold: float = 1e5) -> str | None:
    """Return a fatal-error message when the LD-score design matrix is near-collinear.

    Uses the legacy condition-number threshold on the design matrix itself, and
    names the most correlated annotation pair (from the overlap Gram matrix) as
    the likely culprit. Without ``overlap_matrix`` the pair comes from the Gram
    matrix of the LD-score columns instead. Returns ``None`` when well
    conditioned or under-sized.
    The caller raises on a non-``None`` message: collinear annotations cannot be
    separated, so the partition would be untrustworthy and the run must abort
    rather than emit results the user has to second-guess.
    """
    X = np.asarray(X, dtype=np.float64)
    if X.shape[1] < 2:
        return None
    # cond(X) = sqrt(cond(X.T @ X)); the Gram matrix is small (kxk), so this avoids an
    # SVD of the tall design matrix per fitted model.
    condition_number = float(np.sqrt(np.linalg.cond(X.T @ X)))
    if condition_number <= threshold:
        return None
    O = X.T @ X if overlap_matrix is None else np.asarray(overlap_matrix, dtype=np.float64)
    scale = np.sqrt(np.clip(np.diag(O), 1e-300, None))
    corr = O / np.outer(scale, scale)
    np.fill_diagonal(corr, 0.0)
    i, j = np.unravel_index(int(np.argmax(np.abs(corr))), corr.shape)
    pair = (f"The LD scores of annotations '{columns[i]}' and '{columns[j]}' are almost proportional"
            if overlap_matrix is None else
            f"Annotations '{columns[i]}' and '{columns[j]}' overlap almost entirely")
    return (
        f"LD score regression aborted: the design matrix is collinear (condition number "
        f"{int(condition_number)} > {int(threshold)}). {pair} (r={corr[i, j]:.3f}); their coefficients "
        f"cannot be separated. Drop one of them or coarsen the annotation, then re-run."
    )
