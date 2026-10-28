"""Annotation overlap in partitioned-heritability models.

Assembles the overlap matrix of the annotations in one fitted model,
builds the overlap-aware category table, and detects collinear
annotations before a fit.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from ._kernel import regression as reg
from ._ldscore.overlap_table import LDScoreOverlap


def overlap_aware_category_table(hsq, overlap_matrix, m_annot, m_tot, category_names) -> pd.DataFrame:
    """Legacy overlap-aware annotation table plus LDSC3 augmentation columns.

    The overlap-aware columns are produced verbatim by the ported
    ``Hsq._overlap_output``. This function then renames the z-score column,
    appends the one-sided ``Coefficient_p`` (``H1: coefficient > 0``), the
    conditional ``Category_h2`` / ``Category_h2_std_error`` (from ``hsq.cat`` /
    ``hsq.cat_se``), and a per-model ``overlap_aware`` flag.
    """
    overlap_matrix = np.asarray(overlap_matrix, dtype=np.float64)
    table = reg.Hsq._overlap_output(
        hsq,
        list(category_names),
        overlap_matrix,
        np.asarray(m_annot, dtype=np.float64).reshape(1, -1),
        float(m_tot),
        True,
    ).reset_index(drop=True)
    table = table.rename(columns={"Coefficient_z-score": "Coefficient_z"})
    # Legacy emits the string 'NA' for a category that contains every SNP; use a numeric NaN.
    table["Enrichment_p"] = pd.to_numeric(table["Enrichment_p"], errors="coerce")
    coef = np.ravel(hsq.coef).astype(np.float64)
    coef_se = np.ravel(hsq.coef_se).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.where(coef_se > 0, coef / coef_se, np.nan)
    table["Coefficient_p"] = stats.norm.sf(z)
    table["Category_h2"] = np.ravel(hsq.cat).astype(np.float64)
    table["Category_h2_std_error"] = np.ravel(hsq.cat_se).astype(np.float64)
    off_diagonal = overlap_matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    table["overlap_aware"] = bool(np.any(off_diagonal > 0))
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
    k = len(retained_columns)
    out = np.zeros((k, k), dtype=np.float64)
    for a, ca in enumerate(retained_columns):
        for b, cb in enumerate(retained_columns):
            if ca in baseline_index:
                out[a, b] = float(block.at[ca, cb])
            elif cb in baseline_index:
                out[a, b] = float(block.at[cb, ca])
            else:  # both query -> only the diagonal ca == cb is reachable per model
                out[a, b] = float(qdiag.at[ca])
    return out


def model_collinearity_error(X, columns, overlap_matrix, threshold: float = 1e5) -> str | None:
    """Return a fatal-error message when the LD-score design matrix is near-collinear.

    Uses the legacy condition-number threshold on the design matrix itself, and
    names the most correlated annotation pair (from the overlap Gram matrix) as
    the likely culprit. Returns ``None`` when well conditioned or under-sized.
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
    O = np.asarray(overlap_matrix, dtype=np.float64)
    scale = np.sqrt(np.clip(np.diag(O), 1e-300, None))
    corr = O / np.outer(scale, scale)
    np.fill_diagonal(corr, 0.0)
    i, j = np.unravel_index(int(np.argmax(np.abs(corr))), corr.shape)
    return (
        f"LD score regression aborted: the design matrix is collinear (condition number "
        f"{int(condition_number)} > {int(threshold)}). Annotations '{columns[i]}' and "
        f"'{columns[j]}' overlap almost entirely (r={corr[i, j]:.3f}); their coefficients "
        f"cannot be separated. Drop one of them or coarsen the annotation, then re-run."
    )
