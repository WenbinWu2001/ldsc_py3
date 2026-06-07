"""Query adjusted R² for SNP pairs from index-format reference panels.

Public surface for reading a ``ldsc build-ref-panel`` parquet panel by SNP pair:
the :class:`R2Panel` handle, the one-shot :func:`query_r2` wrapper, and the pure
:func:`unbiased_r2_to_pearson_r` converter. See
``docs/current/ref-panel-r2-query.md`` and the design spec
``docs/superpowers/specs/2026-06-06-ref-panel-r2-query-design.md``.
"""
from __future__ import annotations

import numpy as np


def unbiased_r2_to_pearson_r(r2_adj, n, sign=None):
    """Convert adjusted (unbiased) R² to a signed Pearson correlation r.

    Inverts the unbiased correction ``r2_adj = r2_raw - (1 - r2_raw) / (n - 2)``
    to recover the biased squared correlation, then takes the square root and
    applies ``sign``. Vectorized: ``r2_adj`` and ``sign`` may be scalars or
    array-likes, and ``n`` may be a scalar or array-like reference sample size.

    Parameters
    ----------
    r2_adj : float or array-like
        Adjusted (unbiased) squared correlation as stored in the panel.
    n : int or array-like
        Reference sample size used to estimate R² (``ldsc:n_samples``).
    sign : int, array-like, or None, optional
        Sign of the Pearson correlation (``+1``/``-1``). ``None`` returns the
        non-negative magnitude ``|r|``.

    Returns
    -------
    float or numpy.ndarray
        Signed Pearson r (or ``|r|`` when ``sign`` is ``None``).
    """
    r2_adj = np.asarray(r2_adj, dtype=np.float64)
    n = np.asarray(n, dtype=np.float64)
    r2_raw = (r2_adj * (n - 2.0) + 1.0) / (n - 1.0)
    r2_raw = np.clip(r2_raw, 0.0, 1.0)
    r = np.sqrt(r2_raw)
    if sign is not None:
        r = r * np.asarray(sign, dtype=np.float64)
    return r if r.ndim else r.item()
