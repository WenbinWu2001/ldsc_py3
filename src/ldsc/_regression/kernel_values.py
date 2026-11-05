"""Read the results of the legacy estimator kernel as numbers.

The estimators ported from legacy LDSC (``estimators``) store the string
``'NA'`` where a value is undefined:

- ``intercept_se`` of an ``Hsq`` or ``Gencov`` fit with a fixed intercept;
- ``ratio`` and ``ratio_se`` of an ``Hsq`` fit whose mean chi^2 is at most 1;
- ``rg_ratio``, ``rg``, ``rg_se``, ``z`` and ``p`` of an ``RG`` fit whose
  heritability estimate for either trait is not positive, which the kernel
  flags with ``_negative_hsq`` (legacy LDSC reports "h2 out of bounds");
- ``Enrichment_p`` of a category that contains every SNP, in the
  overlap-aware table of ``Hsq._overlap_output``.

This module is the one place that replaces each sentinel with NaN, where the
package receives the kernel's result; the fits in ``fitting`` and the
overlap-aware table in ``overlap_model`` call it, so nothing downstream sees a
string. An ``RG`` fit also gains the flag ``h2_out_of_bounds``. The kernel
arithmetic is unchanged.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

_KERNEL_NA = "NA"


def _number(value):
    """Return NaN for the kernel's ``'NA'`` sentinel, otherwise ``value`` unchanged."""
    return math.nan if isinstance(value, str) and value == _KERNEL_NA else value


def numeric_ldsc_fit(fit):
    """Replace the ``'NA'`` sentinels of an ``Hsq`` or ``Gencov`` fit with NaN, in place, and return the fit.

    ``intercept_se`` becomes NaN when the intercept was fixed
    (``constrain_intercept``), and an ``Hsq`` fit's ``ratio`` and ``ratio_se``
    become NaN when its mean chi^2 is at most 1. A fixed-intercept ``Hsq`` fit
    has no ``ratio``; that stays so.
    """
    for name in ("intercept_se", "ratio", "ratio_se"):
        if hasattr(fit, name):
            setattr(fit, name, _number(getattr(fit, name)))
    return fit


def numeric_rg_fit(fit):
    """Replace the ``'NA'`` sentinels of an ``RG`` fit and its three component fits with NaN, in place.

    Sets ``h2_out_of_bounds``: True when either trait's heritability estimate
    is not positive, so ``rg_ratio``, ``rg``, ``rg_se``, ``z`` and ``p`` are NaN;
    the component fits ``hsq1``, ``hsq2`` and ``gencov`` keep their estimates.
    Returns the fit.
    """
    for component in (fit.hsq1, fit.hsq2, fit.gencov):
        numeric_ldsc_fit(component)
    for name in ("rg_ratio", "rg", "rg_se", "z", "p"):
        if hasattr(fit, name):
            setattr(fit, name, _number(getattr(fit, name)))
    fit.h2_out_of_bounds = bool(fit._negative_hsq)
    return fit


def numeric_overlap_table(table: pd.DataFrame) -> pd.DataFrame:
    """Return the overlap-aware table with ``Enrichment_p`` as float64, NaN for a category containing every SNP."""
    return table.assign(Enrichment_p=np.asarray([_number(value) for value in table["Enrichment_p"]], dtype=np.float64))
