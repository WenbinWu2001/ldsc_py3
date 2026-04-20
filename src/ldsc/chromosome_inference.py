"""Central chromosome normalization and ordering registry for the package."""

from __future__ import annotations

import logging
import math
import re
import warnings


LOGGER = logging.getLogger("LDSC.chromosomes")
_WARNED_NORMALIZATIONS: set[tuple[str, str, str, str]] = set()
STANDARD_CHROMOSOMES = tuple([str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"])
_CHROMOSOME_ORDER = {chrom: idx for idx, chrom in enumerate(STANDARD_CHROMOSOMES, start=1)}


def _emit_normalization_warning(
    *,
    raw: str,
    canonical: str,
    context: str | None,
    kind: str,
) -> None:
    key = (context or "", raw, canonical, kind)
    if key in _WARNED_NORMALIZATIONS:
        return
    _WARNED_NORMALIZATIONS.add(key)

    suffix = "" if not context else f" in {context}"
    if kind == "sex_code":
        message = (
            f"Inferred canonical chromosome '{canonical}' from numeric sex chromosome code "
            f"'{raw}'{suffix}."
        )
        LOGGER.warning(message)
    else:
        message = f"Inferred canonical chromosome '{canonical}' from input chromosome '{raw}'{suffix}."
        LOGGER.info(message)
    warnings.warn(message, UserWarning, stacklevel=3)


def _normalize_numeric_chromosome(raw: str, text: str, *, context: str | None) -> str:
    try:
        number = float(text)
    except ValueError as exc:
        raise ValueError(
            f"Unsupported chromosome label {text!r}. Allowed chromosome labels are "
            "1-22, X, Y, M, and MT."
        ) from exc

    if not math.isfinite(number) or not number.is_integer():
        raise ValueError(
            f"Unsupported chromosome label {text!r}. Allowed chromosome labels are "
            "1-22, X, Y, M, and MT."
        )

    integer = int(number)
    if 1 <= integer <= 22:
        canonical = str(integer)
    elif integer == 23:
        canonical = "X"
        _emit_normalization_warning(raw=raw, canonical=canonical, context=context, kind="sex_code")
        return canonical
    elif integer == 24:
        canonical = "Y"
        _emit_normalization_warning(raw=raw, canonical=canonical, context=context, kind="sex_code")
        return canonical
    elif integer in {25, 26}:
        raise ValueError(
            f"Numeric chromosome code {text!r} is not supported. If you mean a mitochondrial "
            "chromosome, encode it explicitly as 'M' or 'MT'."
        )
    else:
        raise ValueError(
            f"Unsupported chromosome label {text!r}. Allowed chromosome labels are "
            "1-22, X, Y, M, and MT."
        )

    if raw != canonical:
        _emit_normalization_warning(raw=raw, canonical=canonical, context=context, kind="canonical")
    return canonical


def normalize_chromosome(value: object, *, context: str | None = None) -> str:
    """Normalize one raw chromosome token into the package-wide canonical set."""
    raw = str(value).strip()
    if not raw:
        raise ValueError("Encountered an empty chromosome label.")

    stripped = re.sub(r"^chr", "", raw, flags=re.IGNORECASE).strip()
    if not stripped:
        raise ValueError(f"Encountered an empty chromosome label after removing the 'chr' prefix from {raw!r}.")

    upper = stripped.upper()
    if upper in {"X", "Y", "M", "MT"}:
        canonical = upper
        if raw != canonical:
            _emit_normalization_warning(raw=raw, canonical=canonical, context=context, kind="canonical")
        return canonical

    return _normalize_numeric_chromosome(raw, stripped, context=context)


def chrom_sort_key(chrom: object) -> tuple[int, object]:
    """Return the stable package-wide chromosome ordering key."""
    canonical = normalize_chromosome(chrom)
    order = _CHROMOSOME_ORDER.get(canonical)
    if order is not None:
        return (0, order)
    return (1, canonical)
