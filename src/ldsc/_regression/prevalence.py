"""
prevalence.py

Overview
--------
Parse and validate the binary-trait prevalence inputs (``--samp-prev``,
``--pop-prev``, and the rg ``--prevalence-manifest``) into a normalized per-trait
structure consumed by the regression summary layer for observed-to-liability
scale conversion. Validation is fail-fast: range and consistency are checked
here, before LD scores are loaded.

Key Functions
-------------
resolve_rg_prevalences :
    Dispatch rg inputs (positional comma-list xor manifest) to a per-trait list
    aligned to the resolved trait order, or ``None`` for observed scale.

Design Notes
------------
- Normalized unit is ``TraitPrevalence = (samp_prev | None, pop_prev | None)``;
  ``None`` marks a quantitative trait (no liability scale).
- Range/consistency validation delegates to ``ldsc.config.normalize_prevalence_pair``,
  which ``RegressionConfig`` also applies to the scalar ``h2``/``partitioned-h2``
  pair. Both reach ``ldsc.config.validate_prevalences``, where the (0, 1)-or-NaN
  rule lives in one place.
- Manifest lookup is by exact munged trait name; supersets are allowed, while
  missing resolved traits and (under the manifest scheme) duplicate resolved
  names are hard errors.

Dependencies
------------
- numpy
"""
from __future__ import annotations

from pathlib import Path

from ..config import normalize_prevalence_pair
from ..errors import LDSCUsageError

TraitPrevalence = tuple[float | None, float | None]

_MISSING_TOKENS = {"", "nan", "na", "none"}


def _token_to_value(token: str) -> float | None:
    """Map one comma-list/manifest token to a float or None (quantitative)."""
    text = token.strip()
    if text.lower() in _MISSING_TOKENS:
        return None
    try:
        return float(text)
    except ValueError as exc:
        raise LDSCUsageError(
            f"Could not parse prevalence value {token!r} as a number. "
            "Most likely a malformed `--samp-prev`/`--pop-prev` entry or manifest cell. "
            "Use a probability in (0, 1), or `nan` for a quantitative trait."
        ) from exc


def parse_positional_prevalences(samp_prev: str, pop_prev: str, n_traits: int) -> list[TraitPrevalence]:
    """Parse aligned comma-separated ``--samp-prev``/``--pop-prev`` lists for rg.

    Both strings must split into exactly ``n_traits`` tokens. Returns a list of
    normalized ``(P, K)`` pairs aligned to the resolved trait order.
    """
    if samp_prev is None or pop_prev is None:
        raise LDSCUsageError(
            "rg liability conversion needs both `--samp-prev` and `--pop-prev` (or neither). "
            "Most likely only one positional list was supplied."
        )
    sp = [_token_to_value(t) for t in samp_prev.split(",")]
    pp = [_token_to_value(t) for t in pop_prev.split(",")]
    for name, values in (("--samp-prev", sp), ("--pop-prev", pp)):
        if len(values) != n_traits:
            raise LDSCUsageError(
                f"rg {name} has {len(values)} values but {n_traits} traits resolved from "
                "`--sumstats-sources`. The comma-separated list must have exactly one entry per "
                "resolved trait, in resolved order. Re-check the order printed in the run log, or "
                "use `--prevalence-manifest` to match by trait name."
            )
    return [normalize_prevalence_pair(s, p) for s, p in zip(sp, pp)]


def parse_manifest_prevalences(
    manifest_path: str, trait_names: list[str], trait_paths: list[str]
) -> list[TraitPrevalence]:
    """Look up ``(P, K)`` per resolved trait from a prevalence manifest TSV.

    The manifest is whitespace/tab-delimited with a ``trait_name samp_prev
    pop_prev`` header; ``#`` comment lines are ignored. Extra rows (traits not in
    this run) are allowed; every resolved trait must have a matching row, matched
    by exact name. Returns pairs aligned to ``trait_names``.
    """
    path = Path(manifest_path)
    if not path.is_file():
        raise LDSCUsageError(
            f"rg could not read `--prevalence-manifest` at '{manifest_path}': not a file. "
            "Pass a whitespace- or tab-delimited table with columns trait_name, samp_prev, pop_prev."
        )
    table: dict[str, TraitPrevalence] = {}
    idx: dict[str, int] | None = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        if raw.lstrip().startswith("#") or not raw.strip():
            continue
        fields = raw.split()
        if idx is None:
            header = [f.strip().lower() for f in fields]
            for required in ("trait_name", "samp_prev", "pop_prev"):
                if required not in header:
                    raise LDSCUsageError(
                        f"Prevalence manifest '{manifest_path}' is missing a '{required}' column. "
                        "Required header columns are trait_name, samp_prev, pop_prev."
                    )
            idx = {col: header.index(col) for col in ("trait_name", "samp_prev", "pop_prev")}
            continue
        name = fields[idx["trait_name"]].strip()
        table[name] = normalize_prevalence_pair(
            _token_to_value(fields[idx["samp_prev"]]),
            _token_to_value(fields[idx["pop_prev"]]),
        )
    if idx is None:
        raise LDSCUsageError(
            f"Prevalence manifest '{manifest_path}' has no header row. "
            "Add a header naming trait_name, samp_prev, pop_prev."
        )
    missing = [t for t in trait_names if t not in table]
    if missing:
        annotated = ", ".join(
            f"{t} ({p})" for t, p in zip(trait_names, trait_paths) if t in missing
        )
        raise LDSCUsageError(
            f"Prevalence manifest '{manifest_path}' is missing rows for resolved trait(s): {annotated}. "
            "Add a row per trait (extra rows for other traits are allowed), or use the positional "
            "`--samp-prev`/`--pop-prev` lists instead."
        )
    return [table[t] for t in trait_names]


def _reject_duplicate_names(trait_names: list[str], trait_paths: list[str]) -> None:
    """Raise if resolved munged trait names are not unique (ambiguous by-name lookup)."""
    seen: dict[str, list[str]] = {}
    for name, path in zip(trait_names, trait_paths):
        seen.setdefault(name, []).append(path)
    dupes = {name: paths for name, paths in seen.items() if len(paths) > 1}
    if dupes:
        detail = "; ".join(f"{name}: {', '.join(paths)}" for name, paths in dupes.items())
        raise LDSCUsageError(
            f"rg `--prevalence-manifest` cannot match by name because resolved traits have duplicate "
            f"munged names ({detail}). Either use the positional `--samp-prev`/`--pop-prev` lists "
            "(matched by position), or re-munge so each input has a distinct `trait_name` metadata value."
        )


def resolve_rg_prevalences(
    samp_prev, pop_prev, manifest_path, trait_names: list[str], trait_paths: list[str]
) -> list[TraitPrevalence] | None:
    """Resolve rg prevalence from exactly one supplied scheme, or ``None``.

    Schemes are mutually exclusive: positional comma-lists (``samp_prev`` and
    ``pop_prev``) xor ``manifest_path``. Returns a per-trait list aligned to
    ``trait_names``, or ``None`` when no scheme is supplied (observed scale).
    """
    positional = samp_prev is not None or pop_prev is not None
    if positional and manifest_path is not None:
        raise LDSCUsageError(
            "rg received both `--prevalence-manifest` and `--samp-prev`/`--pop-prev`. "
            "These schemes are mutually exclusive; choose one."
        )
    if manifest_path is not None:
        _reject_duplicate_names(trait_names, trait_paths)
        return parse_manifest_prevalences(manifest_path, trait_names, trait_paths)
    if positional:
        return parse_positional_prevalences(samp_prev, pop_prev, len(trait_names))
    return None
