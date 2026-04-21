"""Genome-build and coordinate-basis inference for ``chr_pos`` inputs."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Iterable

import pandas as pd

from .chromosome_inference import normalize_chromosome
from .column_inference import normalize_genome_build, normalize_snp_identifier_mode


AUTO_GENOME_BUILD = "auto"
DOMINANCE_THRESHOLD = 0.99
MIN_INSPECTED_REFERENCE_SNPS = 200
REFERENCE_RESOURCE_PATH = "data/hm3_chr_pos_reference.tsv.gz"

_HYPOTHESIS_ORDER = (
    ("hg19_1based", "hg19", "1-based", "hg19_POS", 0),
    ("hg19_0based", "hg19", "0-based", "hg19_POS", -1),
    ("hg38_1based", "hg38", "1-based", "hg38_POS", 0),
    ("hg38_0based", "hg38", "0-based", "hg38_POS", -1),
)


@dataclass(frozen=True)
class ChrPosBuildInference:
    """Resolved build/basis decision for one ``chr_pos`` input table."""

    genome_build: str
    coordinate_basis: str
    inspected_snp_count: int
    match_counts: dict[str, int]
    match_fractions: dict[str, float]
    summary_message: str


def is_auto_genome_build(genome_build: str | None) -> bool:
    """Return ``True`` when ``genome_build`` requests automatic inference."""
    return normalize_genome_build(genome_build) == AUTO_GENOME_BUILD


def validate_auto_genome_build_mode(snp_identifier: str, genome_build: str | None) -> None:
    """Reject ``genome_build='auto'`` outside ``chr_pos`` mode."""
    if not is_auto_genome_build(genome_build):
        return
    if normalize_snp_identifier_mode(snp_identifier) != "chr_pos":
        raise ValueError("genome_build='auto' is only supported when snp_identifier='chr_pos'.")


@lru_cache(maxsize=1)
def load_packaged_reference_table() -> pd.DataFrame:
    """Load the packaged HM3 reference subset used for auto inference."""
    resource = resources.files("ldsc").joinpath(REFERENCE_RESOURCE_PATH)
    with resources.as_file(resource) as path:
        df = pd.read_csv(path, sep="\t", compression="gzip")
    return df.loc[:, ["CHR", "hg19_POS", "hg38_POS"]].copy()


def resolve_chr_pos_table(
    df: pd.DataFrame,
    *,
    context: str,
    reference_table: pd.DataFrame | None = None,
    logger=None,
) -> tuple[pd.DataFrame, ChrPosBuildInference]:
    """
    Infer genome build and coordinate basis for a ``CHR``/``POS`` table.

    The returned frame preserves the input columns but normalizes ``CHR`` to the
    package-wide canonical labels. When the inferred basis is ``0-based``, the
    returned ``POS`` values are shifted to canonical 1-based coordinates.
    """
    if not {"CHR", "POS"}.issubset(df.columns):
        raise ValueError(f"{context} must contain CHR and POS columns for auto genome-build inference.")

    reference = load_packaged_reference_table() if reference_table is None else reference_table.copy()
    normalized = df.copy()
    normalized["CHR"] = normalized["CHR"].map(lambda value: normalize_chromosome(value, context=context))
    normalized["POS"] = pd.to_numeric(normalized["POS"], errors="raise").astype("int64")

    inference = infer_chr_pos_build(
        normalized.loc[:, ["CHR", "POS"]],
        context=context,
        reference_table=reference,
    )
    if inference.coordinate_basis == "0-based":
        normalized["POS"] = normalized["POS"] + 1
    if logger is not None and inference.coordinate_basis == "0-based":
        logger.warning(inference.summary_message)
    elif logger is not None:
        logger.info(inference.summary_message)
    return normalized, inference


def infer_chr_pos_build(
    df: pd.DataFrame,
    *,
    context: str,
    reference_table: pd.DataFrame | None = None,
) -> ChrPosBuildInference:
    """Infer build/basis from an already-normalized ``CHR``/``POS`` table."""
    reference = load_packaged_reference_table() if reference_table is None else reference_table.copy()
    keys = _input_keys(df)
    hypothesis_sets = _build_hypothesis_sets(reference)
    match_counts = {name: 0 for name, *_rest in _HYPOTHESIS_ORDER}
    inspected = 0
    for key in keys:
        matched = [name for name, values in hypothesis_sets.items() if key in values]
        if len(matched) != 1:
            continue
        inspected += 1
        match_counts[matched[0]] += 1

    match_fractions = {
        name: (0.0 if inspected == 0 else match_counts[name] / inspected)
        for name, *_rest in _HYPOTHESIS_ORDER
    }
    if inspected < MIN_INSPECTED_REFERENCE_SNPS:
        raise ValueError(
            f"Insufficient evidence to infer genome build for {context}: "
            f"matched {inspected} reference SNPs; need at least {MIN_INSPECTED_REFERENCE_SNPS}."
        )

    best_name, best_fraction = max(
        match_fractions.items(),
        key=lambda item: (item[1], match_counts[item[0]], item[0]),
    )
    if best_fraction < DOMINANCE_THRESHOLD:
        raise ValueError(
            f"Genome build for {context} could not be inferred confidently: "
            f"best hypothesis {best_name} matched {match_counts[best_name]}/{inspected} "
            f"({best_fraction:.1%}); need at least {DOMINANCE_THRESHOLD:.0%}."
        )

    build, basis = _hypothesis_metadata(best_name)
    summary_message = _build_summary_message(build, basis, inspected, match_counts, match_fractions)
    return ChrPosBuildInference(
        genome_build=build,
        coordinate_basis=basis,
        inspected_snp_count=inspected,
        match_counts=match_counts,
        match_fractions=match_fractions,
        summary_message=summary_message,
    )


def _input_keys(df: pd.DataFrame) -> set[tuple[str, int]]:
    """Return the deduplicated, valid input keys used for scoring."""
    clean = df.loc[:, ["CHR", "POS"]].dropna().copy()
    clean["POS"] = pd.to_numeric(clean["POS"], errors="coerce")
    clean = clean.loc[clean["POS"].notna()].copy()
    clean["POS"] = clean["POS"].astype("int64")
    clean = clean.loc[clean["POS"] >= 0].drop_duplicates(subset=["CHR", "POS"])
    return {(str(chrom), int(pos)) for chrom, pos in clean.itertuples(index=False, name=None)}


def _build_hypothesis_sets(reference: pd.DataFrame) -> dict[str, set[tuple[str, int]]]:
    """Convert the packaged reference subset into hypothesis-specific key sets."""
    out: dict[str, set[tuple[str, int]]] = {}
    for name, _build, _basis, pos_column, offset in _HYPOTHESIS_ORDER:
        pos = pd.to_numeric(reference[pos_column], errors="raise").astype("int64") + offset
        valid = pos >= 0
        out[name] = {
            (str(chrom), int(position))
            for chrom, position in zip(reference.loc[valid, "CHR"], pos.loc[valid])
        }
    return out


def _hypothesis_metadata(name: str) -> tuple[str, str]:
    """Return ``(build, basis)`` for one hypothesis name."""
    for hypothesis_name, build, basis, _pos_column, _offset in _HYPOTHESIS_ORDER:
        if hypothesis_name == name:
            return build, basis
    raise KeyError(name)


def _build_summary_message(
    genome_build: str,
    coordinate_basis: str,
    inspected: int,
    match_counts: dict[str, int],
    match_fractions: dict[str, float],
) -> str:
    """Build the user-facing inference message emitted by workflows."""
    parts: list[str] = []
    for name, build, basis, _pos_column, _offset in _HYPOTHESIS_ORDER:
        label = f"{build}/{basis}"
        parts.append(f"{label}={match_counts[name]}/{inspected} ({match_fractions[name]:.1%})")
    basis_note = "converted positions to canonical 1-based coordinates" if coordinate_basis == "0-based" else "positions already matched canonical 1-based coordinates"
    return (
        f"Auto-inferred genome build {genome_build} with {coordinate_basis} coordinates from "
        f"{inspected} reference SNPs; {basis_note}. Hypothesis support: "
        + ", ".join(parts)
        + "."
    )
