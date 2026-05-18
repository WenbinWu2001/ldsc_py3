"""Build the compact HapMap3 chromosome-position inference reference."""

from __future__ import annotations

import argparse
from os import PathLike
from pathlib import Path

import numpy as np
import pandas as pd

from ._kernel.snp_identity import allele_set_series
from .hm3 import _load_hm3_curated_map_from_path, packaged_hm3_curated_map_path

REFERENCE_COLUMNS = ["CHR", "hg19_POS", "hg38_POS"]


def build_hm3_chr_pos_reference(
    curated_map_path: str | PathLike[str],
    output_path: str | PathLike[str],
    *,
    snps_per_chromosome: int = 500,
    min_maf: float = 0.2,
) -> pd.DataFrame:
    """
    Build and write the compact HM3 coordinate reference used for build inference.

    The generated table contains autosomal, common, non-strand-ambiguous SNPs
    that have unique hg19 and hg38 coordinates and differ between builds. Rows
    are evenly sampled by hg38 position within each chromosome.
    """
    if snps_per_chromosome <= 0:
        raise ValueError("snps_per_chromosome must be positive.")
    curated = _load_hm3_curated_map_from_path(curated_map_path, preserve_extra_columns=True)
    filtered = _filter_reference_candidates(curated, min_maf=min_maf)
    sampled = _sample_evenly_by_chromosome(filtered, snps_per_chromosome=snps_per_chromosome)
    reference = sampled.loc[:, REFERENCE_COLUMNS].reset_index(drop=True)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    reference.to_csv(output, sep="\t", index=False, compression="gzip")
    return reference


def _filter_reference_candidates(frame: pd.DataFrame, *, min_maf: float) -> pd.DataFrame:
    """Return HM3 rows eligible for compact genome-build inference."""
    work = frame.copy()
    work = work.loc[work["CHR"].isin([str(chrom) for chrom in range(1, 23)])].copy()
    work["MAF"] = pd.to_numeric(work["MAF"], errors="coerce")
    work = work.loc[work["MAF"] >= min_maf].copy()
    allele_set, allele_reasons = allele_set_series(work, context="HM3 chromosome-position reference build")
    work = work.loc[allele_reasons.isna() & allele_set.notna()].copy()
    work = work.loc[work["hg19_POS"] != work["hg38_POS"]].copy()
    work = _drop_duplicate_coordinate_groups(work, ["CHR", "hg19_POS"])
    work = _drop_duplicate_coordinate_groups(work, ["CHR", "hg38_POS"])
    return work


def _drop_duplicate_coordinate_groups(frame: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Drop all rows in duplicate coordinate groups."""
    duplicated = frame.duplicated(columns, keep=False)
    return frame.loc[~duplicated].copy()


def _sample_evenly_by_chromosome(frame: pd.DataFrame, *, snps_per_chromosome: int) -> pd.DataFrame:
    """Return up to ``snps_per_chromosome`` rows per chromosome, spaced by hg38 position."""
    sampled: list[pd.DataFrame] = []
    for chrom in [str(value) for value in range(1, 23)]:
        group = frame.loc[frame["CHR"] == chrom].sort_values("hg38_POS").reset_index(drop=True)
        if group.empty:
            continue
        if len(group) <= snps_per_chromosome:
            sampled.append(group)
            continue
        indices = np.linspace(0, len(group) - 1, snps_per_chromosome, dtype=int)
        sampled.append(group.iloc[indices].copy())
    if not sampled:
        return pd.DataFrame(columns=frame.columns)
    return pd.concat(sampled, ignore_index=True)


def main(argv: list[str] | None = None) -> int:
    """Run the HM3 reference builder from the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--curated-map",
        default=packaged_hm3_curated_map_path(),
        help="Input curated HM3 map TSV/TSV.GZ. Defaults to the packaged map.",
    )
    parser.add_argument("--output", required=True, help="Output hm3_chr_pos_reference.tsv.gz path.")
    parser.add_argument("--snps-per-chromosome", type=int, default=500)
    parser.add_argument("--min-maf", type=float, default=0.2)
    args = parser.parse_args(argv)
    reference = build_hm3_chr_pos_reference(
        args.curated_map,
        args.output,
        snps_per_chromosome=args.snps_per_chromosome,
        min_maf=args.min_maf,
    )
    print(f"Wrote {len(reference)} SNPs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
