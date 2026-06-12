"""Low-level annotation table and BED projection helpers.

This private module contains concrete-file parsing helpers and BED
intersection primitives used by the public ``ldsc.annotation_builder``
workflow. Public result objects, CLI parsing, path-token orchestration, global
configuration handling, and output policy belong in ``ldsc.annotation_builder``.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from ..chromosome_inference import chrom_sort_key, normalize_chromosome
from ..errors import LDSCDependencyError, LDSCInputError, LDSCInternalError
from ..path_resolution import ensure_output_parent_directory, resolve_scalar_path
from . import regions as kernel_regions


@dataclass(frozen=True)
class _BaselineRow:
    """Minimal metadata for one baseline SNP row during BED projection.

    BED overlap depends only on genomic position; CM (genetic-map distance) is
    population-specific and never used here.
    """

    chrom: str
    pos: int
    snp: str


def _get_pybedtools():
    """Import ``pybedtools`` or raise a user-facing dependency error."""
    try:
        import pybedtools
    except ImportError as exc:  # pragma: no cover - dependency check
        raise LDSCDependencyError(
            "annotate could not project BED intervals because pybedtools is not installed. "
            "Most likely BED-based annotation projection was requested in an environment missing "
            "pybedtools or the bedtools binary. Install pybedtools and bedtools, then retry."
        ) from exc
    return pybedtools


def _read_text_table(path: str) -> pd.DataFrame:
    """Read a whitespace-delimited annotation table with optional gzip compression."""
    compression = "gzip" if str(path).endswith(".gz") else None
    return pd.read_csv(path, sep=r"\s+", compression=compression)


def _annotation_shard_chromosome(path: str | Path) -> str | None:
    """Extract the chromosome token from one LDSC-style shard filename."""
    match = re.match(
        r"^.+\.(?P<chrom>\d+|X|Y|M|MT)\.(?:annot|txt|tsv)(?:\.gz)?$",
        Path(path).name,
        flags=re.IGNORECASE,
    )
    if match is None:
        return None
    return normalize_chromosome(match.group("chrom"))


def _chrom_sort_key(chrom: str) -> tuple[int, str]:
    """Return the stable chromosome ordering used by annotation workflows."""
    return chrom_sort_key(chrom)


def _to_bed_chromosome(chrom: object) -> str:
    """Convert one chromosome label into UCSC-style BED notation."""
    return "chr" + normalize_chromosome(chrom)


def _write_baseline_bed(rows: Sequence[_BaselineRow], path: Path) -> Path:
    """Write baseline SNP rows as single-base BED intervals for overlap queries."""
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            start = row.pos - 1
            if start < 0:
                raise LDSCInputError(
                    f"annotate could not project BED annotations because baseline SNP {row.snp} has invalid POS={row.pos}. "
                    "Most likely the baseline annotation file has non-positive genomic positions. "
                    "Regenerate the baseline annotation file with 1-based positive POS values."
                )
            handle.write(f"{_to_bed_chromosome(row.chrom)}\t{start}\t{row.pos}\t{row.snp}\n")
    return path


def _write_normalized_bed(in_path: Path, out_path: Path, *, bed_padding_bp: int = 0) -> Path:
    """Normalize one BED input into a plain tab-delimited UCSC-style BED file."""
    with (gzip.open(in_path, "rt") if in_path.suffix.lower() == ".gz" else open(in_path, "rt")) as src:
        with out_path.open("w", encoding="utf-8") as dst:
            for row in kernel_regions.parse_bed_text(src.read(), label=str(in_path)):
                chrom = _to_bed_chromosome(row.chrom)
                start = max(0, row.start - bed_padding_bp)
                end = row.end + bed_padding_bp
                dst.write("\t".join([chrom, str(start), str(end), *row.extra_fields]) + "\n")
    return out_path


def _validate_and_convert_intersection(results, baseline_rows: Sequence[_BaselineRow], sanity_mode: str) -> list[bool]:
    """Validate a pybedtools intersection result and convert it into a row-aligned mask."""
    mask: list[bool] = []
    for idx, feature in enumerate(results):
        if idx >= len(baseline_rows):
            raise LDSCInternalError(
                "annotate BED projection returned more intersection rows than baseline SNPs. "
                "Most likely the BED intersection backend changed row cardinality unexpectedly. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            )
        row = baseline_rows[idx]
        fields = feature.fields
        overlap_count = int(fields[-1])
        if fields[3] != row.snp:
            raise LDSCInternalError(
                f"annotate BED projection row order mismatch at index {idx}: expected SNP {row.snp}, got {fields[3]}. "
                "Most likely the BED intersection output no longer matches the baseline template order. "
                "Re-run with `--log-level DEBUG` and report the traceback."
            )
        if sanity_mode == "chr_pos":
            feature_chr = _to_bed_chromosome(fields[0])
            feature_pos = int(fields[2])
            if feature_chr != _to_bed_chromosome(row.chrom) or feature_pos != row.pos:
                raise LDSCInternalError(
                    f"annotate BED projection coordinate mismatch at index {idx}: expected ({row.chrom}, {row.pos}), got ({fields[0]}, {fields[2]}). "
                    "Most likely the BED intersection output desynchronized from the baseline SNP template. "
                    "Re-run with `--log-level DEBUG` and report the traceback."
                )
        mask.append(overlap_count > 0)
    if len(mask) != len(baseline_rows):
        raise LDSCInternalError(
            f"annotate BED projection returned {len(mask)} rows for {len(baseline_rows)} baseline SNPs. "
            "Most likely the BED intersection output omitted or added template rows. "
            "Re-run with `--log-level DEBUG` and report the traceback."
        )
    return mask


def _compute_bed_overlap_mask(baseline_rows: Sequence[_BaselineRow], baseline_bed, bed_path: Path) -> list[bool]:
    """Compute the overlap mask for one annotation BED against the baseline SNP template."""
    results = baseline_bed.intersect(str(bed_path), c=True, wa=True)
    return _validate_and_convert_intersection(results, baseline_rows, sanity_mode="rsid")


def _compute_bed_query_columns(
    metadata: pd.DataFrame,
    bed_paths: Sequence[Path],
    tempdir: Path,
) -> pd.DataFrame:
    """Project BED overlaps onto ``metadata`` without writing ``.annot.gz`` files."""
    pybedtools = _get_pybedtools()
    rows = [
        _BaselineRow(
            chrom=str(row.CHR),
            pos=int(row.POS),
            snp=str(row.SNP),
        )
        for row in metadata.itertuples(index=False)
    ]
    baseline_bed_path = _write_baseline_bed(rows, tempdir / "bed_query_baseline.bed")
    baseline_bed = pybedtools.BedTool(str(baseline_bed_path))
    try:
        return pd.DataFrame(
            {
                path.stem: np.asarray(_compute_bed_overlap_mask(rows, baseline_bed, path), dtype=np.float32)
                for path in bed_paths
            },
            index=metadata.index,
        )
    finally:
        pybedtools.cleanup(remove_all=True)


def _make_single_annotation_file(
    bimfile: str | Path,
    annot_file: str | Path,
    bed_for_annot,
) -> Path:
    """Project one BED-like input onto one BIM file and write a legacy ``.annot`` file."""
    pybedtools = _get_pybedtools()
    bimfile = Path(resolve_scalar_path(bimfile, label="PLINK BIM file"))
    annot_file = ensure_output_parent_directory(annot_file, label="annot_file")

    df_bim = pd.read_csv(
        bimfile,
        sep=r"\s+",
        usecols=[0, 1, 2, 3],
        names=["CHR", "SNP", "CM", "POS"],
        header=None,
    )
    iter_bim = [[_to_bed_chromosome(chrom), int(pos) - 1, int(pos)] for chrom, pos in np.array(df_bim[["CHR", "POS"]])]
    bim_bed = pybedtools.BedTool(iter_bim)
    annot_bed = bim_bed.intersect(bed_for_annot)
    pos = [feature.start + 1 for feature in annot_bed]
    df_int = pd.DataFrame({"POS": pos, "ANNOT": 1})
    df_annot = pd.merge(df_bim, df_int, how="left", on="POS")
    df_annot.fillna(0, inplace=True)
    df_annot = df_annot[["ANNOT"]].astype(int)
    annot_file.parent.mkdir(parents=True, exist_ok=True)
    if str(annot_file).endswith(".gz"):
        with gzip.open(annot_file, "wt") as handle:
            df_annot.to_csv(handle, sep="\t", index=False)
    else:
        df_annot.to_csv(annot_file, sep="\t", index=False)
    return annot_file
