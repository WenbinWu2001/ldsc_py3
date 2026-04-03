"""Create LDSC `.annot` files from one or more BED annotations.

This module converts interval-based BED annotations into chromosome-specific
LDSC annotation files using a directory of baseline `.annot` templates as the
reference SNP universe. Each output row preserves the order and core metadata
(`CHR`, `BP`, `SNP`, `CM`) of the matching baseline template and appends one or
more binary annotation columns derived from BED overlap.

Inputs
- `bed_files`: one or more interval annotations in BED format. Each BED
  basename becomes an output annotation name.
- `baseline_annot_dir`: directory containing chromosome-specific LDSC `.annot`
  or `.annot.gz` files that define the SNP universe and row order.
- `output_dir`: directory where the generated `.annot.gz` files will be
  written.
- `restrict_snps_path` (optional): SNP restriction file used to limit the final
  annotated SNP set before writing output. Default to None.
- `snp_identifier` (optional): how `restrict_snps_path` should be interpreted;
  supported modes are `chr_pos` (default) and `rsID`.
- `batch` (optional): if `True` (default), write one combined output per chromosome with
  one column per BED file; otherwise write one output directory per BED file.

Identifier modes
- `chr_pos`: restriction is matched by chromosome and 1-based position. If the
  restriction file is a tabular SNP list, the script infers chromosome and
  position columns. If it is a BED file, the first three BED columns are used.
- `rsID`: restriction is matched against the baseline `SNP` column. The script
  accepts either a one-column SNP list or a table with an inferred SNP column.

Workflow
1. For each baseline `.annot` file, read the template rows and keep only the
   first four LDSC columns.
2. Convert baseline SNP positions into 1 bp BED intervals for intersection.
3. Build an optional per-chromosome restriction mask on the baseline SNP
   universe.
4. Intersect the baseline SNP BED with each input BED file using
   `pybedtools.intersect(..., c=True, wa=True)`, then convert overlap counts to
   binary masks.
5. Combine BED-overlap masks with the restriction mask and write gzipped LDSC
   `.annot.gz` outputs.

Output modes
- Batch mode: write one chromosome-specific `.annot.gz` per baseline template
  with one annotation column per BED file.
- Per-BED mode: write a separate output directory per BED file, each
  containing chromosome-specific `.annot.gz` files with a single annotation
  column.

TODO: update test file.
"""

from __future__ import annotations

import argparse
import csv
import glob
import gzip
import logging
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence


LOGGER = logging.getLogger("BED2annot")
REQUIRED_ANNOT_COLUMNS = ("CHR", "BP", "SNP", "CM")
CHR_ALIASES = ("chr", "chrom", "chromosome")
BP_ALIASES = ("bp", "pos", "position", "base_pair", "basepair")
SNP_ALIASES = ("snp", "snpid", "rsid", "rs_id", "markername", "marker")


@dataclass(frozen=True)
class BaselineRow:
    # One row from a baseline LDSC template after stripping existing annotations.
    chrom: str
    bp: int
    snp: str
    cm: str


@dataclass(frozen=True)
class RestrictResource:
    # Normalized representation of the optional SNP restriction input.
    # Exactly one of `bed_path` or `snp_ids` is expected depending on `mode`.
    mode: str
    bed_path: Path | None = None
    snp_ids: frozenset[str] | None = None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    # Define the command-line interface for running the conversion as a script.
    parser = argparse.ArgumentParser(
        description="Convert BED annotations into LDSC .annot.gz files.",
    )
    parser.add_argument(
        "--bed-files",
        nargs="+",
        required=True,
        help="BED files, comma-separated lists, or glob patterns.",
    )
    parser.add_argument(
        "--baseline-annot-dir",
        required=True,
        help="Directory containing chromosome-specific baseline .annot files.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Destination directory for generated .annot.gz files.",
    )
    parser.add_argument(
        "--restrict-snps-path",
        default=None,
        help="Optional SNP restriction file matched by --snp-identifier.",
    )
    parser.add_argument(
        "--snp-identifier",
        choices=("chr_pos", "rsID"),
        default="chr_pos",
        help="How to interpret --restrict-snps-path when provided.",
    )
    parser.add_argument(
        "--no-batch",
        dest="batch",
        action="store_false",
        default=True,
        help="Write one output directory per BED file instead of combined output.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Logging verbosity.",
    )
    return parser.parse_args(argv)


def configure_logging(level: str) -> None:
    # Send progress and validation messages to stderr at the requested verbosity.
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(levelname)s: %(message)s",
    )


def get_pybedtools():
    # Import lazily so `--help` and static checks still work when pybedtools is absent.
    try:
        import pybedtools
    except ImportError as exc:  # pragma: no cover - dependency check
        raise SystemExit(
            "pybedtools is required to run BED2annot.py. Install pybedtools and "
            "bedtools, then rerun the command."
        ) from exc
    return pybedtools


def normalize_chromosome(chrom: str) -> str:
    # Normalize chromosome labels to a consistent `chrN` style for intersections.
    chrom = str(chrom).strip()
    if not chrom:
        raise ValueError("Encountered an empty chromosome label.")
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    return f"chr{chrom}"


def split_delimited_line(line: str, delimiter: str | None) -> list[str]:
    # Parse one text row while tolerating tab-, comma-, or whitespace-delimited input.
    line = line.rstrip("\n")
    if delimiter == ",":
        return next(csv.reader([line]))
    if "\t" in line:
        return line.split("\t")
    return re.split(r"\s+", line.strip())


def detect_delimiter(path: Path) -> str | None:
    # Infer the delimiter from the first non-comment data row.
    if path.name.lower().endswith(".csv") or path.name.lower().endswith(".csv.gz"):
        return ","
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "\t" in line:
                return "\t"
            if "," in line:
                return ","
            return None
    raise ValueError(f"{path} does not contain any readable data rows.")


def iter_text_lines(path: Path) -> Iterator[str]:
    # Yield non-empty, non-comment lines from plain text or gzip-compressed inputs.
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                yield line.rstrip("\n")


def infer_column_index(header: Sequence[str], aliases: Sequence[str], label: str, path: Path) -> int:
    # Match expected semantic columns case-insensitively across common aliases.
    normalized = [re.sub(r"[^a-z0-9]+", "", col.lower()) for col in header]
    alias_set = {re.sub(r"[^a-z0-9]+", "", alias.lower()) for alias in aliases}
    for idx, column_name in enumerate(normalized):
        if (
            column_name in alias_set
            or any(column_name.endswith(alias) for alias in alias_set)
            or (
                label == "position"
                and any(alias in column_name for alias in alias_set)
            )
        ):
            LOGGER.info("Using %s column '%s' from %s", label, header[idx], path)
            return idx
    raise ValueError(
        f"Could not infer a {label} column in {path}. "
        f"Header columns were: {', '.join(header)}"
    )


def expand_bed_paths(patterns: Sequence[str]) -> list[Path]:
    # Resolve BED inputs from explicit paths, comma-separated tokens, or globs.
    # Basenames must be unique because they become annotation column names.
    paths: list[Path] = []
    for item in patterns:
        for token in item.split(","):
            token = token.strip()
            if not token:
                continue
            expanded = sorted(Path(match).resolve() for match in glob.glob(token))
            if expanded:
                paths.extend(path.resolve() for path in expanded if path.is_file())
            else:
                path = Path(token).expanduser().resolve()
                if not path.is_file():
                    raise FileNotFoundError(f"BED file not found: {token}")
                paths.append(path)
    unique_paths = []
    seen = set()
    for path in paths:
        if path not in seen:
            unique_paths.append(path)
            seen.add(path)
    if not unique_paths:
        raise ValueError("No BED files were resolved from --bed-files.")
    stems = [path.stem for path in unique_paths]
    duplicates = sorted({stem for stem in stems if stems.count(stem) > 1})
    if duplicates:
        raise ValueError(
            "BED basenames must be unique because they become annotation names. "
            f"Duplicate names: {', '.join(duplicates)}"
        )
    return unique_paths


def list_baseline_annots(baseline_dir: Path) -> list[Path]:
    # Discover chromosome-specific baseline templates that define the SNP universe.
    paths = sorted(
        path
        for path in baseline_dir.iterdir()
        if path.is_file() and (path.name.endswith(".annot") or path.name.endswith(".annot.gz"))
    )
    if not paths:
        raise ValueError(f"No .annot or .annot.gz files found in {baseline_dir}")
    return paths


def read_baseline_annot(path: Path) -> list[BaselineRow]:
    # Load one baseline template and keep only the invariant LDSC columns needed
    # to rebuild custom annotations while preserving baseline row order.
    delimiter = detect_delimiter(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt") as handle:
        reader = csv.reader(handle, delimiter=delimiter) if delimiter else None
        if reader is None:
            header_line = next(handle).rstrip("\n")
            header = re.split(r"\s+", header_line.strip())
            rows_iter = handle
        else:
            header = next(reader)
            rows_iter = reader

        if tuple(header[:4]) != REQUIRED_ANNOT_COLUMNS:
            raise ValueError(
                f"{path} must begin with columns {REQUIRED_ANNOT_COLUMNS}; got {tuple(header[:4])}"
            )

        rows: list[BaselineRow] = []
        if reader is None:
            for line in rows_iter:
                fields = re.split(r"\s+", line.strip())
                if not fields or len(fields) < 4:
                    continue
                rows.append(
                    BaselineRow(
                        chrom=fields[0],
                        bp=int(fields[1]),
                        snp=fields[2],
                        cm=fields[3],
                    )
                )
        else:
            for fields in rows_iter:
                if not fields or len(fields) < 4:
                    continue
                rows.append(
                    BaselineRow(
                        chrom=fields[0],
                        bp=int(fields[1]),
                        snp=fields[2],
                        cm=fields[3],
                    )
                )
    if not rows:
        raise ValueError(f"{path} does not contain any SNP rows.")
    return rows


def write_baseline_bed(rows: Sequence[BaselineRow], path: Path) -> Path:
    # Convert 1-based SNP positions from `.annot` into 0-based 1 bp BED intervals
    # so they can be intersected against interval annotations with bedtools.
    LOGGER.info(
        "Converting baseline BP coordinates from 1-based to 0-based BED intervals for intersection."
    )
    with path.open("w") as handle:
        for row in rows:
            start = row.bp - 1
            if start < 0:
                raise ValueError(f"Invalid BP={row.bp} for SNP {row.snp} in baseline template.")
            handle.write(f"{normalize_chromosome(row.chrom)}\t{start}\t{row.bp}\t{row.snp}\n")
    return path


def write_normalized_bed(in_path: Path, out_path: Path) -> Path:
    # Create a temporary BED copy whose chromosome labels are rewritten into a
    # consistent `chrN` style expected by the internal intersection workflow.
    # This does not change interval coordinates or any trailing annotation
    # columns; it only standardizes labels such as `1` -> `chr1` so that BED
    # inputs and baseline-derived SNP BED intervals use the same chromosome
    # naming convention during pybedtools intersection.
    with (gzip.open(in_path, "rt") if in_path.suffix.lower() == ".gz" else open(in_path, "rt")) as src:
        with out_path.open("w") as dst:
            for line in src:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = split_delimited_line(line, None)
                if len(fields) < 3:
                    raise ValueError(f"BED file {in_path} has a row with fewer than three columns: {line!r}")
                chrom = normalize_chromosome(fields[0])
                dst.write("\t".join([chrom] + fields[1:]) + "\n")
    return out_path


def build_restrict_resource(
    restrict_path: Path | None,
    snp_identifier: str,
    tempdir: Path,
) -> RestrictResource | None:
    # Parse the optional restriction input once, then reuse it across all
    # baseline chromosomes and BED annotations.
    if restrict_path is None:
        return None

    if snp_identifier == "rsID":
        snp_ids = load_restrict_snp_ids(restrict_path)
        LOGGER.info("Loaded %d restriction SNP IDs from %s", len(snp_ids), restrict_path)
        return RestrictResource(mode="rsID", snp_ids=frozenset(snp_ids))

    restrict_bed_path = tempdir / "restrict_snps.normalized.bed"
    if is_bed_path(restrict_path):
        LOGGER.info("Using BED restriction file %s directly for chr_pos matching", restrict_path)
        write_normalized_bed(restrict_path, restrict_bed_path)
    else:
        write_restrict_table_as_bed(restrict_path, restrict_bed_path)
    return RestrictResource(mode="chr_pos", bed_path=restrict_bed_path)


def is_bed_path(path: Path) -> bool:
    # Use the filename extension to distinguish BED-like restriction inputs.
    name = path.name.lower()
    return name.endswith(".bed") or name.endswith(".bed.gz")


def load_restrict_snp_ids(path: Path) -> set[str]:
    # Load a restriction list keyed by SNP ID, accepting either a simple
    # one-column list or a table with an inferred SNP column.
    lines = list(iter_text_lines(path))
    if not lines:
        raise ValueError(f"Restriction file {path} is empty.")

    delimiter = detect_delimiter(path)
    first_fields = split_delimited_line(lines[0], delimiter)
    first_norm = [re.sub(r"[^a-z0-9]+", "", field.lower()) for field in first_fields]
    has_header = any(field in {alias.replace("_", "") for alias in SNP_ALIASES} for field in first_norm)

    if not has_header and len(first_fields) == 1:
        return {line.strip() for line in lines if line.strip()}

    header = first_fields
    snp_idx = infer_column_index(header, SNP_ALIASES, "SNP", path)
    snp_ids = set()
    for line in lines[1:]:
        fields = split_delimited_line(line, delimiter)
        if len(fields) <= snp_idx:
            continue
        value = fields[snp_idx].strip()
        if value:
            snp_ids.add(value)
    return snp_ids


def write_restrict_table_as_bed(in_path: Path, out_path: Path) -> Path:
    # Convert a table of 1-based chromosome/position SNPs into BED format for
    # position-based restriction via pybedtools.
    lines = list(iter_text_lines(in_path))
    if not lines:
        raise ValueError(f"Restriction file {in_path} is empty.")

    delimiter = detect_delimiter(in_path)
    header = split_delimited_line(lines[0], delimiter)
    chr_idx = infer_column_index(header, CHR_ALIASES, "chromosome", in_path)
    bp_idx = infer_column_index(header, BP_ALIASES, "position", in_path)

    LOGGER.info(
        "Converting restriction positions from 1-based to 0-based BED intervals for intersection."
    )
    with out_path.open("w") as handle:
        for line in lines[1:]:
            fields = split_delimited_line(line, delimiter)
            if len(fields) <= max(chr_idx, bp_idx):
                continue
            chrom = normalize_chromosome(fields[chr_idx])
            bp = int(fields[bp_idx])
            start = bp - 1
            if start < 0:
                raise ValueError(f"Invalid restriction BP={bp} in {in_path}")
            handle.write(f"{chrom}\t{start}\t{bp}\n")
    return out_path


def build_restrict_mask(
    baseline_rows: Sequence[BaselineRow],
    baseline_bed,
    restrict_resource: RestrictResource | None,
) -> list[bool]:
    # Project the optional restriction input onto the baseline SNP universe.
    # The result is a reusable per-row boolean mask aligned to `baseline_rows`.
    if restrict_resource is None:
        return [True] * len(baseline_rows)

    if restrict_resource.mode == "rsID":
        assert restrict_resource.snp_ids is not None
        return [row.snp in restrict_resource.snp_ids for row in baseline_rows]

    assert restrict_resource.bed_path is not None
    results = baseline_bed.intersect(str(restrict_resource.bed_path), c=True, wa=True)
    return validate_and_convert_intersection(results, baseline_rows, sanity_mode="chr_pos")


def compute_bed_overlap_mask(baseline_rows: Sequence[BaselineRow], baseline_bed, bed_path: Path) -> list[bool]:
    # Mark which baseline SNPs overlap a single BED annotation.
    results = baseline_bed.intersect(str(bed_path), c=True, wa=True)
    return validate_and_convert_intersection(results, baseline_rows, sanity_mode="rsID")


def validate_and_convert_intersection(results, baseline_rows: Sequence[BaselineRow], sanity_mode: str) -> list[bool]:
    # Convert `bedtools intersect -c` output into a boolean mask and verify that
    # pybedtools preserved baseline row identity and order as expected.
    mask: list[bool] = []
    for idx, feature in enumerate(results):
        if idx >= len(baseline_rows):
            raise ValueError("Intersection returned more rows than the baseline SNP template.")
        row = baseline_rows[idx]
        fields = feature.fields
        overlap_count = int(fields[-1])
        feature_snp = fields[3]
        if feature_snp != row.snp:
            raise ValueError(
                f"Intersection row order mismatch at index {idx}: expected SNP {row.snp}, got {feature_snp}"
            )
        if sanity_mode == "chr_pos":
            feature_chr = normalize_chromosome(fields[0])
            feature_bp = int(fields[2])
            if feature_chr != normalize_chromosome(row.chrom) or feature_bp != row.bp:
                raise ValueError(
                    "Intersection coordinate mismatch at index "
                    f"{idx}: expected ({row.chrom}, {row.bp}), got ({fields[0]}, {fields[2]})"
                )
        mask.append(overlap_count > 0)
    if len(mask) != len(baseline_rows):
        raise ValueError(
            f"Intersection returned {len(mask)} rows for a baseline template with {len(baseline_rows)} rows."
        )
    return mask


def combine_masks(a: Sequence[bool], b: Sequence[bool]) -> list[int]:
    # Final annotations are binary integers: a SNP is annotated only if it
    # overlaps the BED and also survives the optional restriction step.
    if len(a) != len(b):
        raise ValueError("Mask length mismatch while combining BED and restriction masks.")
    return [int(left and right) for left, right in zip(a, b)]


def write_annot_file(
    out_path: Path,
    rows: Sequence[BaselineRow],
    annotation_names: Sequence[str],
    masks: Sequence[Sequence[int]],
) -> None:
    # Write a chromosome-specific LDSC `.annot.gz` file while preserving the
    # baseline template columns and appending the new binary annotation columns.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as handle:
        handle.write("\t".join([*REQUIRED_ANNOT_COLUMNS, *annotation_names]) + "\n")
        for idx, row in enumerate(rows):
            annot_values = [str(mask[idx]) for mask in masks]
            handle.write(
                "\t".join([row.chrom, str(row.bp), row.snp, row.cm, *annot_values]) + "\n"
            )


def process_baseline_file(
    baseline_path: Path,
    bed_paths: Sequence[Path],
    output_dir: Path,
    batch: bool,
    restrict_resource: RestrictResource | None,
    tempdir: Path,
) -> None:
    # Process one chromosome-specific baseline template end-to-end:
    # 1. read template rows,
    # 2. build the baseline SNP BED,
    # 3. compute the reusable restriction mask,
    # 4. compute BED-overlap masks,
    # 5. write either combined or per-BED `.annot.gz` output.
    LOGGER.info("Processing baseline template %s", baseline_path.name)
    rows = read_baseline_annot(baseline_path)
    baseline_bed_path = write_baseline_bed(rows, tempdir / f"{baseline_path.name}.snps.bed")
    pybedtools = get_pybedtools()
    baseline_bed = pybedtools.BedTool(str(baseline_bed_path))

    restrict_mask = build_restrict_mask(rows, baseline_bed, restrict_resource)

    if batch:
        annotation_names = [path.stem for path in bed_paths]
        masks = []
        for bed_path in bed_paths:
            overlap_mask = compute_bed_overlap_mask(rows, baseline_bed, bed_path)
            masks.append(combine_masks(overlap_mask, restrict_mask))
        output_name = baseline_output_name(baseline_path)
        write_annot_file(output_dir / output_name, rows, annotation_names, masks)
        LOGGER.info("Wrote %s", output_dir / output_name)
    else:
        for bed_path in bed_paths:
            annotation_name = bed_path.stem
            overlap_mask = compute_bed_overlap_mask(rows, baseline_bed, bed_path)
            final_mask = combine_masks(overlap_mask, restrict_mask)
            bed_output_dir = output_dir / annotation_name
            output_name = baseline_output_name(baseline_path)
            write_annot_file(
                bed_output_dir / output_name,
                rows,
                [annotation_name],
                [final_mask],
            )
            LOGGER.info("Wrote %s", bed_output_dir / output_name)

    pybedtools.cleanup(remove_all=True)


def baseline_output_name(path: Path) -> str:
    # Standardize outputs to `.annot.gz` regardless of whether the template was
    # plain text or already compressed.
    if path.name.endswith(".annot.gz"):
        stem = path.name[: -len(".annot.gz")]
    elif path.name.endswith(".annot"):
        stem = path.name[: -len(".annot")]
    else:
        stem = path.stem
    return f"{stem}.annot.gz"


def run_bed_to_annot(
    bed_files: Sequence[str | Path],
    baseline_annot_dir: str | Path,
    output_dir: str | Path,
    restrict_snps_path: str | Path | None = None,
    snp_identifier: str = "chr_pos",
    batch: bool = True,
    log_level: str = "INFO",
) -> Path:
    # Package-first entry point:
    # 1. resolve user inputs with sensible defaults,
    # 2. normalize BED inputs into temporary files,
    # 3. build the optional restriction resource,
    # 4. process each chromosome-specific baseline template,
    # 5. return the output directory for downstream use.
    if snp_identifier not in {"chr_pos", "rsID"}:
        raise ValueError("snp_identifier must be either 'chr_pos' or 'rsID'.")

    configure_logging(log_level)

    bed_tokens = [str(path) for path in bed_files]
    bed_paths = expand_bed_paths(bed_tokens)

    baseline_dir = Path(baseline_annot_dir).expanduser().resolve()
    if not baseline_dir.is_dir():
        raise NotADirectoryError(f"Baseline annotation directory not found: {baseline_dir}")
    baseline_paths = list_baseline_annots(baseline_dir)

    output_path = Path(output_dir).expanduser().resolve()
    output_path.mkdir(parents=True, exist_ok=True)

    restrict_path = None
    if restrict_snps_path is not None:
        restrict_path = Path(restrict_snps_path).expanduser().resolve()
        if not restrict_path.is_file():
            raise FileNotFoundError(f"Restriction file not found: {restrict_path}")

    LOGGER.info("Resolved %d BED file(s)", len(bed_paths))
    LOGGER.info("Resolved %d baseline template(s)", len(baseline_paths))

    with tempfile.TemporaryDirectory(prefix="bed2annot_") as tmp:
        tempdir = Path(tmp)
        normalized_beds: list[Path] = []
        for bed_path in bed_paths:
            normalized_path = tempdir / bed_path.name
            write_normalized_bed(bed_path, normalized_path)
            normalized_beds.append(normalized_path)

        restrict_resource = build_restrict_resource(restrict_path, snp_identifier, tempdir)

        for baseline_path in baseline_paths:
            process_baseline_file(
                baseline_path=baseline_path,
                bed_paths=normalized_beds,
                output_dir=output_path,
                batch=batch,
                restrict_resource=restrict_resource,
                tempdir=tempdir,
            )

    LOGGER.info("Finished.")
    return output_path


def main(argv: Sequence[str] | None = None) -> int:
    # Thin CLI wrapper around the package-level `run_bed_to_annot` entry point.
    args = parse_args(argv)
    run_bed_to_annot(
        bed_files=args.bed_files,
        baseline_annot_dir=args.baseline_annot_dir,
        output_dir=args.output_dir,
        restrict_snps_path=args.restrict_snps_path,
        snp_identifier=args.snp_identifier,
        batch=args.batch,
        log_level=args.log_level,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        sys.exit(1)
