"""
Description
-----------
New LD-score estimation core for SNP-level annotation files with support for
either PLINK genotypes or parquet `R2` tables as the reference panel.

This module consumes SNP-level `.annot` / `.annot.gz` inputs only. It does not
perform BED-to-SNP projection; that remains the responsibility of
`utils/BED2annot.py`.

Supported Inputs
----------------
- Query and baseline annotation inputs are provided separately. Each file must
  contain SNP-level rows with required metadata columns `CHR BP SNP CM`, plus
  one or more annotation columns.
- Reference-panel input may be either:
  - PLINK genotype files (`.bed/.bim/.fam`), or
  - parquet `R2` tables containing unordered SNP pairs.
- Optional frequency metadata may be provided to populate `MAF` in outputs and
  to enable `.M_5_50` generation when parquet input is used.
- Optional regression SNP lists may be provided to define the SNP set used for
  `w_ld`. If omitted, the retained reference SNP set is used.

Identifier and Genome-Build Rules
---------------------------------
- Supported SNP identifier modes in v1 are `chr_pos` and `rsID`.
- No allele matching is attempted in v1.
- In `chr_pos` mode, parquet matching uses the selected genome build to choose
  the appropriate position columns from the parquet table. Annotation `BP`
  values must already be in the same build.
- In `rsID` mode, genome build is ignored for matching.

Outputs
-------
- Reference LD-score file: `<out>.l2.ldscore.gz`
- Reference LD-score counts: `<out>.l2.M`
- Common-SNP counts, when MAF is available: `<out>.l2.M_5_50`
- Regression-weight LD-score file: `<out>.w.l2.ldscore.gz`
- Annotation-group manifest: `<out>.annotation_groups.tsv`

Outputs retain `CM` and `MAF` when available. Missing values are written as
`NA`. Multi-chromosome output is aggregated by default; per-chromosome output is
available via a flag.

Example Usage
-------------
Python
    from ldsc_py3_Jerry.ldscore.ldscore_new import build_parser, run_ldscore

    args = build_parser().parse_args([
        "--out", "results/example",
        "--baseline-annot-chr", "data/baseline.@",
        "--query-annot", "data/query.annot.gz",
        "--r2-table-chr", "data/r2/chr@",
        "--snp-identifier", "rsID",
        "--r2-bias-mode", "unbiased",
        "--ld-wind-cm", "1",
    ])
    run_ldscore(args)

CLI
    python -m ldsc_py3_Jerry.ldsc_new \
        --out results/example \
        --baseline-annot-chr data/baseline.@ \
        --query-annot data/query.annot.gz \
        --r2-table-chr data/r2/chr@ \
        --snp-identifier rsID \
        --r2-bias-mode unbiased \
        --ld-wind-cm 1

Computation Overview
--------------------
1. Load query and baseline SNP-level annotation files chromosome by chromosome.
2. Normalize SNP identity using the selected identifier mode and require all
   annotation files for a chromosome to have identical retained SNP rows.
3. Build the canonical reference SNP table from annotation inputs.
4. Intersect annotation SNPs with SNPs actually present in the reference panel.
5. Compute LD scores per chromosome:
   - PLINK mode reuses the original genotype-block implementation.
   - parquet mode constructs a symmetric sparse `R2` matrix directly from the
     filtered pair list.
6. Compute all partitioned LD-score columns in one pass and compute `w_ld` as a
   separate one-column non-partitioned LD score over the regression SNP set.
7. Aggregate chromosome results and write LDSC-compatible outputs.

Raw-vs-Unbiased `R2` Handling
-----------------------------
- parquet input requires the user to declare whether `R2` is `raw` sample
  `r^2` or already `unbiased`.
- If `R2` is raw, off-diagonal entries are transformed using the original LDSC
  correction:

      r2_unbiased = r2_raw - (1 - r2_raw) / (n - 2)

  where `n` is the LD reference sample size.
- If `R2` is already unbiased, it is used as-is.
- The diagonal is always added internally as `r_jj^2 = 1` for every retained
  SNP.

Window Rules
------------
- `--ld-wind-snps`, `--ld-wind-kb`, and `--ld-wind-cm` follow the same semantics
  as the original codebase.
- In parquet mode, window definitions are used to filter pair rows directly;
  the code does not build PLINK-style moving genotype chunks.
- `--ld-wind-cm` requires non-missing CM coordinates for retained SNPs.

MAF and `.M_5_50` Rules
-----------------------
- In PLINK mode, MAF is taken from the genotype data after SNP filtering.
- In parquet mode, MAF is only available when supplied via an optional
  frequency/metadata file.
- If MAF is unavailable, `MAF` is written as `NA` in LD-score outputs and
  `.M_5_50` is not emitted.

Current Limitations
-------------------
- v1 supports only `chr_pos` and `rsID` matching.
- v1 does not attempt allele-aware reconciliation between annotations and
  parquet pairs.
- parquet input requires `pyarrow` at runtime.
- Per-query partitioned LDSC wrappers are intentionally out of scope for this
  module; this module computes LD scores only.

Dependencies
------------
- `numpy`
- `pandas`
- `scipy`
- `pyarrow` for parquet mode
"""

from __future__ import annotations

import argparse
import glob
import gzip
import logging
import os
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import ldscore.parse as legacy_parse
except ImportError:  # pragma: no cover - package import fallback
    from . import parse as legacy_parse


LOGGER = logging.getLogger("LDSC.new")
REQUIRED_ANNOT_COLUMNS = ("CHR", "BP", "SNP", "CM")
ANNOT_META_COLUMNS = ("CHR", "SNP", "BP", "CM", "MAF")
CHROM_ALIASES = ("CHR", "CHROM", "CHROMOSOME")
BP_ALIASES = ("BP", "POS", "POSITION")
SNP_ALIASES = ("SNP", "RSID", "RS_ID", "MARKER", "MARKERNAME")
CM_ALIASES = ("CM", "CMBP", "CENTIMORGAN")
MAF_ALIASES = ("MAF", "FRQ", "FREQ", "FREQUENCY")


@dataclass
class AnnotationBundle:
    metadata: pd.DataFrame
    annotations: pd.DataFrame
    baseline_columns: list[str]
    query_columns: list[str]


@dataclass
class ChromComputationResult:
    chrom: str
    metadata: pd.DataFrame
    ld_scores: np.ndarray
    w_ld: np.ndarray
    M: np.ndarray
    M_5_50: np.ndarray | None
    ldscore_columns: list[str]
    baseline_columns: list[str]
    query_columns: list[str]


# Basic configuration and shared helpers.
def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(levelname)s: %(message)s",
    )


def get_legacy_ld_module():
    try:
        import ldscore.ldscore as legacy_ld
    except ImportError as exc:
        try:
            from . import ldscore as legacy_ld
        except ImportError:
            raise ImportError(
                "PLINK reference-panel mode requires the legacy ldscore backend and its dependencies, "
                "including bitarray."
            ) from exc
    return legacy_ld


def split_arg_list(value: str | None) -> list[str]:
    if not value:
        return []
    out = []
    for token in value.split(","):
        token = os.path.expanduser(os.path.expandvars(token.strip()))
        if token:
            out.append(token)
    return out


def normalize_chromosome(value: object) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError("Encountered an empty chromosome label.")
    if text.lower().startswith("chr"):
        text = text[3:]
    return text.upper()


def chrom_sort_key(chrom: object) -> tuple[int, object]:
    normalized = normalize_chromosome(chrom)
    try:
        return (0, int(normalized))
    except ValueError:
        special = {"X": 23, "Y": 24, "MT": 25, "M": 25}
        return (1, special.get(normalized, normalized))


def find_column(columns: Iterable[str], aliases: Sequence[str]) -> str | None:
    normalized = {str(col).strip().upper(): col for col in columns}
    for alias in aliases:
        if alias in normalized:
            return normalized[alias]
    return None


def get_block_lefts(coords: np.ndarray, max_dist: float) -> np.ndarray:
    M = len(coords)
    j = 0
    block_left = np.zeros(M)
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1
        block_left[i] = j
    return block_left


def identifier_keys(df: pd.DataFrame, mode: str) -> pd.Series:
    if mode == "rsID":
        return df["SNP"].astype(str)
    if mode != "chr_pos":
        raise ValueError(f"Unsupported identifier mode: {mode}")
    chrom = df["CHR"].map(normalize_chromosome)
    bp = df["BP"].astype(np.int64).astype(str)
    return chrom + ":" + bp


def sort_frame_by_genomic_position(df: pd.DataFrame) -> pd.DataFrame:
    sort_df = df.copy()
    sort_df["_chrom_key"] = sort_df["CHR"].map(chrom_sort_key)
    sort_df = sort_df.sort_values(by=["_chrom_key", "BP", "SNP"], kind="mergesort")
    return sort_df.drop(columns="_chrom_key").reset_index(drop=True)


def resolve_prefixed_file(token: str, suffixes: Sequence[str]) -> list[str]:
    if os.path.exists(token):
        return [token]

    matches = sorted({path for path in glob.glob(token)})
    if matches:
        return matches

    for suffix in suffixes:
        candidate = token + suffix
        if os.path.exists(candidate):
            return [candidate]

    raise FileNotFoundError(f"Could not resolve input path or prefix: {token}")


def resolve_annotation_files(spec: str | None) -> list[str]:
    files: list[str] = []
    for token in split_arg_list(spec):
        files.extend(resolve_prefixed_file(token, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz")))
    return files


def resolve_chr_files(spec: str | None, chrom: str, suffixes: Sequence[str]) -> list[str]:
    files: list[str] = []
    for token in split_arg_list(spec):
        chrom_token = legacy_parse.sub_chr(token, chrom)
        files.extend(resolve_prefixed_file(chrom_token, suffixes))
    return files


def resolve_optional_chr_files(spec: str | None, chrom: str, suffixes: Sequence[str]) -> list[str]:
    files: list[str] = []
    for token in split_arg_list(spec):
        chrom_token = legacy_parse.sub_chr(token, chrom)
        try:
            files.extend(resolve_prefixed_file(chrom_token, suffixes))
        except FileNotFoundError:
            continue
    return files


def resolve_parquet_files(args: argparse.Namespace, chrom: str | None = None) -> list[str]:
    if chrom is not None and args.r2_table_chr:
        return resolve_chr_files(args.r2_table_chr, chrom, ("", ".parquet"))
    if args.r2_table:
        files = []
        for token in split_arg_list(args.r2_table):
            files.extend(resolve_prefixed_file(token, ("", ".parquet")))
        return files
    return []


def resolve_bfile_prefix(args: argparse.Namespace, chrom: str | None = None) -> str | None:
    prefix = None
    if chrom is not None and args.bfile_chr:
        prefix = legacy_parse.sub_chr(args.bfile_chr, chrom)
    elif args.bfile:
        prefix = args.bfile
    if prefix is None:
        return None
    if prefix.endswith(".bed") and os.path.exists(prefix):
        prefix = prefix[:-4]
    if os.path.exists(prefix + ".bed"):
        return prefix
    raise FileNotFoundError(f"Could not resolve PLINK prefix: {prefix}")


def resolve_frequency_files(args: argparse.Namespace, chrom: str | None = None) -> list[str]:
    files: list[str] = []
    if chrom is not None and args.frqfile_chr:
        return resolve_optional_chr_files(args.frqfile_chr, chrom, ("", ".frq", ".frq.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
    for token in split_arg_list(args.frqfile):
        files.extend(resolve_prefixed_file(token, ("", ".frq", ".frq.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz")))
    return files


def read_text_table(path: str) -> pd.DataFrame:
    compression = "gzip" if path.endswith(".gz") else None
    return pd.read_csv(path, sep=r"\s+", compression=compression)


# Annotation loading and normalization.
def parse_annotation_file(path: str, chrom: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_text_table(path)
    missing = [col for col in REQUIRED_ANNOT_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required annotation columns: {missing}")

    meta = df.loc[:, list(REQUIRED_ANNOT_COLUMNS)].copy()
    meta["CHR"] = meta["CHR"].map(normalize_chromosome)
    meta["BP"] = pd.to_numeric(meta["BP"], errors="raise").astype(np.int64)
    meta["SNP"] = meta["SNP"].astype(str)
    meta["CM"] = pd.to_numeric(meta["CM"], errors="coerce")
    if "MAF" in df.columns:
        meta["MAF"] = pd.to_numeric(df["MAF"], errors="coerce")

    if chrom is not None:
        keep = meta["CHR"] == normalize_chromosome(chrom)
        meta = meta.loc[keep].reset_index(drop=True)
        df = df.loc[keep].reset_index(drop=True)
    else:
        meta = meta.reset_index(drop=True)

    if len(meta) == 0:
        return meta, pd.DataFrame(index=meta.index)

    meta = sort_frame_by_genomic_position(meta)
    annotation_columns = [col for col in df.columns if col not in list(REQUIRED_ANNOT_COLUMNS) + ["MAF"]]
    if not annotation_columns:
        raise ValueError(f"{path} does not contain any annotation columns.")

    sorted_df = df.copy()
    sorted_df["CHR"] = sorted_df["CHR"].map(normalize_chromosome)
    sorted_df["BP"] = pd.to_numeric(sorted_df["BP"], errors="raise").astype(np.int64)
    sorted_df["SNP"] = sorted_df["SNP"].astype(str)
    sorted_df["_chrom_key"] = sorted_df["CHR"].map(chrom_sort_key)
    sorted_df = sorted_df.sort_values(by=["_chrom_key", "BP", "SNP"], kind="mergesort").drop(columns="_chrom_key")
    annotations = sorted_df.loc[:, annotation_columns].astype(np.float32).reset_index(drop=True)
    meta = meta.reset_index(drop=True)
    if len(meta) != len(annotations):
        raise ValueError(f"{path} metadata and annotation lengths diverged after sorting.")

    return meta, annotations


def combine_annotation_groups(
    baseline_files: Sequence[str],
    query_files: Sequence[str],
    chrom: str,
    identifier_mode: str,
) -> AnnotationBundle | None:
    """
    Read and column-bind baseline/query annotation files for one chromosome.

    Main steps:
    1. Parse each SNP-level annotation file.
    2. Normalize the chosen SNP identifier key.
    3. Require identical SNP rows across all files for that chromosome.
    4. Merge annotation columns into one dense matrix while preserving the
       baseline/query grouping.
    """
    frames: list[pd.DataFrame] = []
    blocks: list[pd.DataFrame] = []
    baseline_columns: list[str] = []
    query_columns: list[str] = []
    seen_columns: set[str] = set()

    for group_name, files in (("baseline", baseline_files), ("query", query_files)):
        for path in files:
            meta, annotations = parse_annotation_file(path, chrom=chrom)
            if len(meta) == 0:
                continue
            meta = meta.copy()
            meta["_key"] = identifier_keys(meta, identifier_mode)
            if not frames:
                frames.append(meta)
            else:
                reference = frames[0]
                if len(reference) != len(meta) or not reference["_key"].equals(meta["_key"]):
                    raise ValueError(
                        f"Annotation SNP rows do not match across files for chromosome {chrom}: {path}"
                    )
                same_core = reference.loc[:, ["CHR", "BP", "SNP"]].equals(meta.loc[:, ["CHR", "BP", "SNP"]])
                if not same_core:
                    raise ValueError(
                        f"Annotation metadata mismatch across files for chromosome {chrom}: {path}"
                    )
                missing_cm = reference["CM"].isna() & meta["CM"].notna()
                if missing_cm.any():
                    frames[0].loc[missing_cm, "CM"] = meta.loc[missing_cm, "CM"].to_numpy()
                if "MAF" in meta.columns:
                    if "MAF" not in frames[0].columns:
                        frames[0]["MAF"] = np.nan
                    missing_maf = frames[0]["MAF"].isna() & meta["MAF"].notna()
                    if missing_maf.any():
                        frames[0].loc[missing_maf, "MAF"] = meta.loc[missing_maf, "MAF"].to_numpy()

            for column in annotations.columns:
                if column in seen_columns:
                    raise ValueError(f"Duplicate annotation column name detected: {column}")
                seen_columns.add(column)

            blocks.append(annotations.reset_index(drop=True))
            if group_name == "baseline":
                baseline_columns.extend(annotations.columns.tolist())
            else:
                query_columns.extend(annotations.columns.tolist())

    if not frames:
        return None

    metadata = frames[0].drop(columns="_key").reset_index(drop=True)
    annotations = pd.concat(blocks, axis=1)
    ordered_columns = baseline_columns + query_columns
    annotations = annotations.loc[:, ordered_columns].reset_index(drop=True)
    return AnnotationBundle(
        metadata=metadata,
        annotations=annotations,
        baseline_columns=baseline_columns,
        query_columns=query_columns,
    )


def read_identifier_list(path: str, mode: str) -> set[str]:
    df = read_text_table(path)
    if df.shape[1] == 1:
        col = df.columns[0]
        values = df[col].astype(str).str.strip()
        if mode == "chr_pos":
            if values.str.contains(":").all():
                return set(values.map(lambda x: normalize_chromosome(x.split(":")[0]) + ":" + str(int(x.split(":")[1]))))
            raise ValueError(f"{path} must contain CHR:BP values in chr_pos mode when only one column is supplied.")
        return set(values)

    snp_col = find_column(df.columns, SNP_ALIASES)
    chr_col = find_column(df.columns, CHROM_ALIASES)
    bp_col = find_column(df.columns, BP_ALIASES)
    if mode == "rsID":
        if snp_col is None:
            raise ValueError(f"{path} must contain a SNP column in rsID mode.")
        return set(df[snp_col].astype(str))
    if chr_col is None or bp_col is None:
        raise ValueError(f"{path} must contain CHR and BP columns in chr_pos mode.")
    chrom = df[chr_col].map(normalize_chromosome)
    bp = pd.to_numeric(df[bp_col], errors="raise").astype(np.int64).astype(str)
    return set(chrom + ":" + bp)


def load_regression_keys(args: argparse.Namespace) -> set[str] | None:
    if not args.regression_snps:
        return None
    return read_identifier_list(args.regression_snps, args.snp_identifier)


def parse_frequency_metadata(path: str, chrom: str | None, identifier_mode: str) -> pd.DataFrame:
    df = read_text_table(path)
    chr_col = find_column(df.columns, CHROM_ALIASES)
    bp_col = find_column(df.columns, BP_ALIASES)
    snp_col = find_column(df.columns, SNP_ALIASES)
    cm_col = find_column(df.columns, CM_ALIASES)
    maf_col = find_column(df.columns, MAF_ALIASES)

    if chrom is not None and chr_col is not None:
        keep = df[chr_col].map(normalize_chromosome) == normalize_chromosome(chrom)
        df = df.loc[keep].reset_index(drop=True)

    if len(df) == 0:
        return pd.DataFrame(columns=["_key", "CM", "MAF"])

    out = pd.DataFrame(index=df.index)
    if identifier_mode == "rsID":
        if snp_col is None:
            raise ValueError(f"{path} must contain a SNP column in rsID mode.")
        out["_key"] = df[snp_col].astype(str)
    else:
        if chr_col is None or bp_col is None:
            raise ValueError(f"{path} must contain CHR and BP columns in chr_pos mode.")
        chrom_norm = df[chr_col].map(normalize_chromosome)
        bp = pd.to_numeric(df[bp_col], errors="raise").astype(np.int64).astype(str)
        out["_key"] = chrom_norm + ":" + bp

    if cm_col is not None:
        out["CM"] = pd.to_numeric(df[cm_col], errors="coerce")
    if maf_col is not None:
        maf = pd.to_numeric(df[maf_col], errors="coerce").astype(float)
        out["MAF"] = np.minimum(maf, 1.0 - maf)

    return out


def merge_frequency_metadata(
    metadata: pd.DataFrame,
    args: argparse.Namespace,
    chrom: str,
    identifier_mode: str,
) -> pd.DataFrame:
    files = resolve_frequency_files(args, chrom=chrom)
    if not files:
        return metadata

    frames = [parse_frequency_metadata(path, chrom=chrom, identifier_mode=identifier_mode) for path in files]
    freq_df = pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame(columns=["_key", "CM", "MAF"])
    if freq_df["_key"].duplicated().any():
        raise ValueError(f"Duplicate SNP metadata keys detected in frequency metadata for chromosome {chrom}.")

    merged = metadata.copy()
    merged["_key"] = identifier_keys(merged, identifier_mode)
    freq_df = freq_df.set_index("_key")
    if "CM" in freq_df.columns:
        if "CM" not in merged.columns:
            merged["CM"] = np.nan
        missing_cm = merged["CM"].isna() & merged["_key"].isin(freq_df.index)
        merged.loc[missing_cm, "CM"] = freq_df.loc[merged.loc[missing_cm, "_key"], "CM"].to_numpy()
    if "MAF" in freq_df.columns:
        if "MAF" not in merged.columns:
            merged["MAF"] = np.nan
        missing_maf = merged["MAF"].isna() & merged["_key"].isin(freq_df.index)
        merged.loc[missing_maf, "MAF"] = freq_df.loc[merged.loc[missing_maf, "_key"], "MAF"].to_numpy()
    return merged.drop(columns="_key")


def apply_maf_filter(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    maf_min: float | None,
    context: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if maf_min is None:
        return metadata, annotations
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        LOGGER.warning("Cannot apply --maf in %s because MAF metadata is unavailable.", context)
        return metadata, annotations
    keep = metadata["MAF"] > maf_min
    removed = int((~keep).sum())
    if removed:
        LOGGER.info("Removed %d SNPs with MAF <= %s in %s.", removed, maf_min, context)
    metadata = metadata.loc[keep].reset_index(drop=True)
    annotations = annotations.loc[keep].reset_index(drop=True)
    return metadata, annotations


def chromosome_set_from_annotation_inputs(args: argparse.Namespace) -> list[str]:
    chromosomes: set[str] = set()
    all_files = (
        resolve_annotation_files(args.query_annot)
        + resolve_annotation_files(args.baseline_annot)
    )
    for path in all_files:
        df = read_text_table(path)
        if "CHR" not in df.columns:
            raise ValueError(f"{path} is missing CHR.")
        chromosomes.update(df["CHR"].map(normalize_chromosome).unique().tolist())

    chr_specs = split_arg_list(args.query_annot_chr) + split_arg_list(args.baseline_annot_chr)
    if chr_specs:
        for chrom in [str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"]:
            files = resolve_optional_chr_files(args.query_annot_chr, chrom, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
            files += resolve_optional_chr_files(args.baseline_annot_chr, chrom, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
            if files:
                chromosomes.add(normalize_chromosome(chrom))

    if not chromosomes:
        raise ValueError("No annotation chromosomes could be resolved from the supplied inputs.")
    return sorted(chromosomes, key=chrom_sort_key)


def build_window_coordinates(metadata: pd.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray, float]:
    selectors = np.array([args.ld_wind_snps is not None, args.ld_wind_kb is not None, args.ld_wind_cm is not None], dtype=bool)
    if selectors.sum() != 1:
        raise ValueError("Must specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.")

    if args.ld_wind_snps is not None:
        return np.arange(len(metadata), dtype=float), float(args.ld_wind_snps)
    if args.ld_wind_kb is not None:
        return metadata["BP"].to_numpy(dtype=float), float(args.ld_wind_kb) * 1000.0
    if metadata["CM"].isna().any():
        raise ValueError("--ld-wind-cm requires non-missing CM values for all retained SNPs.")
    return metadata["CM"].to_numpy(dtype=float), float(args.ld_wind_cm)


def check_whole_chromosome_window(block_left: np.ndarray, args: argparse.Namespace, chrom: str) -> None:
    if len(block_left) == 0:
        return
    if block_left[-1] == 0 and not args.yes_really:
        raise ValueError(
            f"Chromosome {chrom} would use a whole-chromosome LD window. Re-run with --yes-really to allow this."
        )


def get_pyarrow_dataset():
    try:
        import pyarrow.dataset as ds
    except ImportError as exc:
        raise ImportError(
            "pyarrow is required for parquet R2 input. Install pyarrow and retry."
        ) from exc
    return ds


# Parquet R2 adapter.
def read_parquet_pairs(args: argparse.Namespace, chrom: str) -> pd.DataFrame:
    """Read one chromosome of parquet R2 pairs and derive the active SNP keys."""
    files = resolve_parquet_files(args, chrom=chrom)
    if not files:
        raise FileNotFoundError(f"No parquet R2 files resolved for chromosome {chrom}.")

    ds = get_pyarrow_dataset()
    dataset = ds.dataset(files, format="parquet")
    columns = ["chr", "R2"]
    if args.snp_identifier == "rsID":
        columns.extend(["rsID_1", "rsID_2"])
    else:
        if args.genome_build == "hg19":
            columns.extend(["hg19_bp_1", "hg19_bp_2"])
        else:
            columns.extend(["hg38_bp1", "hg38_bp2"])
    table = dataset.to_table(columns=columns)
    pairs = table.to_pandas()
    pairs["chr"] = pairs["chr"].map(normalize_chromosome)
    pairs = pairs.loc[pairs["chr"] == normalize_chromosome(chrom)].reset_index(drop=True)
    if len(pairs) == 0:
        return pairs

    pairs["R2"] = pd.to_numeric(pairs["R2"], errors="raise").astype(np.float32)
    if args.snp_identifier == "rsID":
        pairs["key1"] = pairs["rsID_1"].astype(str)
        pairs["key2"] = pairs["rsID_2"].astype(str)
    else:
        bp1 = "hg19_bp_1" if args.genome_build == "hg19" else "hg38_bp1"
        bp2 = "hg19_bp_2" if args.genome_build == "hg19" else "hg38_bp2"
        pairs["key1"] = pairs["chr"] + ":" + pd.to_numeric(pairs[bp1], errors="raise").astype(np.int64).astype(str)
        pairs["key2"] = pairs["chr"] + ":" + pd.to_numeric(pairs[bp2], errors="raise").astype(np.int64).astype(str)

    return pairs.loc[:, ["key1", "key2", "R2"]]


def filter_reference_to_present_pairs(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    regression_keys: set[str] | None,
    pairs: pd.DataFrame,
    identifier_mode: str,
    chrom: str,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str] | None]:
    metadata = metadata.copy()
    metadata["_key"] = identifier_keys(metadata, identifier_mode)
    present = set(pairs["key1"]).union(set(pairs["key2"]))
    keep = metadata["_key"].isin(present)
    removed = int((~keep).sum())
    if removed:
        LOGGER.warning(
            "Dropping %d annotated SNPs on chromosome %s because they are absent from the parquet R2 table.",
            removed,
            chrom,
        )
    metadata = metadata.loc[keep].reset_index(drop=True)
    annotations = annotations.loc[keep].reset_index(drop=True)
    if regression_keys is not None:
        regression_keys = regression_keys.intersection(set(metadata["_key"]))
    return metadata, annotations, regression_keys


def parquet_to_sparse_r2(
    pairs: pd.DataFrame,
    metadata: pd.DataFrame,
    args: argparse.Namespace,
) -> sparse.csr_matrix:
    """
    Convert filtered parquet pairs into a symmetric sparse R2 matrix.

    Main steps:
    1. Map SNP identifiers onto reference-table row indices.
    2. Drop non-reference pairs and self-pairs from the parquet rows.
    3. Enforce the requested LD window directly on pair coordinates.
    4. Apply raw-to-unbiased correction when requested.
    5. Add the diagonal and symmetrize the matrix.
    """
    if len(metadata) == 0:
        return sparse.csr_matrix((0, 0), dtype=np.float32)

    metadata = metadata.copy()
    keys = metadata["_key"].tolist()
    key_to_index = {key: idx for idx, key in enumerate(keys)}

    pair_df = pairs.copy()
    pair_df["i"] = pair_df["key1"].map(key_to_index)
    pair_df["j"] = pair_df["key2"].map(key_to_index)
    pair_df = pair_df.dropna(subset=["i", "j"]).copy()
    pair_df["i"] = pair_df["i"].astype(np.int64)
    pair_df["j"] = pair_df["j"].astype(np.int64)
    pair_df = pair_df.loc[pair_df["i"] != pair_df["j"]].reset_index(drop=True)
    pair_df["lo"] = np.minimum(pair_df["i"], pair_df["j"])
    pair_df["hi"] = np.maximum(pair_df["i"], pair_df["j"])

    if len(pair_df) == 0:
        diag = sparse.identity(len(metadata), dtype=np.float32, format="csr")
        return diag

    coords, max_dist = build_window_coordinates(metadata.drop(columns="_key"), args)
    i = pair_df["i"].to_numpy(dtype=np.int64)
    j = pair_df["j"].to_numpy(dtype=np.int64)
    if args.ld_wind_snps is not None:
        keep = np.abs(i - j) <= max_dist
    else:
        keep = np.abs(coords[i] - coords[j]) <= max_dist
    pair_df = pair_df.loc[keep].reset_index(drop=True)

    if pair_df.duplicated(subset=["lo", "hi"]).any():
        raise ValueError("Duplicate unordered SNP pairs detected in parquet input after key mapping.")

    data = pair_df["R2"].to_numpy(dtype=np.float32)
    if args.r2_bias_mode == "raw":
        if args.r2_sample_size is None:
            raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
        denom = args.r2_sample_size - 2
        if denom <= 0:
            raise ValueError("--r2-sample-size must be greater than 2 for raw R2 correction.")
        data = data - (1.0 - data) / denom

    rows = pair_df["lo"].to_numpy(dtype=np.int64)
    cols = pair_df["hi"].to_numpy(dtype=np.int64)
    diag_idx = np.arange(len(metadata), dtype=np.int64)
    full_rows = np.concatenate([rows, cols, diag_idx])
    full_cols = np.concatenate([cols, rows, diag_idx])
    full_data = np.concatenate([data, data, np.ones(len(metadata), dtype=np.float32)])
    matrix = sparse.coo_matrix((full_data, (full_rows, full_cols)), shape=(len(metadata), len(metadata)), dtype=np.float32)
    return matrix.tocsr()


def compute_ld_scores_from_sparse(
    r2_matrix: sparse.csr_matrix,
    annotations: pd.DataFrame,
    regression_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=True)
    ld_scores = r2_matrix.dot(annot_matrix)
    w_ld = r2_matrix.dot(regression_mask.astype(np.float32).reshape(-1, 1))
    return np.asarray(ld_scores, dtype=np.float32), np.asarray(w_ld, dtype=np.float32)


def compute_counts(metadata: pd.DataFrame, annotations: pd.DataFrame) -> tuple[np.ndarray, np.ndarray | None]:
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=False)
    M = np.asarray(annot_matrix.sum(axis=0), dtype=np.float64)
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        return M, None
    common = metadata["MAF"] > 0.05
    M_5_50 = np.asarray(annot_matrix[common.to_numpy(), :].sum(axis=0), dtype=np.float64)
    return M, M_5_50


def regression_mask_from_keys(metadata: pd.DataFrame, regression_keys: set[str] | None, identifier_mode: str) -> np.ndarray:
    if regression_keys is None:
        return np.ones(len(metadata), dtype=np.float32)
    keys = identifier_keys(metadata, identifier_mode)
    return keys.isin(regression_keys).to_numpy(dtype=np.float32)


# Per-chromosome compute backends.
def compute_chrom_from_parquet(
    chrom: str,
    bundle: AnnotationBundle,
    args: argparse.Namespace,
    regression_keys: set[str] | None,
) -> ChromComputationResult:
    """
    Compute all LD-score outputs for one chromosome from parquet R2 input.

    Main steps:
    1. Merge optional CM/MAF metadata and apply any MAF filter.
    2. Intersect annotation SNPs with SNPs actually present in the R2 table.
    3. Validate the LD window against the retained SNP coordinates.
    4. Build the sparse R2 matrix and multiply by annotation/regression vectors.
    5. Return chromosome-level LD scores plus M and M_5_50 counts.
    """
    metadata = merge_frequency_metadata(bundle.metadata.copy(), args, chrom=chrom, identifier_mode=args.snp_identifier)
    metadata, annotations = apply_maf_filter(metadata, bundle.annotations.copy(), args.maf, context="parquet mode")
    pairs = read_parquet_pairs(args, chrom=chrom)
    metadata, annotations, regression_keys = filter_reference_to_present_pairs(
        metadata, annotations, regression_keys, pairs, args.snp_identifier, chrom
    )
    if len(metadata) == 0:
        raise ValueError(f"No retained annotation SNPs remain on chromosome {chrom} after parquet intersection.")

    metadata["_key"] = identifier_keys(metadata, args.snp_identifier)
    coords, _ = build_window_coordinates(metadata.drop(columns="_key"), args)
    block_left = get_block_lefts(
        coords,
        args.ld_wind_snps if args.ld_wind_snps is not None else (args.ld_wind_kb * 1000 if args.ld_wind_kb is not None else args.ld_wind_cm),
    )
    check_whole_chromosome_window(block_left, args, chrom)
    regression_mask = regression_mask_from_keys(metadata.drop(columns="_key"), regression_keys, args.snp_identifier)
    r2_matrix = parquet_to_sparse_r2(pairs, metadata, args)
    ld_scores, w_ld = compute_ld_scores_from_sparse(r2_matrix, annotations, regression_mask)
    out_metadata = metadata.drop(columns="_key").reset_index(drop=True)
    if "MAF" not in out_metadata.columns:
        out_metadata["MAF"] = np.nan
    M, M_5_50 = compute_counts(out_metadata, annotations)
    return ChromComputationResult(
        chrom=chrom,
        metadata=out_metadata,
        ld_scores=ld_scores,
        w_ld=w_ld,
        M=M,
        M_5_50=M_5_50,
        ldscore_columns=bundle.baseline_columns + bundle.query_columns,
        baseline_columns=bundle.baseline_columns,
        query_columns=bundle.query_columns,
    )


def compute_chrom_from_plink(
    chrom: str,
    bundle: AnnotationBundle,
    args: argparse.Namespace,
    regression_keys: set[str] | None,
) -> ChromComputationResult:
    """
    Compute all LD-score outputs for one chromosome from a PLINK reference panel.

    Main steps:
    1. Align annotation SNPs to the PLINK BIM table.
    2. Reuse the legacy PLINK genotype reader and LD-score kernel.
    3. Compute partitioned reference LD scores and one-column regression weights.
    4. Return chromosome-level LD scores plus M and M_5_50 counts.
    """
    legacy_ld = get_legacy_ld_module()
    prefix = resolve_bfile_prefix(args, chrom=chrom)
    if prefix is None:
        raise ValueError("PLINK mode requested without --bfile or --bfile-chr.")

    bim = legacy_parse.PlinkBIMFile(prefix + ".bim")
    fam = legacy_parse.PlinkFAMFile(prefix + ".fam")
    panel_df = bim.df.loc[:, ["CHR", "SNP", "CM", "BP"]].copy()
    panel_df["CHR"] = panel_df["CHR"].map(normalize_chromosome)
    panel_df["SNP"] = panel_df["SNP"].astype(str)
    panel_df["BP"] = pd.to_numeric(panel_df["BP"], errors="raise").astype(np.int64)
    panel_df["CM"] = pd.to_numeric(panel_df["CM"], errors="coerce")
    panel_df = panel_df.loc[panel_df["CHR"] == normalize_chromosome(chrom)].copy()
    if len(panel_df) == 0:
        raise ValueError(f"No PLINK SNPs found for chromosome {chrom} in {prefix}.")
    panel_df["_key"] = identifier_keys(panel_df, args.snp_identifier)

    metadata = bundle.metadata.copy()
    metadata = merge_frequency_metadata(metadata, args, chrom=chrom, identifier_mode=args.snp_identifier)
    metadata, annotations = apply_maf_filter(metadata, bundle.annotations.copy(), args.maf, context="PLINK mode")
    metadata["_key"] = identifier_keys(metadata, args.snp_identifier)

    key_to_panel_index = {key: idx for key, idx in zip(panel_df["_key"], panel_df.index)}
    keep = metadata["_key"].isin(key_to_panel_index)
    removed = int((~keep).sum())
    if removed:
        LOGGER.warning(
            "Dropping %d annotated SNPs on chromosome %s because they are absent from the PLINK reference panel.",
            removed,
            chrom,
        )
    metadata = metadata.loc[keep].reset_index(drop=True)
    annotations = annotations.loc[keep].reset_index(drop=True)
    if regression_keys is not None:
        regression_keys = regression_keys.intersection(set(metadata["_key"]))
    if len(metadata) == 0:
        raise ValueError(f"No retained annotation SNPs remain on chromosome {chrom} after PLINK intersection.")

    keep_snps = [key_to_panel_index[key] for key in metadata["_key"]]
    geno = legacy_ld.PlinkBEDFile(prefix + ".bed", len(fam.IDList), bim, keep_snps=keep_snps, mafMin=args.maf)

    geno_meta = pd.DataFrame(geno.df, columns=geno.colnames)
    geno_meta["CHR"] = geno_meta["CHR"].map(normalize_chromosome)
    geno_meta["SNP"] = geno_meta["SNP"].astype(str)
    geno_meta["BP"] = pd.to_numeric(geno_meta["BP"], errors="raise").astype(np.int64)
    geno_meta["CM"] = pd.to_numeric(geno_meta["CM"], errors="coerce")
    geno_meta["MAF"] = pd.to_numeric(geno_meta["MAF"], errors="coerce")
    geno_meta["_key"] = identifier_keys(geno_meta, args.snp_identifier)

    annotation_matrix = annotations.set_index(metadata["_key"]).loc[geno_meta["_key"]]

    coords, max_dist = build_window_coordinates(geno_meta.drop(columns="_key"), args)
    block_left = legacy_ld.getBlockLefts(coords, max_dist)
    check_whole_chromosome_window(block_left, args, chrom)

    ld_scores = geno.ldScoreVarBlocks(block_left, args.chunk_size, annot=annotation_matrix.to_numpy(dtype=np.float32))
    regression_mask = regression_mask_from_keys(geno_meta.drop(columns="_key"), regression_keys, args.snp_identifier)
    w_ld = geno.ldScoreVarBlocks(block_left, args.chunk_size, annot=regression_mask.reshape(-1, 1))
    out_metadata = geno_meta.drop(columns="_key").reset_index(drop=True)
    M, M_5_50 = compute_counts(out_metadata, annotation_matrix.reset_index(drop=True))
    return ChromComputationResult(
        chrom=chrom,
        metadata=out_metadata,
        ld_scores=np.asarray(ld_scores, dtype=np.float32),
        w_ld=np.asarray(w_ld, dtype=np.float32),
        M=M,
        M_5_50=M_5_50,
        ldscore_columns=bundle.baseline_columns + bundle.query_columns,
        baseline_columns=bundle.baseline_columns,
        query_columns=bundle.query_columns,
    )


def result_to_dataframe(result: ChromComputationResult) -> pd.DataFrame:
    df = result.metadata.copy()
    for idx, column in enumerate(result.ldscore_columns):
        df[column + "L2"] = result.ld_scores[:, idx]
    return df


def weight_result_to_dataframe(result: ChromComputationResult) -> pd.DataFrame:
    df = result.metadata.copy()
    df["L2"] = np.ravel(result.w_ld)
    return df


# Output assembly.
def write_ldscore_file(df: pd.DataFrame, path: str) -> None:
    out = df.copy()
    out = out.loc[:, [col for col in ["CHR", "SNP", "BP", "CM", "MAF"] if col in out.columns] + [col for col in out.columns if col not in ANNOT_META_COLUMNS]]
    with gzip.open(path, "wt") as handle:
        out.to_csv(handle, sep="\t", index=False, na_rep="NA", float_format="%.6g")


def write_counts(path: str, counts: np.ndarray) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(str(x) for x in counts))


def write_annotation_groups(path: str, baseline_columns: Sequence[str], query_columns: Sequence[str]) -> None:
    rows = [{"annotation": col, "group": "baseline"} for col in baseline_columns]
    rows.extend({"annotation": col, "group": "query"} for col in query_columns)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def aggregate_results(results: Sequence[ChromComputationResult]) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray | None]:
    ld_frames = [result_to_dataframe(result) for result in results]
    weight_frames = [weight_result_to_dataframe(result) for result in results]
    ld_df = pd.concat(ld_frames, axis=0, ignore_index=True)
    weight_df = pd.concat(weight_frames, axis=0, ignore_index=True)
    ld_df = sort_frame_by_genomic_position(ld_df)
    weight_df = sort_frame_by_genomic_position(weight_df)
    M = np.sum(np.vstack([result.M for result in results]), axis=0)
    all_have_m_5_50 = all(result.M_5_50 is not None for result in results)
    M_5_50 = np.sum(np.vstack([result.M_5_50 for result in results]), axis=0) if all_have_m_5_50 else None
    return ld_df, weight_df, M, M_5_50


def emit_outputs(results: Sequence[ChromComputationResult], args: argparse.Namespace) -> None:
    """Write either per-chromosome outputs or one aggregated LDSC-compatible result set."""
    if not results:
        raise ValueError("No chromosome results were produced.")

    write_annotation_groups(
        args.out + ".annotation_groups.tsv",
        results[0].baseline_columns,
        results[0].query_columns,
    )

    if args.per_chr_output:
        for result in results:
            prefix = f"{args.out}.{result.chrom}"
            write_ldscore_file(result_to_dataframe(result), prefix + ".l2.ldscore.gz")
            write_ldscore_file(weight_result_to_dataframe(result), prefix + ".w.l2.ldscore.gz")
            write_counts(prefix + ".l2.M", result.M)
            if result.M_5_50 is not None:
                write_counts(prefix + ".l2.M_5_50", result.M_5_50)
            else:
                LOGGER.warning("Skipping %s.l2.M_5_50 because MAF is unavailable.", prefix)
        return

    ld_df, weight_df, M, M_5_50 = aggregate_results(results)
    write_ldscore_file(ld_df, args.out + ".l2.ldscore.gz")
    write_ldscore_file(weight_df, args.out + ".w.l2.ldscore.gz")
    write_counts(args.out + ".l2.M", M)
    if M_5_50 is not None:
        write_counts(args.out + ".l2.M_5_50", M_5_50)
    else:
        LOGGER.warning("Skipping %s.l2.M_5_50 because MAF is unavailable.", args.out)


def validate_args(args: argparse.Namespace) -> None:
    if not (args.query_annot or args.query_annot_chr or args.baseline_annot or args.baseline_annot_chr):
        raise ValueError("At least one of query or baseline annotation inputs must be supplied.")
    if bool(args.r2_table or args.r2_table_chr) == bool(args.bfile or args.bfile_chr):
        raise ValueError("Specify exactly one reference-panel mode: parquet or PLINK.")
    if args.r2_table_chr and args.r2_table:
        raise ValueError("Cannot set both --r2-table and --r2-table-chr.")
    if args.bfile_chr and args.bfile:
        raise ValueError("Cannot set both --bfile and --bfile-chr.")
    if args.snp_identifier not in {"chr_pos", "rsID"}:
        raise ValueError("--snp-identifier must be one of: chr_pos, rsID.")
    if args.r2_table or args.r2_table_chr:
        if args.r2_bias_mode is None:
            raise ValueError("--r2-bias-mode is required in parquet mode.")
        if args.r2_bias_mode == "raw" and args.r2_sample_size is None:
            raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
        if args.snp_identifier == "chr_pos" and args.genome_build is None:
            raise ValueError("--genome-build is required in parquet mode when --snp-identifier chr_pos.")
    if args.ld_wind_cm is not None and args.ld_wind_cm <= 0:
        raise ValueError("--ld-wind-cm must be positive.")
    if args.ld_wind_kb is not None and args.ld_wind_kb <= 0:
        raise ValueError("--ld-wind-kb must be positive.")
    if args.ld_wind_snps is not None and args.ld_wind_snps <= 0:
        raise ValueError("--ld-wind-snps must be positive.")
    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be positive.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or parquet R2 input.",
    )
    parser.add_argument("--out", required=True, help="Output prefix.")
    parser.add_argument("--query-annot", default=None, help="Comma-separated SNP-level query annotation files or prefixes.")
    parser.add_argument("--query-annot-chr", default=None, help="Comma-separated chromosome-pattern prefixes for query annotation files.")
    parser.add_argument("--baseline-annot", default=None, help="Comma-separated SNP-level baseline annotation files or prefixes.")
    parser.add_argument("--baseline-annot-chr", default=None, help="Comma-separated chromosome-pattern prefixes for baseline annotation files.")
    parser.add_argument("--bfile", default=None, help="PLINK prefix for the reference panel.")
    parser.add_argument("--bfile-chr", default=None, help="Chromosome-pattern PLINK prefix for the reference panel.")
    parser.add_argument("--r2-table", default=None, help="Comma-separated parquet R2 files or prefixes.")
    parser.add_argument("--r2-table-chr", default=None, help="Comma-separated chromosome-pattern prefixes for parquet R2 files.")
    parser.add_argument("--snp-identifier", choices=("chr_pos", "rsID"), default="chr_pos", help="Identifier mode used to match annotations to the reference panel.")
    parser.add_argument("--genome-build", choices=("hg19", "hg38"), default=None, help="Genome build used for chr_pos parquet matching.")
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default=None, help="Whether parquet R2 values are raw sample r^2 or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
    parser.add_argument("--regression-snps", default=None, help="Optional SNP list defining the regression SNP set for w_ld.")
    parser.add_argument("--frqfile", default=None, help="Optional frequency/metadata file for MAF and CM.")
    parser.add_argument("--frqfile-chr", default=None, help="Chromosome-pattern frequency/metadata file prefix.")
    parser.add_argument("--ld-wind-snps", default=None, type=int, help="LD window size in SNPs.")
    parser.add_argument("--ld-wind-kb", default=None, type=float, help="LD window size in kilobases.")
    parser.add_argument("--ld-wind-cm", default=None, type=float, help="LD window size in centiMorgans.")
    parser.add_argument("--maf", default=None, type=float, help="Optional MAF filter for retained SNPs when MAF is available.")
    parser.add_argument("--chunk-size", default=50, type=int, help="Chunk size for legacy PLINK block computations.")
    parser.add_argument("--per-chr-output", default=False, action="store_true", help="Emit per-chromosome outputs instead of an aggregated output.")
    parser.add_argument("--yes-really", default=False, action="store_true", help="Allow whole-chromosome LD windows.")
    parser.add_argument("--log-level", default="INFO", choices=("DEBUG", "INFO", "WARNING", "ERROR"), help="Logging verbosity.")
    return parser


# Top-level workflow.
def run_ldscore(args: argparse.Namespace) -> list[ChromComputationResult]:
    """
    Execute the end-to-end LD-score workflow.

    Main steps:
    1. Validate CLI arguments and discover the chromosomes to process.
    2. Load and align baseline/query annotations for each chromosome.
    3. Dispatch to either the parquet or PLINK chromosome backend.
    4. Aggregate and write LDSC-compatible outputs.
    """
    validate_args(args)
    chromosomes = chromosome_set_from_annotation_inputs(args)
    regression_keys = load_regression_keys(args)
    results: list[ChromComputationResult] = []

    for chrom in chromosomes:
        # Resolve all annotation files contributing columns for this chromosome.
        baseline_files = resolve_optional_chr_files(args.baseline_annot_chr, chrom, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
        baseline_files += resolve_annotation_files(args.baseline_annot)
        query_files = resolve_optional_chr_files(args.query_annot_chr, chrom, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
        query_files += resolve_annotation_files(args.query_annot)

        bundle = combine_annotation_groups(
            baseline_files=baseline_files,
            query_files=query_files,
            chrom=chrom,
            identifier_mode=args.snp_identifier,
        )
        if bundle is None:
            continue

        LOGGER.info(
            "Processing chromosome %s with %d baseline and %d query annotation columns.",
            chrom,
            len(bundle.baseline_columns),
            len(bundle.query_columns),
        )

        # Compute LD scores from the selected reference-panel backend.
        if args.r2_table or args.r2_table_chr:
            result = compute_chrom_from_parquet(chrom, bundle, args, regression_keys)
        else:
            result = compute_chrom_from_plink(chrom, bundle, args, regression_keys)
        results.append(result)

    emit_outputs(results, args)
    return results


def main(argv: Sequence[str] | None = None) -> list[ChromComputationResult]:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return run_ldscore(args)


if __name__ == "__main__":
    main()
