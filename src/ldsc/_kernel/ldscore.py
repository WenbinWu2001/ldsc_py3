"""
Description
-----------
New LD-score estimation core for SNP-level annotation files with support for
either PLINK genotypes or sorted per-chromosome parquet `R2` tables as the
reference panel.

This module consumes SNP-level `.annot` / `.annot.gz` inputs only. It does not
perform BED-to-SNP projection; that remains the responsibility of
`utils/BED2annot.py`.

Supported Inputs
----------------
- Query and baseline annotation inputs are provided separately. Each file must
  contain SNP-level rows with required metadata columns `CHR POS SNP CM`, plus
  one or more annotation columns.
- Reference-panel input may be either:
  - PLINK genotype files (`.bed/.bim/.fam`), or
  - one per-chromosome sorted parquet `R2` file.
- Optional frequency metadata may be provided to populate `MAF` in outputs and
  to enable `.M_5_50` generation when parquet input is used.
- Optional regression SNP lists may be provided to define the SNP set used for
  `w_ld`. If omitted, the retained reference SNP set is used.

Sorted Parquet `R2` Format
--------------------------
- Runtime parquet input may be either:
  - normalized sorted parquet with helper columns `chr`, `pos_1`, `pos_2`, or
  - raw pairwise parquet carrying the legacy source columns listed below.
- The source table is expected to contain the pairwise columns:
  - `chr`
  - `rsID_1`, `rsID_2`
  - `hg38_pos_1`, `hg38_pos_2`
  - `hg19_pos_1`, `hg19_pos_2`
  - `hg38_Uniq_ID_1`, `hg38_Uniq_ID_2`
  - `hg19_Uniq_ID_1`, `hg19_Uniq_ID_2`
  - `R2`, `Dprime`, `+/-corr`
- The normalized sorted parquet created by this module preserves those source
  columns and adds helper columns:
  - `pos_1`
  - `pos_2`
- Pair orientation is canonicalized by the selected build so `pos_1 <= pos_2`.
- Rows are sorted by `(pos_1, pos_2)`.
- Exact duplicate normalized pairs with the same `R2` are dropped. Duplicate
  normalized pairs with different `R2` raise an error.
- In `chr_pos` mode, retained reference SNP rows must have unique chromosome
  positions. If two retained SNPs share the same `CHR` and `POS`, matching to
  the `R2` table is ambiguous and the code fails fast.

Identifier and Genome-Build Rules
---------------------------------
- Supported SNP identifier modes in v1 are `chr_pos` and `rsid`.
- No allele matching is attempted in v1.
- In `chr_pos` mode, retained reference SNP rows are matched by chromosome and
  position and must be unique within each chromosome.
- The sorted parquet helper is build-specific at normalization time. When the
  sorted parquet is used at runtime, the selected `genome_build` is assumed to
  match the build intended for the LD-score run.

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
    from ldsc_py3_Jerry.ldscore.ldscore_new import run_ldscore

    run_ldscore(
        out="results/example",
        baseline_annot_chr="data/baseline.@",
        query_annot="data/query.annot.gz",
        r2_table_chr="data/r2/chr@_sorted",
        snp_identifier="rsid",
        r2_bias_mode="unbiased",
        ld_wind_cm=1,
    )

CLI
    python -m ldsc_py3_Jerry.ldsc_new \
        --out results/example \
        --baseline-annot-chr data/baseline.@ \
        --query-annot data/query.annot.gz \
        --r2-table-chr data/r2/chr@_sorted \
        --snp-identifier rsid \
        --r2-bias-mode unbiased \
        --ld-wind-cm 1

Preprocessing helper
--------------------
- Use `convert_r2_table_to_sorted_parquet(source_path, genome_build, output_path)`
  to convert a common tabular `R2` file into the normalized sorted parquet
  format required by `run_ldscore(...)`.

Computation Overview
--------------------
1. Load query and baseline SNP-level annotation files chromosome by chromosome.
2. Normalize SNP identity using the selected identifier mode and require all
   annotation files for a chromosome to have identical retained SNP rows.
3. Build the canonical reference SNP table from annotation inputs.
4. Intersect annotation SNPs with SNPs actually present in the reference panel.
5. Compute LD scores per chromosome:
   - PLINK mode reuses the original genotype-block implementation.
   - parquet mode follows the old LDSC sliding-block accumulation workflow and
     fills block-local dense `R2` matrices from the sorted parquet file.
6. Compute all partitioned LD-score columns and `w_ld` from the same block
   traversal.
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
- Each unordered SNP pair contributes to two symmetric matrix entries in the
  block-local dense matrices.
- The diagonal is always added internally as `r_jj^2 = 1`.

Window Rules
------------
- `--ld-wind-snps`, `--ld-wind-kb`, and `--ld-wind-cm` follow the same semantics
  as the original codebase.
- In parquet mode, `block_left` and the old sliding-block iteration determine
  the query bounds and accumulation order.
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
- v1 supports only `chr_pos` and `rsid` matching.
- v1 does not attempt allele-aware reconciliation between annotations and
  parquet pairs.
- Runtime parquet input may be either normalized sorted parquet or raw parquet
  with the legacy pairwise schema.
- parquet input requires `pyarrow` at runtime.
- Per-query partitioned LDSC wrappers are intentionally out of scope for this
  module; this module computes LD scores only.

Dependencies
------------
- `numpy`
- `pandas`
- `pyarrow` for sorted parquet conversion and runtime queries
"""

from __future__ import annotations

import argparse
import glob
import gzip
import logging
import os
import sys
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from ..column_inference import (
    ANNOTATION_METADATA_SPEC_MAP,
    CHR_COLUMN_ALIASES,
    CHR_COLUMN_SPEC,
    CM_COLUMN_ALIASES,
    CM_COLUMN_SPEC,
    MAF_COLUMN_ALIASES,
    MAF_COLUMN_SPEC,
    POS_COLUMN_ALIASES,
    POS_COLUMN_SPEC,
    REFERENCE_METADATA_SPEC_MAP,
    R2_HELPER_COLUMN_SPECS,
    RESTRICTION_CHRPOS_SPEC_MAP,
    RESTRICTION_RSID_SPEC_MAP,
    R2_SOURCE_COLUMN_SPECS,
    SNP_COLUMN_ALIASES,
    SNP_COLUMN_SPEC,
    normalize_genome_build,
    normalize_snp_identifier_mode,
    resolve_optional_column,
    resolve_required_column,
    resolve_required_columns,
)
from ..chromosome_inference import STANDARD_CHROMOSOMES, chrom_sort_key, normalize_chromosome
from ..genome_build_inference import (
    load_packaged_reference_table,
    resolve_chr_pos_table,
    validate_auto_genome_build_mode,
)
from . import formats as legacy_parse
from .identifiers import build_snp_id_series, read_global_snp_restriction

try:  # pragma: no cover - optional dependency
    import bitarray as ba
except ImportError:  # pragma: no cover - optional dependency
    ba = None


LOGGER = logging.getLogger("LDSC.new")
REQUIRED_ANNOT_COLUMNS = ("CHR", "POS", "SNP", "CM")
ANNOT_META_COLUMNS = ("CHR", "SNP", "POS", "CM", "MAF")
CHROM_ALIASES = CHR_COLUMN_ALIASES
POS_ALIASES = POS_COLUMN_ALIASES
SNP_ALIASES = SNP_COLUMN_ALIASES
CM_ALIASES = CM_COLUMN_ALIASES
MAF_ALIASES = MAF_COLUMN_ALIASES
R2_CANONICAL_SOURCE_COLUMNS = tuple(spec.canonical for spec in R2_SOURCE_COLUMN_SPECS)
R2_HELPER_COLUMNS = tuple(spec.canonical for spec in R2_HELPER_COLUMN_SPECS)
SORTED_R2_REQUIRED_COLUMNS = R2_CANONICAL_SOURCE_COLUMNS + R2_HELPER_COLUMNS


@dataclass
class AnnotationBundle:
    """Per-chromosome annotation payload consumed by the LD-score backends."""
    metadata: pd.DataFrame
    annotations: pd.DataFrame
    baseline_columns: list[str]
    query_columns: list[str]


@dataclass
class ChromComputationResult:
    """Chromosome-level LD-score outputs before cross-chromosome aggregation."""
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
    """Configure the module logger for CLI-style execution."""
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(levelname)s: %(message)s",
    )


def get_legacy_ld_module():
    """Return the bitarray-backed LD kernel or raise a dependency error."""
    if ba is None:
        raise ImportError(
            "PLINK reference-panel mode requires the internal bitarray-backed LD kernel. "
            "Install bitarray and retry."
        )
    return sys.modules[__name__]


def split_arg_list(value: str | None) -> list[str]:
    """Split a comma-delimited CLI argument into normalized path tokens."""
    if not value:
        return []
    out = []
    for token in value.split(","):
        token = os.path.expanduser(os.path.expandvars(token.strip()))
        if token:
            out.append(token)
    return out


def find_column(columns: Iterable[str], aliases: Sequence[str]) -> str | None:
    """Return the first column whose normalized name matches one of ``aliases``."""
    spec = {
        tuple(CHROM_ALIASES): CHR_COLUMN_SPEC,
        tuple(POS_ALIASES): POS_COLUMN_SPEC,
        tuple(SNP_ALIASES): SNP_COLUMN_SPEC,
        tuple(CM_ALIASES): CM_COLUMN_SPEC,
        tuple(MAF_ALIASES): MAF_COLUMN_SPEC,
    }.get(tuple(aliases))
    if spec is not None:
        return resolve_optional_column(columns, spec)
    return None


def get_block_lefts(coords: np.ndarray, max_dist: float) -> np.ndarray:
    """
    Compute the standard LDSC ``block_left`` window start array.

    ``block_left[i]`` stores the first index retained in the LD window centered
    on SNP ``i`` for the chosen distance metric.
    """
    M = len(coords)
    j = 0
    block_left = np.zeros(M)
    for i in range(M):
        while j < M and abs(coords[j] - coords[i]) > max_dist:
            j += 1
        block_left[i] = j
    return block_left


def getBlockLefts(coords, max_dist):
    """Backward-compatible alias for :func:`get_block_lefts`."""
    return get_block_lefts(coords, max_dist)


def block_left_to_right(block_left):
    """Convert block-left coordinates to block-right coordinates."""
    M = len(block_left)
    j = 0
    block_right = np.zeros(M)
    for i in range(M):
        while j < M and block_left[j] <= i:
            j += 1
        block_right[i] = j
    return block_right


class __GenotypeArrayInMemory__(object):
    """Parent class for in-memory genotype matrices."""

    def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
        """Load, filter, and normalize one in-memory genotype matrix."""
        self.m = len(snp_list.IDList)
        self.n = n
        self.keep_snps = keep_snps
        self.keep_indivs = keep_indivs
        self.df = np.array(snp_list.df[["CHR", "SNP", "BP", "CM"]])
        self.colnames = ["CHR", "SNP", "POS", "CM"]
        self.mafMin = mafMin if mafMin is not None else 0
        self._currentSNP = 0
        (self.nru, self.geno) = self.__read__(fname, self.m, n)
        if keep_indivs is not None:
            keep_indivs = np.array(keep_indivs, dtype="int")
            if np.any(keep_indivs > self.n):
                raise ValueError("keep_indivs indices out of bounds")
            (self.geno, self.m, self.n) = self.__filter_indivs__(self.geno, keep_indivs, self.m, self.n)
            if self.n <= 0:
                raise ValueError("After filtering, no individuals remain")
        if keep_snps is not None:
            keep_snps = np.array(keep_snps, dtype="int")
            if np.any(keep_snps > self.m):
                raise ValueError("keep_snps indices out of bounds")
        (self.geno, self.m, self.n, self.kept_snps, self.freq) = self.__filter_snps_maf__(
            self.geno, self.m, self.n, self.mafMin, keep_snps
        )
        if self.m <= 0:
            raise ValueError("After filtering, no SNPs remain")
        self.df = self.df[self.kept_snps, :]
        self.maf = np.minimum(self.freq, np.ones(self.m) - self.freq)
        self.sqrtpq = np.sqrt(self.freq * (np.ones(self.m) - self.freq))
        self.df = np.c_[self.df, self.maf]
        self.colnames.append("MAF")

    def __read__(self, fname, m, n):
        """Read the backend-specific genotype representation into memory."""
        raise NotImplementedError

    def __filter_indivs__(geno, keep_indivs, m, n):
        """Apply backend-specific sample filtering to the genotype matrix."""
        raise NotImplementedError

    def __filter_maf_(geno, m, n, maf):
        """Apply backend-specific SNP and MAF filtering to the genotype matrix."""
        raise NotImplementedError

    def ldScoreVarBlocks(self, block_left, c, annot=None):
        """Compute LD-score block sums using the unbiased :math:`r^2` transform."""
        func = lambda x: self.__l2_unbiased__(x, self.n)
        snp_getter = self.nextSNPs
        return self.__corSumVarBlocks__(block_left, c, func, snp_getter, annot)

    def ldScoreBlockJackknife(self, block_left, c, annot=None, jN=10):
        """Compute block-jackknife LD-score summaries using squared correlations."""
        func = lambda x: np.square(x)
        snp_getter = self.nextSNPs
        return self.__corSumBlockJackknife__(block_left, c, func, snp_getter, annot, jN)

    def __l2_unbiased__(self, x, n):
        """Convert correlation values to the unbiased LD-score contribution."""
        denom = n - 2 if n > 2 else n
        sq = np.square(x)
        return sq - (1 - sq) / denom

    def __corSumVarBlocks__(self, block_left, c, func, snp_getter, annot=None):
        """Accumulate transformed correlation sums over LDSC-style LD blocks."""
        m, n = self.m, self.n
        block_sizes = np.array(np.arange(m) - block_left)
        block_sizes = np.ceil(block_sizes / c) * c
        if annot is None:
            annot = np.ones((m, 1))
        else:
            annot_m = annot.shape[0]
            if annot_m != self.m:
                raise ValueError("Incorrect number of SNPs in annot")

        n_a = annot.shape[1]
        cor_sum = np.zeros((m, n_a))
        b = np.nonzero(block_left > 0)
        if np.any(b):
            b = b[0][0]
        else:
            b = m
        b = int(np.ceil(b / c) * c)
        if b > m:
            c = 1
            b = m
        l_A = 0
        A = snp_getter(b)
        rfuncAB = np.zeros((b, c))
        rfuncBB = np.zeros((c, c))
        for l_B in range(0, b, c):
            B = A[:, l_B:l_B + c]
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + c, :])
        b0 = b
        md = int(c * np.floor(m / c))
        end = md + 1 if md != m else md
        for l_B in range(b0, end, c):
            old_b = b
            b = int(block_sizes[l_B])
            if l_B > b0 and b > 0:
                A = np.hstack((A[:, old_b - b + c:old_b], B))
                l_A += old_b - b + c
            elif l_B == b0 and b > 0:
                A = A[:, b0 - b:b0]
                l_A = b0 - b
            elif b == 0:
                A = np.array(()).reshape((n, 0))
                l_A = l_B
            if l_B == md:
                c = m - md
                rfuncAB = np.zeros((b, c))
                rfuncBB = np.zeros((c, c))
            if b != old_b:
                rfuncAB = np.zeros((b, c))
            B = snp_getter(c)
            p1 = np.all(annot[l_A:l_A + b, :] == 0)
            p2 = np.all(annot[l_B:l_B + c, :] == 0)
            if p1 and p2:
                continue
            np.dot(A.T, B / n, out=rfuncAB)
            rfuncAB = func(rfuncAB)
            cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + c, :])
            cor_sum[l_B:l_B + c, :] += np.dot(annot[l_A:l_A + b, :].T, rfuncAB).T
            np.dot(B.T, B / n, out=rfuncBB)
            rfuncBB = func(rfuncBB)
            cor_sum[l_B:l_B + c, :] += np.dot(rfuncBB, annot[l_B:l_B + c, :])
        return cor_sum


if ba is not None:
    class PlinkBEDFile(__GenotypeArrayInMemory__):
        """Interface for PLINK .bed format."""

        def __init__(self, fname, n, snp_list, keep_snps=None, keep_indivs=None, mafMin=None):
            """Initialize the PLINK reader and configure BED bit-pattern decoding."""
            self._bedcode = {
                2: ba.bitarray("11"),
                9: ba.bitarray("10"),
                1: ba.bitarray("01"),
                0: ba.bitarray("00"),
            }
            __GenotypeArrayInMemory__.__init__(
                self,
                fname,
                n,
                snp_list,
                keep_snps=keep_snps,
                keep_indivs=keep_indivs,
                mafMin=mafMin,
            )

        def __read__(self, fname, m, n):
            """Read the raw PLINK BED payload and validate the file header."""
            if not fname.endswith(".bed"):
                raise ValueError(".bed filename must end in .bed")
            fh = open(fname, "rb")
            magicNumber = ba.bitarray(endian="little")
            magicNumber.fromfile(fh, 2)
            bedMode = ba.bitarray(endian="little")
            bedMode.fromfile(fh, 1)
            e = (4 - n % 4) if n % 4 != 0 else 0
            nru = n + e
            self.nru = nru
            if magicNumber != ba.bitarray("0011011011011000"):
                raise IOError("Magic number from Plink .bed file not recognized")
            if bedMode != ba.bitarray("10000000"):
                raise IOError("Plink .bed file must be in default SNP-major mode")
            self.geno = ba.bitarray(endian="little")
            self.geno.fromfile(fh)
            self.__test_length__(self.geno, self.m, self.nru)
            return (self.nru, self.geno)

        def __test_length__(self, geno, m, nru):
            """Validate that the BED payload length matches the expected shape."""
            exp_len = 2 * m * nru
            real_len = len(geno)
            if real_len != exp_len:
                raise IOError("Plink .bed file has {n1} bits, expected {n2}".format(n1=real_len, n2=exp_len))

        def __filter_indivs__(self, geno, keep_indivs, m, n):
            """Subset the BED bitarray to the requested individuals."""
            n_new = len(keep_indivs)
            e = (4 - n_new % 4) if n_new % 4 != 0 else 0
            nru_new = n_new + e
            nru = self.nru
            z = ba.bitarray(m * 2 * nru_new, endian="little")
            z.setall(0)
            for e, i in enumerate(keep_indivs):
                z[2 * e::2 * nru_new] = geno[2 * i::2 * nru]
                z[2 * e + 1::2 * nru_new] = geno[2 * i + 1::2 * nru]
            self.nru = nru_new
            return (z, m, n_new)

        def __filter_snps_maf__(self, geno, m, n, mafMin, keep_snps):
            """Filter SNPs by explicit keep list and minor-allele frequency."""
            nru = self.nru
            m_poly = 0
            y = ba.bitarray()
            if keep_snps is None:
                keep_snps = range(m)
            kept_snps = []
            freq = []
            for j in keep_snps:
                z = geno[2 * nru * j:2 * nru * (j + 1)]
                A = z[0::2]
                a = A.count()
                B = z[1::2]
                b = B.count()
                c = (A & B).count()
                major_ct = b + c
                n_nomiss = n - a + c
                f = major_ct / (2 * n_nomiss) if n_nomiss > 0 else 0
                het_miss_ct = a + b - 2 * c
                if np.minimum(f, 1 - f) > mafMin and het_miss_ct < n:
                    freq.append(f)
                    y += z
                    m_poly += 1
                    kept_snps.append(j)
            return (y, m_poly, n, kept_snps, freq)

        def nextSNPs(self, b, minorRef=None):
            """Return the next ``b`` standardized SNP columns from the BED stream."""
            try:
                b = int(b)
                if b <= 0:
                    raise ValueError("b must be > 0")
            except TypeError:
                raise TypeError("b must be an integer")
            if self._currentSNP + b > self.m:
                raise ValueError("{b} SNPs requested, {k} SNPs remain".format(b=b, k=(self.m - self._currentSNP)))
            c = self._currentSNP
            n = self.n
            nru = self.nru
            slice = self.geno[2 * c * nru:2 * (c + b) * nru]
            X = np.array(slice.decode(self._bedcode), dtype="float64").reshape((b, nru)).T
            X = X[0:n, :]
            Y = np.zeros(X.shape)
            for j in range(0, b):
                newsnp = X[:, j]
                ii = newsnp != 9
                avg = np.mean(newsnp[ii])
                newsnp[np.logical_not(ii)] = avg
                denom = np.std(newsnp)
                if denom == 0:
                    denom = 1
                if minorRef is not None and self.freq[self._currentSNP + j] > 0.5:
                    denom = denom * -1
                Y[:, j] = (newsnp - avg) / denom
            self._currentSNP += b
            return Y
else:
    class PlinkBEDFile:  # pragma: no cover - dependency-gated fallback
        """Fallback PLINK reader that raises when `bitarray` is unavailable."""
        def __init__(self, *args, **kwargs):
            """Raise an informative import error for dependency-gated PLINK support."""
            raise ImportError("PLINK LD-score support requires the optional dependency 'bitarray'.")


def identifier_keys(df: pd.DataFrame, mode: str) -> pd.Series:
    """Build the canonical SNP identifier series used for matching within the kernel."""
    return build_snp_id_series(df, normalize_snp_identifier_mode(mode))


def sort_frame_by_genomic_position(df: pd.DataFrame) -> pd.DataFrame:
    """Sort a metadata-like frame by chromosome, position, and SNP name."""
    pos_col = find_column(df.columns, POS_ALIASES)
    if pos_col is None:
        raise KeyError("No POS-like column available for genomic sorting.")
    sort_df = df.copy()
    sort_df["_chrom_key"] = sort_df["CHR"].map(chrom_sort_key)
    sort_df = sort_df.sort_values(by=["_chrom_key", pos_col, "SNP"], kind="mergesort")
    return sort_df.drop(columns="_chrom_key").reset_index(drop=True)


def resolve_prefixed_file(token: str, suffixes: Sequence[str]) -> list[str]:
    """Resolve one explicit path, glob, or filename prefix into concrete files."""
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
    """Resolve comma-delimited annotation inputs into concrete file paths."""
    files: list[str] = []
    for token in split_arg_list(spec):
        files.extend(resolve_prefixed_file(token, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz")))
    return files


def resolve_chr_files(spec: str | None, chrom: str, suffixes: Sequence[str]) -> list[str]:
    """Resolve chromosome-pattern file prefixes for one chromosome."""
    files: list[str] = []
    for token in split_arg_list(spec):
        chrom_token = legacy_parse.sub_chr(token, chrom)
        files.extend(resolve_prefixed_file(chrom_token, suffixes))
    return files


def resolve_optional_chr_files(spec: str | None, chrom: str, suffixes: Sequence[str]) -> list[str]:
    """Resolve chromosome-pattern file prefixes while ignoring missing chromosomes."""
    files: list[str] = []
    for token in split_arg_list(spec):
        chrom_token = legacy_parse.sub_chr(token, chrom)
        try:
            files.extend(resolve_prefixed_file(chrom_token, suffixes))
        except FileNotFoundError:
            continue
    return files


def resolve_parquet_files(args: argparse.Namespace, chrom: str | None = None) -> list[str]:
    """Resolve the sorted parquet R2 files participating in one LD-score run."""
    if chrom is not None and args.r2_table_chr:
        return resolve_chr_files(args.r2_table_chr, chrom, ("", ".parquet"))
    if args.r2_table:
        files = []
        for token in split_arg_list(args.r2_table):
            files.extend(resolve_prefixed_file(token, ("", ".parquet")))
        return files
    return []


def resolve_bfile_prefix(args: argparse.Namespace, chrom: str | None = None) -> str | None:
    """Resolve the PLINK prefix for the requested chromosome, if any."""
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
    """Resolve optional frequency or metadata files for one chromosome."""
    files: list[str] = []
    if chrom is not None and args.frqfile_chr:
        return resolve_optional_chr_files(args.frqfile_chr, chrom, ("", ".frq", ".frq.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
    for token in split_arg_list(args.frqfile):
        files.extend(resolve_prefixed_file(token, ("", ".frq", ".frq.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz")))
    return files


def read_text_table(path: str) -> pd.DataFrame:
    """Read a whitespace-delimited kernel input table with optional gzip compression."""
    compression = "gzip" if path.endswith(".gz") else None
    return pd.read_csv(path, sep=r"\s+", compression=compression)


def _resolve_r2_source_columns(schema_names: Iterable[str], context: str | None = None) -> dict[str, str]:
    """Resolve raw R2 source columns from an on-disk schema into canonical slots."""
    return resolve_required_columns(schema_names, R2_SOURCE_COLUMN_SPECS, context=context)


def _resolve_r2_source_subset(
    schema_names: Iterable[str],
    canonicals: Sequence[str],
    context: str | None = None,
) -> dict[str, str]:
    """Resolve a required subset of raw R2 source columns from ``schema_names``."""
    spec_map = {spec.canonical: spec for spec in R2_SOURCE_COLUMN_SPECS}
    return {
        canonical: resolve_required_column(schema_names, spec_map[canonical], context=context)
        for canonical in canonicals
    }


def normalize_r2_source_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename raw pairwise-R2 source columns to the canonical POS-based schema."""
    rename_map = {
        actual: canonical
        for canonical, actual in _resolve_r2_source_columns(df.columns).items()
        if actual != canonical
    }
    return df.rename(columns=rename_map)


def _resolve_helper_columns(schema_names: Iterable[str], context: str | None = None) -> tuple[str, str] | None:
    """Resolve normalized parquet helper position columns from schema aliases."""
    chr_col = resolve_optional_column(schema_names, R2_SOURCE_COLUMN_SPECS[0], context=context)
    pos_1 = resolve_optional_column(schema_names, R2_HELPER_COLUMN_SPECS[0], context=context)
    pos_2 = resolve_optional_column(schema_names, R2_HELPER_COLUMN_SPECS[1], context=context)
    if chr_col is None or pos_1 is None or pos_2 is None:
        return None
    return pos_1, pos_2


def get_r2_build_columns(genome_build: str, schema_names: Iterable[str] | None = None) -> tuple[str, str]:
    """Return the source R2-table coordinate columns for the selected build."""
    genome_build = _require_runtime_genome_build(genome_build)
    if genome_build == "hg19":
        canonical = ("hg19_pos_1", "hg19_pos_2")
    else:
        canonical = ("hg38_pos_1", "hg38_pos_2")
    if schema_names is None:
        return canonical
    mapping = _resolve_r2_source_subset(schema_names, canonical)
    return mapping[canonical[0]], mapping[canonical[1]]


def get_pyarrow_modules():
    """Import the pyarrow dataset module or raise a user-facing dependency error."""
    try:
        import pyarrow.dataset as ds
    except ImportError as exc:
        raise ImportError(
            "pyarrow is required for sorted parquet R2 input. Install pyarrow and retry."
        ) from exc
    return ds


def _parquet_schema_layout(schema_names: Sequence[str]) -> str:
    """Classify a runtime parquet schema as normalized, raw, or unsupported."""
    if _resolve_helper_columns(schema_names) is not None:
        return "normalized"
    try:
        _resolve_r2_source_columns(schema_names)
    except Exception:
        return "unsupported"
    return "raw"


def _require_runtime_genome_build(genome_build: str | None) -> str:
    """Require a build whenever raw parquet rows must be canonicalized at runtime."""
    normalized = normalize_genome_build(genome_build)
    if normalized in {"hg19", "hg38"}:
        return normalized
    raise ValueError(
        "Raw parquet R2 input requires genome_build='hg19'/'hg37'/'GRCh37' or 'hg38'/'GRCh38' so runtime normalization can choose the coordinate columns."
    )


def read_common_tabular_r2(path: str) -> pd.DataFrame:
    """Read one pairwise-R2 source file from parquet, CSV, TSV, or inferred text."""
    lower = path.lower()
    if lower.endswith(".parquet"):
        try:
            return pd.read_parquet(path)
        except ImportError as exc:
            raise ImportError(
                "Reading parquet R2 input requires pyarrow or fastparquet."
            ) from exc
    if lower.endswith(".csv") or lower.endswith(".csv.gz"):
        return pd.read_csv(path)
    if lower.endswith(".tsv") or lower.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep=None, engine="python")


def validate_r2_source_columns(df: pd.DataFrame, path: str) -> None:
    """Validate that a source R2 table contains the required legacy columns."""
    try:
        _resolve_r2_source_columns(df.columns, context=path)
    except ValueError as exc:
        raise ValueError(f"{path} is missing required R2 columns.") from exc


def canonicalize_r2_pairs(df: pd.DataFrame, genome_build: str) -> pd.DataFrame:
    """Orient pairwise-R2 rows so the selected build always satisfies ``pos_1 <= pos_2``."""
    df = normalize_r2_source_columns(df.copy())
    left_pos_col, right_pos_col = get_r2_build_columns(genome_build, df.columns)
    df["chr"] = df["chr"].map(lambda value: normalize_chromosome(value, context="R2 chromosome column"))

    paired_columns = (
        ("rsID_1", "rsID_2"),
        ("hg38_pos_1", "hg38_pos_2"),
        ("hg19_pos_1", "hg19_pos_2"),
        ("hg38_Uniq_ID_1", "hg38_Uniq_ID_2"),
        ("hg19_Uniq_ID_1", "hg19_Uniq_ID_2"),
    )
    swap = pd.to_numeric(df[left_pos_col], errors="raise") > pd.to_numeric(df[right_pos_col], errors="raise")
    for left_col, right_col in paired_columns:
        left_values = df.loc[swap, left_col].copy()
        df.loc[swap, left_col] = df.loc[swap, right_col].to_numpy()
        df.loc[swap, right_col] = left_values.to_numpy()

    df["pos_1"] = pd.to_numeric(df[left_pos_col], errors="raise").astype(np.int64)
    df["pos_2"] = pd.to_numeric(df[right_pos_col], errors="raise").astype(np.int64)
    return df


def deduplicate_normalized_r2_pairs(df: pd.DataFrame) -> pd.DataFrame:
    """Drop exact duplicate normalized pairs and reject conflicting duplicate R2 rows."""
    pair_key_columns = [
        "chr",
        "rsID_1",
        "rsID_2",
        "hg38_pos_1",
        "hg38_pos_2",
        "hg19_pos_1",
        "hg19_pos_2",
        "hg38_Uniq_ID_1",
        "hg38_Uniq_ID_2",
        "hg19_Uniq_ID_1",
        "hg19_Uniq_ID_2",
    ]
    duplicated = df.duplicated(subset=pair_key_columns, keep=False)
    if not duplicated.any():
        return df

    dup_df = df.loc[duplicated, pair_key_columns + ["R2"]].copy()
    nunique = dup_df.groupby(pair_key_columns, dropna=False)["R2"].nunique(dropna=False)
    inconsistent = nunique[nunique > 1]
    if len(inconsistent) > 0:
        raise ValueError("Duplicate normalized R2 pairs detected with conflicting R2 values.")

    return df.drop_duplicates(subset=pair_key_columns, keep="first").reset_index(drop=True)


def convert_r2_table_to_sorted_parquet(source_path: str, genome_build: str, output_path: str) -> str:
    """
    Convert a common tabular pairwise-R2 file into the normalized sorted parquet
    format required by the runtime parquet LD-score path.
    """
    source_path = os.path.expanduser(os.path.expandvars(source_path))
    output_path = os.path.expanduser(os.path.expandvars(output_path))
    df = read_common_tabular_r2(source_path)
    validate_r2_source_columns(df, source_path)
    df = canonicalize_r2_pairs(df, genome_build)
    df = deduplicate_normalized_r2_pairs(df)
    chroms = df["chr"].dropna().map(lambda value: normalize_chromosome(value, context=source_path)).unique().tolist()
    if len(chroms) != 1:
        raise ValueError("Each sorted parquet R2 file must contain exactly one chromosome.")
    df = df.sort_values(by=["pos_1", "pos_2"], kind="mergesort").reset_index(drop=True)
    parent = os.path.dirname(output_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    try:
        df.to_parquet(output_path, index=False)
    except ImportError as exc:
        raise ImportError(
            "Writing sorted parquet R2 files requires pyarrow or fastparquet."
        ) from exc
    return output_path


def validate_retained_identifier_uniqueness(metadata: pd.DataFrame, identifier_mode: str, chrom: str) -> None:
    """Reject ambiguous retained SNP identifiers before parquet matching begins."""
    if identifier_mode == "chr_pos":
        duplicated = metadata.duplicated(subset=["CHR", "POS"], keep=False)
        if duplicated.any():
            raise ValueError(
                f"Chromosome {chrom} has duplicate retained SNP positions. "
                "This is ambiguous in chr_pos mode."
            )
        return

    duplicated = metadata["SNP"].duplicated(keep=False)
    if duplicated.any():
        raise ValueError(
            f"Chromosome {chrom} has duplicate retained SNP IDs. "
            "This is ambiguous in rsid mode."
        )


def read_sorted_r2_presence(args: argparse.Namespace, chrom: str) -> set[object]:
    """Read the SNP identifiers or positions present in one parquet chromosome source."""
    files = resolve_parquet_files(args, chrom=chrom)
    if not files:
        raise FileNotFoundError(f"No sorted parquet R2 files resolved for chromosome {chrom}.")

    ds = get_pyarrow_modules()
    dataset = ds.dataset(files, format="parquet")
    layout = _parquet_schema_layout(dataset.schema.names)
    normalized_chrom = normalize_chromosome(chrom)
    identifier_mode = normalize_snp_identifier_mode(args.snp_identifier)
    if layout == "normalized":
        helper_pos1, helper_pos2 = _resolve_helper_columns(dataset.schema.names, context="normalized parquet R2 schema")
        columns = ["chr"]
        if identifier_mode == "rsid":
            columns.extend(["rsID_1", "rsID_2"])
        else:
            columns.extend([helper_pos1, helper_pos2])
        filter_expr = ds.field("chr") == normalized_chrom
        table = dataset.to_table(columns=columns, filter=filter_expr)
    elif layout == "raw":
        genome_build = _require_runtime_genome_build(getattr(args, "genome_build", None))
        left_pos_col, right_pos_col = get_r2_build_columns(genome_build, dataset.schema.names)
        columns = ["chr"]
        if identifier_mode == "rsid":
            columns.extend(["rsID_1", "rsID_2"])
        else:
            columns.extend([left_pos_col, right_pos_col])
        filter_expr = ds.field("chr") == normalized_chrom
        table = dataset.to_table(columns=columns, filter=filter_expr)
    else:
        raise ValueError(
            "Parquet R2 input must contain either normalized helper columns "
            "`chr`, `pos_1`, `pos_2` or the raw legacy pairwise columns."
        )
    pairs = table.to_pandas()
    if layout == "raw" and len(pairs) > 0:
        pairs = pairs.loc[
            pairs["chr"].map(lambda value: normalize_chromosome(value, context="raw parquet R2 schema")) == normalized_chrom
        ].reset_index(drop=True)
    if identifier_mode == "rsid":
        return set(pairs["rsID_1"].astype(str)).union(set(pairs["rsID_2"].astype(str)))
    if layout == "normalized":
        left_series = pd.to_numeric(pairs[helper_pos1], errors="raise").astype(np.int64)
        right_series = pd.to_numeric(pairs[helper_pos2], errors="raise").astype(np.int64)
    else:
        genome_build = _require_runtime_genome_build(getattr(args, "genome_build", None))
        left_pos_col, right_pos_col = get_r2_build_columns(genome_build, pairs.columns)
        left_series = pd.to_numeric(pairs[left_pos_col], errors="raise").astype(np.int64)
        right_series = pd.to_numeric(pairs[right_pos_col], errors="raise").astype(np.int64)
    return set(left_series).union(set(right_series))


def filter_reference_to_present_r2(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    regression_keys: set[str] | None,
    present_values: set[object],
    identifier_mode: str,
    chrom: str,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str] | None]:
    """Restrict annotation SNPs to those that are actually present in the R2 source."""
    metadata = metadata.copy()
    if normalize_snp_identifier_mode(identifier_mode) == "rsid":
        keep = metadata["SNP"].astype(str).isin(present_values)
    else:
        keep = metadata["POS"].astype(np.int64).isin(present_values)
    removed = int((~keep).sum())
    if removed:
        LOGGER.warning(
            "Dropping %d annotated SNPs on chromosome %s because they are absent from the R2 table.",
            removed,
            chrom,
        )
    metadata = metadata.loc[keep].reset_index(drop=True)
    annotations = annotations.loc[keep].reset_index(drop=True)
    if regression_keys is not None:
        regression_keys = regression_keys.intersection(set(identifier_keys(metadata, identifier_mode)))
    return metadata, annotations, regression_keys


# Annotation loading and normalization.
def parse_annotation_file(path: str, chrom: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse one SNP-level annotation table into normalized metadata and values."""
    df = read_text_table(path)
    context = path
    chr_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CHR"], context=context)
    pos_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["POS"], context=context)
    snp_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["SNP"], context=context)
    cm_col = resolve_required_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["CM"], context=context)
    maf_col = resolve_optional_column(df.columns, ANNOTATION_METADATA_SPEC_MAP["MAF"], context=context)

    meta = pd.DataFrame(
        {
            "CHR": df[chr_col],
            "POS": df[pos_col],
            "SNP": df[snp_col],
            "CM": df[cm_col],
        }
    )
    meta["CHR"] = meta["CHR"].map(lambda value: normalize_chromosome(value, context=path))
    meta["POS"] = pd.to_numeric(meta["POS"], errors="raise").astype(np.int64)
    meta["SNP"] = meta["SNP"].astype(str)
    meta["CM"] = pd.to_numeric(meta["CM"], errors="coerce")
    if maf_col is not None:
        meta["MAF"] = pd.to_numeric(df[maf_col], errors="coerce")

    if chrom is not None:
        keep = meta["CHR"] == normalize_chromosome(chrom, context=path)
        meta = meta.loc[keep].reset_index(drop=True)
        df = df.loc[keep].reset_index(drop=True)
    else:
        meta = meta.reset_index(drop=True)

    if len(meta) == 0:
        return meta, pd.DataFrame(index=meta.index)

    meta = sort_frame_by_genomic_position(meta)
    annotation_columns = [col for col in df.columns if col not in {chr_col, pos_col, snp_col, cm_col, maf_col}]
    if not annotation_columns:
        raise ValueError(f"{path} does not contain any annotation columns.")

    sorted_df = df.copy()
    sorted_df["_CHR"] = sorted_df[chr_col].map(lambda value: normalize_chromosome(value, context=path))
    sorted_df["_POS"] = pd.to_numeric(sorted_df[pos_col], errors="raise").astype(np.int64)
    sorted_df["_SNP"] = sorted_df[snp_col].astype(str)
    sorted_df["_chrom_key"] = sorted_df["_CHR"].map(chrom_sort_key)
    sorted_df = sorted_df.sort_values(by=["_chrom_key", "_POS", "_SNP"], kind="mergesort").drop(columns="_chrom_key")
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
                same_core = reference.loc[:, ["CHR", "POS", "SNP"]].equals(meta.loc[:, ["CHR", "POS", "SNP"]])
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
    """Read a SNP list file into the canonical identifier set for ``mode``."""
    return read_global_snp_restriction(path, normalize_snp_identifier_mode(mode))


def load_regression_keys(args: argparse.Namespace) -> set[str] | None:
    """Load the optional regression SNP universe from CLI arguments."""
    if not args.regression_snps:
        return None
    return read_identifier_list(args.regression_snps, args.snp_identifier)


def parse_frequency_metadata(path: str, chrom: str | None, identifier_mode: str) -> pd.DataFrame:
    """Parse an optional frequency metadata file into ``_key``/``CM``/``MAF`` columns."""
    df = read_text_table(path)
    identifier_mode = normalize_snp_identifier_mode(identifier_mode)
    context = path
    chr_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=context)
    pos_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=context)
    snp_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=context)
    cm_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=context)
    maf_col = resolve_optional_column(df.columns, REFERENCE_METADATA_SPEC_MAP["MAF"], context=context)

    if chrom is not None and chr_col is not None:
        keep = df[chr_col].map(lambda value: normalize_chromosome(value, context=path)) == normalize_chromosome(chrom, context=path)
        df = df.loc[keep].reset_index(drop=True)

    if len(df) == 0:
        return pd.DataFrame(columns=["_key", "CM", "MAF"])

    out = pd.DataFrame(index=df.index)
    if identifier_mode == "rsid":
        if snp_col is None:
            raise ValueError(f"{path} must contain a SNP column in rsid mode.")
        out["_key"] = df[snp_col].astype(str)
    else:
        if chr_col is None or pos_col is None:
            raise ValueError(f"{path} must contain CHR and POS columns in chr_pos mode.")
        chrom_norm = df[chr_col].map(lambda value: normalize_chromosome(value, context=path))
        pos = pd.to_numeric(df[pos_col], errors="raise").astype(np.int64).astype(str)
        out["_key"] = chrom_norm + ":" + pos

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
    """Merge optional CM and MAF metadata onto the retained reference SNP table."""
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
    """Apply the optional MAF threshold to metadata and annotation rows together."""
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
    """Discover the chromosome set implied by the supplied annotation inputs."""
    chromosomes: set[str] = set()
    all_files = (
        resolve_annotation_files(args.query_annot)
        + resolve_annotation_files(args.baseline_annot)
    )
    for path in all_files:
        df = read_text_table(path)
        if "CHR" not in df.columns:
            raise ValueError(f"{path} is missing CHR.")
        chromosomes.update(df["CHR"].map(lambda value: normalize_chromosome(value, context=path)).unique().tolist())

    chr_specs = split_arg_list(args.query_annot_chr) + split_arg_list(args.baseline_annot_chr)
    if chr_specs:
        for chrom in STANDARD_CHROMOSOMES:
            files = resolve_optional_chr_files(args.query_annot_chr, chrom, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
            files += resolve_optional_chr_files(args.baseline_annot_chr, chrom, ("", ".annot", ".annot.gz", ".txt", ".txt.gz", ".tsv", ".tsv.gz"))
            if files:
                chromosomes.add(normalize_chromosome(chrom))

    if not chromosomes:
        raise ValueError("No annotation chromosomes could be resolved from the supplied inputs.")
    return sorted(chromosomes, key=chrom_sort_key)


def build_window_coordinates(metadata: pd.DataFrame, args: argparse.Namespace) -> tuple[np.ndarray, float]:
    """Build the coordinate array and maximum distance for the active LD window mode."""
    selectors = np.array([args.ld_wind_snps is not None, args.ld_wind_kb is not None, args.ld_wind_cm is not None], dtype=bool)
    if selectors.sum() != 1:
        raise ValueError("Must specify exactly one of --ld-wind-snps, --ld-wind-kb, or --ld-wind-cm.")

    if args.ld_wind_snps is not None:
        return np.arange(len(metadata), dtype=float), float(args.ld_wind_snps)
    if args.ld_wind_kb is not None:
        return metadata["POS"].to_numpy(dtype=float), float(args.ld_wind_kb) * 1000.0
    if metadata["CM"].isna().any():
        raise ValueError("--ld-wind-cm requires non-missing CM values for all retained SNPs.")
    return metadata["CM"].to_numpy(dtype=float), float(args.ld_wind_cm)


def check_whole_chromosome_window(block_left: np.ndarray, args: argparse.Namespace, chrom: str) -> None:
    """Guard against accidental whole-chromosome windows unless explicitly allowed."""
    if len(block_left) == 0:
        return
    if block_left[-1] == 0 and not args.yes_really:
        raise ValueError(
            f"Chromosome {chrom} would use a whole-chromosome LD window. Re-run with --yes-really to allow this."
        )


# Parquet R2 adapter.
class SortedR2BlockReader:
    """
    Query block-local dense R2 matrices from a normalized sorted per-chromosome
    parquet table.

    The runtime parquet contract is path-based. This reader assumes the file
    already contains `chr`, `pos_1`, and `pos_2`, and that pair orientation has
    been canonicalized so `pos_1 <= pos_2`.
    """

    def __init__(
        self,
        paths: Sequence[str],
        chrom: str,
        metadata: pd.DataFrame,
        identifier_mode: str,
        r2_bias_mode: str,
        r2_sample_size: float | None,
        genome_build: str | None = None,
    ) -> None:
        """Open one chromosome's sorted parquet R2 tables and build index maps."""
        if not paths:
            raise FileNotFoundError(f"No sorted parquet R2 files resolved for chromosome {chrom}.")
        ds = get_pyarrow_modules()
        self.dataset = ds.dataset(list(paths), format="parquet")
        self.ds = ds
        self.chrom = normalize_chromosome(chrom)
        self.identifier_mode = identifier_mode
        self.r2_bias_mode = r2_bias_mode
        self.r2_sample_size = r2_sample_size
        self.genome_build = genome_build
        metadata = metadata.copy()
        metadata_context = f"SortedR2BlockReader[{self.chrom}] metadata"
        renamed = {
            resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["CHR"], context=metadata_context): "CHR",
            resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["POS"], context=metadata_context): "POS",
            resolve_required_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["SNP"], context=metadata_context): "SNP",
        }
        optional_cm = resolve_optional_column(metadata.columns, REFERENCE_METADATA_SPEC_MAP["CM"], context=metadata_context)
        if optional_cm is not None:
            renamed[optional_cm] = "CM"
        metadata = metadata.rename(columns=renamed)
        if self.genome_build == "auto":
            validate_auto_genome_build_mode(self.identifier_mode, self.genome_build)
            metadata, inference = resolve_chr_pos_table(
                metadata,
                context=f"SortedR2BlockReader[{self.chrom}]",
                reference_table=load_packaged_reference_table(),
                logger=LOGGER,
            )
            self.genome_build = inference.genome_build
        validate_retained_identifier_uniqueness(metadata, identifier_mode, chrom)
        self.pos = metadata["POS"].to_numpy(dtype=np.int64)
        self.m = len(metadata)
        self.query_columns = ["chr", "pos_1", "pos_2", "R2"]
        self._runtime_layout = _parquet_schema_layout(self.dataset.schema.names)
        self._raw_pos_columns: tuple[str, str] | None = None
        self._raw_query_columns: list[str] | None = None
        if self.identifier_mode == "rsid":
            self.query_columns.extend(["rsID_1", "rsID_2"])
            self.index_map = {str(snp): idx for idx, snp in enumerate(metadata["SNP"].astype(str))}
        else:
            self.index_map = {int(pos): idx for idx, pos in enumerate(metadata["POS"].astype(np.int64))}

        if self._runtime_layout == "normalized":
            missing_columns = sorted(set(self.query_columns) - set(self.dataset.schema.names))
            if missing_columns:
                raise ValueError(
                    "Sorted parquet R2 input is missing required normalized columns: "
                    + ", ".join(missing_columns)
                )
        elif self._runtime_layout == "raw":
            self.genome_build = _require_runtime_genome_build(self.genome_build)
            raw_mapping = _resolve_r2_source_columns(self.dataset.schema.names, context="raw parquet R2 schema")
            self._raw_pos_columns = get_r2_build_columns(self.genome_build, self.dataset.schema.names)
            self._raw_query_columns = [raw_mapping[canonical] for canonical in R2_CANONICAL_SOURCE_COLUMNS]
        else:
            raise ValueError(
                "Parquet R2 input must contain either normalized helper columns "
                "`chr`, `pos_1`, `pos_2` or the raw legacy pairwise columns."
            )

        self._last_query_key: tuple[int, int] | None = None
        self._last_query_rows: pd.DataFrame | None = None

    def _transform_r2(self, values: np.ndarray) -> np.ndarray:
        """Apply the configured raw-to-unbiased R2 correction when required."""
        values = values.astype(np.float32, copy=False)
        if self.r2_bias_mode == "raw":
            if self.r2_sample_size is None:
                raise ValueError("--r2-sample-size is required when --r2-bias-mode raw.")
            denom = self.r2_sample_size - 2
            if denom <= 0:
                raise ValueError("--r2-sample-size must be greater than 2 for raw R2 correction.")
            values = values - (1.0 - values) / denom
        return values

    def _query_union_rows(self, pos_min: int, pos_max: int) -> pd.DataFrame:
        """Query cached or on-disk pair rows spanning a union genomic window."""
        key = (int(pos_min), int(pos_max))
        if self._last_query_key == key and self._last_query_rows is not None:
            return self._last_query_rows.copy()

        if self._runtime_layout == "normalized":
            helper_pos1, helper_pos2 = _resolve_helper_columns(self.dataset.schema.names, context="normalized parquet R2 schema")
            filter_expr = (
                (self.ds.field("chr") == self.chrom)
                & (self.ds.field(helper_pos1) >= int(pos_min))
                & (self.ds.field(helper_pos2) <= int(pos_max))
            )
            table = self.dataset.to_table(columns=self.query_columns, filter=filter_expr)
        else:
            left_pos_col, right_pos_col = self._raw_pos_columns
            filter_expr = (
                (self.ds.field("chr") == self.chrom)
                & (self.ds.field(left_pos_col) >= int(pos_min))
                & (self.ds.field(left_pos_col) <= int(pos_max))
                & (self.ds.field(right_pos_col) >= int(pos_min))
                & (self.ds.field(right_pos_col) <= int(pos_max))
            )
            table = self.dataset.to_table(columns=self._raw_query_columns, filter=filter_expr)
        rows = table.to_pandas()
        if len(rows) == 0:
            base_columns = ["chr", "pos_1", "pos_2", "R2"]
            if self.identifier_mode == "rsid":
                base_columns.extend(["rsID_1", "rsID_2"])
            rows = pd.DataFrame(columns=base_columns + ["i", "j"])
        else:
            if self._runtime_layout == "raw":
                left_pos_col, right_pos_col = self._raw_pos_columns
                rows = rows.loc[
                    rows["chr"].map(lambda value: normalize_chromosome(value, context="raw parquet R2 query")) == self.chrom
                ].copy()
                keep = (
                    pd.to_numeric(rows[left_pos_col], errors="raise").between(int(pos_min), int(pos_max))
                    & pd.to_numeric(rows[right_pos_col], errors="raise").between(int(pos_min), int(pos_max))
                )
                rows = rows.loc[keep].reset_index(drop=True)
                if len(rows) > 0:
                    rows = canonicalize_r2_pairs(rows, self.genome_build)
            rows["R2"] = self._transform_r2(pd.to_numeric(rows["R2"], errors="raise").to_numpy(dtype=np.float32))
            if self.identifier_mode == "rsid":
                rows["i"] = rows["rsID_1"].astype(str).map(self.index_map)
                rows["j"] = rows["rsID_2"].astype(str).map(self.index_map)
            else:
                rows["i"] = pd.to_numeric(rows["pos_1"], errors="raise").astype(np.int64).map(self.index_map)
                rows["j"] = pd.to_numeric(rows["pos_2"], errors="raise").astype(np.int64).map(self.index_map)
            rows = rows.dropna(subset=["i", "j"]).copy()
            rows["i"] = rows["i"].astype(np.int64)
            rows["j"] = rows["j"].astype(np.int64)

        self._last_query_key = key
        self._last_query_rows = rows.copy()
        return rows

    def _deduplicate_pairs(self, pair_rows: pd.DataFrame, context: str) -> pd.DataFrame:
        """Drop duplicate unordered SNP pairs after verifying matching R2 values."""
        if len(pair_rows) == 0:
            return pair_rows
        pair_rows = pair_rows.copy()
        pair_rows["lo"] = np.minimum(pair_rows["i"], pair_rows["j"])
        pair_rows["hi"] = np.maximum(pair_rows["i"], pair_rows["j"])
        duplicated = pair_rows.duplicated(subset=["lo", "hi"], keep=False)
        if not duplicated.any():
            return pair_rows

        dup_df = pair_rows.loc[duplicated, ["lo", "hi", "R2"]]
        nunique = dup_df.groupby(["lo", "hi"], dropna=False)["R2"].nunique(dropna=False)
        if (nunique > 1).any():
            raise ValueError(f"Duplicate unordered SNP pairs with conflicting R2 values detected in {context}.")
        return pair_rows.drop_duplicates(subset=["lo", "hi"], keep="first").reset_index(drop=True)

    def cross_block_matrix(self, l_A: int, b: int, l_B: int, c: int) -> np.ndarray:
        """Build the dense cross-block R2 matrix used by the parquet backend."""
        if b <= 0 or c <= 0:
            return np.zeros((max(b, 0), max(c, 0)), dtype=np.float32)

        a_start, a_stop = l_A, l_A + b
        b_start, b_stop = l_B, l_B + c
        union_min = int(self.pos[min(a_start, b_start)])
        union_max = int(self.pos[max(a_stop - 1, b_stop - 1)])
        rows = self._deduplicate_pairs(
            self._query_union_rows(union_min, union_max),
            context=f"cross-block query {self.chrom}:{l_A}:{b}:{l_B}:{c}",
        )

        matrix = np.zeros((b, c), dtype=np.float32)
        if len(rows) > 0:
            i = rows["i"].to_numpy(dtype=np.int64)
            j = rows["j"].to_numpy(dtype=np.int64)
            values = rows["R2"].to_numpy(dtype=np.float32)

            left_i = (a_start <= i) & (i < a_stop) & (b_start <= j) & (j < b_stop)
            if left_i.any():
                matrix[i[left_i] - a_start, j[left_i] - b_start] = values[left_i]

            left_j = (a_start <= j) & (j < a_stop) & (b_start <= i) & (i < b_stop)
            if left_j.any():
                matrix[j[left_j] - a_start, i[left_j] - b_start] = values[left_j]

        overlap_start = max(a_start, b_start)
        overlap_stop = min(a_stop, b_stop)
        for idx in range(overlap_start, overlap_stop):
            matrix[idx - a_start, idx - b_start] = 1.0

        return matrix

    def within_block_matrix(self, l_B: int, c: int) -> np.ndarray:
        """Build the dense within-block R2 matrix used by the parquet backend."""
        if c <= 0:
            return np.zeros((0, 0), dtype=np.float32)

        b_start, b_stop = l_B, l_B + c
        pos_min = int(self.pos[b_start])
        pos_max = int(self.pos[b_stop - 1])
        rows = self._deduplicate_pairs(
            self._query_union_rows(pos_min, pos_max),
            context=f"within-block query {self.chrom}:{l_B}:{c}",
        )

        matrix = np.zeros((c, c), dtype=np.float32)
        if len(rows) > 0:
            i = rows["i"].to_numpy(dtype=np.int64)
            j = rows["j"].to_numpy(dtype=np.int64)
            values = rows["R2"].to_numpy(dtype=np.float32)
            keep = (b_start <= i) & (i < b_stop) & (b_start <= j) & (j < b_stop)
            if keep.any():
                i = i[keep] - b_start
                j = j[keep] - b_start
                values = values[keep]
                matrix[i, j] = values
                matrix[j, i] = values

        np.fill_diagonal(matrix, 1.0)
        return matrix


def ld_score_var_blocks_from_r2_reader(
    block_left: np.ndarray,
    c: int,
    annot: np.ndarray,
    block_reader: SortedR2BlockReader,
) -> np.ndarray:
    """
    Mirror the old LDSC sliding-block accumulation while sourcing block-local
    dense R2 matrices from a sorted parquet reader instead of genotype blocks.
    """
    m = annot.shape[0]
    n_a = annot.shape[1]
    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / c) * c
    cor_sum = np.zeros((m, n_a), dtype=np.float64)

    b = np.nonzero(block_left > 0)
    if np.any(b):
        b = b[0][0]
    else:
        b = m
    b = int(np.ceil(b / c) * c)
    if b > m:
        c = 1
        b = m

    l_A = 0
    for l_B in range(0, b, c):
        rfuncAB = block_reader.cross_block_matrix(l_A, b, l_B, c)
        cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + c, :])

    b0 = b
    md = int(c * np.floor(m / c))
    end = md + 1 if md != m else md
    for l_B in range(b0, end, c):
        old_b = b
        b = int(block_sizes[l_B])
        if l_B > b0 and b > 0:
            l_A += old_b - b + c
        elif l_B == b0 and b > 0:
            l_A = b0 - b
        elif b == 0:
            l_A = l_B

        chunk_width = c
        if l_B == md:
            chunk_width = m - md

        p1 = np.all(annot[l_A:l_A + b, :] == 0)
        p2 = np.all(annot[l_B:l_B + chunk_width, :] == 0)
        if p1 and p2:
            continue

        rfuncAB = block_reader.cross_block_matrix(l_A, b, l_B, chunk_width)
        cor_sum[l_A:l_A + b, :] += np.dot(rfuncAB, annot[l_B:l_B + chunk_width, :])
        cor_sum[l_B:l_B + chunk_width, :] += np.dot(annot[l_A:l_A + b, :].T, rfuncAB).T

        rfuncBB = block_reader.within_block_matrix(l_B, chunk_width)
        cor_sum[l_B:l_B + chunk_width, :] += np.dot(rfuncBB, annot[l_B:l_B + chunk_width, :])

    return np.asarray(cor_sum, dtype=np.float32)


def compute_counts(metadata: pd.DataFrame, annotations: pd.DataFrame) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Compute LDSC-style annotation count vectors ``M`` and optional ``M_5_50``.

    ``M`` is the column-wise sum over the retained reference SNP universe.
    ``M_5_50`` is the same sum restricted to rows with ``MAF > 0.05`` when MAF
    metadata is available.
    """
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=False)
    M = np.asarray(annot_matrix.sum(axis=0), dtype=np.float64)
    if "MAF" not in metadata.columns or metadata["MAF"].isna().all():
        return M, None
    common = metadata["MAF"] > 0.05
    M_5_50 = np.asarray(annot_matrix[common.to_numpy(), :].sum(axis=0), dtype=np.float64)
    return M, M_5_50


def regression_mask_from_keys(metadata: pd.DataFrame, regression_keys: set[str] | None, identifier_mode: str) -> np.ndarray:
    """Build the binary mask column used to compute regression-weight LD scores."""
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
    Compute all LD-score outputs for one chromosome from sorted parquet R2 input.

    Main steps:
    1. Merge optional CM/MAF metadata and apply any MAF filter.
    2. Intersect annotation SNPs with SNPs actually present in the R2 table.
    3. Validate the LD window against the retained SNP coordinates.
    4. Query block-local dense R2 matrices from the sorted parquet file while
       following the old LDSC sliding-block traversal.
    5. Return chromosome-level LD scores plus M and M_5_50 counts.
    """
    metadata = merge_frequency_metadata(bundle.metadata.copy(), args, chrom=chrom, identifier_mode=args.snp_identifier)
    if args.snp_identifier == "chr_pos" and getattr(args, "genome_build", None) == "auto":
        metadata, inference = resolve_chr_pos_table(
            metadata,
            context=f"parquet metadata chromosome {chrom}",
            reference_table=load_packaged_reference_table(),
            logger=LOGGER,
        )
        args.genome_build = inference.genome_build
    metadata, annotations = apply_maf_filter(metadata, bundle.annotations.copy(), args.maf, context="parquet mode")
    present_values = read_sorted_r2_presence(args, chrom=chrom)
    metadata, annotations, regression_keys = filter_reference_to_present_r2(
        metadata, annotations, regression_keys, present_values, args.snp_identifier, chrom
    )
    if len(metadata) == 0:
        raise ValueError(f"No retained annotation SNPs remain on chromosome {chrom} after parquet intersection.")

    validate_retained_identifier_uniqueness(metadata, args.snp_identifier, chrom)
    coords, max_dist = build_window_coordinates(metadata, args)
    block_left = get_block_lefts(
        coords,
        max_dist,
    )
    check_whole_chromosome_window(block_left, args, chrom)
    regression_mask = regression_mask_from_keys(metadata, regression_keys, args.snp_identifier).reshape(-1, 1)
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=True)
    combined_annot = np.c_[annot_matrix, regression_mask]
    block_reader = SortedR2BlockReader(
        paths=resolve_parquet_files(args, chrom=chrom),
        chrom=chrom,
        metadata=metadata,
        identifier_mode=args.snp_identifier,
        r2_bias_mode=args.r2_bias_mode,
        r2_sample_size=args.r2_sample_size,
        genome_build=args.genome_build,
    )
    combined_scores = ld_score_var_blocks_from_r2_reader(
        block_left=block_left,
        c=args.chunk_size,
        annot=combined_annot,
        block_reader=block_reader,
    )
    ld_scores = combined_scores[:, :annot_matrix.shape[1]]
    w_ld = combined_scores[:, annot_matrix.shape[1]:]
    out_metadata = metadata.reset_index(drop=True)
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
    panel_df = bim.df.loc[:, ["CHR", "SNP", "CM", "BP"]].copy().rename(columns={"BP": "POS"})
    panel_df["CHR"] = panel_df["CHR"].map(lambda value: normalize_chromosome(value, context=prefix + ".bim"))
    panel_df["SNP"] = panel_df["SNP"].astype(str)
    panel_df["POS"] = pd.to_numeric(panel_df["POS"], errors="raise").astype(np.int64)
    panel_df["CM"] = pd.to_numeric(panel_df["CM"], errors="coerce")
    panel_df = panel_df.loc[panel_df["CHR"] == normalize_chromosome(chrom, context=prefix + ".bim")].copy()
    if len(panel_df) == 0:
        raise ValueError(f"No PLINK SNPs found for chromosome {chrom} in {prefix}.")
    panel_df["_key"] = identifier_keys(panel_df, args.snp_identifier)

    metadata = bundle.metadata.copy()
    annotations = bundle.annotations.copy()
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

    keep_indivs = resolve_keep_individuals(getattr(args, "keep", None), fam)
    keep_snps = [key_to_panel_index[key] for key in metadata["_key"]]
    geno = legacy_ld.PlinkBEDFile(
        prefix + ".bed",
        len(fam.IDList),
        bim,
        keep_snps=keep_snps,
        keep_indivs=keep_indivs,
        mafMin=args.maf,
    )

    geno_meta = pd.DataFrame(geno.df, columns=geno.colnames)
    if "BP" in geno_meta.columns:
        geno_meta = geno_meta.rename(columns={"BP": "POS"})
    geno_meta["CHR"] = geno_meta["CHR"].map(lambda value: normalize_chromosome(value, context=prefix + ".bim"))
    geno_meta["SNP"] = geno_meta["SNP"].astype(str)
    geno_meta["POS"] = pd.to_numeric(geno_meta["POS"], errors="raise").astype(np.int64)
    geno_meta["CM"] = pd.to_numeric(geno_meta["CM"], errors="coerce")
    geno_meta["MAF"] = pd.to_numeric(geno_meta["MAF"], errors="coerce")
    geno_meta["_key"] = identifier_keys(geno_meta, args.snp_identifier)

    annotation_matrix = annotations.set_index(metadata["_key"]).loc[geno_meta["_key"]]

    coords, max_dist = build_window_coordinates(geno_meta.drop(columns="_key"), args)
    block_left = legacy_ld.getBlockLefts(coords, max_dist)
    check_whole_chromosome_window(block_left, args, chrom)

    ld_scores = geno.ldScoreVarBlocks(block_left, args.chunk_size, annot=annotation_matrix.to_numpy(dtype=np.float32))
    regression_mask = regression_mask_from_keys(geno_meta.drop(columns="_key"), regression_keys, args.snp_identifier)
    geno._currentSNP = 0
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


def resolve_keep_individuals(keep_path: str | None, fam) -> list[int] | None:
    """Return FAM row indices retained by ``keep_path`` using IID matching."""
    if not keep_path:
        return None
    keep_indivs = fam.loj(legacy_parse.FilterFile(keep_path).IDList)
    if len(keep_indivs) == 0:
        raise ValueError("No individuals retained for analysis")
    return keep_indivs.tolist()


def result_to_dataframe(result: ChromComputationResult) -> pd.DataFrame:
    """Materialize one chromosome result as a standard ``.l2.ldscore`` table."""
    df = result.metadata.copy()
    for idx, column in enumerate(result.ldscore_columns):
        df[column + "L2"] = result.ld_scores[:, idx]
    return df


def weight_result_to_dataframe(result: ChromComputationResult) -> pd.DataFrame:
    """Materialize one chromosome regression-weight result as a one-column table."""
    df = result.metadata.copy()
    df["L2"] = np.ravel(result.w_ld)
    return df


# Output assembly.
def write_ldscore_file(df: pd.DataFrame, path: str) -> None:
    """Write one LDSC-compatible LD-score table, preserving metadata columns first."""
    out = df.copy()
    out = out.loc[:, [col for col in ["CHR", "SNP", "POS", "CM", "MAF"] if col in out.columns] + [col for col in out.columns if col not in ANNOT_META_COLUMNS]]
    with gzip.open(path, "wt") as handle:
        out.to_csv(handle, sep="\t", index=False, na_rep="NA", float_format="%.6g")


def write_counts(path: str, counts: np.ndarray) -> None:
    """Write one LDSC ``.M``-style count vector to disk."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(str(x) for x in counts))


def write_annotation_groups(path: str, baseline_columns: Sequence[str], query_columns: Sequence[str]) -> None:
    """Write the baseline/query annotation-group manifest used by downstream tools."""
    rows = [{"annotation": col, "group": "baseline"} for col in baseline_columns]
    rows.extend({"annotation": col, "group": "query"} for col in query_columns)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def aggregate_results(results: Sequence[ChromComputationResult]) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray, np.ndarray | None]:
    """
    Combine chromosome-level results into one genome-wide LD-score output set.

    Returns the aggregated reference LD-score table, regression-weight table,
    and the chromosome-summed ``M`` / ``M_5_50`` vectors.
    """
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
    """Validate top-level LD-score CLI arguments before any heavy work starts."""
    args.snp_identifier = normalize_snp_identifier_mode(args.snp_identifier)
    args.genome_build = normalize_genome_build(args.genome_build)
    validate_auto_genome_build_mode(args.snp_identifier, args.genome_build)
    keep = getattr(args, "keep", None)
    if not (args.query_annot or args.query_annot_chr or args.baseline_annot or args.baseline_annot_chr):
        raise ValueError("At least one of query or baseline annotation inputs must be supplied.")
    if bool(args.r2_table or args.r2_table_chr) == bool(args.bfile or args.bfile_chr):
        raise ValueError("Specify exactly one reference-panel mode: parquet or PLINK.")
    if args.r2_table_chr and args.r2_table:
        raise ValueError("Cannot set both --r2-table and --r2-table-chr.")
    if args.bfile_chr and args.bfile:
        raise ValueError("Cannot set both --bfile and --bfile-chr.")
    if args.r2_table or args.r2_table_chr:
        if keep:
            raise ValueError("--keep is only supported in PLINK mode.")
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
    """Build the standalone LD-score kernel CLI parser."""
    parser = argparse.ArgumentParser(
        description="Estimate LDSC-compatible LD scores from SNP-level annotation files using PLINK or sorted parquet R2 input.",
    )
    parser.add_argument("--out", required=True, help="Output prefix.")
    parser.add_argument("--query-annot", default=None, help="Comma-separated SNP-level query annotation files or prefixes.")
    parser.add_argument("--query-annot-chr", default=None, help="Comma-separated chromosome-pattern prefixes for query annotation files.")
    parser.add_argument("--baseline-annot", default=None, help="Comma-separated SNP-level baseline annotation files or prefixes.")
    parser.add_argument("--baseline-annot-chr", default=None, help="Comma-separated chromosome-pattern prefixes for baseline annotation files.")
    parser.add_argument("--bfile", default=None, help="PLINK prefix for the reference panel.")
    parser.add_argument("--bfile-chr", default=None, help="Chromosome-pattern PLINK prefix for the reference panel.")
    parser.add_argument("--r2-table", default=None, help="Comma-separated sorted parquet R2 files or prefixes.")
    parser.add_argument("--r2-table-chr", default=None, help="Comma-separated chromosome-pattern prefixes for sorted parquet R2 files.")
    parser.add_argument("--snp-identifier", default="chr_pos", help="Identifier mode used to match annotations to the reference panel.")
    parser.add_argument(
        "--genome-build",
        default=None,
        choices=("auto", "hg19", "hg37", "GRCh37", "hg38", "GRCh38"),
        help="Genome build assumed for the sorted parquet R2 file and chr_pos matching. Use 'auto' to infer hg19/hg38 and 0-based/1-based coordinates.",
    )
    parser.add_argument("--r2-bias-mode", choices=("raw", "unbiased"), default=None, help="Whether sorted parquet R2 values are raw sample r^2 or already unbiased.")
    parser.add_argument("--r2-sample-size", default=None, type=float, help="LD reference sample size used to correct raw parquet R2 values.")
    parser.add_argument("--regression-snps", default=None, help="Optional SNP list defining the regression SNP set for w_ld.")
    parser.add_argument("--frqfile", default=None, help="Optional frequency/metadata file for MAF and CM.")
    parser.add_argument("--frqfile-chr", default=None, help="Chromosome-pattern frequency/metadata file prefix.")
    parser.add_argument("--keep", default=None, help="File with individuals to include in LD Score estimation. The file should contain one IID per row.")
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
def run_ldscore_from_args(args: argparse.Namespace) -> list[ChromComputationResult]:
    """
    Execute the end-to-end LD-score workflow from parsed CLI-style arguments.

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


def run_ldscore(
    *,
    out: str,
    query_annot: str | None = None,
    query_annot_chr: str | None = None,
    baseline_annot: str | None = None,
    baseline_annot_chr: str | None = None,
    bfile: str | None = None,
    bfile_chr: str | None = None,
    r2_table: str | None = None,
    r2_table_chr: str | None = None,
    snp_identifier: str = "chr_pos",
    genome_build: str | None = None,
    r2_bias_mode: str | None = None,
    r2_sample_size: float | None = None,
    regression_snps: str | None = None,
    frqfile: str | None = None,
    frqfile_chr: str | None = None,
    keep: str | None = None,
    ld_wind_snps: int | None = None,
    ld_wind_kb: float | None = None,
    ld_wind_cm: float | None = None,
    maf: float | None = None,
    chunk_size: int = 50,
    per_chr_output: bool = False,
    yes_really: bool = False,
    log_level: str = "INFO",
) -> list[ChromComputationResult]:
    """
    Execute the end-to-end LD-score workflow through a normal Python API.

    This is the package-level function entrypoint. It mirrors the CLI options
    with keyword arguments and writes the same outputs as the command-line path.
    """
    args = argparse.Namespace(
        out=out,
        query_annot=query_annot,
        query_annot_chr=query_annot_chr,
        baseline_annot=baseline_annot,
        baseline_annot_chr=baseline_annot_chr,
        bfile=bfile,
        bfile_chr=bfile_chr,
        r2_table=r2_table,
        r2_table_chr=r2_table_chr,
        snp_identifier=snp_identifier,
        genome_build=genome_build,
        r2_bias_mode=r2_bias_mode,
        r2_sample_size=r2_sample_size,
        regression_snps=regression_snps,
        frqfile=frqfile,
        frqfile_chr=frqfile_chr,
        keep=keep,
        ld_wind_snps=ld_wind_snps,
        ld_wind_kb=ld_wind_kb,
        ld_wind_cm=ld_wind_cm,
        maf=maf,
        chunk_size=chunk_size,
        per_chr_output=per_chr_output,
        yes_really=yes_really,
        log_level=log_level,
    )
    configure_logging(log_level)
    return run_ldscore_from_args(args)


def main(argv: Sequence[str] | None = None) -> list[ChromComputationResult]:
    """CLI entrypoint for the standalone LD-score kernel module."""
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.log_level)
    return run_ldscore_from_args(args)


if __name__ == "__main__":
    main()
