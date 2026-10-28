# Includes primitives derived from legacy LDSC ldscore/ldscore.py.
# GPLv3 (LICENSE); upstream authors and Python 3 port credit are in NOTICE.
# Modified 2026 for prepared reference panels, streaming R2, and projection;
# the PLINK LD-score block sums moved here from the genotype reader.
# See NOTICE for attribution and dated modifications.

"""LD-score projection primitives and streaming R2 panel reader.

``RefPanel.prepare_chromosome`` owns source resolution, reference filtering,
annotation alignment, LD-window validation, and reader lifetime. This module's
``compute_chromosome`` consumes that ``PreparedChromosome`` state and projects
baseline/query annotations and regression weights in one traversal. It also
provides annotation parsing, window geometry, and count primitives used by the
workflow adapters. BED-to-SNP projection belongs to ``ldsc._annotation.builder``.

Scientific contracts
--------------------
Prepared metadata, selected annotation rows, and window bounds describe the
same retained reference SNPs in genomic order. Binary annotations can remain
packed until selected reads; continuous values retain float32 normalization.
Projection accumulates in float64. Reference-panel CM and MAF are
authoritative. PLINK sample selection precedes genotype QC and inclusive MAF
filtering; chromosome LD scores retain the existing float32 result conversion.
Parquet projection accumulates in float64. Common-SNP counts and overlap use
``MAF >= common_maf_min`` (0.05 by default), independently of panel filtering.

R2 panel Parquet files contain sidecar-row indices ``IDX_1``/``IDX_2``,
``R2``, and ``SIGN_R``. The required sidecar is bound by its identity digest and
row count. Each chromosome builds one identity remap; pairs stream through
``iter_all_pairs`` without a dense chromosome matrix or row-group cache. Every
file must match the ``build-r2-panel`` metadata contract checked by
``read_r2_parquet_contract``: int16 ``R2`` quantized with scale 32767 and already
bias-corrected, so decoding divides by the scale and applies no further
correction. Each unordered off-diagonal pair contributes in both directions,
and the diagonal is added internally as one.

This module returns in-memory results and writes no artifacts. Public result
normalization and aggregation live in ``ldsc._ldscore.calculator``; canonical
output layout belongs to ``ldsc.outputs``. See
``docs/reference/formats/r2-panels.md`` for the file contract.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .._core.columns import normalize_snp_identifier_mode
from ..errors import LDSCInternalError
from .._annotation.overlap import OverlapContribution, annotation_statistics

from .._core.snp_identity import (
    RestrictionIdentityKeys,
    empty_identity_drop_frame,
    restriction_membership_mask,
)
from .._core.regions import RegionIntervals, region_exclusion_keep_mask
from .._core.identifiers import identifier_keys
from .._refpanel.prepared import ArrayAnnotations, PreparedChromosome
from .._refpanel.r2_format import SortedR2BlockReader


LOGGER = logging.getLogger("LDSC.ldscore")
# K bounds buffered pairs; CSR construction and projection add working arrays.
_CSR_CHUNK_PAIRS = 16_000_000


@dataclass
class ChromComputationResult:
    """Chromosome-level LD-score outputs before cross-chromosome aggregation."""
    chrom: str
    metadata: pd.DataFrame
    ld_scores: np.ndarray
    w_ld: np.ndarray
    M: np.ndarray
    M_5_50: np.ndarray
    ldscore_columns: list[str]
    baseline_columns: list[str]
    query_columns: list[str]
    reference_snp_count: int = 0
    regression_selected_snp_count: int = 0
    regression_region_removed_snp_count: int = 0
    annotation_types: dict = field(default_factory=dict)
    overlap: OverlapContribution | None = None
    identity_drops: pd.DataFrame = field(default_factory=empty_identity_drop_frame)


# Shared computational helpers.


def ld_score_streaming_from_r2_reader(
    block_left: np.ndarray,
    annot,
    block_reader: SortedR2BlockReader,
    chunk_pairs: int = _CSR_CHUNK_PAIRS,
    *, output_rows=None, n_baseline=None, query_batch_size=100, weight_mask=None,
) -> np.ndarray:
    """Project selected annotation batches through each stored R2 chunk once.

    Only output SNPs own score buffers, but both pair directions retain every
    eligible reference SNP. Diagonals and the separate regression
    weight mask use the same float64 accumulator; direct results are float32.
    """
    block_left = np.asarray(block_left, dtype=np.int64)
    if isinstance(annot, np.ndarray):
        annot = ArrayAnnotations(annot, tuple(str(i) for i in range(annot.shape[1])))
    accumulator = ProjectionAccumulator(annot,
        n_baseline=annot.shape[1] if n_baseline is None else n_baseline,
        output_rows=output_rows, query_batch_size=query_batch_size, weight_mask=weight_mask)
    accumulator.add_diagonal()
    buf_i, buf_j, buf_r2 = [], [], []
    buffered = 0

    def flush():
        nonlocal buffered
        if buffered:
            accumulator.add_pairs(np.concatenate(buf_i), np.concatenate(buf_j),
                                  np.concatenate(buf_r2), block_left)
            buf_i.clear()
            buf_j.clear()
            buf_r2.clear()
            buffered = 0

    for i, j, r2 in block_reader.iter_all_pairs():
        start = 0
        while start < len(i):
            stop = min(len(i), start + chunk_pairs-buffered)
            buf_i.append(i[start:stop])
            buf_j.append(j[start:stop])
            buf_r2.append(r2[start:stop])
            buffered += stop-start
            start = stop
            if buffered == chunk_pairs:
                flush()
    flush()
    return np.asarray(accumulator.values, dtype=np.float32)


def ld_scores_from_genotypes(reader, block_left, chunk_size, annot=None, *, output_rows=None,
                             n_baseline=None, query_batch_size=100, weight_mask=None):
    """Project shared and batched annotations in one traversal of PLINK genotypes.

    ``reader`` is a PLINK genotype reader positioned at its first SNP; its
    ``nextSNPs`` supplies standardized genotypes for the ``m`` retained SNPs
    and ``n`` individuals. ``annot`` is a selected-read source or a
    caller-owned numerical array. Only ``output_rows`` receive float64 score
    buffers; all aligned reference SNPs still contribute. The optional weight
    mask selects its own SNP set. LD-score sums use exactly the pairs inside
    ``block_left``; unlike legacy LDSC, the chunk size does not widen the
    window. Legacy LDSC computes the same sums in the genotype reader's
    ``ldScoreVarBlocks`` method (``ldscore/ldscore.py``).
    """
    n = reader.n
    return _cor_sum_var_blocks(reader.m, n, block_left, chunk_size, lambda x: _l2_unbiased(x, n),
        reader.nextSNPs, annot, output_rows=output_rows, n_baseline=n_baseline,
        query_batch_size=query_batch_size, weight_mask=weight_mask)


def _l2_unbiased(x, n):
    """Convert correlation values to the unbiased LD-score contribution."""
    denom = n - 2 if n > 2 else n
    sq = np.square(x)
    return sq - (1 - sq) / denom


def _cor_sum_var_blocks(m, n, block_left, c, func, snp_getter, annot=None, *,
                        output_rows=None, n_baseline=None, query_batch_size=100,
                        weight_mask=None):
    """Accumulate transformed correlation sums over the exact LD window.

    The legacy LDSC traversal below rounds each block up to a multiple of ``c``
    and shares the chunk start's block across the chunk, so it computes a
    superset of the window pairs, each exactly once. Pairs outside the
    window are zeroed before accumulation, which makes the result
    independent of ``c``.
    """
    window_left = np.asarray(block_left, dtype=np.int64)
    block_sizes = np.array(np.arange(m) - block_left)
    block_sizes = np.ceil(block_sizes / c) * c
    if annot is None:
        annot = np.ones((m, 1))
    else:
        annot_m = annot.shape[0]
        if annot_m != m:
            raise LDSCInternalError(
                "PLINK LD-score calculation received an annotation matrix with the "
                f"wrong SNP count ({annot_m} rows for {m} retained SNPs). Most "
                "likely annotation/reference alignment desynchronized before PLINK "
                "computation. Re-run with DEBUG logging and report the traceback."
            )

    if isinstance(annot, np.ndarray):
        annot = ArrayAnnotations(annot, tuple(str(i) for i in range(annot.shape[1])))
    accumulator = ProjectionAccumulator(annot,
        n_baseline=annot.shape[1] if n_baseline is None else n_baseline,
        output_rows=output_rows, query_batch_size=query_batch_size, weight_mask=weight_mask)
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
        _zero_pairs_outside_window(rfuncAB, l_A, l_B, window_left)
        accumulator.add_dense(np.arange(l_A, l_A+b), np.arange(l_B, l_B+c), rfuncAB)
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
        np.dot(A.T, B / n, out=rfuncAB)
        rfuncAB = func(rfuncAB)
        _zero_pairs_outside_window(rfuncAB, l_A, l_B, window_left)
        accumulator.add_dense(np.arange(l_A, l_A+b), np.arange(l_B, l_B+c), rfuncAB)
        accumulator.add_dense(np.arange(l_B, l_B+c), np.arange(l_A, l_A+b), rfuncAB.T)
        np.dot(B.T, B / n, out=rfuncBB)
        rfuncBB = func(rfuncBB)
        _zero_pairs_outside_window(rfuncBB, l_B, l_B, window_left)
        accumulator.add_dense(np.arange(l_B, l_B+c), np.arange(l_B, l_B+c), rfuncBB)
    return accumulator.values


def _zero_pairs_outside_window(values, row_start, col_start, block_left):
    """Zero, in place, the entries of one dense block outside the exact LD window.

    ``values[p, q]`` is the pair ``(row_start + p, col_start + q)``, which lies
    in the window iff ``min(i, j) >= block_left[max(i, j)]``. Because
    ``block_left`` is non-decreasing, each column keeps one contiguous row
    range, so only per-column bounds are computed (no block-sized mask).
    """
    n_rows, n_cols = values.shape
    if n_rows == 0 or n_cols == 0:
        return
    cols = np.arange(col_start, col_start + n_cols)
    first = block_left[cols] - row_start
    stop = np.searchsorted(block_left, cols, side="right") - row_start
    for q in np.flatnonzero((first > 0) | (stop < n_rows)):
        values[:max(first[q], 0), q] = 0.0
        values[max(stop[q], 0):, q] = 0.0


def compute_counts(
    metadata: pd.DataFrame,
    annotations: pd.DataFrame,
    common_maf_min: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute LDSC-style annotation count vectors ``M`` and common-SNP counts.

    ``M`` is the column-wise sum over the retained reference SNPs.
    The common-SNP count vector is the same sum restricted to rows with
    ``MAF >= common_maf_min`` (inclusive; deviates from legacy LDSC's strict
    ``0.05 < FRQ < 0.95`` on canonical folded MAF). ``metadata`` must carry
    MAF for every row; both reference-panel adapters guarantee it.
    """
    if "MAF" not in metadata.columns or metadata["MAF"].isna().any():
        raise LDSCInternalError(
            "LD-score counts received reference metadata without MAF for every SNP, so the "
            "common-SNP counts cannot be computed. Most likely a custom reference-panel adapter "
            "dropped the MAF column or returned missing values. Supply MAF for every retained "
            "SNP; for the built-in panels, re-run with DEBUG logging and report the traceback."
        )
    annot_matrix = annotations.to_numpy(dtype=np.float32, copy=False)
    M = np.asarray(annot_matrix.sum(axis=0), dtype=np.float64)
    common = metadata["MAF"] >= common_maf_min
    M_5_50 = np.asarray(annot_matrix[common.to_numpy(), :].sum(axis=0), dtype=np.float64)
    return M, M_5_50


def regression_mask_from_keys(
    metadata: pd.DataFrame,
    regression_keys: set[str] | RestrictionIdentityKeys | None,
    identifier_mode: str,
    region_intervals: RegionIntervals | None = None,
) -> np.ndarray:
    """Build the selected-regression-set mask after named region subtraction.

    This is deliberately separate from the reference SNPs: callers use the
    returned mask only for written regression SNPs and the SNPs that contribute to ``w_ld``.
    """
    if regression_keys is None:
        keep = np.ones(len(metadata), dtype=bool)
    elif isinstance(regression_keys, RestrictionIdentityKeys):
        keep = restriction_membership_mask(
            metadata,
            regression_keys,
            identifier_mode,
            context="regression-SNP restriction matching for LD scores",
        ).to_numpy(dtype=bool)
    else:
        keys = identifier_keys(metadata, identifier_mode)
        keep = keys.isin(regression_keys).to_numpy(dtype=bool)
    if region_intervals is not None and region_intervals.intervals:
        keep &= np.asarray(region_exclusion_keep_mask(metadata, region_intervals), dtype=bool)
    return keep.astype(np.float32)


def compute_chromosome(
    chrom: str, prepared: PreparedChromosome, *, snp_identifier: str,
    snp_batch_size: int, common_maf_min: float = 0.05, query_batch_size: int = 100,
    regression_keys: set[str] | RestrictionIdentityKeys | None = None,
    regression_regions: RegionIntervals | None = None,
) -> ChromComputationResult:
    """Project annotations and regression weights through one prepared chromosome.

    Preparation owns identity, sample/MAF selection, reader policy, and window
    validation. This function only projects the aligned matrices and derives
    counts and overlap; it never resolves paths or reopens a reference panel.
    """
    metadata = prepared.metadata
    regression_selected = regression_mask_from_keys(metadata, regression_keys, snp_identifier)
    regression_mask = regression_mask_from_keys(
        metadata, regression_keys, snp_identifier, region_intervals=regression_regions,
    )
    output_rows = np.flatnonzero(regression_mask)
    options = dict(output_rows=output_rows, n_baseline=len(prepared.baseline_columns),
                   query_batch_size=query_batch_size, weight_mask=regression_mask)
    if prepared.backend == "plink":
        prepared.reader._currentSNP = 0
        combined_scores = np.asarray(ld_scores_from_genotypes(
            prepared.reader, prepared.block_left, snp_batch_size, annot=prepared.annotations, **options,
        ), dtype=np.float32)
    else:
        combined_scores = ld_score_streaming_from_r2_reader(
            block_left=prepared.block_left, annot=prepared.annotations, block_reader=prepared.reader, **options,
        )
    n_columns = prepared.annotations.shape[1]
    M, M_5_50, overlap, kinds = annotation_statistics(metadata, prepared.annotations,
        len(prepared.baseline_columns), common_maf_min=common_maf_min, query_batch_size=query_batch_size)
    return ChromComputationResult(
        chrom=chrom, metadata=metadata.iloc[output_rows].reset_index(drop=True),
        ld_scores=combined_scores[:, :n_columns], w_ld=combined_scores[:, n_columns:],
        M=M, M_5_50=M_5_50,
        ldscore_columns=prepared.baseline_columns + prepared.query_columns,
        baseline_columns=prepared.baseline_columns, query_columns=prepared.query_columns,
        reference_snp_count=len(metadata), regression_selected_snp_count=int(regression_selected.sum()),
        regression_region_removed_snp_count=int(regression_selected.sum()-len(output_rows)),
        annotation_types=kinds, overlap=overlap, identity_drops=prepared.identity_drops,
    )


class ProjectionAccumulator:
    """Accumulate selected output rows while preserving all reference SNPs.

    ``annotations`` is a borrowed chromosome accessor with shape
    ``(n_reference_rows, n_annotation_columns)`` and ordered ``columns``.
    ``output_rows`` selects accumulator rows only; it never filters the
    reference SNPs whose annotations contribute. ``values`` owns a float64 array with shape
    ``(len(output_rows), n_annotation_columns + (weight_mask is not None))``.

    Shared baseline columns and the optional weight are projected together;
    ``query_batch_size`` limits query columns per operand read, not the total
    width of ``values``. The workflow supplies only its active query batch.
    Each float64 reference-SNP tile is reused across output blocks.
    ``read_budget_bytes`` targets one operand tile (at least one row), not
    total RAM: sparse operators, slices, gathers, and products also allocate.
    The caller retains ownership of the accessor and controls public output
    conversion after float64 accumulation.
    """

    def __init__(self, annotations, *, n_baseline, output_rows=None, weight_mask=None, query_batch_size=100, read_budget_bytes=16*1024*1024):
        if isinstance(query_batch_size, bool) or not isinstance(query_batch_size, (int, np.integer)) or query_batch_size < 1:
            raise ValueError('query_batch_size must be a positive integer.')
        self.annotations = annotations
        self.n_baseline = n_baseline
        self.query_batch_size = int(query_batch_size)
        self.read_budget_bytes = read_budget_bytes
        self.weight_mask = weight_mask
        self.output_rows = np.arange(annotations.shape[0]) if output_rows is None else np.asarray(output_rows, dtype=np.int64)
        self.output_index = np.full(annotations.shape[0], -1, dtype=np.int64)
        self.output_index[self.output_rows] = np.arange(len(self.output_rows))
        self.values = np.zeros((len(self.output_rows), annotations.shape[1] + (weight_mask is not None)), dtype=np.float64)

    def batches(self):
        shared = list(range(self.n_baseline))
        if self.weight_mask is not None:
            shared.append(self.annotations.shape[1])
        if shared:
            yield shared
        for start in range(self.n_baseline, self.annotations.shape[1], self.query_batch_size):
            yield list(range(start, min(start+self.query_batch_size, self.annotations.shape[1])))

    def _read(self, rows, columns):
        names = [self.annotations.columns[index] for index in columns if index < self.annotations.shape[1]]
        result = np.empty((len(rows), len(columns)), dtype=np.float64)
        if names:
            result[:, :len(names)] = self.annotations.read(rows=rows, columns=names)
        if len(names) != len(columns):
            result[:, -1] = self.weight_mask[rows]
        return result

    def add_diagonal(self):
        """Seed R2=1 contributions on output rows, including selected weights."""
        for columns in self.batches():
            tile = max(1, self.read_budget_bytes // (8*len(columns)))
            for start in range(0, len(self.output_rows), tile):
                stop = min(start+tile, len(self.output_rows))
                self.values[np.ix_(np.arange(start, stop), columns)] += self._read(self.output_rows[start:stop], columns)

    def _accumulate(self, operator, output_indices, contributors):
        if not len(output_indices) or not len(contributors):
            return
        sparse = hasattr(operator, 'nnz')
        for columns in self.batches():
            tile = max(1, self.read_budget_bytes // (8*len(columns)))
            for start in range(0, len(contributors), tile):
                stop = min(start+tile, len(contributors))
                values = self._read(contributors[start:stop], columns)
                for first in range(0, len(output_indices), tile):
                    last = min(first+tile, len(output_indices))
                    block = operator[first:last, start:stop]
                    if sparse and not block.nnz:
                        continue
                    self.values[np.ix_(output_indices[first:last], columns)] += block @ values

    def add_dense(self, targets, contributors, correlations):
        """Project one already-computed dense block without narrowing ``contributors``."""
        targets, contributors = np.asarray(targets, dtype=np.int64), np.asarray(contributors, dtype=np.int64)
        indices = self.output_index[targets]
        keep = indices >= 0
        self._accumulate(correlations[keep], indices[keep], contributors)

    def add_pairs(self, i, j, r2, block_left):
        """Project both directions of upper-triangle pairs through one sparse block."""
        from scipy import sparse

        keep = i >= block_left[j]
        i, j, r2 = i[keep], j[keep], r2[keep]
        left, right = self.output_index[i], self.output_index[j]
        first, second = left >= 0, right >= 0
        targets = np.concatenate([left[first], right[second]])
        contributors = np.concatenate([j[first], i[second]])
        if not len(targets):
            return
        output_indices, rows = np.unique(targets, return_inverse=True)
        source_indices, columns = np.unique(contributors, return_inverse=True)
        data = np.concatenate([r2[first], r2[second]]).astype(np.float64)
        operator = sparse.csr_matrix((data, (rows, columns)), shape=(len(output_indices), len(source_indices)))
        self._accumulate(operator, output_indices, source_indices)
