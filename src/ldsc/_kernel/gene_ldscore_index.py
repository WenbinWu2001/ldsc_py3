"""Pure sparse primitives for exact disjoint-atom gene LD-score indexes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence

import numpy as np
from scipy import sparse

from ..errors import LDSCInternalError


@dataclass(frozen=True)
class ChromosomeAtomModel:
    """Deterministic chromosome-local atom geometry and gene membership."""

    chromosome: str
    starts: np.ndarray
    ends: np.ndarray
    gene_to_atom: sparse.csr_matrix

    @property
    def n_atoms(self) -> int:
        """Return the number of disjoint atoms."""
        return int(len(self.starts))

    def validate(self) -> None:
        """Validate geometry, Boolean CSR membership, and excluded empty rows."""
        starts = np.asarray(self.starts)
        ends = np.asarray(self.ends)
        matrix = self.gene_to_atom
        if starts.dtype != np.int64 or ends.dtype != np.int64:
            raise LDSCInternalError("Gene LD-score atom coordinates must use int64 dtype.")
        if starts.shape != ends.shape or starts.ndim != 1:
            raise LDSCInternalError("Gene LD-score atom coordinate arrays must be aligned one-dimensional arrays.")
        if len(starts) and (np.any(starts < 0) or np.any(ends <= starts)):
            raise LDSCInternalError("Gene LD-score atoms must be nonempty non-negative half-open intervals.")
        if len(starts) > 1 and (np.any(starts[1:] < ends[:-1]) or np.any(starts[1:] < starts[:-1])):
            raise LDSCInternalError("Gene LD-score atoms must be sorted and nonoverlapping.")
        if not sparse.isspmatrix_csr(matrix) or matrix.shape[1] != len(starts):
            raise LDSCInternalError("Gene-to-atom membership must be CSR with one column per atom.")
        if matrix.dtype != np.bool_ or np.any(matrix.data != 1):
            raise LDSCInternalError("Gene-to-atom membership must contain Boolean true values only.")
        if matrix.indices.dtype != np.int32 or matrix.indptr.dtype != np.int32:
            raise LDSCInternalError("Gene-to-atom CSR indices must use int32 dtype.")
        if not matrix.has_sorted_indices:
            raise LDSCInternalError("Gene-to-atom CSR column indices must be sorted within rows.")
        if not matrix.has_canonical_format:
            raise LDSCInternalError("Gene-to-atom CSR must not contain duplicate entries.")


@dataclass(frozen=True)
class AtomStatistics:
    """All/common atom counts and ordered baseline-to-atom overlaps."""

    atom_count_all: np.ndarray
    atom_count_common: np.ndarray
    baseline_atom_overlap_all: np.ndarray
    baseline_atom_overlap_common: np.ndarray


@dataclass(frozen=True)
class SelectedAtomStatistics:
    """Counts and overlaps assembled for one Boolean selected-atom union."""

    count_all: int
    count_common: int
    baseline_overlap_all: np.ndarray
    baseline_overlap_common: np.ndarray
    control_overlap_all: int | None
    control_overlap_common: int | None


def build_disjoint_atoms(
    chromosome: str,
    intervals: np.ndarray,
    *,
    included: np.ndarray | None = None,
    padding_bp: int = 0,
) -> ChromosomeAtomModel:
    """Build maximal nonempty half-open atoms with constant covering-gene sets."""
    raw = np.asarray(intervals, dtype=np.int64)
    if raw.ndim != 2 or raw.shape[1] != 2:
        raise LDSCInternalError("Gene intervals must have shape (n_genes, 2).")
    n_genes = int(raw.shape[0])
    include = np.ones(n_genes, dtype=bool) if included is None else np.asarray(included, dtype=bool)
    if include.shape != (n_genes,):
        raise LDSCInternalError("Gene inclusion mask must have one value per interval.")
    if isinstance(padding_bp, bool) or int(padding_bp) != padding_bp or padding_bp < 0:
        raise LDSCInternalError("Gene atom padding must be a non-negative integer.")
    if n_genes and (np.any(raw[:, 0] < 0) or np.any(raw[:, 1] <= raw[:, 0])):
        raise LDSCInternalError("Gene intervals must be nonempty non-negative half-open intervals.")

    additions: dict[int, list[int]] = {}
    removals: dict[int, list[int]] = {}
    for gene_row in np.flatnonzero(include):
        start = max(0, int(raw[gene_row, 0]) - int(padding_bp))
        end = int(raw[gene_row, 1]) + int(padding_bp)
        additions.setdefault(start, []).append(int(gene_row))
        removals.setdefault(end, []).append(int(gene_row))

    boundaries = sorted(set(additions) | set(removals))
    starts: list[int] = []
    ends: list[int] = []
    atom_gene_rows: list[tuple[int, ...]] = []
    active: set[int] = set()
    for boundary_index, start in enumerate(boundaries[:-1]):
        active.difference_update(removals.get(start, ()))
        active.update(additions.get(start, ()))
        end = boundaries[boundary_index + 1]
        if active and end > start:
            starts.append(start)
            ends.append(end)
            atom_gene_rows.append(tuple(sorted(active)))

    row_indices: list[int] = []
    col_indices: list[int] = []
    for atom_id, gene_rows in enumerate(atom_gene_rows):
        row_indices.extend(gene_rows)
        col_indices.extend([atom_id] * len(gene_rows))
    membership = sparse.csr_matrix(
        (
            np.ones(len(row_indices), dtype=bool),
            (np.asarray(row_indices, dtype=np.int32), np.asarray(col_indices, dtype=np.int32)),
        ),
        shape=(n_genes, len(starts)),
        dtype=bool,
    )
    membership.sort_indices()
    model = ChromosomeAtomModel(
        chromosome=str(chromosome),
        starts=np.asarray(starts, dtype=np.int64),
        ends=np.asarray(ends, dtype=np.int64),
        gene_to_atom=membership,
    )
    model.validate()
    return model


def map_snps_to_atoms(snp_pos0: np.ndarray, model: ChromosomeAtomModel) -> np.ndarray:
    """Map 0-based SNP point positions to atom IDs, using -1 outside all atoms."""
    positions = np.asarray(snp_pos0, dtype=np.int64)
    if positions.ndim != 1:
        raise LDSCInternalError("SNP positions for atom mapping must be one-dimensional.")
    if model.n_atoms == 0:
        return np.full(len(positions), -1, dtype=np.int32)
    candidate = np.searchsorted(model.starts, positions, side="right") - 1
    mapped = np.full(len(positions), -1, dtype=np.int32)
    valid = candidate >= 0
    valid_indices = np.flatnonzero(valid)
    inside = positions[valid_indices] < model.ends[candidate[valid_indices]]
    mapped[valid_indices[inside]] = candidate[valid_indices[inside]].astype(np.int32)
    return mapped


def assemble_atom_selector(
    gene_to_atom: sparse.csr_matrix,
    selected_gene_rows: Sequence[int] | np.ndarray,
) -> np.ndarray:
    """Return the Boolean union of atoms covered by selected catalog rows."""
    rows = np.asarray(selected_gene_rows, dtype=np.int64)
    if rows.size == 0:
        return np.zeros(gene_to_atom.shape[1], dtype=bool)
    if np.any(rows < 0) or np.any(rows >= gene_to_atom.shape[0]):
        raise LDSCInternalError("Selected gene row is outside the gene-to-atom matrix.")
    return np.asarray(gene_to_atom[rows].max(axis=0).toarray()).reshape(-1).astype(bool)


def assemble_snp_annotation(snp_atom_indices: np.ndarray, atom_selector: np.ndarray) -> np.ndarray:
    """Scatter one selected-atom vector to the retained reference-SNP grid."""
    mapped = np.asarray(snp_atom_indices, dtype=np.int32)
    selector = np.asarray(atom_selector, dtype=bool)
    annotation = np.zeros(len(mapped), dtype=bool)
    inside = mapped >= 0
    if np.any(mapped[inside] >= len(selector)):
        raise LDSCInternalError("SNP-to-atom mapping references an undeclared atom.")
    annotation[inside] = selector[mapped[inside]]
    return annotation


def iter_snp_atom_blocks(
    snp_atom_indices: np.ndarray,
    n_atoms: int,
    batch_size: int,
) -> Iterator[tuple[int, int, np.ndarray]]:
    """Yield bounded dense Boolean SNP-by-atom blocks for offline construction."""
    if batch_size <= 0:
        raise LDSCInternalError("Atom batch size must be positive.")
    mapped = np.asarray(snp_atom_indices, dtype=np.int32)
    for start in range(0, int(n_atoms), int(batch_size)):
        end = min(start + int(batch_size), int(n_atoms))
        block = np.zeros((len(mapped), end - start), dtype=bool)
        selected = (mapped >= start) & (mapped < end)
        rows = np.flatnonzero(selected)
        block[rows, mapped[rows] - start] = True
        yield start, end, block


def compute_atom_statistics(
    snp_atom_indices: np.ndarray,
    n_atoms: int,
    baseline_annotations: np.ndarray,
    common_mask: np.ndarray,
) -> AtomStatistics:
    """Compute exact all/common atom sufficient statistics without dense H."""
    mapped = np.asarray(snp_atom_indices, dtype=np.int32)
    baseline = np.asarray(baseline_annotations, dtype=np.float64)
    common = np.asarray(common_mask, dtype=bool)
    if baseline.ndim != 2 or baseline.shape[0] != len(mapped) or common.shape != (len(mapped),):
        raise LDSCInternalError("Atom statistic inputs have incompatible SNP dimensions.")
    if np.any(mapped >= n_atoms):
        raise LDSCInternalError("SNP-to-atom mapping references an undeclared atom.")
    assigned = mapped >= 0
    atom_all = np.bincount(mapped[assigned], minlength=n_atoms).astype(np.int64)
    assigned_common = assigned & common
    atom_common = np.bincount(mapped[assigned_common], minlength=n_atoms).astype(np.int64)
    membership = sparse.csr_matrix(
        (
            np.ones(int(assigned.sum()), dtype=np.float64),
            (np.flatnonzero(assigned).astype(np.int32), mapped[assigned]),
        ),
        shape=(len(mapped), n_atoms),
        dtype=np.float64,
    )
    overlap_all = np.asarray(baseline.T @ membership, dtype=np.float64)
    overlap_common = np.asarray(baseline[common].T @ membership[common], dtype=np.float64)
    return AtomStatistics(atom_all, atom_common, overlap_all, overlap_common)


def validate_ldscore_operator(
    operator: sparse.csr_matrix,
    *,
    n_rows: int,
    n_atoms: int,
) -> None:
    """Validate the exact persisted-row-by-atom LD-score operator."""
    if not sparse.isspmatrix_csr(operator) or operator.shape != (n_rows, n_atoms):
        raise LDSCInternalError("Gene LD-score operator must be CSR with the declared dimensions.")
    if operator.dtype != np.float64:
        raise LDSCInternalError("Gene LD-score operator values must use float64 dtype.")
    if operator.indices.dtype != np.int32 or operator.indptr.dtype != np.int32:
        raise LDSCInternalError("Gene LD-score operator CSR indices must use int32 dtype.")
    if not operator.has_sorted_indices:
        raise LDSCInternalError("Gene LD-score operator column indices must be sorted within rows.")
    if not operator.has_canonical_format:
        raise LDSCInternalError("Gene LD-score operator CSR must not contain duplicate entries.")
    if not np.isfinite(operator.data).all():
        raise LDSCInternalError("Gene LD-score operator contains non-finite values.")


def assemble_selected_atom_statistics(
    statistics: AtomStatistics,
    atom_selector: np.ndarray,
    *,
    control_selector: np.ndarray | None = None,
) -> SelectedAtomStatistics:
    """Assemble exact counts and overlaps for one Boolean atom selector."""
    selector = np.asarray(atom_selector, dtype=bool)
    n_atoms = len(statistics.atom_count_all)
    if selector.shape != (n_atoms,) or np.asarray(statistics.atom_count_common).shape != (n_atoms,):
        raise LDSCInternalError("Selected atom statistics have incompatible atom dimensions.")
    overlap_all = np.asarray(statistics.baseline_atom_overlap_all, dtype=np.float64)
    overlap_common = np.asarray(statistics.baseline_atom_overlap_common, dtype=np.float64)
    if overlap_all.ndim != 2 or overlap_common.shape != overlap_all.shape or overlap_all.shape[1] != n_atoms:
        raise LDSCInternalError("Baseline-to-atom overlap arrays have incompatible dimensions.")
    control_all = control_common = None
    if control_selector is not None:
        control = np.asarray(control_selector, dtype=bool)
        if control.shape != (n_atoms,):
            raise LDSCInternalError("Control atom selector has incompatible dimensions.")
        intersection = selector & control
        control_all = int(np.asarray(statistics.atom_count_all, dtype=np.int64) @ intersection.astype(np.int64))
        control_common = int(
            np.asarray(statistics.atom_count_common, dtype=np.int64) @ intersection.astype(np.int64)
        )
    return SelectedAtomStatistics(
        count_all=int(np.asarray(statistics.atom_count_all, dtype=np.int64) @ selector.astype(np.int64)),
        count_common=int(np.asarray(statistics.atom_count_common, dtype=np.int64) @ selector.astype(np.int64)),
        baseline_overlap_all=overlap_all @ selector.astype(np.float64),
        baseline_overlap_common=overlap_common @ selector.astype(np.float64),
        control_overlap_all=control_all,
        control_overlap_common=control_common,
    )


def assemble_indexed_ld_scores(
    operator: sparse.csr_matrix,
    atom_selectors: np.ndarray,
) -> np.ndarray:
    """Multiply float64 ``Y`` by one or more Boolean selected-atom vectors."""
    selectors = np.asarray(atom_selectors, dtype=bool)
    if selectors.ndim not in {1, 2} or selectors.shape[0] != operator.shape[1]:
        raise LDSCInternalError("Atom selector dimensions do not match the LD-score operator.")
    return np.asarray(operator @ selectors.astype(np.float64), dtype=np.float64)
