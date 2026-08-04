from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from ldsc.errors import LDSCInternalError

from ldsc._kernel.gene_ldscore_index import (
    assemble_atom_selector,
    assemble_indexed_ld_scores,
    assemble_selected_atom_statistics,
    assemble_snp_annotation,
    build_disjoint_atoms,
    compute_atom_statistics,
    iter_snp_atom_blocks,
    map_snps_to_atoms,
    validate_ldscore_operator,
)


def test_disjoint_atoms_reproduce_boolean_interval_unions():
    intervals = np.array(
        [
            [0, 10],
            [5, 15],
            [6, 8],
            [20, 30],
            [15, 20],
        ],
        dtype=np.int64,
    )
    included = np.array([True, True, True, False, True])

    model = build_disjoint_atoms("1", intervals, included=included)

    np.testing.assert_array_equal(model.starts, [0, 5, 6, 8, 10, 15])
    np.testing.assert_array_equal(model.ends, [5, 6, 8, 10, 15, 20])
    np.testing.assert_array_equal(
        model.gene_to_atom.toarray(),
        np.array(
            [
                [1, 1, 1, 1, 0, 0],
                [0, 1, 1, 1, 1, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1],
            ],
            dtype=bool,
        ),
    )

    snp_pos0 = np.array([0, 4, 5, 6, 7, 8, 9, 10, 14, 15, 19, 20, 25])
    snp_atoms = map_snps_to_atoms(snp_pos0, model)
    selector = assemble_atom_selector(model.gene_to_atom, [0, 1, 0])
    actual = assemble_snp_annotation(snp_atoms, selector)
    expected = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0], dtype=bool)
    np.testing.assert_array_equal(actual, expected)
    assert actual.max(initial=False) <= 1


def test_atom_statistics_match_known_dense_counts_and_overlaps():
    model = build_disjoint_atoms(
        "1",
        np.array([[0, 5], [5, 10]], dtype=np.int64),
        included=np.array([True, True]),
    )
    snp_atoms = map_snps_to_atoms(np.array([0, 4, 5, 7, 10]), model)
    baseline = np.array(
        [
            [1.0, 0.5],
            [0.0, 1.5],
            [2.0, 0.0],
            [3.0, 2.0],
            [4.0, 4.0],
        ],
        dtype=np.float64,
    )
    common = np.array([True, False, True, False, True])

    stats = compute_atom_statistics(snp_atoms, 2, baseline, common)

    np.testing.assert_array_equal(stats.atom_count_all, [2, 2])
    np.testing.assert_array_equal(stats.atom_count_common, [1, 1])
    np.testing.assert_allclose(stats.baseline_atom_overlap_all, [[1.0, 5.0], [2.0, 2.0]])
    np.testing.assert_allclose(stats.baseline_atom_overlap_common, [[1.0, 2.0], [0.5, 0.0]])

    selected = assemble_selected_atom_statistics(
        stats,
        np.array([True, True]),
        control_selector=np.array([False, True]),
    )
    assert selected.count_all == 4
    assert selected.count_common == 2
    np.testing.assert_allclose(selected.baseline_overlap_all, [6.0, 4.0])
    np.testing.assert_allclose(selected.baseline_overlap_common, [3.0, 0.5])
    assert selected.control_overlap_all == 2
    assert selected.control_overlap_common == 1


def test_empty_atom_model_is_structurally_valid():
    model = build_disjoint_atoms(
        "22",
        np.empty((0, 2), dtype=np.int64),
        included=np.empty(0, dtype=bool),
    )

    assert model.gene_to_atom.shape == (0, 0)
    assert map_snps_to_atoms(np.array([1, 2]), model).tolist() == [-1, -1]


def test_padding_clips_at_chromosome_start_and_atom_batches_reconstruct_h():
    model = build_disjoint_atoms(
        "1",
        np.array([[2, 4], [10, 12]], dtype=np.int64),
        padding_bp=5,
    )
    np.testing.assert_array_equal(model.starts, [0, 5, 9])
    np.testing.assert_array_equal(model.ends, [5, 9, 17])
    snp_atoms = map_snps_to_atoms(np.array([0, 4, 5, 6, 7, 16, 17]), model)

    blocks = [block for _start, _end, block in iter_snp_atom_blocks(snp_atoms, model.n_atoms, 2)]
    h = np.concatenate(blocks, axis=1)

    np.testing.assert_array_equal(h.sum(axis=1), [1, 1, 1, 1, 1, 1, 0])
    assert all(block.shape[1] <= 2 for block in blocks)


def test_genes_with_no_retained_atom_snps_have_zero_statistics():
    model = build_disjoint_atoms("22", np.array([[100, 200]], dtype=np.int64))
    snp_atoms = map_snps_to_atoms(np.array([1, 2]), model)
    stats = compute_atom_statistics(
        snp_atoms,
        model.n_atoms,
        np.ones((2, 1), dtype=np.float64),
        np.array([True, False]),
    )

    np.testing.assert_array_equal(stats.atom_count_all, [0])
    np.testing.assert_array_equal(stats.atom_count_common, [0])
    np.testing.assert_array_equal(stats.baseline_atom_overlap_all, [[0.0]])


def test_index_operator_preserves_exact_linear_identity_and_negative_values():
    # SNPs 0/1 share atom 0; SNP 2 belongs to atom 1; SNP 3 is outside all genes.
    h = np.array(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]
    )
    r = np.array(
        [
            [1.0, -0.25, 0.5, 0.1],
            [-0.25, 1.0, -0.8, 0.2],
            [0.5, -0.8, 1.0, -0.3],
            [0.1, 0.2, -0.3, 1.0],
        ],
        dtype=np.float64,
    )
    persisted_rows = np.array([True, False, True, False])
    y = sparse.csr_matrix((r @ h)[persisted_rows], dtype=np.float64)
    y.sort_indices()
    z = np.array([True, True])

    validate_ldscore_operator(y, n_rows=2, n_atoms=2)
    actual = assemble_indexed_ld_scores(y, z)
    direct = (r @ (h @ z.astype(np.float64)))[persisted_rows]

    np.testing.assert_allclose(actual, direct, rtol=0, atol=0)
    assert y.data.min() < 0


def test_index_operator_rejects_duplicate_csr_entries():
    operator = sparse.csr_matrix(
        (
            np.array([1.0, -0.25], dtype=np.float64),
            np.array([0, 0], dtype=np.int32),
            np.array([0, 2], dtype=np.int32),
        ),
        shape=(1, 1),
    )
    assert operator.has_sorted_indices
    assert not operator.has_canonical_format

    with pytest.raises(LDSCInternalError, match="duplicate"):
        validate_ldscore_operator(operator, n_rows=1, n_atoms=1)
