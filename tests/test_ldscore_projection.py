"""Independent dense oracles for batched output-row LD projection."""

import numpy as np
import pytest


@pytest.mark.parametrize('batch_size', [1, 2, 3, 1000])
@pytest.mark.parametrize('sparse_pairs', [False, True])
def test_projection_preserves_nonoutput_contributors_and_weight_universe(batch_size, sparse_pairs):
    from ldsc._kernel.ldscore_projection import ArrayAnnotations, ProjectionAccumulator

    annotations = np.array([[1, 1, 0, -2], [1, 0, 3, 4], [1, 2, 1, 0], [1, 0, 0, 5]], dtype=np.float32)
    r2 = np.array([[1, .25, .5, 0], [.25, 1, .125, -.05], [.5, .125, 1, .75], [0, -.05, .75, 1]], dtype=np.float64)
    output = np.array([2, 0])
    weights = np.array([1, 0, 1, 0], dtype=np.float32)
    source = ArrayAnnotations(annotations, ('base', 'a', 'b', 'c'))
    accumulator = ProjectionAccumulator(source, n_baseline=1, output_rows=output, weight_mask=weights, query_batch_size=batch_size, read_budget_bytes=32)
    if sparse_pairs:
        accumulator.add_diagonal()
        i, j = np.triu_indices(4, 1)
        for start in range(0, len(i), 2):
            accumulator.add_pairs(i[start:start+2], j[start:start+2], r2[i[start:start+2],j[start:start+2]], np.zeros(4,dtype=int))
    else:
        accumulator.add_dense(np.arange(4), np.arange(4), r2)
    assert accumulator.values.shape == (2, 5)
    assert accumulator.values.dtype == np.float64
    np.testing.assert_allclose(accumulator.values[:,:4], (r2 @ annotations.astype(np.float64))[output], rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(accumulator.values[:,-1], (r2 @ weights)[output], rtol=1e-14, atol=1e-14)


def test_pair_window_filter_and_no_output_rows(tmp_path):
    from ldsc._kernel.ldscore_projection import ArrayAnnotations, ProjectionAccumulator

    source = ArrayAnnotations(np.array([[1,2],[1,4],[1,8]],dtype=np.float32), ('base','q'))
    accumulator=ProjectionAccumulator(source, n_baseline=1, output_rows=[2], query_batch_size=1)
    accumulator.add_diagonal()
    accumulator.add_pairs(np.array([0,1]),np.array([2,2]),np.array([.5,.25]),np.array([0,0,1]))
    np.testing.assert_array_equal(accumulator.values, [[1.25,9]])
    empty=ProjectionAccumulator(source, n_baseline=1, output_rows=[], query_batch_size=1)
    empty.add_diagonal()
    empty.add_dense(np.arange(3),np.arange(3),np.eye(3))
    assert empty.values.shape == (0,2)


def test_accessor_reads_only_query_batches_and_bounded_contributor_tiles():
    from ldsc._kernel.ldscore_projection import ArrayAnnotations, ProjectionAccumulator

    reads=[]
    class RecordingAnnotations(ArrayAnnotations):
        def read(self, *, rows=None, columns=None):
            result=super().read(rows=rows,columns=columns)
            reads.append((tuple(columns),result.shape))
            return result
    source=RecordingAnnotations(np.ones((20,8),dtype=np.float32),tuple(['base',*[f'q{i}' for i in range(7)]]))
    accumulator=ProjectionAccumulator(source,n_baseline=1,output_rows=[0,2],query_batch_size=3,read_budget_bytes=48)
    accumulator.add_dense(np.arange(20),np.arange(20),np.ones((20,20)))
    np.testing.assert_array_equal(accumulator.values,np.full((2,8),20.))
    assert all(len(columns)<=3 and np.prod(shape)*8<=48 for columns,shape in reads)


@pytest.mark.parametrize('batch_size', [1, 2, 1000])
def test_plink_projection_reuses_genotype_traversal(batch_size):
    from ldsc._kernel.plink_bed import __GenotypeArrayInMemory__
    from ldsc._kernel.ldscore import ld_scores_from_genotypes
    from ldsc._kernel.ldscore_projection import ArrayAnnotations

    genotypes = np.array([[0,1,2,0], [1,2,0,1], [2,0,1,0], [0,0,2,2], [1,1,0,2]], dtype=float)
    genotypes = (genotypes-genotypes.mean(axis=0))/genotypes.std(axis=0)
    annotations = np.array([[1,0,2,3], [1,4,0,2], [1,1,-1,0], [1,0,3,1]], dtype=np.float32)
    reader = object.__new__(__GenotypeArrayInMemory__)
    reader.m, reader.n = 4, 5
    seen = []
    def next_snps(count):
        start = sum(seen)
        seen.append(count)
        return genotypes[:, start:start+count]
    reader.nextSNPs = next_snps
    scores = ld_scores_from_genotypes(reader, np.array([0,0,0,1]), 1,
        annot=ArrayAnnotations(annotations, ('base','q1','q2','q3')),
        output_rows=[2,0], n_baseline=1, query_batch_size=batch_size,
        weight_mask=np.array([1,0,1,0], dtype=np.float32))
    r2 = (genotypes.T @ genotypes / 5)**2
    r2 -= (1-r2)/3
    r2[0,3] = r2[3,0] = 0
    expected = r2 @ np.column_stack([annotations, [1,0,1,0]])
    np.testing.assert_allclose(scores, expected[[2,0]], rtol=1e-14, atol=1e-14)
    assert sum(seen) == 4
    assert scores.dtype == np.float64


def _window_coordinates(name, m):
    """Coordinates and maximum distance defining one exact LD window."""
    if name == 'gap':
        # A 100-unit gap leaves the first SNP after it with no left partners.
        return np.r_[np.arange(10), 110 + np.arange(m - 10)].astype(float), 3.0
    if name == 'whole':
        return np.arange(m, dtype=float), float(m)
    return np.arange(m, dtype=float), float(name.removeprefix('snps'))


@pytest.mark.parametrize('chunk', [1, 2, 3, 4, 7, 50])
@pytest.mark.parametrize('window', ['snps1', 'snps5', 'gap', 'whole'])
def test_plink_ld_scores_use_the_exact_window_for_every_chunk_size(window, chunk):
    """LD scores sum adjusted r2 over the exact window, whatever the batch size.

    The oracle is a dense adjusted-r2 matrix zeroed outside the window rule
    ``min(i, j) >= block_left[max(i, j)]`` (the parquet reader's rule). m=23 is
    prime, so batches 2-7 end with a partial chunk; batch 50 exceeds m and
    exercises the single initial block. ``gap`` reaches the empty-left-block path.
    """
    from ldsc._kernel.ldscore import get_block_lefts, ld_scores_from_genotypes
    from ldsc._kernel.ldscore_projection import ArrayAnnotations
    from ldsc._kernel.plink_bed import __GenotypeArrayInMemory__

    rng = np.random.default_rng(20261027)
    n, m = 40, 23
    genotypes = rng.integers(0, 3, size=(n, m)).astype(float)
    genotypes = (genotypes - genotypes.mean(axis=0)) / genotypes.std(axis=0)
    annotations = np.column_stack([np.ones(m), rng.normal(size=(m, 2))]).astype(np.float32)
    weights = rng.integers(0, 2, size=m).astype(np.float32)
    output_rows = np.array([20, 3, 0, 11, 22])
    coords, max_dist = _window_coordinates(window, m)
    block_left = get_block_lefts(coords, max_dist)

    reader = object.__new__(__GenotypeArrayInMemory__)
    reader.m, reader.n = m, n
    cursor = [0]
    def next_snps(count):
        start = cursor[0]
        cursor[0] += count
        return genotypes[:, start:start + count]
    reader.nextSNPs = next_snps
    scores = ld_scores_from_genotypes(reader, block_left, chunk,
        annot=ArrayAnnotations(annotations, ('base', 'q1', 'q2')),
        output_rows=output_rows, n_baseline=1, query_batch_size=2, weight_mask=weights)

    r = genotypes.T @ genotypes / n
    r2 = r**2 - (1 - r**2) / (n - 2)
    i, j = np.meshgrid(np.arange(m), np.arange(m), indexing='ij')
    inside = np.minimum(i, j) >= block_left.astype(int)[np.maximum(i, j)]
    expected = np.where(inside, r2, 0.0) @ np.column_stack([annotations, weights]).astype(np.float64)
    assert cursor[0] == m
    assert scores.shape == (len(output_rows), 4)
    # Only summation order differs from the oracle, so float64 roundoff bounds the error.
    np.testing.assert_allclose(scores, expected[output_rows], rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize('batch_size', [1, 2, 1000])
def test_r2_projection_reuses_pair_traversal(batch_size):
    from ldsc._kernel.ldscore import ld_score_streaming_from_r2_reader
    from ldsc._kernel.ldscore_projection import ArrayAnnotations

    class PairReader:
        traversals = 0
        def iter_all_pairs(self):
            self.traversals += 1
            yield np.array([0,0,1]), np.array([1,2,2]), np.array([.5,.25,.75], dtype=np.float32)
    reader = PairReader()
    annotations = np.array([[1,2,0,1], [1,0,4,2], [1,8,2,-1]], dtype=np.float32)
    r2 = np.array([[1,.5,.25],[.5,1,.75],[.25,.75,1]])
    actual = ld_score_streaming_from_r2_reader(np.zeros(3,dtype=int),
        ArrayAnnotations(annotations, ('base','q1','q2','q3')), reader,
        chunk_pairs=2, output_rows=[2,0], n_baseline=1, query_batch_size=batch_size,
        weight_mask=np.array([1,0,1],dtype=np.float32))
    np.testing.assert_array_equal(actual, (r2 @ np.column_stack([annotations,[1,0,1]]))[[2,0]].astype(np.float32))
    assert reader.traversals == 1
