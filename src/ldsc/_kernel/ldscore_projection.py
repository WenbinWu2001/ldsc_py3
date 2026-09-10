"""Batched LD projection with full contributors and output-row-only scores.

The numerical traversal supplies each correlation block or sparse pair chunk
once. This module reuses it for shared columns/weights and every focal batch.
Contributors may be absent from the output SNP set. Accumulation is float64;
the calling workflow controls its established artifact precision conversion.
"""

from dataclasses import dataclass

import numpy as np


@dataclass
class ArrayAnnotations:
    """Selected reads over a caller-owned, explicitly supplied numerical array."""

    values: np.ndarray
    columns: tuple[str, ...]

    @property
    def shape(self):
        return self.values.shape

    def read(self, *, rows=None, columns=None):
        indices = np.arange(self.shape[0]) if rows is None else np.arange(self.shape[0])[rows] if isinstance(rows, slice) else np.asarray(rows, dtype=int)
        selected = range(self.shape[1]) if columns is None else [self.columns.index(name) for name in columns]
        return self.values[np.ix_(indices, selected)]


@dataclass
class MappedAnnotations:
    """Map reference row requests to a resolved source's row/column reader.

    ``source`` supplies ``read(rows=..., columns=...)`` and contains only a
    chromosome descriptor. ``row_indices`` maps each retained reference SNP
    to its annotation row. No genome-wide matrix or source cache is created.
    """

    source: object
    row_indices: np.ndarray
    columns: tuple[str, ...]

    @property
    def shape(self):
        return len(self.row_indices), len(self.columns)

    def read(self, *, rows=None, columns=None):
        selected = self.row_indices if rows is None else self.row_indices[rows]
        return self.source.read(rows=selected, columns=self.columns if columns is None else columns)


class ProjectionAccumulator:
    """Own float64 output-SNP-by-annotation scores and bounded read workspace."""

    def __init__(self, annotations, *, n_baseline, output_rows=None, weight_mask=None, query_batch_size=1000, read_budget_bytes=16*1024*1024):
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
        """Seed R²=1 contributions on output rows, including selected weights."""
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
        """Project one already-computed dense block without narrowing contributors."""
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
