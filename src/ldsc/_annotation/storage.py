"""Detached, column-selective annotation I/O.

Private values use column-major NumPy files, with binary SNP values packed into
bytes. Reads seek and decode bounded runs without a memory map or value cache.
Descriptors carry logical dimensions and encoding across chromosome workers.
The files live in a run's ``ScratchWorkspace`` (``ldsc._runtime.scratch``).
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

BINARY_BITORDER = "little"


@dataclass(frozen=True)
class ColumnStore:
    """Descriptor for a seekable SNP-by-column matrix with no resident values.

    Boolean stores pack eight SNP rows per uint8, with SNP zero in the low bit
    and zero padding in the final byte. ``n_rows`` and ``dtype`` describe logical
    values; the NPY header describes physical bytes. Continuous annotations use
    float32 and external quantile targets retain float64. Reads are detached.
    """

    path: Path
    n_rows: int
    columns: tuple[str, ...]
    dtype: np.dtype
    offset: int
    bitorder: str | None = None

    @property
    def stored_rows(self):
        return (self.n_rows + 7) // 8 if self.bitorder else self.n_rows

    @classmethod
    def create(cls, path, n_rows, columns, *, dtype=np.float32):
        """Create a private matrix with logical ``n_rows`` and ordered columns.

        Boolean dtype selects packed storage; other dtypes use dense storage.
        Fill every retained row before consumption. Newly created padding bits
        are zero and subsequent partial writes preserve them.
        """
        path = Path(path)
        dtype = np.dtype(dtype)
        columns = tuple(columns)
        if n_rows < 0 or len(set(columns)) != len(columns):
            raise ValueError("Matrix dimensions and column names must be valid and unique.")
        bitorder = BINARY_BITORDER if dtype == np.dtype(bool) else None
        stored_rows = (n_rows + 7) // 8 if bitorder else n_rows
        stored_dtype = np.dtype(np.uint8) if bitorder else dtype
        with path.open("wb") as stream:
            np.lib.format.write_array_header_2_0(stream, {
                "descr": np.lib.format.dtype_to_descr(stored_dtype),
                "fortran_order": True,
                "shape": (stored_rows, len(columns)),
            })
            offset = stream.tell()
            stream.truncate(offset + stored_rows * len(columns) * stored_dtype.itemsize)
        return cls(path, n_rows, columns, dtype, offset, bitorder)

    def write(self, start: int, values: np.ndarray, *, columns: Sequence[str] | None = None) -> None:
        """Write consecutive rows of all or selected columns without a mapping."""
        selected = self._column_indices(columns)
        values = np.asarray(values)
        if values.ndim != 2 or values.shape[1] != len(selected):
            raise ValueError("Values must have one column per selected annotation.")
        if start < 0 or start + len(values) > self.n_rows:
            raise IndexError("Annotation write exceeds the shard's row count.")
        if self.bitorder and values.dtype != bool and not np.isin(values, (0, 1)).all():
            raise ValueError("Packed annotation values must be exactly zero or one.")
        if not len(values):
            return
        # A byte-aligned write covering whole bytes, or ending at the last row
        # (whose padding bits are always zero), replaces its bytes without reading them.
        whole_bytes = start % 8 == 0 and (len(values) % 8 == 0 or start + len(values) == self.n_rows)
        with self.path.open("r+b") as stream:
            for j, column in enumerate(selected):
                if self.bitorder and whole_bytes:
                    stream.seek(self.offset + column * self.stored_rows + start // 8)
                    stream.write(np.packbits(values[:, j].astype(bool, copy=False), bitorder=self.bitorder).tobytes())
                elif self.bitorder:
                    first, stop = start // 8, (start + len(values) + 7) // 8
                    offset = self.offset + column * self.stored_rows + first
                    stream.seek(offset)
                    # Preserve neighboring SNPs when writes share a boundary byte.
                    packed = np.frombuffer(stream.read(stop - first), dtype=np.uint8)
                    if len(packed) != stop - first:
                        raise OSError(f"Truncated annotation shard: {self.path}")
                    bits = np.unpackbits(packed, bitorder=self.bitorder)
                    bits[start % 8:start % 8 + len(values)] = values[:, j]
                    stream.seek(offset)
                    stream.write(np.packbits(bits, bitorder=self.bitorder).tobytes())
                else:
                    stream.seek(self.offset + (column * self.n_rows + start) * self.dtype.itemsize)
                    stream.write(np.asarray(values[:, j], dtype=self.dtype).tobytes())

    def _column_indices(self, columns):
        if columns is None:
            return list(range(len(self.columns)))
        lookup = {name: index for index, name in enumerate(self.columns)}
        return [lookup[name] for name in columns]

    def read(self, *, rows=None, columns=None, max_read_rows=65536) -> np.ndarray:
        """Read selected rows/columns, preserving order and repeated row indices.

        ``max_read_rows`` bounds each row run, including coalesced gaps. Packed
        runs read covering bytes, decoding at most fourteen extra boundary bits.
        The returned array itself has the explicitly requested shape. Callers
        must select a query batch and an LD block to bound that shape.
        """
        if max_read_rows < 1:
            raise ValueError("max_read_rows must be positive.")
        selected = self._column_indices(columns)
        if rows is None:
            indices = np.arange(self.n_rows, dtype=np.int64)
        elif isinstance(rows, slice):
            indices = np.arange(*rows.indices(self.n_rows), dtype=np.int64)
        else:
            indices = np.asarray(rows)
            if indices.dtype == bool:
                if indices.shape != (self.n_rows,):
                    raise IndexError("Boolean annotation selector has the wrong length.")
                indices = np.flatnonzero(indices)
            indices = indices.astype(np.int64, copy=False)
        if indices.ndim != 1 or np.any(indices < 0) or np.any(indices >= self.n_rows):
            raise IndexError("Annotation row index is outside the shard.")
        result = np.empty((len(indices), len(selected)), dtype=self.dtype)
        if not len(indices) or not selected:
            return result
        order = np.argsort(indices, kind="stable")
        sorted_rows = indices[order]
        # The row runs are the same for every column, so plan them once.
        runs, cursor = [], 0
        while cursor < len(order):
            start = int(sorted_rows[cursor])
            stop = int(np.searchsorted(sorted_rows, start + max_read_rows, side="left"))
            # Coalesce only nearby requested rows, avoiding large reads
            # for sparse selections while retaining sequential I/O.
            gaps = np.flatnonzero(np.diff(sorted_rows[cursor:stop]) > 64)
            if len(gaps):
                stop = cursor + int(gaps[0]) + 1
            runs.append((cursor, stop, start, int(sorted_rows[stop - 1]) - start + 1))
            cursor = stop
        with self.path.open("rb") as stream:
            for j, column in enumerate(selected):
                for cursor, stop, start, count in runs:
                    if self.bitorder:
                        first, last = start // 8, (start + count + 7) // 8
                        stream.seek(self.offset + column * self.stored_rows + first)
                        packed = np.frombuffer(stream.read(last - first), dtype=np.uint8)
                        if len(packed) != last - first:
                            raise OSError(f"Truncated annotation shard: {self.path}")
                        values = np.unpackbits(packed, bitorder=self.bitorder)[start % 8:start % 8 + count]
                    else:
                        stream.seek(self.offset + (column * self.n_rows + start) * self.dtype.itemsize)
                        values = np.frombuffer(stream.read(count * self.dtype.itemsize), dtype=self.dtype)
                    if len(values) != count:
                        raise OSError(f"Truncated annotation shard: {self.path}")
                    result[order[cursor:stop], j] = values[sorted_rows[cursor:stop] - start]
        return result


@dataclass(frozen=True)
class AnnotationShard:
    """Chromosome descriptor with separate row metadata and column value files."""

    metadata_path: Path
    stores: tuple[ColumnStore, ...]
    n_rows: int

    @property
    def columns(self):
        return tuple(name for store in self.stores for name in store.columns)

    def metadata(self, columns=None):
        """Load this chromosome's row metadata, or only ``columns``; the caller owns the frame."""
        return pd.read_parquet(self.metadata_path, columns=columns)

    def read(self, *, columns, rows=None, max_read_rows=65536):
        """Read separate artifacts in the required explicit ``columns`` order.

        Binary-only selections return Boolean values; mixed selections use the
        stores' common NumPy dtype. The public bundle converts reads to float32
        and supplies its logical default column order.
        """
        columns = tuple(columns)
        lookup = {column: store for store in self.stores for column in store.columns}
        if not columns:
            return self.stores[0].read(rows=rows, columns=(), max_read_rows=max_read_rows)
        groups = {}
        for index, name in enumerate(columns):
            groups.setdefault(lookup[name], []).append((index, name))
        if len(groups) == 1:
            return lookup[columns[0]].read(rows=rows, columns=columns, max_read_rows=max_read_rows)
        result = None
        for store, selection in groups.items():
            values = store.read(rows=rows, columns=[name for _, name in selection], max_read_rows=max_read_rows)
            if result is None:
                result = np.empty((len(values), len(columns)), dtype=np.result_type(*[item.dtype for item in groups]))
            result[:, [index for index, _ in selection]] = values
        return result
