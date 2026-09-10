"""Output-contained scratch and detached, column-selective annotation I/O.

Private value files use the NumPy column-major format. Reads explicitly seek
and copy bounded runs; no memory map or chromosome cache keeps unused columns
resident. Metadata and value descriptors are safe to pass to chromosome workers.
"""

from dataclasses import dataclass
from pathlib import Path
import shutil
import tempfile
from typing import Sequence

import numpy as np
import pandas as pd


class AnnotationWorkspace:
    """Own a unique scratch directory below an explicit workflow destination.

    Closing is idempotent and removes only this owner's private directory.
    Public artifacts and source files are never owned by this object. Borrowers
    must not close their caller's workspace.
    """

    def __init__(self, output_dir: str | Path):
        destination = Path(output_dir)
        destination.mkdir(parents=True, exist_ok=True)
        self.path = Path(tempfile.mkdtemp(prefix=".ldsc-annotation-", dir=destination))
        self.closed = False

    def require_open(self) -> None:
        if self.closed:
            raise ValueError("Annotation workspace is closed.")

    def close(self) -> None:
        if not self.closed:
            shutil.rmtree(self.path)
            self.closed = True

    def __enter__(self):
        self.require_open()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


@dataclass(frozen=True)
class ColumnStore:
    """Descriptor for a seekable SNP-by-column matrix with no resident values.

    Ordinary annotations use float32; exact external quantile targets use
    float64. Requested arrays are detached and belong to the caller. An active
    chromosome consumer controls their lifetime.
    """

    path: Path
    n_rows: int
    columns: tuple[str, ...]
    dtype: np.dtype
    offset: int

    @classmethod
    def create(cls, path, n_rows, columns, *, dtype=np.float32):
        """Create an empty matrix; fill every retained row before consumption."""
        path = Path(path)
        dtype = np.dtype(dtype)
        columns = tuple(columns)
        if n_rows < 0 or len(set(columns)) != len(columns):
            raise ValueError("Matrix dimensions and column names must be valid and unique.")
        with path.open("wb") as stream:
            np.lib.format.write_array_header_2_0(stream, {
                "descr": np.lib.format.dtype_to_descr(dtype),
                "fortran_order": True,
                "shape": (n_rows, len(columns)),
            })
            offset = stream.tell()
            stream.truncate(offset + n_rows * len(columns) * dtype.itemsize)
        return cls(path, n_rows, columns, dtype, offset)

    def write(self, start: int, values: np.ndarray, *, columns: Sequence[str] | None = None) -> None:
        """Write consecutive rows of all or selected columns without a mapping."""
        selected = self._column_indices(columns)
        values = np.asarray(values)
        if values.ndim != 2 or values.shape[1] != len(selected):
            raise ValueError("Values must have one column per selected annotation.")
        if start < 0 or start + len(values) > self.n_rows:
            raise IndexError("Annotation write exceeds the shard's row count.")
        with self.path.open("r+b") as stream:
            for j, column in enumerate(selected):
                stream.seek(self.offset + (column * self.n_rows + start) * self.dtype.itemsize)
                stream.write(np.asarray(values[:, j], dtype=self.dtype).tobytes())

    def _column_indices(self, columns):
        if columns is None:
            return list(range(len(self.columns)))
        lookup = {name: index for index, name in enumerate(self.columns)}
        return [lookup[name] for name in columns]

    def read(self, *, rows=None, columns=None, max_read_rows=65536) -> np.ndarray:
        """Read selected rows/columns, preserving order and repeated row indices.

        ``max_read_rows`` bounds each physical read, including coalesced gaps.
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
        with self.path.open("rb") as stream:
            for j, column in enumerate(selected):
                cursor = 0
                while cursor < len(order):
                    start = int(sorted_rows[cursor])
                    stop = int(np.searchsorted(sorted_rows, start + max_read_rows, side="left"))
                    # Coalesce only nearby requested rows, avoiding large reads
                    # for sparse selections while retaining sequential I/O.
                    gaps = np.flatnonzero(np.diff(sorted_rows[cursor:stop]) > 64)
                    if len(gaps):
                        stop = cursor + int(gaps[0]) + 1
                    count = int(sorted_rows[stop - 1]) - start + 1
                    stream.seek(self.offset + (column * self.n_rows + start) * self.dtype.itemsize)
                    values = np.frombuffer(stream.read(count * self.dtype.itemsize), dtype=self.dtype)
                    if len(values) != count:
                        raise OSError(f"Truncated annotation shard: {self.path}")
                    result[order[cursor:stop], j] = values[sorted_rows[cursor:stop] - start]
                    cursor = stop
        return result


@dataclass
class FrameSpool:
    """Append bounded private frames and replay them without retaining records.

    The pickles are internal, locally generated scratch, never user inputs.
    Only a path, row count, and part count are retained in memory.
    """

    path: Path
    n_parts: int = 0
    n_rows: int = 0

    def append(self, frame):
        if frame.empty:
            return
        self.path.mkdir(parents=True, exist_ok=True)
        frame.to_pickle(self.path / f"{self.n_parts}.pkl")
        self.n_parts += 1
        self.n_rows += len(frame)

    def frames(self):
        for index in range(self.n_parts):
            yield pd.read_pickle(self.path / f"{index}.pkl")

    def write_tsv(self, path, *, columns=(), na_rep=""):
        """Publish all records, retaining only one append chunk at a time."""
        import gzip

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "wt") as stream:
            header = True
            for frame in self.frames():
                frame.to_csv(stream, sep="\t", index=False, header=header, na_rep=na_rep)
                header = False
            if header:
                pd.DataFrame(columns=columns).to_csv(stream, sep="\t", index=False)


@dataclass(frozen=True)
class AnnotationShard:
    """Chromosome descriptor with separate row metadata and column value files."""

    metadata_path: Path
    stores: tuple[ColumnStore, ...]
    n_rows: int

    @property
    def columns(self):
        return tuple(name for store in self.stores for name in store.columns)

    def metadata(self):
        """Load this chromosome's row metadata; the caller owns the frame."""
        return pd.read_parquet(self.metadata_path)

    def read(self, *, rows=None, columns=None, max_read_rows=65536):
        """Read selected annotations in requested order from separate artifacts."""
        columns = self.columns if columns is None else tuple(columns)
        lookup = {column: store for store in self.stores for column in store.columns}
        if not columns:
            return self.stores[0].read(rows=rows, columns=(), max_read_rows=max_read_rows)
        first = lookup[columns[0]].read(rows=rows, columns=[columns[0]], max_read_rows=max_read_rows)
        result = np.empty((len(first), len(columns)), dtype=first.dtype)
        result[:, 0] = first[:, 0]
        for j, column in enumerate(columns[1:], 1):
            result[:, j] = lookup[column].read(rows=rows, columns=[column], max_read_rows=max_read_rows)[:, 0]
        return result
