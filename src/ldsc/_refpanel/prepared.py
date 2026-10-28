"""The prepared-chromosome contract between reference panels and LD-score computation.

``RefPanel.prepare_chromosome`` returns a ``PreparedChromosome``: retained
reference-SNP metadata, an annotation accessor that maps each retained SNP
to its stored annotation row, LD-window bounds and one owned reader.
``ChromosomeAnnotations`` is the chromosome's annotation input to that
preparation. LD-score computation reads the accessor only through its
``shape``, ``columns`` and ``read()``.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .._core.snp_identity import empty_identity_drop_frame


@dataclass
class ChromosomeAnnotations:
    """Chromosome metadata and a selected-read annotation source.

    Explicitly supplied DataFrames are numerical inputs owned by the caller;
    workflow preparation supplies a shard descriptor instead.
    """
    metadata: pd.DataFrame
    annotations: object
    baseline_columns: list[str]
    query_columns: list[str]

    def __post_init__(self):
        if isinstance(self.annotations, pd.DataFrame):
            self.annotations = ArrayAnnotations(
                self.annotations.to_numpy(dtype=np.float32, copy=False),
                tuple(self.baseline_columns + self.query_columns))


@dataclass
class PreparedChromosome:
    """Aligned reference state whose reader belongs to one chromosome calculation.

    ``annotations`` maps selected reads to source rows in ``metadata`` order.
    ``block_left`` describes all of those reference SNPs. Use as a context
    manager; closing releases the reader and chromosome working state.
    """

    backend: str
    reader: object
    metadata: pd.DataFrame
    annotations: object
    block_left: np.ndarray
    baseline_columns: list[str]
    query_columns: list[str]
    reference_rows_before_genotype_qc: int = 0
    genotype_qc_removed: int = 0
    maf_removed: int = 0
    selected_individual_count: int = 0
    cm_source: str = "bim_cm"
    identity_drops: pd.DataFrame = field(default_factory=empty_identity_drop_frame)

    def close(self):
        """Release chromosome state once, including when reader closure fails."""
        reader = self.reader
        self.reader = None
        self.annotations = None
        self.metadata = None
        self.block_left = None
        self.identity_drops = None
        if reader is not None:
            reader.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


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
