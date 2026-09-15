"""Prepared annotations and complete diagnostics for small file-free calculations."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .annotation_semantics import require_unique_annotation_names
from .chromosome_inference import normalize_chromosome
from .config import get_global_config
from .errors import LDSCInputError
from ._annotation_storage import BINARY_BITORDER
from ._kernel.snp_identity import clean_identity_artifact_table, identity_base_mode, identity_mode_family


@dataclass(frozen=True)
class MemoryDiagnostics:
    """Complete diagnostic rows owned by an in-memory result."""

    frame: pd.DataFrame

    def frames(self, *, chunk_rows=65536):
        for start in range(0, len(self.frame), chunk_rows):
            yield self.frame.iloc[start:start + chunk_rows]

    def write_to(self, path):
        """Serialize only when the caller explicitly chooses a writing workflow."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.frame.to_csv(path, sep="\t", index=False, na_rep="")


@dataclass(frozen=True)
class MemoryAnnotationShard:
    """Own chromosome metadata, packed binary values, and dense float32 values."""

    frame: pd.DataFrame
    values: np.ndarray
    continuous_columns: tuple[str, ...]
    binary_values: np.ndarray
    binary_columns: tuple[str, ...]
    bitorder: str = BINARY_BITORDER

    @property
    def columns(self):
        return self.binary_columns + self.continuous_columns

    @property
    def n_rows(self):
        return len(self.frame)

    def metadata(self):
        return self.frame.copy()

    def read(self, *, columns, rows=None, max_read_rows=65536):
        """Decode selected SNPs in the requested column order as float32."""
        indices = np.arange(self.n_rows) if rows is None else np.arange(self.n_rows)[rows]
        binary = {name: index for index, name in enumerate(self.binary_columns)}
        continuous = {name: index for index, name in enumerate(self.continuous_columns)}
        result = np.empty((len(indices), len(columns)), dtype=np.float32)
        shifts = indices % 8 if self.bitorder == "little" else 7 - indices % 8
        for j, name in enumerate(columns):
            if name in binary:
                result[:, j] = (self.binary_values[indices // 8, binary[name]] >> shifts) & 1
            else:
                result[:, j] = self.values[indices, continuous[name]]
        return result


def bundle_from_frames(cls, metadata, baseline_annotations, query_annotations, config_snapshot):
    """Apply the global annotation identity policy without creating a workspace."""
    config = config_snapshot or get_global_config()
    required = {"CHR", "POS"} | ({"SNP"} if identity_mode_family(config.snp_identifier) == "rsid" else set())
    if not required.issubset(metadata):
        raise LDSCInputError(f"Prepared SNP metadata is missing required columns: {sorted(required.difference(metadata))}.")
    if ("A1" in metadata) != ("A2" in metadata):
        raise LDSCInputError("Prepared SNP metadata must provide both allele columns or neither.")
    baseline = list(baseline_annotations.columns)
    queries = [] if query_annotations is None else list(query_annotations.columns)
    require_unique_annotation_names(baseline, queries)
    names = [*baseline, *queries]
    if not baseline or len(baseline_annotations) != len(metadata) or (query_annotations is not None and len(query_annotations) != len(metadata)):
        raise LDSCInputError("Prepared annotations require baseline columns and the same row count as metadata.")
    if set(names).intersection(metadata.columns):
        raise LDSCInputError("Prepared annotation names must not overlap SNP metadata column names.")
    values = pd.concat([baseline_annotations.reset_index(drop=True),
                        query_annotations.reset_index(drop=True) if query_annotations is not None else pd.DataFrame()], axis=1).apply(pd.to_numeric)
    binary = values.isin((0, 1)).all(axis=0)
    binary_columns, continuous_columns = tuple(values.columns[binary]), tuple(values.columns[~binary])
    continuous = values.loc[:, list(continuous_columns)].to_numpy(dtype=np.float32)
    if not np.isfinite(continuous).all():
        raise LDSCInputError("Prepared annotation values must be finite numeric values.")
    frame = metadata.reset_index(drop=True).copy()
    frame["CHR"] = frame["CHR"].map(normalize_chromosome)
    frame["POS"] = pd.to_numeric(frame["POS"], errors="raise").astype(np.int64)
    frame["_memory_row"] = np.arange(len(frame))
    mode = config.snp_identifier if {"A1", "A2"}.issubset(frame) else identity_base_mode(config.snp_identifier)
    clean = clean_identity_artifact_table(frame, mode, context="prepared annotations", stage="annotation_identity_cleanup")
    shards = {}
    for chrom, rows in clean.cleaned.groupby("CHR", sort=False):
        rows = rows.sort_values("POS", kind="stable")
        indices = rows["_memory_row"].to_numpy(dtype=np.int64)
        packed = np.empty(((len(rows) + 7) // 8, len(binary_columns)), dtype=np.uint8, order="F")
        for j, name in enumerate(binary_columns):
            packed[:, j] = np.packbits(values[name].to_numpy()[indices].astype(bool, copy=False), bitorder=BINARY_BITORDER)
        shards[chrom] = MemoryAnnotationShard(rows.drop(columns="_memory_row").reset_index(drop=True),
                                             continuous[indices], continuous_columns, packed, binary_columns)
    bundle = cls(shards, baseline, queries, None, config_snapshot=config,
                 identity_drops=MemoryDiagnostics(clean.dropped), diagnostics_in_memory=True)
    bundle.validate()
    return bundle
