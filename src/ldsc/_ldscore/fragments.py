"""Write a direct or indexed batch from owned chromosome fragments.

The producer owns the fragment directory until the writer finishes the batch.
Only one chromosome query frame is read, validated, and narrowed at a time;
the shared baseline and compact scientific summaries remain in memory.
"""

from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ..errors import LDSCInputError, LDSCInternalError
from .outputs import _cast_parquet_floats


@dataclass(frozen=True)
class StagedLDScoreBatch:
    """Metadata and ordered query paths for one private LD-score batch.

    ``query_fragments`` retain direct post-kernel float32 scores or indexed
    float64 scores before public narrowing, in chromosome/SNP order.
    ``query_columns`` are finalized using variation summaries in the respective
    LD-score producer's precision; reads select only these retained columns.
    This internal descriptor is neither a materialized ``LDScoreResult`` nor
    a persistent ``LDScoreSource`` and must not escape its producer lifetime.
    """

    baseline_table: pd.DataFrame
    query_fragments: tuple[Path, ...]
    baseline_columns: list[str]
    query_columns: list[str]
    count_records: list[dict]
    count_config: dict
    config_snapshot: Any
    overlap: Any
    query_statuses: tuple
    gene_list_batch: Any
    chromosome_scope: dict
    index_provenance: dict
    annotation_types: dict = field(default_factory=dict)
    identity_drops_by_chrom: dict = field(default_factory=dict)
    snp_universe_policy: dict | None = None
    ld_window_snp_counts: pd.DataFrame | None = None
    reference_snp_funnel: pd.DataFrame | None = None


def write_query_fragments(result: StagedLDScoreBatch, baseline_table: pd.DataFrame, path: Path, compression):
    """Write finalized queries against the run's already-validated shared baseline.

    ``baseline_table`` is the first batch's baseline. Read one chromosome at a
    time, require its identity columns to equal the matching baseline rows
    exactly, narrow its values to float32, and delete the consumed fragment
    after releasing its frame. The producer keeps the fragment directory alive
    until consumption or worker shutdown.
    """
    identity = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in baseline_table]
    columns = [*identity, *result.query_columns]
    row_groups, offset = [], 0
    parquet_writer = None
    with ExitStack() as stack:
        for fragment in result.query_fragments:
            frame = pd.read_parquet(fragment, columns=columns, use_threads=False)
            _require_baseline_rows(baseline_table.iloc[offset:offset + len(frame)], frame, identity, offset)
            if not frame.empty:
                if frame.CHR.astype(str).nunique() != 1:
                    raise LDSCInternalError("A query fragment must contain exactly one chromosome.")
                narrowed = _cast_parquet_floats(frame)
                table = pa.Table.from_pandas(narrowed, preserve_index=False)
                if parquet_writer is None:
                    parquet_writer = stack.enter_context(pq.ParquetWriter(path, table.schema, compression=compression))
                parquet_writer.write_table(table, row_group_size=len(frame))
                row_groups.append({"chrom": str(frame.CHR.iloc[0]), "row_group_index": len(row_groups),
                                   "row_offset": offset, "n_rows": len(frame)})
                offset += len(frame)
                del table, narrowed
            del frame
            fragment.unlink()
    if offset != len(baseline_table):
        raise LDSCInternalError("Query fragments do not cover every baseline SNP row.")
    return row_groups


def _require_baseline_rows(baseline: pd.DataFrame, frame: pd.DataFrame, identity: list[str], offset: int) -> None:
    """Raise unless a query fragment repeats the baseline rows' identity columns exactly, in order."""
    for column in identity:
        expected, actual = baseline[column].to_numpy(), frame[column].to_numpy()
        same = np.zeros(len(actual), dtype=bool) if len(expected) != len(actual) else (
            (expected == actual) | (pd.isna(expected) & pd.isna(actual)))
        if not same.all():
            row = int(np.flatnonzero(~same)[0]) if len(same) else 0
            found = actual[row] if row < len(actual) else "no row"
            wanted = expected[row] if row < len(expected) else "no row"
            raise LDSCInputError(
                f"LD-score query batch rows must equal the baseline rows: SNP row mismatch at row {offset + row} "
                f"in {column} (baseline {wanted!r}, query batch {found!r}). Most likely an input file "
                "(reference panel, annotations or SNP restriction) changed while the run was in progress, so a "
                "later query batch selected different SNPs. Rerun with unchanged inputs."
            )
