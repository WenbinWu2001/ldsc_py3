"""Write an indexed batch from owned float64 chromosome fragments.

The producer owns the fragment directory until the writer finishes the batch.
Only one chromosome query frame is read, validated, and narrowed at a time;
the shared baseline and compact scientific summaries remain in memory.
"""

from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ._row_alignment import assert_same_snp_rows
from .errors import LDSCInternalError
from .ldscore_source import validate_ldscore_schemas
from .outputs import _cast_parquet_floats


@dataclass(frozen=True)
class StagedLDScoreBatch:
    """Metadata and ordered query paths for one private indexed batch.

    ``query_fragments`` contain float64 scores before public narrowing, in
    chromosome/SNP order. ``query_columns`` are finalized using float64
    variation summaries; fragment reads select only these retained columns.
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


def validate_fragment_baseline(writer, result):
    """Validate the shared baseline independently of deferred query reads."""
    writer._validate_tables(SimpleNamespace(**(vars(result) | {
        "query_table": None, "query_columns": [],
    })))


def write_query_fragments(writer, result, path, compression):
    """Validate and write one chromosome per row group without concatenation."""
    identity = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in result.baseline_table]
    columns = [*identity, *result.query_columns]
    row_groups, offset = [], 0
    parquet_writer = None
    with ExitStack() as stack:
        for fragment in result.query_fragments:
            frame = pd.read_parquet(fragment, columns=columns, use_threads=False)
            baseline = result.baseline_table.iloc[offset:offset + len(frame)]
            part = SimpleNamespace(baseline_table=baseline, query_table=frame,
                                   baseline_columns=result.baseline_columns,
                                   query_columns=result.query_columns,
                                   config_snapshot=result.config_snapshot)
            writer._validate_tables(part)
            validate_ldscore_schemas(result.baseline_columns, result.query_columns,
                                    baseline.columns, frame.columns)
            assert_same_snp_rows(baseline, frame, context="indexed query rows must match baseline rows on CHR/SNP/POS",
                                snp_identifier=result.config_snapshot.snp_identifier)
            if not frame.empty:
                if frame.CHR.astype(str).nunique() != 1:
                    raise LDSCInternalError("An indexed query fragment must contain exactly one chromosome.")
                narrowed = _cast_parquet_floats(frame)
                table = pa.Table.from_pandas(narrowed, preserve_index=False)
                if parquet_writer is None:
                    parquet_writer = stack.enter_context(pq.ParquetWriter(path, table.schema, compression=compression))
                parquet_writer.write_table(table, row_group_size=len(frame))
                row_groups.append({"chrom": str(frame.CHR.iloc[0]), "row_group_index": len(row_groups),
                                   "row_offset": offset, "n_rows": len(frame)})
                offset += len(frame)
                del table, narrowed
            del part, baseline, frame
            fragment.unlink()
    if offset != len(result.baseline_table):
        raise LDSCInternalError("Indexed query fragments do not cover every baseline SNP row.")
    return row_groups
