"""Selective access to canonical aggregate LD-score directories.

Shared baseline values and SNP metadata stay resident for regression alignment.
Query reads return detached batches; the source never caches query values.
The original directory must remain available and unchanged while it is used.
"""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .annotation_semantics import require_unique_annotation_names
from .config import GlobalConfig
from .errors import LDSCInputError, LDSCInternalError


def validate_ldscore_schemas(baseline_columns, query_columns, baseline_schema, query_schema):
    """Check every declared column without loading its values."""
    require_unique_annotation_names(baseline_columns, query_columns)
    required = {"CHR", "SNP", "POS", "regression_ld_scores", *baseline_columns}
    missing = required.difference(baseline_schema)
    if missing:
        raise LDSCInternalError(f"LD-score baseline_table is missing required columns {sorted(missing)}. Regenerate the complete LD-score directory.")
    if bool(query_columns) != (query_schema is not None):
        raise LDSCInternalError("LD-score query columns and query_table presence disagree. Regenerate the complete LD-score directory.")
    if query_schema is not None:
        missing = {"CHR", "SNP", "POS", *query_columns}.difference(query_schema)
        if missing:
            raise LDSCInternalError(f"LD-score query_table is missing required columns {sorted(missing)}. Regenerate the complete LD-score directory.")


@dataclass(frozen=True)
class LDScoreSource:
    """Validated shared LD metadata and baseline values, with query paths.

    ``read_queries(columns)`` loads only explicit query columns and returns
    their SNP identities in stored order. Alignment with baseline rows remains
    a fit-time check. Complete artifact schemas and allele metadata are checked
    when the source is opened, including columns that no model will select.
    """

    baseline_table: pd.DataFrame
    query_metadata: pd.DataFrame | None
    query_path: Path | None
    baseline_columns: list[str]
    query_columns: list[str]
    count_records: list[dict]
    config_snapshot: GlobalConfig | None
    count_config: dict
    overlap: object | None
    output_paths: dict[str, str]
    annotation_types: dict[str, str]
    ld_regression_snps: frozenset[str]
    chromosome_results: tuple = ()
    ld_reference_snps: frozenset[str] = frozenset()

    def read_queries(self, columns):
        """Return one detached query batch, including its stored SNP metadata."""
        columns = list(columns)
        missing = set(columns).difference(self.query_columns)
        if missing:
            raise LDSCInputError(f"Unknown query LD-score columns: {sorted(missing)}.")
        if self.query_path is None:
            if columns:
                raise LDSCInputError("LD-score directory has no query table.")
            return None
        values = pd.read_parquet(self.query_path, columns=columns)
        return pd.concat([self.query_metadata, values], axis=1)
