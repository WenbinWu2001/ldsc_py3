"""The staged batch a direct or indexed LD-score producer hands the writer.

The LD-score producer owns the batch's chromosome fragment directory until the
writer (``_ldscore.outputs``) has consumed it, reading, validating, and
narrowing one chromosome query frame at a time; the shared baseline and
compact scientific summaries remain in memory.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


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
    ld_window: dict | None = None
    ld_window_snp_counts: pd.DataFrame | None = None
    reference_snp_funnel: pd.DataFrame | None = None
