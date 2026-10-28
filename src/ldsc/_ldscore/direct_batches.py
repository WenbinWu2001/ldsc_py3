"""Stage direct post-kernel query scores before worker-result transport.

The numerical workflow still computes complete chromosomes. Writing workers
return their baseline and compact scientific/diagnostic summaries, while their
float32 query scores stay in private Parquet until the shared writer consumes
them. The caller owns the batch directory through consumption or pool shutdown.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .query_variation import ScoreVariation, summarize_variation


@dataclass(frozen=True)
class DirectChromosomeFragment:
    """Direct chromosome statistics with a private query path instead of scores."""

    chrom: str
    baseline_table: pd.DataFrame
    query_fragment: Path | None
    baseline_columns: list[str]
    query_columns: list[str]
    snp_count_totals: dict
    config_snapshot: Any
    overlap: Any
    reference_snp_count: int
    regression_selected_snp_count: int
    regression_region_removed_snp_count: int
    identity_drops: Any
    annotation_types: dict
    variation: ScoreVariation
    reference_snp_funnel: pd.DataFrame | None = None


def stage_direct_chromosome(result, directory):
    """Persist existing float32 scores without recomputation or further casting."""
    frame = result.query_table
    variation = summarize_variation(
        frame[name].to_numpy(copy=False) for name in result.query_columns)
    path = None
    if frame is not None:
        path = Path(directory) / f"chr{result.chrom}.parquet"
        frame.to_parquet(path, index=False, compression="zstd", use_dictionary=False)
    return DirectChromosomeFragment(
        chrom=result.chrom, baseline_table=result.baseline_table, query_fragment=path,
        baseline_columns=result.baseline_columns, query_columns=result.query_columns,
        snp_count_totals=result.snp_count_totals, config_snapshot=result.config_snapshot,
        overlap=result.overlap, reference_snp_count=result.reference_snp_count,
        regression_selected_snp_count=result.regression_selected_snp_count,
        regression_region_removed_snp_count=result.regression_region_removed_snp_count,
        identity_drops=result.identity_drops, annotation_types=result.annotation_types,
        variation=variation, reference_snp_funnel=result.reference_snp_funnel,
    )
