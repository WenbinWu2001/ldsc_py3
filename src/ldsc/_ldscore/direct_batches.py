"""Stage direct post-kernel query scores before worker-result transport.

Writing workers return compact scientific/diagnostic summaries, while their
float32 query scores stay in private Parquet until the shared writer consumes
them. Only the first query batch returns the chromosome's baseline table; later
batches return ``baseline_table=None`` and their query rows are checked against
the shared first-batch baseline. The caller owns the batch directory through
consumption or pool shutdown.
"""

from pathlib import Path

from .query_variation import summarize_variation
from .results import DirectChromosomeFragment


def stage_query_table(chrom, frame, query_columns, directory):
    """Persist existing float32 query scores and summarize their variation."""
    variation = summarize_variation(frame[name].to_numpy(copy=False) for name in query_columns)
    path = None
    if frame is not None:
        path = Path(directory) / f"chr{chrom}.parquet"
        frame.to_parquet(path, index=False, compression="zstd", use_dictionary=False)
    return path, variation


def stage_direct_chromosome(result, directory):
    """Persist existing float32 scores without recomputation or further casting."""
    path, variation = stage_query_table(result.chrom, result.query_table, result.query_columns, directory)
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
