"""Compute one indexed query batch before advancing to the next.

A persistent pool reloads one chromosome operator per task and writes lossless
float64 Parquet in batch-owned scratch. Only compact measurements return to the
parent. The output writer consumes chromosome fragments in canonical order;
no genome-wide query table or future-batch score fragments are retained.
"""

from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass
import logging
import multiprocessing
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ..._runtime.progress import PhaseProgress
from ..._runtime.parallelism import _resolve_worker_count, _run_chromosome_tasks
from ..diagnostics import chromosome_log, report_chromosome_failure
from ..query_variation import ScoreVariation, summarize_variation, combine_variation
from .atoms import assemble_indexed_ld_scores, assemble_selected_atom_statistics

LOGGER = logging.getLogger("LDSC.gene_ldscore_index")


@dataclass(frozen=True)
class _ChromosomeBatchSummary:
    """Small float64 statistics; score arrays stay in the worker's fragment."""

    counts: np.ndarray
    blocks: np.ndarray
    variation: ScoreVariation
    baseline_statistics: tuple | None


def _write_baseline(record, control_selector, directory):
    """Stage shared baseline scores once and return all/common statistics."""
    supplied = len(record.baseline_count_all)
    n_baseline = supplied + int(control_selector is not None)
    counts = np.zeros((2, n_baseline), dtype=np.float64)
    blocks = np.zeros((2, n_baseline, n_baseline), dtype=np.float64)
    counts[:, :supplied] = [record.baseline_count_all, record.baseline_count_common]
    blocks[:, :supplied, :supplied] = [record.baseline_overlap_all, record.baseline_overlap_common]
    baseline = record.baseline_rows
    if control_selector is not None:
        baseline = baseline.copy()
        baseline["gene_control"] = assemble_indexed_ld_scores(record.operator, control_selector)
        control = assemble_selected_atom_statistics(record.atom_statistics, control_selector)
        counts[:, supplied] = [control.count_all, control.count_common]
        cross = np.array([control.baseline_overlap_all, control.baseline_overlap_common])
        blocks[:, :supplied, supplied] = cross
        blocks[:, supplied, :supplied] = cross
        blocks[:, supplied, supplied] = counts[:, supplied]
    baseline.to_parquet(directory / "baseline.parquet", index=False, compression="zstd")
    totals = np.array([record.total_reference_snps_all, record.total_reference_snps_common], dtype=np.float64)
    return counts, blocks, totals


def _write_chromosome(index, chrom, batch, declarations, control_resolution, ordinal, total_batches, root, log_path=None):
    """Load, compute, and release one chromosome for the active query batch."""
    from .scoring import _selector_for_resolution

    with chromosome_log(log_path, chrom, f"indexed batch {ordinal + 1}/{total_batches}"):
        directory = Path(root) / f"chr{chrom}"
        directory.mkdir()
        record = index.load_chromosome(chrom)
        control_selector = (None if control_resolution is None else _selector_for_resolution(
            index.gene_catalog, record.atom_model, chrom, control_resolution))
        baseline_statistics = _write_baseline(record, control_selector, directory) if ordinal == 0 else None
        supplied = len(index.baseline_columns)
        n_baseline = supplied + int(control_selector is not None)
        selectors = [_selector_for_resolution(
            index.gene_catalog, record.atom_model, chrom,
            batch.selection(item["input_role"], item["source_ordinal"])) for item in declarations]
        query_counts = np.empty((2, len(declarations)), dtype=np.float64)
        query_blocks = np.empty((2, n_baseline, len(declarations)), dtype=np.float64)
        for offset, selector in enumerate(selectors):
            selected = assemble_selected_atom_statistics(
                record.atom_statistics, selector, control_selector=control_selector)
            query_counts[:, offset] = [selected.count_all, selected.count_common]
            query_blocks[:, :supplied, offset] = [selected.baseline_overlap_all, selected.baseline_overlap_common]
            if control_selector is not None:
                query_blocks[:, supplied, offset] = [selected.control_overlap_all, selected.control_overlap_common]
        scores = assemble_indexed_ld_scores(record.operator, np.column_stack(selectors))
        del selectors, selector, selected
        variation = summarize_variation(scores.T)
        metadata = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in record.baseline_rows]
        identity = pa.Table.from_pandas(record.baseline_rows.loc[:, metadata], preserve_index=False)
        table = pa.Table.from_arrays([*identity.columns, *[pa.array(column) for column in scores.T]],
                                    names=[*metadata, *[item["query"] for item in declarations]])
        # Arrow owns the columns now; compression does not need the operator
        # or the original row-major score matrix alongside its column buffers.
        del scores, identity, record, control_selector
        pq.write_table(table, directory / "query.parquet", compression="zstd", use_dictionary=False)
        return _ChromosomeBatchSummary(query_counts, query_blocks, variation, baseline_statistics)


def _write_chromosomes(executor, index, batch, declarations, control, ordinal, total_batches, root, workers, log_paths):
    """Collect compact measurements in index order using the shared scheduler."""
    summaries = {}
    stage = f"indexed batch {ordinal + 1}/{total_batches}"

    def submit(chrom):
        LOGGER.info("Starting chromosome %s indexed LD-score computation; batch %d/%d.", chrom, ordinal + 1, total_batches)
        return executor.submit(_write_chromosome, index, chrom, batch, declarations, control,
                               ordinal, total_batches, root, log_paths.get(chrom))

    def collect(chrom, result):
        summaries[chrom] = result
        LOGGER.info("Finished chromosome %s indexed LD-score computation; batch %d/%d.", chrom, ordinal + 1, total_batches)

    def on_failure(error, failed_chrom, outstanding):
        report_chromosome_failure(error, outstanding if isinstance(error, BrokenProcessPool)
                                  else [failed_chrom] if failed_chrom is not None else outstanding,
                                  log_paths, stage, LOGGER)

    _run_chromosome_tasks(index.chromosomes, workers, submit, collect, on_failure)
    return [summaries[chrom] for chrom in index.chromosomes]


def indexed_results(index, batch, declarations, control, width, root, workers, statuses, scope, log_paths=None):
    """Yield one file-backed batch; release its fragments before submitting more.

    Workers retain one chromosome's float64 arrays. The parent retains the
    shared baseline and compact summaries, while the writer reads one query
    chromosome at a time. Index loading repeats for each batch. A failed task
    must finish pool shutdown before its batch workspace is removed.
    """
    from ...config import GlobalConfig
    from ..calculator import _init_worker
    from ..fragments import StagedLDScoreBatch
    from ..overlap_table import LDScoreOverlap
    from ..._annotation.query_outcomes import finalize_query_measurements

    root = Path(root) / "indexed-batches"
    root.mkdir()
    workers = _resolve_worker_count(workers, len(index.chromosomes))
    baseline_columns = [*index.baseline_columns, *(["gene_control"] if control is not None else [])]
    config = GlobalConfig(snp_identifier="rsid") if index.snp_identifier == "rsid" else GlobalConfig(
        snp_identifier="chr_pos", genome_build=index.genome_build)
    total_batches = (len(declarations) + width - 1) // width
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                             initializer=_init_worker,
                             initargs=(None, None, logging.getLogger("LDSC").getEffectiveLevel())) as executor:
        for ordinal, start in enumerate(range(0, len(declarations), width)):
            active = declarations[start:start + width]
            with TemporaryDirectory(prefix=f"batch-{ordinal + 1:05d}-", dir=root) as temporary:
                directory = Path(temporary)
                try:
                    summaries = _write_chromosomes(executor, index, batch, active, control, ordinal,
                                                   total_batches, directory, workers, log_paths or {})
                except BaseException:
                    executor.shutdown(wait=True)
                    raise
                if ordinal == 0:
                    baseline = pd.concat([pd.read_parquet(directory / f"chr{chrom}" / "baseline.parquet",
                                                         use_threads=False) for chrom in index.chromosomes], ignore_index=True)
                    baseline_counts, baseline_blocks, totals = [sum(item.baseline_statistics[pos] for item in summaries)
                                                               for pos in range(3)]
                with PhaseProgress(LOGGER, "batch", f"indexed LD-score assembly {ordinal + 1}/{total_batches}",
                                   len(active)) as progress:
                    query_columns = [item["query"] for item in active]
                    query_counts = sum(item.counts for item in summaries)
                    query_blocks = sum(item.blocks for item in summaries)
                    varying = dict(zip(query_columns, combine_variation(
                        (item.variation for item in summaries), len(active)), strict=True))
                    active_statuses = tuple(status for status in statuses if status.query in varying or (
                        ordinal == 0 and status.status not in {"ok", "warning"}))
                    finalized, retained = finalize_query_measurements(
                        active_statuses, dict(zip(query_columns, query_counts[0], strict=True)), varying)
                    block = np.concatenate([baseline_blocks, query_blocks], axis=2)
                    all_columns = [*baseline_columns, *query_columns]
                    overlap = LDScoreOverlap(
                        baseline_block_all=pd.DataFrame(block[0], index=baseline_columns, columns=all_columns),
                        baseline_block_common=pd.DataFrame(block[1], index=baseline_columns, columns=all_columns),
                        query_diagonal_all=pd.Series(query_counts[0], index=query_columns),
                        query_diagonal_common=pd.Series(query_counts[1], index=query_columns),
                        total_all_reference_snps=float(totals[0]), total_common_reference_snps=float(totals[1]))
                    counts = [{"group": group, "column": name, "all_reference_snp_count": float(values[0, pos]),
                               "common_reference_snp_count": float(values[1, pos])}
                              for group, columns, values in [("baseline", baseline_columns, baseline_counts),
                                                             ("query", query_columns, query_counts)]
                              for pos, name in enumerate(columns) if group == "baseline" or name in retained]
                    result = StagedLDScoreBatch(
                        baseline_table=baseline,
                        query_fragments=tuple(directory / f"chr{chrom}" / "query.parquet" for chrom in index.chromosomes),
                        count_records=counts, baseline_columns=baseline_columns, query_columns=retained,
                        count_config={"common_reference_snp_maf_min": 0.05, "common_reference_snp_maf_operator": ">="},
                        config_snapshot=config, overlap=overlap.select_queries(retained), query_statuses=finalized,
                        gene_list_batch=batch, chromosome_scope=scope,
                        index_provenance={"index_id": index.index_id, "index_snp_identifier": index.snp_identifier,
                                          "index_genome_build": index.genome_build})
                    yield result
                    del result
                    progress.completed = len(active)
