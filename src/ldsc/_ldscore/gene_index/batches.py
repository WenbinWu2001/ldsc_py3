"""The indexed LD-score producer: one stored chromosome operator per task.

``_ldscore.batch_execution`` runs one query batch at a time through a
persistent pool. Each task reloads one chromosome operator and writes lossless
float64 Parquet in batch-owned scratch; only compact measurements return to the
parent. The output writer consumes chromosome fragments in canonical order;
no genome-wide query table or future-batch score fragments are retained.
"""

from contextlib import contextmanager
from dataclasses import dataclass
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from ..._runtime.parallelism import chromosome_pool, resolve_worker_count
from ..._runtime.worker_runtime import initialize_worker
from ..batch_execution import run_query_batches
from ..diagnostics import chromosome_log
from ..query_variation import ScoreVariation, summarize_variation
from .atoms import assemble_atom_selector, assemble_indexed_ld_scores, assemble_selected_atom_statistics

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


_SCORE_ROW_BLOCK = 65536


def _write_chromosome(loader, chrom, query_names, query_rows, control_rows, ordinal, total_batches, root,
                      log_path=None):
    """Load, compute, and release one chromosome for the active query batch.

    ``query_rows`` and ``control_rows`` are this chromosome's included gene
    rows of each query and of the control, resolved by the parent; a task
    carries neither the index catalog nor the gene-list batch.
    """
    with chromosome_log(log_path, chrom, f"indexed batch {ordinal + 1}/{total_batches}"):
        directory = Path(root) / f"chr{chrom}"
        directory.mkdir()
        record = loader.load(chrom)
        gene_to_atom = record.atom_model.gene_to_atom
        control_selector = None if control_rows is None else assemble_atom_selector(gene_to_atom, control_rows)
        baseline_statistics = _write_baseline(record, control_selector, directory) if ordinal == 0 else None
        supplied = loader.n_baseline_columns
        n_baseline = supplied + int(control_selector is not None)
        selectors = [assemble_atom_selector(gene_to_atom, rows) for rows in query_rows]
        del gene_to_atom
        query_counts = np.empty((2, len(query_names)), dtype=np.float64)
        query_blocks = np.empty((2, n_baseline, len(query_names)), dtype=np.float64)
        for offset, selector in enumerate(selectors):
            selected = assemble_selected_atom_statistics(
                record.atom_statistics, selector, control_selector=control_selector)
            query_counts[:, offset] = [selected.count_all, selected.count_common]
            query_blocks[:, :supplied, offset] = [selected.baseline_overlap_all, selected.baseline_overlap_common]
            if control_selector is not None:
                query_blocks[:, supplied, offset] = [selected.control_overlap_all, selected.control_overlap_common]
        # Column-major scores, filled by row blocks of the same product (each
        # row's arithmetic is unchanged), give Arrow contiguous columns to
        # share instead of copying.
        selector_matrix = np.column_stack(selectors)
        del selectors, selector, selected
        scores = np.empty((record.operator.shape[0], selector_matrix.shape[1]), dtype=np.float64, order="F")
        for start in range(0, len(scores), _SCORE_ROW_BLOCK):
            stop = start + _SCORE_ROW_BLOCK
            scores[start:stop] = assemble_indexed_ld_scores(record.operator[start:stop], selector_matrix)
        del selector_matrix
        variation = summarize_variation(scores.T)
        metadata = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in record.baseline_rows]
        identity = pa.Table.from_pandas(record.baseline_rows.loc[:, metadata], preserve_index=False)
        table = pa.Table.from_arrays([*identity.columns, *[pa.array(column) for column in scores.T]],
                                    names=[*metadata, *query_names])
        # Arrow shares the score columns now; compression does not need the
        # operator alongside them.
        del scores, identity, record, control_selector
        pq.write_table(table, directory / "query.parquet", compression="zstd", use_dictionary=False)
        return _ChromosomeBatchSummary(query_counts, query_blocks, variation, baseline_statistics)


class _IndexedProducer:
    """Indexed LD-score producer: each task projects one chromosome's stored operator.

    The first batch also stages and combines the shared baseline and its
    statistics, which later batches reuse. ``reduce`` combines the compact
    float64 summaries in index order; query scores stay in the workers'
    fragments.
    """

    batch_name, phase, task = "indexed batch", "indexed LD-score assembly", "indexed LD-score computation"

    def __init__(self, index, gene_batch, declarations, control, width, scope, log_paths):
        from .scoring import _gene_rows_by_chromosome

        self.index = index
        self.loader = index.chromosome_loader()
        self.gene_batch = gene_batch
        self.declarations = declarations
        self.control = control
        self.control_rows = (None if control is None else
                             _gene_rows_by_chromosome(index.gene_catalog, control, index.chromosomes))
        self.scope = scope
        self.log_paths = log_paths or {}
        self.logger = LOGGER
        self.chromosomes = index.chromosomes
        self.total_batches = (len(declarations) + width - 1) // width
        self.baseline_columns = [*index.baseline_columns, *(["gene_control"] if control is not None else [])]
        self.baseline = self.baseline_statistics = None

    @contextmanager
    def open_batch(self, ordinal, query_columns, statuses, scratch):
        """Select the batch's gene-list declarations and resolve their gene rows per chromosome.

        Workers write into ``scratch``.
        """
        from .scoring import _gene_rows_by_chromosome

        selected = set(query_columns)
        declarations = [item for item in self.declarations if item["query"] in selected]
        rows = {chrom: [] for chrom in self.chromosomes}
        for item in declarations:
            by_chromosome = _gene_rows_by_chromosome(
                self.index.gene_catalog, self.gene_batch.selection(item["input_role"], item["source_ordinal"]),
                self.chromosomes)
            for chrom in self.chromosomes:
                rows[chrom].append(by_chromosome[chrom])
        yield SimpleNamespace(ordinal=ordinal, root=scratch, declarations=declarations, rows=rows)

    def submit(self, pool, chrom, batch, label):
        """Submit one chromosome's projection for the batch's declarations."""
        return pool.submit(_write_chromosome, self.loader, chrom, [item["query"] for item in batch.declarations],
                           batch.rows[chrom], None if self.control_rows is None else self.control_rows[chrom],
                           batch.ordinal, self.total_batches, batch.root, self.log_paths.get(chrom))

    def reduce(self, results, batch):
        """Sum chromosome counts and overlap blocks into one staged batch."""
        from ..fragments import StagedLDScoreBatch
        from ..overlap_table import LDScoreOverlap

        index, directory = self.index, batch.root
        summaries = [results[chrom] for chrom in self.chromosomes]
        if batch.ordinal == 0:
            self.baseline = pd.concat([pd.read_parquet(directory / f"chr{chrom}" / "baseline.parquet", use_threads=False)
                                       for chrom in self.chromosomes], ignore_index=True)
            self.baseline_statistics = [sum(item.baseline_statistics[pos] for item in summaries) for pos in range(3)]
        baseline_counts, baseline_blocks, totals = self.baseline_statistics
        baseline_columns = self.baseline_columns
        query_columns = [item["query"] for item in batch.declarations]
        query_counts = sum(item.counts for item in summaries)
        query_blocks = sum(item.blocks for item in summaries)
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
                  for pos, name in enumerate(columns)]
        staged = StagedLDScoreBatch(
            baseline_table=self.baseline,
            query_fragments=tuple(directory / f"chr{chrom}" / "query.parquet" for chrom in self.chromosomes),
            count_records=counts, baseline_columns=baseline_columns, query_columns=query_columns,
            count_config=index.count_config, config_snapshot=index.global_config, overlap=overlap,
            query_statuses=(), gene_list_batch=self.gene_batch, chromosome_scope=self.scope,
            index_provenance={"index_id": index.index_id, "index_snp_identifier": index.snp_identifier,
                              "index_genome_build": index.genome_build},
            reference_snp_funnel=index.reference_snp_funnel)
        return staged, [item.variation for item in summaries]


def indexed_results(index, batch, declarations, control, width, root, workers, statuses, scope, log_paths=None):
    """Yield one file-backed batch; release its fragments before submitting more.

    Workers retain one chromosome's float64 arrays. The parent retains the
    shared baseline and compact summaries, while the writer reads one query
    chromosome at a time. Index loading repeats for each batch. Batch
    scratch lives under ``root``.
    """
    workers = resolve_worker_count(workers, len(index.chromosomes))
    producer = _IndexedProducer(index, batch, declarations, control, width, scope, log_paths)
    with chromosome_pool(workers, operation="ldscore chromosome computation", initializer=initialize_worker,
                         initargs=(logging.getLogger("LDSC").getEffectiveLevel(),)) as pool:
        yield from run_query_batches(producer, [item["query"] for item in declarations], statuses, width=width,
                                     pool=pool, workers=workers, scratch_root=root)
