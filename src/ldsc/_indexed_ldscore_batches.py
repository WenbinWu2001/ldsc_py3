"""Assemble immutable-index queries with bounded chromosome and batch ownership.

Workers retain one chromosome operator and serialize float64 batch fragments in
run-owned scratch. The coordinator reads one genome-wide batch in chromosome
order. These private fragments are never accepted as public input artifacts.
"""

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import multiprocessing
from pathlib import Path
import pickle

import numpy as np
import pandas as pd


def _save(path, value):
    with path.open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)


def _read(path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def _write_chromosome(index, chrom, batch, declarations, control_resolution, width, root):
    from . import gene_ldscore_index as core

    directory = Path(root) / f"chr{chrom}"
    directory.mkdir()
    record = index.load_chromosome(chrom)
    control_selector = (None if control_resolution is None else core._selector_for_resolution(
        index.gene_catalog, record.atom_model, chrom, control_resolution))
    baseline = record.baseline_rows.copy()
    supplied = len(index.baseline_columns)
    n_baseline = supplied + int(control_selector is not None)
    counts = np.zeros((2, n_baseline), dtype=np.float64)
    blocks = np.zeros((2, n_baseline, n_baseline), dtype=np.float64)
    counts[:, :supplied] = [record.baseline_count_all, record.baseline_count_common]
    blocks[:, :supplied, :supplied] = [record.baseline_overlap_all, record.baseline_overlap_common]
    if control_selector is not None:
        baseline["gene_control"] = core.assemble_indexed_ld_scores(record.operator, control_selector)
        control = core.assemble_selected_atom_statistics(record.atom_statistics, control_selector)
        counts[:, supplied] = [control.count_all, control.count_common]
        cross = np.array([control.baseline_overlap_all, control.baseline_overlap_common])
        blocks[:, :supplied, supplied] = cross
        blocks[:, supplied, :supplied] = cross
        blocks[:, supplied, supplied] = counts[:, supplied]
    totals = np.array([record.total_reference_snps_all, record.total_reference_snps_common], dtype=np.float64)
    _save(directory / "baseline.pkl", (baseline, counts, blocks, totals))
    identity = baseline.loc[:, [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in baseline]]
    for ordinal, start in enumerate(range(0, len(declarations), width)):
        active = declarations[start:start + width]
        selectors = [core._selector_for_resolution(
            index.gene_catalog, record.atom_model, chrom,
            batch.selection(item["input_role"], item["source_ordinal"])) for item in active]
        scores = core.assemble_indexed_ld_scores(record.operator, np.column_stack(selectors))
        query_counts = np.empty((2, len(active)), dtype=np.float64)
        query_blocks = np.empty((2, n_baseline, len(active)), dtype=np.float64)
        for offset, selector in enumerate(selectors):
            selected = core.assemble_selected_atom_statistics(
                record.atom_statistics, selector, control_selector=control_selector)
            query_counts[:, offset] = [selected.count_all, selected.count_common]
            query_blocks[:, :supplied, offset] = [selected.baseline_overlap_all, selected.baseline_overlap_common]
            if control_selector is not None:
                query_blocks[:, supplied, offset] = [selected.control_overlap_all, selected.control_overlap_common]
        table = pd.concat([identity, pd.DataFrame(scores, columns=[item["query"] for item in active])], axis=1)
        _save(directory / f"query-{ordinal}.pkl", (table, query_counts, query_blocks))
        del table, scores, selectors, selector, selected, query_counts, query_blocks


def _write_chromosomes(index, batch, declarations, control, width, root, threads):
    from .ldscore_calculator import _init_worker
    from ._parallelism import _resolve_worker_count

    workers = _resolve_worker_count(threads, len(index.chromosomes))
    def arguments(chrom):
        return index, chrom, batch, declarations, control, width, root
    if workers == 1:
        for chrom in index.chromosomes:
            _write_chromosome(*arguments(chrom))
        return
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                             initializer=_init_worker, initargs=(None,)) as executor:
        remaining = iter(index.chromosomes)
        pending = {executor.submit(_write_chromosome, *arguments(next(remaining))) for _ in range(workers)}
        try:
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    future.result()
                for _ in done:
                    chrom = next(remaining, None)
                    if chrom is not None:
                        pending.add(executor.submit(_write_chromosome, *arguments(chrom)))
        except BaseException:
            for future in pending:
                future.cancel()
            raise


def indexed_results(index, batch, declarations, control, width, root, threads, statuses, scope):
    """Yield finalized float64 genome-wide results one query batch at a time."""
    from .config import GlobalConfig
    from .ldscore_calculator import LDScoreResult
    from .overlap_matrix import LDScoreOverlap
    from .query_annotations import finalize_query_statuses, _log_query_annotation_statuses

    root = Path(root) / "indexed-batches"
    root.mkdir()
    _write_chromosomes(index, batch, declarations, control, width, root, threads)
    baselines = [_read(root / f"chr{chrom}" / "baseline.pkl") for chrom in index.chromosomes]
    baseline = pd.concat([item[0] for item in baselines], ignore_index=True)
    baseline_counts = sum((item[1] for item in baselines))
    baseline_blocks = sum((item[2] for item in baselines))
    totals = sum((item[3] for item in baselines))
    del baselines
    baseline_columns = [*index.baseline_columns, *(["gene_control"] if control is not None else [])]
    config = GlobalConfig(snp_identifier="rsid") if index.snp_identifier == "rsid" else GlobalConfig(
        snp_identifier="chr_pos", genome_build=index.genome_build)
    for ordinal, start in enumerate(range(0, len(declarations), width)):
        query_columns = [item["query"] for item in declarations[start:start + width]]
        fragments = [_read(root / f"chr{chrom}" / f"query-{ordinal}.pkl") for chrom in index.chromosomes]
        query_table = pd.concat([item[0] for item in fragments], ignore_index=True)
        query_counts = sum((item[1] for item in fragments))
        query_blocks = sum((item[2] for item in fragments))
        del fragments
        all_columns = [*baseline_columns, *query_columns]
        block = np.concatenate([baseline_blocks, query_blocks], axis=2)
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
        active_statuses = tuple(status for status in statuses if status.query in query_columns or (
            ordinal == 0 and status.status not in {"ok", "warning"}))
        result = LDScoreResult(
            baseline_table=baseline, query_table=query_table, count_records=counts,
            baseline_columns=baseline_columns, query_columns=query_columns,
            ld_reference_snps=frozenset(), ld_regression_snps=frozenset(), chromosome_results=[],
            count_config={"common_reference_snp_maf_min": 0.05, "common_reference_snp_maf_operator": ">="},
            config_snapshot=config, overlap=overlap, query_statuses=active_statuses,
            gene_list_batch=batch, chromosome_scope=scope,
            index_provenance={"index_id": index.index_id, "index_snp_identifier": index.snp_identifier,
                              "index_genome_build": index.genome_build})
        result = finalize_query_statuses(result, active_statuses)
        _log_query_annotation_statuses(result.query_statuses)
        result.validate()
        del query_table
        yield result
        del result
        for chrom in index.chromosomes:
            (root / f"chr{chrom}" / f"query-{ordinal}.pkl").unlink()
