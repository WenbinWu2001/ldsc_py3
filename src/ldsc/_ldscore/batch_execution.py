"""Run LD-score query batches for direct and indexed producers.

An LD-score producer computes one chromosome per task for the active query
batch and reduces the chromosome results to one staged batch; this module owns
everything around that. Query columns are split into batches of at most
``width``. Only one batch is in flight: its chromosome tasks run in the
caller's persistent worker pool, the parent reduces their compact results and
finalizes the batch's query statuses, and the writer consumes the staged batch
before the next batch starts. A failed task shuts the pool down, waiting for
running workers, before the batch's scratch is removed.
"""

from __future__ import annotations

from concurrent.futures import Future
from concurrent.futures.process import BrokenProcessPool
from dataclasses import replace
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, ContextManager, Iterator, Protocol, Sequence

from .._annotation.query_outcomes import finalize_query_measurements, finalize_query_statuses
from .._runtime.parallelism import run_chromosome_tasks
from .._runtime.progress import PhaseProgress
from .diagnostics import report_chromosome_failure
from .fragments import StagedLDScoreBatch
from .query_variation import ScoreVariation, combine_variation


class LDScoreProducer(Protocol):
    """Direct or indexed source of the chromosome results for one query batch.

    ``batch_name`` and ``phase`` label batches in chromosome and progress logs;
    ``task`` names one chromosome's work in the start and finish log lines.
    """

    chromosomes: Sequence[str]
    log_paths: dict
    logger: logging.Logger
    batch_name: str
    phase: str
    task: str

    def open_batch(self, ordinal: int, query_columns: list[str], statuses: tuple,
                   scratch: Path) -> ContextManager[Any]:
        """Prepare one batch's inputs; they stay live until the batch is consumed."""

    def submit(self, pool, chrom: str, batch: Any, label: str) -> Future:
        """Submit one chromosome task for ``batch`` and return its future."""

    def reduce(self, results: dict[str, Any], batch: Any) -> tuple[Any, list[ScoreVariation]]:
        """Combine chromosome results into one staged batch and its variation summaries."""


def run_chromosome_batch(producer: LDScoreProducer, pool, workers: int, batch: Any, label: str) -> dict[str, Any]:
    """Run one batch's chromosome tasks and return their results by chromosome.

    At most ``workers`` tasks are outstanding. On any failure the remaining
    tasks are cancelled, the failure is reported to the workflow log, and the
    pool shuts down, waiting for running workers, before the error propagates
    to the caller that owns the batch's scratch.
    """
    results = {}

    def submit(chrom):
        producer.logger.info("Starting chromosome %s %s; %s; log: %s.", chrom, producer.task, label,
                             producer.log_paths.get(chrom, "no file output"))
        return producer.submit(pool, chrom, batch, label)

    def collect(chrom, result):
        results[chrom] = result
        producer.logger.info("Finished chromosome %s %s (%d/%d); %s.", chrom, producer.task,
                             len(results), len(producer.chromosomes), label)

    def on_failure(error, failed_chrom, outstanding):
        affected = (outstanding if isinstance(error, BrokenProcessPool)
                    else [failed_chrom] if failed_chrom is not None else outstanding)
        report_chromosome_failure(error, affected, producer.log_paths, label, producer.logger)

    try:
        run_chromosome_tasks(producer.chromosomes, workers, submit, collect, on_failure)
    except BaseException:
        pool.shutdown(wait=True)
        raise
    return results


def run_query_batches(producer: LDScoreProducer, query_columns: Sequence[str], statuses: tuple, *,
                      width: int, pool, workers: int, scratch_root) -> Iterator[StagedLDScoreBatch]:
    """Yield one finalized staged batch at a time, in query order.

    The first batch also carries the statuses of queries skipped before
    calculation. With no query columns, one baseline-only batch runs. Each
    batch's scratch directory is removed after the consumer resumes the
    generator, so the consumer reads the batch's fragments before they go.
    """
    names = list(query_columns)
    groups = [names[start:start + width] for start in range(0, len(names), width)] or [[]]
    for ordinal, columns in enumerate(groups):
        label = f"{producer.batch_name} {ordinal + 1}/{len(groups)}"
        with PhaseProgress(producer.logger, "batch", f"{producer.phase} {ordinal + 1}/{len(groups)}",
                           len(columns)) as progress:
            selected = set(columns)
            active = tuple(status for status in statuses if status.query in selected
                           or (ordinal == 0 and status.status not in {"ok", "warning"}))
            with TemporaryDirectory(prefix=".ldsc-fragments-", dir=scratch_root) as scratch, \
                    producer.open_batch(ordinal, columns, active, Path(scratch)) as batch:
                results = run_chromosome_batch(producer, pool, workers, batch, label)
                staged, variation = producer.reduce(results, batch)
                del results
                yield finalize_batch(staged, active, variation)
                del staged, batch
            progress.completed = len(columns)


def finalize_batch(staged, statuses: tuple, variation: Sequence[ScoreVariation]):
    """Apply final query statuses from genome-wide counts and score variation.

    Queries with no annotation SNPs or constant LD scores are dropped from the
    batch's columns, count records, overlap, and annotation types. A batch
    without query statuses (prebuilt query annotations) is returned unchanged.
    A materialized ``LDScoreResult``, from the no-output direct path, is
    finalized from its own query table and needs no ``variation``.
    """
    if not statuses:
        return staged
    if not isinstance(staged, StagedLDScoreBatch):
        return finalize_query_statuses(staged, statuses)
    varying = dict(zip(staged.query_columns, combine_variation(variation, len(staged.query_columns)), strict=True))
    counts = {record["column"]: float(record["all_reference_snp_count"])
              for record in staged.count_records if record["group"] == "query"}
    finalized, retained = finalize_query_measurements(statuses, counts, varying)
    kept = {*staged.baseline_columns, *retained}
    return replace(
        staged, query_columns=retained, query_statuses=finalized,
        count_records=[record for record in staged.count_records
                       if record["group"] != "query" or record["column"] in retained],
        overlap=None if staged.overlap is None else staged.overlap.select_queries(retained),
        annotation_types={name: kind for name, kind in staged.annotation_types.items() if name in kept},
    )
