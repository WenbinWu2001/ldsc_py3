"""Query batches run one at a time through the caller's pool for any producer."""

from concurrent.futures import Future, ThreadPoolExecutor
from contextlib import contextmanager
import logging
from pathlib import Path
from threading import Event
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ldsc._annotation.query_outcomes import QueryAnnotationStatus
from ldsc._ldscore.batch_execution import run_query_batches
from ldsc._ldscore.fragments import StagedLDScoreBatch
from ldsc._ldscore.query_variation import summarize_variation


class FakeProducer:
    """Each chromosome scores every query 1, 2, 3 except ``flat``, scored 1, 1, 1."""

    batch_name, phase, task = "batch", "fake LD-score calculation", "fake LD-score computation"

    def __init__(self, task=None, chromosomes=("1", "2", "3")):
        self.chromosomes = chromosomes
        self.log_paths = {}
        self.logger = logging.getLogger("LDSC.test_batch_execution")
        self.task_function = task or (lambda chrom, batch: chrom)
        self.events, self.pools = [], []

    @contextmanager
    def open_batch(self, ordinal, query_columns, statuses, scratch):
        self.events.append(("open", ordinal, tuple(query_columns), tuple(s.query for s in statuses)))
        (scratch / "borrowed.txt").write_text("borrowed by tasks")
        try:
            yield SimpleNamespace(ordinal=ordinal, query_columns=list(query_columns), scratch=scratch)
        finally:
            self.events.append(("close", ordinal))

    def submit(self, pool, chrom, batch, label):
        self.pools.append(pool)
        return pool.submit(self.task_function, chrom, batch)

    def reduce(self, results, batch):
        scores = {name: np.array([1., 1., 1.]) if name == "flat" else np.array([1., 2., 3.])
                  for name in batch.query_columns}
        variation = [summarize_variation(scores[name] for name in batch.query_columns) for _ in results]
        counts = [{"group": "query", "column": name, "all_reference_snp_count": 3.0,
                   "common_reference_snp_count": 3.0} for name in batch.query_columns]
        staged = StagedLDScoreBatch(
            baseline_table=pd.DataFrame({"base": [1.]}), query_fragments=(), baseline_columns=["base"],
            query_columns=list(batch.query_columns), count_records=counts, count_config={},
            config_snapshot=None, overlap=None, query_statuses=(), gene_list_batch=None,
            chromosome_scope={}, index_provenance={},
            annotation_types={"base": "binary", **{name: "binary" for name in batch.query_columns}})
        return staged, variation


def ok(query):
    return QueryAnnotationStatus(query, f"{query}.txt", "gene_list", "ok")


def test_batches_share_the_callers_pool_and_run_one_at_a_time(tmp_path):
    producer = FakeProducer()
    queries = ["q1", "q2", "q3", "q4", "q5"]
    with ThreadPoolExecutor(max_workers=2) as pool:
        for staged in run_query_batches(producer, queries, tuple(map(ok, queries)), width=2, pool=pool,
                                        workers=2, scratch_root=tmp_path):
            open_batches = [event for event in producer.events if event[0] == "open"]
            assert producer.events[-1] == open_batches[-1], "a batch closed or opened before consumption"
            assert len(list(tmp_path.iterdir())) == 1, "an earlier batch's scratch remains"
            producer.events.append(("consume", tuple(staged.query_columns)))
    assert [event[:3] for event in producer.events if event[0] != "close"] == [
        ("open", 0, ("q1", "q2")), ("consume", ("q1", "q2")),
        ("open", 1, ("q3", "q4")), ("consume", ("q3", "q4")),
        ("open", 2, ("q5",)), ("consume", ("q5",))]
    assert len(producer.pools) == 9 and all(seen is pool for seen in producer.pools)
    assert not list(tmp_path.iterdir())


def test_statuses_are_filtered_per_batch_and_finalized(tmp_path):
    producer = FakeProducer()
    skipped = QueryAnnotationStatus("empty", "empty.txt", "gene_list", "skipped", "empty_gene_list")
    statuses = (ok("q1"), ok("flat"), skipped, ok("q3"))
    with ThreadPoolExecutor(max_workers=1) as pool:
        batches = list(run_query_batches(producer, ["q1", "flat", "q3"], statuses, width=2, pool=pool,
                                         workers=1, scratch_root=tmp_path))
    assert [event[3] for event in producer.events if event[0] == "open"] == [
        ("q1", "flat", "empty"), ("q3",)]
    assert [batch.query_columns for batch in batches] == [["q1"], ["q3"]]
    assert [(s.query, s.status, s.reason) for s in batches[0].query_statuses] == [
        ("q1", "ok", ""), ("flat", "skipped", "zero_variance_ld_scores"), ("empty", "skipped", "empty_gene_list")]
    assert [record["column"] for record in batches[0].count_records] == ["q1"]
    assert batches[0].annotation_types == {"base": "binary", "q1": "binary"}


def test_failure_cancels_pending_tasks_and_submits_no_more(tmp_path):
    submitted, shutdown = [], []
    failed, pending = Future(), Future()
    failed.set_exception(ValueError("task failure"))

    class Pool:
        def submit(self, function, chrom, batch):
            submitted.append(chrom)
            return failed if chrom == "1" else pending

        def shutdown(self, *, wait):
            shutdown.append(wait)

    with pytest.raises(ValueError, match="task failure"):
        list(run_query_batches(FakeProducer(), ["q1"], (ok("q1"),), width=1, pool=Pool(), workers=2,
                               scratch_root=tmp_path))
    assert submitted == ["1", "2"]
    assert pending.cancelled()
    assert shutdown == [True]


def test_failure_waits_for_running_workers_before_removing_scratch(tmp_path):
    running, shutting_down, finished = Event(), Event(), Event()
    order = []

    class Pool(ThreadPoolExecutor):
        def shutdown(self, wait=True, **kwargs):
            shutting_down.set()
            super().shutdown(wait=wait, **kwargs)
            order.append("workers finished")

    def task(chrom, batch):
        if chrom == "1":
            assert running.wait(5)
            raise OSError("injected worker failure")
        running.set()
        assert shutting_down.wait(5)
        assert (batch.scratch / "borrowed.txt").read_text() == "borrowed by tasks"
        finished.set()

    producer = FakeProducer(task, chromosomes=("1", "2"))
    with Pool(max_workers=2) as pool:
        with pytest.raises(OSError, match="injected worker failure"):
            list(run_query_batches(producer, ["q1"], (ok("q1"),), width=1, pool=pool, workers=2,
                                   scratch_root=tmp_path))
        order.append(producer.events[-1][0])
    assert finished.is_set()
    assert order[:2] == ["workers finished", "close"]
    assert not list(tmp_path.iterdir())
