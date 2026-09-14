"""Evaluate annotation preparation concurrency without changing package code.

Run in fresh processes. Reference snapshots and full-value comparisons happen
after measurement; all package scratch remains owned by AnnotationWorkspace.
See docs/audits/annotation-memory/preparation-parallelism.md for the decision.
"""

import argparse
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
from contextlib import contextmanager
import gzip
import json
import multiprocessing
import os
from pathlib import Path
import platform
import resource
import shutil
import sqlite3
import threading
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pyarrow as pa

from annotation_memory import tree_rss
from ldsc import _annotation_sources as sources
from ldsc._annotation_identity import DiskIdentityIndex, IdentityDropSpool
from ldsc._annotation_storage import AnnotationWorkspace
from ldsc.ldscore_calculator import _resolve_worker_count


SERIAL_SELECT = sources._select_rows
SERIAL_WRITE = sources._write_values


class ReadOnlyIdentity(DiskIdentityIndex):
    """Benchmark adapter reusing the real selection policy on a committed DB."""

    def __init__(self, path, mode):
        self.mode = mode
        self._pending_rows = 0
        self.connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
        self.connection.execute("PRAGMA cache_size=-8192")
        self.connection.execute("PRAGMA temp_store=MEMORY")
        self.connection.execute("PRAGMA query_only=ON")


def initialize_worker():
    pa.set_cpu_count(1)
    pa.set_io_thread_count(1)


def finish_group(task):
    ordinal, group, chunk_rows, root, database, identity_mode, mode, chrom, columns, fail_group, events = task
    owner = SimpleNamespace(path=root / f"group-{ordinal}")
    owner.path.mkdir()
    started = time.perf_counter()
    if events is not None:
        (events / f"started-{ordinal}").touch()
    try:
        drops = IdentityDropSpool(owner.path / "drops")
        with ReadOnlyIdentity(database, identity_mode) as identity:
            selections, counts = SERIAL_SELECT([group], [chunk_rows], owner, identity, drops, mode, chrom)
        selected = time.perf_counter()
        if ordinal == fail_group:
            raise OSError("injected chromosome worker failure")
        shards = SERIAL_WRITE([group], selections, counts, owner, columns)
        return {"ordinal": ordinal, "shards": shards, "drops": drops,
                "selection_seconds": selected - started,
                "numeric_seconds": time.perf_counter() - selected}
    finally:
        if events is not None:
            (events / f"finished-{ordinal}").write_text(str(database.exists()))


def run_groups(tasks, workers, executor):
    options = {"max_workers": workers}
    factory = ThreadPoolExecutor
    if executor == "process":
        factory = ProcessPoolExecutor
        options.update(mp_context=multiprocessing.get_context("spawn"), initializer=initialize_worker)
    remaining = iter(tasks)
    with factory(**options) as pool:
        pending = {pool.submit(finish_group, next(remaining)) for _ in range(workers)}
        completed = []
        try:
            while pending:
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                # Check the entire completed set before scheduling another task.
                completed.extend(future.result() for future in done)
                for _ in done:
                    task = next(remaining, None)
                    if task is not None:
                        pending.add(pool.submit(finish_group, task))
        except BaseException:
            for future in pending:
                future.cancel()
            raise
    return sorted(completed, key=lambda record: record["ordinal"])


@contextmanager
def parallel_prototype(workers, executor, fail_group=None, events=None):
    """Replace only coordinator calls during this benchmark; restore afterward."""
    report = {"effective_workers": 1, "group_phases": [], "parallel_used": False}
    ready = None

    def select(groups, group_rows, workspace, identity, drops, mode, chrom):
        nonlocal ready
        effective = _resolve_worker_count(workers, len(groups))
        if effective == 1 or len(groups) == 1:
            return SERIAL_SELECT(groups, group_rows, workspace, identity, drops, mode, chrom)
        identity._finish_transaction()
        assert not identity.connection.in_transaction
        columns = tuple(c for source in groups[0] for c in source.columns)
        tasks = [(i, group, rows, workspace.path, workspace.path / "identity.sqlite", identity.mode,
                  mode, chrom, columns, fail_group, events) for i, (group, rows) in enumerate(zip(groups, group_rows))]
        completed = run_groups(tasks, effective, executor)
        ready = {}
        for record in completed:
            ready.update(record["shards"])
            for frame in record["drops"].frames():
                drops.append(frame)
            if record["drops"].path.exists():
                shutil.rmtree(record["drops"].path)
        report.update(effective_workers=effective, parallel_used=True,
                      group_phases=[{k: r[k] for k in ("ordinal", "selection_seconds", "numeric_seconds")} for r in completed])
        return [], {}

    def write(*args, **kwargs):
        return ready if ready is not None else SERIAL_WRITE(*args, **kwargs)

    original_select, original_write = sources._select_rows, sources._write_values
    sources._select_rows, sources._write_values = select, write
    try:
        yield report
    finally:
        sources._select_rows, sources._write_values = original_select, original_write


@contextmanager
def phase_times():
    totals, counts, input_rows, saved = Counter(), Counter(), {}, []

    def wrap(owner, name):
        original = getattr(owner, name)
        saved.append((owner, name, original))

        def measured(*args, **kwargs):
            if threading.current_thread() is not threading.main_thread():
                return original(*args, **kwargs)
            started = time.perf_counter()
            try:
                result = original(*args, **kwargs)
                if name == "_scan":
                    input_rows[str(result.path)] = result.spool.n_rows
                return result
            finally:
                totals[name] += time.perf_counter() - started
                counts[name] += 1

        setattr(owner, name, measured)

    for name in ("_scan", "_select_rows", "_write_values"):
        wrap(sources, name)
    for name in ("add", "select"):
        wrap(DiskIdentityIndex, name)
    try:
        yield totals, counts, input_rows
    finally:
        for owner, name, original in saved:
            setattr(owner, name, original)


class Monitor:
    """Sample simultaneous process-tree RSS and owned scratch every 50 ms."""

    def __init__(self, root):
        self.root = root
        self.stop = threading.Event()
        self.peak = {"rss_bytes": 0, "scratch_bytes": 0, "scratch_files": 0, "processes": 0}
        self.error = None
        self.thread = threading.Thread(target=self.sample, daemon=True)

    def sample(self):
        try:
            while not self.stop.is_set():
                rss, processes = tree_rss(os.getpid())
                count = size = 0
                for directory, _, names in os.walk(self.root):
                    for name in names:
                        try:
                            size += os.stat(os.path.join(directory, name)).st_size
                            count += 1
                        except FileNotFoundError:
                            pass
                for key, value in (("rss_bytes", rss), ("scratch_bytes", size),
                                   ("scratch_files", count), ("processes", processes)):
                    self.peak[key] = max(self.peak[key], value)
                self.stop.wait(.05)
        except Exception as error:
            self.error = error

    def close(self):
        self.stop.set()
        self.thread.join()
        if self.error is not None:
            raise self.error


def snapshot(prepared, destination):
    destination.mkdir()
    record = {"baseline_columns": prepared.baseline_columns,
              "query_columns": prepared.query_columns,
              "scope": prepared.scope_chromosomes, "shards": {}}
    for chrom, shard in prepared.shards.items():
        root = destination / chrom
        root.mkdir()
        shutil.copyfile(shard.metadata_path, root / "metadata.parquet")
        shutil.copyfile(shard.stores[0].path, root / "values.npy")
        record["shards"][chrom] = shard.n_rows
    prepared.drops.write_tsv(destination / "drops.tsv.gz")
    (destination / "reference.json").write_text(json.dumps(record, indent=2) + "\n")


def compare(prepared, reference, output):
    expected = json.loads((reference / "reference.json").read_text())
    assert list(prepared.baseline_columns) == expected["baseline_columns"]
    assert list(prepared.query_columns) == expected["query_columns"]
    assert list(prepared.scope_chromosomes) == expected["scope"]
    assert list(prepared.shards) == list(expected["shards"])
    for chrom, shard in prepared.shards.items():
        assert shard.n_rows == expected["shards"][chrom]
        pd.testing.assert_frame_equal(shard.metadata(), pd.read_parquet(reference / chrom / "metadata.parquet"))
        values = np.load(reference / chrom / "values.npy", mmap_mode="r")
        for start in range(0, shard.n_rows, 8192):
            actual = shard.read(rows=slice(start, start + 8192))
            np.testing.assert_array_equal(actual, values[start:start + 8192])
            assert actual.dtype == values.dtype == np.float32
        del values
    prepared.drops.write_tsv(output / "drops.tsv.gz")
    with gzip.open(reference / "drops.tsv.gz", "rb") as left, gzip.open(output / "drops.tsv.gz", "rb") as right:
        while True:
            block = left.read(1024 * 1024)
            assert block == right.read(1024 * 1024)
            if not block:
                break
    (output / "drops.tsv.gz").unlink()


def measure(args):
    initialize_worker()
    args.output.mkdir(parents=True, exist_ok=False)
    paths = [Path(p) for p in args.baseline]
    queries = [Path(p) for p in args.query]
    with AnnotationWorkspace(args.output) as workspace, parallel_prototype(args.workers, args.executor) as prototype, phase_times() as timing:
        monitor = Monitor(workspace.path)
        monitor.thread.start()
        started = time.perf_counter()
        try:
            prepared = sources.prepare_annotation_sources(workspace, paths, queries, mode=args.identity,
                                                          chunk_rows=args.chunk_rows, chrom=args.chrom)
            elapsed = time.perf_counter() - started
            parent_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        finally:
            monitor.close()
        totals, counts, input_rows = timing
        result = {"variant": "serial" if args.workers == 1 else args.executor, "workers": args.workers,
                  "prototype": prototype, "elapsed_seconds": elapsed,
                  "phase_seconds": dict(totals), "phase_calls": dict(counts), "peak_sampled": monitor.peak,
                  "parent_peak_rss_bytes": parent_rss if platform.system() == "Darwin" else parent_rss * 1024,
                  "input_rows": input_rows, "retained_rows": {c: s.n_rows for c, s in prepared.shards.items()},
                  "baseline_columns": len(prepared.baseline_columns), "query_columns": len(prepared.query_columns),
                  "identity_drops": prepared.drops.n_rows, "input_compressed_bytes": sum(p.stat().st_size for p in paths + queries),
                  "platform": platform.platform(), "versions": {m.__name__: m.__version__ for m in (np, pd, pa)},
                  "identity": args.identity, "chrom": args.chrom}
        if args.reference:
            compare(prepared, args.reference, args.output)
            result["verification"] = "all metadata, values, dtypes, order, and diagnostics equal the serial reference"
        if args.save_reference:
            snapshot(prepared, args.output / "reference")
    assert not list(args.output.glob(".ldsc-annotation-*"))
    result["owned_scratch_cleanup"] = "passed"
    (args.output / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", nargs="+", required=True)
    parser.add_argument("--query", nargs="*", default=[])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--save-reference", action="store_true")
    parser.add_argument("--identity", default="rsid")
    parser.add_argument("--chrom")
    parser.add_argument("--chunk-rows", type=int)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--executor", choices=["thread", "process"], default="thread")
    measure(parser.parse_args())


if __name__ == "__main__":
    main()
