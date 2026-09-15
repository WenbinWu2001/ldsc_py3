"""Bounded whole-query processes and immutable, output-owned regression maps.

The parent owns every file and publishes results. Workers borrow mapped numeric
columns, keep row labels/model metadata once per process, and stage one fit at a
time. Dedicated pipes carry only query descriptors and compact outcomes. Public
process sentinels detect worker death; all children exit before scratch removal.
Statistical preparation and fitting remain in ``regression_runner``.
"""

from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, replace
import logging
import mmap
import multiprocessing as mp
from multiprocessing.connection import wait
import os
from pathlib import Path
import pickle
import signal
import sys
import traceback

import numpy as np
import pandas as pd

from .errors import LDSCInternalError, LDSCUsageError


@dataclass(frozen=True)
class _MappedColumn:
    path: Path
    offset: int
    rows: int
    dtype: np.dtype

    @contextmanager
    def open(self):
        """Borrow only this column's pages, preserving dtype and row order."""
        start = self.offset // mmap.ALLOCATIONGRANULARITY * mmap.ALLOCATIONGRANULARITY
        offset = self.offset - start
        with self.path.open("rb") as stream:
            with mmap.mmap(stream.fileno(), offset + self.rows * self.dtype.itemsize,
                           access=mmap.ACCESS_READ, offset=start) as mapping:
                yield np.ndarray((self.rows,), dtype=self.dtype, buffer=mapping, offset=offset)


@dataclass(frozen=True)
class _MappedFrame:
    columns: tuple[str, ...]
    numeric: dict[str, _MappedColumn]
    objects_path: Path

    @contextmanager
    def open(self):
        """Load labels once and borrow numeric columns without pandas copies."""
        with ExitStack() as stack:
            objects = pd.read_pickle(self.objects_path)
            values = {name: objects[name] for name in objects}
            values.update({name: stack.enter_context(column.open()) for name, column in self.numeric.items()})
            yield pd.DataFrame(values, columns=self.columns, copy=False)


def _stage_frame(frame, directory):
    """Write dtype-preserving columns once; never pickle a numeric matrix."""
    directory = Path(directory)
    directory.mkdir()
    numeric, objects = {}, []
    with ExitStack() as stack:
        streams = {}
        for name in frame:
            dtype = frame[name].dtype
            if not isinstance(dtype, np.dtype) or dtype.hasobject:
                objects.append(name)
                continue
            if dtype not in streams:
                path = directory / f"values-{len(streams)}.bin"
                streams[dtype] = (path, stack.enter_context(path.open("wb")))
            path, stream = streams[dtype]
            numeric[name] = _MappedColumn(path, stream.tell(), len(frame), dtype)
            np.ascontiguousarray(frame[name].to_numpy(copy=False)).tofile(stream)
    objects_path = directory / "labels.pkl"
    frame.loc[:, objects].to_pickle(objects_path)
    return _MappedFrame(tuple(frame.columns), numeric, objects_path)


@dataclass(frozen=True)
class _PickledQueryColumn:
    """A non-NumPy column keeps its dtype and original model-error behavior."""

    path: Path

    @contextmanager
    def open(self):
        yield pd.read_pickle(self.path)


def _stage_query_batch(frame, directory):
    """Map numeric queries and snapshot other columns separately, once each."""
    snapshot = _stage_frame(frame, directory)
    columns = dict(snapshot.numeric)
    for index, name in enumerate(frame):
        if name not in columns:
            path = Path(directory) / f"query-{index}.pkl"
            frame[name].to_pickle(path)
            columns[name] = _PickledQueryColumn(path)
    return columns


class _WorkerLog(logging.Handler):
    """Send bounded per-model log text to the parent's existing handlers."""

    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append((record.name, record.levelno, self.format(record)))


def _set_worker_accelerate_threads():
    """Use Apple's per-thread BLAS/LAPACK control, absent from threadpoolctl.

    Called only in a child before its first model. macOS 15 added the public
    BLASSetThreading API; older systems need the environment limit at launch.
    """
    if sys.platform != "darwin":
        return
    import ctypes

    library = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
    setter = getattr(library, "BLASSetThreading", None)
    if setter is None:
        if os.environ.get("VECLIB_MAXIMUM_THREADS") != "1":
            raise LDSCUsageError(
                "Parallel query workers on macOS before 15 require VECLIB_MAXIMUM_THREADS=1 "
                "before starting Python. Set it when launching ldsc, or use --threads 1."
            )
        return
    setter.argtypes, setter.restype = [ctypes.c_uint], ctypes.c_int
    # BLAS_THREADING_SINGLE_THREADED = 1 in Accelerate's public thread_api.h.
    if setter(1) != 0:
        raise LDSCUsageError("Accelerate could not enable single-threaded BLAS/LAPACK; use --threads 1.")


def _query_worker(connection, snapshot, inputs_path):
    """Own one complete model at a time; let the parent handle interrupts."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    logs = _WorkerLog()
    try:
        from threadpoolctl import threadpool_limits
        import pyarrow as pa
        from .regression_runner import RegressionRunner

        _set_worker_accelerate_threads()
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS"):
            os.environ[key] = "1"
        pa.set_cpu_count(1)
        pa.set_io_thread_count(1)
        with inputs_path.open("rb") as stream:
            inputs, global_config, config, log_level = pickle.load(stream)
        logger = logging.getLogger("LDSC")
        logger.handlers = [logs]
        logger.setLevel(log_level)
        logger.propagate = False
        # Apply after NumPy/SciPy imports: environment settings alone are too late.
        with threadpool_limits(limits=1), snapshot.open() as frame:
            inputs = replace(inputs, dataset=replace(inputs.dataset, merged=frame))
            runner = RegressionRunner(global_config, config)
            while (task := connection.recv()) is not None:
                query, column, directory = task
                logs.records.clear()
                with column.open() as values:
                    result = runner._fit_partitioned_query(
                        inputs, query, pd.DataFrame({query: values}, copy=False), directory, config=config,
                    )
                connection.send((result, logs.records, None))
    except BaseException as error:
        # Output/transport errors and process-control exceptions are run failures.
        detail = "".join(traceback.format_exception(error))
        error.__traceback__ = error.__cause__ = error.__context__ = None
        try:
            connection.send((None, logs.records, (error, detail)))
        except Exception:
            pass  # The parent also watches the process sentinel and pipe EOF.
    finally:
        connection.close()


class _QueryWorkers:
    """One task per spawned process, reused across bounded loading batches."""

    def __init__(self, workspace, inputs, runner, config, workers):
        self.processes = []
        self.connections = []
        self.snapshot = _stage_frame(inputs.dataset.merged, workspace / "shared")
        inputs_path = workspace / "model.pkl"
        with inputs_path.open("wb") as stream:
            pickle.dump((replace(inputs, dataset=replace(inputs.dataset, merged=pd.DataFrame())),
                         runner.global_config, config, logging.getLogger("LDSC").getEffectiveLevel()), stream)
        context = mp.get_context("spawn")
        try:
            for _ in range(workers):
                parent, child = context.Pipe()
                process = context.Process(target=_query_worker, args=(child, self.snapshot, inputs_path))
                try:
                    process.start()
                except BaseException:
                    parent.close()
                    raise
                finally:
                    child.close()
                self.processes.append(process)
                self.connections.append(parent)
        except BaseException:
            self.close(abort=True)
            raise

    def map(self, tasks):
        """Return compact outcomes in input order, replenishing ready workers."""
        pending = iter(enumerate(tasks))
        active, results = {}, {}

        def submit(connection):
            item = next(pending, None)
            if item is not None:
                index, task = item
                active[connection] = (index, task[0])
                connection.send(task)

        try:
            for connection in self.connections:
                submit(connection)
            sentinels = [process.sentinel for process in self.processes]
            while active:
                ready = wait([*active, *sentinels])
                for connection in list(active):
                    if connection not in ready:
                        continue
                    index, query = active[connection]
                    result, records, fatal = connection.recv()
                    for name, level, message in records:
                        logging.getLogger(name).log(level, "Query %r: %s", query, message)
                    if fatal is not None:
                        error, detail = fatal
                        logging.getLogger("LDSC.regression_runner").error("Query %r worker failed:\n%s", query, detail)
                        error.add_note(f"partitioned-h2 query {query!r} failed in a worker.")
                        raise error
                    results[index] = result
                    del active[connection]
                    submit(connection)
                if any(sentinel in ready for sentinel in sentinels):
                    raise LDSCInternalError("A partitioned-h2 worker exited unexpectedly.")
        except (EOFError, BrokenPipeError, ConnectionResetError, LDSCInternalError) as error:
            self.close(abort=True)
            queries = [query for _, query in active.values()]
            raise LDSCInternalError(
                f"partitioned-h2 worker communication failed for in-flight queries {queries!r}: {error} "
                "Check memory limits or retry with fewer --threads. No results published."
            ) from error
        except BaseException as error:
            # Reap borrowers before the caller unwinds its batch-directory owner.
            self.close(abort=True)
            error.add_note(f"partitioned-h2 in-flight queries: {[query for _, query in active.values()]!r}.")
            raise
        return [results[index] for index in range(len(results))]

    def close(self, *, abort=False):
        """Stop and reap every child before the parent releases mapped files."""
        if abort:
            for process in self.processes:
                if process.is_alive():
                    process.terminate()
        else:
            for connection in self.connections:
                try:
                    connection.send(None)
                except (BrokenPipeError, EOFError, OSError):
                    pass
        for process in self.processes:
            process.join(timeout=5)
            if process.is_alive():
                process.kill()
                process.join()
            process.close()
        for connection in self.connections:
            connection.close()
        self.processes.clear()
        self.connections.clear()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close(abort=exc_type is not None)
