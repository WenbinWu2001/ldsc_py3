"""A controlled child for numerical entry points without a task-worker option.

Calls made within an existing compute child execute there directly. This keeps
standalone fitting and preparation under the same native policy without nested
pools or changes to the caller's numerical-library configuration.
"""

from concurrent.futures import ProcessPoolExecutor
from functools import wraps
import importlib
import logging
import multiprocessing as mp
import traceback
import warnings

from . import worker_runtime as _worker_runtime
from .workflow_logging import last_workflow_log_path, set_workflow_log_path


class _LogRecords(logging.Handler):
    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append((record.name, record.levelno, self.format(record)))


def _invoke(module, name, args, kwargs, log_level):
    logger = logging.getLogger("LDSC")
    handler = _LogRecords()
    logger.handlers = [handler]
    logger.propagate = False
    _worker_runtime.initialize_worker(log_level)
    _worker_runtime.log_native_thread_state(logger)
    function = importlib.import_module(module)
    for part in name.split("."):
        function = getattr(function, part)
    with warnings.catch_warnings(record=True) as caught:
        try:
            result = function(*args, **kwargs)
            error = None
        except Exception as exc:
            # A note keeps the child traceback for the parent's workflow-log
            # footer without adding an ERROR record, which the CLI would echo.
            exc.add_note("Calculation child traceback:\n" + traceback.format_exc().rstrip())
            exc.__traceback__ = exc.__context__ = exc.__cause__ = None
            result, error = None, exc
    notices = [(str(w.message), w.category, w.filename, w.lineno) for w in caught]
    return result, handler.records, notices, error, last_workflow_log_path()


def serial_calculation(function):
    """Run a complete serial calculation in one child, or reuse its caller child.

    A workflow log the child opens becomes the caller's most recent log, so a
    failure marker or error message written in the caller names it.
    """
    @wraps(function)
    def controlled(*args, **kwargs):
        if _worker_runtime._INITIALIZED:
            return function(*args, **kwargs)
        with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn")) as pool:
            result, records, notices, error, log_path = pool.submit(
                _invoke, function.__module__, function.__qualname__, args, kwargs,
                logging.getLogger("LDSC").getEffectiveLevel()).result()
        for name, level, message in records:
            logging.getLogger(name).log(level, message)
        for message, category, filename, lineno in notices:
            warnings.warn_explicit(message, category, filename, lineno)
        if log_path is not None:
            set_workflow_log_path(log_path)
        if error is not None:
            raise error
        return result
    return controlled
