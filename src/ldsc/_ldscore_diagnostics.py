"""Run-owned LD-score logs, independent of disposable scientific scratch."""

from contextlib import contextmanager, nullcontext
import logging
from pathlib import Path

from ._logging import workflow_logging
from ._worker_runtime import log_native_thread_state

LOGGER = logging.getLogger("LDSC.ldscore")


@contextmanager
def ldscore_logs(output_dir, chromosomes, *, overwrite=False, log_level="INFO"):
    """Preflight/reset owned logs once per run; chromosome tasks append later."""
    from .outputs import LDScoreDirectoryWriter
    if output_dir is None:
        yield {}
        return
    root = Path(output_dir)
    main = root / "diagnostics" / "ldscore.log"
    active = any(isinstance(handler, logging.FileHandler) and Path(handler.baseFilename).resolve() == main.resolve()
                 for handler in logging.getLogger("LDSC").handlers)
    if not active:
        LDScoreDirectoryWriter.artifact_family(root, include_logs=True).preflight(overwrite=overwrite)
    directory = main.parent / "chromosome-logs"
    directory.mkdir(parents=True, exist_ok=True)
    if overwrite:
        for path in directory.glob("chr*.log"):
            path.unlink()
    paths = {chrom: directory / f"chr{chrom}.log" for chrom in chromosomes}
    for path in paths.values():
        path.touch(exist_ok=False)
    with (nullcontext() if active else workflow_logging("ldscore", main, log_level=log_level)):
        for chrom, path in paths.items():
            LOGGER.info("Chromosome %s diagnostic log: %s", chrom, path)
        yield paths


@contextmanager
def chromosome_log(path, chrom, batch):
    """Attach one child-local handler and retain full task exception details."""
    logger = logging.getLogger("LDSC")
    handler = None
    if path is not None:
        handler = logging.FileHandler(path, mode="a", encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
        logger.addHandler(handler)
    try:
        LOGGER.info("Starting chromosome %s; %s.", chrom, batch)
        log_native_thread_state(LOGGER)
        yield
        LOGGER.info("Finished chromosome %s; %s.", chrom, batch)
    except Exception:
        LOGGER.exception("Chromosome %s failed; %s.", chrom, batch)
        raise
    finally:
        if handler is not None:
            logger.removeHandler(handler)
            handler.close()


def report_chromosome_failure(error, pending, paths, batch, logger=LOGGER):
    """Keep a broken pool's affected tasks distinct from an identified failure."""
    from concurrent.futures.process import BrokenProcessPool
    tasks = "; ".join(f"chromosome {chrom}: {paths.get(chrom, 'no file output')}" for chrom in pending)
    if isinstance(error, BrokenProcessPool):
        logger.error("LD-score worker terminated abnormally; %s. Affected in-flight tasks: %s. "
                     "The chromosome that caused termination is unconfirmed. Error: %s", batch, tasks, error)
    else:
        logger.error("LD-score task failed; %s; %s; %s: %s. %s", batch, tasks, type(error).__name__, error,
                     " ".join(getattr(error, "__notes__", [])))
