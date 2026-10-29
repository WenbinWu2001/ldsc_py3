"""Run-owned LD-score logs and the shared stop-with-diagnostics sequences.

Logs are independent of disposable scientific scratch. Direct and indexed
scoring stop on rejected gene lists and unusable queries through the same
write-diagnostics, log, raise sequences defined here.
"""

from contextlib import contextmanager, nullcontext
import logging
from pathlib import Path

from .._runtime.workflow_logging import workflow_logging
from .._runtime.worker_runtime import log_native_thread_state

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


def fail_gene_list_resolution(batch, output_config):
    """Stop on rejected gene-list identifiers or query names (Gate A).

    Writes the gene-list audit and resolution summary, logs the rejections,
    and raises the consolidated ``LDSCInputError``.
    """
    from .._annotation.query_outcomes import gene_list_gate_a_message, log_gene_list_rejections
    from ..errors import LDSCInputError
    from .outputs import LDScoreDirectoryWriter

    LDScoreDirectoryWriter().write_gene_list_preflight(batch, output_config)
    log_gene_list_rejections(batch)
    raise LDSCInputError(gene_list_gate_a_message(batch))


def fail_unusable_queries(diagnostics, message, output_config):
    """Stop when the queries or the gene control cannot be used.

    ``diagnostics`` (an annotation bundle or ``QueryDiagnostics``) supplies the
    query statuses, gene-list batch, and chromosome scope to write. The written
    paths, gene-list rejections and SNP support, and query statuses go to the
    workflow log before ``LDSCInputError(message)`` is raised.
    """
    from .._runtime.workflow_logging import log_outputs
    from ..errors import LDSCInputError
    from .outputs import LDScoreDirectoryWriter

    log_outputs(**LDScoreDirectoryWriter().write_query_diagnostics(diagnostics, output_config))
    log_query_outcomes(diagnostics.gene_list_batch, diagnostics.query_statuses)
    raise LDSCInputError(message)


def log_query_outcomes(gene_list_batch, statuses) -> None:
    """Log gene-list rejections and SNP support, then query statuses, after their diagnostics are written."""
    from .._annotation.query_outcomes import (
        log_gene_list_rejections, log_gene_list_snp_support, log_query_annotation_statuses,
    )

    if gene_list_batch is not None:
        log_gene_list_rejections(gene_list_batch)
        log_gene_list_snp_support(gene_list_batch)
    log_query_annotation_statuses(statuses)


def reference_snp_funnel_path(output_dir) -> str:
    """Return the published reference-SNP funnel path of an LD-score output directory."""
    return str(Path(output_dir) / "diagnostics" / "reference_snp_funnel.tsv")


def log_reference_snp_funnel_warning(funnel, diagnostics, *, command="ldscore") -> None:
    """Log one run-level WARNING when the ``annotation_snp_grid`` step excluded panel SNPs."""
    from .._refpanel.funnel import reference_snp_funnel_warning

    message = reference_snp_funnel_warning(funnel, command=command, diagnostics=str(diagnostics))
    if message is not None:
        LOGGER.warning(message)


def emit_reference_snp_funnel_notice(funnel, diagnostics, *, command="ldscore") -> None:
    """Print the funnel WARNING to stderr after a successful command-line run.

    The workflow log already holds the same line; the notice exists because
    batch-scheduler users rarely read logs after a successful run.
    """
    import sys

    from .._refpanel.funnel import reference_snp_funnel_warning

    message = reference_snp_funnel_warning(funnel, command=command, diagnostics=str(diagnostics))
    if message is not None:
        print(f"WARNING: {message}", file=sys.stderr)
