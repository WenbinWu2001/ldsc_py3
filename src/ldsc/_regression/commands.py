"""Command-line and keyword entry points of ``ldsc h2``, ``ldsc partitioned-h2`` and ``ldsc rg``.

Each command has an argument registrar (``add_*_arguments``), a
``run_*_from_args`` function that ``ldsc`` dispatches to, and a Python form
(``run_h2``, ``run_partitioned_h2``, ``run_rg``) that takes the command's
option names as keyword arguments and writes the same result directory.
A run turns its options into a ``RegressionConfig``, preflights its owned
outputs and workflow log under ``--output-dir``, checks its declared input
paths, loads the munged summary statistics and the LD-score directory, and
then calls the shared alignment and fits (h2), the partitioned batch
(``partitioned``) or the pair run (``rg``), and writes the result through
its directory writer. The registered ``GlobalConfig`` is only the SNP
identifier fallback for an input without recorded provenance and the source
of the workflow-log level.
"""

from __future__ import annotations

import logging
import numbers
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from typing import Sequence

from .._runtime.cli_help import CLIHelpFormatter, SCALAR_PATH_HELP, add_output_dir_argument, add_log_level_argument
from .._runtime.workflow_logging import args_output_scope
from ..config import GlobalConfig, RegressionConfig, get_global_config, print_global_config_banner
from .._runtime.paths import ensure_output_directory, normalize_path_token
from .._runtime.workflow_logging import log_inputs, log_outputs, materializing_overwrite_guard, workflow_logging
from . import partitioned, rg
from .._runtime.parallelism import (
    WORKERS_HELP, parse_query_batch_size, parse_workers, validate_query_batch_size, validate_workers,
)
from .._ldscore.source import load_ldscore_from_dir
from .outputs import (
    H2DirectoryWriter,
    H2OutputConfig,
    PartitionedH2DirectoryWriter,
    RgDirectoryWriter,
    RgOutputConfig,
)
from .results import H2Result, PartitionedH2BatchResult, RgResultFamily
from .alignment import _log_regression_identity, _prepare_h2_dataset, _trait_label
from .fitting import _fit_h2_dataset, _h2_metadata, _ldscore_dir_of, _log_phase_timing, summarize_total_h2
from .partitioned import (
    PARTITIONED_H2_SUMMARY_SORT_COLUMNS,
    _log_partitioned_h2_regime,
    _log_quantitative_annotation_interpretation,
)
from .rg import _disambiguate_trait_names
from .._sumstats.table import SumstatsTable, load_sumstats
from ..errors import LDSCUsageError, LDSCUserError
from .._runtime.progress import PhaseProgress
from .prevalence import resolve_rg_prevalences
from .._runtime.preflight import InputGate, inspect_artifact_paths, inspect_declared_inputs
from .._runtime.workflow_kwargs import command_parser, run_from_kwargs


LOGGER = logging.getLogger("LDSC.regression_runner")


def _add_scalar_prevalence_arguments(parser) -> None:
    """Add paired case fractions for optional liability-scale estimates."""
    parser.add_argument(
        '--samp-prev', type=float, default=None, metavar='VALUE',
        help=(
            'Sample case fraction, strictly between 0 and 1, for liability-scale estimates. Requires '
            '--pop-prev; omit both, or pass nan for both, to report observed-scale estimates only.'
        ),
    )
    parser.add_argument(
        '--pop-prev', type=float, default=None, metavar='VALUE',
        help=(
            'Population case fraction, strictly between 0 and 1, for liability-scale estimates. Requires '
            '--samp-prev; omit both, or pass nan for both, to report observed-scale estimates only.'
        ),
    )


def add_h2_arguments(parser) -> None:
    """Register grouped h2 command-line options."""
    parser.prog = "ldsc h2"
    parser.formatter_class = CLIHelpFormatter
    parser.description = "Estimate heritability from munged summary statistics and LD scores."
    sections = _add_regression_sections(parser, include_h2_intercept=True)
    sections.inputs.add_argument(
        '--trait-name', default=None,
        help=(
            'Trait label used in summaries. If omitted, use the stored trait name or derive one from the '
            'input filename.'
        ),
    )
    _add_regression_output_arguments(sections.output)
    _add_scalar_prevalence_arguments(sections.scale)
    _add_regression_model_arguments(sections, include_h2_intercept=True)
    add_log_level_argument(sections.logs)


def add_partitioned_h2_arguments(parser) -> None:
    """Register grouped partitioned-h2 command-line options."""
    parser.prog = "ldsc partitioned-h2"
    parser.formatter_class = CLIHelpFormatter
    parser.description = "Estimate partitioned heritability with independent query models."
    sections = _add_regression_sections(parser, include_h2_intercept=True)
    sections.inputs.add_argument(
        '--trait-name', default=None,
        help=(
            'Trait label used in summaries. If omitted, use the stored trait name or derive one from the '
            'input filename.'
        ),
    )
    _add_regression_output_arguments(sections.output)
    sections.output.add_argument(
        '--continue-on-query-error', action='store_true', default=False,
        help=(
            'Skip and record queries whose model preparation, regression, or result calculation fails; '
            'publish successful fits and write diagnostics/query_status.tsv plus error tracebacks in the log. '
            'Default: off; any query error prevents publication. Shared input/output failures and runs '
            'with no successful queries still fail.'
        ),
    )
    sections.output.add_argument(
        '--summary-sort-by', default='auto', choices=('auto', *PARTITIONED_H2_SUMMARY_SORT_COLUMNS),
        help=(
            'Order partitioned_h2.tsv by the selected column. Default: auto, coefficient-p for per-query '
            'fits and category for baseline-only fits. P-values sort ascending; other numeric values sort '
            'descending.'
        ),
    )
    _add_scalar_prevalence_arguments(sections.scale)
    _add_regression_model_arguments(sections, include_h2_intercept=True)
    sections.computing.add_argument(
        '--workers', type=parse_workers, default=1, metavar='N',
        help='Maximum complete query-plus-baseline models running simultaneously; also capped by query batch size. ' + WORKERS_HELP,
    )
    sections.computing.add_argument(
        '--query-batch-size', type=parse_query_batch_size, default=100, metavar='N',
        help=(
            'Maximum number of query annotations loaded at once; each query is fitted separately with the '
            'shared baseline. Default: 100; a smaller positive value reduces memory and caps active query workers.'
        ),
    )
    add_log_level_argument(sections.logs)


def add_rg_arguments(parser) -> None:
    """Register grouped rg command-line options."""
    parser.prog = "ldsc rg"
    parser.formatter_class = CLIHelpFormatter
    parser.description = "Estimate genetic correlation."
    sections = _add_regression_sections(parser, include_h2_intercept=False)
    sections.inputs.add_argument(
        '--anchor-trait', default=None,
        help=(
            'Trait name or input file path to compare against every other trait. Must identify one of '
            '--sumstats-sources. If omitted, estimate every unique trait pair.'
        ),
    )
    _add_regression_output_arguments(sections.output)
    sections.output.add_argument(
        '--write-per-pair-detail', action='store_true', default=False,
        help=(
            'Also write detailed results for each tested pair under diagnostics/pairs/. Default: off; write '
            'aggregate pair and trait tables only.'
        ),
    )
    sections.scale.add_argument(
        '--samp-prev', default=None,
        help=(
            'Comma-separated sample case fractions in resolved --sumstats-sources order, one per trait; '
            'each in (0, 1) or nan for a quantitative trait. Requires --pop-prev; cannot be combined with '
            '--prevalence-manifest. If neither form is supplied, report observed-scale estimates only.'
        ),
    )
    sections.scale.add_argument(
        '--pop-prev', default=None,
        help=(
            'Comma-separated population case fractions in resolved --sumstats-sources order, one per trait; '
            'each in (0, 1) or nan. Requires --samp-prev; cannot be combined with --prevalence-manifest. If '
            'neither form is supplied, report observed-scale estimates only.'
        ),
    )
    sections.scale.add_argument(
        '--prevalence-manifest', default=None,
        help=(
            'TSV with trait_name, samp_prev, and pop_prev columns, matched to exact input trait names. '
            'Cannot be combined with --samp-prev or --pop-prev. If omitted, use paired prevalence lists '
            'when supplied; otherwise report observed-scale estimates only.'
        ),
    )
    _add_regression_model_arguments(sections, include_h2_intercept=False)
    add_log_level_argument(sections.logs)


@materializing_overwrite_guard(
    args_output_scope,
    command="run_h2_from_args(...)",
)
def run_h2_from_args(args):
    """Run single-trait heritability estimation from parsed CLI arguments.

    The regression settings become a ``RegressionConfig`` before the output
    preflight. The workflow requires ``args.output_dir`` and preflights
    ``h2.tsv``, the LD score regression-bin diagnostic, metadata, and
    ``diagnostics/h2.log`` before loading inputs. A successful overwrite
    removes default plots and liability-scale post-processing derived from the
    superseded h2 result. Returns the :class:`H2Result` of the run.
    """
    global_config, config = _configs_from_args(args, "h2")
    return _run_h2(
        global_config, config, sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir,
        trait_name=args.trait_name, output_dir=args.output_dir, overwrite=args.overwrite, log_level=args.log_level,
    )


def _run_h2(
    global_config: GlobalConfig,
    config: RegressionConfig,
    *,
    sumstats_file: str,
    ldscore_dir: str,
    trait_name: str | None,
    output_dir: str,
    overwrite: bool,
    log_level: str,
) -> H2Result:
    """Estimate SNP heritability from typed settings and write the ``h2`` result directory.

    ``global_config`` is the registered configuration, the SNP identifier
    fallback for an input without recorded provenance; the caller opens the
    failure-marker scope.
    """
    output_dir, log_path = _preflight_regression_outputs(output_dir, overwrite, "h2", H2DirectoryWriter)
    with workflow_logging("h2", log_path, log_level=log_level):
        _inspect_regression_inputs(output_dir, sumstats=sumstats_file, ldscore_dir=ldscore_dir)
        print_global_config_banner("run_h2_from_args", global_config)
        log_inputs(sumstats_file=sumstats_file, ldscore_dir=ldscore_dir, output_dir=output_dir)
        LOGGER.info(f"Starting h2 regression for '{sumstats_file}' using LD-score directory '{ldscore_dir}'.")
        sumstats_table = _load_sumstats_table(sumstats_file, trait_name)
        ldscore_result = load_ldscore_from_dir(ldscore_dir)
        prepared = _prepare_h2_dataset(sumstats_table, ldscore_result, config, global_config)
        dataset, legacy_drops, nonfinite_rows = (prepared.dataset, prepared.legacy_drops,
                                                 prepared.dropped_nonfinite_statistic_rows)
        _log_regression_identity([sumstats_table], ldscore_result, global_config,
                                 [(dataset.effective_snp_identifier, dataset.identity_downgrade_applied)])
        prepared = None  # the fitted dataset below replaces the unfiltered one
        outcome = _fit_h2_dataset(dataset, config)
        hsq, dataset = outcome.estimator, outcome.dataset
        n_snps_used = outcome.n_snps
        summary = summarize_total_h2(
            hsq, dataset, trait_name=_trait_label(sumstats_table), n_snps_used=n_snps_used,
            samp_prev=config.samp_prev, pop_prev=config.pop_prev,
        )
        written = H2DirectoryWriter().write(
            summary,
            H2OutputConfig(
                output_dir=output_dir,
                overwrite=overwrite,
            ),
            metadata=_h2_metadata(
                sumstats_table, _ldscore_dir_of(ldscore_result), outcome, config,
                dropped_nonfinite_statistic_rows=nonfinite_rows,
            ),
            diagnostic_bins=outcome.diagnostic_bins,
            legacy_sumstats_drops=legacy_drops,
        )
        log_outputs(**written)
        LOGGER.info(f"Finished h2 regression with {n_snps_used} regression SNPs.")
    return H2Result(summary, written, legacy_drops)


@materializing_overwrite_guard(
    args_output_scope,
    command="run_partitioned_h2_from_args(...)",
)
def run_partitioned_h2_from_args(args):
    """Run batch partitioned heritability from parsed CLI arguments.

    The regression settings become a ``RegressionConfig`` before the output
    preflight; the argparse types of ``--workers`` and ``--query-batch-size``,
    or :func:`run_partitioned_h2`, have already checked those two options.
    The workflow requires ``args.output_dir`` and preflights
    ``partitioned_h2.tsv``, ``diagnostics/query_annotations/``, and
    ``diagnostics/partitioned-h2.log`` before loading inputs. Query-annotation
    runs always write per-query results; baseline-only runs keep the complete
    fitted-model artifacts at the result root. It returns the
    ``PartitionedH2BatchResult``, whose ``output_paths`` also name the legacy
    summary-statistics audit when one is written. The saved root TSV appends
    ``query_result_dir`` relative to ``args.output_dir`` (``NaN`` for baseline-only
    runs); the returned summary and per-query tables retain scientific columns
    only. A successful overwrite removes
    the default plot root derived from the superseded result.
    ``args.continue_on_query_error`` explicitly permits query regression
    failures: successful models are published and every attempted query is
    recorded in ``diagnostics/query_status.tsv``. Strict mode, the default,
    collects query errors and publishes only the diagnostic ledger on failure.
    ``args.workers`` uses the shared LD-score/index worker-count policy, with
    query count and ``args.query_batch_size`` as the work limits. The default
    and every effective count of one use one worker. Each worker fits
    one complete query model with one native numerical thread; baseline-only
    input remains one joint fit in one worker. Worker death and transport/output
    failures abort under either query-error policy. Requested/effective counts
    and phase timings are recorded in the workflow log.
    """
    global_config, config = _configs_from_args(args, "partitioned-h2")
    return _run_partitioned_h2(
        global_config, config, sumstats_file=args.sumstats_file, ldscore_dir=args.ldscore_dir,
        trait_name=args.trait_name, output_dir=args.output_dir, overwrite=args.overwrite, log_level=args.log_level,
        workers=args.workers, query_batch_size=args.query_batch_size,
        continue_on_query_error=args.continue_on_query_error, summary_sort_by=args.summary_sort_by,
    )


def _run_partitioned_h2(
    global_config: GlobalConfig,
    config: RegressionConfig,
    *,
    sumstats_file: str,
    ldscore_dir: str,
    trait_name: str | None,
    output_dir: str,
    overwrite: bool,
    log_level: str,
    workers: int,
    query_batch_size: int,
    continue_on_query_error: bool,
    summary_sort_by: str,
) -> PartitionedH2BatchResult:
    """Fit partitioned heritability from typed settings and write the ``partitioned-h2`` result directory.

    ``global_config`` is the SNP identifier fallback, as in :func:`_run_h2`;
    the caller opens the failure-marker scope.
    """
    output_dir, log_path = _preflight_regression_outputs(
        output_dir, overwrite, "partitioned-h2", PartitionedH2DirectoryWriter,
    )
    with workflow_logging("partitioned-h2", log_path, log_level=log_level):
        _inspect_regression_inputs(output_dir, sumstats=sumstats_file, ldscore_dir=ldscore_dir)
        print_global_config_banner("run_partitioned_h2_from_args", global_config)
        log_inputs(sumstats_file=sumstats_file, ldscore_dir=ldscore_dir, output_dir=output_dir)
        LOGGER.info(
            f"Starting partitioned-h2 regression for '{sumstats_file}' using LD-score directory '{ldscore_dir}'."
        )
        with _log_phase_timing("sumstats loading"):
            sumstats_table = _load_sumstats_table(sumstats_file, trait_name)
        with _log_phase_timing("LD-score loading"):
            ldscore_result = load_ldscore_from_dir(ldscore_dir)
        has_queries = bool(ldscore_result.query_columns)
        _log_partitioned_h2_regime(ldscore_result, has_queries)
        _log_quantitative_annotation_interpretation(ldscore_result)
        result = partitioned.estimate_partitioned_h2_batch(
            sumstats_table, ldscore_result, config, global_config, output_dir=output_dir,
            query_batch_size=query_batch_size,
            workers=workers,
            continue_on_query_error=continue_on_query_error,
            overwrite=overwrite,
            summary_sort_by=summary_sort_by,
        )
        log_outputs(**result.output_paths)
        LOGGER.info(
            f"Finished partitioned-h2 regression for {len(ldscore_result.query_columns)} query annotations "
            f"and {len(result.summary)} summary rows."
        )
    return result


@materializing_overwrite_guard(
    args_output_scope,
    command="run_rg_from_args(...)",
)
def run_rg_from_args(args):
    """Run multi-trait genetic-correlation estimation from parsed CLI args.

    The regression settings become a ``RegressionConfig`` before the output
    preflight. The workflow requires ``--output-dir`` and writes the rg result
    family: ``rg.tsv``, ``rg_full.tsv``, ``h2_per_trait.tsv``, optional
    ``diagnostics/pairs/``, and workflow-owned ``diagnostics/rg.log``. A
    successful overwrite removes the default plot root derived from the
    superseded result.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed ``rg`` options holding every option of
        :func:`add_rg_arguments`, as ``ldsc rg`` and :func:`run_rg` supply
        them; a missing option raises ``AttributeError``. ``sumstats_sources``,
        ``ldscore_dir`` and ``output_dir`` must be set.

    Returns
    -------
    RgResultFamily
        Complete in-memory result family, with ``output_paths`` naming the
        files written. From the command line, a run whose harmonization
        dropped SNPs also prints one WARNING with the counts to stderr; the
        Python wrapper ``run_rg`` never prints.

    Raises
    ------
    LDSCConfigError
        Before the output preflight, if a regression setting such as
        ``num_blocks`` or ``chisq_max`` is invalid.
    LDSCUsageError
        If ``output_dir`` is unset, or if the prevalence lists, the prevalence
        manifest or ``anchor_trait`` do not match the resolved traits.
    LDSCUserError
        If fewer than two summary-statistics paths resolve. A declared input
        that is missing or unreadable raises ``LDSCInputError``.
    """
    global_config, config = _configs_from_args(args, "rg")
    return _run_rg(
        global_config, config, sumstats_sources=args.sumstats_sources, ldscore_dir=args.ldscore_dir,
        anchor_trait=args.anchor_trait, samp_prev=args.samp_prev, pop_prev=args.pop_prev,
        prevalence_manifest=args.prevalence_manifest, output_dir=args.output_dir, overwrite=args.overwrite,
        write_per_pair_detail=args.write_per_pair_detail, log_level=args.log_level,
    )


def _run_rg(
    global_config: GlobalConfig,
    config: RegressionConfig,
    *,
    sumstats_sources: Sequence[str],
    ldscore_dir: str,
    anchor_trait: str | None,
    samp_prev: str | None,
    pop_prev: str | None,
    prevalence_manifest: str | None,
    output_dir: str,
    overwrite: bool,
    write_per_pair_detail: bool,
    log_level: str,
) -> RgResultFamily:
    """Estimate genetic correlations from typed settings and write the ``rg`` result family.

    ``samp_prev`` and ``pop_prev`` are the comma-separated per-trait lists,
    resolved here against the loaded traits; ``config`` carries no
    prevalences. ``global_config`` is the SNP identifier fallback, as in
    :func:`_run_h2`; the caller opens the failure-marker scope.
    """
    output_dir, log_path = _preflight_regression_outputs(output_dir, overwrite, "rg", RgDirectoryWriter)
    with workflow_logging("rg", log_path, log_level=log_level):
        declared_inputs = _inspect_regression_inputs(
            output_dir, sumstats=sumstats_sources, ldscore_dir=ldscore_dir, prevalence_manifest=prevalence_manifest,
        )
        sumstats_paths = declared_inputs.files['sumstats']
        if len(sumstats_paths) < 2:
            raise LDSCUserError(
                "rg requires at least two summary-statistics inputs. Most likely `--sumstats-sources` "
                "resolved to fewer than two files. Pass two or more munged `.sumstats` artifacts."
            )
        print_global_config_banner("run_rg_from_args", global_config)
        log_inputs(
            sumstats_sources=f"{len(sumstats_paths)} files; source paths are recorded in diagnostics/metadata.json",
            anchor_trait=anchor_trait or "none",
            ldscore_dir=ldscore_dir,
            output_dir=output_dir,
        )
        LOGGER.info(
            f"Starting rg regression for {len(sumstats_paths)} summary-statistics files "
            f"using LD-score directory '{ldscore_dir}'."
        )
        gate = InputGate('Summary-statistic content validation', Path(output_dir)/'diagnostics/input_issues.tsv')
        sumstats_tables = []
        with PhaseProgress(LOGGER, 'validation', 'summary-statistic content', len(sumstats_paths)) as progress:
            for path in sumstats_paths:
                table = gate.check('sumstats', path, lambda path=path: _load_sumstats_table(str(path), None))
                if table is not None:
                    sumstats_tables.append(table)
                progress.advance(object=str(path))
            gate.finish()
        # Resolve prevalence against the original munged names so the manifest
        # duplicate-name guard fires before disambiguation renames collisions.
        prevalences = resolve_rg_prevalences(
            samp_prev=samp_prev,
            pop_prev=pop_prev,
            manifest_path=prevalence_manifest,
            trait_names=[_trait_label(table) for table in sumstats_tables],
            trait_paths=[str(path) for path in sumstats_paths],
        )
        sumstats_tables = _disambiguate_trait_names(sumstats_tables)
        if prevalences is not None:
            LOGGER.info("Liability-scale prevalences applied to %d traits; values are recorded in diagnostics/metadata.json.", len(sumstats_tables))
            if LOGGER.isEnabledFor(logging.DEBUG):
                LOGGER.debug(
                    "Liability-scale prevalences applied: "
                    + ", ".join(
                        f"{_trait_label(table)}=(P={sp}, K={pp})"
                        for table, (sp, pp) in zip(sumstats_tables, prevalences)
                    )
                )
        anchor_index = _resolve_anchor_index(anchor_trait, sumstats_paths, sumstats_tables)
        ldscore_result = load_ldscore_from_dir(ldscore_dir)
        result = rg.estimate_rg_pairs(
            sumstats_tables,
            ldscore_result,
            config,
            global_config,
            anchor_index=anchor_index,
            prevalences=prevalences,
        )
        # The pair run fits in a child; its per-pair metadata records each pair's realized mode.
        resolved = [pair for pair in result.per_pair_metadata if "effective_snp_identifier" in pair]
        _log_regression_identity(
            sumstats_tables, ldscore_result, global_config,
            [(pair["effective_snp_identifier"], pair["identity_downgrade_applied"]) for pair in resolved],
            unresolved_pairs=len(result.per_pair_metadata) - len(resolved),
        )
        written = RgDirectoryWriter().write(
            result,
            RgOutputConfig(
                output_dir=output_dir,
                overwrite=overwrite,
                write_per_pair_detail=write_per_pair_detail,
            ),
        )
        log_outputs(**written)
        LOGGER.info(
            f"Finished rg regression for {len(result.rg)} trait pairs "
            f"and {len(result.h2_per_trait)} traits."
        )
    return replace(result, output_paths=written)


def emit_rg_notices(result: RgResultFamily) -> None:
    """Print one stderr WARNING when harmonization dropped SNPs from any rg pair.

    ``ldsc rg`` calls this after a successful run (``ldsc.cli._COMMANDS``);
    Python callers read the counts from the result. The workflow log already
    holds each pair's line; the notice exists because batch-scheduler users
    rarely read logs after a successful run.
    """
    unknown = sum(result.dropped_orientation_unknown_rows or [])
    incompatible = sum(result.dropped_incompatible_allele_rows or [])
    if not unknown and not incompatible:
        return
    print(
        f"WARNING: rg could not align {unknown + incompatible:,} SNPs between traits and dropped them: "
        f"{unknown:,} A/T or C/G SNPs whose strand, and so the sign of their statistic, is unknown (expected), "
        f"and {incompatible:,} SNPs whose two traits list different allele pairs (check the inputs if nonzero). "
        "Per-pair counts are dropped_orientation_unknown_rows and dropped_incompatible_allele_rows in "
        "diagnostics/metadata.json.",
        file=sys.stderr,
    )


def run_h2(**kwargs) -> H2Result:
    """Estimate SNP heritability and save the ``h2`` result directory.

    This is the Python form of ``ldsc h2`` and writes the same files. Keyword
    arguments are the command's option names without the leading ``--``:
    ``sumstats_file``, ``ldscore_dir`` and ``output_dir`` are required;
    ``trait_name``, ``count_kind``, ``intercept_h2``, ``two_step_cutoff``,
    ``chisq_max``, ``num_blocks``, ``allow_identity_downgrade``, ``samp_prev``,
    ``pop_prev`` and ``overwrite`` are optional and default to the command's
    defaults. Run the call under ``if __name__ == "__main__":`` because the
    fit runs in a spawned worker.

    The SNP identifier mode and genome build come from the munged summary
    statistics and the LD-score directory; the registered ``GlobalConfig``
    supplies only the level of ``diagnostics/h2.log``.

    Returns
    -------
    H2Result
        ``summary``, the one-row table saved as ``h2.tsv``; ``output_paths``,
        the files written; and ``legacy_sumstats_drops`` for a legacy input.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is unknown, a required one
        is missing, a value is outside the command's choices, a flag is not a
        bool, or ``log_level``, ``snp_identifier`` or ``genome_build`` is passed.
    LDSCConfigError
        Before the output preflight, when a regression setting is invalid,
        such as ``num_blocks`` below 2, half a prevalence pair, or a
        prevalence that is not a number in ``(0, 1)``; a keyword call skips
        the command line's float conversion, so a string such as ``"0.5"``
        is rejected too.

    Notes
    -----
    The result directory holds ``h2.tsv``, the regression-bin diagnostic and
    the metadata that ``plot_result()`` and ``convert_h2_scale()`` read.
    Existing outputs stop the run unless ``overwrite=True``; a failed
    authorized overwrite leaves ``RUN_FAILED.txt`` naming this call.
    """
    return run_from_kwargs(
        "run_h2", command_parser(add_h2_arguments), run_h2_from_args, kwargs,
        registered=get_global_config(),
    )


def run_partitioned_h2(**kwargs) -> PartitionedH2BatchResult:
    """Estimate partitioned heritability and save the ``partitioned-h2`` result directory.

    This is the Python form of ``ldsc partitioned-h2`` and writes the same
    files. Keyword arguments are the command's option names without the
    leading ``--``: ``sumstats_file``, ``ldscore_dir`` and ``output_dir`` are
    required; ``trait_name``, ``count_kind``, ``intercept_h2``,
    ``two_step_cutoff``, ``chisq_max``, ``num_blocks``,
    ``allow_identity_downgrade``, ``samp_prev``, ``pop_prev``,
    ``summary_sort_by``, ``workers``, ``query_batch_size``,
    ``continue_on_query_error`` and ``overwrite`` are optional and default to
    the command's defaults. Run the call under ``if __name__ == "__main__":``
    because fits run in spawned workers.

    An LD-score directory with query annotations gives one baseline-plus-query
    model per query annotation; one without them gives a single joint model of
    the baseline annotations. The SNP identifier mode and genome build come
    from the munged summary statistics and the LD-score directory; the
    registered ``GlobalConfig`` supplies only the level of
    ``diagnostics/partitioned-h2.log``.

    Returns
    -------
    PartitionedH2BatchResult
        ``summary`` with the scientific columns of the saved
        ``partitioned_h2.tsv``, ``query_status`` with one row per attempted
        query (empty for a baseline-only run), ``per_query_artifacts`` and
        ``output_paths``.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is unknown, a required one
        is missing, a value is outside the command's choices, a flag is not a
        bool, ``query_batch_size`` is not a positive integer, or
        ``log_level``, ``snp_identifier`` or ``genome_build`` is passed.
    ValueError
        Before anything is written, when ``workers`` is not a nonzero integer
        or exceeds the detected available CPU count.
    LDSCConfigError
        Before the output preflight, when a regression setting is invalid,
        such as ``num_blocks`` below 2, half a prevalence pair, or a
        prevalence that is not a number in ``(0, 1)``; a keyword call skips
        the command line's float conversion, so a string such as ``"0.5"``
        is rejected too.

    Notes
    -----
    Query errors stop publication unless ``continue_on_query_error=True``.
    Existing outputs stop the run unless ``overwrite=True``; a failed
    authorized overwrite leaves ``RUN_FAILED.txt`` naming this call.
    """
    # The argparse types check these two options on the command line, before any run starts. A keyword
    # call is checked here, before run_from_kwargs opens the failure-marker scope, so it writes nothing.
    if "workers" in kwargs:
        validate_workers(kwargs["workers"])
    if "query_batch_size" in kwargs:
        validate_query_batch_size(kwargs["query_batch_size"], error_type=LDSCUsageError)
    return run_from_kwargs(
        "run_partitioned_h2", command_parser(add_partitioned_h2_arguments), run_partitioned_h2_from_args, kwargs,
        registered=get_global_config(),
    )


def run_rg(**kwargs) -> RgResultFamily:
    """Estimate genetic correlation and save the ``rg`` result directory.

    This is the Python form of ``ldsc rg`` and writes the same files. Keyword
    arguments are the command's option names without the leading ``--``:
    ``sumstats_sources`` (two or more munged summary-statistics paths or
    ``*`` patterns), ``ldscore_dir`` and ``output_dir`` are required;
    ``anchor_trait``, ``count_kind``, ``intercept_h2``, ``intercept_gencov``,
    ``two_step_cutoff``, ``chisq_max``, ``num_blocks``,
    ``allow_identity_downgrade``, ``samp_prev``, ``pop_prev``,
    ``prevalence_manifest``, ``write_per_pair_detail`` and ``overwrite`` are
    optional and default to the command's defaults. Run the call under
    ``if __name__ == "__main__":`` because fits run in a spawned worker.

    ``samp_prev`` and ``pop_prev`` take one value per trait in resolved
    ``sumstats_sources`` order, as a list such as ``[0.5, float("nan")]``
    (``nan`` for a quantitative trait) or as the command's comma-separated
    string; ``prevalence_manifest`` matches prevalences by trait name instead.
    Every trait needs ``A1`` (the effect allele, the allele the signed
    statistic refers to) and ``A2``. Traits that share a label are
    relabeled ``label@<path suffix>``. The SNP identifier mode and genome
    build come from the inputs; the registered ``GlobalConfig`` supplies only
    the level of ``diagnostics/rg.log``.

    Returns
    -------
    RgResultFamily
        ``rg``, ``rg_full`` and ``h2_per_trait`` with the contents of the
        saved ``rg.tsv``, ``rg_full.tsv`` and ``h2_per_trait.tsv``, one
        ``per_pair_metadata`` record per pair, and ``output_paths`` naming
        the files written.

    Raises
    ------
    LDSCUsageError
        Before anything is written, when a keyword is unknown, a required one
        is missing, a value is outside the command's choices, a flag is not a
        bool, ``samp_prev`` or ``pop_prev`` is a single number, or
        ``log_level``, ``snp_identifier`` or ``genome_build`` is passed.

    Notes
    -----
    A pair whose fit fails keeps its row with ``status="failed"``, or
    ``status="unestimable"`` for a singular jackknife deletion or a
    heritability estimate that is not positive; the other pairs are still
    fitted. Existing outputs stop the run unless
    ``overwrite=True``; a failed authorized overwrite leaves
    ``RUN_FAILED.txt`` naming this call.
    """
    for name in ("samp_prev", "pop_prev"):
        value = kwargs.get(name)
        if isinstance(value, numbers.Real):
            raise LDSCUsageError(
                f"run_rg() cannot run: {name}={value!r} must give one value per trait, "
                f"for example {name}=[0.5, float('nan')]."
            )
        if value is not None and not isinstance(value, str):
            kwargs[name] = ",".join(str(item) for item in value)
    return run_from_kwargs(
        "run_rg", command_parser(add_rg_arguments), run_rg_from_args, kwargs,
        registered=get_global_config(),
    )


def _add_regression_sections(parser, include_h2_intercept: bool) -> SimpleNamespace:
    """Create the regression help sections in display order and add the shared inputs.

    Callers add the remaining options section by section, in display order,
    so the usage line follows the help sections. Sections that receive no
    options are omitted from the help.
    """
    parser.allow_abbrev = False
    sections = SimpleNamespace(
        inputs=parser.add_argument_group("Step 1: Munged summary statistics and LD scores"),
        output=parser.add_argument_group("Step 2: Output"),
        scale=parser.add_argument_group("Optional: liability scale"),
        model=parser.add_argument_group("Optional: model settings"),
        advanced=parser.add_argument_group("Optional: advanced estimation and SNP matching"),
        computing=parser.add_argument_group("Computing resources"),
        logs=parser.add_argument_group("Logging"),
    )
    if include_h2_intercept:
        sections.inputs.add_argument(
            '--sumstats-file', required=True, metavar='FILE',
            help=(
            'Required munged summary-statistics Parquet or legacy .sumstats(.gz) file. SNP identity must be '
            'compatible with --ldscore-dir. '
            + SCALAR_PATH_HELP
        ),
        )
    else:
        sections.inputs.add_argument(
            '--sumstats-sources', nargs='+', required=True, metavar='SOURCES',
            help=(
            "Required munged summary-statistics Parquet or legacy .sumstats(.gz) files; accepts exact paths or quoted '*' patterns "
            'matching filename text and requires at least two files. SNP identity must be compatible with '
            "--ldscore-dir. '@' is not expanded."
        ),
        )
    sections.inputs.add_argument(
        '--ldscore-dir', required=True, metavar='DIR',
        help=(
            'Required LD-score directory from ldscore or convert-legacy-ldscores, containing '
            'predictors, regression weights, and SNP counts.'
        ),
    )
    return sections


def _add_regression_output_arguments(output) -> None:
    """Add the regression output directory and overwrite policy."""
    add_output_dir_argument(
        output,
        help=(
            'Required destination for regression estimates and diagnostics.'
        ),
    )
    output.add_argument(
        '--overwrite', action='store_true', default=False,
        help=(
            "Replace this command's existing output files and remove obsolete outputs from an earlier run. "
            'Default: off; stop if output files already exist.'
        ),
    )


def _add_regression_model_arguments(sections: SimpleNamespace, include_h2_intercept: bool) -> None:
    """Add SNP-count, intercept, chi-square, jackknife, and SNP-matching settings."""
    sections.model.add_argument(
        '--count-kind', choices=('common', 'all'), default='common',
        help=(
            'Reference SNP counts used to scale regression estimates: common uses common-SNP counts, '
            'all uses all-SNP counts. Default: common. The requested counts must be in the LD-score '
            'directory; the other kind is never substituted.'
        ),
    )
    if include_h2_intercept:
        sections.advanced.add_argument(
            '--intercept-h2', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the heritability intercept to this value; use 1 for the standard fixed intercept. '
                'Cannot be combined with an explicit --two-step-cutoff. If omitted, estimate the intercept.'
            ),
        )
    else:
        sections.advanced.add_argument(
            '--intercept-h2', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the heritability intercept to the same value for both traits in every pair; use 1 for '
                'the standard fixed intercept. Cannot be combined with an explicit --two-step-cutoff. '
                'If omitted, estimate the intercepts.'
            ),
        )
        sections.advanced.add_argument(
            '--intercept-gencov', type=float, default=None, metavar='VALUE',
            help=(
                'Fix the genetic-covariance intercept to the same value for every trait pair; use 0 for '
                'the standard fixed intercept. Cannot be combined with --two-step-cutoff. Single-annotation '
                'fits also require --intercept-h2 because automatic two-step estimation otherwise conflicts '
                'with this fixed value. If omitted, estimate the covariance intercepts.'
            ),
        )
    if include_h2_intercept:
        sections.advanced.add_argument(
            '--two-step-cutoff', type=float, default=None, metavar='VALUE',
            help=(
                'First-step chi-square cutoff, including the boundary. Requires a single-annotation model with '
                'a free intercept; cannot be combined with --intercept-h2. If omitted, use 30 '
                'for such models; otherwise disable two-step estimation.'
            ),
        )
        sections.advanced.add_argument(
            '--chisq-max', type=float, default=None, metavar='VALUE',
            help=(
                'Keep SNPs with chi-square at or below this value. If omitted, use max(0.001 * maximum N, 80) '
                'for multi-annotation models; single-annotation models have no chi-square cap.'
            ),
        )
    else:
        sections.advanced.add_argument(
            '--two-step-cutoff', type=float, default=None, metavar='VALUE',
            help=(
                'First-step chi-square cutoff, including the boundary. Requires a single-annotation model with '
                'free intercepts; cannot be combined with --intercept-h2 or '
                '--intercept-gencov. If omitted, use 30 with a free heritability intercept in single-annotation fits; '
                'otherwise disable it.'
            ),
        )
        sections.advanced.add_argument(
            '--chisq-max', type=float, default=None, metavar='VALUE',
            help=(
                'Keep SNPs satisfying Z1^2 * Z2^2 <= VALUE^2 for each trait pair. If omitted, apply no '
                'chi-square product cap. This differs from the single-trait chi-square filter.'
            ),
        )
    sections.advanced.add_argument(
        '--num-blocks', type=int, default=200, metavar='N',
        help=(
            'Number of SNP blocks for jackknife standard errors; must be at least 2. Default: 200, capped '
            'at the number of retained SNPs.'
        ),
    )
    sections.advanced.add_argument(
        '--allow-identity-downgrade', action='store_true', default=False,
        help=(
            'Allow allele-aware and allele-unaware inputs of the same mode family by matching without '
            'alleles. rsID-family and coordinate-family modes still cannot be mixed. Default: off; require '
            'compatible SNP identifier modes.'
        ),
    )


def _preflight_regression_outputs(
    output_dir: str | None,
    overwrite: bool,
    workflow_name: str,
    writer,
) -> tuple[str, Path]:
    """Preflight regression outputs and return normalized output dir plus log path."""
    if not output_dir:
        raise LDSCUsageError(
            f"{workflow_name} requires `--output-dir`. Most likely a programmatic namespace omitted output_dir. "
            "Pass an explicit result directory."
        )
    output_path = ensure_output_directory(output_dir, label="output directory")
    log_path = output_path / "diagnostics" / f"{workflow_name}.log"
    writer.artifact_family(output_path).preflight(
        overwrite=overwrite,
        additional_paths=[log_path],
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)
    return str(output_path), log_path


def _inspect_regression_inputs(output_dir, *, sumstats, ldscore_dir, prevalence_manifest=None):
    """Check the declared summary statistics, prevalence manifest and LD-score directory before any fit.

    ``sumstats`` is the ``--sumstats-file`` path or the ``--sumstats-sources``
    tokens. Every issue found is recorded in ``diagnostics/input_issues.tsv``
    before the run stops; otherwise returns the resolved inputs.
    """
    return inspect_declared_inputs(
        files=[('sumstats', sumstats), ('prevalence manifest', prevalence_manifest)],
        checks=[('LD scores', ldscore_dir, lambda: inspect_artifact_paths(ldscore_dir))],
        issues_path=Path(output_dir)/'diagnostics/input_issues.tsv')


def _resolve_anchor_index(
    anchor_trait: str | None,
    sumstats_paths: Sequence[str],
    sumstats_tables: Sequence[SumstatsTable],
) -> int | None:
    """Resolve ``--anchor-trait`` by trait label first, then source path."""
    if not anchor_trait:
        return None
    token = normalize_path_token(anchor_trait)
    trait_matches = {idx for idx, table in enumerate(sumstats_tables) if table.trait_name == token}
    if len(trait_matches) == 1:
        return next(iter(trait_matches))
    if len(trait_matches) > 1:
        available = [table.trait_name or Path(path).name for table, path in zip(sumstats_tables, sumstats_paths)]
        raise LDSCUsageError(
            f"rg could not resolve `--anchor-trait={anchor_trait}` because it matched "
            f"{len(trait_matches)} trait names. Most likely multiple inputs have the same trait label. "
            f"Use a unique input path or one of these resolved trait names: {available}."
        )

    path_matches: set[int] = set()
    try:
        anchor_path = Path(token).resolve(strict=False)
        for idx, path in enumerate(sumstats_paths):
            if Path(path).resolve(strict=False) == anchor_path:
                path_matches.add(idx)
    except OSError:
        pass  # The actionable unknown-anchor error below includes valid inputs.
    if len(path_matches) != 1:
        available = [table.trait_name or Path(path).name for table, path in zip(sumstats_tables, sumstats_paths)]
        raise LDSCUsageError(
            f"rg could not resolve `--anchor-trait={anchor_trait}` because it matched "
            f"{len(path_matches)} input paths. Most likely the anchor token is misspelled or not part "
            f"of `--sumstats-sources`. Use one of these trait names or exact input paths: {available}."
        )
    return next(iter(path_matches))


def _configs_from_args(args, command: str) -> tuple[GlobalConfig, RegressionConfig]:
    """Return the ``GlobalConfig`` and ``RegressionConfig`` of parsed ``command`` arguments.

    ``command`` is ``"h2"``, ``"partitioned-h2"`` or ``"rg"``. Regression
    reads the SNP identifier mode and genome build from its saved inputs; the
    registered ``GlobalConfig`` is the fallback for an input without recorded
    provenance. The scalar ``--samp-prev``/``--pop-prev`` of ``h2`` and
    ``partitioned-h2`` go on the config, which validates them. The per-trait
    lists of ``rg`` are resolved against the loaded traits in :func:`_run_rg`,
    so its config has no prevalences.
    """
    is_rg = command == "rg"
    config = RegressionConfig(
        num_blocks=args.num_blocks,
        use_common_counts=(args.count_kind == "common"),
        intercept_h2=args.intercept_h2,
        intercept_gencov=args.intercept_gencov if is_rg else None,
        two_step_cutoff=args.two_step_cutoff,
        chisq_max=args.chisq_max,
        samp_prev=None if is_rg else args.samp_prev,
        pop_prev=None if is_rg else args.pop_prev,
        allow_identity_downgrade=args.allow_identity_downgrade,
    )
    return get_global_config(), config


def _load_sumstats_table(path: str, trait_name: str | None) -> SumstatsTable:
    """Load one munged summary-statistics file through the public workflow helper."""
    return load_sumstats(path, trait_name=trait_name)
