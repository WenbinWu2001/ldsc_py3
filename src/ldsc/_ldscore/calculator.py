"""The public LD-score calculation service ``LDScoreCalculator``.

``LDScoreCalculator.run`` computes LD scores from prepared annotations and a
resolved reference panel; the ``ldscore`` command and ``run_ldscore`` in
``ldsc._ldscore.workflow`` resolve user inputs first and call the same
calculation. The result records it returns live in
``ldsc._ldscore.results``, and the direct LD-score producer that runs each
chromosome lives in ``ldsc._ldscore.direct_batches``. Each chromosome task asks
``ref_panel.prepare_chromosome`` for aligned reference data, annotations,
window bounds, and an owned reader. The numerical kernel consumes that
prepared state without reopening the panel.

Writing workflows create ``diagnostics/ldscore.log`` and detailed
``diagnostics/chromosome-logs/chrN.log`` files after output preflight.
``LDScoreCalculator.run(..., output_config=None)`` creates no files.
Every chromosome executes in a worker with one BLAS thread, including
serial execution; workers are reused across query batches.
"""

from __future__ import annotations

from contextlib import ExitStack
from copy import copy
from tempfile import TemporaryDirectory
from dataclasses import dataclass, replace as dataclass_replace
import logging
from typing import Any

import pandas as pd

from .diagnostics import ldscore_logs, log_reference_snp_funnel_warning
from .batch_execution import finalize_batch, run_chromosome_batch, run_query_batches
from .direct_batches import _DirectProducer, _init_worker, _warn_ld_window_snp_subset
from .._runtime.parallelism import chromosome_pool, resolve_worker_count
from .._runtime.progress import report_phase
from .._annotation.saved_results import AnnotationResult
from .._core.snp_identity import identity_mode_family
from ..config import GlobalConfig, LDScoreConfig
from .outputs import LDScoreDirectoryWriter, LDScoreOutputConfig
from ..errors import LDSCUsageError
from .._genome.build_policy import GenomeBuildCheck


LOGGER = logging.getLogger("LDSC.ldscore_calculator")


def _log_annotation_classification(result):
    """Report whole-run type counts; full fitted column identities stay in metadata."""
    counts = pd.Series(list(result.annotation_types.values()), dtype="string").value_counts()
    LOGGER.info("Annotation classification: binary=%d; quantitative=%d. Classification does not change fitted coefficients; see annotation_types in LD-score metadata.",
                counts.get("binary", 0), counts.get("quantitative", 0))
    if LOGGER.isEnabledFor(logging.DEBUG):
        LOGGER.debug("Annotation types: %s", result.annotation_types)


@dataclass(frozen=True)
class _RunContext:
    """What a direct LD-score run settled before ``LDScoreCalculator`` computes it.

    The command-line workflow settles everything first and says so;
    ``LDScoreCalculator.run`` settles nothing, so the calculation does each
    remaining step once and repeats none.

    Attributes
    ----------
    global_config : GlobalConfig
        The run's SNP identity, genome build and log level.
    ref_panel : RefPanel
        The reference-panel adapter, built with ``global_config``.
    main_log_open : bool
        The outputs were preflighted and the workflow log
        ``diagnostics/ldscore.log`` is open; otherwise the calculation does
        both itself.
    build_resolved : bool
        The genome build in ``global_config`` and ``ref_panel`` was resolved
        against every input.
    reference_checked : bool
        ``ref_panel.check_reference`` already checked every chromosome's
        reference files.
    """

    global_config: GlobalConfig
    ref_panel: Any
    main_log_open: bool = False
    build_resolved: bool = False
    reference_checked: bool = False


class LDScoreCalculator:
    """Orchestrate chromosome-wise LD-score calculation.

    This service assembles annotation and reference-panel inputs, delegates the
    heavy computation to the internal LD-score kernel, aggregates chromosome
    outputs, and optionally hands the result to the output layer. The calculator
    delegates reference filtering, annotation alignment, LD windows, and reader
    policy to ``ref_panel.prepare_chromosome``. Reference-panel CM and MAF are
    authoritative. The prepared reader is closed at the chromosome boundary,
    including when numerical computation fails.
    """

    def run(
        self, annotation_bundle, ref_panel, ldscore_config: LDScoreConfig,
        global_config: GlobalConfig, output_config: LDScoreOutputConfig | None = None,
        regression_snps=None, regression_regions=None,
    ):
        """Calculate sequential query batches over one fixed set of reference SNPs.

        Parameters
        ----------
        annotation_bundle : AnnotationBundle or AnnotationResult
            Prepared baseline and query inputs. ``AnnotationBundle.from_frames``
            supplies a small in-memory dataset without filesystem staging. An
            annotation result (``AnnotationResult``) from ``run_annotate``
            needs ``output_config``; it is prepared inside the output
            directory and the preparation is removed when the run ends.
        ref_panel : RefPanel
            Resolved PLINK or R2 matrix adapter; input files are read without mutation.
        ldscore_config : LDScoreConfig
            Scientific settings, query execution width, and chromosome worker
            count. Chromosome workers run within each active query batch.
        global_config : GlobalConfig
            Shared identifier, genome-build, and logging configuration.
        output_config : LDScoreOutputConfig or None, optional
            Required for multiple query batches. If supplied, all batch files
            stay private until success. If omitted, one materialized batch and
            complete diagnostics are returned without filesystem writes.
        regression_snps : set of str or RestrictionIdentityKeys or None, optional
            Already-resolved regression SNP candidates. ``None`` selects all
            retained reference rows before region exclusions. This low-level
            method does not load ``ldscore_config.regr_snps_file`` or apply the
            public workflow's packaged HapMap3 default.
        regression_regions : RegionIntervals or None, optional
            Already-resolved region exclusions for the regression SNPs.
            ``None`` applies no region exclusions. Neither this nor
            ``regression_snps`` changes the reference SNPs used for
            annotation counts and LD scores.

        Returns
        -------
        LDScoreSource or LDScoreResult
            Saved-artifact handle when writing; otherwise an in-memory result.
            Explicit ``read_queries(names)`` on a saved handle has no width cap.

        Raises
        ------
        LDSCUsageError
            Output is omitted for multiple batches or an annotation result
            (``AnnotationResult``), or reference-metadata export is requested
            without an output directory.
        LDSCInputError
            An R2 matrix has an invalid chromosome footer or sidecar.
            Every chromosome is checked before any is computed, and one error
            lists all issues. Also raised, naming every such chromosome, when
            no reference-panel SNP left on a chromosome has an annotation row.

        Notes
        -----
        The annotation rows define the candidate SNPs: reference-panel SNPs
        without an annotation row are excluded, counted in the result's
        ``reference_snp_funnel``, and reported once as a WARNING in the log.
        A chromosome emptied by an explicit restriction or by genotype QC is
        skipped with a WARNING.
        Each direct batch repeats reference traversal. Writing workers stage
        their post-kernel float32 query scores in private chromosome fragments;
        the parent combines scientific summaries and streams those fragments
        without constructing a genome-wide query matrix. Query workspace scales
        with active batch width and chromosome workers. Counts and LD
        accumulation use float64; canonical output retains its usual precision.
        The calculator borrows inputs and never closes the caller's bundle.
        Each worker uses 1 BLAS thread, including
        serial execution. Writing runs retain chromosome logs under
        ``diagnostics/chromosome-logs``; no-output calls create no files.
        """
        return self._run(annotation_bundle, ldscore_config, output_config, regression_snps, regression_regions,
                         context=_RunContext(global_config, ref_panel))

    @report_phase(LOGGER, "computation", "LD-score calculation")
    def _run(self, annotation_bundle, ldscore_config, output_config=None, regression_snps=None,
             regression_regions=None, *, context, skipped_funnels=(), deferred_queries=None):
        """Implement :meth:`run`; ``context`` (``_RunContext``) says what the caller already settled.

        An annotation result (``AnnotationResult``) is prepared inside this
        run's output directory, so the annotate output directory, an input
        here, is never written; the preparation is removed when the run ends.
        ``skipped_funnels`` are the reference-SNP funnels of chromosomes the
        synthetic ``base`` preparation already skipped. ``deferred_queries``
        (``DeferredQueries``), which only the ``ldscore`` workflow passes with
        an output directory, projects the bundle's BED or gene-list queries
        one query batch at a time.
        """
        arguments = (ldscore_config, output_config, regression_snps, regression_regions)
        options = dict(context=context, skipped_funnels=skipped_funnels, deferred_queries=deferred_queries)
        if not isinstance(annotation_bundle, AnnotationResult):
            return self._run_prepared(annotation_bundle, *arguments, **options)
        if output_config is None:
            raise LDSCUsageError(
                "LDScoreCalculator.run() without output_config cannot read an annotation result (AnnotationResult) "
                "from run_annotate directly: it reads no files until you prepare it. Call result.prepared_in(directory), "
                "which prepares the exported annotation files in a private folder of a directory you choose, and pass "
                "the returned AnnotationBundle; close that bundle when you are done. Alternatively pass output_config, "
                "and run() prepares it inside the output directory.")
        with annotation_bundle.prepared_in(output_config.output_dir) as bundle:
            return self._run_prepared(bundle, *arguments, **options)

    def _run_prepared(self, annotation_bundle, ldscore_config, output_config, regression_snps, regression_regions,
                      *, context, skipped_funnels, deferred_queries):
        """Run the LD-score calculation for an annotation bundle that needs no further preparation.

        Resolves the genome build, checks the R2 matrix and opens the workflow
        log only where ``context`` says the caller has not.
        """
        ref_panel = context.ref_panel
        if ldscore_config.export_ref_metadata and output_config is None:
            raise LDSCUsageError(
                "LDScoreCalculator.run() requires output_config when export_ref_metadata=True. "
                "Reference metadata is a filesystem artifact, so pass LDScoreOutputConfig(output_dir=...)."
            )
        if ldscore_config.export_ref_metadata and ref_panel.reference_metadata_sources(ldscore_config)[0] == "parquet_sidecar":
            raise LDSCUsageError(
                "ldscore cannot export reference metadata from an R2 matrix: --export-ref-metadata "
                "(LDScoreConfig.export_ref_metadata) applies to a PLINK reference panel (--plink-prefix). "
                "The R2 matrix sidecars chr<N>_meta.tsv.gz already hold the same CHR, POS, SNP, A1, A2, CM "
                "and MAF metadata; use them directly and drop the option."
            )
        global_config = context.global_config
        if not context.build_resolved:
            global_config, ref_panel = _resolve_calculator_build(annotation_bundle, ref_panel, global_config)
        chromosomes = list(annotation_bundle.chromosomes)
        if not context.reference_checked:
            ref_panel.check_reference(chromosomes, ldscore_config)
        with ExitStack() as stack:
            output_dir = None if output_config is None else output_config.output_dir
            log_paths = stack.enter_context(ldscore_logs(
                output_dir, chromosomes, overwrite=bool(output_config and output_config.overwrite),
                log_level=global_config.log_level, main_log_open=context.main_log_open))
            worker_count = resolve_worker_count(ldscore_config.workers, len(chromosomes))
            export_root = (stack.enter_context(TemporaryDirectory(prefix=".ldsc-ref-", dir=output_dir))
                           if ldscore_config.export_ref_metadata else None)
            pool = stack.enter_context(chromosome_pool(
                worker_count, operation="ldscore chromosome computation", initializer=_init_worker,
                initargs=(regression_snps, regression_regions, global_config.log_level)))
            LOGGER.info(_format_ldscore_start_message(annotation_bundle, len(annotation_bundle.chromosomes)))
            cm_source, maf_source = ref_panel.reference_metadata_sources(ldscore_config)
            LOGGER.info("Reference-panel CM source: %s; MAF source: %s (annotation CM/MAF are ignored).",
                        cm_source, maf_source)
            producer = _DirectProducer(annotation_bundle, ref_panel, ldscore_config, global_config, regression_regions,
                                       log_paths, output_dir=output_dir, export_dir=export_root,
                                       skipped_funnels=skipped_funnels, deferred_queries=deferred_queries)
            if output_config is None:
                if len(annotation_bundle.query_columns) > ldscore_config.query_batch_size:
                    raise LDSCUsageError("Multiple query execution batches require output_config with an output directory. Use one batch for a small in-memory calculation.")
                batch = producer.prepare(dataclass_replace(annotation_bundle, diagnostics_in_memory=True))
                outcomes = run_chromosome_batch(producer, pool, worker_count, batch, "batch 1/1")
                result, _ = producer.reduce(outcomes, batch)
                result = finalize_batch(result, tuple(batch.bundle.query_statuses), ())
                _warn_ld_window_snp_subset(result.ld_window_snp_counts, "LDScoreResult.ld_window_snp_counts")
                log_reference_snp_funnel_warning(result.reference_snp_funnel, "LDScoreResult.reference_snp_funnel")
                _log_annotation_classification(result)
                return result

            results = run_query_batches(
                producer, annotation_bundle.query_columns, tuple(annotation_bundle.query_statuses),
                width=ldscore_config.query_batch_size, pool=pool, workers=worker_count,
                scratch_root=output_config.output_dir)
            result = LDScoreDirectoryWriter().write_batches(results, output_config,
                query_status_order=[status.query for status in annotation_bundle.query_statuses],
                ref_metadata_dir=export_root)
            _log_annotation_classification(result)
            return result


def _format_ldscore_start_message(annotation_bundle, n_chromosomes: int) -> str:
    """Return the workflow-level LD-score start message for a bundle."""
    source_summary = getattr(annotation_bundle, "source_summary", {}) or {}
    baseline_source = str(source_summary.get("baseline", ""))
    if (
        list(getattr(annotation_bundle, "baseline_columns", [])) == ["base"]
        and not list(getattr(annotation_bundle, "query_columns", []))
        and baseline_source.startswith("synthetic all-ones base annotation")
    ):
        return (
            f"Starting LD-score calculation for {n_chromosomes} chromosomes "
            "with synthetic base annotation and no query annotations."
        )
    return (
        f"Starting LD-score calculation for {n_chromosomes} chromosomes "
        f"with {len(annotation_bundle.baseline_columns)} baseline columns "
        f"and {len(annotation_bundle.query_columns)} query columns."
    )


def _resolve_calculator_build(annotation_bundle, ref_panel, global_config: GlobalConfig):
    """Resolve the build of a prepared-input calculation against its annotations and reference panel.

    ``LDScoreCalculator.run`` callers supply inputs the workflow has not
    checked. Returns the ``GlobalConfig`` and the reference panel, both with
    a concrete build in coordinate-family modes.
    """
    if identity_mode_family(global_config.snp_identifier) == "chr_pos":
        check = GenomeBuildCheck(global_config.genome_build, context='LDScoreCalculator', logger=LOGGER)
        for role, obj in (("annotations", annotation_bundle), ("reference panel", ref_panel)):
            snapshot = getattr(obj, 'config_snapshot', None) or getattr(obj, 'global_config', None)
            build = getattr(snapshot, 'genome_build', None)
            if build in {'hg19', 'hg38'}:
                check.declared(role, role, build)
            if role == 'annotations':
                for chrom in obj.chromosomes:
                    check.frames(role, f'annotations [chr{chrom}]', [obj.metadata_for_chromosome(chrom)],
                                 infer_unlabeled=build not in {'hg19', 'hg38'})
            if role == 'reference panel':
                if obj.spec.backend == 'plink':
                    check.plink('--plink-prefix', obj.spec.plink_prefix)
                else:
                    check.r2('--r2-matrix-dir', obj.spec.r2_matrix_dir)
                check.paths('--ref-panel-snps-file', obj.spec.ref_panel_snps_file, restriction_mode=global_config.snp_identifier)
        global_config = dataclass_replace(global_config, genome_build=check.resolve())
    if hasattr(ref_panel, 'global_config') and global_config.genome_build in {'hg19', 'hg38'}:
        ref_panel = copy(ref_panel)
        ref_panel.global_config = dataclass_replace(ref_panel.global_config, genome_build=global_config.genome_build)
    return global_config, ref_panel
