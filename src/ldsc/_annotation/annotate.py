"""Standalone BED/gene annotation over a bounded BIM or annotation grid.

This workflow resolves source and catalog gates before projection, measures
support on the annotation rows, writes canonical chromosome shards
incrementally, and returns a record of the written files. Reference-SNP and
regression-SNP selection are deliberately absent from standalone annotation
generation.
"""

import logging
from dataclasses import replace

import numpy as np
import pandas as pd

from .bim_grid import resolve_bim_grid_inputs, prepare_bim_grid
from .interval_queries import prepare_bed_queries, build_query_shards, require_readable_bed_sources
from .sources import prepare_annotation_sources
from .._runtime.scratch import ScratchWorkspace
from .gene_query_storage import resolve_gene_lists_staged
from .._core.snp_identity import identity_mode_family
from .._runtime.workflow_logging import workflow_logging, materializing_overwrite_guard, log_inputs, log_outputs
from .preflight import (
    annotation_columns, gene_list_gate_a_error, inspect_annotation_inputs, require_usable_bed_query_names,
    resolve_projection_build,
)
from ..config import AnnotationBuildConfig, GlobalConfig, get_global_config, print_global_config_banner
from ..errors import LDSCInputError, LDSCUsageError
from .gene_lists import GeneCatalog, GeneCatalogValidationError
from .outputs import AnnotationDirectoryWriter
from .saved_results import annotation_result
from .._runtime.paths import ensure_output_directory, split_cli_path_tokens
from .query_outcomes import (
    gene_query_statuses, gene_viability_errors, all_query_annotations_skipped_message,
    bed_query_notice, gene_gate_b_notice, log_notices, resolved_only_notice,
    log_gene_catalog_issues,
    log_gene_list_rejections,
    log_gene_list_snp_support,
    log_query_annotation_statuses,
)

LOGGER = logging.getLogger('LDSC.annotation')


class _AnnotationRowSupport:
    """Measure each selected gene's annotation-row SNP support during query projection.

    ``build_query_shards`` calls ``measure`` with each chromosome's projector,
    so support needs no second pass over the grid, and then ``finalize``
    with the projected per-query SNP totals. Finalization records the
    support on the batch and recomputes the gene query statuses from it; a
    usable query whose projection holds no SNP is skipped with
    ``zero_annotation_snps``.
    """

    def __init__(self, batch):
        frame = batch.catalog.frame
        self.batch = batch
        self.selected = np.zeros(len(frame), dtype=bool)
        for selection in batch.selections:
            self.selected[list(selection.catalog_indices)] = True
        self.counts = pd.Series(pd.NA, index=frame.index, dtype='Int64')
        self.counts.loc[self.selected] = 0

    def measure(self, chrom, projector, padding_bp):
        frame = self.batch.catalog.frame
        genes = frame.loc[self.selected & frame.chrom.eq(chrom).to_numpy()]
        self.counts.loc[genes.index] = projector.counts(genes.start0, genes.end, padding_bp=padding_bp)

    def finalize(self, bundle, totals):
        batch = self.batch.with_snp_support(self.counts, support_kind='annotation')
        finalized = []
        for status in gene_query_statuses(batch):
            if status.status in {'ok', 'warning'}:
                count = totals[status.query]
                if count == 0:
                    status = status.updated(status='skipped', reason='zero_annotation_snps')
                status = status.updated(n_annotation_snps=float(count))
            elif status.reason == 'zero_annotation_snps':
                status = status.updated(n_annotation_snps=0.0)
            finalized.append(status)
        return replace(bundle, gene_list_batch=batch, query_statuses=tuple(finalized),
                       query_columns=[item.query for item in finalized if item.status in {'ok', 'warning'}])


@materializing_overwrite_guard(lambda *args, **kwargs: (kwargs.get('output_dir'), kwargs.get('overwrite', False)) if kwargs.get('output_dir') else None, command='run_annotate(...)')
def run_annotate(*, output_dir, ref_snp_annot_sources=None, ref_snp_bim_sources=None,
                 query_annot_bed_sources=None, query_annot_gene_list_sources=None,
                 gene_coordinate_file=None, padding_bp=None, gene_list_resolution_policy=None, gene_exclude_regions=None,
                 genome_build=None, global_config=None, overwrite=False):
    """Write BED or gene-list query annotations and return a record of the written files.

    Parameters
    ----------
    ref_snp_annot_sources, ref_snp_bim_sources : str or sequence of path-like, optional
        Exactly one grid source. Annotation files supply SNP metadata and real
        baseline values; BIM files supply only SNP metadata and need no BED/FAM
        companions. Exact paths, globs, comma-separated tokens, and @
        chromosome suites are accepted. BIM @ suites require chromosomes 1-22.
        Existing identity cleanup determines the retained grid and row order.
    output_dir : path-like
        Required destination for canonical outputs and owned private staging.
    query_annot_bed_sources, query_annot_gene_list_sources : str or sequence of path-like, optional
        Exactly one query route. Gene lists use exact catalog identifiers.
    gene_coordinate_file : path-like, optional
        Required one-build gene coordinate catalog for gene mode; no packaged default.
    padding_bp : int, optional
        Nonnegative interval padding. Gene mode requires an explicit value;
        omitted BED padding is zero. Gene starts are clipped at zero.
    gene_list_resolution_policy : {'strict', 'resolved-only'} or None, optional
        Gene-only rejection policy. Omitted or ``None`` means ``'strict'``:
        rejected identifiers stop preparation after the batch audit.
        ``'resolved-only'`` omits allowlisted rejected rows and records the
        resolved subset in diagnostics. Naming, structural input, and
        chromosome-coverage errors remain fatal. The CLI selects this policy with the value-free
        ``--allow-unresolved-genes`` flag; the Python parameter remains a string.
    gene_exclude_regions : {'none', 'mhc'}, optional
        Gene-only explicit exclusion before padding, default none.
    genome_build : str, optional
        Projection build assertion, kept separate from build-independent
        rsID identity. A concrete build must agree with applicable input
        declarations and supplies the build for unlabeled inputs without
        inference or conversion. Omitted gene mode uses auto reconciliation.
    global_config : GlobalConfig, optional
        Shared identity and logging configuration; defaults to registered globals.
    overwrite : bool, optional
        Permit replacement of the owned artifact family, default false.

    Returns
    -------
    AnnotationResult
        The exported query files and, for annotation-input mode, aligned
        baseline files, with query statuses, gene-list diagnostics and output
        paths. BIM mode has no baseline columns. Construction scratch is
        released before returning, and the result owns no files. To read
        values, call ``result.prepared_in(directory)``, which prepares the
        exported files in a directory you choose; original inputs are not
        needed. No genotype QC or LD scoring occurs.

    Raises
    ------
    LDSCUsageError
        Both or neither grid source, both or neither query route, no
        ``output_dir``, gene-list options without gene lists, or gene lists
        without ``padding_bp``. No output is written, but with
        ``overwrite=True`` the failed call still leaves ``RUN_FAILED.txt`` in
        ``output_dir``.
    LDSCInputError
        Invalid inputs, query names, gene identifiers or chromosome coverage.

    Notes
    -----
    After input-header and catalog validation, gene identifier and query-name
    failures are collected together before build inference, annotation SNP
    grid preparation, or projection. Query names come from source filenames and
    are checked even for empty or unresolved lists. Failed preflight writes the
    complete gene audit and source summary; no new scientific annotation files
    are published.
    """
    bed = split_cli_path_tokens(query_annot_bed_sources)
    genes = split_cli_path_tokens(query_annot_gene_list_sources)
    baseline = split_cli_path_tokens(ref_snp_annot_sources)
    bim = split_cli_path_tokens(ref_snp_bim_sources)
    if bool(baseline) == bool(bim):
        raise LDSCUsageError('annotate requires exactly one of --ref-snp-annot-sources or --ref-snp-bim-sources.')
    if bool(bed) == bool(genes):
        raise LDSCUsageError('annotate requires exactly one of --query-annot-bed-sources or --query-annot-gene-list-sources.')
    if output_dir is None:
        raise LDSCUsageError('annotate requires --output-dir.')
    if genes and padding_bp is None:
        raise LDSCUsageError('Gene-list annotate requires explicit nonnegative --padding-bp; use 0 for gene bodies.')
    if not genes and (gene_coordinate_file is not None or gene_list_resolution_policy is not None or gene_exclude_regions is not None):
        raise LDSCUsageError('Gene coordinate catalog, resolution, and exclusion options require --query-annot-gene-list-sources.')
    config = global_config or get_global_config()
    requested = genome_build if genome_build is not None else config.genome_build
    print_global_config_banner('run_annotate', config)
    spec = AnnotationBuildConfig(baseline_annot_sources=baseline, query_annot_bed_sources=bed, query_annot_gene_list_sources=genes,
                                 gene_coordinate_file=gene_coordinate_file, padding_bp=0 if padding_bp is None else padding_bp,
                                 gene_list_resolution_policy=gene_list_resolution_policy or 'strict', gene_exclude_regions=gene_exclude_regions or 'none')
    output = ensure_output_directory(output_dir, label='output directory')
    writer = AnnotationDirectoryWriter()
    family = writer.artifact_family(output, gene_lists=bool(genes))
    family.preflight(overwrite=overwrite)
    workspace = ScratchWorkspace(output)
    batch = catalog = bundle = None
    try:
        with workflow_logging('annotate', family.paths['log'], log_level=config.log_level):
            log_inputs(grid_kind='bim' if bim else 'annotation', grid_sources=f"{len(bim or baseline)} input declarations",
                       query_sources=f"{len(genes or bed)} input declarations; see diagnostics/metadata.json for source provenance",
                       padding_bp=spec.padding_bp, output_dir=str(output))
            declared_inputs = inspect_annotation_inputs(baseline=baseline,
                files=[('BED', bed), ('gene lists', genes)], scalar_files=[('gene catalog', spec.gene_coordinate_file)],
                mode=config.snp_identifier, strict_annotations=bool(genes),
                issues_path=output/'diagnostics/input_issues.tsv',
                checks=[('BIM grid', ', '.join(bim), lambda: resolve_bim_grid_inputs(bim))] if bim else ())
            if genes:
                try:
                    catalog = GeneCatalog.load(spec.gene_coordinate_file)
                except GeneCatalogValidationError as exc:
                    writer.write_diagnostics(output, catalog_issues=exc.issues)
                    raise
                # Query names depend on declarations and validated headers, not resolved genes or grid rows.
                batch = resolve_gene_lists_staged(genes, catalog, workspace, existing_columns=annotation_columns(declared_inputs),
                                                  resolution_policy=spec.gene_list_resolution_policy,
                                                  gene_exclude_regions=spec.gene_exclude_regions)
                writer.write_diagnostics(output, batch=batch, catalog_issues=catalog.issues)
                log_gene_catalog_issues(catalog.issues)
                log_gene_list_rejections(batch)
                if batch.has_fatal_gate_a_issues:
                    raise gene_list_gate_a_error(batch)
            else:
                # Names and BED contents are cheap to check; both must pass before the grid is read.
                require_usable_bed_query_names(declared_inputs.inputs.files['BED'], annotation_columns(declared_inputs),
                                               output/'diagnostics/input_issues.tsv')
                bed_sources = prepare_bed_queries(declared_inputs.inputs.files['BED'], workspace)
                require_readable_bed_sources(bed_sources, output/'diagnostics/input_issues.tsv')
            projection_build = resolve_projection_build(
                requested, config.snp_identifier, catalog, context='annotate', baseline_option='--ref-snp-annot-sources',
                baseline_tokens=baseline, bim_tokens=bim, bed_tokens=bed)
            prepared = None
            if bim:
                paths, declared = declared_inputs.inputs.checks['BIM grid']
                prepared = prepare_bim_grid(workspace, paths, declared, mode=config.snp_identifier, autosomes_only=bool(genes))
            else:
                paths, declared = declared_inputs.baseline, declared_inputs.declarations
            resolved_config = GlobalConfig(snp_identifier=config.snp_identifier,
                                           genome_build=projection_build if identity_mode_family(config.snp_identifier) == 'chr_pos' else None,
                                           log_level=config.log_level)
            if prepared is None:
                prepared = prepare_annotation_sources(workspace, paths, [], mode=config.snp_identifier, declared_chromosomes=declared,
                                                       autosomes_only=bool(genes), header_widths=declared_inputs.widths)
            bundle = prepared.bundle(workspace, resolved_config, gene_list_batch=batch)
            scope = {'chromosomes': list(prepared.scope_chromosomes), 'grid_kind': 'bim' if bim else 'annotation',
                     'requires_all_autosomes': bool(declared), 'selection': 'validated_input_contents', 'validation_status': 'passed'}
            if baseline:
                scope['baseline_chromosomes'] = list(prepared.scope_chromosomes)
            if batch is not None:
                batch, coverage_errors = batch.with_coverage(prepared.scope_chromosomes)
                bundle = replace(bundle, gene_list_batch=batch)
                if coverage_errors:
                    scope['validation_status'] = 'failed'
                    writer.write_diagnostics(output, batch=batch, statuses=gene_query_statuses(batch), drops=prepared.drops, scope=scope)
                    raise LDSCInputError('Annotation chromosome coverage preflight failed: ' + '; '.join(coverage_errors[:10]) + '. See diagnostics/chromosome_scope.json and gene-list audit/summary.')
                bundle = build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp,
                                            gene_support=_AnnotationRowSupport(batch))
                writer.write_diagnostics(output, batch=bundle.gene_list_batch, statuses=bundle.query_statuses, drops=prepared.drops, scope=scope)
                errors = gene_viability_errors(bundle.gene_list_batch, bundle.query_statuses)
                if errors:
                    raise LDSCInputError('; '.join(errors))
            else:
                bundle = build_query_shards(bundle, bed_sources=bed_sources, padding_bp=spec.padding_bp)
                if not bundle.query_columns:
                    writer.write_diagnostics(output, statuses=bundle.query_statuses, drops=prepared.drops, scope=scope)
                    raise LDSCInputError(all_query_annotations_skipped_message(bundle.query_statuses))
            provenance = {'query_source_kind': 'gene_list' if genes else 'bed', 'query_sources': genes or bed,
                          'grid_kind': 'bim' if bim else 'annotation', 'ref_snp_annot_sources': baseline, 'ref_snp_bim_sources': bim,
                          'query_annot_bed_sources': bed, 'query_annot_gene_list_sources': genes,
                          'gene_coordinate_file': spec.gene_coordinate_file, 'gene_catalog_build': None if catalog is None else catalog.genome_build,
                          'projection_genome_build': projection_build, 'padding_bp': spec.padding_bp, 'snp_rows': prepared.snp_rows,
                          'gene_exclude_regions': spec.gene_exclude_regions, 'gene_list_resolution_policy': spec.gene_list_resolution_policy}
            bundle = replace(bundle, source_summary=provenance)
            bundle = replace(bundle, output_paths=writer.write(
                bundle, output, overwrite=overwrite, provenance=provenance, scope=scope,
                catalog_issues=None if catalog is None else catalog.issues, preflighted=True))
            if bundle.gene_list_batch is not None:
                log_gene_list_snp_support(bundle.gene_list_batch)
            log_query_annotation_statuses(bundle.query_statuses)
            log_notices(query_notices(bundle))
            log_outputs(**bundle.output_paths)
            counts = prepared.snp_rows
            LOGGER.info("%s SNP rows: input=%s, output=%s, dropped=%s.",
                        'BIM grid' if bim else 'Baseline', counts['input'], counts['output'], counts['dropped'])
            if baseline:
                LOGGER.info("For downstream ldscore, use the exported baseline files with --baseline-annot-sources "
                            "and the exported query files with --query-annot-sources; do not use the original baseline inputs. "
                            "Exported paths are listed above.")
            return annotation_result(bundle)
    except BaseException as exc:
        try:
            writer.write_diagnostics(output, batch=None, drops=getattr(exc, 'annotation_drops', None), input_issues=getattr(exc, 'input_issues', None))
        finally:
            workspace.close()
        raise


def query_notices(result) -> list:
    """Return an annotate run's gene-list and BED query notices in console order, ``None`` where one does not apply.

    ``result`` is the run's annotation bundle or the returned ``AnnotationResult``;
    the workflow logs these notices and ``ldsc annotate`` also prints them.
    """
    return [resolved_only_notice(result.gene_list_batch),
            gene_gate_b_notice(result.query_statuses, result.gene_list_batch, command='annotate'),
            bed_query_notice(result.query_statuses, command='annotate')]
