"""Standalone BED/gene annotation over a bounded BIM or annotation grid.

This workflow resolves source and catalog gates before projection, measures
annotation-grid support, writes canonical chromosome shards incrementally,
and returns a resource-owning dataset handle. Reference-SNP and regression-SNP
selection are deliberately absent from standalone annotation generation.
"""

import logging
from dataclasses import replace
from pathlib import Path

from ._annotation_bundle import AnnotationBundle
from ._annotation_bim import resolve_bim_grid_inputs, prepare_bim_grid
from ._annotation_queries import prepare_bed_queries, build_query_shards, require_readable_bed_sources
from ._annotation_sources import prepare_annotation_sources
from ._annotation_storage import AnnotationWorkspace
from ._gene_query_storage import resolve_gene_lists_staged
from ._core.snp_identity import identity_mode_family
from ._logging import workflow_logging, materializing_overwrite_guard, log_inputs, log_outputs
from .annotation_semantics import require_unique_annotation_names
from ._core.columns import normalize_genome_build
from .config import AnnotationBuildConfig, GlobalConfig, get_global_config, print_global_config_banner
from .errors import LDSCInputError
from .gene_list_resolver import GeneCatalog, GeneCatalogValidationError
from ._genome_build_policy import GenomeBuildCheck
from .outputs import AnnotationDirectoryWriter
from .path_resolution import ensure_output_directory, resolve_file_group, split_cli_path_tokens
from .query_annotations import (
    gene_query_statuses, gene_viability_errors, _all_query_annotations_skipped_message,
    _emit_bed_query_notice, _emit_gene_gate_b_notice, _emit_resolved_only_notice,
)

LOGGER = logging.getLogger('LDSC.annotation')


def _projection_build(requested, mode, baseline_tokens, catalog, *, bim_tokens=(), bed_tokens=(), query_tokens=()):
    """Reconcile grid/catalog declarations before resolving unlabeled inputs."""
    requested = normalize_genome_build(requested)
    if catalog is None and identity_mode_family(mode) == 'rsid' and requested is None:
        return None
    if catalog is not None and requested is None:
        requested = 'auto'
    if requested is None:
        raise LDSCInputError('annotate requires --genome-build for coordinate-family SNP identifier modes. Supply auto, hg19, or hg38.')
    check = GenomeBuildCheck(requested, context='annotate', logger=LOGGER)
    check.paths('--ref-snp-annot-sources', baseline_tokens)
    check.paths('--query-annot-sources', query_tokens)
    if bim_tokens:
        check.paths('--ref-snp-bim-sources', bim_tokens, bim=True)
    if catalog is not None:
        check.declared('--gene-coordinate-file', getattr(catalog, 'source', 'gene coordinate catalog'), catalog.genome_build)
    check.unchecked('--query-annot-bed-sources', bed_tokens, 'ordinary BED; inherits workflow build')
    return check.resolve()


def _gate_a_error(batch, naming_error=None):
    samples, rejected = [], 0
    for frame in batch.audit_frames():
        invalid = frame.loc[frame.disposition.eq('rejected')]
        rejected += len(invalid)
        for row in invalid.head(max(0, 10-len(samples))).itertuples(index=False):
            samples.append(f'{row.source} line {row.line}: {row.input_gene!r} ({row.reason})')
    sources = batch.summary.loc[batch.summary.source_status.eq('error')]
    issues = [f'{row.source}: {row.source_reasons}' for row in sources.head(10).itertuples(index=False)]
    issues.extend(samples)
    if naming_error:
        issues.append(str(naming_error))
    return LDSCInputError(f'Gene-list catalog preflight rejected {rejected} submitted row(s). ' + '; '.join(issues) + '. Complete diagnostics: diagnostics/gene_list_audit.tsv.gz and diagnostics/gene_list_resolution_summary.tsv. Repair the listed inputs and rerun. Other causes & fixes: docs/user-guide/troubleshooting.md#annotate-gene-list-preflight')


@materializing_overwrite_guard(lambda *args, **kwargs: (kwargs.get('output_dir'), kwargs.get('overwrite', False), 'RUN_FAILED.txt') if kwargs.get('output_dir') else None, command='run_annotate(...)')
def run_annotate(*, output_dir, ref_snp_annot_sources=None, ref_snp_bim_sources=None,
                 query_annot_bed_sources=None, query_annot_gene_list_sources=None,
                 gene_coordinate_file=None, padding_bp=None, gene_list_resolution_policy=None, gene_exclude_regions=None,
                 genome_build=None, global_config=None, overwrite=False, _console_notices=False):
    """Write BED or query gene-list annotations and return a live shard handle.

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
    AnnotationBundle
        Persistent query outputs and, for annotation-input mode, aligned baseline
        outputs. Original inputs are not needed for later reads. BIM mode has no
        baseline columns. Construction scratch is released before returning;
        explicit float32 reads prepare packed-binary/dense-float32 private stores
        on demand under the output directory. Use a context manager or ``close()``;
        closure preserves public outputs. No genotype QC or LD scoring occurs.

    Notes
    -----
    After input-header and catalog validation, gene identifier and query-name
    failures are collected together before build inference, SNP-grid preparation,
    or projection. Query names come from source filenames and are checked even
    for empty or unresolved lists. Failed preflight writes the complete gene
    audit and source summary; no new scientific annotation files are published.
    """
    bed = split_cli_path_tokens(query_annot_bed_sources)
    genes = split_cli_path_tokens(query_annot_gene_list_sources)
    baseline = split_cli_path_tokens(ref_snp_annot_sources)
    bim = split_cli_path_tokens(ref_snp_bim_sources)
    if bool(baseline) == bool(bim):
        raise LDSCInputError('annotate requires exactly one of --ref-snp-annot-sources or --ref-snp-bim-sources.')
    if bool(bed) == bool(genes):
        raise LDSCInputError('annotate requires exactly one of --query-annot-bed-sources or --query-annot-gene-list-sources.')
    if output_dir is None:
        raise LDSCInputError('annotate requires --output-dir.')
    if genes and padding_bp is None:
        raise LDSCInputError('Gene-list annotate requires explicit nonnegative --padding-bp; use 0 for gene bodies.')
    if not genes and (gene_coordinate_file is not None or gene_list_resolution_policy is not None or gene_exclude_regions is not None):
        raise LDSCInputError('Gene coordinate catalog, resolution, and exclusion options require --query-annot-gene-list-sources.')
    config = global_config or get_global_config()
    requested = genome_build if genome_build is not None else config.genome_build
    print_global_config_banner('run_annotate', config)
    spec = AnnotationBuildConfig(baseline_annot_sources=baseline, query_annot_bed_sources=bed, query_annot_gene_list_sources=genes,
                                 gene_coordinate_file=gene_coordinate_file, padding_bp=0 if padding_bp is None else padding_bp,
                                 gene_list_resolution_policy=gene_list_resolution_policy or 'strict', gene_exclude_regions=gene_exclude_regions or 'none')
    output = ensure_output_directory(output_dir, label='output directory')
    writer = AnnotationDirectoryWriter()
    writer.artifact_family(output, gene_lists=bool(genes)).preflight(overwrite=overwrite)
    workspace = AnnotationWorkspace(output)
    batch = catalog = bundle = None
    try:
        with workflow_logging('annotate', output/'diagnostics'/'annotate.log', log_level=config.log_level):
            log_inputs(grid_kind='bim' if bim else 'annotation', grid_sources=f"{len(bim or baseline)} input declarations",
                       query_sources=f"{len(genes or bed)} input declarations; see diagnostics/metadata.json for source provenance",
                       padding_bp=spec.padding_bp, output_dir=str(output))
            from ._annotation_preflight import inspect_annotation_inputs
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
                batch = resolve_gene_lists_staged(genes, catalog, workspace, resolution_policy=spec.gene_list_resolution_policy, gene_exclude_regions=spec.gene_exclude_regions)
                # Names depend on declarations and validated headers, not resolved genes or grid rows.
                names = [declaration['query'] for declaration in batch.declarations]
                baseline_columns = {column for path in declared_inputs.baseline for column in declared_inputs.columns[path]}
                reserved = {'gene_control', 'CHR', 'POS', 'BP', 'SNP', 'CM', 'A1', 'A2', 'MAF'}
                naming_errors = []
                try:
                    require_unique_annotation_names(sorted(baseline_columns), names)
                except LDSCInputError as exc:
                    naming_errors.append(str(exc))
                reserved_queries = sorted(reserved.intersection(names))
                if reserved_queries:
                    naming_errors.append('Gene-list query names collide with reserved annotation metadata/control names: ' + ', '.join(reserved_queries) + '.')
                naming_error = None
                if naming_errors:
                    naming_error = LDSCInputError(' '.join(naming_errors) + ' Rename conflicting query source files to change their query names.')
                    summary = batch.summary.copy()
                    # Duplicate-query reasons are already recorded by the resolver.
                    invalid_names = summary['query'].isin(baseline_columns | reserved)
                    summary.loc[invalid_names, 'source_status'] = 'error'
                    summary.loc[invalid_names, 'source_reasons'] = summary.loc[invalid_names, 'source_reasons'].fillna('').map(lambda text: ';'.join(filter(None, (text, 'annotation_name_collision'))))
                    batch = replace(batch, summary=summary, has_fatal_gate_a_issues=True)
                writer.write_diagnostics(output, batch=batch, catalog_issues=catalog.issues)
                from .query_annotations import _log_gene_catalog_issues, _log_gene_list_rejections
                _log_gene_catalog_issues(catalog.issues)
                _log_gene_list_rejections(batch)
                if batch.has_fatal_gate_a_issues:
                    raise _gate_a_error(batch, naming_error)
            else:
                # BED parsing is cheap; it must pass before the grid is read.
                bed_sources = prepare_bed_queries(declared_inputs.files['BED'], workspace)
                require_readable_bed_sources(bed_sources, output/'diagnostics/input_issues.tsv')
            projection_build = _projection_build(requested, config.snp_identifier, baseline, catalog,
                                                 bim_tokens=bim, bed_tokens=bed)
            prepared = None
            if bim:
                paths, declared = declared_inputs.files['BIM grid']
                prepared = prepare_bim_grid(workspace, paths, declared, mode=config.snp_identifier, autosomes_only=bool(genes))
            else:
                paths, declared = declared_inputs.baseline, declared_inputs.declarations
            resolved_config = GlobalConfig(snp_identifier=config.snp_identifier,
                                           genome_build=projection_build if identity_mode_family(config.snp_identifier) == 'chr_pos' else None,
                                           log_level=config.log_level)
            if prepared is None:
                prepared = prepare_annotation_sources(workspace, paths, [], mode=config.snp_identifier, declared_chromosomes=declared,
                                                       autosomes_only=bool(genes), header_widths=declared_inputs.widths)
            bundle = AnnotationBundle(prepared.shards, list(prepared.baseline_columns), [], workspace,
                                      config_snapshot=resolved_config, identity_drops=prepared.drops, gene_list_batch=batch)
            if batch is None:
                names = [Path(path).stem for path in resolve_file_group(bed, label='BED file')]
                require_unique_annotation_names(bundle.baseline_columns, names)
            scope = {'chromosomes': list(prepared.scope_chromosomes), 'grid_kind': 'bim' if bim else 'annotation',
                     'requires_all_autosomes': bool(declared), 'selection': 'validated_input_contents', 'validation_status': 'passed'}
            if baseline:
                scope['baseline_chromosomes'] = list(prepared.scope_chromosomes)
            if batch is not None:
                batch, coverage_errors = batch.with_coverage(prepared.scope_chromosomes)
                bundle.gene_list_batch = batch
                if coverage_errors:
                    scope['validation_status'] = 'failed'
                    writer.write_diagnostics(output, batch=batch, statuses=gene_query_statuses(batch), drops=prepared.drops, scope=scope)
                    raise LDSCInputError('Annotation chromosome coverage preflight failed: ' + '; '.join(coverage_errors[:10]) + '. See diagnostics/chromosome_scope.json and gene-list audit/summary.')
                build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp)
                writer.write_diagnostics(output, batch=bundle.gene_list_batch, statuses=bundle.query_statuses, drops=prepared.drops, scope=scope)
                errors = gene_viability_errors(bundle.gene_list_batch, bundle.query_statuses)
                if errors:
                    raise LDSCInputError('; '.join(errors))
            else:
                build_query_shards(bundle, bed_sources=bed_sources, padding_bp=spec.padding_bp)
                if not bundle.query_columns:
                    writer.write_diagnostics(output, statuses=bundle.query_statuses, drops=prepared.drops, scope=scope)
                    raise LDSCInputError(_all_query_annotations_skipped_message(bundle.query_statuses))
            provenance = {'query_source_kind': 'gene_list' if genes else 'bed', 'query_sources': genes or bed,
                          'grid_kind': 'bim' if bim else 'annotation', 'ref_snp_annot_sources': baseline, 'ref_snp_bim_sources': bim,
                          'query_annot_bed_sources': bed, 'query_annot_gene_list_sources': genes,
                          'gene_coordinate_file': spec.gene_coordinate_file, 'gene_catalog_build': None if catalog is None else catalog.genome_build,
                          'projection_genome_build': projection_build, 'padding_bp': spec.padding_bp, 'snp_rows': prepared.snp_rows,
                          'gene_exclude_regions': spec.gene_exclude_regions, 'gene_list_resolution_policy': spec.gene_list_resolution_policy}
            bundle.source_summary = provenance
            writer.write(bundle, output, overwrite=overwrite, provenance=provenance, scope=scope,
                         catalog_issues=None if catalog is None else catalog.issues, preflighted=True)
            from .query_annotations import _log_gene_list_snp_support, _log_query_annotation_statuses
            if bundle.gene_list_batch is not None:
                _log_gene_list_snp_support(bundle.gene_list_batch)
            _log_query_annotation_statuses(bundle.query_statuses)
            if _console_notices:
                _emit_resolved_only_notice(bundle.gene_list_batch)
                _emit_gene_gate_b_notice(bundle.query_statuses, bundle.gene_list_batch, command='annotate')
                _emit_bed_query_notice(bundle.query_statuses, command='annotate')
            log_outputs(**bundle.output_paths)
            counts = prepared.snp_rows
            LOGGER.info("%s SNP rows: input=%s, output=%s, dropped=%s.",
                        'BIM grid' if bim else 'Baseline', counts['input'], counts['output'], counts['dropped'])
            if baseline:
                LOGGER.info("For downstream ldscore, use the exported baseline files with --baseline-annot-sources "
                            "and the exported query files with --query-annot-sources; do not use the original baseline inputs. "
                            "Exported paths are listed above.")
            from ._annotation_outputs import persistent_annotation_bundle
            return persistent_annotation_bundle(bundle, output)
    except BaseException as exc:
        try:
            writer.write_diagnostics(output, batch=None, drops=getattr(exc, 'annotation_drops', None), input_issues=getattr(exc, 'input_issues', None))
        finally:
            workspace.close()
        raise
