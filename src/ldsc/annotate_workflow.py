"""Standalone BED/gene annotation over a bounded, validated baseline grid.

This workflow resolves source and catalog gates before projection, measures
annotation-grid support, writes canonical chromosome shards incrementally,
and returns a resource-owning dataset handle. Reference-panel and regression
universes are deliberately absent from standalone annotation generation.
"""

import logging
from dataclasses import replace
from pathlib import Path

from ._annotation_bundle import AnnotationBundle
from ._annotation_queries import prepare_bed_queries, build_query_shards
from ._annotation_sources import prepare_annotation_sources
from ._annotation_storage import AnnotationWorkspace
from ._chr_sampler import sample_frame_from_chr_pattern
from ._gene_query_storage import resolve_gene_lists_staged
from ._kernel.snp_identity import identity_mode_family
from ._logging import workflow_logging, materializing_overwrite_guard, log_inputs, log_outputs
from .annotation_semantics import require_unique_annotation_names
from .column_inference import normalize_genome_build
from .config import AnnotationBuildConfig, GlobalConfig, get_global_config, print_global_config_banner
from .errors import LDSCInputError
from .gene_list_resolver import GeneCatalog, GeneCatalogValidationError
from .genome_build_inference import resolve_genome_build
from .outputs import AnnotationDirectoryWriter
from .path_resolution import ensure_output_directory, resolve_file_group, split_cli_path_tokens
from .query_annotations import gene_query_statuses, gene_viability_errors

LOGGER = logging.getLogger('LDSC.annotation')


def _projection_build(requested, mode, baseline_tokens, catalog):
    requested = normalize_genome_build(requested)
    if catalog is None and identity_mode_family(mode) == 'rsid':
        return None
    if catalog is not None and requested is None:
        requested = 'auto'
    if requested is None:
        raise LDSCInputError('annotate requires --genome-build for chr_pos-family identifiers. Supply auto, hg19, or hg38.')
    if requested == 'auto':
        sample, _ = sample_frame_from_chr_pattern(baseline_tokens, context='annotation inputs')
        resolved = resolve_genome_build('auto', 'chr_pos', sample, context='annotation inputs', logger=LOGGER)
    else:
        resolved = requested
    if catalog is not None and resolved != catalog.genome_build:
        raise LDSCInputError(f'Gene-coordinate catalog build does not match annotation build: catalog={catalog.genome_build}, annotation={resolved}. Supply matching inputs; no implicit liftover is performed.')
    return resolved


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
    return LDSCInputError(f'Gene-list catalog preflight rejected {rejected} submitted row(s). ' + '; '.join(issues) + '. Complete diagnostics: diagnostics/gene_list_audit.tsv.gz and diagnostics/gene_list_resolution_summary.tsv. Repair the listed inputs and rerun. Other causes & fixes: docs/troubleshooting.md#annotate-gene-list-preflight')


@materializing_overwrite_guard(lambda *args, **kwargs: (kwargs.get('output_dir'), kwargs.get('overwrite', False), 'RUN_FAILED.txt') if kwargs.get('output_dir') else None, command='run_annotate(...)')
def run_annotate(*, baseline_annot_sources, output_dir, query_annot_bed_sources=None, query_annot_gene_list_sources=None,
                 gene_coordinate_file=None, padding_bp=None, gene_list_resolution_policy=None, gene_exclude_regions=None,
                 genome_build=None, global_config=None, overwrite=False):
    """Write BED or focal gene-list annotations and return a live shard handle.

    Parameters
    ----------
    baseline_annot_sources : path-like or sequence
        Existing annotation files defining SNP identities, order, and baseline
        columns. Exact paths, globs, and explicit chromosome suites are accepted.
    output_dir : path-like
        Required destination for canonical outputs and owned private staging.
    query_annot_bed_sources, query_annot_gene_list_sources : path-like or sequence, optional
        Exactly one focal query route. Gene lists use exact catalog identifiers.
    gene_coordinate_file : path-like, optional
        Required one-build coordinate catalog for gene mode; no packaged default.
    padding_bp : int, optional
        Nonnegative interval padding. Gene mode requires an explicit value;
        omitted BED padding is zero. Gene starts are clipped at zero.
    gene_list_resolution_policy : {'strict', 'resolved-only'}, optional
        Gene-only rejection policy, default strict.
    gene_exclude_regions : {'none', 'mhc'}, optional
        Gene-only explicit exclusion before padding, default none.
    genome_build : str, optional
        Projection build override, kept separate from build-independent rsID
        identity. Omitted gene mode uses auto, following existing inference rules.
    global_config : GlobalConfig, optional
        Shared identity and logging configuration; defaults to registered globals.
    overwrite : bool, optional
        Permit replacement of the owned artifact family, default false.

    Returns
    -------
    AnnotationBundle
        Persistent query outputs plus original baseline dependencies. Construction
        scratch is released before returning; explicit reads prepare private
        shards on demand under the output directory. Use a context manager or
        call ``close()`` after consumption. Closure preserves public outputs.
    """
    bed = split_cli_path_tokens(query_annot_bed_sources)
    genes = split_cli_path_tokens(query_annot_gene_list_sources)
    baseline = split_cli_path_tokens(baseline_annot_sources)
    if bool(bed) == bool(genes):
        raise LDSCInputError('annotate requires exactly one of --query-annot-bed-sources or --query-annot-gene-list-sources.')
    if not baseline or output_dir is None:
        raise LDSCInputError('annotate requires --baseline-annot-sources and --output-dir.')
    if genes and padding_bp is None:
        raise LDSCInputError('Gene-list annotate requires explicit nonnegative --padding-bp; use 0 for gene bodies.')
    if not genes and (gene_coordinate_file is not None or gene_list_resolution_policy is not None or gene_exclude_regions is not None):
        raise LDSCInputError('Gene catalog, resolution, and exclusion options require --query-annot-gene-list-sources.')
    config = global_config or get_global_config()
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
            log_inputs(baseline_annot_sources=baseline, query_sources=genes or bed, padding_bp=spec.padding_bp, output_dir=str(output))
            from ._input_preflight import inspect_declared_inputs
            declared_inputs = inspect_declared_inputs(baseline=baseline,
                files=[('BED', bed), ('gene lists', genes)], scalar_files=[('gene catalog', spec.gene_coordinate_file)],
                mode=config.snp_identifier, strict_annotations=bool(genes),
                issues_path=output/'diagnostics/input_issues.tsv')
            if genes:
                try:
                    catalog = GeneCatalog.load(spec.gene_coordinate_file)
                except GeneCatalogValidationError as exc:
                    writer.write_diagnostics(output, catalog_issues=exc.issues)
                    raise
                batch = resolve_gene_lists_staged(genes, catalog, workspace, resolution_policy=spec.gene_list_resolution_policy, gene_exclude_regions=spec.gene_exclude_regions)
                writer.write_diagnostics(output, batch=batch, catalog_issues=catalog.issues)
                if batch.has_fatal_gate_a_issues:
                    raise _gate_a_error(batch)
            requested = genome_build if genome_build is not None else config.genome_build
            projection_build = _projection_build(requested, config.snp_identifier, baseline, catalog)
            resolved_config = GlobalConfig(snp_identifier=config.snp_identifier,
                                           genome_build=projection_build if identity_mode_family(config.snp_identifier) == 'chr_pos' else None,
                                           log_level=config.log_level)
            paths, declared, issues = declared_inputs.baseline, declared_inputs.declarations, []
            prepared = prepare_annotation_sources(workspace, paths, [], mode=config.snp_identifier, declared_chromosomes=declared, input_issues=issues, autosomes_only=bool(genes), header_widths=declared_inputs.widths)
            bundle = AnnotationBundle(prepared.shards, list(prepared.baseline_columns), [], workspace,
                                      config_snapshot=resolved_config, identity_drops=prepared.drops, gene_list_batch=batch)
            names = [d['query'] for d in batch.declarations] if batch is not None else [Path(path).stem for path in resolve_file_group(bed, label='BED file')]
            naming_error = None
            try:
                require_unique_annotation_names(bundle.baseline_columns, names)
                if batch is not None and any(name in {'gene_control', 'CHR', 'POS', 'BP', 'SNP', 'CM', 'A1', 'A2', 'MAF'} for name in names):
                    raise LDSCInputError('Gene-list query names collide with reserved annotation metadata/control names. Rename the focal source files.')
            except LDSCInputError as exc:
                naming_error = exc
            if batch is not None and (batch.has_fatal_gate_a_issues or naming_error):
                if naming_error:
                    summary = batch.summary.copy()
                    forbidden = set(bundle.baseline_columns) | {'gene_control', 'CHR', 'POS', 'BP', 'SNP', 'CM', 'A1', 'A2', 'MAF'}
                    invalid_names = summary['query'].isin(forbidden) | summary['query'].duplicated(keep=False)
                    summary.loc[invalid_names, 'source_status'] = 'error'
                    summary.loc[invalid_names, 'source_reasons'] = summary.loc[invalid_names, 'source_reasons'].fillna('').map(lambda text: ';'.join(filter(None, (text, 'annotation_name_collision'))))
                    batch = replace(batch, summary=summary)
                    writer.write_diagnostics(output, batch=batch)
                raise _gate_a_error(batch, naming_error)
            if naming_error:
                raise naming_error
            scope = {'chromosomes': list(prepared.scope_chromosomes), 'baseline_chromosomes': list(prepared.scope_chromosomes),
                     'requires_all_autosomes': bool(declared), 'selection': 'validated_input_contents', 'validation_status': 'passed'}
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
                sources = prepare_bed_queries(declared_inputs.files['BED'], workspace)
                build_query_shards(bundle, bed_sources=sources, padding_bp=spec.padding_bp)
            provenance = {'query_source_kind': 'gene_list' if genes else 'bed', 'query_sources': genes or bed,
                          'baseline_annot_sources': baseline, 'query_annot_bed_sources': bed, 'query_annot_gene_list_sources': genes,
                          'gene_coordinate_file': spec.gene_coordinate_file, 'gene_catalog_build': None if catalog is None else catalog.genome_build,
                          'projection_genome_build': projection_build, 'padding_bp': spec.padding_bp,
                          'gene_exclude_regions': spec.gene_exclude_regions, 'gene_list_resolution_policy': spec.gene_list_resolution_policy}
            bundle.source_summary = provenance
            writer.write(bundle, output, overwrite=overwrite, provenance=provenance, scope=scope,
                         catalog_issues=None if catalog is None else catalog.issues, preflighted=True)
            log_outputs(**bundle.output_paths)
            from ._annotation_outputs import persistent_annotation_bundle
            return persistent_annotation_bundle(bundle, paths, output)
    except BaseException as exc:
        try:
            writer.write_diagnostics(output, batch=None, drops=getattr(exc, 'annotation_drops', None), input_issues=getattr(exc, 'input_issues', None))
        finally:
            workspace.close()
        raise
