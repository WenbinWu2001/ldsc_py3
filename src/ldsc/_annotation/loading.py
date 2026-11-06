"""Build a complete source-backed annotation dataset in a borrowed workspace."""

from dataclasses import replace

from .gene_lists import GeneCatalog
from .interval_queries import build_query_shards, prepare_bed_queries, require_readable_bed_sources
from .preflight import (
    annotation_columns, gene_list_gate_a_error, inspect_annotation_inputs, require_usable_bed_query_names,
    resolve_projection_build,
)
from .sources import prepare_annotation_sources
from .gene_query_storage import resolve_gene_lists_staged
from .._core.snp_identity import identity_mode_family
from ..errors import LDSCInputError


def build_annotation_shards(spec, config, workspace, *, projection_genome_build=None, issues_path=None):
    """Prepare aligned sources and construct optional interval queries.

    The caller owns ``workspace`` and the returned bundle borrows that same
    owner. This preparation API writes no canonical artifacts. Standalone
    annotate and direct LD-score workflows own their additional validation
    gates and publication policies.

    When ``issues_path`` is given, unusable BED query names are written there
    as ``annotation_name_collision`` rows before the error is raised, and a
    passing check removes the file. Either way, the error carries the rows as
    ``input_issues`` for callers that write their own diagnostics.
    """
    declared = inspect_annotation_inputs(baseline=spec.baseline_annot_sources, query=spec.query_annot_sources,
        files=[('BED', spec.query_annot_bed_sources), ('gene lists', spec.query_annot_gene_list_sources)],
        scalar_files=[('gene catalog', spec.gene_coordinate_file), ('control genes', spec.control_gene_list_file)],
        mode=config.snp_identifier, strict_annotations=False)
    baseline, query = declared.baseline, declared.query
    require_usable_bed_query_names(declared.inputs.files['BED'], annotation_columns(declared), issues_path)

    catalog = GeneCatalog.load(spec.gene_coordinate_file) if spec.query_annot_gene_list_sources else None
    requested = projection_genome_build or config.genome_build
    projection_build = resolve_projection_build(
        requested, config.snp_identifier, catalog, context='annotation preparation',
        baseline_option='--baseline-annot-sources', baseline_tokens=spec.baseline_annot_sources,
        query_tokens=spec.query_annot_sources, bed_tokens=spec.query_annot_bed_sources)
    if identity_mode_family(config.snp_identifier) == "chr_pos":
        config = replace(config, genome_build=projection_build)
    if spec.query_annot_gene_list_sources:
        batch = resolve_gene_lists_staged(spec.query_annot_gene_list_sources, catalog, workspace,
                                         existing_columns=annotation_columns(declared),
                                         control_path=spec.control_gene_list_file,
                                         resolution_policy=spec.gene_list_resolution_policy,
                                         gene_exclude_regions=spec.gene_exclude_regions)
        if batch.has_fatal_gate_a_issues:
            raise gene_list_gate_a_error(batch)
    else:
        batch = None
    beds = prepare_bed_queries(declared.inputs.files['BED'], workspace) if spec.query_annot_bed_sources else []
    require_readable_bed_sources(beds)
    prepared = prepare_annotation_sources(workspace, baseline, query, mode=config.snp_identifier,
                                          header_widths=declared.widths)
    source_summary = {"baseline_annot_sources": baseline, "query_annot_sources": query, "padding_bp": spec.padding_bp}
    if batch is not None:
        source_summary.update(query_annot_gene_list_sources=list(spec.query_annot_gene_list_sources),
                              gene_catalog_build=catalog.genome_build, projection_genome_build=projection_build)
    elif spec.query_annot_bed_sources:
        source_summary['query_annot_bed_sources'] = declared.inputs.files['BED']
    bundle = prepared.bundle(workspace, config, source_summary=source_summary)
    if batch is not None:
        bundle = replace(build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp),
                         gene_list_batch=batch)
    elif spec.query_annot_bed_sources:
        bundle = build_query_shards(bundle, bed_sources=beds, padding_bp=spec.padding_bp)
    if not bundle.baseline_columns:
        raise LDSCInputError("Annotation sources must supply at least one baseline annotation column.")
    bundle.validate()
    return bundle
