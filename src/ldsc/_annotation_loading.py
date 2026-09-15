"""Build a complete source-backed annotation dataset in a borrowed workspace."""

from ._annotation_bundle import AnnotationBundle
from ._annotation_queries import build_query_shards, prepare_bed_queries
from ._annotation_sources import prepare_annotation_sources
from ._gene_query_storage import resolve_gene_lists_staged
from .errors import LDSCInputError


def build_annotation_shards(spec, config, workspace, *, projection_genome_build=None):
    """Prepare aligned sources and construct optional interval queries.

    The caller owns ``workspace`` and the returned bundle borrows that same
    owner. This preparation API writes no canonical artifacts. Standalone
    annotate and direct LD-score workflows own their additional validation
    gates and publication policies.
    """
    from ._input_preflight import inspect_declared_inputs
    declared = inspect_declared_inputs(baseline=spec.baseline_annot_sources, query=spec.query_annot_sources,
        files=[('BED', spec.query_annot_bed_sources), ('gene lists', spec.query_annot_gene_list_sources)],
        scalar_files=[('gene catalog', spec.gene_coordinate_file), ('control genes', spec.control_gene_list_file)],
        mode=config.snp_identifier, strict_annotations=False)
    baseline, query = declared.baseline, declared.query
    if spec.query_annot_gene_list_sources:
        from .annotate_workflow import _projection_build, _gate_a_error
        from .gene_list_resolver import GeneCatalog

        catalog = GeneCatalog.load(spec.gene_coordinate_file)
        projection_build = _projection_build(projection_genome_build or config.genome_build,
                                            config.snp_identifier, spec.baseline_annot_sources, catalog)
        batch = resolve_gene_lists_staged(spec.query_annot_gene_list_sources, catalog, workspace,
                                         control_path=spec.control_gene_list_file,
                                         resolution_policy=spec.gene_list_resolution_policy,
                                         gene_exclude_regions=spec.gene_exclude_regions)
        if batch.has_fatal_gate_a_issues:
            error = _gate_a_error(batch)
            error.gene_list_batch = batch
            raise error
    else:
        batch = None
    beds = prepare_bed_queries(declared.files['BED'], workspace) if spec.query_annot_bed_sources else []
    prepared = prepare_annotation_sources(workspace, baseline, query, mode=config.snp_identifier,
                                          header_widths=declared.widths)
    bundle = AnnotationBundle(prepared.shards, list(prepared.baseline_columns), list(prepared.query_columns), workspace,
                              config_snapshot=config, identity_drops=prepared.drops,
                              source_summary={"baseline_annot_sources": baseline, "query_annot_sources": query,
                                              "padding_bp": spec.padding_bp})
    if batch is not None:
        build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp, evaluate_support=False)
        bundle.gene_list_batch = batch
        bundle.source_summary.update(query_annot_gene_list_sources=list(spec.query_annot_gene_list_sources),
                                     gene_catalog_build=catalog.genome_build, projection_genome_build=projection_build)
    elif spec.query_annot_bed_sources:
        build_query_shards(bundle, bed_sources=beds, padding_bp=spec.padding_bp, evaluate_support=False)
        bundle.source_summary['query_annot_bed_sources'] = declared.files['BED']
    if not bundle.baseline_columns:
        raise LDSCInputError("Annotation sources must supply at least one baseline annotation column.")
    bundle.validate()
    return bundle
