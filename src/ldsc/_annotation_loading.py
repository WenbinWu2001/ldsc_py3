"""Build a complete source-backed annotation dataset in a borrowed workspace."""

from ._annotation_bundle import AnnotationBundle
from ._annotation_queries import build_query_shards, prepare_bed_queries
from ._annotation_sources import prepare_annotation_sources
from ._gene_query_storage import resolve_gene_lists_staged
from .errors import LDSCInputError
from .path_resolution import ANNOTATION_SUFFIXES, resolve_file_group


def build_annotation_shards(spec, config, workspace, *, projection_genome_build=None):
    """Normalize aligned sources once and construct optional interval queries.

    The caller owns ``workspace`` and the returned bundle borrows that same
    owner. This preparation API writes no canonical artifacts. Standalone
    annotate and direct LD-score workflows own their additional validation
    gates and publication policies.
    """
    baseline = resolve_file_group(spec.baseline_annot_sources, suffixes=ANNOTATION_SUFFIXES,
                                  label="baseline annotation", allow_chromosome_suite=True)
    query = resolve_file_group(spec.query_annot_sources, suffixes=ANNOTATION_SUFFIXES,
                               label="query annotation", allow_chromosome_suite=True) if spec.query_annot_sources else []
    prepared = prepare_annotation_sources(workspace, baseline, query, mode=config.snp_identifier)
    bundle = AnnotationBundle(prepared.shards, list(prepared.baseline_columns), list(prepared.query_columns), workspace,
                              config_snapshot=config, identity_drops=prepared.drops,
                              source_summary={"baseline_annot_sources": baseline, "query_annot_sources": query,
                                              "padding_bp": spec.padding_bp})
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
        build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp, evaluate_support=False)
        bundle.gene_list_batch = batch
        bundle.source_summary.update(query_annot_gene_list_sources=list(spec.query_annot_gene_list_sources),
                                     gene_catalog_build=catalog.genome_build, projection_genome_build=projection_build)
    elif spec.query_annot_bed_sources:
        paths = resolve_file_group(spec.query_annot_bed_sources, label="BED file")
        beds = prepare_bed_queries(paths, workspace)
        build_query_shards(bundle, bed_sources=beds, padding_bp=spec.padding_bp, evaluate_support=False)
        bundle.source_summary["query_annot_bed_sources"] = paths
    if not bundle.baseline_columns:
        raise LDSCInputError("Annotation sources must supply at least one baseline annotation column.")
    bundle.validate()
    return bundle
