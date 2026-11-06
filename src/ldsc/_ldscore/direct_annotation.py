"""Prepare direct LD annotations and scope from one bounded source scan.

The workflow owns the workspace. This module borrows it for global annotation
validation and separate chromosome metadata/value artifacts; reference input
integrity is checked before preparation, and that evidence is reused when
content-derived annotation scope becomes available.
"""

from dataclasses import replace
import logging
from pathlib import Path

import numpy as np

from .._annotation.bundle import AnnotationBundle
from .._runtime.preflight import input_issue
from .._annotation.interval_queries import DeferredQueries, build_query_shards, prepare_bed_queries, query_source_statuses
from .._annotation.sources import prepare_annotation_sources
from .._annotation.storage import AnnotationShard, ColumnStore
from .preflight import validate_direct_scope
from ..errors import EmptyReferenceSNPs, LDSCInputError

LOGGER = logging.getLogger("LDSC.ldscore_calculator")


def prepare_direct_annotations(spec, config, workspace, output_config, *, declared_inputs, reference_inputs, batch=None,
                               gene_catalog_build=None):
    """Validate and construct one shard-backed dataset and its deferred queries.

    ``spec`` is the run's ``AnnotationBuildConfig`` and ``config`` its
    ``GlobalConfig`` with the resolved genome build; ``gene_catalog_build``
    is the gene-list projection build the run resolved for ``batch``.

    Returns the bundle, a ``DeferredQueries`` record or ``None``, and the
    validated chromosome scope. For BED or gene-list queries the bundle
    declares the query columns and statuses while its shards hold only the
    baseline (and gene control) columns: the LD-score producer projects the
    queries one batch at a time from the ``DeferredQueries`` record. Without
    such queries the record is ``None`` and the shards hold every column.
    """
    if batch is not None and batch.catalog.genome_build != gene_catalog_build:
        raise LDSCInputError(f"Gene coordinate catalog build does not match the resolved analysis build: catalog={batch.catalog.genome_build}, analysis={gene_catalog_build}. Use a matching catalog; no implicit liftover is performed.")
    baseline, queries = declared_inputs.baseline, declared_inputs.query
    declarations, issues = declared_inputs.declarations, []
    bed_sources = prepare_bed_queries(declared_inputs.files['query_annot_bed_sources'],workspace) if spec.query_annot_bed_sources else []
    has_queries = bool(spec.query_annot_sources or spec.query_annot_bed_sources or spec.query_annot_gene_list_sources)
    try:
        prepared = prepare_annotation_sources(workspace,baseline,queries,mode=config.snp_identifier,
            declared_chromosomes=declarations,input_issues=issues,autosomes_only=has_queries,
            header_widths=declared_inputs.widths)
    except LDSCInputError as exc:
        if has_queries:
            failure_issues = getattr(exc,'input_issues',None)
            records = failure_issues.to_dict('records') if failure_issues is not None else [input_issue('alignment','annotation sources','','invalid_required_input',exc)]
            validate_direct_scope(spec,batch,output_config,annotation_issues=records,bed_sources=bed_sources,reference_inputs=reference_inputs)
        raise
    scope = {}
    if has_queries:
        batch, scope = validate_direct_scope(spec,batch,output_config,
            annotation_sources=prepared,bed_sources=bed_sources,reference_inputs=reference_inputs)
    bundle = prepared.bundle(workspace, config, gene_list_batch=batch,
        source_summary={
            'baseline_annot_sources':list(spec.baseline_annot_sources),
            'query_annot_sources':list(spec.query_annot_sources),
            'query_annot_bed_sources':list(spec.query_annot_bed_sources),
            'query_annot_gene_list_sources':[Path(p).name for p in spec.query_annot_gene_list_sources],
            'gene_coordinate_file':None if spec.gene_coordinate_file is None else Path(spec.gene_coordinate_file).name,
            'control_gene_list_file':None if spec.control_gene_list_file is None else Path(spec.control_gene_list_file).name,
            'gene_exclude_regions':spec.gene_exclude_regions,'gene_list_resolution_policy':spec.gene_list_resolution_policy,
            'padding_bp':spec.padding_bp,'chromosome_scope':scope,
        })
    if batch is None and not bed_sources:
        bundle.validate()
        return bundle, None, scope
    statuses = query_source_statuses(bed_sources, batch)
    if batch is not None and any(item['input_role'] == 'control' for item in batch.declarations):
        bundle = build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp, statuses=())
    bundle.validate()
    bundle = replace(bundle, query_statuses=statuses,
                     query_columns=[item.query for item in statuses if item.status in {'ok', 'warning'}])
    return bundle, DeferredQueries(tuple(bed_sources), spec.padding_bp), scope


def prepare_synthetic_base(panel, config, workspace):
    """Stage one all-ones chromosome at a time from retained panel metadata.

    A chromosome that the panel's filters leave without SNPs is skipped with a
    WARNING, as in runs with annotation files. Returns the bundle and the
    skipped chromosomes' reference-SNP funnels, which the run publishes.
    """
    shards, skipped_funnels = {}, []
    for chrom in panel.available_chromosomes():
        try:
            metadata = panel.load_metadata(chrom)
        except EmptyReferenceSNPs as exc:
            LOGGER.warning("Skipping chromosome %s: %s", chrom, exc)
            skipped_funnels.append(exc.reference_snp_funnel)
            continue
        root = workspace.path/f'base-{chrom}'
        root.mkdir()
        metadata = metadata.rename(columns={'BP':'POS'}, copy=False)
        columns = ['CHR','SNP','CM','POS',*[c for c in ('A1','A2') if c in metadata]]
        metadata.loc[:,columns].to_parquet(root/'metadata.parquet',index=False)
        values = ColumnStore.create(root/'values.npy',len(metadata),['base'],dtype=bool)
        for start in range(0,len(metadata),65536):
            values.write(start,np.ones((min(65536,len(metadata)-start),1),dtype=bool))
        shards[str(chrom)] = AnnotationShard(root/'metadata.parquet',(values,),len(metadata))
        del metadata
    if not shards:
        raise LDSCInputError('ldscore could not build the synthetic `base` annotation: the reference-panel filters left no SNP on any chromosome, and each skipped chromosome is a WARNING in the log. Check the panel, --ref-panel-snps-file and --maf-min.')
    bundle = AnnotationBundle(shards,['base'],[],workspace,config_snapshot=config,
        source_summary={'baseline':'synthetic all-ones base annotation from retained reference-panel metadata'})
    return bundle, skipped_funnels
