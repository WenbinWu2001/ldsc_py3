"""Prepare direct LD annotations and scope from one bounded source scan.

The workflow owns the workspace. This module borrows it for global annotation
validation and separate chromosome metadata/value artifacts; reference input
integrity is checked before preparation, and that evidence is reused when
content-derived annotation scope becomes available.
"""

from dataclasses import replace
from pathlib import Path

import numpy as np

from ._annotation_bundle import AnnotationBundle
from ._annotation_preflight import input_issue
from ._annotation_queries import build_query_shards, prepare_bed_queries, query_source_statuses, ProjectedQueries
from ._annotation_sources import prepare_annotation_sources
from ._annotation_storage import AnnotationShard, ColumnStore
from ._ldscore_preflight import validate_direct_scope
from .annotation_semantics import require_unique_annotation_names
from .errors import LDSCInputError
from .outputs import LDScoreDirectoryWriter
from .query_annotations import _gene_list_gate_a_message


def prepare_direct_annotations(args, config, spec, workspace, output_config, *, declared_inputs, reference_inputs, batch=None):
    """Validate, construct, and return one complete shard-backed dataset."""
    if batch is not None and batch.catalog.genome_build != getattr(args,'gene_catalog_build',None):
        raise LDSCInputError(f"Gene-coordinate catalog build does not match the resolved analysis build: catalog={batch.catalog.genome_build}, analysis={getattr(args,'gene_catalog_build',None)}. Use a matching catalog; no implicit liftover is performed.")
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
            validate_direct_scope(args,config,batch,output_config,annotation_issues=records,bed_sources=bed_sources,reference_inputs=reference_inputs)
        raise
    if batch is not None:
        batch = validate_gene_query_names(batch, prepared.baseline_columns)
    if batch is not None and batch.has_fatal_gate_a_issues:
        LDScoreDirectoryWriter().write_gene_list_preflight(batch,output_config)
        raise LDSCInputError(_gene_list_gate_a_message(batch))
    scope = {}
    if has_queries:
        batch, scope = validate_direct_scope(args,config,batch,output_config,
            annotation_sources=prepared,bed_sources=bed_sources,reference_inputs=reference_inputs)
    bundle = AnnotationBundle(prepared.shards,list(prepared.baseline_columns),list(prepared.query_columns),workspace,
        config_snapshot=config,identity_drops=prepared.drops,gene_list_batch=batch,
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
    if batch is not None or bed_sources:
        names = [d['query'] for d in batch.declarations] if batch is not None else [b.query for b in bed_sources]
        require_unique_annotation_names(bundle.baseline_columns,names)
        statuses = query_source_statuses(bed_sources, batch)
        if batch is not None and any(item['input_role'] == 'control' for item in batch.declarations):
            build_query_shards(bundle, gene_batch=batch, padding_bp=spec.padding_bp,
                               evaluate_support=False, query_columns=[])
        bundle.query_preparation = ProjectedQueries(tuple(bed_sources), spec.padding_bp)
        bundle.query_statuses = statuses
        bundle.query_columns = [item.query for item in statuses if item.status in {'ok', 'warning'}]
    bundle.validate()
    return bundle, scope


def prepare_synthetic_base(panel, config, workspace):
    """Stage one all-ones chromosome at a time from retained panel metadata."""
    shards = {}
    for chrom in panel.available_chromosomes():
        metadata = panel.load_metadata(chrom)
        if metadata.empty:
            continue
        root = workspace.path/f'base-{chrom}'
        root.mkdir()
        metadata = metadata.rename(columns={'BP':'POS'})
        columns = ['CHR','SNP','CM','POS',*[c for c in ('A1','A2') if c in metadata]]
        metadata.loc[:,columns].to_parquet(root/'metadata.parquet',index=False)
        values = ColumnStore.create(root/'values.npy',len(metadata),['base'])
        for start in range(0,len(metadata),65536):
            values.write(start,np.ones((min(65536,len(metadata)-start),1),dtype=np.float32))
        shards[str(chrom)] = AnnotationShard(root/'metadata.parquet',(values,),len(metadata))
        del metadata
    if not shards:
        raise LDSCInputError('ldscore could not build the synthetic `base` annotation: no retained reference-panel SNP metadata rows were available. Check the panel and SNP restrictions.')
    return AnnotationBundle(shards,['base'],[],workspace,config_snapshot=config,
        source_summary={'baseline':'synthetic all-ones base annotation from retained reference-panel metadata'})


def validate_gene_query_names(batch, baseline_columns):
    """Record header-detectable gene query collisions alongside identifier issues."""
    summary = batch.summary.copy()
    reserved = {'gene_control', 'CHR', 'SNP', 'POS', 'BP', 'CM', 'A1', 'A2', 'MAF'}
    invalid = summary['query'].isin(baseline_columns) | (
        summary.input_role.eq('focal') & summary['query'].isin(reserved))
    if invalid.any():
        summary.loc[invalid, 'source_status'] = 'error'
        summary.loc[invalid, 'source_reasons'] = summary.loc[invalid, 'source_reasons'].fillna('').map(
            lambda reason: ';'.join(filter(None, (reason, 'annotation_name_collision'))))
        batch = replace(batch, summary=summary, has_fatal_gate_a_issues=True)
    return batch
