"""Prepare direct LD annotations and scope from one bounded source scan.

The workflow owns the workspace. This module borrows it for global annotation
validation and separate chromosome metadata/value artifacts; reference input
integrity is checked by the existing exhaustive scope gate.
"""

from dataclasses import replace
from pathlib import Path

import numpy as np

from ._annotation_bundle import AnnotationBundle
from ._annotation_preflight import input_issue, resolve_annotation_inputs
from ._annotation_queries import build_query_shards, prepare_bed_queries
from ._annotation_sources import prepare_annotation_sources
from ._annotation_storage import AnnotationShard, ColumnStore
from ._ldscore_preflight import validate_direct_scope
from .annotation_semantics import require_unique_annotation_names
from .errors import LDSCInputError
from .outputs import LDScoreDirectoryWriter
from .path_resolution import resolve_file_group
from .query_annotations import _gene_list_gate_a_message


def prepare_direct_annotations(args, config, spec, workspace, output_config, *, batch=None):
    """Validate, construct, and return one complete shard-backed dataset."""
    if batch is not None and batch.catalog.genome_build != getattr(args,'gene_catalog_build',None):
        raise LDSCInputError(f"Gene-coordinate catalog build does not match the resolved analysis build: catalog={batch.catalog.genome_build}, analysis={getattr(args,'gene_catalog_build',None)}. Use a matching catalog; no implicit liftover is performed.")
    baseline, declarations, issues = resolve_annotation_inputs(spec.baseline_annot_sources)
    queries, query_declarations, query_issues = resolve_annotation_inputs(spec.query_annot_sources,role='query')
    declarations.update(query_declarations)
    issues.extend(query_issues)
    bed_sources = prepare_bed_queries(resolve_file_group(spec.query_annot_bed_sources,label='BED file'),workspace) if spec.query_annot_bed_sources else []
    has_queries = bool(spec.query_annot_sources or spec.query_annot_bed_sources or spec.query_annot_gene_list_sources)
    try:
        prepared = prepare_annotation_sources(workspace,baseline,queries,mode=config.snp_identifier,
            declared_chromosomes=declarations,input_issues=issues,autosomes_only=has_queries)
    except LDSCInputError as exc:
        if has_queries:
            failure_issues = getattr(exc,'input_issues',None)
            records = failure_issues.to_dict('records') if failure_issues is not None else [input_issue('alignment','annotation sources','','invalid_required_input',exc)]
            validate_direct_scope(args,config,batch,output_config,annotation_issues=records,bed_sources=bed_sources)
        raise
    if batch is not None:
        summary = batch.summary.copy()
        reserved = {'gene_control', 'CHR', 'SNP', 'POS', 'BP', 'CM', 'A1', 'A2', 'MAF'}
        invalid = summary['query'].isin(prepared.baseline_columns) | (
            summary.input_role.eq('focal') & summary['query'].isin(reserved))
        if invalid.any():
            summary.loc[invalid, 'source_status'] = 'error'
            summary.loc[invalid, 'source_reasons'] = summary.loc[invalid, 'source_reasons'].fillna('').map(
                lambda reason: ';'.join(filter(None, (reason, 'annotation_name_collision'))))
            batch = replace(batch, summary=summary, has_fatal_gate_a_issues=True)
    if batch is not None and batch.has_fatal_gate_a_issues:
        LDScoreDirectoryWriter().write_gene_list_preflight(batch,output_config)
        raise LDSCInputError(_gene_list_gate_a_message(batch))
    scope = {}
    if has_queries:
        batch, scope = validate_direct_scope(args,config,batch,output_config,
            annotation_sources=prepared,bed_sources=bed_sources)
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
        build_query_shards(bundle,bed_sources=bed_sources,gene_batch=batch,padding_bp=spec.padding_bp,evaluate_support=False)
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
