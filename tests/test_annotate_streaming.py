"""Standalone annotation expectations independent of the projection code."""

import json
import pandas as pd
import pytest

from ldsc.config import GlobalConfig
from ldsc.errors import LDSCInputError


def inputs(tmp_path, sharded=False):
    catalog = tmp_path / 'catalog.tsv'
    pd.DataFrame([['G1', 'ONE', '1', 101, 110, 'hg19'], ['G2', 'TWO', '2', 201, 210, 'hg19']],
                 columns=['gene_id', 'gene_name', 'chrom', 'start', 'end', 'genome_build']).to_csv(catalog, sep='\t', index=False)
    genes = tmp_path / 'pathway.txt'
    genes.write_text('G1\nONE\n')
    bed = tmp_path / 'pathway.bed'
    bed.write_text('chr1\t100\t110\n')
    frame = pd.DataFrame({'CHR': [1,1,1,1,2], 'SNP': ['a','b','c','d','e'], 'POS': [100,101,110,111,201], 'base': 1})
    paths=[]
    for chrom, rows in (frame.groupby('CHR') if sharded else [('all', frame)]):
        path=tmp_path / f'baseline.{chrom}.annot.gz'
        rows.to_csv(path, sep='\t', index=False)
        paths.append(path)
    return paths, genes, bed, catalog


@pytest.mark.parametrize('sharded', [False, True])
def test_gene_and_bed_membership_boundaries_reload_and_lifetime(tmp_path, sharded):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, bed, catalog = inputs(tmp_path, sharded)
    config=GlobalConfig(snp_identifier='rsid')
    with run_annotate(baseline_annot_sources=baseline, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                      padding_bp=0, genome_build='hg19', global_config=config, output_dir=tmp_path/'gene') as bundle:
        assert bundle.chromosomes == ['1', '2']
        assert bundle.read('1', columns=['pathway'])[:,0].tolist() == [0,1,1,0]
        assert bundle.read('2', columns=['pathway'])[:,0].tolist() == [0]
        assert bundle.query_statuses[0].n_annotation_snps == 2
        private = bundle.workspace.path
    assert not private.exists()
    with run_annotate(baseline_annot_sources=baseline, query_annot_bed_sources=[bed], global_config=config, output_dir=tmp_path/'bed'):
        pass
    for chrom in ['1','2']:
        gene_frame = pd.read_csv(tmp_path/'gene'/f'query.{chrom}.annot.gz', sep=r'\s+')
        bed_frame = pd.read_csv(tmp_path/'bed'/f'query.{chrom}.annot.gz', sep=r'\s+')
        pd.testing.assert_frame_equal(gene_frame, bed_frame)
    audit=pd.read_csv(tmp_path/'gene'/'diagnostics'/'gene_list_audit.tsv.gz', sep='\t')
    assert audit['annotation_snp_count'].tolist() == [2,2]
    assert audit.reference_snp_count.isna().all()
    metadata=json.loads((tmp_path/'gene'/'diagnostics'/'metadata.json').read_text())
    assert metadata['projection_genome_build'] == 'hg19'
    assert metadata['genome_build'] is None


def test_gene_input_requires_explicit_padding(tmp_path):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog=inputs(tmp_path)
    with pytest.raises(LDSCInputError, match='padding-bp'):
        run_annotate(baseline_annot_sources=baseline, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                     genome_build='hg19', global_config=GlobalConfig(snp_identifier='rsid'), output_dir=tmp_path/'out')


def test_returned_annotation_uses_persistent_outputs_after_build_scratch_is_released(tmp_path):
    from ldsc import run_annotate

    baseline, genes, _, catalog = inputs(tmp_path)
    output = tmp_path / 'out'
    output.mkdir()
    marker = output / 'RUN_FAILED.txt'
    marker.write_text('earlier failed overwrite')
    bundle = run_annotate(baseline_annot_sources=baseline, query_annot_gene_list_sources=[genes],
                          gene_coordinate_file=catalog, padding_bp=0, genome_build='hg19',
                          global_config=GlobalConfig(snp_identifier='rsid'), output_dir=output)
    assert not marker.exists()
    assert not list(output.glob('.ldsc-annotation-*'))
    assert bundle.n_rows == 5
    assert bundle.gene_list_batch.audit_path.is_file()
    genes.unlink()
    catalog.unlink()
    with bundle:
        assert bundle.read('1', columns=['base', 'pathway']).tolist() == [[1, 0], [1, 1], [1, 1], [1, 0]]
        assert bundle.read('2', columns=['pathway']).tolist() == [[0]]
    assert not list(output.glob('.ldsc-annotation-*'))
    assert (output / 'query.1.annot.gz').is_file()


@pytest.mark.parametrize('policy', ['strict', 'resolved-only'])
def test_empty_and_unsupported_sources_are_skipped_but_supported_query_survives(tmp_path, policy):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog = inputs(tmp_path)
    empty, unsupported = tmp_path/'empty.txt', tmp_path/'unsupported.txt'
    empty.write_text('')
    unsupported.write_text('G2\n')
    frame = pd.read_csv(baseline[0], sep='\t')
    frame.loc[frame.CHR.eq(2), 'POS'] = 500
    frame.to_csv(baseline[0], sep='\t', index=False)
    with run_annotate(baseline_annot_sources=baseline, query_annot_gene_list_sources=[empty, genes, unsupported], gene_coordinate_file=catalog,
                      padding_bp=0, genome_build='hg19', gene_list_resolution_policy=policy,
                      global_config=GlobalConfig(snp_identifier='rsid'), output_dir=tmp_path/'out') as bundle:
        assert bundle.query_columns == ['pathway']
        assert [(s.query,s.status,s.reason) for s in bundle.query_statuses] == [('empty','skipped','empty_gene_list'), ('pathway','ok',''), ('unsupported','skipped','zero_annotation_snps')]
        summary = bundle.gene_list_batch.summary.set_index('query')
        assert summary.loc['unsupported','genes_with_snp_support'] == 0


def test_all_focal_skipped_writes_diagnostics_and_no_scientific_family(tmp_path):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog = inputs(tmp_path)
    genes.write_text('')
    output = tmp_path/'out'
    with pytest.raises(LDSCInputError, match='skipped'):
        run_annotate(baseline_annot_sources=baseline, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build='hg19', global_config=GlobalConfig(snp_identifier='rsid'), output_dir=output)
    assert (output/'diagnostics'/'query_annotation_status.tsv').is_file()
    assert not list(output.glob('query.*.annot.gz'))
    assert not list(output.glob('.ldsc-annotation-*'))


@pytest.mark.parametrize('policy', ['strict', 'resolved-only'])
def test_coverage_failure_is_not_treated_as_zero_support(tmp_path, policy):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog = inputs(tmp_path, True)
    genes.write_text('G2\n')
    output = tmp_path/'out'
    with pytest.raises(LDSCInputError, match='coverage'):
        run_annotate(baseline_annot_sources=baseline[:1], query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build='hg19', gene_list_resolution_policy=policy,
                     global_config=GlobalConfig(snp_identifier='rsid'), output_dir=output)
    audit = pd.read_csv(output/'diagnostics'/'gene_list_audit.tsv.gz', sep='\t')
    assert audit.coverage_status.tolist() == ['uncovered']
    assert audit.reference_snp_count.isna().all()
    assert 'annotation_snp_count' not in audit


def test_gate_a_collects_independent_sources_and_omits_later_statuses(tmp_path):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog = inputs(tmp_path)
    genes.write_text('UNKNOWN\nG1\textra\n')
    missing = tmp_path/'missing.txt'
    missing.write_text('G1\textra\n')
    output = tmp_path/'out'
    with pytest.raises(LDSCInputError, match='preflight'):
        run_annotate(baseline_annot_sources=baseline, query_annot_gene_list_sources=[genes,missing], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build='hg19', gene_list_resolution_policy='resolved-only',
                     global_config=GlobalConfig(snp_identifier='rsid'), output_dir=output)
    summary = pd.read_csv(output/'diagnostics'/'gene_list_resolution_summary.tsv', sep='\t').set_index('query')
    assert summary.loc['missing','nonblank_input_rows'] == 1
    audit = pd.read_csv(output/'diagnostics'/'gene_list_audit.tsv.gz', sep='\t')
    assert set(audit.reason) == {'unmatched_identifier','malformed_input'}
    assert not (output/'diagnostics'/'query_annotation_status.tsv').exists()
    assert not list(output.glob('.ldsc-annotation-*'))


def test_explicit_suite_collects_missing_autosomes_and_bad_members(tmp_path):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog = inputs(tmp_path, True)
    (tmp_path/'baseline.3.annot.gz').write_bytes(b'not gzip')
    output=tmp_path/'out'
    with pytest.raises(LDSCInputError, match='[Ii]nput preflight'):
        run_annotate(baseline_annot_sources=[tmp_path/'baseline.@.annot.gz'], query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                     padding_bp=0, genome_build='hg19', global_config=GlobalConfig(snp_identifier='rsid'), output_dir=output)
    issues=pd.read_csv(output/'diagnostics'/'input_issues.tsv', sep='\t')
    assert len(issues) == 20
    assert set(issues.reason) == {'missing_required_input','invalid_required_input'}


def test_overwrite_failure_marker_and_persistent_family(tmp_path):
    from ldsc.annotation_builder import run_annotate

    baseline, genes, _, catalog=inputs(tmp_path)
    output=tmp_path/'out'
    kwargs=dict(baseline_annot_sources=baseline, query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog,
                padding_bp=0, genome_build='hg19', global_config=GlobalConfig(snp_identifier='rsid'), output_dir=output)
    with run_annotate(**kwargs):
        pass
    with pytest.raises(FileExistsError):
        run_annotate(**kwargs)
    unrelated=output/'notes.txt'
    unrelated.write_text('keep')
    genes.write_text('UNKNOWN\n')
    with pytest.raises(LDSCInputError):
        run_annotate(**kwargs, overwrite=True)
    assert (output/'RUN_FAILED.txt').is_file()
    assert unrelated.read_text() == 'keep'
    genes.write_text('G1\n')
    with run_annotate(**kwargs, overwrite=True):
        pass
    assert not (output/'RUN_FAILED.txt').exists()
    assert not list(output.glob('.ldsc-annotation-*'))


def test_real_module_cli_writes_reloadable_outputs_and_releases_scratch(tmp_path):
    import subprocess
    import sys

    baseline, genes, _, catalog=inputs(tmp_path)
    output=tmp_path/'out'
    run=subprocess.run([sys.executable,'-m','ldsc','annotate','--baseline-annot-sources',str(baseline[0]),
                        '--query-annot-gene-list-sources',str(genes),'--gene-coordinate-file',str(catalog),
                        '--padding-bp','0','--snp-identifier','rsid','--genome-build','hg19','--output-dir',str(output)], capture_output=True,text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    frame=pd.read_csv(output/'query.1.annot.gz',sep=r'\s+')
    assert frame.pathway.tolist() == [0,1,1,0]
    assert frame.CM.isna().all()
    assert not list(output.glob('.ldsc-annotation-*'))
