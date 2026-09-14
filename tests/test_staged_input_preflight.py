"""Input gates finish independent checks before staging or numerical work."""

from pathlib import Path

import pandas as pd
import pytest

from ldsc import gene_ldscore_index as index_module, ldscore_calculator
from ldsc.errors import LDSCInputError
from tests.test_plink_workflow_resolution import write_inputs, index_args


def forbid_work(*args, **kwargs):
    pytest.fail("expensive work ran before the input gate completed")


def test_direct_reports_late_reference_and_annotation_defects_before_staging(tmp_path, monkeypatch):
    from ldsc import _annotation_sources
    write_inputs(tmp_path)
    (tmp_path / "1000G.EUR.QC.22.bed").unlink()
    bad = tmp_path / "bad.annot"
    bad.write_text("SNP value\nrs1 1\n")
    monkeypatch.setattr(_annotation_sources, "_scan", forbid_work)
    from ldsc import GlobalConfig, set_global_config
    set_global_config(GlobalConfig(snp_identifier="chr_pos", genome_build="hg19"))
    with pytest.raises(LDSCInputError) as failure:
        ldscore_calculator.run_ldscore(
            baseline_annot_sources=[str(tmp_path / "baseline.annot"), str(bad)],
            plink_prefix=str(tmp_path / "1000G.EUR.QC.*"),
            ld_wind_cm=.05,
            regr_snps_file=str(tmp_path / "regression.tsv"), regr_snps_exclude_regions="none",
            output_dir=tmp_path / "out",
        )
    assert "22.bed" in str(failure.value) and "bad.annot" in str(failure.value)
    issues = pd.read_csv(tmp_path / "out/diagnostics/input_issues.tsv", sep="\t")
    assert len(issues) >= 2


def test_index_missing_late_component_precedes_any_operator_load(tmp_path, monkeypatch):
    from tests.test_gene_index_streaming import two_chromosome_index
    path = two_chromosome_index(tmp_path)
    (path / "chromosomes/chr22/ldscore_operator.npz").unlink()
    monkeypatch.setattr(index_module, "_load_index_chromosome", forbid_work)
    with pytest.raises(LDSCInputError, match="ldscore_operator"):
        index_module._load_gene_ldscore_index(path, _allow_partial_for_tests=True)


def test_builder_reports_all_declared_inputs_before_staging(tmp_path, monkeypatch):
    from ldsc import _annotation_sources
    write_inputs(tmp_path)
    (tmp_path / "1000G.EUR.QC.22.fam").unlink()
    args = index_args(tmp_path, str(tmp_path / "1000G.EUR.QC.*"), tmp_path / "index")
    args.baseline_annot_sources = [str(tmp_path / "missing.annot")]
    args.keep_indivs_file = str(tmp_path / "missing.keep")
    monkeypatch.setattr(_annotation_sources, "_scan", forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        index_module.run_build_gene_ldscore_index_from_args(args)
    message = str(failure.value)
    assert all(name in message for name in ("22.fam", "missing.annot", "missing.keep"))


def test_annotate_reports_baseline_and_bed_paths_before_staging(tmp_path, monkeypatch):
    from ldsc import run_annotate, GlobalConfig, _annotation_sources
    monkeypatch.setattr(_annotation_sources, "_scan", forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        run_annotate(baseline_annot_sources=[tmp_path / "missing.annot"],
                     query_annot_bed_sources=[tmp_path / "missing.bed"], output_dir=tmp_path / "out",
                     global_config=GlobalConfig(snp_identifier="rsid"))
    assert "missing.annot" in str(failure.value) and "missing.bed" in str(failure.value)


def test_progress_has_phase_boundaries_and_periodic_counts(caplog, monkeypatch):
    import logging
    from ldsc import _progress
    time = [0.0]
    monkeypatch.setattr(_progress, 'perf_counter', lambda: time[0])
    with caplog.at_level(logging.INFO, logger='LDSC'):
        with _progress.PhaseProgress(logging.getLogger('LDSC.test'), 'validation', 'objects', 3) as phase:
            phase.advance(object='first')
            time[0] = 31
            phase.advance(object='second')
            phase.advance(object='third')
    records = [r for r in caplog.records if hasattr(r, 'phase_event')]
    assert [r.phase_event for r in records] == ['start', 'progress', 'complete']
    assert records[1].completed == 2 and records[1].total == 3
    assert records[1].current_object == 'second'
    assert records[-1].completed == 3


def test_header_gate_checks_all_sources_without_scanning(tmp_path, monkeypatch):
    from ldsc import _annotation_sources
    from ldsc._annotation_storage import AnnotationWorkspace
    first, last = tmp_path/'first.annot', tmp_path/'last.annot'
    first.write_text('CHR SNP POS value\n21 rs1 100 1\n')
    last.write_text('CHR SNP wrong value\n22 rs2 200 1\n')
    monkeypatch.setattr(_annotation_sources, '_scan', forbid_work)
    with AnnotationWorkspace(tmp_path) as workspace:
        with pytest.raises(LDSCInputError, match='last.annot'):
            _annotation_sources.prepare_annotation_sources(workspace, [str(first), str(last)], [], mode='rsid')


def test_r2_builder_aggregates_late_trio_and_keep_defects(tmp_path, monkeypatch):
    from ldsc import ref_panel_builder
    write_inputs(tmp_path)
    (tmp_path/'1000G.EUR.QC.22.bed').unlink()
    args = ref_panel_builder.build_parser().parse_args([
        '--plink-prefix', str(tmp_path/'1000G.EUR.QC.*'), '--source-genome-build', 'hg19',
        '--snp-identifier', 'chr_pos', '--ld-wind-kb', '1', '--keep-indivs-file', str(tmp_path/'missing.keep'),
        '--output-dir', str(tmp_path/'out')])
    monkeypatch.setattr(ref_panel_builder.ReferencePanelBuilder, '_build_chromosome', forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        ref_panel_builder.run_build_ref_panel_from_args(args)
    assert '22.bed' in str(failure.value) and 'missing.keep' in str(failure.value)


def test_index_late_metadata_and_early_header_errors_precede_operator_loading(tmp_path, monkeypatch):
    import json
    from tests.test_gene_index_streaming import two_chromosome_index
    path = two_chromosome_index(tmp_path)
    component = path/'chromosomes/chr22/metadata.json'
    metadata = json.loads(component.read_text())
    metadata['chromosome'] = '21'
    component.write_text(json.dumps(metadata))
    (path/'chromosomes/chr21/atoms.parquet').write_bytes(b'not parquet')
    monkeypatch.setattr(index_module, '_load_index_chromosome', forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        index_module._load_gene_ldscore_index(path, _allow_partial_for_tests=True)
    assert len(failure.value.input_issues) == 2
    assert set(failure.value.input_issues.chrom) == {'21', '22'}


def test_builder_aggregates_reference_and_restriction_content_errors(tmp_path, monkeypatch):
    from ldsc import _annotation_sources
    write_inputs(tmp_path)
    (tmp_path/'1000G.EUR.QC.22.bed').write_bytes(b'bad')
    (tmp_path/'regression.tsv').write_text('nonsense\nvalue\n')
    args = index_args(tmp_path, str(tmp_path/'1000G.EUR.QC.*'), tmp_path/'index')
    monkeypatch.setattr(_annotation_sources, '_scan', forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        index_module.run_build_gene_ldscore_index_from_args(args)
    assert {'regression SNPs', 'reference'} <= set(failure.value.input_issues.input_role)


@pytest.mark.parametrize('overwrite', [False, True])
def test_early_failure_preserves_science_and_output_contract(tmp_path, monkeypatch, overwrite):
    from ldsc import GlobalConfig, set_global_config, _annotation_sources
    from ldsc.errors import LDSCUserError
    output = tmp_path/'out'
    output.mkdir()
    artifact = output/'ldscore.baseline.parquet'
    artifact.write_bytes(b'existing scientific output')
    set_global_config(GlobalConfig(snp_identifier='chr_pos', genome_build='hg19'))
    monkeypatch.setattr(_annotation_sources, '_scan', forbid_work)
    with pytest.raises(LDSCUserError if overwrite else FileExistsError) as failure:
        ldscore_calculator.run_ldscore(baseline_annot_sources=[str(tmp_path/'missing.annot')],
            plink_prefix=str(tmp_path/'missing_reference'), ld_wind_kb=1,
            output_dir=output, overwrite=overwrite)
    assert artifact.read_bytes() == b'existing scientific output'
    assert (output/'RUN_FAILED.txt').exists() == overwrite
    assert (output/'diagnostics/input_issues.tsv').exists() == overwrite
    if overwrite:
        assert 'missing.annot' in str(failure.value)
    else:
        assert 'already exists' in str(failure.value)


def test_direct_reuses_discovery_and_reference_validation_with_phase_records(tmp_path, monkeypatch, caplog):
    import logging
    from unittest.mock import Mock
    import numpy as np
    from ldsc import GlobalConfig, set_global_config, path_resolution, _ldscore_preflight
    write_inputs(tmp_path)
    discover = Mock(wraps=path_resolution._discover_plink_members)
    inspect = Mock(wraps=_ldscore_preflight.inspect_plink_inputs)
    # The declaration owner imports discovery once; the content owner receives its result.
    from ldsc import _input_preflight
    monkeypatch.setattr(_input_preflight, '_discover_plink_members', discover)
    monkeypatch.setattr(path_resolution, '_discover_plink_members', discover)
    monkeypatch.setattr(_ldscore_preflight, 'inspect_plink_inputs', inspect)
    set_global_config(GlobalConfig(snp_identifier='chr_pos', genome_build='hg19'))
    with caplog.at_level(logging.INFO, logger='LDSC'):
        result = ldscore_calculator.run_ldscore(
            baseline_annot_sources=[str(tmp_path/'baseline.annot')], query_annot_sources=[str(tmp_path/'query.annot')],
            plink_prefix=str(tmp_path/'1000G.EUR.QC.*'), ld_wind_cm=.05, snp_batch_size=1,
            regr_snps_file=str(tmp_path/'regression.tsv'), regr_snps_exclude_regions='none', output_dir=tmp_path/'out')
    assert discover.call_count == inspect.call_count == 1
    np.testing.assert_allclose(pd.read_parquet(tmp_path/'out/ldscore.baseline.parquet')['base'], 1)
    phases = [(r.phase, r.phase_event) for r in caplog.records if hasattr(r, 'phase_event')]
    assert ('validation', 'complete') in phases and ('computation', 'complete') in phases
    assert phases.index(('validation', 'complete')) < phases.index(('computation', 'start'))


def test_r2_builder_aggregates_map_header_and_late_companion(tmp_path, monkeypatch):
    from ldsc import ref_panel_builder
    write_inputs(tmp_path)
    bad_map = tmp_path/'map.tsv'
    bad_map.write_text('CHR POS wrong\n21 100 0.1\n')
    (tmp_path/'1000G.EUR.QC.22.fam').unlink()
    args = ref_panel_builder.build_parser().parse_args([
        '--plink-prefix', str(tmp_path/'1000G.EUR.QC.*'), '--source-genome-build', 'hg19',
        '--snp-identifier', 'chr_pos', '--ld-wind-kb', '1', '--genetic-map-hg19-sources', str(bad_map),
        '--output-dir', str(tmp_path/'out')])
    monkeypatch.setattr(ref_panel_builder.ReferencePanelBuilder, '_build_chromosome', forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        ref_panel_builder.run_build_ref_panel_from_args(args)
    assert 'map.tsv' in str(failure.value) and '22.fam' in str(failure.value)


def test_empty_annotation_header_stops_before_other_source_scans(tmp_path, monkeypatch):
    from ldsc import _annotation_sources
    from ldsc._annotation_storage import AnnotationWorkspace
    good, empty = tmp_path/'good.annot', tmp_path/'empty.annot'
    good.write_text('CHR SNP POS base\n1 rs1 100 1\n')
    empty.write_text('CHR SNP POS extra\n')
    monkeypatch.setattr(_annotation_sources, '_scan', forbid_work)
    with AnnotationWorkspace(tmp_path/'out') as workspace:
        with pytest.raises(LDSCInputError, match='no SNP rows'):
            _annotation_sources.prepare_annotation_sources(workspace, [good, empty], [], mode='rsid')


def test_declared_header_chromosome_and_columns_fail_before_scans(tmp_path, monkeypatch):
    from ldsc import _annotation_sources
    from ldsc._input_preflight import inspect_declared_inputs
    for chrom in range(1, 23):
        column = 'wrong' if chrom == 21 else 'base'
        member = 1 if chrom == 22 else chrom
        (tmp_path/f'base.{chrom}.annot').write_text(f'CHR SNP POS {column}\n{member} rs{chrom} 100 1\n')
    monkeypatch.setattr(_annotation_sources, '_scan', forbid_work)
    with pytest.raises(LDSCInputError) as failure:
        inspect_declared_inputs(baseline=[str(tmp_path/'base.@.annot')], mode='rsid')
    assert {'21', '22'} <= set(failure.value.input_issues.chrom)


@pytest.mark.parametrize('manifest', ['[]', '{"files": []}'])
def test_malformed_manifest_aggregates_with_other_missing_inputs(tmp_path, manifest):
    from ldsc._input_preflight import inspect_declared_inputs, inspect_artifact_paths
    model = tmp_path/'model'
    model.mkdir()
    (model/'metadata.json').write_text(manifest)
    with pytest.raises(LDSCInputError) as failure:
        inspect_declared_inputs(files=[('sumstats', tmp_path/'missing.tsv')],
            checks=[('model', model, lambda: inspect_artifact_paths(model))])
    assert set(failure.value.input_issues.input_role) == {'sumstats', 'model'}


def test_r2_successful_retry_clears_stale_input_repair_table(tmp_path):
    from ldsc import ref_panel_builder
    write_inputs(tmp_path)
    bed = tmp_path/'1000G.EUR.QC.22.bed'
    original = bed.read_bytes()
    bed.unlink()
    args = ref_panel_builder.build_parser().parse_args([
        '--plink-prefix', str(tmp_path/'1000G.EUR.QC.*'), '--source-genome-build', 'hg19',
        '--snp-identifier', 'chr_pos', '--ld-wind-kb', '.05', '--output-dir', str(tmp_path/'out')])
    with pytest.raises(LDSCInputError):
        ref_panel_builder.run_build_ref_panel_from_args(args)
    issue_path = tmp_path/'out/diagnostics/input_issues.tsv'
    assert issue_path.exists()
    bed.write_bytes(original)
    result = ref_panel_builder.run_build_ref_panel_from_args(args)
    assert result.chromosomes == ['21', '22']
    assert not issue_path.exists()
    assert (tmp_path/'out/hg19/chr22_r2.parquet').exists()
