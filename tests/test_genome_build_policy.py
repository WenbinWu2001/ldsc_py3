"""Workflow build policy exercised with independently labeled HM3 coordinates."""

from argparse import Namespace

import pandas as pd
import pytest

from ldsc.errors import LDSCInputError
from ldsc._genome.build_inference import load_packaged_reference_table
from ldsc.ldscore_calculator import _resolve_ldscore_chr_pos_genome_build


def coordinates(chrom, build, n=500):
    reference = load_packaged_reference_table()
    rows = reference.loc[reference.CHR.astype(str).eq(str(chrom))].head(n)
    return pd.DataFrame({"CHR": rows.CHR, "POS": rows[f"{build}_POS"]})


def test_ldscore_collects_all_build_evidence_before_error(tmp_path):
    paths = []
    for chrom, build in [(1, "hg19"), (2, "hg38")]:
        path = tmp_path / f"baseline.{chrom}.annot"
        coordinates(chrom, build).assign(base=1).to_csv(path, sep="\t", index=False)
        paths.append(str(path))
    keep = tmp_path / "small.tsv"
    coordinates(1, "hg19", n=2).to_csv(keep, sep="\t", index=False)
    bed = tmp_path / "regions.bed"
    bed.write_text("1\t0\t100\n")
    args = Namespace(baseline_annot_sources=paths, r2_dir=None, plink_prefix=None,
                     regr_snps_file=str(keep), query_annot_bed_sources=[str(bed)])
    with pytest.raises(LDSCInputError) as caught:
        _resolve_ldscore_chr_pos_genome_build(args, "auto")
    message = str(caught.value)
    for expected in ("hg19:", "hg38:", "unresolved:", "unchecked:",
                     "baseline.1.annot", "baseline.2.annot", "small.tsv", "regions.bed",
                     "--baseline-annot-sources", "--regr-snps-file", "--genome-build hg19"):
        assert expected in message


def test_annotate_exact_chr2_grid_infers_and_publishes_build(tmp_path):
    from ldsc import GlobalConfig, run_annotate

    grid, bed, output = tmp_path / "grid.annot", tmp_path / "region.bed", tmp_path / "out"
    coordinates(2, "hg38").assign(base=1).to_csv(grid, sep="\t", index=False)
    bed.write_text("2\t0\t1000000000\n")
    with run_annotate(ref_snp_annot_sources=[grid], query_annot_bed_sources=[bed],
                      global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="auto"),
                      output_dir=output) as result:
        assert result.config_snapshot.genome_build == "hg38"
        emitted = pd.read_csv(output / "query.2.annot.gz", sep="\t")
        assert emitted.BP.tolist() == coordinates(2, "hg38").POS.tolist()


def test_annotate_checks_every_bim_before_publication(tmp_path):
    from ldsc import GlobalConfig, run_annotate

    bims = []
    for chrom, build in [(1, "hg19"), (2, "hg38")]:
        path = tmp_path / f"chr{chrom}.bim"
        frame = coordinates(chrom, build).assign(SNP=lambda f: "rs" + f.POS.astype(str), CM=0, A1="A", A2="C")
        frame[["CHR", "SNP", "CM", "POS", "A1", "A2"]].to_csv(path, sep="\t", header=False, index=False)
        bims.append(path)
    bed = tmp_path / "region.bed"
    bed.write_text("1\t0\t1000000000\n2\t0\t1000000000\n")
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="Genome-build consistency") as caught:
        run_annotate(ref_snp_bim_sources=bims, query_annot_bed_sources=[bed],
                     global_config=GlobalConfig(snp_identifier="chr_pos", genome_build="auto"), output_dir=output)
    assert "chr1.bim" in str(caught.value) and "chr2.bim" in str(caught.value)
    assert not list(output.glob("*.annot.gz"))


def test_restriction_prefers_matching_build_column(tmp_path):
    from ldsc._core.identifiers import read_snp_restriction_keys
    path = tmp_path / 'keep.tsv'
    pd.DataFrame({'CHR': [1], 'POS': [999], 'hg19_POS': [123], 'hg38_POS': [456]}).to_csv(path, sep='\t', index=False)
    keys = read_snp_restriction_keys(path, 'chr_pos', genome_build='hg19')
    assert keys.keys == {'1:123'}


def test_panel_cm_auto_waits_for_inference(tmp_path):
    from ldsc import ReferencePanelBuildConfig
    config = ReferencePanelBuildConfig(plink_prefix='panel', output_dir=tmp_path,
        source_genome_build='auto', ld_wind_cm=1,
        genetic_map_hg19_sources='hg19.map', genetic_map_hg38_sources='hg38.map')
    assert config.source_genome_build == 'auto'


@pytest.mark.parametrize('requested, fails', [('auto', True), ('hg19', False)])
def test_panel_keep_list_requires_evidence_only_in_auto(tmp_path, requested, fails):
    from ldsc import ReferencePanelBuildConfig, ReferencePanelBuilder, GlobalConfig
    prefix = tmp_path / 'panel'
    frame = coordinates(1, 'hg19').assign(SNP=lambda f: 'rs' + f.POS.astype(str), CM=0, A1='A', A2='C')
    frame[['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2']].to_csv(str(prefix)+'.bim', sep='\t', header=False, index=False)
    keep = tmp_path / 'keep.tsv'
    coordinates(1, 'hg38', 2).to_csv(keep, sep='\t', index=False)
    config = ReferencePanelBuildConfig(plink_prefix=str(prefix), source_genome_build=requested,
        output_dir=tmp_path / 'out', ld_wind_kb=100, ref_panel_snps_file=keep)
    builder = ReferencePanelBuilder(GlobalConfig(snp_identifier='chr_pos'))
    if fails:
        with pytest.raises(LDSCInputError, match='--source-genome-build hg19') as caught:
            builder._resolve_source_genome_build(config, [str(prefix)])
        assert 'keep.tsv' in str(caught.value)
    else:
        concrete = builder._resolve_source_genome_build(config, [str(prefix)])
        state = builder._prepare_build_state(concrete)
        assert state is not None


@pytest.mark.parametrize('requested, fails', [('auto', True), ('hg19', False)])
def test_munge_checks_small_coordinate_keep_list(tmp_path, requested, fails):
    from ldsc import GlobalConfig
    from ldsc._sumstats_input import _resolve_source_coordinates
    source, keep = tmp_path / 'raw.tsv', tmp_path / 'keep.tsv'
    coordinates(1, 'hg19').to_csv(source, sep='\t', index=False)
    coordinates(1, 'hg38', 2).to_csv(keep, sep='\t', index=False)
    config = GlobalConfig(snp_identifier='chr_pos', genome_build=requested)
    if fails:
        with pytest.raises(LDSCInputError, match='--source-genome-build hg19') as caught:
            _resolve_source_coordinates(source, {'CHR': 'CHR', 'POS': 'POS'}, config, None, 0, restriction_path=keep)
        assert 'keep.tsv' in str(caught.value) and 'raw.tsv' in str(caught.value)
    else:
        result = _resolve_source_coordinates(source, {'CHR': 'CHR', 'POS': 'POS'}, config, None, 0, restriction_path=keep)
        assert result['genome_build'] == 'hg19'


def test_munge_trusts_the_selected_build_specific_source_column(tmp_path):
    from ldsc import GlobalConfig
    from ldsc._sumstats_input import _resolve_source_coordinates
    source = tmp_path / 'raw.tsv'
    coordinates(1, 'hg38', 2).rename(columns={'POS': 'hg19_POS'}).to_csv(source, sep='\t', index=False)
    result = _resolve_source_coordinates(source, {'CHR': 'CHR', 'hg19_POS': 'POS'},
        GlobalConfig(snp_identifier='chr_pos', genome_build='auto'), None, 0)
    assert result['genome_build'] == 'hg19'
    assert result['genome_build_inferred'] is False


def test_munge_preview_trusts_the_selected_build_specific_source_column(tmp_path):
    from ldsc.sumstats_munger import main
    source = tmp_path / 'raw.tsv'
    source.write_text('CHR hg19_POS SNP A1 A2 P BETA N\n1 123 rs1 A G 0.05 0.1 1000\n')
    result = main(['--raw-sumstats-file', str(source), '--output-dir', str(tmp_path / 'unused'),
                   '--pos', 'hg19_POS', '--source-genome-build', 'auto', '--output-genome-build', 'hg19',
                   '--infer-only', '--no-snp-restriction'])
    assert result.runnable
    assert result.source_genome_build == 'hg19'
    assert not (tmp_path / 'unused').exists()


def test_converter_explicit_build_is_authoritative():
    from ldsc.legacy_ldscore_converter import _resolve_output_coordinates
    original = coordinates(1, 'hg38')
    emitted, build, provenance = _resolve_output_coordinates(original, 'chr_pos', 'hg19', context='reference')
    assert build == 'hg19'
    pd.testing.assert_frame_equal(emitted, original)
    assert provenance['inferred_build'] is None


@pytest.mark.parametrize('requested', ['auto', 'hg19'])
def test_in_memory_build_labels_are_declarations(requested):
    from ldsc._genome.build_policy import GenomeBuildCheck
    frame = coordinates(1, 'hg19', 2).assign(genome_build='GRCh38')
    check = GenomeBuildCheck(requested, context='in-memory annotations')
    check.declared('reference', 'panel', 'hg19')
    check.frames('annotations', 'frame', [frame])
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence'):
        check.resolve()


def test_munge_preview_rejects_conflicting_source_declaration(tmp_path):
    from ldsc.sumstats_munger import main
    source = tmp_path / 'raw.tsv'
    source.write_text('CHR hg38_POS SNP A1 A2 P BETA N\n1 123 rs1 A G 0.05 0.1 1000\n')
    result = main(['--raw-sumstats-file', str(source), '--output-dir', str(tmp_path / 'unused'),
                   '--pos', 'hg38_POS', '--source-genome-build', 'hg19', '--output-genome-build', 'hg19',
                   '--infer-only', '--no-snp-restriction'])
    assert not result.runnable
    assert 'Conflicting genome-build evidence' in ' '.join(result.notes)
    assert not (tmp_path / 'unused').exists()


def test_converter_validates_weight_build_declaration_before_conversion(tmp_path, monkeypatch):
    from tests.test_legacy_ldscore_converter import _write_unpartitioned_suite, _write_text
    from ldsc.legacy_ldscore_converter import convert_ldsc2_ldscores
    monkeypatch.setattr('ldsc._runtime.worker_runtime._INITIALIZED', True)
    reference, weights = tmp_path / 'reference', tmp_path / 'weights'
    _write_unpartitioned_suite(reference, weights)
    _write_text(weights / 'weights.2.l2.ldscore.gz', '# ldsc:genome_build=hg38\nCHR SNP BP L2\n2 rs2 201 2\n')
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence') as caught:
        convert_ldsc2_ldscores(legacy_reference_dir=reference, legacy_weight_dir=weights,
                              output_dir=tmp_path / 'out', snp_identifier='chr_pos', genome_build='hg19')
    assert 'weights.2.l2.ldscore.gz' in str(caught.value)
    assert not (tmp_path / 'out' / 'metadata.json').exists()


@pytest.mark.parametrize('validator', ['_validate_r2_identity_matches_config', '_validate_metadata_identity_matches_config'])
def test_low_level_reference_validation_does_not_replace_input_build(validator):
    from ldsc import GlobalConfig
    from ldsc._refpanel import adapters as ref_panel
    metadata = {'snp_identifier': 'chr_pos', 'genome_build': 'GRCh38'}
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence'):
        getattr(ref_panel, validator)(metadata, GlobalConfig(snp_identifier='chr_pos', genome_build='hg19'))
    assert metadata['genome_build'] == 'GRCh38'


def test_regression_reports_all_traits_and_panel_without_override():
    from ldsc import GlobalConfig
    from ldsc.regression_runner import _regression_genome_build
    panel = Namespace(config_snapshot=GlobalConfig(snp_identifier='chr_pos', genome_build='hg19'), output_paths={})
    tables = [Namespace(config_snapshot=GlobalConfig(snp_identifier='chr_pos', genome_build=build),
                        source_path=label, trait_name=label) for label, build in [('first', 'hg38'), ('second', 'hg19')]]
    with pytest.raises(LDSCInputError) as caught:
        _regression_genome_build(tables, panel, 'chr_pos', context='rg')
    message = str(caught.value)
    assert all(label in message for label in ('first', 'second', '--ldscore-dir', 'hg19:', 'hg38:'))
    assert '--genome-build' not in message


def test_keep_list_generic_fallback_only_when_matching_build_column_absent(tmp_path):
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'keep.csv'
    frame = coordinates(1, 'hg38').assign(hg19_POS=lambda f: f.POS + 100)
    frame.to_csv(path, index=False)
    check = GenomeBuildCheck('auto', context='test')
    check.declared('--r2-dir', 'panel', 'hg38')
    check.paths('--regr-snps-file', path, restriction_mode='chr_pos')
    assert check.resolve() == 'hg38'


def test_metadata_and_alternative_columns_are_trusted_without_inference(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from ldsc._genome.build_policy import GenomeBuildCheck
    reference, keep = tmp_path / 'panel.parquet', tmp_path / 'keep.tsv'
    table = pa.Table.from_pandas(coordinates(1, 'hg38', 2)).replace_schema_metadata({b'ldsc:genome_build': b'hg19'})
    pq.write_table(table, reference)
    pd.DataFrame({'CHR': [1], 'hg19_POS': [123], 'hg38_POS': [456], 'POS': [999]}).to_csv(keep, sep='\t', index=False)
    check = GenomeBuildCheck('auto', context='test')
    check.table('--r2-dir', reference)
    check.paths('--regr-snps-file', keep, restriction_mode='chr_pos')
    assert check.resolve() == 'hg19'
    assert check.inferences == []


def test_panel_auto_failure_is_written_to_workflow_log(tmp_path, monkeypatch):
    from ldsc._refpanel import builder as ref_panel_builder
    monkeypatch.setattr('ldsc._runtime.worker_runtime._INITIALIZED', True)
    prefix = tmp_path / 'panel'
    frame = coordinates(1, 'hg19').assign(SNP=lambda f: 'rs' + f.POS.astype(str), CM=0, A1='A', A2='C')
    frame[['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2']].to_csv(str(prefix)+'.bim', sep='\t', header=False, index=False)
    (tmp_path / 'panel.fam').write_text('f1 i1 0 0 1 -9\nf2 i2 0 0 2 -9\n')
    (tmp_path / 'panel.bed').write_bytes(b'\x6c\x1b\x01' + bytes(len(frame)))
    keep = tmp_path / 'keep.tsv'
    coordinates(1, 'hg19', 2).to_csv(keep, sep='\t', index=False)
    output = tmp_path / 'out'
    args = ref_panel_builder.build_parser().parse_args([
        '--plink-prefix', str(prefix), '--output-dir', str(output), '--source-genome-build', 'auto',
        '--snp-identifier', 'chr_pos', '--ref-panel-snps-file', str(keep), '--ld-wind-kb', '1'])
    with pytest.raises(LDSCInputError, match='Genome-build consistency'):
        ref_panel_builder.run_build_ref_panel_from_args(args)
    logs = list((output / 'diagnostics').glob('build-r2-panel*.log'))
    assert len(logs) == 1
    report = logs[0].read_text()
    assert 'keep.tsv' in report and 'unresolved:' in report and '--source-genome-build hg19' in report
    assert not list(output.glob('*/chr*_r2.parquet'))


def test_ldscore_plink_only_auto_resolves_from_bim(tmp_path):
    prefix = tmp_path / 'panel'
    frame = coordinates(2, 'hg38').assign(SNP=lambda f: 'rs' + f.POS.astype(str), CM=0, A1='A', A2='C')
    frame[['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2']].to_csv(str(prefix)+'.bim', sep='\t', header=False, index=False)
    (tmp_path / 'panel.bed').write_bytes(b'\x6c\x1b\x01' + bytes(len(frame)))
    (tmp_path / 'panel.fam').write_text('f i 0 0 1 -9\n')
    args = Namespace(plink_prefix=str(prefix), r2_dir=None, snp_identifier='chr_pos')
    assert _resolve_ldscore_chr_pos_genome_build(args, 'auto') == 'hg38'


def test_build_sample_preserves_existing_invalid_coordinate_drop_policy():
    from ldsc._genome.build_policy import GenomeBuildCheck
    raw = pd.concat([pd.DataFrame({'CHR': ['chrUn', None], 'POS': ['bad', 123]}), coordinates(1, 'hg19')], ignore_index=True)
    check = GenomeBuildCheck('auto', context='test')
    check.frames('--raw-sumstats-file', 'raw', [raw])
    assert check.resolve() == 'hg19'


def test_conflicting_r2_footer_declarations_are_both_reported(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'chr1_r2.parquet'
    table = pa.table({'IDX_1': [0]}).replace_schema_metadata({b'ldsc:genome_build': b'hg19', b'ldsc:sorted_by_build': b'hg38'})
    pq.write_table(table, path)
    check = GenomeBuildCheck('auto', context='test')
    check.table('--r2-dir', path)
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence') as caught:
        check.resolve()
    assert 'ldsc:genome_build' in str(caught.value) and 'ldsc:sorted_by_build' in str(caught.value)


def test_r2_sidecar_declaration_participates_with_an_explicit_build(tmp_path):
    import gzip
    import pyarrow as pa
    import pyarrow.parquet as pq
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'chr1_r2.parquet'
    pq.write_table(pa.table({'IDX_1': [0]}).replace_schema_metadata(
        {b'ldsc:genome_build': b'hg19'}), path)
    with gzip.open(tmp_path / 'chr1_meta.tsv.gz', 'wt') as stream:
        stream.write('# ldsc:genome_build=hg38\nCHR\tPOS\n1\t123\n')
    check = GenomeBuildCheck('auto', context='test')
    check.table('--r2-dir', path)
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence') as caught:
        check.resolve()
    assert 'chr1_r2.parquet' in str(caught.value) and 'chr1_meta.tsv.gz' in str(caught.value)
    explicit = GenomeBuildCheck('hg19', context='test')
    explicit.table('--r2-dir', path)
    with pytest.raises(LDSCInputError, match='hg38') as caught:
        explicit.resolve()
    assert 'chr1_meta.tsv.gz' in str(caught.value)
    assert 'hg19' in str(caught.value)


@pytest.mark.parametrize('requested, declared', [('hg19', 'GRCh38'), ('GRCh37', 'hg38')])
def test_explicit_build_rejects_conflicting_declarations(requested, declared):
    from ldsc._genome.build_policy import GenomeBuildCheck
    check = GenomeBuildCheck(requested, context='test')
    check.declared('gene catalog', 'genes.tsv', declared)
    with pytest.raises(LDSCInputError) as caught:
        check.resolve()
    message = str(caught.value)
    assert all(value in message for value in ('genes.tsv', 'hg19', 'hg38', 'No coordinates or build labels were changed'))


def test_explicit_build_accepts_aliases_and_unlabeled_coordinates_without_inference(tmp_path, monkeypatch):
    from ldsc._genome.build_policy import GenomeBuildCheck
    monkeypatch.setattr('ldsc._genome.build_policy.infer_chr_pos_build',
                        lambda *a, **k: pytest.fail('Explicit builds must not infer coordinates'))
    declared, unlabeled = tmp_path / 'declared.tsv', tmp_path / 'unlabeled.tsv'
    declared.write_text('# ldsc:genome_build=GRCh37\nCHR\tPOS\n1\t123\n')
    unlabeled.write_text('CHR\tPOS\n1\t456\n')
    check = GenomeBuildCheck('hg19', context='test')
    check.table('annotations', declared)
    check.table('annotations', unlabeled)
    assert check.resolve() == 'hg19'
    assert check.inferences == []


@pytest.mark.parametrize('declared, fails', [('GRCh37', False), ('hg38', True)])
def test_munge_preserves_bzip2_build_validation(tmp_path, declared, fails):
    import bz2
    from ldsc import GlobalConfig
    from ldsc._sumstats_input import _resolve_source_coordinates
    path = tmp_path / 'raw.tsv.bz2'
    with bz2.open(path, 'wt') as stream:
        stream.write(f'CHR\tPOS\tgenome_build\n1\t123\t{declared}\n')
    config = GlobalConfig(snp_identifier='chr_pos', genome_build='hg19')
    def resolve():
        return _resolve_source_coordinates(str(path), {'CHR': 'CHR', 'POS': 'POS'}, config, 'bz2', 0)
    if fails:
        with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence') as caught:
            resolve()
        assert str(path) in str(caught.value)
    else:
        assert resolve()['genome_build'] == 'hg19'


def test_annotation_sidecar_declarations_must_agree(tmp_path):
    import json
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'baseline.annot'
    path.write_text('CHR\tPOS\tbase\n1\t123\t1\n')
    (tmp_path / 'diagnostics').mkdir()
    (tmp_path / 'diagnostics' / 'metadata.json').write_text(json.dumps({
        'artifact_type': 'annotation_projection', 'files': {'baseline_annotations': [path.name]},
        'projection_genome_build': 'hg19', 'genome_build': 'hg38'}))
    check = GenomeBuildCheck('hg19', context='test')
    check.table('baseline annotations', path)
    with pytest.raises(LDSCInputError, match='Conflicting genome-build evidence'):
        check.resolve()


def test_mixed_catalog_labels_are_conflicts_not_coordinate_alternatives(tmp_path):
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'mixed.tsv'
    path.write_text('CHR\tPOS\tgenome_build\n1\t123\thg19\n2\t456\thg38\n')
    check = GenomeBuildCheck('hg19', context='test')
    check.table('catalog', path)
    with pytest.raises(LDSCInputError, match='mixed.tsv'):
        check.resolve()


@pytest.mark.parametrize('build', ['hg19', 'hg38'])
def test_explicit_build_selects_coordinate_alternatives(tmp_path, build):
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'alternatives.tsv'
    path.write_text('CHR\thg19_POS\thg38_POS\n1\t123\t456\n')
    check = GenomeBuildCheck(build, context='test')
    check.table('restriction', path)
    assert check.resolve() == build


def test_ldscore_explicit_build_checks_annotation_metadata(tmp_path):
    path = tmp_path / 'baseline.annot'
    path.write_text('# ldsc:genome_build=GRCh38\nCHR\tPOS\tbase\n1\t123\t1\n')
    args = Namespace(baseline_annot_sources=[path], snp_identifier='chr_pos')
    with pytest.raises(LDSCInputError, match='baseline.annot'):
        _resolve_ldscore_chr_pos_genome_build(args, 'hg19')


def test_munge_explicit_source_build_checks_selected_column_and_metadata(tmp_path):
    from ldsc import GlobalConfig
    from ldsc._sumstats_input import _resolve_source_coordinates
    source = tmp_path / 'raw.tsv'
    source.write_text('CHR\thg38_POS\n1\t123\n')
    with pytest.raises(LDSCInputError, match='raw.tsv'):
        _resolve_source_coordinates(source, {'CHR': 'CHR', 'hg38_POS': 'POS'},
                                    GlobalConfig(snp_identifier='chr_pos', genome_build='hg19'), None, 0)


def test_panel_explicit_source_build_checks_keep_list_declaration(tmp_path):
    from ldsc import GlobalConfig, ReferencePanelBuilder, ReferencePanelBuildConfig
    keep = tmp_path / 'keep.tsv'
    keep.write_text('# ldsc:genome_build=hg38\nCHR\tPOS\n1\t123\n')
    config = ReferencePanelBuildConfig(plink_prefix='panel', output_dir=tmp_path / 'out',
        source_genome_build='hg19', ref_panel_snps_file=keep, ld_wind_kb=1)
    with pytest.raises(LDSCInputError, match='keep.tsv'):
        ReferencePanelBuilder(GlobalConfig(snp_identifier='chr_pos'))._resolve_source_genome_build(config, [])


def test_declared_footer_and_column_must_agree(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from ldsc._genome.build_policy import GenomeBuildCheck
    path = tmp_path / 'conflicting.parquet'
    pq.write_table(pa.table({'CHR': [1], 'POS': [123], 'genome_build': ['hg38']})
                   .replace_schema_metadata({b'ldsc:genome_build': b'hg19'}), path)
    check = GenomeBuildCheck('hg19', context='test')
    check.table('annotations', path)
    with pytest.raises(LDSCInputError, match='conflicting.parquet'):
        check.resolve()


@pytest.mark.parametrize('requested', ['hg19', 'auto'])
def test_ldscore_selects_one_available_panel_build(tmp_path, requested):
    import pyarrow as pa
    import pyarrow.parquet as pq
    root = tmp_path / 'panel'
    for build in ('hg19', 'hg38'):
        child = root / build
        child.mkdir(parents=True)
        pq.write_table(pa.table({'IDX_1': [0]}).replace_schema_metadata(
            {b'ldsc:genome_build': build.encode()}), child / 'chr1_r2.parquet')
    baseline = tmp_path / 'baseline.annot'
    baseline.write_text('# ldsc:genome_build=hg19\nCHR\tPOS\tbase\n1\t123\t1\n')
    args = Namespace(r2_dir=root, baseline_annot_sources=[baseline], snp_identifier='chr_pos')
    assert _resolve_ldscore_chr_pos_genome_build(args, requested) == 'hg19'
