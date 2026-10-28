"""Standalone grid/output contracts, with independently expected SNP rows."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import GlobalConfig, run_annotate
from ldsc._annotation.sources import prepare_annotation_sources
from ldsc._annotation.storage import AnnotationWorkspace
from ldsc.errors import LDSCInputError


def test_exported_baseline_query_pair_survives_cleanup_and_input_removal(tmp_path):
    baseline = tmp_path / "input.annot"
    baseline.write_text(
        "CHR SNP BP A1 A2 binary continuous\n"
        "1 duplicate 100 A C 1 0.25\n"
        "1 rs1 110 A C 1 0.5\n"
        "2 duplicate 120 A C 0 0.75\n"
        "2 rs2 130 A C 0 1.25\n"
    )
    bed = tmp_path / "focal.bed"
    bed.write_text("1\t109\t110\n")
    output = tmp_path / "out"
    with run_annotate(
        ref_snp_annot_sources=[baseline], query_annot_bed_sources=[bed],
        global_config=GlobalConfig(snp_identifier="rsid"), output_dir=output,
    ) as bundle:
        baselines = list(map(Path, bundle.output_paths["baseline_annotations"]))
        queries = list(map(Path, bundle.output_paths["query_annotations"]))
        assert len(baselines) == len(queries) == 2
        for path, query in zip(baselines, queries):
            left, right = (pd.read_csv(p, sep=r"\s+") for p in (path, query))
            pd.testing.assert_frame_equal(left.drop(columns=bundle.baseline_columns),
                                          right.drop(columns=bundle.query_columns))
        assert baseline.read_text().count("duplicate") == 2
        baseline.unlink()
        bed.unlink()
        np.testing.assert_array_equal(bundle.read("1"), [[1, .5, 1]])
        np.testing.assert_array_equal(bundle.read("2"), [[0, 1.25, 0]])
        with AnnotationWorkspace(tmp_path / "downstream") as workspace:
            prepared = prepare_annotation_sources(workspace, baselines, queries, mode="rsid")
            assert prepared.shards["1"].metadata().SNP.tolist() == ["rs1"]
            np.testing.assert_array_equal(
                prepared.shards["2"].read(columns=["binary", "continuous", "focal"]),
                [[0, 1.25, 0]],
            )
        metadata = json.loads((output / "diagnostics/metadata.json").read_text())
        assert metadata["snp_rows"] == {
            "input": 4, "output": 2, "dropped": 2,
            "by_chromosome": {
                "1": {"input": 2, "output": 1, "dropped": 1},
                "2": {"input": 2, "output": 1, "dropped": 1},
            },
        }
    log = (output / "diagnostics/annotate.log").read_text()
    assert "input=4, output=2, dropped=2" in log
    assert "--baseline-annot-sources" in log and "--query-annot-sources" in log
    assert "original baseline" in log


def test_bim_only_grid_is_query_only_and_ignores_cm_without_companions(tmp_path):
    bim = tmp_path / "grid.bim"
    bim.write_text(
        "1 duplicate ignored 100 A C\n"
        "1 rs1 ignored 110 A C\n"
        "2 duplicate ignored 120 A C\n"
        "2 rs2 ignored 130 A C\n"
    )
    bed = tmp_path / "focal.bed"
    bed.write_text("1\t109\t110\n")
    with run_annotate(ref_snp_bim_sources=[bim], query_annot_bed_sources=[bed],
                      global_config=GlobalConfig(snp_identifier="rsid"),
                      output_dir=tmp_path / "out") as bundle:
        assert bundle.baseline_columns == []
        assert bundle.query_columns == ["focal"]
        assert "baseline_annotations" not in bundle.output_paths
        assert not list((tmp_path / "out").glob("baseline.*.annot.gz"))
        bim.unlink()
        np.testing.assert_array_equal(bundle.read("1"), [[1]])
        np.testing.assert_array_equal(bundle.read("2"), [[0]])
        assert bundle.read("1").dtype == np.float32
        assert bundle.metadata_for_chromosome("1").CM.isna().all()
        assert bundle.shard("1").stores[0].bitorder == "little"
        assert bundle.summary()["n_rows"] == 2


@pytest.mark.parametrize("sources", [{}, {
    "ref_snp_bim_sources": ["grid.bim"],
    "ref_snp_annot_sources": ["grid.annot"],
}])
def test_grid_sources_are_required_and_mutually_exclusive(tmp_path, sources):
    with pytest.raises(LDSCInputError, match="exactly one"):
        run_annotate(**sources, query_annot_bed_sources=["focal.bed"], output_dir=tmp_path / "out")
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize('mode', ['rsid', 'rsid_allele_aware', 'chr_pos', 'chr_pos_allele_aware'])
def test_bim_gene_and_bed_projection_share_identity_and_boundaries(tmp_path, mode):
    from tests.annotation.test_annotate_streaming import inputs
    _, genes, bed, catalog = inputs(tmp_path)
    bim = tmp_path / 'grid.bim.gz'
    import gzip
    with gzip.open(bim, 'wt') as stream:
        stream.write(''.join(f'1 rs{i} 0 {pos} A C\n' for i, pos in enumerate([100, 101, 110, 111])))
        stream.write('1 bad 0 105 A N\n')
    config = GlobalConfig(snp_identifier=mode, genome_build='hg19' if mode.startswith('chr_pos') else None)
    for route in ('gene', 'bed'):
        options = dict(query_annot_gene_list_sources=[genes], gene_coordinate_file=catalog) if route == 'gene' else dict(query_annot_bed_sources=[bed])
        with run_annotate(ref_snp_bim_sources=[bim], **options, padding_bp=0, genome_build='hg19',
                          global_config=config, output_dir=tmp_path / route) as bundle:
            expected = [0, 1, 1, 0] if mode.endswith('allele_aware') else [0, 1, 1, 0, 1]
            assert bundle.read('1')[:, 0].tolist() == expected
            assert bundle.source_summary['snp_rows']['dropped'] == (1 if mode.endswith('allele_aware') else 0)
    pd.testing.assert_frame_equal(*(pd.read_csv(tmp_path / route / 'query.1.annot.gz', sep='\t') for route in ('gene', 'bed')))


@pytest.mark.parametrize('content, message', [
    ('1 rs1 0 100 A\n', 'six headerless columns'),
    ('1 rs1 0 0 A C\n', 'positive integer'),
    ('1 rs1 0 10.5 A C\n', 'positive integer'),
])
def test_bim_bad_structure_has_actionable_failure(tmp_path, content, message):
    bim, bed = tmp_path / 'grid.bim', tmp_path / 'focal.bed'
    bim.write_text(content)
    bed.write_text('1\t0\t100\n')
    with pytest.raises(LDSCInputError, match=message):
        run_annotate(ref_snp_bim_sources=[bim], query_annot_bed_sources=[bed],
                     global_config=GlobalConfig(snp_identifier='rsid'), output_dir=tmp_path / 'out')
    assert not list((tmp_path / 'out').glob('*.annot.gz'))
    assert (tmp_path / 'out/diagnostics/input_issues.tsv').is_file()
    assert not list((tmp_path / 'out').glob('.ldsc-annotation-*'))


def test_bim_auto_build_uses_packaged_evidence(tmp_path):
    from ldsc._genome.build_inference import load_packaged_reference_table
    reference = load_packaged_reference_table().head(500)
    bim, bed = tmp_path / 'grid.bim', tmp_path / 'focal.bed'
    bim.write_text(''.join(f'{row.CHR} rs{i} 0 {row.hg19_POS} A C\n' for i, row in enumerate(reference.itertuples(index=False))))
    bed.write_text(f'{reference.iloc[0].CHR}\t0\t1000000000\n')
    with run_annotate(ref_snp_bim_sources=[bim], query_annot_bed_sources=[bed],
                     global_config=GlobalConfig(snp_identifier='chr_pos', genome_build='auto'), output_dir=tmp_path / 'out') as bundle:
        assert bundle.source_summary['projection_genome_build'] == 'hg19'


def test_cli_switch_to_bim_cleans_stale_baselines_and_preserves_unrelated_file(tmp_path):
    from ldsc._annotation.builder import main
    baseline, bim, bed = tmp_path / 'grid.annot', tmp_path / 'grid.bim', tmp_path / 'focal.bed'
    baseline.write_text('CHR SNP BP base\n1 rs1 100 1\n')
    bim.write_text('1 rs1 0 100 A C\n')
    bed.write_text('1\t99\t100\n')
    output = tmp_path / 'out'
    common = ['--query-annot-bed-sources', str(bed), '--snp-identifier', 'rsid', '--output-dir', str(output)]
    with main(['--ref-snp-annot-sources', str(baseline), *common]):
        pass
    marker = output / 'notes.txt'
    marker.write_text('keep')
    with main(['--ref-snp-bim-sources', str(bim), *common, '--overwrite']) as bundle:
        assert bundle.baseline_columns == []
    assert not (output / 'baseline.1.annot.gz').exists()
    assert marker.read_text() == 'keep'
    assert 'input=1, output=1, dropped=0' in (output / 'diagnostics/annotate.log').read_text()


def test_exported_pair_computes_known_ld_scores_after_cleanup(tmp_path, monkeypatch):
    from ldsc import config, run_ldscore
    from tests.ldscore.test_ldscore_chromosome_scope import write_r2_inputs
    monkeypatch.setattr(config, '_GLOBAL_CONFIG', GlobalConfig(snp_identifier='rsid'))
    args = write_r2_inputs(tmp_path)
    baseline = Path(args.baseline_annot_sources)
    original = pd.read_csv(baseline, sep='\t')
    duplicates = original.iloc[[0, 0]].assign(SNP='duplicate', POS=[40, 50])
    pd.concat([original, duplicates], ignore_index=True).to_csv(baseline, sep='\t', index=False)
    bed = tmp_path / 'focal.bed'
    bed.write_text('22\t19\t20\n')
    with run_annotate(ref_snp_annot_sources=[baseline], query_annot_bed_sources=[bed], output_dir=tmp_path / 'annot') as bundle:
        baseline.unlink()
        regression = tmp_path / 'regression.txt'
        regression.write_text('SNP\nrs1\nrs2\nrs3\n')
        result = run_ldscore(baseline_annot_sources=bundle.output_paths['baseline_annotations'],
                            query_annot_sources=bundle.output_paths['query_annotations'], r2_dir=args.r2_dir,
                            regr_snps_file=regression, regr_snps_exclude_regions='none',
                            ld_wind_snps=2, yes_really=True, output_dir=tmp_path / 'ldscore')
    assert result.read_queries(['focal']).SNP.tolist() == ['rs1', 'rs2', 'rs3']
    np.testing.assert_allclose(result.read_queries(['focal']).focal, [.5, 1, .25], atol=2e-5, rtol=0)


def test_aligned_column_sources_count_logical_rows_once(tmp_path):
    left, right, bed = tmp_path / 'left.annot', tmp_path / 'right.annot', tmp_path / 'focal.bed'
    left.write_text('CHR SNP BP base\n1 rs1 100 1\n2 rs2 200 1\n')
    right.write_text('CHR SNP BP continuous\n1 rs1 100 0.125\n2 rs2 200 2.25\n')
    bed.write_text('1\t99\t100\n')
    with run_annotate(ref_snp_annot_sources=[left, right], query_annot_bed_sources=[bed],
                      global_config=GlobalConfig(snp_identifier='rsid'), output_dir=tmp_path / 'out') as bundle:
        assert bundle.source_summary['snp_rows']['input'] == 2
        assert bundle.source_summary['snp_rows']['dropped'] == 0
        np.testing.assert_array_equal(bundle.read('2'), [[1, 2.25, 0]])


def test_bim_suite_requires_all_declared_chromosomes_and_matching_content(tmp_path):
    from ldsc._annotation.bim_grid import prepare_bim_grid, resolve_bim_grid_inputs
    token = str(tmp_path / 'grid.@.bim')
    with pytest.raises(LDSCInputError, match='missing'):
        resolve_bim_grid_inputs([token])
    for chrom in range(1, 23):
        (tmp_path / f'grid.{chrom}.bim').write_text(f'{chrom} rs{chrom} 0 100 A C\n')
    paths, declared = resolve_bim_grid_inputs([token])
    (tmp_path / 'grid.22.bim').write_text('21 rs22 0 100 A C\n')
    with AnnotationWorkspace(tmp_path / 'out') as workspace:
        with pytest.raises(LDSCInputError, match='declared @ chromosome 22'):
            prepare_bim_grid(workspace, paths, declared, mode='rsid', chunk_rows=1)


def test_cli_grid_selection_rejects_both_and_retired_flag():
    from ldsc._annotation.builder import parse_annotate_args
    common = ['--query-annot-bed-sources', 'focal.bed', '--output-dir', 'out']
    for sources in [[], ['--ref-snp-bim-sources', 'grid.bim', '--ref-snp-annot-sources', 'grid.annot'],
                    ['--baseline-annot-sources', 'grid.annot']]:
        with pytest.raises(SystemExit):
            parse_annotate_args([*common, *sources])
