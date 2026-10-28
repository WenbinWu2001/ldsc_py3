"""Unusable BED queries stop annotate; skipped queries reach CLI stderr."""

import gzip
from pathlib import Path

import pandas as pd
import pytest

from ldsc import cli
from ldsc._annotation.builder import run_annotate
from ldsc.config import GlobalConfig
from ldsc.errors import LDSCInputError
from tests.test_annotate_streaming import inputs

RSID = GlobalConfig(snp_identifier='rsid')
PLINK = Path(__file__).parent / 'fixtures/minimal_external_resources/plink/hm3_chr22_subset'


def test_annotate_reports_every_unusable_bed_before_reading_the_grid(tmp_path, monkeypatch):
    from ldsc._annotation import annotate as annotate_workflow

    baseline, _, bed, _ = inputs(tmp_path)
    bad_start = tmp_path / 'bad_start.bed'
    bad_start.write_text('chr1\tbad\t110\n')
    reversed_interval = tmp_path / 'reversed.bed'
    reversed_interval.write_text('chr1\t110\t100\n')
    truncated = tmp_path / 'truncated.bed.gz'
    truncated.write_bytes(gzip.compress(b'chr1\t100\t110\n' * 1000)[:-30])

    def grid_must_wait(*args, **kwargs):
        raise AssertionError('the annotation grid was read before BED validation finished')

    monkeypatch.setattr(annotate_workflow, 'prepare_annotation_sources', grid_must_wait)
    output = tmp_path / 'out'
    with pytest.raises(LDSCInputError) as error:
        run_annotate(ref_snp_annot_sources=baseline, global_config=RSID, output_dir=output,
                     query_annot_bed_sources=[bed, bad_start, reversed_interval, truncated])
    message = str(error.value)
    for name in ('bad_start.bed', 'reversed.bed', 'truncated.bed.gz'):
        assert name in message
    assert 'pathway.bed' not in message
    issues = pd.read_csv(output / 'diagnostics/input_issues.tsv', sep='\t')
    assert sorted(issues.source) == ['bad_start.bed', 'reversed.bed', 'truncated.bed.gz']
    assert not list(output.glob('*.annot.gz'))


def test_annotate_fails_when_every_bed_query_is_empty(tmp_path):
    baseline, _, _, _ = inputs(tmp_path)
    empty_paths = [tmp_path / 'first.bed', tmp_path / 'second.bed']
    for path in empty_paths:
        path.write_text('')
    output = tmp_path / 'out'
    with pytest.raises(LDSCInputError, match='all 2 requested query annotations were skipped'):
        run_annotate(ref_snp_annot_sources=baseline, query_annot_bed_sources=empty_paths,
                     global_config=RSID, output_dir=output)
    status = pd.read_csv(output / 'diagnostics/query_annotation_status.tsv', sep='\t')
    assert status[['query', 'status', 'reason']].values.tolist() == [
        ['first', 'skipped', 'empty_input'], ['second', 'skipped', 'empty_input']]
    assert not list(output.glob('query.*.annot.gz'))


def test_annotate_cli_names_skipped_bed_queries_and_python_api_stays_quiet(tmp_path, capsys):
    baseline, _, bed, _ = inputs(tmp_path)
    empty = tmp_path / 'empty.bed'
    empty.write_text('')
    cli.main(['annotate', '--ref-snp-annot-sources', str(baseline[0]), '--snp-identifier', 'rsid',
              '--query-annot-bed-sources', f'{bed},{empty}', '--output-dir', str(tmp_path / 'cli')])
    notice = capsys.readouterr().err
    assert '1 of 2 BED query annotation(s) were skipped and are omitted from annotation outputs' in notice
    assert 'empty=empty_input' in notice
    assert 'the BED file contains no intervals' in notice
    assert 'diagnostics/query_annotation_status.tsv' in notice

    with run_annotate(ref_snp_annot_sources=baseline, query_annot_bed_sources=[bed, empty],
                      global_config=RSID, output_dir=tmp_path / 'api'):
        pass
    assert capsys.readouterr().err == ''


def test_annotate_cli_reports_gene_list_omissions_and_support(tmp_path, capsys):
    baseline, genes, _, catalog = inputs(tmp_path)
    frame = pd.read_csv(baseline[0], sep='\t')
    frame.loc[frame.CHR.eq(2), 'POS'] = 500  # G2 now covers no grid SNP.
    frame.to_csv(baseline[0], sep='\t', index=False)
    genes.write_text('G1\nUNKNOWN\n')
    unsupported = tmp_path / 'unsupported.txt'
    unsupported.write_text('G2\n')
    cli.main(['annotate', '--ref-snp-annot-sources', str(baseline[0]), '--snp-identifier', 'rsid',
              '--query-annot-gene-list-sources', f'{genes},{unsupported}', '--gene-coordinate-file', str(catalog),
              '--padding-bp', '0', '--genome-build', 'hg19', '--allow-unresolved-genes',
              '--output-dir', str(tmp_path / 'out')])
    notice = capsys.readouterr().err
    assert '--allow-unresolved-genes completed with 1 of 3 submitted row(s) omitted' in notice
    assert 'Gene-list SNP-support checks' in notice
    assert 'unsupported=zero_annotation_snps' in notice
    assert 'Skipped queries are omitted from annotation outputs.' in notice
    assert 'LD-score' not in notice


def test_ldscore_cli_names_skipped_bed_queries(tmp_path, capsys):
    bim = pd.read_csv(f'{PLINK}.bim', sep=r'\s+', header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    baseline = tmp_path / 'baseline.22.annot'
    bim.assign(base=1)[['CHR', 'POS', 'SNP', 'CM', 'base']].to_csv(baseline, sep='\t', index=False)
    half = tmp_path / 'half.bed'
    half.write_text(f'chr22\t{bim.POS.iloc[0] - 1}\t{bim.POS.iloc[15]}\n')
    empty = tmp_path / 'empty.bed'
    empty.write_text('')
    regression_snps = tmp_path / 'regression_snps.tsv'  # Fixture IDs are not packaged HM3 rsIDs.
    bim[['SNP']].to_csv(regression_snps, sep='\t', index=False)
    cli.main(['ldscore', '--baseline-annot-sources', str(baseline), '--query-annot-bed-sources', f'{half},{empty}',
              '--plink-prefix', str(PLINK), '--snp-identifier', 'rsid', '--ld-wind-snps', '10',
              '--regr-snps-file', str(regression_snps), '--regr-snps-exclude-regions', 'none',
              '--output-dir', str(tmp_path / 'ld')])
    notice = capsys.readouterr().err
    assert '1 of 2 BED query annotation(s) were skipped and are omitted from LD-score outputs' in notice
    assert 'empty=empty_input' in notice
