"""Regression coverage for no-write munging diagnostics."""
import re
import shlex

import pytest

from ldsc._sumstats import munger as workflow
from ldsc.errors import LDSCUsageError, LDSCUserError


def test_infer_only_rejects_blank_output_before_reading(capsys):
    with pytest.raises(SystemExit) as caught:
        workflow.main(['--raw-sumstats-file', 'absent.tsv', '--infer-only', '--output-dir', ' '])
    assert caught.value.code == 2
    assert 'argument --output-dir: must name a directory' in capsys.readouterr().err


def test_suggested_command_preserves_overrides(tmp_path, capsys):
    raw = tmp_path / 'raw file.tsv'
    # FRQ backs --maf-min, so the report is runnable and the replayed command would run.
    raw.write_text('chromosome_custom SNP BP A1 A2 P effect FRQ\n1 rs1 100 A G 0.05 0.2 0.3\n')
    options = ['--chr', 'chromosome_custom', '--signed-sumstats', 'effect,0',
               '--N-cas', '78576', '--N-con', '77189', '--maf-min', '0.02',
               '--trait-name', 'trait with spaces', '--output-format', 'both',
               '--ignore', 'unused', '--a1', 'A1']
    inference = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
                               '--snp-identifier', 'rsid', '--infer-only', *options])
    assert inference.runnable, inference.missing_fields
    command = capsys.readouterr().out.split('Suggested command:\n')[1]
    tokens = shlex.split(command.replace('\\\n', ''))
    replay = workflow.build_parser().parse_args(tokens[2:])
    expected = workflow.build_parser().parse_args([
        '--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
        '--snp-identifier', 'rsid', *options])
    assert vars(replay) == {**vars(expected), 'sumstats_format': 'plain'}
    assert not (tmp_path / 'out').exists()


def test_suggested_command_lists_flags_in_help_order(tmp_path, capsys):
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P effect\nrs1 A G 0.05 0.2\n')
    # One option from each help section, passed in reverse help order.
    workflow.main(['--log-level', 'ERROR', '--chunksize', '1', '--signed-sumstats', 'effect,0',
                   '--no-snp-restriction', '--output-format', 'both', '--output-dir', str(tmp_path / 'out'),
                   '--N', '1000', '--trait-name', 'trait', '--raw-sumstats-file', str(raw),
                   '--snp-identifier', 'rsid', '--infer-only'])
    command = capsys.readouterr().out.split('Suggested command:\n')[1]
    flags = [token for token in shlex.split(command.replace('\\\n', '')) if token.startswith('--')]
    help_flags = re.findall(r'^  (--[\w-]+)', workflow.build_parser().format_help(), re.M)
    assert flags == sorted(flags, key=help_flags.index)
    assert flags[0] == '--snp-identifier' and flags[-1] == '--log-level'


def test_infer_only_reports_explicit_qc_threshold_without_column(tmp_path, capsys):
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n')
    inference = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
                               '--snp-identifier', 'rsid', '--infer-only', '--maf-min', '0.02'])
    assert not inference.runnable
    assert any('--maf-min' in field for field in inference.missing_fields)
    assert not any('--info-min' in field for field in inference.missing_fields)
    assert '--frq COLUMN' in capsys.readouterr().out
    assert not (tmp_path / 'out').exists()


@pytest.mark.parametrize('compressed', [False, True])
@pytest.mark.parametrize('extra,blocked', [([], True), (['--ignore', 'N', '--N', '100'], False)])
def test_infer_only_checks_entire_columns(tmp_path, capsys, compressed, extra, blocked):
    import gzip
    raw = tmp_path / ('raw.tsv.gz' if compressed else 'raw.tsv')
    contents = 'SNP P Z N unused\nrs1 NA 1 NA .\nrs2 0.05 2 NaN NA\n'
    if compressed:
        with gzip.open(raw, 'wt') as handle:
            handle.write(contents)
    else:
        raw.write_text(contents)
    report = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
                            '--snp-identifier', 'rsid', '--infer-only', '--chunksize', '1', *extra])
    assert report.runnable is not blocked
    output = capsys.readouterr().out
    assert ('All-missing columns: N' in output) is blocked
    assert 'unused' not in output
    assert ('all-missing column: N' in report.missing_fields) is blocked
    assert not (tmp_path / 'out').exists()


@pytest.mark.parametrize('column', ['P', 'Z', 'SNP', 'N', 'MAF', 'INFO'])
def test_infer_only_checks_required_and_filter_columns(tmp_path, column):
    raw = tmp_path / 'raw.tsv'
    headers = ['SNP', 'P', 'Z', 'N', 'MAF', 'INFO']
    values = dict(zip(headers, ['rs1', '0.05', '2', '100', '0.2', '0.9']))
    values[column] = 'NA'
    raw.write_text(' '.join(headers) + '\n' + ' '.join(values.values()) + '\n')
    report = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
                            '--snp-identifier', 'rsid', '--infer-only'])
    assert f'all-missing column: {column}' in report.missing_fields


def test_infer_only_excludes_replaced_alias_and_checks_info_lists(tmp_path):
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP P Z N MAF frequency INFO\nrs1 0.05 2 100 NA 0.2 NA,.,NaN\n')
    report = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
                            '--snp-identifier', 'rsid', '--infer-only', '--frq', 'frequency'])
    assert report.missing_fields == ('all-missing column: INFO',)


def test_infer_only_requires_raw_sumstats_file(tmp_path, capsys):
    with pytest.raises(SystemExit) as error:
        workflow.main(['--infer-only', '--output-dir', str(tmp_path / 'out')])
    assert error.value.code == 2
    assert 'required: --raw-sumstats-file' in capsys.readouterr().err
    assert not (tmp_path / 'out').exists()


@pytest.mark.parametrize('value', ['', ' '])
def test_infer_only_rejects_blank_raw_sumstats_file(tmp_path, value):
    with pytest.raises(LDSCUsageError, match='nonempty --raw-sumstats-file'):
        workflow.main(['--infer-only', '--raw-sumstats-file', value,
                       '--output-dir', str(tmp_path / 'out')])
    assert not (tmp_path / 'out').exists()


_AUDIT_COMMON_OPTIONS = [
    '--trait-name', 'audit trait', '--output-format', 'both',
    '--snp-identifier', 'chr_pos', '--source-genome-build', 'hg19',
    '--output-genome-build', 'hg38', '--N-cas', '40', '--N-con', '60',
    '--N-min', '25', '--nstudy-min', '2', '--info-min', '0.8', '--maf-min', '0.02',
    '--snp', 'variant', '--chr', 'chrom_custom', '--pos', 'base_custom',
    '--a1', 'effect_allele_custom', '--a2', 'other_allele_custom', '--p', 'p_custom',
    '--frq', 'frequency_custom', '--nstudy', 'studies_custom', '--ignore', 'unused_custom',
    '--input-format', 'plain', '--chunksize', '1', '--overwrite', '--log-level', 'ERROR',
]
_AUDIT_VARIANTS = [
    ['--N-col', 'total_custom', '--N', '100', '--info', 'quality_custom',
     '--signed-sumstats', 'effect_custom,0', '--no-snp-restriction'],
    ['--N-cas-col', 'cases_custom', '--N-con-col', 'controls_custom',
     '--info-list', 'list_custom', '--a1-inc'],
]


def test_suggested_command_audit_covers_every_parser_flag():
    exercised = {token for token in [*_AUDIT_COMMON_OPTIONS, *sum(_AUDIT_VARIANTS, [])]
                 if token.startswith('--')}
    exercised.update({'--raw-sumstats-file', '--output-dir', '--infer-only',
                      '--liftover-chain-file', '--sumstats-snps-file', '--help'})
    parser_flags = {flag for action in workflow.build_parser()._actions
                    for flag in action.option_strings if flag.startswith('--')}
    assert exercised == parser_flags


@pytest.mark.parametrize('variant', _AUDIT_VARIANTS)
def test_suggested_command_round_trips_all_flag_families(tmp_path, capsys, variant):
    raw = tmp_path / 'raw file.tsv'
    raw.write_text(
        'variant chrom_custom base_custom effect_allele_custom other_allele_custom '
        'p_custom frequency_custom studies_custom unused_custom total_custom '
        'cases_custom controls_custom quality_custom list_custom effect_custom\n'
        'rs1 1 100 A G 0.05 0.2 3 NA 100 40 60 0.99 0.98,NA 0.2\n'
    )
    chain = tmp_path / 'example chain'
    chain.write_text('')  # Infer-only validates the path without parsing the chain.
    keep = tmp_path / 'keep list.tsv'
    keep.write_text('CHR POS\n1 100\n')
    options = ['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'),
               '--liftover-chain-file', str(chain), *_AUDIT_COMMON_OPTIONS, *variant]
    if '--no-snp-restriction' not in variant:
        options.extend(['--sumstats-snps-file', str(keep)])
    report = workflow.main([*options, '--infer-only'])
    assert report.runnable
    command = capsys.readouterr().out.split('Suggested command:\n')[1]
    tokens = shlex.split(command.replace('\\\n', ''))
    parser = workflow.build_parser()
    assert vars(parser.parse_args(tokens[2:])) == vars(parser.parse_args(options))
    assert not (tmp_path / 'out').exists()


# Each header problem stops the run; --infer-only must report it instead of "Runnable: yes".
_HEADER_PROBLEMS = {
    'two signed statistics': ('SNP A1 A2 P BETA Z N\nrs1 A G 0.05 0.1 1.0 1000\n', [], 'conflict: signed statistic'),
    'two p-value columns': ('SNP A1 A2 P PVAL BETA N\nrs1 A G 0.05 0.05 0.1 1000\n', [], 'conflict: P'),
    'two sample-size strategies': ('SNP A1 A2 P BETA N N_CAS N_CON\nrs1 A G 0.05 0.1 1000 400 600\n', [],
                                   'conflict: sample size'),
    'case column without control column': ('SNP A1 A2 P BETA NCAS NCON\nrs1 A G 0.05 0.1 400 600\n',
                                           ['--N-cas-col', 'NCAS'], 'sample-size options'),
    'forced old DANER without counts': ('SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n', ['--input-format', 'daner-old'],
                                        'DANER case/control counts'),
}


@pytest.mark.parametrize('case', list(_HEADER_PROBLEMS))
def test_infer_only_reports_the_header_problem_the_run_stops_on(tmp_path, capsys, case):
    text, extra, label = _HEADER_PROBLEMS[case]
    raw = tmp_path / 'raw.tsv'
    raw.write_text(text)
    base = ['--raw-sumstats-file', str(raw), '--snp-identifier', 'rsid', '--no-snp-restriction', *extra]
    report = workflow.main([*base, '--output-dir', str(tmp_path / 'infer'), '--infer-only'])
    assert not report.runnable
    assert label in report.missing_fields
    assert 'Runnable: no' in capsys.readouterr().out
    with pytest.raises(LDSCUserError) as caught:
        workflow.main([*base, '--output-dir', str(tmp_path / 'run')])
    assert str(caught.value) in report.notes


def test_infer_only_honours_ignore_for_the_signed_statistic(tmp_path):
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n')
    report = workflow.main(['--raw-sumstats-file', str(raw), '--snp-identifier', 'rsid', '--no-snp-restriction',
                            '--ignore', 'BETA', '--output-dir', str(tmp_path / 'infer'), '--infer-only'])
    assert report.missing_fields == ('signed statistic',)
    assert any('BETA is passed to --ignore' in note for note in report.notes)


# Before, --infer-only named these fields without saying why the run would stop.
_MISSING_FIELDS = {
    'no signed statistic': ('SNP A1 A2 P N\nrs1 A G 0.05 1000\n', 'rsid', [], 'signed statistic'),
    'ignored signed statistic': ('SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n', 'rsid', ['--ignore', 'BETA'],
                                 'signed statistic'),
    'no sample size': ('SNP A1 A2 P BETA\nrs1 A G 0.05 0.1\n', 'rsid', [], 'N'),
    'no alleles in an allele-aware mode': ('SNP P BETA N\nrs1 0.05 0.1 1000\n', 'rsid_allele_aware', [], 'A1/A2'),
}


@pytest.mark.parametrize('case', list(_MISSING_FIELDS))
def test_infer_only_explains_every_missing_field_with_the_runs_message(tmp_path, case):
    text, mode, extra, label = _MISSING_FIELDS[case]
    raw = tmp_path / 'raw.tsv'
    raw.write_text(text)
    base = ['--raw-sumstats-file', str(raw), '--snp-identifier', mode, '--no-snp-restriction', *extra]
    report = workflow.main([*base, '--output-dir', str(tmp_path / 'infer'), '--infer-only'])
    assert report.missing_fields == (label,)
    with pytest.raises(LDSCUserError) as caught:
        workflow.main([*base, '--output-dir', str(tmp_path / 'run')])
    assert str(caught.value) in report.notes


def test_python_inference_reports_header_problems(tmp_path):
    from ldsc import GlobalConfig, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P BETA Z N\nrs1 A G 0.05 0.1 1.0 1000\n')
    report = infer_raw_sumstats(raw, global_config=GlobalConfig(snp_identifier='rsid'))
    assert report.missing_fields == ('conflict: signed statistic',)


def test_run_lists_every_header_problem_in_one_error(tmp_path):
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P effect_size\nrs1 A G 0.05 0.1\n')
    with pytest.raises(LDSCUserError) as caught:
        workflow.main(['--raw-sumstats-file', str(raw), '--snp-identifier', 'rsid', '--no-snp-restriction',
                       '--output-dir', str(tmp_path / 'out')])
    message = str(caught.value)
    assert message.startswith(f"munge-sumstats found 2 problems with the columns of '{raw}'")
    assert '\n  1. munge-sumstats could not determine sample size (N)' in message
    assert '\n  2. munge-sumstats could not find a signed summary statistic column' in message


def test_infer_only_does_not_turn_an_internal_build_error_into_a_note(tmp_path, monkeypatch):
    from ldsc._genome.build_policy import GenomeBuildCheck

    def broken(self):
        raise RuntimeError('internal failure')

    monkeypatch.setattr(GenomeBuildCheck, 'resolve', broken)
    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP CHR POS A1 A2 P BETA N\nrs1 1 100 A G 0.05 0.1 1000\n')
    with pytest.raises(RuntimeError, match='internal failure'):
        workflow.main(['--raw-sumstats-file', str(raw), '--snp-identifier', 'chr_pos', '--source-genome-build', 'hg19',
                       '--output-genome-build', 'hg19', '--no-snp-restriction',
                       '--output-dir', str(tmp_path / 'out'), '--infer-only'])


def test_loading_does_not_relabel_an_internal_footer_error(tmp_path, monkeypatch):
    from ldsc import load_sumstats
    from ldsc._sumstats import table

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n')
    workflow.main(['--raw-sumstats-file', str(raw), '--snp-identifier', 'rsid', '--no-snp-restriction',
                   '--output-dir', str(tmp_path / 'out')])

    def broken(*args, **kwargs):
        raise RuntimeError('internal failure')

    monkeypatch.setattr(table, 'validate_identity_artifact_metadata', broken)
    with pytest.raises(RuntimeError, match='internal failure'):
        load_sumstats(tmp_path / 'out' / 'sumstats.parquet')
