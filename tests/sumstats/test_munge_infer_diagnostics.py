"""Regression coverage for no-write munging diagnostics."""
import logging
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
    from ldsc import GlobalConfig, MungeConfig, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P BETA Z N\nrs1 A G 0.05 0.1 1.0 1000\n')
    report = infer_raw_sumstats(MungeConfig(raw_sumstats_file=raw), global_config=GlobalConfig(snp_identifier='rsid'))
    assert report.missing_fields == ('conflict: signed statistic',)


def test_run_and_python_inference_read_every_option_from_one_config(tmp_path):
    # Audit D9: the documented one-config preview ignored info_min, which the run rejects without an INFO column.
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP A1 A2 P BETA N\nrs1 A G 0.05 0.1 1000\n')
    config = MungeConfig(raw_sumstats_file=raw, info_min=0.5, no_snp_restriction=True, output_dir=tmp_path / 'out')
    rsid = GlobalConfig(snp_identifier='rsid')
    report = infer_raw_sumstats(config, global_config=rsid)
    assert report.missing_fields == ('column required by --info-min',)
    with pytest.raises(LDSCUserError) as caught:
        SumstatsMunger().run(config, global_config=rsid)
    assert str(caught.value) in report.notes


def test_run_infer_only_and_python_inference_share_one_build_plan(tmp_path, capsys):
    # Audit D9: for 50 chromosome-1 rows that no HM3 SNP supports, --infer-only reported the source build as
    # missing while infer_raw_sumstats() reported the file runnable; the run stops on the same build check.
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP CHR POS A1 A2 P BETA N\n'
                   + ''.join(f'rs{i} 1 {1000 * (i + 1)} A G 0.5 0.1 1000\n' for i in range(50)))
    cli = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'unused'),
                         '--snp-identifier', 'chr_pos', '--output-genome-build', 'hg38', '--infer-only'])
    config = MungeConfig(raw_sumstats_file=raw, output_genome_build='hg38', output_dir=tmp_path / 'out')
    chr_pos = GlobalConfig(snp_identifier='chr_pos')
    python = infer_raw_sumstats(config, global_config=chr_pos)
    assert not cli.runnable and not python.runnable
    assert python.missing_fields == cli.missing_fields == ('source_genome_build',)
    assert python.notes == cli.notes
    assert (python.source_genome_build, python.output_genome_build) == (None, 'hg38')
    assert (cli.source_genome_build, cli.output_genome_build) == (None, 'hg38')
    with pytest.raises(LDSCUserError, match='Genome-build consistency check failed') as caught:
        SumstatsMunger().run(config, global_config=chr_pos)
    assert any(str(caught.value) in note for note in python.notes)


def test_python_inference_suggests_the_options_of_its_config_as_infer_only_does(tmp_path):
    # Final review of audit D9: infer_raw_sumstats() left the request's own options (QC thresholds and others) out of
    # suggested_args, which --infer-only keeps, so a command rebuilt from the Python report ran with other QC.
    from ldsc import GlobalConfig, MungeConfig, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP CHR POS A1 A2 Z P N INFO FRQ\n'
                   + ''.join(f'rs{i} 1 {1000 * (i + 1)} A G 0.1 0.5 1000 0.99 0.3\n' for i in range(20)))
    options = {'--info-min': '0.95', '--maf-min': '0.02', '--N-min': '100.0', '--chunksize': '5',
               '--trait-name': 'T', '--output-format': 'both', '--info': 'INFO', '--log-level': 'WARNING'}
    cli = workflow.main(['--raw-sumstats-file', str(raw), '--output-dir', str(tmp_path / 'out'), '--infer-only',
                         '--snp-identifier', 'chr_pos', '--source-genome-build', 'hg19',
                         '--output-genome-build', 'hg19', '--no-snp-restriction',
                         *[token for option in options.items() for token in option]])
    config = MungeConfig(raw_sumstats_file=raw, output_dir=tmp_path / 'out', source_genome_build='hg19',
                         output_genome_build='hg19', no_snp_restriction=True, info_min=0.95, maf_min=0.02,
                         N_min=100.0, chunk_size=5, trait_name='T', output_format='both',
                         column_hints={'info': 'INFO'})
    python = infer_raw_sumstats(config, global_config=GlobalConfig(snp_identifier='chr_pos', log_level='WARNING'))
    assert cli.runnable, cli.missing_fields
    assert set(options.items()) <= set(zip(cli.suggested_args, cli.suggested_args[1:]))
    assert python == cli
    assert not (tmp_path / 'out').exists()


def _unresolvable_chr1_rows(n_column=True):
    """Return 50 chromosome-1 rows that no HM3 SNP supports, so ``source_genome_build='auto'`` cannot resolve."""
    header = 'SNP CHR POS A1 A2 P BETA' + (' N' if n_column else '')
    return header + '\n' + ''.join(
        f'rs{i} 1 {1000 * (i + 1)} A G 0.5 0.1' + (' 1000' if n_column else '') + '\n' for i in range(50))


def test_header_problems_stop_the_run_and_the_previews_before_the_build_check(tmp_path):
    # Review of audit D9: the build check ran before the header problems were raised, so a run whose header lacked N
    # and whose build could not be inferred scanned every coordinate row and logged the build ERROR before the N error.
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text(_unresolvable_chr1_rows(n_column=False))
    out = tmp_path / 'out'
    config = MungeConfig(raw_sumstats_file=raw, output_genome_build='hg19', no_snp_restriction=True, output_dir=out)
    chr_pos = GlobalConfig(snp_identifier='chr_pos')
    preview = infer_raw_sumstats(config, global_config=chr_pos)
    with pytest.raises(LDSCUserError, match='could not determine sample size') as caught:
        SumstatsMunger().run(config, global_config=chr_pos)
    assert preview.missing_fields == ('N',)
    assert str(caught.value) in preview.notes
    assert 'Genome-build consistency check failed' not in (out / 'diagnostics' / 'sumstats.log').read_text()


def test_the_run_raises_the_build_check_error_with_its_evidence(tmp_path):
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger

    raw = tmp_path / 'raw.tsv'
    raw.write_text(_unresolvable_chr1_rows())
    with pytest.raises(LDSCUserError, match='Genome-build consistency check failed') as caught:
        SumstatsMunger().run(MungeConfig(raw_sumstats_file=raw, output_genome_build='hg38', no_snp_restriction=True,
                                         output_dir=tmp_path / 'out'), global_config=GlobalConfig(snp_identifier='chr_pos'))
    # GenomeBuildCheck attaches its evidence to the error it raises; the run raises that error, not a copy.
    assert [entry.role for entry in caught.value.build_evidence] == ['--raw-sumstats-file']


@pytest.mark.parametrize('mode', ['rsid', 'chr_pos'])
def test_previews_report_a_keep_list_that_does_not_resolve(tmp_path, mode):
    # Review of audit D9: the plan did not resolve the keep-list, so in rsID modes the previews reported the file
    # runnable, and in chr_pos modes they blamed the source build, while the run stops on the missing keep-list.
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP CHR POS A1 A2 P BETA N\nrs1 1 100 A G 0.05 0.1 1000\n')
    builds = {} if mode == 'rsid' else {'source_genome_build': 'hg19', 'output_genome_build': 'hg19'}
    config = MungeConfig(raw_sumstats_file=raw, sumstats_snps_file=tmp_path / 'missing_keep.tsv',
                         output_dir=tmp_path / 'out', **builds)
    global_config = GlobalConfig(snp_identifier=mode)
    report = infer_raw_sumstats(config, global_config=global_config)
    assert report.missing_fields == ('sumstats_snps_file',)
    assert any('missing_keep.tsv' in note and 'summary-statistics SNPs file' in note for note in report.notes)
    with pytest.raises(LDSCUserError, match='missing_keep.tsv'):
        SumstatsMunger().run(config, global_config=global_config)


def test_previews_report_increasing_allele_columns_that_the_signed_statistic_contradicts(tmp_path):
    # Review of audit D9: only the run checked INC_ALLELE/DEC_ALLELE against the signed statistic.
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text('SNP INC_ALLELE DEC_ALLELE P BETA N\nrs1 A G 0.05 0.1 1000\nrs2 C T 0.05 -0.2 1000\n')
    config = MungeConfig(raw_sumstats_file=raw, no_snp_restriction=True, output_dir=tmp_path / 'out')
    rsid = GlobalConfig(snp_identifier='rsid')
    report = infer_raw_sumstats(config, global_config=rsid)
    assert report.missing_fields == ('A1/A2 orientation',)
    with pytest.raises(LDSCUserError, match='INC_ALLELE as A1') as caught:
        SumstatsMunger().run(config, global_config=rsid)
    assert str(caught.value) in report.notes


def test_the_python_preview_logs_no_warning_for_a_problem_it_reports(tmp_path, caplog):
    # Review of audit D9: the Python preview logged the grouped build-check ERROR, which reached stderr through
    # logging's last-resort handler, besides returning it in its notes.
    from ldsc import GlobalConfig, MungeConfig, infer_raw_sumstats

    raw = tmp_path / 'raw.tsv'
    raw.write_text(_unresolvable_chr1_rows())
    with caplog.at_level(logging.WARNING):
        report = infer_raw_sumstats(MungeConfig(raw_sumstats_file=raw, output_genome_build='hg38', no_snp_restriction=True),
                                    global_config=GlobalConfig(snp_identifier='chr_pos'))
    assert report.missing_fields == ('source_genome_build',)
    assert [record.getMessage() for record in caplog.records if record.levelno >= logging.WARNING] == []


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
