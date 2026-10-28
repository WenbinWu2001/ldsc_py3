"""User error guidance retains complete evidence and scientific failure policy."""

import logging

import pandas as pd
import pytest

from ldsc.errors import LDSCInputError


def test_saved_input_table_bounds_primary_error_without_losing_causes(tmp_path):
    from ldsc._runtime.preflight import InputGate

    table = tmp_path / "diagnostics/input_issues.tsv"
    gate = InputGate("Input validation", table)
    inputs = [tmp_path / f"missing-{i:03d}.tsv" for i in range(21)]
    for source in inputs:
        gate.check("baseline annotations", source, source.read_bytes)
    with pytest.raises(LDSCInputError) as failure:
        gate.finish()

    error = failure.value
    message = str(error)
    assert "21 input issues" in message
    assert str(inputs[0]) in message
    assert str(inputs[-1]) not in message
    assert str(table) in message
    saved = pd.read_csv(table, sep="\t", keep_default_na=False)
    assert saved.source.tolist() == [str(path) for path in inputs]
    assert len(error.input_issues) == 21
    assert isinstance(error.__cause__, ExceptionGroup)
    assert len(error.__cause__.exceptions) == 21


def test_without_written_input_table_primary_error_keeps_each_input(tmp_path):
    from ldsc._runtime.preflight import InputGate

    gate = InputGate("Input validation")
    inputs = [tmp_path / f"missing-{i:03d}.tsv" for i in range(21)]
    for source in inputs:
        gate.check("baseline annotations", source, source.read_bytes)
    with pytest.raises(LDSCInputError) as failure:
        gate.finish()
    assert all(str(path) in str(failure.value) for path in inputs)
    assert "Complete diagnostics:" not in str(failure.value)


def test_missing_saved_result_file_identifies_artifact_and_repair(tmp_path):
    from ldsc._runtime.result_files import declared_result_file

    with pytest.raises(LDSCInputError) as failure:
        declared_result_file(tmp_path, {"files": {"summary": "h2.tsv"}},
                             "summary", context="Saved h2 result")
    message = str(failure.value)
    assert "h2.tsv" in message and "Saved h2 result" in message
    assert "restore" in message.lower() and "same run" in message.lower()
    assert not list(tmp_path.iterdir())


def test_all_qc_rows_lost_reports_counts_and_current_options_without_sidecar(tmp_path):
    from ldsc.config import GlobalConfig, MungeConfig
    from ldsc.sumstats_munger import SumstatsMunger

    source = tmp_path / "raw.tsv"
    source.write_text("SNP A1 A2 P BETA N INFO\n"
                      "synthetic1 A C 0.05 0.1 1000 0.1\n"
                      "synthetic2 A G 0.10 0.2 1000 0.2\n")
    output = tmp_path / "result"
    config = MungeConfig(raw_sumstats_file=str(source), output_dir=str(output),
                         no_snp_restriction=True)
    with pytest.raises(LDSCInputError) as failure:
        SumstatsMunger().run(config, global_config=GlobalConfig(snp_identifier="rsid"))
    message = str(failure.value)
    assert "0 of 2" in message
    assert "--info-min" in message and "QC" in message
    assert "not individual rows" in message
    assert not (output / "diagnostics/dropped_snps/dropped.tsv.gz").exists()
    logfile = output / "diagnostics/sumstats.log"
    assert "Traceback" in logfile.read_text()
    assert message in logfile.read_text()


@pytest.mark.parametrize("reason,consequence", [
    ("zero_annotation_snps", "no SNP"),
    ("zero_variance_ld_scores", "constant"),
])
def test_default_query_warning_explains_skip_and_preserves_reason(caplog, reason, consequence):
    from ldsc.query_annotations import QueryAnnotationStatus, _log_query_annotation_statuses

    statuses = [QueryAnnotationStatus("query", "query.txt", "gene_list", "skipped", reason, 0.0)]
    with caplog.at_level(logging.INFO, logger="LDSC"):
        _log_query_annotation_statuses(statuses)
    record = next(record for record in caplog.records if record.levelno == logging.WARNING)
    message = record.getMessage()
    assert reason in message and "queries=1" in message
    assert consequence.lower() in message.lower()
    assert "omitted" in message.lower()
    assert "diagnostics/query_annotation_status.tsv" in message


def test_all_queries_skipped_error_explains_distinct_causes_and_repair():
    from ldsc.query_annotations import QueryAnnotationStatus, _all_query_annotations_skipped_message

    statuses = [
        QueryAnnotationStatus("empty_support", "query1.txt", "gene_list", "skipped", "zero_annotation_snps"),
        QueryAnnotationStatus("constant_scores", "query2.txt", "gene_list", "skipped", "zero_variance_ld_scores"),
    ]
    message = _all_query_annotations_skipped_message(statuses)
    assert "all 2" in message
    assert "empty_support" in message and "constant_scores" in message
    assert "no SNP" in message and "constant" in message
    assert "retained SNP support" in message
    assert "diagnostics/query_annotation_status.tsv" in message
