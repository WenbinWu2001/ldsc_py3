"""Batch logs expose aggregate progress while diagnostic artifacts retain detail."""

import logging
from pathlib import Path

import pandas as pd
import pytest

from ldsc import cli
from tests.regression.test_regression_streaming import cli_inputs
from tests.ldscore.test_gene_query_integration import inputs, direct_args, build_index


@pytest.mark.parametrize("workers", [1, 2])
@pytest.mark.parametrize("level", ["INFO", "DEBUG"])
def test_partitioned_log_has_numbered_batches_and_aggregate_outcomes(tmp_path, workers, level):
    args, _ = cli_inputs(tmp_path)
    result = cli.main([*args, "--workers", str(workers), "--query-batch-size", "2", "--log-level", level]).summary
    log = (tmp_path / "out/diagnostics/partitioned-h2.log").read_text()
    for ordinal in (1, 2):
        assert f"Phase start: batch; query fitting {ordinal}/2" in log
        assert f"Phase complete: batch; query fitting {ordinal}/2" in log
    assert "completed=3/3; success=3" in log
    assert ("Phase start: query model preparation." in log) == (level == "DEBUG")
    assert ("Phase start: query output staging." in log) == (level == "DEBUG")
    status = pd.read_csv(tmp_path / "out/diagnostics/query_status.tsv", sep="\t")
    assert status.status.tolist() == ["success"] * 3
    assert set(result.category) == {"first", "second", "third"}


def test_progress_counts_outcomes_and_retains_object_only_at_debug(caplog, monkeypatch):
    from ldsc._runtime import progress as _progress

    now = [0.0]
    monkeypatch.setattr(_progress, "perf_counter", lambda: now[0])
    with caplog.at_level(logging.INFO, logger="LDSC"):
        with _progress.PhaseProgress(logging.getLogger("LDSC.test"), "computation", "queries", 3) as progress:
            progress.advance(outcome="success", object="long-query-name")
            now[0] = 31
            progress.advance(outcome="unestimable", object="another-query-name")
            progress.advance(outcome="success")
    records = [r for r in caplog.records if hasattr(r, "phase_event")]
    assert [r.phase_event for r in records] == ["start", "progress", "complete"]
    assert "completed=2/3; success=1; unestimable=1" in records[1].getMessage()
    assert "completed=3/3; success=2; unestimable=1" in records[-1].getMessage()
    assert records[1].current_object == "another-query-name"
    assert "query-name" not in caplog.text


@pytest.mark.parametrize("indexed", [False, True])
def test_ldscore_batches_log_boundaries_and_one_global_status(tmp_path, inputs, indexed):
    from ldsc._runtime.workflow_logging import workflow_logging
    from ldsc._ldscore.calculator import run_ldscore_from_args
    from ldsc._ldscore.gene_index.scoring import run_indexed_ldscore

    sources = [tmp_path / f"{name}.txt" for name in ("first", "second", "third")]
    for source, gene in zip(sources, ("G1", "G2", "G1")):
        source.write_text(gene + "\n")
    if indexed:
        index = build_index(tmp_path, inputs)
        output = tmp_path / "indexed"
        with workflow_logging("ldscore", output / "diagnostics/ldscore.log"):
            result = run_indexed_ldscore(index, query_gene_list_sources=sources, query_batch_size=2,
                                        output_dir=output)
        label = "indexed LD-score assembly"
    else:
        output = tmp_path / "direct"
        args = direct_args(tmp_path, inputs, sources)
        args.query_batch_size = 2
        result = run_ldscore_from_args(args)
        label = "LD-score calculation"
    log = (output / "diagnostics/ldscore.log").read_text()
    for ordinal in (1, 2):
        assert f"Phase start: batch; {label} {ordinal}/2" in log
        assert f"Phase complete: batch; {label} {ordinal}/2" in log
    assert log.count("Query annotation status summary:") == 1
    assert "total=3, ok=3, warning=0, skipped=0" in log
    assert "binary=base," not in log
    assert result.query_columns == ["first", "second", "third"]


@pytest.mark.parametrize("n_queries", [3, 13])
def test_constant_gene_query_notice_explains_skips_without_claiming_zero_snp_support(
    tmp_path, capsys, n_queries
):
    from types import SimpleNamespace
    from ldsc._runtime.workflow_logging import workflow_logging
    from ldsc._annotation.query_outcomes import QueryAnnotationStatus, emit_gene_gate_b_notice

    batch = SimpleNamespace(audit_frames=lambda: iter([
        pd.DataFrame({"disposition": ["selected"]})
    ]))
    statuses = tuple(
        QueryAnnotationStatus(f"query-{i}", f"query-{i}.txt", "gene_list", "skipped",
                              "zero_variance_ld_scores", 4.0)
        for i in range(n_queries)
    )
    logfile = tmp_path / "ldscore.log"
    with workflow_logging("ldscore", logfile):
        emit_gene_gate_b_notice(statuses, batch, command="ldscore")

    console = capsys.readouterr()
    notice = console.err
    assert console.out == ""
    assert f"{n_queries} query annotation(s) skipped" in notice
    assert "LD scores are constant on regression SNP rows (zero_variance_ld_scores)" in notice
    assert "Skipped queries are omitted from LD-score outputs." in notice
    assert "zero retained-SNP support for 0" not in notice
    assert "preflight" not in notice
    assert "diagnostics/query_annotation_status.tsv" in notice
    assert notice.removeprefix("WARNING: ").strip() in logfile.read_text()
    if n_queries == 13:
        assert "query-9=zero_variance_ld_scores" in notice
        assert "query-10=" not in notice
        assert "plus 3 more affected query(s)" in notice


def test_many_annotation_and_bed_sources_have_collection_progress(tmp_path, caplog):
    from ldsc._annotation.sources import prepare_annotation_sources
    from ldsc._annotation.interval_queries import prepare_bed_queries
    from ldsc._annotation.storage import AnnotationWorkspace

    paths = [tmp_path / f"source-{index}.annot" for index in range(40)]
    beds = [tmp_path / f"source-{index}.bed" for index in range(40)]
    for index, (path, bed) in enumerate(zip(paths, beds)):
        path.write_text(f"CHR SNP BP CM MAF col{index}\n1 rs1 10 0 .2 1\n2 rs2 20 0 .3 0\n")
        bed.write_text("1\t0\t10\n")
    with caplog.at_level(logging.INFO, logger="LDSC"), AnnotationWorkspace(tmp_path / "out") as workspace:
        prepared = prepare_annotation_sources(workspace, paths[:1], paths[1:], mode="rsid")
        sources = prepare_bed_queries(beds, workspace)
        assert prepared.shards["1"].read(columns=[f"col{index}" for index in range(40)]).shape == (1, 40)
        assert sum(source.n_intervals for source in sources) == 40
    messages = [record.getMessage() for record in caplog.records if record.levelno == logging.INFO]
    assert len(messages) < 25
    assert any("annotation input files; completed=40/40" in message for message in messages)
    assert not any("source-" in message for message in messages)


def test_munging_warning_counts_do_not_repeat_per_chunk_or_leak_between_runs(tmp_path):
    from ldsc import GlobalConfig, MungeConfig, SumstatsMunger

    raw = tmp_path / "raw.tsv"
    raw.write_text("SNP P BETA N INFO FRQ\n" + "".join(
        f"bad{i} 0 0 1000 3 2\nkeep{i} .05 .05 1000 1 .2\n" for i in range(12)))
    for chunk_size in (2, 100):
        output = tmp_path / f"out{chunk_size}"
        munger = SumstatsMunger()
        table = munger.run(MungeConfig(raw_sumstats_file=raw, output_dir=output, no_snp_restriction=True,
                                      chunk_size=chunk_size), global_config=GlobalConfig(snp_identifier="rsid"))
        log = (output / "diagnostics/sumstats.log").read_text()
        assert log.count("Invalid INFO values") == 1
        assert log.count("Invalid FRQ values") == 1
        assert log.count("Invalid P values") == 1
        assert "Invalid-value totals observed: FRQ=12, INFO=12, P=12" in log
        assert "[0,2]" in log
        assert table.data.SNP.tolist() == [f"keep{i}" for i in range(12)]
        assert munger.build_run_summary().drop_counts["FRQ"] == 12


def test_rg_reports_trait_and_pair_totals_without_repeating_input_provenance(tmp_path, caplog):
    from dataclasses import replace
    from ldsc import load_sumstats, RegressionRunner, RegressionConfig

    _, source = cli_inputs(tmp_path)
    table = load_sumstats(tmp_path / "trait.parquet")
    # rg aligns signed statistics by allele, so each trait needs A1/A2.
    table = replace(table, data=table.data.assign(A1="A", A2="G"), has_alleles=True)
    traits = [replace(table, trait_name=f"trait-{i}") for i in range(3)]
    with caplog.at_level(logging.INFO, logger="LDSC"):
        result = RegressionRunner(regression_config=RegressionConfig(num_blocks=6, use_intercept=False)).estimate_rg_pairs(traits, source)
    assert len(result.h2_per_trait) == len(result.rg_full) == 3
    assert "Phase complete: computation; per-trait h2; completed=3/3" in caplog.text
    assert "Phase complete: computation; rg pairs; completed=3/3" in caplog.text
    assert "Regression genome build for sumstats" not in caplog.text


def test_quantile_log_points_to_bin_table():
    from tests.regression.test_quantile_h2 import QuantileH2WorkflowTest
    from ldsc._regression.quantile_h2 import run_quantile_h2_from_args

    case = QuantileH2WorkflowTest()
    case.setUp()
    try:
        result = run_quantile_h2_from_args(case.args)
        text = (Path(case.args.output_dir) / "diagnostics/quantile-h2.log").read_text()
        assert "Quantile summary: bins=2;" in text and "See quantile_h2.tsv" in text
        assert "Quantile 1:" not in text
        assert result.quantile_h2.n_snps.tolist() == [3, 1]
    finally:
        case.doCleanups()


@pytest.mark.parametrize("invalid_query", [False, True])
def test_unused_catalog_defects_are_saved_with_one_warning(tmp_path, inputs, invalid_query):
    from ldsc.errors import LDSCInputError
    from ldsc._ldscore.calculator import run_ldscore_from_args
    import json

    with inputs[3].open("a") as stream:
        stream.writelines(f"BAD{i}\tBAD{i}\t1\tinvalid\t20\thg19\n" for i in range(20))
    focal = tmp_path / "focal.txt"
    focal.write_text("UNKNOWN\n" if invalid_query else "G1\n")
    args = direct_args(tmp_path, inputs, [focal])
    if invalid_query:
        with pytest.raises(LDSCInputError):
            run_ldscore_from_args(args)
    else:
        result = run_ldscore_from_args(args)
        assert result.query_columns == ["focal"]
        metadata = json.loads((tmp_path / "direct/metadata.json").read_text())
        assert metadata["query_diagnostics"]["gene_catalog_issues"] == "diagnostics/gene_catalog_issues.tsv"
    issues = pd.read_csv(tmp_path / "direct/diagnostics/gene_catalog_issues.tsv", sep="\t")
    assert len(issues) == 20
    log = (tmp_path / "direct/diagnostics/ldscore.log").read_text()
    assert log.count("Gene coordinate catalog defects:") == 1
    assert "records=20" in log and "See diagnostics/gene_catalog_issues.tsv" in log
    assert "BAD19" not in log
    if not invalid_query:
        inputs[3].write_text("\n".join(line for line in inputs[3].read_text().splitlines() if not line.startswith("BAD")) + "\n")
        args.overwrite = True
        cleaned = run_ldscore_from_args(args)
        pd.testing.assert_frame_equal(result.baseline_table, cleaned.baseline_table)
        assert not (tmp_path / "direct/diagnostics/gene_catalog_issues.tsv").exists()


def test_output_audit_summarizes_query_files_with_manifest_pointer(tmp_path):
    from ldsc._runtime.workflow_logging import workflow_logging, log_outputs

    logfile = tmp_path / "run.log"
    paths = {f"query_batch{i:05d}": str(tmp_path / f"ldscore.query.batch{i:05d}.parquet") for i in range(1, 101)}
    paths.update(metadata=str(tmp_path / "metadata.json"), query_status=str(tmp_path / "diagnostics/query_annotation_status.tsv"))
    with workflow_logging("ldscore", logfile):
        log_outputs(**paths)
    text = logfile.read_text()
    assert "query_batches" in text and "100 files; see " + paths["metadata"] in text
    assert paths["query_status"] in text
    assert "query_batch00100" not in text
    assert len(paths) == 102


def test_rg_duplicate_trait_names_are_summarized(tmp_path, caplog):
    from dataclasses import replace
    from ldsc import load_sumstats
    from ldsc._regression.runner import _disambiguate_trait_names

    cli_inputs(tmp_path)
    table = load_sumstats(tmp_path / "trait.parquet")
    tables = [replace(table, trait_name="same", source_path=str(tmp_path / str(i) / "trait.parquet")) for i in range(40)]
    with caplog.at_level(logging.INFO, logger="LDSC"):
        renamed = _disambiguate_trait_names(tables)
    assert len({table.trait_name for table in renamed}) == 40
    assert len(caplog.records) == 1
    assert "Renamed 40 duplicate rg trait labels" in caplog.text
    assert "diagnostics/metadata.json" in caplog.text


def _rg_tables(tmp_path, named_paths):
    from dataclasses import replace
    from ldsc import load_sumstats

    if not (tmp_path / "trait.parquet").exists():
        cli_inputs(tmp_path)
    table = load_sumstats(tmp_path / "trait.parquet")
    return [replace(table, trait_name=name, source_path=str(tmp_path / path)) for name, path in named_paths]


def test_rg_label_collisions_add_path_parts_instead_of_hashes(tmp_path):
    from ldsc._regression.runner import _disambiguate_trait_names

    named_paths = [
        ("SCZ", "pgc/v3/trait.parquet"), ("SCZ", "ukb/v3/trait.parquet"), ("BIP", "pgc/v3/bip.parquet"),
        ("ADHD", "same/adhd1.parquet"), ("ADHD", "same/adhd2.parquet"),
        ("MDD", "a/x/mdd1.parquet"), ("MDD", "a/x/mdd2.parquet"), ("MDD", "b/x/mdd1.parquet"),
    ]
    # Parent directories are tried first; same-directory files fall back to the
    # shortest distinguishing path that ends in the file name.
    expected = ["SCZ@pgc/v3", "SCZ@ukb/v3", "BIP", "ADHD@adhd1.parquet", "ADHD@adhd2.parquet",
                "MDD@a/x/mdd1.parquet", "MDD@a/x/mdd2.parquet", "MDD@b/x/mdd1.parquet"]
    labels = [table.trait_name for table in _disambiguate_trait_names(_rg_tables(tmp_path, named_paths))]
    reversed_labels = [table.trait_name for table in
                       _disambiguate_trait_names(_rg_tables(tmp_path, named_paths[::-1]))]
    assert labels == expected
    assert reversed_labels == expected[::-1]


def test_rg_rejects_trait_name_that_matches_a_tie_breaker_label(tmp_path):
    from ldsc.errors import LDSCInputError
    from ldsc._regression.runner import _disambiguate_trait_names

    tables = _rg_tables(tmp_path, [("SCZ", "pgc/trait.parquet"), ("SCZ", "ukb/trait.parquet"),
                                   ("SCZ@pgc", "other/trait.parquet")])
    with pytest.raises(LDSCInputError, match=r"'SCZ@pgc'.*--trait-name"):
        _disambiguate_trait_names(tables)
