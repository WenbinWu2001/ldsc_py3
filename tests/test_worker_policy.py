"""Shared chromosome/query-worker validation and execution contracts."""

import json

import pytest

from ldsc.config import LDScoreConfig, GeneLDScoreIndexBuildConfig
from ldsc.errors import LDSCConfigError
from ldsc.gene_ldscore_index import run_indexed_ldscore


@pytest.mark.parametrize("threads", [0, True, False, 1.5, 1.0, "2", None])
@pytest.mark.parametrize("workflow", ["direct", "build", "indexed", "partitioned"])
def test_invalid_workers_rejected_before_io(tmp_path, threads, workflow):
    error = ValueError if workflow in {"indexed", "partitioned"} else LDSCConfigError
    with pytest.raises(error, match="threads"):
        if workflow == "direct":
            LDScoreConfig(ld_wind_cm=1, threads=threads)
        elif workflow == "build":
            GeneLDScoreIndexBuildConfig(
                baseline_annot_sources=("missing",), plink_prefix="missing",
                gene_coordinate_file="missing", output_dir=str(tmp_path / "out"),
                genome_build="hg19", snp_identifier="rsid", threads=threads,
            )
        elif workflow == "indexed":
            run_indexed_ldscore("missing", query_gene_list_sources=[],
                                output_dir=tmp_path / "out", threads=threads)
        else:
            from ldsc import RegressionRunner
            RegressionRunner().estimate_partitioned_h2_batch(
                None, None, output_dir=tmp_path / "out", threads=threads,
            )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("module", ["ldscore_calculator", "gene_ldscore_index"])
@pytest.mark.parametrize("value", ["0", "True", "False", "1.5", "1.0"])
def test_cli_rejects_invalid_workers(module, value, capsys):
    from importlib import import_module
    parser = import_module(f"ldsc.{module}").build_parser()
    with pytest.raises(SystemExit) as error:
        parser.parse_args(["--output-dir", "out", "--threads", value])
    assert error.value.code == 2
    assert "argument --threads:" in capsys.readouterr().err


@pytest.mark.parametrize("cpus,threads,expected", [
    (1, 1, 1), (1, 4, 2), (1, -1, 1), (1, -2, 1),
    (3, -1, 2), (3, -2, 2), (3, -3, 1), (3, -20, 1),
])
def test_workflows_share_affinity_resolution_and_inline_execution(tmp_path, monkeypatch, cpus, threads, expected):
    import os
    import numpy as np
    from ldsc import gene_ldscore_index as builder, ldscore_calculator as direct
    from ldsc import _indexed_ldscore_batches as indexed
    from ldsc import RegressionRunner, RegressionConfig
    from tests.test_plink_workflow_resolution import write_inputs, index_args
    from tests.test_regression_streaming import batch_inputs

    monkeypatch.setattr(os, "cpu_count", lambda: 8)
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(range(cpus)), raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    calls = {"build": [], "direct": [], "indexed": []}

    def observe(module, attribute, workflow):
        executor = getattr(module, attribute)

        def construct(*args, **kwargs):
            assert expected > 1, f"{workflow} constructed an executor for one effective worker"
            calls[workflow].append(kwargs["max_workers"])
            return executor(*args, **kwargs)

        monkeypatch.setattr(module, attribute, construct)

    observe(builder, "ThreadPoolExecutor", "build")
    observe(direct, "ProcessPoolExecutor", "direct")
    observe(indexed, "ProcessPoolExecutor", "indexed")
    write_inputs(tmp_path)
    token = str(tmp_path / "1000G.EUR.QC.*")
    args = index_args(tmp_path, token, tmp_path / "index")
    args.threads = threads
    path = builder.run_build_gene_ldscore_index_from_args(args)
    index = builder._load_gene_ldscore_index(path, _allow_partial_for_tests=True)
    np.testing.assert_array_equal(index.gene_support, [4, 4])

    args = direct.build_parser().parse_args([
        "--baseline-annot-sources", str(tmp_path / "baseline.annot"),
        "--plink-prefix", token, "--snp-identifier", "chr_pos", "--genome-build", "hg19",
        "--ld-wind-cm", ".05", "--snp-batch-size", "1",
        "--regr-snps-exclude-regions", "none", "--regr-snps-file", str(tmp_path / "regression.tsv"),
        "--output-dir", str(tmp_path / "direct"), "--threads", str(threads),
    ])
    result = direct.run_ldscore_from_args(args)
    np.testing.assert_allclose(result.baseline_table["base"], np.ones(8), atol=1e-7)

    query = tmp_path / "query.txt"
    query.write_text("G21\n")
    result = builder.run_indexed_ldscore(path, query_gene_list_sources=[query], threads=threads,
                                        output_dir=tmp_path / "indexed", _allow_partial_for_tests=True)
    np.testing.assert_allclose(result.baseline_table["base"], np.ones(8), atol=1e-7)
    np.testing.assert_allclose(result.read_queries(["query"])["query"], [1.] * 4 + [0.] * 4, atol=1e-7)
    assert calls == {name: ([] if expected == 1 else [expected]) for name in calls}

    table, source = batch_inputs(tmp_path / "query-inputs")
    result = RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6)).estimate_partitioned_h2_batch(
        table, source, output_dir=tmp_path / "partitioned", threads=threads,
        query_columns=["first", "second"],
    )
    metadata = json.loads((tmp_path / "partitioned/diagnostics/metadata.json").read_text())
    assert metadata["query_workers_requested"] == threads
    assert metadata["query_workers_effective"] == expected
    assert result.query_status.status.tolist() == ["success", "success"]


@pytest.mark.parametrize("affinity", ["missing", "unavailable", "empty"])
@pytest.mark.parametrize("cpu_count,expected", [(4, 3), (None, 1)])
def test_worker_cpu_discovery_fallback(monkeypatch, affinity, cpu_count, expected):
    from ldsc import _parallelism as policy
    monkeypatch.setattr(policy.os, "cpu_count", lambda: cpu_count)
    if affinity == "missing":
        monkeypatch.delattr(policy.os, "sched_getaffinity", raising=False)
    else:
        def getaffinity(pid):
            if affinity == "unavailable":
                raise OSError("affinity unavailable")
            return set()
        monkeypatch.setattr(policy.os, "sched_getaffinity", getaffinity, raising=False)
    assert policy._resolve_worker_count(-1, 3) == (1 if affinity == "empty" else expected)


@pytest.mark.parametrize("threads", [0, True, False, 1.0, 1.5, "2", None])
def test_direct_python_wrapper_rejects_invalid_workers(tmp_path, threads):
    from ldsc import run_ldscore
    with pytest.raises(LDSCConfigError, match="threads"):
        run_ldscore(r2_dir="missing", ld_wind_cm=1, output_dir=tmp_path / "out", threads=threads)
    assert not (tmp_path / "out").exists()
