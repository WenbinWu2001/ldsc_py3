"""Shared chromosome/query-worker validation and execution contracts."""

import json

import pytest

from ldsc.config import LDScoreConfig, GeneLDScoreIndexBuildConfig
from ldsc.errors import LDSCConfigError
from ldsc._ldscore.gene_index.scoring import run_indexed_ldscore


@pytest.mark.parametrize("workers", [0, True, False, 1.5, 1.0, "2", None])
@pytest.mark.parametrize("workflow", ["direct", "build", "indexed", "partitioned"])
def test_invalid_workers_rejected_before_io(tmp_path, workers, workflow):
    error = ValueError if workflow in {"indexed", "partitioned"} else LDSCConfigError
    with pytest.raises(error, match="workers"):
        if workflow == "direct":
            LDScoreConfig(ld_wind_cm=1, workers=workers)
        elif workflow == "build":
            GeneLDScoreIndexBuildConfig(ld_wind_cm=1.0,
                baseline_annot_sources=("missing",), plink_prefix="missing",
                gene_coordinate_file="missing", output_dir=str(tmp_path / "out"),
                genome_build="hg19", snp_identifier="rsid", workers=workers,
            )
        elif workflow == "indexed":
            run_indexed_ldscore("missing", query_gene_list_sources=[],
                                output_dir=tmp_path / "out", workers=workers)
        else:
            from ldsc import RegressionRunner
            RegressionRunner().estimate_partitioned_h2_batch(
                None, None, output_dir=tmp_path / "out", workers=workers,
            )
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("module", ["_ldscore.calculator", "_ldscore.gene_index.build"])
@pytest.mark.parametrize("value", ["0", "True", "False", "1.5", "1.0"])
def test_cli_rejects_invalid_workers(module, value, capsys):
    from importlib import import_module
    parser = import_module(f"ldsc.{module}").build_parser()
    with pytest.raises(SystemExit) as error:
        parser.parse_args(["--output-dir", "out", "--workers", value])
    assert error.value.code == 2
    assert "argument --workers:" in capsys.readouterr().err


@pytest.mark.parametrize("cpus,workers,expected", [
    (1, 1, 1), (4, 4, 2), (1, -1, 1), (1, -2, 1),
    (3, -1, 2), (3, -2, 2), (3, -3, 1), (3, -20, 1),
])
def test_workflows_share_detected_budget_and_managed_processes(tmp_path, monkeypatch, cpus, workers, expected):
    import os
    import numpy as np
    from ldsc._ldscore.gene_index import build as builder
    from ldsc._ldscore import calculator as direct
    from ldsc._ldscore.gene_index import scoring as _gene_index_scoring
    from ldsc._ldscore.gene_index.scoring import _load_gene_ldscore_index
    from ldsc._ldscore.gene_index import batches as indexed
    from ldsc import RegressionRunner, RegressionConfig
    from tests.ldscore.test_plink_workflow_resolution import write_inputs, index_args
    from tests.regression.test_regression_streaming import batch_inputs

    monkeypatch.setattr("joblib.cpu_count", lambda: cpus)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    calls = {"build": [], "direct": [], "indexed": []}

    def observe(module, attribute, workflow):
        executor = getattr(module, attribute)

        def construct(*args, **kwargs):
            calls[workflow].append(kwargs["max_workers"])
            return executor(*args, **kwargs)

        monkeypatch.setattr(module, attribute, construct)

    observe(builder, "ProcessPoolExecutor", "build")
    observe(direct, "ProcessPoolExecutor", "direct")
    observe(indexed, "ProcessPoolExecutor", "indexed")
    write_inputs(tmp_path)
    token = str(tmp_path / "1000G.EUR.QC.*")
    args = index_args(tmp_path, token, tmp_path / "index", workers=workers)
    args.workers = workers
    path = builder.run_build_gene_ldscore_index_from_args(args)
    index = _load_gene_ldscore_index(path, _allow_partial_for_tests=True)
    np.testing.assert_array_equal(index.gene_support, [4, 4])

    args = direct.build_parser().parse_args([
        "--baseline-annot-sources", str(tmp_path / "baseline.annot"),
        "--plink-prefix", token, "--snp-identifier", "chr_pos", "--genome-build", "hg19",
        "--ld-wind-cm", ".05", "--snp-batch-size", "1",
        "--regr-snps-exclude-regions", "none", "--regr-snps-file", str(tmp_path / "regression.tsv"),
        "--output-dir", str(tmp_path / "direct"), "--workers", str(workers),
    ])
    result = direct.run_ldscore_from_args(args)
    np.testing.assert_allclose(result.baseline_table["base"], np.ones(8), atol=1e-7)

    query = tmp_path / "query.txt"
    query.write_text("G21\n")
    result = _gene_index_scoring.run_indexed_ldscore(path, query_gene_list_sources=[query], workers=workers,
                                        output_dir=tmp_path / "indexed", _allow_partial_for_tests=True)
    np.testing.assert_allclose(result.baseline_table["base"], np.ones(8), atol=1e-7)
    np.testing.assert_allclose(result.read_queries(["query"])["query"], [1.] * 4 + [0.] * 4, atol=1e-7)
    assert calls == {name: [expected] for name in calls}

    table, source = batch_inputs(tmp_path / "query-inputs")
    result = RegressionRunner(source.config_snapshot, RegressionConfig(n_blocks=6)).estimate_partitioned_h2_batch(
        table, source, output_dir=tmp_path / "partitioned", workers=workers,
        query_columns=["first", "second"],
    )
    metadata = json.loads((tmp_path / "partitioned/diagnostics/metadata.json").read_text())
    assert metadata["query_workers_requested"] == workers
    assert metadata["query_workers_effective"] == expected
    assert result.query_status.status.tolist() == ["success", "success"]


@pytest.mark.parametrize("cpu_count,expected", [(4, 3), (1, 1)])
def test_worker_cpu_discovery_delegates_to_joblib(monkeypatch, cpu_count, expected):
    from ldsc._runtime import parallelism as policy
    monkeypatch.setattr(policy.joblib, "cpu_count", lambda: cpu_count)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "1")
    assert policy.resolve_worker_count(-1, 3) == expected


@pytest.mark.parametrize("workers", [0, True, False, 1.0, 1.5, "2", None])
def test_direct_python_wrapper_rejects_invalid_workers(tmp_path, workers):
    from ldsc import run_ldscore
    with pytest.raises(LDSCConfigError, match="workers"):
        run_ldscore(r2_dir="missing", ld_wind_cm=1, output_dir=tmp_path / "out", workers=workers)
    assert not (tmp_path / "out").exists()
