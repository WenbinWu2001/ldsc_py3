"""Python regression workflows save the same result directories as their commands."""

import json
from pathlib import Path
import re

import pandas as pd
import pytest

import ldsc
from ldsc import cli
from ldsc.errors import LDSCUsageError
from tests import test_regression_workflow as fixtures
from tests.test_regression_streaming import batch_inputs


def _inputs(tmp_path):
    """Write one LD-score directory and one munged trait; return their paths."""
    table, _ = batch_inputs(tmp_path / "ld")
    sumstats = tmp_path / "trait.parquet"
    fixtures.RegressionWorkflowTest().write_footer_sumstats_parquet(sumstats, frame=table.data, trait_name="MDD")
    return str(tmp_path / "ld"), str(sumstats)


def _table(path):
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def _metadata(result_dir):
    return json.loads((result_dir / "diagnostics" / "metadata.json").read_text())


# Required inputs of each workflow; validation fails before any of them is read.
_REQUIRED = {
    "run_h2": {"sumstats_file": "trait.parquet", "ldscore_dir": "ld"},
    "run_partitioned_h2": {"sumstats_file": "trait.parquet", "ldscore_dir": "ld"},
    "run_rg": {"sumstats_sources": ["a.parquet", "b.parquet"], "ldscore_dir": "ld"},
    "run_quantile_h2": {
        "partitioned_h2_result_dir": "fit", "target_annot_sources": ["target.tsv"], "target_annotation": "target",
        "baseline_annot_sources": ["baseline.tsv"], "ref_metadata_sources": ["ref.tsv"],
    },
}
_REGRESSION_CHOICE = ({"count_kind": "Common"}, "count_kind='Common' is not one of 'common', 'all'")
# Each case would otherwise run with a default or a truthy string in place of the intended value.
_INVALID_KEYWORDS = [
    *[(function, {"ld_score_dir": "ld"}, "unknown keyword argument(s) ld_score_dir") for function in _REQUIRED],
    *[(function, {name: "x"}, f"{name} is not a keyword argument")
      for function in _REQUIRED for name in ("log_level", "snp_identifier", "genome_build")],
    *[(function, {"overwrite": "False"}, "overwrite='False' must be False or True") for function in _REQUIRED],
    ("run_h2", *_REGRESSION_CHOICE),
    ("run_partitioned_h2", *_REGRESSION_CHOICE),
    ("run_rg", *_REGRESSION_CHOICE),
    ("run_partitioned_h2", {"continue_on_query_error": "no"}, "continue_on_query_error='no' must be False or True"),
    ("run_quantile_h2", {"gene_exclude_regions": "MHC"}, "gene_exclude_regions='MHC' is not one of 'none', 'mhc'"),
    ("run_quantile_h2", {"gene_list_resolution_policy": "resolved_only"},
     "gene_list_resolution_policy='resolved_only' must be 'strict' or 'resolved-only'"),
]


@pytest.mark.parametrize(("function", "extra", "message"), _INVALID_KEYWORDS)
def test_python_workflows_reject_invalid_keywords_before_writing(tmp_path, function, extra, message):
    output_dir = tmp_path / "out"
    with pytest.raises(LDSCUsageError, match=re.escape(f"{function}() cannot run: ") + ".*" + re.escape(message)):
        getattr(ldsc, function)(**_REQUIRED[function], output_dir=str(output_dir), **extra)
    assert not output_dir.exists()


@pytest.mark.parametrize("function", _REQUIRED)
def test_python_workflows_list_every_missing_required_keyword(function):
    with pytest.raises(LDSCUsageError, match=re.escape("missing required keyword argument(s)")) as raised:
        getattr(ldsc, function)()
    listed = str(raised.value).split("missing required keyword argument(s) ")[1].rstrip(".").split(", ")
    assert sorted(listed) == sorted([*_REQUIRED[function], "output_dir"])


def test_python_workflows_are_exported_and_their_namespace_entry_points_are_not():
    for function in _REQUIRED:
        assert function in ldsc.__all__
    assert not {f"{function}_from_args" for function in _REQUIRED} & set(ldsc.__all__)


def test_a_failed_python_overwrite_leaves_a_marker_naming_the_python_function(tmp_path):
    ldscore_dir, _ = _inputs(tmp_path)
    output_dir = tmp_path / "out"
    with pytest.raises(ldsc.LDSCInputError):
        ldsc.run_h2(sumstats_file=str(tmp_path / "missing.parquet"), ldscore_dir=ldscore_dir,
                    output_dir=str(output_dir), overwrite=True)
    assert "Command: Python API: run_h2(...)" in (output_dir / "RUN_FAILED.txt").read_text()


def test_run_h2_saves_the_result_the_h2_command_saves(tmp_path):
    ldscore_dir, sumstats = _inputs(tmp_path)
    cli.main(["h2", "--ldscore-dir", ldscore_dir, "--sumstats-file", sumstats,
              "--output-dir", str(tmp_path / "cli"), "--n-blocks", "6"])

    summary = ldsc.run_h2(ldscore_dir=ldscore_dir, sumstats_file=sumstats, output_dir=str(tmp_path / "py"), n_blocks=6)

    for name in ("h2.tsv", "diagnostics/ld_score_regression_bins.tsv"):
        pd.testing.assert_frame_equal(_table(tmp_path / "py" / name), _table(tmp_path / "cli" / name))
    assert _metadata(tmp_path / "py") == _metadata(tmp_path / "cli")
    assert summary.loc[0, "total_h2_obs"] == _table(tmp_path / "py" / "h2.tsv").loc[0, "total_h2_obs"]
    assert (tmp_path / "py" / "diagnostics" / "h2.log").is_file()


def test_a_python_h2_result_can_be_converted_and_plotted(tmp_path):
    ldscore_dir, sumstats = _inputs(tmp_path)
    result_dir = tmp_path / "h2"
    ldsc.run_h2(ldscore_dir=ldscore_dir, sumstats_file=sumstats, output_dir=str(result_dir), n_blocks=6)

    conversion = ldsc.convert_h2_scale(result_dir, samp_prev=0.5, pop_prev=0.01)
    plot = ldsc.plot_result(result_dir)

    assert Path(conversion.table_path).resolve().is_relative_to((result_dir / "postprocessing" / "liability-scale").resolve())
    assert plot.path.resolve().is_relative_to((result_dir / "plots").resolve())


def test_run_partitioned_h2_saves_the_command_result_and_returns_every_query_status(tmp_path):
    ldscore_dir, sumstats = _inputs(tmp_path)
    cli.main(["partitioned-h2", "--ldscore-dir", ldscore_dir, "--sumstats-file", sumstats,
              "--output-dir", str(tmp_path / "cli"), "--n-blocks", "6", "--query-batch-size", "2"])

    result = ldsc.run_partitioned_h2(
        ldscore_dir=ldscore_dir, sumstats_file=sumstats, output_dir=str(tmp_path / "py"), n_blocks=6, query_batch_size=2,
    )

    pd.testing.assert_frame_equal(_table(tmp_path / "py" / "partitioned_h2.tsv"), _table(tmp_path / "cli" / "partitioned_h2.tsv"))
    assert isinstance(result, ldsc.PartitionedH2BatchResult)
    assert result.query_status["query_annotation"].tolist() == ["first", "second", "third"]
    assert set(result.query_status["status"]) == {"success"}
    assert Path(result.output_paths["summary"]).resolve() == (tmp_path / "py" / "partitioned_h2.tsv").resolve()


def _rg_inputs(tmp_path):
    """Write one LD-score directory and two munged traits with alleles; return their paths."""
    table, _ = batch_inputs(tmp_path / "ld")
    paths = [tmp_path / "SCZ.parquet", tmp_path / "BMI.parquet"]
    for path, z_scale in zip(paths, (1.0, 0.5)):
        frame = table.data.assign(A1="A", A2="G", Z=table.data.Z * z_scale)
        fixtures.RegressionWorkflowTest().write_footer_sumstats_parquet(path, frame=frame, trait_name=path.stem)
    return str(tmp_path / "ld"), [str(path) for path in paths]


def test_run_rg_takes_prevalence_lists_and_saves_the_command_result(tmp_path):
    ldscore_dir, sumstats = _rg_inputs(tmp_path)
    cli.main(["rg", "--ldscore-dir", ldscore_dir, "--sumstats-sources", *sumstats, "--output-dir", str(tmp_path / "cli"),
              "--n-blocks", "6", "--samp-prev", "0.5,nan", "--pop-prev", "0.01,nan"])

    result = ldsc.run_rg(
        ldscore_dir=ldscore_dir, sumstats_sources=sumstats, output_dir=str(tmp_path / "py"), n_blocks=6,
        samp_prev=[0.5, float("nan")], pop_prev=[0.01, float("nan")],
    )

    for name in ("rg.tsv", "rg_full.tsv", "h2_per_trait.tsv"):
        pd.testing.assert_frame_equal(_table(tmp_path / "py" / name), _table(tmp_path / "cli" / name))
    assert _metadata(tmp_path / "py") == _metadata(tmp_path / "cli")
    assert isinstance(result, ldsc.RgResultFamily)
    assert result.h2_per_trait["trait_name"].tolist() == ["SCZ", "BMI"]
    assert result.h2_per_trait["samp_prev"].iloc[0] == 0.5


def test_run_rg_rejects_a_single_prevalence_number(tmp_path):
    with pytest.raises(LDSCUsageError, match="one value per trait"):
        ldsc.run_rg(ldscore_dir="ld", sumstats_sources=["a", "b"], output_dir=str(tmp_path / "out"),
                    samp_prev=0.5, pop_prev=0.01)
    assert not (tmp_path / "out").exists()
