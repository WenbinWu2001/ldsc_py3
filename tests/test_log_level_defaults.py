"""Workflow functions that take ``log_level`` default to the registered GlobalConfig level."""

from pathlib import Path

import pytest

import ldsc
from ldsc._ldscore.gene_index import build as gene_ldscore_index
from ldsc._regression import h2_scale
from ldsc._plotting import dispatch as plotting
from ldsc._ldscore.legacy_converter import LegacyLDScoreConverter
from tests.test_python_regression_workflows import _inputs


class _Requested(Exception):
    """Stop a workflow once it has requested its log level."""


@pytest.fixture(scope="module")
def h2_result(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("h2")
    ldscore_dir, sumstats = _inputs(tmp_path)
    ldsc.run_h2(ldscore_dir=ldscore_dir, sumstats_file=sumstats, output_dir=str(tmp_path / "h2"), n_blocks=6)
    return tmp_path / "h2"


@pytest.fixture
def debug_session():
    ldsc.set_global_config(ldsc.GlobalConfig(log_level="DEBUG"))
    yield
    ldsc.reset_global_config()


def _stop_with_level(monkeypatch, owner, name, level_of):
    """Replace ``owner.name`` with a stub that raises the level it was asked for."""
    def stub(*args, **kwargs):
        raise _Requested(level_of(*args, **kwargs))
    monkeypatch.setattr(owner, name, stub)


def _plot(h2_result, monkeypatch, **level):
    # Plotting and conversion open their workflow log in this process.
    _stop_with_level(monkeypatch, plotting, "workflow_logging", lambda *args, log_level: log_level)
    ldsc.plot_result(h2_result, **level)


def _convert_scale(h2_result, monkeypatch, **level):
    _stop_with_level(monkeypatch, h2_scale, "workflow_logging", lambda *args, log_level: log_level)
    ldsc.convert_h2_scale(h2_result, samp_prev=0.5, pop_prev=0.01, **level)


def _convert_ldscores(h2_result, monkeypatch, **level):
    # The conversion body runs in a spawned child with its own default config,
    # so the level must already be resolved when the parent hands it over.
    _stop_with_level(monkeypatch, LegacyLDScoreConverter, "_convert", lambda self, **kwargs: kwargs["log_level"])
    ldsc.convert_ldsc2_ldscores(legacy_reference_dir="ref", legacy_weight_dir="w", output_dir="out", **level)


def _converter_method(h2_result, monkeypatch, **level):
    _stop_with_level(monkeypatch, LegacyLDScoreConverter, "_convert", lambda self, **kwargs: kwargs["log_level"])
    LegacyLDScoreConverter().convert(legacy_reference_dir="ref", legacy_weight_dir="w", output_dir="out", **level)


def _build_index(h2_result, monkeypatch, **level):
    _stop_with_level(monkeypatch, gene_ldscore_index, "run_build_gene_ldscore_index_from_args", lambda args: args.log_level)
    config = ldsc.GeneLDScoreIndexBuildConfig(
        baseline_annot_sources=("baseline",), plink_prefix="panel", output_dir="index",
        gene_coordinate_file="genes", genome_build="hg19", snp_identifier="rsid", ld_wind_kb=100,
    )
    ldsc.build_gene_ldscore_index(config, **level)


_CALLS = [_plot, _convert_scale, _convert_ldscores, _converter_method, _build_index]


@pytest.mark.parametrize("call", _CALLS, ids=lambda call: call.__name__.lstrip("_"))
def test_an_omitted_log_level_follows_the_registered_global_config(call, h2_result, monkeypatch, debug_session):
    with pytest.raises(_Requested, match="^DEBUG$"):
        call(h2_result, monkeypatch)


@pytest.mark.parametrize("call", _CALLS, ids=lambda call: call.__name__.lstrip("_"))
def test_an_explicit_log_level_overrides_the_global_config(call, h2_result, monkeypatch, debug_session):
    with pytest.raises(_Requested, match="^ERROR$"):
        call(h2_result, monkeypatch, log_level="ERROR")
