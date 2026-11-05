from dataclasses import replace

import pytest

from ldsc._regression import prevalence
from ldsc._regression import runner as workflow
from ldsc.config import RegressionConfig
from ldsc.errors import LDSCConfigError, LDSCUsageError
from tests.regression.test_trait_preparation import munged, panel
from tests.worker_fixtures import calculation_body


def test_config_without_prevalence_is_observed_scale():
    config = RegressionConfig()
    assert (config.samp_prev, config.pop_prev) == (None, None)


def test_config_keeps_a_prevalence_pair():
    config = RegressionConfig(samp_prev=0.5, pop_prev=0.01)
    assert (config.samp_prev, config.pop_prev) == (0.5, 0.01)


def test_config_reads_nan_as_a_quantitative_trait():
    config = RegressionConfig(samp_prev=float("nan"), pop_prev=float("nan"))
    assert (config.samp_prev, config.pop_prev) == (None, None)


@pytest.mark.parametrize("samp_prev, pop_prev", [(0.5, None), (None, 0.01), (0.5, float("nan"))])
def test_config_rejects_half_a_prevalence_pair(samp_prev, pop_prev):
    # Before, the config kept the half pair and every fit silently stayed on the observed scale.
    with pytest.raises(LDSCConfigError, match="exactly one of sample/population prevalence"):
        RegressionConfig(samp_prev=samp_prev, pop_prev=pop_prev)


@pytest.mark.parametrize("samp_prev, pop_prev", [(0.5, 1.0), (0.0, 0.01), (1.5, 0.01)])
def test_config_rejects_a_prevalence_outside_the_open_unit_interval(samp_prev, pop_prev):
    with pytest.raises(LDSCConfigError, match=r"outside the open interval \(0, 1\)"):
        RegressionConfig(samp_prev=samp_prev, pop_prev=pop_prev)


def test_config_rejects_a_prevalence_that_is_not_a_number():
    with pytest.raises(LDSCConfigError, match="not a number"):
        RegressionConfig(samp_prev="0.5", pop_prev=0.01)


@pytest.mark.parametrize("prevalences, error", [
    ([(0.5, None), (None, None)], LDSCConfigError),
    ([(0.5, 1.2), (None, None)], LDSCConfigError),
    ([(0.5, 0.01)], LDSCUsageError),
])
def test_rg_pairs_reject_bad_prevalences_before_any_fit(monkeypatch, prevalences, error):
    def no_fit(*args, **kwargs):
        raise AssertionError("a fit ran before the prevalences were checked")

    monkeypatch.setattr(workflow, "_fit_h2_dataset", no_fit)
    tables = [munged("A", "chr_pos_allele_aware", 1), munged("B", "chr_pos_allele_aware", 2)]
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))
    # In-process, so the replaced fit is the one the run would call.
    with calculation_body(), pytest.raises(error):
        runner.estimate_rg_pairs(tables, panel(), prevalences=prevalences)


def test_positional_list_aligns_and_marks_qt():
    out = prevalence.parse_positional_prevalences("0.5,0.5,nan", "0.01,0.02,nan", n_traits=3)
    assert out == [(0.5, 0.01), (0.5, 0.02), (None, None)]


def test_positional_empty_token_is_quantitative():
    out = prevalence.parse_positional_prevalences("0.5,,0.3", "0.01,,0.02", n_traits=3)
    assert out == [(0.5, 0.01), (None, None), (0.3, 0.02)]


def test_positional_length_mismatch_raises():
    with pytest.raises(LDSCUsageError):
        prevalence.parse_positional_prevalences("0.5,0.5", "0.01,0.01", n_traits=3)


def test_positional_half_specified_trait_raises():
    with pytest.raises((LDSCUsageError, LDSCConfigError)):
        prevalence.parse_positional_prevalences("0.5,nan", "0.01,0.01", n_traits=2)


def test_manifest_superset_lookup_by_name(tmp_path):
    manifest = tmp_path / "prev.tsv"
    manifest.write_text(
        "# standing repository\n"
        "trait_name\tsamp_prev\tpop_prev\n"
        "scz2022\t0.43\t0.01\n"
        "scz2014\t0.40\t0.01\n"
        "height\tnan\tnan\n",
        encoding="utf-8",
    )
    out = prevalence.parse_manifest_prevalences(
        str(manifest), trait_names=["scz2022", "height"], trait_paths=["a", "b"]
    )
    assert out == [(0.43, 0.01), (None, None)]


def test_manifest_whitespace_delimiter(tmp_path):
    manifest = tmp_path / "prev.txt"
    manifest.write_text("trait_name samp_prev pop_prev\nscz2022   0.43   0.01\n", encoding="utf-8")
    out = prevalence.parse_manifest_prevalences(
        str(manifest), trait_names=["scz2022"], trait_paths=["a"]
    )
    assert out == [(0.43, 0.01)]


def test_manifest_missing_trait_raises(tmp_path):
    manifest = tmp_path / "prev.tsv"
    manifest.write_text("trait_name\tsamp_prev\tpop_prev\nscz2022\t0.43\t0.01\n", encoding="utf-8")
    with pytest.raises(LDSCUsageError):
        prevalence.parse_manifest_prevalences(
            str(manifest), trait_names=["scz2022", "bip"], trait_paths=["a", "b"]
        )


def test_resolve_rg_rejects_both_schemes():
    with pytest.raises(LDSCUsageError):
        prevalence.resolve_rg_prevalences(
            samp_prev="0.5,0.5",
            pop_prev="0.01,0.01",
            manifest_path="x.tsv",
            trait_names=["a", "b"],
            trait_paths=["a", "b"],
        )


def test_resolve_rg_none_returns_none():
    assert (
        prevalence.resolve_rg_prevalences(
            samp_prev=None, pop_prev=None, manifest_path=None,
            trait_names=["a", "b"], trait_paths=["a", "b"],
        )
        is None
    )


def test_resolve_rg_manifest_duplicate_name_raises(tmp_path):
    manifest = tmp_path / "prev.tsv"
    manifest.write_text("trait_name\tsamp_prev\tpop_prev\nscz\t0.43\t0.01\n", encoding="utf-8")
    with pytest.raises(LDSCUsageError):
        prevalence.resolve_rg_prevalences(
            samp_prev=None, pop_prev=None, manifest_path=str(manifest),
            trait_names=["scz", "scz"], trait_paths=["/data/a.sumstats", "/data/b.sumstats"],
        )


def test_a_failed_rg_pair_keeps_the_prevalences_of_its_traits():
    # The two halves of one trait share no SNP, so their pair fails; the pair
    # used to show NaN prevalences and the observed scale although the run
    # applied liability-scale prevalences to both traits.
    whole = munged("A", "chr_pos_allele_aware", 1)
    halves = [replace(whole, data=rows.reset_index(drop=True), trait_name=name, source_path=f"{name}.parquet")
              for name, rows in (("first", whole.data.iloc[:80]), ("second", whole.data.iloc[80:]))]
    runner = workflow.RegressionRunner(regression_config=RegressionConfig(num_blocks=10))

    result = runner.estimate_rg_pairs(halves, panel(), prevalences=[(0.5, 0.01), (0.3, 0.02)])

    row = result.rg_full.iloc[0]
    assert row["status"] == "failed"
    assert (row["samp_prev_1"], row["pop_prev_1"], row["samp_prev_2"], row["pop_prev_2"]) == (0.5, 0.01, 0.3, 0.02)
    metadata = result.per_pair_metadata[0]
    assert (metadata["samp_prev_1"], metadata["pop_prev_1"]) == (0.5, 0.01)
    assert (metadata["samp_prev_2"], metadata["pop_prev_2"]) == (0.3, 0.02)
    assert metadata["scale"] == "liability"
