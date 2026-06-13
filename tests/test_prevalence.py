import pytest

from ldsc import prevalence
from ldsc.errors import LDSCConfigError, LDSCUsageError


def test_scalar_both_none_is_quantitative():
    assert prevalence.parse_scalar_prevalence(None, None) == (None, None)


def test_scalar_pair():
    assert prevalence.parse_scalar_prevalence(0.5, 0.01) == (0.5, 0.01)


def test_scalar_half_specified_raises():
    with pytest.raises(LDSCUsageError):
        prevalence.parse_scalar_prevalence(0.5, None)


def test_scalar_out_of_range_raises():
    with pytest.raises(LDSCConfigError):
        prevalence.parse_scalar_prevalence(0.5, 1.0)


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
