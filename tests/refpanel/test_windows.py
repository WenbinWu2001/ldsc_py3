"""LD-window geometry in ``ldsc._refpanel.windows`` and the reference panel's window preparation.

``window_coordinates`` and ``get_block_lefts`` in every window unit, the
usable-CM and sorted-position checks, and the missing-CM error of
``ldsc._refpanel.adapters._prepare_window``.
"""
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from ldsc.config import LDWindow
from ldsc.errors import LDSCInputError, LDSCInternalError
from ldsc._refpanel import adapters as ref_panel_adapters, windows as ld_windows
from ldsc._refpanel.windows import (
    assert_cm_usable,
    validate_window_positions_sorted,
    window_coordinates,
)


def test_window_coordinates_select_each_window_mode():
    """One window builder serves ldscore and build-r2-matrix for every mode."""
    positions = np.array([100, 250, 900, 2_000])
    cm = np.array([0.0, 0.4, 0.9, 1.6])

    coords, dist = window_coordinates(positions, cm, LDWindow("cm", 1.0))
    np.testing.assert_array_equal(coords, cm)
    assert dist == 1.0
    np.testing.assert_array_equal(ld_windows.get_block_lefts(coords, dist), [0, 0, 0, 2])

    coords, dist = window_coordinates(positions, None, LDWindow("kb", 1.0))
    np.testing.assert_array_equal(coords, positions.astype(float))
    assert dist == 1000.0

    coords, dist = window_coordinates(positions, None, LDWindow("snps", 2))
    np.testing.assert_array_equal(coords, [0.0, 1.0, 2.0, 3.0])
    assert dist == 2.0


def test_window_coordinates_count_panel_rows_when_given():
    """``--ld-wind-snps`` counts R2 matrix sidecar rows, not retained rows."""
    coords, dist = window_coordinates(
        np.array([10, 20, 30]), None, LDWindow("snps", 1), snp_positions=np.array([0, 5, 6]),
    )
    np.testing.assert_array_equal(coords, [0.0, 5.0, 6.0])
    np.testing.assert_array_equal(ld_windows.get_block_lefts(coords, dist), [0, 1, 1])


def test_ldscore_window_rejects_missing_cm():
    metadata = pd.DataFrame({"CHR": ["22"] * 3, "POS": [1, 2, 3], "SNP": ["rs1", "rs2", "rs3"],
                             "CM": [0.0, np.nan, 0.2]})
    config = SimpleNamespace(ld_wind_snps=None, ld_wind_kb=None, ld_wind_cm=1.0, whole_chromosome_ok=True)
    with pytest.raises(LDSCInputError, match="ldscore cannot use `--ld-wind-cm`.*missing CM"):
        ref_panel_adapters._prepare_window(metadata, "22", config, "rsid")


def test_assert_cm_usable_rejects_all_zero():
    with pytest.raises(LDSCInputError, match="all zero|uninformative"):
        assert_cm_usable(pd.Series([0.0, 0.0, 0.0]), chrom="1")


def test_assert_cm_usable_rejects_constant():
    with pytest.raises(LDSCInputError):
        assert_cm_usable(pd.Series([2.5, 2.5, 2.5]), chrom="1")


def test_assert_cm_usable_accepts_two_distinct_values():
    assert_cm_usable(pd.Series([0.0, 0.1, 0.1]), chrom="1")  # no raise


def test_validate_window_positions_sorted_accepts_non_decreasing():
    metadata = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 100, 250]})
    # Non-decreasing (including a tie) must pass without raising.
    validate_window_positions_sorted(metadata, chrom="1")


def test_validate_window_positions_sorted_rejects_unsorted_metadata():
    metadata = pd.DataFrame({"CHR": ["1", "1", "1"], "POS": [100, 300, 250]})
    with pytest.raises(LDSCInternalError, match="non-decreasing order"):
        validate_window_positions_sorted(metadata, chrom="1")


def test_validate_window_positions_sorted_ignores_chromosome_boundary_drop():
    # A position drop across a chromosome boundary is expected, not an error:
    # the window is computed per chromosome.
    metadata = pd.DataFrame({"CHR": ["1", "1", "2"], "POS": [100, 999, 1]})
    validate_window_positions_sorted(metadata, chrom="1")
