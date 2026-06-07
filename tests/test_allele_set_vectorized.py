"""Equivalence tests: vectorized allele-set normalization must match the scalar oracle.

The vectorized ``_allele_set_series_vectorized`` replaces the per-row Python loop in
``allele_set_series`` for performance. Its output must be bit-identical to the
trusted scalar reference (``normalize_allele_set`` / ``_allele_failure_reason``)
for every input, because the resulting keys define allele-aware SNP identity used
to align panels, restriction lists, sumstats, and annotations.
"""

import random

import numpy as np
import pandas as pd

import ldsc._kernel.snp_identity as si


def _scalar_reference(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Compute (values, reasons) row-by-row with the scalar oracle functions."""
    values: list[object] = []
    reasons: list[object] = []
    for a1, a2 in zip(frame["A1"], frame["A2"]):
        reason = si._allele_failure_reason(a1, a2)
        if reason is None:
            values.append(si.normalize_allele_set(a1, a2))
            reasons.append(pd.NA)
        else:
            values.append(pd.NA)
            reasons.append(reason)
    return (
        pd.Series(values, index=frame.index, dtype=object),
        pd.Series(reasons, index=frame.index, dtype=object),
    )


def _assert_matches_oracle(frame: pd.DataFrame) -> None:
    exp_values, exp_reasons = _scalar_reference(frame)
    got_values, got_reasons = si._allele_set_series_vectorized(frame, context="test")

    # Index preserved exactly.
    pd.testing.assert_index_equal(got_values.index, frame.index)
    pd.testing.assert_index_equal(got_reasons.index, frame.index)
    assert got_values.dtype == object
    assert got_reasons.dtype == object

    # Element-by-element equality, treating NA as a distinct sentinel.
    for idx in frame.index:
        ev, er = exp_values.loc[idx], exp_reasons.loc[idx]
        gv, gr = got_values.loc[idx], got_reasons.loc[idx]
        assert pd.isna(ev) == pd.isna(gv), f"value-NA mismatch at {idx}: {ev!r} vs {gv!r}"
        if not pd.isna(ev):
            assert ev == gv, f"value mismatch at {idx}: {ev!r} vs {gv!r}"
        assert pd.isna(er) == pd.isna(gr), f"reason-NA mismatch at {idx}: {er!r} vs {gr!r}"
        if not pd.isna(er):
            assert er == gr, f"reason mismatch at {idx}: {er!r} vs {gr!r}"


def test_handcrafted_matrix_matches_scalar_oracle():
    rows = [
        # valid, both orders and both strands -> canonical "A:C"
        ("A", "C"), ("C", "A"), ("T", "G"), ("G", "T"),
        # valid A/G family
        ("A", "G"), ("G", "A"), ("T", "C"), ("C", "T"),
        # lowercase + surrounding whitespace must clean before canonicalizing
        ("a", "c"), (" A ", "c"), ("g", " t "),
        # missing
        (None, "A"), ("", "A"), (np.nan, "A"), ("A", pd.NA),
        # identical
        ("A", "A"), ("g", "G"),
        # multi-character / invalid base tokens
        ("AC", "T"), ("A", "N"), ("", ""), ("X", "Y"),
        # strand-ambiguous
        ("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"),
        # non-string objects fall through str()/isna handling
        (1, "A"), ("A", 2),
    ]
    frame = pd.DataFrame(rows, columns=["A1", "A2"])
    _assert_matches_oracle(frame)


def test_non_default_index_is_preserved():
    frame = pd.DataFrame(
        {"A1": ["A", "G", "A"], "A2": ["C", "T", "T"]},
        index=[5, 99, 7],
    )
    _assert_matches_oracle(frame)


def test_empty_frame_returns_empty_object_series():
    frame = pd.DataFrame({"A1": pd.Series([], dtype=object), "A2": pd.Series([], dtype=object)})
    values, reasons = si._allele_set_series_vectorized(frame, context="test")
    assert len(values) == 0 and len(reasons) == 0
    assert values.dtype == object and reasons.dtype == object


def test_seeded_fuzz_matches_scalar_oracle():
    rng = random.Random(20260607)
    alphabet = ["A", "C", "G", "T", "a", "c", "g", "t", " A", "C ", "N", "AC", "", "X", "1"]
    specials = [None, np.nan, pd.NA]
    pool = alphabet + specials
    rows = [(rng.choice(pool), rng.choice(pool)) for _ in range(4000)]
    frame = pd.DataFrame(rows, columns=["A1", "A2"])
    _assert_matches_oracle(frame)
