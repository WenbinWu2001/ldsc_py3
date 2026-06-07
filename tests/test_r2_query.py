"""Tests for the reference-panel R² pair-query API and R²→r converter."""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import ldsc._kernel.ref_panel_builder as kernel_builder
from ldsc._kernel.snp_identity import sidecar_identity_sha256
from ldsc.r2_query import unbiased_r2_to_pearson_r

pytest.importorskip("pyarrow")


def build_test_panel(
    panel_dir: Path,
    *,
    chrom: str = "1",
    snp_identifier: str = "chr_pos_allele_aware",
    genome_build: str = "hg38",
    n_samples: int = 1000,
):
    """Write one binding-valid index panel (chrN_r2.parquet + chrN_meta.tsv.gz).

    Sidecar rows (0-based index → identity):
        0: 1:100 A/G   1: 1:200 C/T   2: 1:300 A/C   3: 1:400 G/T
    Stored pairs (canonical i<j), R2 (adjusted) and SIGN (Pearson r >= 0):
        (0,1) R2=0.64 sign=+   (0,2) R2=0.04 sign=-   (1,2) R2=0.25 sign=+
    Pair (0,3) is intentionally absent (out-of-window placeholder).
    """
    build_dir = panel_dir / genome_build
    build_dir.mkdir(parents=True, exist_ok=True)
    sidecar = pd.DataFrame(
        {
            "CHR": [chrom] * 4,
            "POS": [100, 200, 300, 400],
            "SNP": ["rs1", "rs2", "rs3", "rs4"],
            "A1": ["A", "C", "A", "G"],
            "A2": ["G", "T", "C", "T"],
            "CM": [0.0, 0.1, 0.2, 0.3],
            "MAF": [0.2, 0.3, 0.25, 0.4],
        }
    )
    meta_path = build_dir / f"chr{chrom}_meta.tsv.gz"
    kernel_builder.write_runtime_metadata_sidecar(
        sidecar, meta_path, genome_build=genome_build, snp_identifier=snp_identifier
    )
    identity_hash = sidecar_identity_sha256(sidecar)
    i = np.array([0, 0, 1], dtype=np.int64)
    j = np.array([1, 2, 2], dtype=np.int64)
    r2 = np.array([0.64, 0.04, 0.25], dtype=np.float32)
    sign = np.array([1, -1, 1], dtype=np.int8)
    r2_path = build_dir / f"chr{chrom}_r2.parquet"
    kernel_builder.write_r2_parquet(
        pair_chunks=[(i, j, r2, sign)],
        path=r2_path,
        genome_build=genome_build,
        n_samples=n_samples,
        snp_identifier=snp_identifier,
        n_snps=len(sidecar),
        sidecar_identity_sha256=identity_hash,
        min_r2=0.0,
    )
    return build_dir, meta_path, r2_path


class TestUnbiasedR2ToPearsonR:
    def test_inverts_unbiased_correction_then_roots(self):
        # Build adjusted R2 from a known raw R2, then recover signed r.
        n = 1000
        r2_raw = 0.49
        r2_adj = r2_raw - (1.0 - r2_raw) / (n - 2)
        r = unbiased_r2_to_pearson_r(r2_adj, n, sign=-1)
        assert r == pytest.approx(-0.7, abs=1e-3)

    def test_negative_adjusted_value_still_yields_real_r(self):
        # Near-zero LD: adjusted R2 can be slightly negative; r must be real (not NaN).
        r = unbiased_r2_to_pearson_r(-0.0005, 1000, sign=1)
        assert np.isfinite(r)
        assert r >= 0.0

    def test_sign_none_returns_magnitude(self):
        r = unbiased_r2_to_pearson_r(0.25, 1000, sign=None)
        assert r >= 0.0

    def test_vectorized_arrays(self):
        r2 = np.array([0.25, 0.04], dtype=np.float64)
        sign = np.array([1, -1], dtype=np.int8)
        out = unbiased_r2_to_pearson_r(r2, 1000, sign=sign)
        assert out.shape == (2,)
        assert out[0] > 0 and out[1] < 0


class TestLookupPairsInParquet:
    def _open(self, r2_path):
        import pyarrow.parquet as pq

        return pq.ParquetFile(str(r2_path))

    def test_returns_stored_r2_and_sign(self, tmp_path):
        from ldsc._kernel.r2_query import lookup_pairs_in_parquet

        _, _, r2_path = build_test_panel(tmp_path)
        pf = self._open(r2_path)
        i = np.array([0, 0, 1], dtype=np.int64)
        j = np.array([1, 2, 2], dtype=np.int64)
        r2, sign = lookup_pairs_in_parquet(
            pf, i, j, n_snps=4, r2_scale=32767.0, strategy="stream"
        )
        np.testing.assert_allclose(r2, [0.64, 0.04, 0.25], atol=1.5e-5)
        np.testing.assert_array_equal(sign, [1, -1, 1])

    def test_absent_pair_is_nan_with_zero_sign(self, tmp_path):
        from ldsc._kernel.r2_query import lookup_pairs_in_parquet

        _, _, r2_path = build_test_panel(tmp_path)
        pf = self._open(r2_path)
        i = np.array([0], dtype=np.int64)
        j = np.array([3], dtype=np.int64)  # (0,3) not stored
        r2, sign = lookup_pairs_in_parquet(pf, i, j, n_snps=4, r2_scale=32767.0)
        assert np.isnan(r2[0])
        assert sign[0] == 0

    def test_random_and_stream_strategies_match(self, tmp_path):
        from ldsc._kernel.r2_query import lookup_pairs_in_parquet

        _, _, r2_path = build_test_panel(tmp_path)
        pf = self._open(r2_path)
        i = np.array([0, 0, 1, 0], dtype=np.int64)
        j = np.array([1, 2, 2, 3], dtype=np.int64)
        rg_r2, rg_sign = lookup_pairs_in_parquet(
            pf, i, j, n_snps=4, r2_scale=32767.0, strategy="random"
        )
        st_r2, st_sign = lookup_pairs_in_parquet(
            pf, i, j, n_snps=4, r2_scale=32767.0, strategy="stream"
        )
        np.testing.assert_array_equal(np.nan_to_num(rg_r2, nan=-1), np.nan_to_num(st_r2, nan=-1))
        np.testing.assert_array_equal(rg_sign, st_sign)


class TestR2PanelOpen:
    def test_open_panel_dir_exposes_metadata(self, tmp_path):
        from ldsc.r2_query import R2Panel

        build_test_panel(tmp_path, snp_identifier="chr_pos_allele_aware", n_samples=1234)
        panel = R2Panel.open(tmp_path, genome_build="hg38")
        assert panel.chromosomes == ["1"]
        assert panel.snp_identifier == "chr_pos_allele_aware"
        assert panel.n_samples == 1234

    def test_open_explicit_meta_and_parquet(self, tmp_path):
        from ldsc.r2_query import R2Panel

        _, meta_path, r2_path = build_test_panel(tmp_path)
        panel = R2Panel.open(meta_path=meta_path, parquet_path=r2_path)
        assert panel.chromosomes == ["1"]

    def test_open_requires_exactly_one_input_mode(self, tmp_path):
        from ldsc.errors import LDSCUsageError
        from ldsc.r2_query import R2Panel

        with pytest.raises(LDSCUsageError):
            R2Panel.open()  # neither
        _, meta_path, r2_path = build_test_panel(tmp_path)
        with pytest.raises(LDSCUsageError):
            R2Panel.open(tmp_path, meta_path=meta_path, parquet_path=r2_path)  # both

    def test_binding_mismatch_is_hard_error(self, tmp_path):
        from ldsc.errors import LDSCInputError
        from ldsc.r2_query import R2Panel

        build_dir, meta_path, _ = build_test_panel(tmp_path)
        # Corrupt the sidecar identity by rewriting one allele (breaks the hash).
        with gzip.open(meta_path, "rt") as handle:
            lines = handle.readlines()
        body = [ln for ln in lines if not ln.startswith("#")]
        header = [ln for ln in lines if ln.startswith("#")]
        body[1] = body[1].replace("\tA\t", "\tT\t", 1)  # first data row allele edit
        with gzip.open(meta_path, "wt") as handle:
            handle.writelines(header + body)
        panel = R2Panel.open(tmp_path, genome_build="hg38")
        with pytest.raises(LDSCInputError):
            panel.query_pairs(pd.DataFrame({
                "CHR_1": [1], "POS_1": [100], "A1_1": ["A"], "A2_1": ["G"],
                "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"],
            }))


class TestQueryPairs:
    def _panel(self, tmp_path, mode="chr_pos_allele_aware"):
        from ldsc.r2_query import R2Panel

        build_test_panel(tmp_path, snp_identifier=mode)
        return R2Panel.open(tmp_path, genome_build="hg38")

    def test_stored_diagonal_absent_cross_and_missing(self, tmp_path):
        from ldsc.r2_query import R2Panel

        # Two-chromosome panel so the cross-chromosome row resolves both endpoints.
        build_test_panel(tmp_path, chrom="1", snp_identifier="chr_pos_allele_aware")
        build_test_panel(tmp_path, chrom="2", snp_identifier="chr_pos_allele_aware")
        panel = R2Panel.open(tmp_path, genome_build="hg38")
        pairs = pd.DataFrame(
            {
                # stored (0,1); diagonal (0,0); absent (0,3); not-in-panel; cross-chr
                "CHR_1": [1, 1, 1, 1, 1],
                "POS_1": [100, 100, 100, 100, 100],
                "A1_1": ["A", "A", "A", "A", "A"],
                "A2_1": ["G", "G", "G", "G", "G"],
                "CHR_2": [1, 1, 1, 1, 2],
                "POS_2": [200, 100, 400, 999, 200],
                "A1_2": ["C", "A", "G", "X", "C"],
                "A2_2": ["T", "G", "T", "Y", "T"],
            }
        )
        out = panel.query_pairs(pairs)
        assert out["r2"].iloc[0] == pytest.approx(0.64, abs=1.5e-5)
        assert out["status"].iloc[0] == ""
        assert out["r2"].iloc[1] == 1.0           # diagonal
        assert out["status"].iloc[1] == ""
        assert np.isnan(out["r2"].iloc[2]) and out["status"].iloc[2] == "absent"
        assert np.isnan(out["r2"].iloc[3]) and out["status"].iloc[3] == "not_in_panel"
        assert np.isnan(out["r2"].iloc[4]) and out["status"].iloc[4] == "cross_chromosome"

    def test_sign_harmonized_in_allele_aware_mode(self, tmp_path):
        panel = self._panel(tmp_path)
        # Panel pair (0,1)=SNP1 A/G, SNP2 C/T, stored SIGN=+ (r>=0).
        pairs = pd.DataFrame(
            {
                "CHR_1": [1, 1, 1],
                "POS_1": [100, 100, 100],
                "A1_1": ["A", "G", "G"],   # aligned / swapped / swapped
                "A2_1": ["G", "A", "A"],
                "CHR_2": [1, 1, 1],
                "POS_2": [200, 200, 200],
                "A1_2": ["C", "C", "T"],   # aligned / aligned / swapped
                "A2_2": ["T", "T", "C"],
            }
        )
        out = panel.query_pairs(pairs)
        assert out["sign"].tolist() == [1, -1, 1]  # 0 / 1 / 2 swaps -> +,-,+

    def test_base_mode_ignores_alleles_and_sign_is_na(self, tmp_path):
        panel = self._panel(tmp_path, mode="chr_pos")
        # With alleles supplied vs not supplied: identical r2, sign always NA.
        with_alleles = pd.DataFrame(
            {"CHR_1": [1], "POS_1": [100], "A1_1": ["A"], "A2_1": ["G"],
             "CHR_2": [1], "POS_2": [200], "A1_2": ["C"], "A2_2": ["T"]}
        )
        without = pd.DataFrame({"CHR_1": [1], "POS_1": [100], "CHR_2": [1], "POS_2": [200]})
        out_a = panel.query_pairs(with_alleles)
        out_b = panel.query_pairs(without)
        assert out_a["r2"].iloc[0] == pytest.approx(out_b["r2"].iloc[0], abs=1.5e-5)
        assert pd.isna(out_a["sign"].iloc[0]) and pd.isna(out_b["sign"].iloc[0])
