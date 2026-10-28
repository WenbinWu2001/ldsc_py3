"""PLINK LD scores use the exact requested window, whatever ``--snp-batch-size``.

Fixture: the deterministic 32-SNP chromosome-22 1KG-30x subset (3202 samples,
hg38 positions; see ``tests/fixtures/minimal_external_resources/README.md``).
Its SNPs are 30-900 kb apart, so ``--ld-wind-snps 3`` and ``--ld-wind-kb 200``
cover only a few neighbours, while legacy LDSC chunk rounding at a batch of 128
would have widened every window to the whole chromosome.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc import cli

pytest.importorskip("bitarray")
pytest.importorskip("pyarrow")

PLINK = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"
IDENTITY = ["--snp-identifier", "chr_pos"]
WINDOWS = [["--ld-wind-snps", "3"], ["--ld-wind-kb", "200"]]
SCORES = ["base", "regression_ld_scores"]


def _ldscore(output_dir, *args):
    status = cli.run_cli(["ldscore", "--output-dir", str(output_dir), *IDENTITY, "--genome-build", "hg38", *args])
    assert status == 0
    return pd.read_parquet(output_dir / "ldscore.baseline.parquet").set_index(["CHR", "POS"])


@pytest.mark.parametrize("window", WINDOWS, ids=lambda w: "".join(w))
def test_plink_ld_scores_do_not_depend_on_snp_batch_size(tmp_path, window):
    scores = {
        batch: _ldscore(tmp_path / f"batch{batch}", "--plink-prefix", str(PLINK), *window,
                        "--snp-batch-size", str(batch))
        for batch in (1, 3, 128)
    }
    for batch in (3, 128):
        pd.testing.assert_index_equal(scores[batch].index, scores[1].index)
        # Values are float64 sums written as float32; only summation order differs.
        np.testing.assert_allclose(scores[batch][SCORES], scores[1][SCORES], rtol=1e-6, atol=1e-6)


@pytest.mark.parametrize("window", WINDOWS, ids=lambda w: "".join(w))
def test_plink_and_parquet_backends_agree_on_the_same_window(tmp_path, window):
    panel = tmp_path / "panel"
    assert cli.run_cli(["build-r2-matrix", "--plink-prefix", str(PLINK), "--output-dir", str(panel), *IDENTITY,
                        "--source-genome-build", "hg38", *window]) == 0
    from_plink = _ldscore(tmp_path / "plink", "--plink-prefix", str(PLINK), *window)
    from_parquet = _ldscore(tmp_path / "parquet", "--r2-matrix-dir", str(panel), *window)

    pd.testing.assert_index_equal(from_parquet.index, from_plink.index)
    # int16 R2 storage rounds each pair by at most 0.5/32767 (1.5e-5); a SNP
    # here has at most 31 partners, so scores agree within 31 * 1.5e-5.
    np.testing.assert_allclose(from_parquet[SCORES], from_plink[SCORES], rtol=0, atol=5e-4)


@pytest.mark.parametrize("window", WINDOWS, ids=lambda w: "".join(w))
def test_build_r2_matrix_pairs_do_not_depend_on_snp_batch_size(tmp_path, window):
    """The builder walks the same chunks but emits only in-window pairs.

    Pairs, their order, and signs are identical. The builder correlates float32
    genotypes, so the batch shape changes summation order inside the matrix
    product; a value at an int16 rounding boundary may then move by one stored
    unit (1/32767 in R2).
    """
    import pyarrow.parquet as pq

    tables = {}
    for batch in (2, 3, 512):
        out = tmp_path / f"batch{batch}"
        assert cli.run_cli(["build-r2-matrix", "--plink-prefix", str(PLINK), "--output-dir", str(out), *IDENTITY,
                            "--source-genome-build", "hg38", *window, "--snp-batch-size", str(batch)]) == 0
        tables[batch] = pq.read_table(out / "hg38" / "chr22_r2.parquet").to_pandas()
    for batch in (3, 512):
        pd.testing.assert_frame_equal(tables[batch][["IDX_1", "IDX_2", "SIGN_R"]], tables[2][["IDX_1", "IDX_2", "SIGN_R"]])
        stored_units = (tables[batch]["R2"].astype(int) - tables[2]["R2"].astype(int)).abs()
        assert stored_units.max() <= 1
