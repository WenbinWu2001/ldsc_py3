"""ldscore rejects R2 matrices built with a positive ``--min-r2``; query-r2-matrix still reads them.

Fixture: the deterministic 32-SNP chromosome-22 1KG-30x subset (hg38; see
``tests/fixtures/minimal_external_resources/README.md``). A positive
``--min-r2`` omits low-R2 pairs, which LD-score sums would count as zero.
"""
from pathlib import Path

import pandas as pd
import pytest

from ldsc import cli

pytest.importorskip("bitarray")
pytest.importorskip("pyarrow")

PLINK = Path(__file__).resolve().parents[1] / "fixtures" / "minimal_external_resources" / "plink" / "hm3_chr22_subset"


@pytest.fixture(scope="module")
def truncated_panel(tmp_path_factory):
    panel = tmp_path_factory.mktemp("panel") / "kb200_min_r2"
    assert cli.run_cli(["build-r2-matrix", "--plink-prefix", str(PLINK), "--output-dir", str(panel),
                        "--snp-identifier", "chr_pos", "--source-genome-build", "hg38",
                        "--ld-wind-kb", "200", "--min-r2", "0.05"]) == 0
    return panel


def test_ldscore_rejects_a_truncated_panel_before_any_chromosome_is_computed(tmp_path, truncated_panel):
    output = tmp_path / "out"
    assert cli.run_cli(["ldscore", "--output-dir", str(output), "--snp-identifier", "chr_pos",
                        "--genome-build", "hg38", "--r2-matrix-dir", str(truncated_panel), "--ld-wind-kb", "200"]) != 0

    issues = pd.read_csv(output / "diagnostics" / "input_issues.tsv", sep="\t")
    assert issues["chrom"].astype(str).tolist() == ["22"]
    assert "omits SNP pairs with unbiased R2 below 0.05" in issues["details"].iloc[0]
    assert "ldsc query-r2-matrix" in issues["details"].iloc[0]
    assert not (output / "diagnostics" / "chromosome-logs").exists()
    assert not (output / "ldscore.baseline.parquet").exists()


def test_query_r2_matrix_still_reads_a_truncated_panel(tmp_path, truncated_panel):
    sidecar = pd.read_csv(truncated_panel / "hg38" / "chr22_meta.tsv.gz", sep="\t", comment="#")
    pairs = tmp_path / "pairs.tsv"
    pd.DataFrame({"CHR_1": ["22"], "POS_1": [sidecar["POS"].iloc[0]],
                  "CHR_2": ["22"], "POS_2": [sidecar["POS"].iloc[1]]}).to_csv(pairs, sep="\t", index=False)
    assert cli.run_cli(["query-r2-matrix", "--r2-matrix-dir", str(truncated_panel), "--genome-build", "hg38",
                        "--snp-identifier", "chr_pos", "--pairs", str(pairs),
                        "--output-dir", str(tmp_path / "query")]) == 0
    assert len(pd.read_csv(tmp_path / "query" / "query_r2_matrix.tsv", sep="\t")) == 1
