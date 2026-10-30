"""Row order and alignment of direct LD-score query batches against the batch-1 baseline.

Only the first execution batch projects and returns the baseline. Later batches
return query scores alone, so their rows must follow the shared baseline's rows
exactly; these tests pin that contract on the PLINK and R2 matrix paths.
"""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

from ldsc import AnnotationBundle, GlobalConfig, LDScoreCalculator, LDScoreConfig, PlinkRefPanel, RefPanelConfig
from ldsc._core.row_alignment import concat_chromosome_tables, sort_frame_by_genomic_position
from ldsc._ldscore import calculator as ldscore_calculator
from ldsc._ldscore import fragments
from ldsc._ldscore.outputs import LDScoreOutputConfig
from ldsc._refpanel.adapters import R2MatrixRefPanel
from ldsc.errors import LDSCInputError, LDSCInternalError
from tests.ldscore.test_ldscore_parallelism import _write_index_panel

PLINK_FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/minimal_external_resources/plink/hm3_chr22_subset"


def _r2_case(tmp_path):
    """Chromosomes 10, 2 and 1 in that input order; their labels sort differently as text."""
    def chromosome(label, positions):
        frame = pd.DataFrame({"CHR": label, "POS": positions, "SNP": [f"rs{label}_{p}" for p in positions],
                              "A1": "A", "A2": "C", "MAF": np.linspace(0.1, 0.45, len(positions))})
        pairs = sorted([(i, i + 1, 0.2 + 0.1 * (i % 3)) for i in range(len(positions) - 1)] + [(0, 2, -0.05)])
        return frame, pairs
    chroms = {"10": chromosome("10", [150, 250, 350, 450]), "2": chromosome("2", [100, 200, 300]),
              "1": chromosome("1", [120, 220, 320, 420, 520])}
    _write_index_panel(tmp_path / "panel", chroms, snp_identifier="chr_pos")
    metadata = pd.concat([frame[["CHR", "POS", "SNP"]] for frame, _ in chroms.values()], ignore_index=True)
    rng = np.random.default_rng(3)
    baseline = pd.DataFrame({"base": np.ones(len(metadata)), "half": (np.arange(len(metadata)) % 2).astype(float)})
    queries = pd.DataFrame({**{f"q{i}": rng.choice([0., 1.], len(metadata)) for i in range(4)},
                            "q4": rng.gamma(2., 1., len(metadata))})
    config = GlobalConfig(snp_identifier="chr_pos", genome_build="hg19")
    panel = R2MatrixRefPanel(config, RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(tmp_path / "panel")))
    return metadata, baseline, queries, config, panel, LDScoreConfig(ld_wind_kb=1, whole_chromosome_ok=True)


def _plink_case(tmp_path):
    bim = pd.read_csv(PLINK_FIXTURE.with_suffix(".bim"), sep=r"\s+", header=None,
                      names=["CHR", "SNP", "CM", "POS", "A1", "A2"])
    metadata = bim.loc[:, ["CHR", "SNP", "CM", "POS"]].assign(CHR=bim.CHR.astype(str))
    rng = np.random.default_rng(5)
    baseline = pd.DataFrame({"base": np.ones(len(metadata)), "third": (np.arange(len(metadata)) % 3 == 0) * 1.})
    queries = pd.DataFrame({**{f"q{i}": rng.choice([0., 1.], len(metadata)) for i in range(4)},
                            "q4": rng.gamma(2., 1., len(metadata))})
    config = GlobalConfig(snp_identifier="chr_pos", genome_build="hg38")
    panel = PlinkRefPanel(config, RefPanelConfig(backend="plink", plink_prefix=str(PLINK_FIXTURE)))
    return metadata, baseline, queries, config, panel, LDScoreConfig(ld_wind_snps=10)


CASES = {"r2_matrix": _r2_case, "plink": _plink_case}


def _run(case, output, width, workers=1):
    metadata, baseline, queries, config, panel, numerical = case
    with AnnotationBundle.from_frames(metadata, baseline, queries, config_snapshot=config) as bundle:
        return LDScoreCalculator().run(bundle, panel, replace(numerical, query_batch_size=width, workers=workers),
                                       config, output_config=LDScoreOutputConfig(output))


def _identity(path):
    columns = [name for name in ("CHR", "SNP", "POS", "A1", "A2") if name in pq.read_schema(path).names]
    return pd.read_parquet(path, columns=columns)


def _row_groups(path):
    footer = pq.ParquetFile(path).metadata
    return [footer.row_group(index).num_rows for index in range(footer.num_row_groups)]


def test_later_batches_carry_no_baseline_table(tmp_path, monkeypatch):
    seen = []
    original = ldscore_calculator._DirectProducer.reduce

    def reduce(self, outcomes, batch):
        tables = [outcome.result.baseline_table is None for outcome in outcomes.values() if not outcome.skipped]
        staged, variation = original(self, outcomes, batch)
        seen.append((batch.ordinal, tables, staged.baseline_table is None))
        return staged, variation

    monkeypatch.setattr(ldscore_calculator._DirectProducer, "reduce", reduce)
    _run(_r2_case(tmp_path), tmp_path / "out", width=2)
    assert [ordinal for ordinal, _, _ in seen] == [0, 1, 2]
    assert seen[0][1] == [False, False, False] and not seen[0][2]
    assert all(tables == [True, True, True] and staged for _, tables, staged in seen[1:])


@pytest.mark.parametrize("backend", sorted(CASES))
@pytest.mark.parametrize("width,workers", [(1, 1), (2, 2)])
def test_multi_batch_outputs_match_one_batch(tmp_path, backend, width, workers):
    case = CASES[backend](tmp_path)
    single = _run(case, tmp_path / "single", width=100)
    multi = _run(case, tmp_path / "multi", width=width, workers=workers)
    queries = list(case[2].columns)
    assert len(multi.query_batches) == -(-len(queries) // width)
    for name in ("ldscore.baseline.parquet", "ldscore.overlap.parquet"):
        assert (tmp_path / "multi" / name).read_bytes() == (tmp_path / "single" / name).read_bytes()
    assert multi.count_records == single.count_records
    pd.testing.assert_frame_equal(multi.read_queries(queries), single.read_queries(queries), check_exact=True)
    metadata = json.loads((tmp_path / "multi" / "metadata.json").read_text())
    baseline_path = tmp_path / "multi" / "ldscore.baseline.parquet"
    baseline_identity, baseline_groups = _identity(baseline_path), _row_groups(baseline_path)
    for batch in metadata["query_batches"]:
        path = tmp_path / "multi" / batch["file"]
        pd.testing.assert_frame_equal(_identity(path), baseline_identity, check_exact=True)
        assert _row_groups(path) == baseline_groups
        assert [group["n_rows"] for group in batch["row_groups"]] == baseline_groups


def test_published_rows_follow_genomic_chromosome_order(tmp_path):
    _run(_r2_case(tmp_path), tmp_path / "out", width=2)
    identity = _identity(tmp_path / "out" / "ldscore.baseline.parquet")
    assert identity.CHR.drop_duplicates().tolist() == ["1", "2", "10"]
    assert identity.groupby("CHR", sort=False).POS.apply(lambda pos: pos.is_monotonic_increasing).all()


def _corrupt_second_batch(monkeypatch, change):
    original = fragments.write_query_fragments
    calls = []

    def write(result, *args, **kwargs):
        calls.append(result)
        if len(calls) == 2:
            path = next(path for path in result.query_fragments if len(pd.read_parquet(path)) > 1)
            frame = pd.read_parquet(path)
            change(frame).to_parquet(path, index=False)
        return original(result, *args, **kwargs)

    monkeypatch.setattr(fragments, "write_query_fragments", write)


@pytest.mark.parametrize("change", [
    pytest.param(lambda frame: frame.iloc[[1, 0, *range(2, len(frame))]].reset_index(drop=True), id="swapped-rows"),
    pytest.param(lambda frame: frame.assign(SNP=["renamed", *frame.SNP.iloc[1:]]), id="other-snp-label"),
])
def test_misaligned_later_batch_fails_without_publishing(tmp_path, monkeypatch, change):
    _corrupt_second_batch(monkeypatch, change)
    output = tmp_path / "out"
    with pytest.raises(LDSCInputError, match="SNP row mismatch"):
        _run(_r2_case(tmp_path), output, width=2)
    assert not (output / "metadata.json").exists()
    assert not list(output.glob("ldscore.*.parquet"))
    assert not list(output.glob(".ldsc-*"))


def _piece(chrom, positions, snps):
    return pd.DataFrame({"CHR": chrom, "SNP": snps, "POS": positions, "value": np.arange(len(positions), dtype=float)})


def test_chromosome_tables_concatenate_in_genomic_order():
    pieces = [("10", _piece("10", [5, 7], ["a", "b"])), ("2", _piece("2", [3, 3, 9], ["x", "y", "a"])),
              ("1", _piece("1", [1, 4], ["c", "d"]))]
    expected = sort_frame_by_genomic_position(pd.concat([table for _, table in pieces], ignore_index=True))
    pd.testing.assert_frame_equal(concat_chromosome_tables(pieces), expected, check_exact=True)


@pytest.mark.parametrize("pieces", [
    pytest.param([("1", _piece("1", [4, 1], ["d", "c"]))], id="unsorted-piece"),
    pytest.param([("1", _piece("1", [3, 3], ["y", "x"]))], id="unsorted-tie"),
    pytest.param([("1", _piece("2", [1, 4], ["c", "d"]))], id="other-chromosome"),
])
def test_chromosome_tables_reject_pieces_out_of_order(pieces):
    with pytest.raises(LDSCInternalError):
        concat_chromosome_tables(pieces)
