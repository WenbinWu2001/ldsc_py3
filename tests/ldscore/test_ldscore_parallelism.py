"""Tests for cross-chromosome LD-score parallelism (`workers`).

The fixtures here build a minimal but faithful two-chromosome canonical
index-format reference panel (the same `IDX_1/IDX_2/R2/SIGN_R` parquet + sidecar
that `build-r2-matrix` emits, written via the real builder helpers) so the
parallel path is exercised end-to-end through `SortedR2BlockReader`, not a mock.
"""

from __future__ import annotations
from tests.annotation_fixtures import make_annotation_bundle

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from ldsc.errors import EmptyReferenceSNPs, LDSCInputError

_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None
pytestmark = pytest.mark.skipif(not _HAS_PYARROW, reason="pyarrow is required for R2 matrices")


# --- Fixture builders -------------------------------------------------------

# Per-chromosome: (sidecar panel rows, in-window pairs as (i, j, unbiased_r2)).
_PANEL_CHROMS: dict[str, tuple[pd.DataFrame, list[tuple[int, int, float]]]] = {
    "1": (
        pd.DataFrame(
            {
                "CHR": ["1", "1", "1", "1"],
                "POS": [100, 200, 300, 400],
                "SNP": ["rs1", "rs2", "rs3", "rs4"],
                "A1": ["A", "A", "A", "A"],
                "A2": ["C", "C", "C", "C"],
                "MAF": [0.3, 0.25, 0.4, 0.35],
            }
        ),
        [(0, 1, 0.5), (0, 2, 0.1), (1, 2, 0.4), (2, 3, 0.6)],
    ),
    "2": (
        pd.DataFrame(
            {
                "CHR": ["2", "2", "2"],
                "POS": [100, 200, 300],
                "SNP": ["rs5", "rs6", "rs7"],
                "A1": ["A", "A", "A"],
                "A2": ["C", "C", "C"],
                "MAF": [0.3, 0.45, 0.2],
            }
        ),
        [(0, 1, 0.5), (1, 2, 0.3)],
    ),
}


def _write_index_panel(
    build_dir: Path,
    chroms: dict[str, tuple[pd.DataFrame, list[tuple[int, int, float]]]],
    *,
    snp_identifier: str = "rsid",
    genome_build: str = "hg19",
    n_samples: int = 200,
) -> None:
    """Write canonical index R2 parquet + sidecar for each chromosome."""
    from ldsc._refpanel import build_kernel as kb
    from ldsc._core.snp_identity import sidecar_identity_sha256
    from ldsc._refpanel.r2_format import (
        build_runtime_metadata_table,
        write_r2_parquet,
        write_runtime_metadata_sidecar,
    )

    build_dir.mkdir(parents=True, exist_ok=True)
    for chrom, (panel, pairs) in chroms.items():
        positions = panel["POS"].to_numpy(dtype=np.int64)
        runtime_meta = build_runtime_metadata_table(metadata=panel, positions=positions, cm_values=None)
        identity_hash = sidecar_identity_sha256(runtime_meta)
        i = np.array([p[0] for p in pairs], dtype=np.int64)
        j = np.array([p[1] for p in pairs], dtype=np.int64)
        r2 = np.array([p[2] for p in pairs], dtype=np.float32)
        sign = np.ones(len(i), dtype=np.int8)
        write_r2_parquet(
            pair_chunks=[(i, j, r2, sign)],
            path=build_dir / f"chr{chrom}_r2.parquet",
            genome_build=genome_build,
            n_samples=n_samples,
            snp_identifier=snp_identifier,
            n_snps=len(panel),
            sidecar_identity_sha256=identity_hash,
            ld_window_mode="kb",
            ld_window_value=1000.0,
        )
        write_runtime_metadata_sidecar(
            runtime_meta,
            build_dir / f"chr{chrom}_meta.tsv.gz",
            genome_build=genome_build,
            snp_identifier=snp_identifier,
        )


def _run_minimal_ldscore(build_dir: Path, out_dir: Path, workers: int) -> Path:
    """Run a synthetic-base LD-score job over the index panel; return out_dir."""
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc._ldscore.calculator import run_ldscore

    regression_snps = build_dir / "regression_snps.tsv"
    pd.concat(
        [
            pd.read_csv(path, sep="\t", comment="#", usecols=["SNP"])
            for path in sorted(build_dir.glob("chr*_meta.tsv.gz"))
        ],
        ignore_index=True,
    ).to_csv(regression_snps, sep="\t", index=False)
    set_global_config(GlobalConfig(snp_identifier="rsid"))
    run_ldscore(
        r2_matrix_dir=str(build_dir),
        output_dir=str(out_dir),
        regr_snps_file=str(regression_snps),
        ld_wind_kb=1.0,
        snp_batch_size=2,
        yes_really=True,
        workers=workers,
        regr_snps_exclude_regions="none",  # rsid smoke test; region exclusion is out of scope
    )
    return out_dir


@pytest.fixture
def two_chrom_panel(tmp_path: Path) -> Path:
    build_dir = tmp_path / "panel" / "hg19"
    _write_index_panel(build_dir, _PANEL_CHROMS)
    return build_dir


@pytest.fixture
def run_minimal_ldscore(two_chrom_panel: Path, tmp_path: Path):
    def _runner(workers: int, out_subdir: str) -> Path:
        return _run_minimal_ldscore(two_chrom_panel, tmp_path / out_subdir, workers)

    return _runner


# --- Fixture smoke test -----------------------------------------------------


def test_fixture_builds_and_runs_sequential(run_minimal_ldscore):
    out = run_minimal_ldscore(workers=1, out_subdir="seq_smoke")
    baseline = pd.read_parquet(out / "ldscore.baseline.parquet")
    assert baseline["CHR"].astype(str).tolist() == ["1", "1", "1", "1", "2", "2", "2"]
    assert baseline["SNP"].tolist() == ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7"]
    assert "base" in baseline.columns


# --- LDScoreConfig.workers --------------------------------------------------


def _cfg(**kw):
    from ldsc.config import LDScoreConfig

    base = dict(ld_wind_cm=1.0)
    base.update(kw)
    return LDScoreConfig(**base)


def test_threads_defaults_to_one_sequential():
    assert _cfg().workers == 1


def test_threads_negative_and_positive_allowed(monkeypatch):
    # Positive requests are checked against the detected CPU count.
    monkeypatch.setattr("joblib.cpu_count", lambda: 8)
    assert _cfg(workers=-1).workers == -1
    assert _cfg(workers=-2).workers == -2
    assert _cfg(workers=4).workers == 4


# --- Task 3: module-level chromosome worker ---------------------------------


def _chrom1_inputs(build_dir: Path):
    """Build (chrom, bundle, spec, ldscore_config, global_config, regression) for chr1."""
    from ldsc.config import GlobalConfig, LDScoreConfig, RefPanelConfig

    panel, _pairs = _PANEL_CHROMS["1"]
    metadata = pd.DataFrame(
        {
            "CHR": panel["CHR"].astype(str),
            "SNP": panel["SNP"].astype(str),
            "POS": panel["POS"].astype(int),
            "CM": [np.nan] * len(panel),
        }
    )
    bundle = make_annotation_bundle(
        metadata=metadata,
        baseline_annotations=pd.DataFrame({"base": np.ones(len(panel), dtype=np.float32)}),
        query_annotations=pd.DataFrame(index=metadata.index),
        baseline_columns=["base"],
        query_columns=[],
        chromosomes=["1"],
        source_summary={},
        config_snapshot=GlobalConfig(snp_identifier="rsid"),
    )
    spec = RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(build_dir))
    ldscore_config = LDScoreConfig(ld_wind_kb=1.0, snp_batch_size=2, whole_chromosome_ok=True)
    return "1", bundle, spec, ldscore_config, GlobalConfig(snp_identifier="rsid"), None


def test_compute_one_chromosome_returns_success_outcome(two_chrom_panel):
    from ldsc._ldscore.calculator import _compute_one_chromosome

    chrom, bundle, spec, ldcfg, gcfg, reg = _chrom1_inputs(two_chrom_panel)
    outcome = _compute_one_chromosome(chrom, bundle, spec, ldcfg, gcfg, regression_snps=reg)
    assert outcome.chrom == "1"
    assert outcome.skipped is False
    assert outcome.result is not None
    assert outcome.result.baseline_table["SNP"].tolist() == ["rs1", "rs2", "rs3", "rs4"]


@pytest.mark.parametrize("error_type", [ValueError, LDSCInputError])
@pytest.mark.parametrize("message", [
    "No retained annotation SNPs remain on chromosome 1 after parquet intersection.",
    "ldscore retained no annotation SNPs on chromosome 1 after intersecting with the PLINK panel.",
])
def test_worker_propagates_untyped_errors_even_with_empty_intersection_wording(
    two_chrom_panel, monkeypatch, error_type, message,
):
    from ldsc._ldscore.calculator import LDScoreCalculator, _compute_one_chromosome

    error = error_type(message)

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(LDScoreCalculator, "compute_chromosome", fail)
    chrom, bundle, spec, ldcfg, gcfg, reg = _chrom1_inputs(two_chrom_panel)
    with pytest.raises(error_type) as caught:
        _compute_one_chromosome(chrom, bundle, spec, ldcfg, gcfg, regression_snps=reg)
    assert caught.value is error


def test_worker_skips_typed_empty_reference_independent_of_message(two_chrom_panel, monkeypatch):
    from ldsc._ldscore.calculator import LDScoreCalculator, _compute_one_chromosome

    def empty(*args, **kwargs):
        raise EmptyReferenceSNPs("The requested SNP filters retained no computational rows.")

    monkeypatch.setattr(LDScoreCalculator, "compute_chromosome", empty)
    chrom, bundle, spec, ldcfg, gcfg, reg = _chrom1_inputs(two_chrom_panel)
    outcome = _compute_one_chromosome(chrom, bundle, spec, ldcfg, gcfg, regression_snps=reg)
    assert outcome.chrom == "1"
    assert outcome.skipped is True
    assert outcome.result is None
    assert outcome.skip_message == "The requested SNP filters retained no computational rows."


# --- Task 4: pool initializer -----------------------------------------------


def _inspect_ld_worker():
    from ldsc._ldscore.calculator import _WORKER_STATE
    from ldsc._runtime.worker_runtime import native_thread_state
    return _WORKER_STATE, native_thread_state(), dict((key, os.environ[key]) for key in
                                                   ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"))


def test_init_worker_overrides_inherited_limits_in_child(monkeypatch):
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    from ldsc._ldscore.calculator import _init_worker
    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    with ProcessPoolExecutor(max_workers=1, mp_context=mp.get_context("spawn"), initializer=_init_worker,
                             initargs=({"rs1", "rs2"}, None, "WARNING")) as pool:
        state, native, environment = pool.submit(_inspect_ld_worker).result()
    assert state["regression_snps"] == {"rs1", "rs2"}
    assert native["blas_target"] == 1
    assert set(environment.values()) == {"1"}
    assert os.environ["OMP_NUM_THREADS"] == "4"


# --- Task 5: sequential-vs-pool equivalence ---------------------------------


def test_parallel_output_matches_sequential(run_minimal_ldscore):
    seq = run_minimal_ldscore(workers=1, out_subdir="seq")
    par = run_minimal_ldscore(workers=3, out_subdir="par")
    seq_base = pd.read_parquet(seq / "ldscore.baseline.parquet")
    par_base = pd.read_parquet(par / "ldscore.baseline.parquet")
    pd.testing.assert_frame_equal(seq_base, par_base)


def test_partitioned_parallel_preserves_scores_counts_overlap_and_classification(two_chrom_panel, tmp_path):
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc._ldscore.calculator import run_ldscore
    from ldsc import load_ldscore_from_dir

    metadata = pd.concat([panel for panel, _ in _PANEL_CHROMS.values()], ignore_index=True)
    baseline = metadata.loc[:, ["CHR", "POS", "SNP"]].assign(
        base=1, cont=[0, 0.5, 1, 1.5, 0, 2, 3]
    )
    query = metadata.loc[:, ["CHR", "POS", "SNP"]].assign(query=[1, 0, 1, 0, 1, 0, 1])
    baseline_path = tmp_path / "baseline.annot"
    query_path = tmp_path / "query.annot"
    baseline.to_csv(baseline_path, sep="\t", index=False)
    query.to_csv(query_path, sep="\t", index=False)
    regression_path = tmp_path / "regression.tsv"
    regression_path.write_text("SNP\nrs6\nrs1\nrs5\nrs3\n", encoding="utf-8")

    outputs = []
    for workers in (1, 2):
        set_global_config(GlobalConfig(snp_identifier="rsid"))
        output = tmp_path / f"threads_{workers}"
        result = run_ldscore(
            r2_matrix_dir=str(two_chrom_panel),
            baseline_annot_sources=[str(baseline_path)],
            query_annot_sources=[str(query_path)],
            output_dir=str(output),
            regr_snps_file=str(regression_path),
            ld_wind_kb=1.0,
            common_maf_min=0.3,
            snp_batch_size=2,
            yes_really=True,
            workers=workers,
            regr_snps_exclude_regions="none",
        )
        assert not hasattr(result, "annotation_fingerprints")
        assert all(not hasattr(chrom, "common_annotation_values") for chrom in result.chromosome_results)
        persisted = json.loads((output / "metadata.json").read_text(encoding="utf-8"))
        assert "annotation_fingerprints" not in persisted
        loaded = load_ldscore_from_dir(output)
        assert loaded.annotation_types == {"base": "binary", "cont": "quantitative", "query": "binary"}
        assert loaded.baseline_table["SNP"].tolist() == ["rs1", "rs3", "rs5", "rs6"]
        np.testing.assert_allclose(
            loaded.baseline_table[["base", "cont", "regression_ld_scores"]],
            [[1.6, 0.35, 1.1], [2.1, 2.1, 1.1], [1.5, 1.0, 1.5], [1.8, 2.9, 1.5]],
            rtol=0, atol=1e-4,  # The fixture's R2 values use int16 parquet quantization.
        )
        np.testing.assert_allclose(loaded.read_queries(["query"])["query"], [1.1, 1.1, 1.0, 0.8], rtol=0, atol=1e-4)
        assert [record["all_reference_snp_count"] for record in loaded.count_records] == [7, 8, 4]
        assert [record["common_reference_snp_count"] for record in loaded.count_records] == [5, 4.5, 3]
        np.testing.assert_array_equal(loaded.overlap.baseline_block_all, [[7, 8, 4], [8, 16.5, 4]])
        np.testing.assert_array_equal(loaded.overlap.baseline_block_common, [[5, 4.5, 3], [4.5, 7.25, 1]])
        np.testing.assert_array_equal(loaded.overlap.query_diagonal_all, [4])
        np.testing.assert_array_equal(loaded.overlap.query_diagonal_common, [3])
        assert loaded.overlap.total_all_reference_snps == 7
        assert loaded.overlap.total_common_reference_snps == 5
        outputs.append(output)

    for filename in ("ldscore.baseline.parquet", "ldscore.query.parquet", "ldscore.overlap.parquet"):
        pd.testing.assert_frame_equal(
            pd.read_parquet(outputs[0] / filename), pd.read_parquet(outputs[1] / filename)
        )
    sequential, parallel = [json.loads((output / "metadata.json").read_text()) for output in outputs]
    for field in ("counts", "count_config", "overlap_config", "annotation_types", "chromosomes"):
        assert sequential[field] == parallel[field]


def test_all_cores_output_matches_sequential(run_minimal_ldscore):
    # workers=-1 (all cores) exercises the affinity-resolved pool path end to end.
    seq = run_minimal_ldscore(workers=1, out_subdir="seq_allcores")
    par = run_minimal_ldscore(workers=-1, out_subdir="allcores")
    pd.testing.assert_frame_equal(
        pd.read_parquet(seq / "ldscore.baseline.parquet"),
        pd.read_parquet(par / "ldscore.baseline.parquet"),
    )


def test_parallel_run_is_deterministic(run_minimal_ldscore):
    a = run_minimal_ldscore(workers=3, out_subdir="a")
    b = run_minimal_ldscore(workers=3, out_subdir="b")
    pd.testing.assert_frame_equal(
        pd.read_parquet(a / "ldscore.baseline.parquet"),
        pd.read_parquet(b / "ldscore.baseline.parquet"),
    )


# --- Task 6: empty-intersection skip parity under the pool ------------------


def _run_with_one_empty_chrom(build_dir: Path, out_dir: Path, workers: int) -> Path:
    """Run with a panel restriction that leaves chr2 no SNPs (forces a skip).

    A restriction empties chr2 before the ``annotation_snp_grid`` step, so the
    chromosome is skipped; a grid that misses a whole chromosome fails instead.
    """
    import warnings as _warnings

    from ldsc._refpanel.adapters import RefPanelLoader
    from ldsc.config import GlobalConfig, LDScoreConfig, RefPanelConfig
    from ldsc._ldscore.calculator import LDScoreCalculator
    from ldsc._ldscore.outputs import LDScoreOutputConfig

    metadata = pd.DataFrame(
        {
            "CHR": ["1", "1", "1", "1", "2", "2"],
            "SNP": ["rs1", "rs2", "rs3", "rs4", "rs_nomatch1", "rs_nomatch2"],
            "POS": [100, 200, 300, 400, 100, 200],
            "CM": [np.nan] * 6,
        }
    )
    gcfg = GlobalConfig(snp_identifier="rsid")
    bundle = make_annotation_bundle(
        metadata=metadata,
        baseline_annotations=pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float32)}),
        query_annotations=pd.DataFrame(index=metadata.index),
        baseline_columns=["base"],
        query_columns=[],
        chromosomes=["1", "2"],
        source_summary={},
        config_snapshot=gcfg,
    )
    restriction = out_dir.parent / f"{out_dir.name}_chr1_only.tsv"
    restriction.parent.mkdir(parents=True, exist_ok=True)
    restriction.write_text("SNP\nrs1\nrs2\nrs3\nrs4\n")
    spec = RefPanelConfig(backend="r2_matrix", r2_matrix_dir=str(build_dir), ref_panel_snps_file=str(restriction))
    ref_panel = RefPanelLoader(gcfg).load(spec)
    ldcfg = LDScoreConfig(ld_wind_kb=1.0, snp_batch_size=2, whole_chromosome_ok=True, workers=workers)
    with _warnings.catch_warnings():
        _warnings.simplefilter("ignore")
        LDScoreCalculator().run(
            annotation_bundle=bundle,
            ref_panel=ref_panel,
            ldscore_config=ldcfg,
            global_config=gcfg,
            output_config=LDScoreOutputConfig(output_dir=str(out_dir)),
        )
    return out_dir


def test_skip_path_matches_under_pool(two_chrom_panel, tmp_path):
    seq = _run_with_one_empty_chrom(two_chrom_panel, tmp_path / "skip_seq", workers=1)
    par = _run_with_one_empty_chrom(two_chrom_panel, tmp_path / "skip_par", workers=3)
    seq_base = pd.read_parquet(seq / "ldscore.baseline.parquet")
    par_base = pd.read_parquet(par / "ldscore.baseline.parquet")
    # Only the matching chromosome 1 survives in both paths.
    assert seq_base["CHR"].astype(str).tolist() == ["1", "1", "1", "1"]
    pd.testing.assert_frame_equal(seq_base, par_base)


# --- --workers CLI flag -----------------------------------------------------


def test_cli_threads_reaches_config():
    from ldsc._ldscore.calculator import _ldscore_config_from_args, build_parser

    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1", "--workers", "2"])
    assert _ldscore_config_from_args(args).workers == 2


def test_cli_threads_accepts_negative():
    from ldsc._ldscore.calculator import _ldscore_config_from_args, build_parser

    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1", "--workers", "-1"])
    assert _ldscore_config_from_args(args).workers == -1


def test_cli_threads_defaults_to_one():
    from ldsc._ldscore.calculator import _ldscore_config_from_args, build_parser

    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1"])
    assert _ldscore_config_from_args(args).workers == 1


def test_retired_threads_setting_is_rejected_by_cli_and_python():
    from ldsc.errors import LDSCUsageError
    from ldsc._ldscore.calculator import build_parser, run_ldscore

    with pytest.raises(SystemExit) as caught:
        build_parser().parse_args(["--output-dir", "unused", "--threads", "2"])
    assert caught.value.code == 2
    with pytest.raises(LDSCUsageError, match="threads"):
        run_ldscore(threads=2)
