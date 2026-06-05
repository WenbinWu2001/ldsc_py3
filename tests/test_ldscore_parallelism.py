"""Tests for cross-chromosome LD-score parallelism (`num_workers`).

The fixtures here build a minimal but faithful two-chromosome canonical
index-format reference panel (the same `IDX_1/IDX_2/R2/SIGN` parquet + sidecar
that `build-ref-panel` emits, written via the real builder helpers) so the
parallel path is exercised end-to-end through `SortedR2BlockReader`, not a mock.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_HAS_PYARROW = importlib.util.find_spec("pyarrow") is not None
pytestmark = pytest.mark.skipif(not _HAS_PYARROW, reason="pyarrow is required for index parquet panels")


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
    from ldsc._kernel import ref_panel_builder as kb

    build_dir.mkdir(parents=True, exist_ok=True)
    for chrom, (panel, pairs) in chroms.items():
        positions = panel["POS"].to_numpy(dtype=np.int64)
        runtime_meta = kb.build_runtime_metadata_table(metadata=panel, positions=positions, cm_values=None)
        identity_hash = kb.sidecar_identity_sha256(runtime_meta)
        i = np.array([p[0] for p in pairs], dtype=np.int64)
        j = np.array([p[1] for p in pairs], dtype=np.int64)
        r2 = np.array([p[2] for p in pairs], dtype=np.float32)
        sign = np.ones(len(i), dtype=np.int8)
        kb.write_r2_parquet(
            pair_chunks=[(i, j, r2, sign)],
            path=build_dir / f"chr{chrom}_r2.parquet",
            genome_build=genome_build,
            n_samples=n_samples,
            snp_identifier=snp_identifier,
            n_snps=len(panel),
            sidecar_identity_sha256=identity_hash,
        )
        kb.write_runtime_metadata_sidecar(
            runtime_meta,
            build_dir / f"chr{chrom}_meta.tsv.gz",
            genome_build=genome_build,
            snp_identifier=snp_identifier,
        )


def _run_minimal_ldscore(build_dir: Path, out_dir: Path, num_workers: int, parallel: bool = True) -> Path:
    """Run a synthetic-base LD-score job over the index panel; return out_dir."""
    from ldsc.config import GlobalConfig, set_global_config
    from ldsc.ldscore_calculator import run_ldscore

    set_global_config(GlobalConfig(snp_identifier="rsid"))
    run_ldscore(
        r2_dir=str(build_dir),
        output_dir=str(out_dir),
        ld_wind_kb=1.0,
        snp_batch_size=2,
        yes_really=True,
        parallel=parallel,
        num_workers=num_workers,
    )
    return out_dir


@pytest.fixture
def two_chrom_panel(tmp_path: Path) -> Path:
    build_dir = tmp_path / "panel" / "hg19"
    _write_index_panel(build_dir, _PANEL_CHROMS)
    return build_dir


@pytest.fixture
def run_minimal_ldscore(two_chrom_panel: Path, tmp_path: Path):
    def _runner(num_workers: int, out_subdir: str, parallel: bool = True) -> Path:
        return _run_minimal_ldscore(two_chrom_panel, tmp_path / out_subdir, num_workers, parallel=parallel)

    return _runner


# --- Fixture smoke test -----------------------------------------------------


def test_fixture_builds_and_runs_sequential(run_minimal_ldscore):
    out = run_minimal_ldscore(num_workers=1, out_subdir="seq_smoke")
    baseline = pd.read_parquet(out / "ldscore.baseline.parquet")
    assert baseline["CHR"].astype(str).tolist() == ["1", "1", "1", "1", "2", "2", "2"]
    assert baseline["SNP"].tolist() == ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7"]
    assert "base" in baseline.columns


# --- Task 1: LDScoreConfig.num_workers --------------------------------------


def _cfg(**kw):
    from ldsc.config import LDScoreConfig

    base = dict(ld_wind_cm=1.0)
    base.update(kw)
    return LDScoreConfig(**base)


def test_parallel_defaults_to_true():
    assert _cfg().parallel is True


def test_num_workers_defaults_to_zero_auto():
    assert _cfg().num_workers == 0


def test_num_workers_negative_rejected():
    from ldsc.errors import LDSCConfigError

    with pytest.raises(LDSCConfigError):
        _cfg(num_workers=-1)


def test_num_workers_explicit_value_preserved():
    assert _cfg(num_workers=4).num_workers == 4


# --- Task 2: worker-count resolution ----------------------------------------


def test_resolve_worker_count():
    from ldsc.ldscore_calculator import _resolve_worker_count

    assert _resolve_worker_count(1, n_chromosomes=22) == 1
    assert _resolve_worker_count(8, n_chromosomes=3) == 3
    assert _resolve_worker_count(4, n_chromosomes=22) == 4
    auto = _resolve_worker_count(0, n_chromosomes=2)
    assert 1 <= auto <= 2
    assert _resolve_worker_count(8, n_chromosomes=1) == 1
    assert _resolve_worker_count(8, n_chromosomes=0) == 1


def test_resolve_worker_count_parallel_disabled_forces_sequential():
    from ldsc.ldscore_calculator import _resolve_worker_count

    # parallel=False forces sequential regardless of num_workers (explicit or auto).
    assert _resolve_worker_count(8, n_chromosomes=22, parallel=False) == 1
    assert _resolve_worker_count(0, n_chromosomes=22, parallel=False) == 1
    # parallel=True with auto (0) requests up to the chromosome count.
    assert _resolve_worker_count(0, n_chromosomes=22, parallel=True) >= 1


# --- Task 3: module-level chromosome worker ---------------------------------


def _chrom1_inputs(build_dir: Path):
    """Build (chrom, bundle, spec, ldscore_config, global_config, regression) for chr1."""
    from ldsc.annotation_builder import AnnotationBundle
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
    bundle = AnnotationBundle(
        metadata=metadata,
        baseline_annotations=pd.DataFrame({"base": np.ones(len(panel), dtype=np.float32)}),
        query_annotations=pd.DataFrame(index=metadata.index),
        baseline_columns=["base"],
        query_columns=[],
        chromosomes=["1"],
        source_summary={},
        config_snapshot=GlobalConfig(snp_identifier="rsid"),
    )
    spec = RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir))
    ldscore_config = LDScoreConfig(ld_wind_kb=1.0, snp_batch_size=2, whole_chromosome_ok=True)
    return "1", bundle, spec, ldscore_config, GlobalConfig(snp_identifier="rsid"), None


def test_compute_one_chromosome_returns_success_outcome(two_chrom_panel):
    from ldsc.ldscore_calculator import _compute_one_chromosome

    chrom, bundle, spec, ldcfg, gcfg, reg = _chrom1_inputs(two_chrom_panel)
    outcome = _compute_one_chromosome(chrom, bundle, spec, ldcfg, gcfg, regression_snps=reg)
    assert outcome.chrom == "1"
    assert outcome.skipped is False
    assert outcome.result is not None
    assert outcome.result.baseline_table["SNP"].tolist() == ["rs1", "rs2", "rs3", "rs4"]


# --- Task 4: pool initializer -----------------------------------------------


def test_init_worker_sets_state_and_blas_env(monkeypatch):
    from ldsc.ldscore_calculator import _WORKER_STATE, _init_worker

    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.delenv("OPENBLAS_NUM_THREADS", raising=False)
    monkeypatch.delenv("MKL_NUM_THREADS", raising=False)
    _init_worker(regression_snps={"rs1", "rs2"}, log_level="WARNING")
    assert _WORKER_STATE["regression_snps"] == {"rs1", "rs2"}
    assert os.environ["OMP_NUM_THREADS"] == "1"
    assert os.environ["OPENBLAS_NUM_THREADS"] == "1"
    assert os.environ["MKL_NUM_THREADS"] == "1"


def test_init_worker_respects_user_blas_env(monkeypatch):
    from ldsc.ldscore_calculator import _init_worker

    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    _init_worker(regression_snps=None, log_level="WARNING")
    assert os.environ["OMP_NUM_THREADS"] == "4"


# --- Task 5: sequential-vs-pool equivalence ---------------------------------


def test_run_chromosomes_dispatch_exists():
    # Guards that the sequential-vs-pool dispatch seam exists, so the equivalence
    # test below actually exercises the pool path rather than passing trivially.
    from ldsc.ldscore_calculator import LDScoreCalculator

    assert hasattr(LDScoreCalculator, "_run_chromosomes")


def test_pool_path_executes_for_multiple_workers(run_minimal_ldscore, monkeypatch):
    # Confirm a ProcessPoolExecutor is actually constructed when worker_count>1.
    import ldsc.ldscore_calculator as mod

    calls = {"n": 0}
    real_executor = mod.ProcessPoolExecutor

    def _spy(*args, **kwargs):
        calls["n"] += 1
        return real_executor(*args, **kwargs)

    monkeypatch.setattr(mod, "ProcessPoolExecutor", _spy)
    run_minimal_ldscore(num_workers=2, out_subdir="pool_spy")
    assert calls["n"] == 1


def test_parallel_output_matches_sequential(run_minimal_ldscore):
    seq = run_minimal_ldscore(num_workers=1, out_subdir="seq")
    par = run_minimal_ldscore(num_workers=3, out_subdir="par")
    seq_base = pd.read_parquet(seq / "ldscore.baseline.parquet")
    par_base = pd.read_parquet(par / "ldscore.baseline.parquet")
    pd.testing.assert_frame_equal(seq_base, par_base)


def test_parallel_run_is_deterministic(run_minimal_ldscore):
    a = run_minimal_ldscore(num_workers=3, out_subdir="a")
    b = run_minimal_ldscore(num_workers=3, out_subdir="b")
    pd.testing.assert_frame_equal(
        pd.read_parquet(a / "ldscore.baseline.parquet"),
        pd.read_parquet(b / "ldscore.baseline.parquet"),
    )


# --- Task 6: empty-intersection skip parity under the pool ------------------


def _run_with_one_empty_chrom(build_dir: Path, out_dir: Path, num_workers: int) -> Path:
    """Run with a bundle whose chr2 SNPs do not match the panel (forces a skip)."""
    import warnings as _warnings

    from ldsc._kernel.ref_panel import RefPanelLoader
    from ldsc.annotation_builder import AnnotationBundle
    from ldsc.config import GlobalConfig, LDScoreConfig, RefPanelConfig
    from ldsc.ldscore_calculator import LDScoreCalculator
    from ldsc.outputs import LDScoreOutputConfig

    metadata = pd.DataFrame(
        {
            "CHR": ["1", "1", "1", "1", "2", "2"],
            "SNP": ["rs1", "rs2", "rs3", "rs4", "rs_nomatch1", "rs_nomatch2"],
            "POS": [100, 200, 300, 400, 100, 200],
            "CM": [np.nan] * 6,
        }
    )
    gcfg = GlobalConfig(snp_identifier="rsid")
    bundle = AnnotationBundle(
        metadata=metadata,
        baseline_annotations=pd.DataFrame({"base": np.ones(len(metadata), dtype=np.float32)}),
        query_annotations=pd.DataFrame(index=metadata.index),
        baseline_columns=["base"],
        query_columns=[],
        chromosomes=["1", "2"],
        source_summary={},
        config_snapshot=gcfg,
    )
    spec = RefPanelConfig(backend="parquet_r2", r2_dir=str(build_dir))
    ref_panel = RefPanelLoader(gcfg).load(spec)
    ldcfg = LDScoreConfig(ld_wind_kb=1.0, snp_batch_size=2, whole_chromosome_ok=True, num_workers=num_workers)
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
    seq = _run_with_one_empty_chrom(two_chrom_panel, tmp_path / "skip_seq", num_workers=1)
    par = _run_with_one_empty_chrom(two_chrom_panel, tmp_path / "skip_par", num_workers=3)
    seq_base = pd.read_parquet(seq / "ldscore.baseline.parquet")
    par_base = pd.read_parquet(par / "ldscore.baseline.parquet")
    # Only the matching chromosome 1 survives in both paths.
    assert seq_base["CHR"].astype(str).tolist() == ["1", "1", "1", "1"]
    pd.testing.assert_frame_equal(seq_base, par_base)


# --- Task 7: --num-workers CLI flag -----------------------------------------


def test_cli_num_workers_reaches_config():
    from ldsc.ldscore_calculator import _ldscore_config_from_args, build_parser

    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1", "--num-workers", "2"])
    assert _ldscore_config_from_args(args).num_workers == 2


def test_cli_defaults_are_parallel_auto():
    from ldsc.ldscore_calculator import _ldscore_config_from_args, build_parser

    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1"])
    config = _ldscore_config_from_args(args)
    assert config.parallel is True
    assert config.num_workers == 0


def test_cli_no_parallel_disables():
    from ldsc.ldscore_calculator import _ldscore_config_from_args, build_parser

    parser = build_parser()
    args = parser.parse_args(["--output-dir", "x", "--ld-wind-cm", "1", "--no-parallel"])
    assert _ldscore_config_from_args(args).parallel is False


# --- Binary parallel switch: end-to-end behavior ----------------------------


def test_no_parallel_output_matches_parallel(run_minimal_ldscore):
    seq = run_minimal_ldscore(num_workers=0, out_subdir="noparallel", parallel=False)
    par = run_minimal_ldscore(num_workers=3, out_subdir="parallel", parallel=True)
    pd.testing.assert_frame_equal(
        pd.read_parquet(seq / "ldscore.baseline.parquet"),
        pd.read_parquet(par / "ldscore.baseline.parquet"),
    )


def test_no_parallel_does_not_construct_pool(run_minimal_ldscore, monkeypatch):
    import ldsc.ldscore_calculator as mod

    calls = {"n": 0}
    real_executor = mod.ProcessPoolExecutor

    def _spy(*args, **kwargs):
        calls["n"] += 1
        return real_executor(*args, **kwargs)

    monkeypatch.setattr(mod, "ProcessPoolExecutor", _spy)
    # parallel=False must stay sequential even though num_workers would request a pool.
    run_minimal_ldscore(num_workers=4, out_subdir="noparallel_spy", parallel=False)
    assert calls["n"] == 0
