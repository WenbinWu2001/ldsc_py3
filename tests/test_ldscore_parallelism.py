"""Tests for cross-chromosome LD-score parallelism (`num_workers`).

The fixtures here build a minimal but faithful two-chromosome canonical
index-format reference panel (the same `IDX_1/IDX_2/R2/SIGN` parquet + sidecar
that `build-ref-panel` emits, written via the real builder helpers) so the
parallel path is exercised end-to-end through `SortedR2BlockReader`, not a mock.
"""

from __future__ import annotations

import importlib.util
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


def _run_minimal_ldscore(build_dir: Path, out_dir: Path, num_workers: int) -> Path:
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
    def _runner(num_workers: int, out_subdir: str) -> Path:
        return _run_minimal_ldscore(two_chrom_panel, tmp_path / out_subdir, num_workers)

    return _runner


# --- Fixture smoke test -----------------------------------------------------


def test_fixture_builds_and_runs_sequential(run_minimal_ldscore):
    out = run_minimal_ldscore(num_workers=1, out_subdir="seq_smoke")
    baseline = pd.read_parquet(out / "ldscore.baseline.parquet")
    assert baseline["CHR"].astype(str).tolist() == ["1", "1", "1", "1", "2", "2", "2"]
    assert baseline["SNP"].tolist() == ["rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7"]
    assert "base" in baseline.columns
