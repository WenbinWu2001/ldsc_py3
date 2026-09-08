from __future__ import annotations

import importlib.util
import json
import subprocess
import sys

import pandas as pd
import pytest


HAS_MATPLOTLIB = importlib.util.find_spec("matplotlib") is not None


def test_root_import_and_lazy_public_exports_do_not_import_matplotlib():
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; import ldsc; "
                "assert 'matplotlib' not in sys.modules; "
                "from ldsc import plot_result, convert_h2_scale; "
                "assert 'matplotlib' not in sys.modules"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr


@pytest.mark.skipif(HAS_MATPLOTLIB, reason="exercises the core-only dependency failure")
def test_plot_and_range_conversion_fail_before_output_without_matplotlib(tmp_path):
    from ldsc.errors import LDSCDependencyError
    from ldsc.h2_scale import convert_h2_scale
    from ldsc.plotting import plot_result

    result_dir = tmp_path / "h2"
    diagnostics = result_dir / "diagnostics"
    diagnostics.mkdir(parents=True)
    pd.DataFrame(
        [{"trait_name": "Trait", "total_h2_obs": 0.2, "total_h2_obs_se": 0.03}]
    ).to_csv(result_dir / "h2.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "bin": 1,
                "n_snps": 1,
                "ld_score_min": 1.0,
                "ld_score_max": 1.0,
                "mean_ld_score": 1.0,
                "mean_chi_square": 1.1,
                "sd_chi_square": float("nan"),
                "mean_sample_size": 1000,
                "mean_fitted_chi_square": 1.1,
                "mean_regression_weight": 0.5,
            }
        ]
    ).to_csv(diagnostics / "ld_score_regression_bins.tsv", sep="\t", index=False)
    (diagnostics / "metadata.json").write_text(
        json.dumps(
            {
                "artifact_type": "h2_result",
                "files": {
                    "summary": "h2.tsv",
                    "ld_score_regression_bins": "diagnostics/ld_score_regression_bins.tsv",
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(LDSCDependencyError, match=r"ldsc\[plot\]"):
        plot_result(result_dir)
    with pytest.raises(LDSCDependencyError, match=r"ldsc\[plot\]"):
        convert_h2_scale(
            result_dir,
            samp_prev=0.5,
            pop_prev_range=(0.01, 0.02),
            num_points=2,
        )

    assert not (result_dir / "plots").exists()
    assert not (result_dir / "postprocessing").exists()
