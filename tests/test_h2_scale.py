from __future__ import annotations

import builtins
import json
import sys

import numpy as np
import pandas as pd
import pytest

from ldsc._kernel.regression import liability_conversion_factor
from ldsc.errors import LDSCInputError

def _write_h2_result(tmp_path):
    result_dir = tmp_path / "h2-result"
    diagnostics = result_dir / "diagnostics"
    diagnostics.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "trait_name": "Example disease",
                "total_h2_obs": 0.2,
                "total_h2_obs_se": 0.03,
                "total_h2_liab": 0.9,
                "total_h2_liab_se": 0.4,
            }
        ]
    ).to_csv(result_dir / "h2.tsv", sep="\t", index=False)
    (diagnostics / "metadata.json").write_text(
        json.dumps(
            {
                "artifact_type": "h2_result",
                "files": {"summary": "h2.tsv"},
            }
        ),
        encoding="utf-8",
    )
    return result_dir


def test_exact_conversion_uses_observed_fields_and_default_destination(tmp_path, monkeypatch):
    from ldsc.h2_scale import convert_h2_scale

    result_dir = _write_h2_result(tmp_path)
    original_import = builtins.__import__

    def reject_matplotlib(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise AssertionError("exact conversion imported Matplotlib")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", reject_matplotlib)
    artifact = convert_h2_scale(result_dir, samp_prev=0.5, pop_prev=0.01)

    expected_dir = result_dir / "postprocessing" / "liability-scale"
    assert artifact.table_path == expected_dir / "h2_scale_conversion.tsv"
    assert artifact.plot_path is None
    table = pd.read_csv(artifact.table_path, sep="\t")
    factor = float(liability_conversion_factor(0.5, 0.01))
    assert table.loc[0, "conversion_factor"] == pytest.approx(factor)
    assert table.loc[0, "total_h2_liab"] == pytest.approx(0.2 * factor)
    assert table.loc[0, "total_h2_liab_se"] == pytest.approx(0.03 * factor)
    metadata = json.loads(artifact.metadata_path.read_text(encoding="utf-8"))
    assert metadata["artifact_type"] == "h2_scale_conversion_result"
    assert metadata["source_artifact_type"] == "h2_result"
    assert metadata["mode"] == "exact"


def test_range_conversion_is_inclusive_and_uses_python_output_override(tmp_path):
    from ldsc.h2_scale import convert_h2_scale

    result_dir = _write_h2_result(tmp_path)
    output_dir = tmp_path / "custom"
    artifact = convert_h2_scale(
        result_dir,
        output_dir=output_dir,
        samp_prev=0.5,
        pop_prev_range=(0.01, 0.03),
        num_points=3,
    )

    assert artifact.plot_path == output_dir / "h2_prevalence_sensitivity.png"
    assert artifact.plot_path.is_file()
    table = pd.read_csv(artifact.table_path, sep="\t")
    np.testing.assert_allclose(table["pop_prev"], [0.01, 0.02, 0.03])
    np.testing.assert_allclose(
        table["conversion_factor"],
        liability_conversion_factor(0.5, np.array([0.01, 0.02, 0.03])),
    )


def test_exact_overwrite_removes_stale_sensitivity_plot(tmp_path):
    from ldsc.h2_scale import convert_h2_scale

    result_dir = _write_h2_result(tmp_path)
    ranged = convert_h2_scale(
        result_dir,
        samp_prev=0.5,
        pop_prev_range=(0.01, 0.03),
        num_points=3,
    )
    assert ranged.plot_path.is_file()

    exact = convert_h2_scale(
        result_dir,
        samp_prev=0.5,
        pop_prev=0.02,
        overwrite=True,
    )

    assert exact.plot_path is None
    assert not ranged.plot_path.exists()


def test_conversion_requires_declared_summary_file(tmp_path):
    from ldsc.h2_scale import convert_h2_scale

    result_dir = _write_h2_result(tmp_path)
    metadata_path = result_dir / "diagnostics" / "metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["files"] = {}
    metadata_path.write_text(json.dumps(metadata), encoding="utf-8")

    with pytest.raises(LDSCInputError, match="files.summary"):
        convert_h2_scale(result_dir, samp_prev=0.5, pop_prev=0.01)


def test_conversion_cli_has_no_output_dir_and_writes_default_family(tmp_path):
    from ldsc import cli

    result_dir = _write_h2_result(tmp_path)
    result = cli.main(
        [
            "convert-h2-scale",
            "--h2-result-dir",
            str(result_dir),
            "--samp-prev",
            "0.5",
            "--pop-prev",
            "0.01",
        ]
    )

    assert result.table_path.is_file()
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "convert-h2-scale",
                "--h2-result-dir",
                str(result_dir),
                "--samp-prev",
                "0.5",
                "--pop-prev",
                "0.01",
                "--output-dir",
                str(tmp_path / "bad"),
            ]
        )
