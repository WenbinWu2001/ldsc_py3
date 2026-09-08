from __future__ import annotations

import argparse
from pathlib import Path
from unittest import mock

import pytest

from ldsc import cli
from ldsc._logging import overwrite_failure_marker
from ldsc import ref_panel_builder


def test_authorized_overwrite_failure_writes_durable_marker_without_rollback(tmp_path):
    existing = tmp_path / "result.tsv"
    existing.write_text("partially replaced", encoding="utf-8")

    with pytest.raises(RuntimeError, match="publication failed"):
        with overwrite_failure_marker(
            tmp_path,
            overwrite=True,
            command="ldsc example --output-dir result --overwrite",
            log_path=tmp_path / "diagnostics" / "example.log",
        ):
            raise RuntimeError("publication failed")

    marker = tmp_path / "RUN_FAILED.txt"
    text = marker.read_text(encoding="utf-8")
    assert "publication failed" in text
    assert "incomplete or mixed" in text
    assert "not restored" in text
    assert "diagnostics/example.log" in text
    assert existing.read_text(encoding="utf-8") == "partially replaced"


def test_successful_retry_removes_marker_and_no_overwrite_failure_does_not_create_one(tmp_path):
    marker = tmp_path / "RUN_FAILED.txt"
    marker.write_text("old failure", encoding="utf-8")

    with overwrite_failure_marker(tmp_path, overwrite=True, command="retry"):
        pass
    assert not marker.exists()

    with pytest.raises(RuntimeError):
        with overwrite_failure_marker(tmp_path, overwrite=False, command="fresh run"):
            raise RuntimeError("collision")
    assert not marker.exists()


def test_cli_failure_marker_uses_nested_plot_scope(tmp_path):
    result_dir = tmp_path / "result"
    result_dir.mkdir()

    with mock.patch.object(cli, "main", side_effect=RuntimeError("plot failed")):
        status = cli.run_cli(["plot", "--result-dir", str(result_dir), "--overwrite"])

    assert status == 2
    marker = result_dir / "plots" / "RUN_FAILED.txt"
    assert marker.is_file()
    assert "ldsc plot" in marker.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("command", "source_flag", "nested"),
    [
        ("convert-h2-scale", "--h2-result-dir", "postprocessing/liability-scale"),
        ("h2", "--output-dir", ""),
        ("partitioned-h2", "--output-dir", ""),
        ("quantile-h2", "--output-dir", ""),
        ("rg", "--output-dir", ""),
        ("annotate", "--output-dir", ""),
        ("ldscore", "--output-dir", ""),
        ("build-ref-panel", "--output-dir", ""),
        ("build-gene-ldscore-index", "--output-dir", ""),
        ("convert-ldsc2-ldscores", "--output-dir", ""),
        ("munge-sumstats", "--output-dir", ""),
        ("query-r2", "--output-dir", ""),
    ],
)
def test_cli_materializing_commands_share_marker_policy(tmp_path, command, source_flag, nested):
    root = tmp_path / command
    with mock.patch.object(cli, "main", side_effect=RuntimeError("failed")):
        status = cli.run_cli([command, source_flag, str(root), "--overwrite"])

    assert status == 2
    assert (root / nested / "RUN_FAILED.txt").is_file()


def test_cli_infer_only_munge_failure_writes_no_marker(tmp_path):
    with mock.patch.object(cli, "main", side_effect=RuntimeError("failed")):
        cli.run_cli(
            [
                "munge-sumstats",
                "--output-dir",
                str(tmp_path / "out"),
                "--overwrite",
                "--infer-only",
            ]
        )

    assert not (tmp_path / "out" / "RUN_FAILED.txt").exists()


def test_unknown_cli_command_does_not_claim_an_output_scope(tmp_path):
    output_dir = tmp_path / "out"
    with mock.patch.object(cli, "main", side_effect=RuntimeError("failed")):
        cli.run_cli(["not-a-command", "--output-dir", str(output_dir), "--overwrite"])

    assert not (output_dir / "RUN_FAILED.txt").exists()


def test_cli_concrete_reference_panel_failure_uses_chromosome_marker(tmp_path):
    output_dir = tmp_path / "panel"
    with mock.patch.object(cli, "main", side_effect=RuntimeError("failed")):
        cli.run_cli(
            [
                "build-ref-panel",
                "--plink-prefix",
                str(tmp_path / "reference.chr22"),
                "--output-dir",
                str(output_dir),
                "--overwrite",
            ]
        )

    assert (output_dir / "RUN_FAILED.chr22.txt").is_file()
    assert not (output_dir / "RUN_FAILED.txt").exists()


def test_python_reference_panel_failure_uses_chromosome_marker(tmp_path):
    output_dir = tmp_path / "panel"
    args = argparse.Namespace(
        output_dir=str(output_dir),
        overwrite=True,
        plink_prefix=str(tmp_path / "reference.chr22"),
    )

    with mock.patch.object(
        ref_panel_builder,
        "config_from_args",
        side_effect=RuntimeError("failed before config"),
    ):
        with pytest.raises(RuntimeError, match="failed before config"):
            ref_panel_builder.run_build_ref_panel_from_args(args)

    assert (output_dir / "RUN_FAILED.chr22.txt").is_file()
    assert not (output_dir / "RUN_FAILED.txt").exists()
