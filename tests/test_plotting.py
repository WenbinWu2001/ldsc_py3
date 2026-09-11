from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

def _write_result(tmp_path, *, name, metadata, filename, rows):
    result_dir = tmp_path / name
    (result_dir / "diagnostics").mkdir(parents=True)
    pd.DataFrame(rows).to_csv(result_dir / filename, sep="\t", index=False)
    (result_dir / "diagnostics" / "metadata.json").write_text(
        json.dumps(metadata), encoding="utf-8"
    )
    return result_dir


def test_h2_dispatch_plots_saved_bins_without_refitting(tmp_path):
    from ldsc.plotting import plot_result

    result_dir = _write_result(
        tmp_path,
        name="h2",
        filename="diagnostics/ld_score_regression_bins.tsv",
        metadata={
            "artifact_type": "h2_result",
            "trait_name": "Example trait",
            "files": {"ld_score_regression_bins": "diagnostics/ld_score_regression_bins.tsv"},
        },
        rows=[
            {
                "bin": 1,
                "n_snps": 10,
                "ld_score_min": 1.0,
                "ld_score_max": 2.0,
                "mean_ld_score": 1.5,
                "mean_chi_square": 1.1,
                "sd_chi_square": 0.2,
                "mean_sample_size": 1000,
                "mean_fitted_chi_square": 1.12,
                "mean_regression_weight": 0.25,
            },
            {
                "bin": 2,
                "n_snps": 10,
                "ld_score_min": 2.0,
                "ld_score_max": 4.0,
                "mean_ld_score": 3.0,
                "mean_chi_square": 1.3,
                "sd_chi_square": 0.3,
                "mean_sample_size": 1000,
                "mean_fitted_chi_square": 1.28,
                "mean_regression_weight": 0.5,
            },
        ],
    )

    artifact = plot_result(result_dir)

    assert artifact.kind == "ld_score_regression"
    assert artifact.path == result_dir / "plots" / "ld_score_regression.png"
    assert artifact.path.is_file()
    offsets = artifact.axes.collections[0].get_offsets()
    np.testing.assert_allclose(offsets, [[1.5, 1.1], [3.0, 1.3]])
    np.testing.assert_allclose(artifact.axes.lines[0].get_xdata(), [1.5, 3.0])
    np.testing.assert_allclose(artifact.axes.lines[0].get_ydata(), [1.12, 1.28])
    metadata = json.loads((result_dir / "plots" / "diagnostics" / "metadata.json").read_text())
    assert metadata["artifact_type"] == "plot_result"
    assert metadata["plot_kind"] == "ld_score_regression"
    artifact.figure.clf()


def test_old_h2_result_without_bins_hard_fails_with_rerun_guidance(tmp_path):
    from ldsc.errors import LDSCInputError
    from ldsc.plotting import plot_result

    result_dir = tmp_path / "old-h2"
    (result_dir / "diagnostics").mkdir(parents=True)
    (result_dir / "diagnostics" / "metadata.json").write_text(
        json.dumps({"artifact_type": "h2_result", "files": {"summary": "h2.tsv"}}),
        encoding="utf-8",
    )

    with pytest.raises(LDSCInputError, match="ld_score_regression_bins.tsv.*rerun.*h2"):
        plot_result(result_dir)


@pytest.mark.parametrize(
    ("pair_kind", "filename", "kind"),
    [
        ("all_pairs", "rg_heatmap.png", "rg_heatmap"),
        ("anchor", "rg_anchor_forest.png", "rg_anchor_forest"),
    ],
)
def test_rg_dispatch_selects_pair_kind_and_never_adds_heatmap_colorbar(
    tmp_path, pair_kind, filename, kind
):
    from ldsc.plotting import plot_result

    rows = [
        {"trait_1": "Anchor", "trait_2": "Trait B", "rg": 0.42, "rg_se": 0.08},
        {"trait_1": "Anchor", "trait_2": "Trait C", "rg": -0.31, "rg_se": 0.11},
    ]
    if pair_kind == "all_pairs":
        rows.append({"trait_1": "Trait B", "trait_2": "Trait C", "rg": 0.15, "rg_se": 0.07})
    result_dir = _write_result(
        tmp_path,
        name=f"rg-{pair_kind}",
        filename="rg.tsv",
        metadata={
            "artifact_type": "rg_result",
            "pair_kind": pair_kind,
            "trait_names": ["Anchor", "Trait B", "Trait C"],
            "files": {"rg": "rg.tsv"},
        },
        rows=rows,
    )

    artifact = plot_result(result_dir, output_dir=tmp_path / "custom-plots")

    assert artifact.kind == kind
    assert artifact.path.name == filename
    assert artifact.path.is_file()
    if pair_kind == "all_pairs":
        assert len(artifact.figure.axes) == 1
        assert "0.42\n(0.08)" in {text.get_text() for text in artifact.axes.texts}
    else:
        assert [tick.get_text() for tick in artifact.axes.get_yticklabels()] == ["Trait B", "Trait C"]
        assert "0.42 (0.08)" in {text.get_text() for text in artifact.axes.texts}
    artifact.figure.clf()


def _write_rg_with_h2(tmp_path, *, pair_kind="all_pairs"):
    result_dir = _write_result(
        tmp_path,
        name="rg-with-h2",
        filename="rg.tsv",
        metadata={
            "artifact_type": "rg_result",
            "pair_kind": pair_kind,
            "trait_names": ["Anchor", "Trait B", "Trait C"],
            "files": {"rg": "rg.tsv", "h2_per_trait": "heritabilities.tsv"},
        },
        rows=[
            {"trait_1": "Anchor", "trait_2": "Trait B", "rg": 0.42, "rg_se": 0.08},
            {"trait_1": "Anchor", "trait_2": "Trait C", "rg": -0.31, "rg_se": 0.11},
        ],
    )
    pd.DataFrame({
        "trait_name": ["Trait C", "Anchor", "Unused", "Trait B"],
        "total_h2_obs": [1.23, 0.25, 0.99, -0.12],
        "total_h2_obs_se": [0.04, 0.03, 0.01, 0.0],
        "total_h2_liab": [2.46, 0.50, 1.98, -0.24],
        "total_h2_liab_se": [0.08, 0.06, 0.02, 0.0],
    }).to_csv(result_dir / "heritabilities.tsv", sep="\t", index=False)
    return result_dir


def test_rg_diagonal_uses_observed_single_trait_h2_in_metadata_order(tmp_path):
    from matplotlib import pyplot as plt
    from matplotlib.colors import to_rgba
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path)
    artifact = plot_result(result_dir)
    try:
        labels = {text.get_position(): text.get_text() for text in artifact.axes.texts}
        assert labels[(0.5, 0.5)] == "0.25\n(0.03)"
        assert labels[(1.5, 1.5)] == "−0.12\n(0.00)"
        assert labels[(2.5, 2.5)] == "1.23\n(0.04)"
        assert labels[(0.5, 1.5)] == "0.42\n(0.08)"
        assert labels[(0.5, 2.5)] == "−0.31\n(0.11)"
        assert labels[(1.5, 2.5)] == "failed"
        assert all(x <= y for x, y in labels)
        assert len(artifact.figure.axes) == 1
        diagonal = [patch for patch in artifact.axes.patches if patch.get_x() == patch.get_y()]
        assert len(diagonal) == 3
        assert all(patch.get_facecolor() == to_rgba("#E0E0E0") for patch in diagonal)
        assert "observed" in artifact.axes.get_title().lower()
        metadata = json.loads((result_dir / "plots/diagnostics/metadata.json").read_text())
        assert metadata["source_table"] == "rg.tsv"
        assert metadata["heritability_annotations"] == {
            "source_table": "heritabilities.tsv",
            "source_available": True,
            "scale": "observed",
            "estimate_column": "total_h2_obs",
            "se_column": "total_h2_obs_se",
            "uncertainty": "block_jackknife_standard_error",
            "missing_label": "failed",
        }
    finally:
        plt.close(artifact.figure)


def test_anchor_h2_column_and_subtitle_preserve_correlations(tmp_path):
    from matplotlib import pyplot as plt
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path, pair_kind="anchor")
    artifact = plot_result(result_dir)
    try:
        assert artifact.path.name == "rg_anchor_forest.png"
        assert [tick.get_text() for tick in artifact.axes.get_yticklabels()] == ["Trait B", "Trait C"]
        texts = {text.get_text(): text for text in artifact.axes.texts}
        assert "0.25 (0.03)" in artifact.axes.get_title()
        assert "Observed" in artifact.axes.get_title()
        assert texts["−0.12 (0.00)"].get_position()[1] == 0
        assert texts["1.23 (0.04)"].get_position()[1] == 1
        assert "0.42 (0.08)" in texts
        assert "−0.31 (0.11)" in texts
        intervals = [container for container in artifact.axes.containers if hasattr(container, "has_xerr")]
        np.testing.assert_allclose(intervals[0].lines[2][0].get_segments(), [[[0.34, 0], [0.50, 0]]])
        np.testing.assert_allclose(intervals[1].lines[2][0].get_segments(), [[[-0.42, 1], [-0.20, 1]]])
        artifact.figure.canvas.draw()
        renderer = artifact.figure.canvas.get_renderer()
        for tick, label in zip(artifact.axes.get_yticklabels(), ["−0.12 (0.00)", "1.23 (0.04)"]):
            h2_box = texts[label].get_window_extent(renderer)
            assert tick.get_window_extent(renderer).x1 < h2_box.x0
            assert h2_box.x1 < artifact.axes.get_window_extent(renderer).x0
    finally:
        plt.close(artifact.figure)


@pytest.mark.parametrize("pair_kind", ["all_pairs", "anchor"])
@pytest.mark.parametrize("missing", ["declaration", "file", "rows"])
def test_rg_missing_heritabilities_are_failed(tmp_path, pair_kind, missing):
    from matplotlib import pyplot as plt
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path, pair_kind=pair_kind)
    h2_path = result_dir / "heritabilities.tsv"
    if missing == "declaration":
        path = result_dir / "diagnostics/metadata.json"
        metadata = json.loads(path.read_text())
        del metadata["files"]["h2_per_trait"]
        path.write_text(json.dumps(metadata))
    elif missing == "file":
        h2_path.unlink()
    else:
        pd.read_csv(h2_path, sep="\t").iloc[:0].to_csv(h2_path, sep="\t", index=False)
    artifact = plot_result(result_dir)
    try:
        failed = [text for text in artifact.axes.texts if text.get_text() == "failed"]
        if pair_kind == "all_pairs":
            assert {text.get_position() for text in failed} >= {(0.5, 0.5), (1.5, 1.5), (2.5, 2.5)}
        else:
            assert "failed" in artifact.axes.get_title()
            assert len(failed) == 2
        metadata = json.loads((result_dir / "plots/diagnostics/metadata.json").read_text())
        annotation = metadata["heritability_annotations"]
        assert annotation["source_table"] == (None if missing == "declaration" else "heritabilities.tsv")
        assert annotation["source_available"] == (missing == "rows")
    finally:
        plt.close(artifact.figure)


def test_anchor_failed_last_row_stays_above_axis_spine(tmp_path):
    from matplotlib import pyplot as plt
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path, pair_kind="anchor")
    path = result_dir / "rg.tsv"
    table = pd.read_csv(path, sep="\t")
    table.loc[1, ["rg", "rg_se"]] = np.nan
    table.to_csv(path, sep="\t", index=False)
    artifact = plot_result(result_dir)
    try:
        bottom, top = artifact.axes.get_ylim()
        assert bottom > 1.2
        assert top < -0.2
        assert "failed" in {text.get_text() for text in artifact.axes.texts}
    finally:
        plt.close(artifact.figure)


@pytest.mark.parametrize("pair_kind", ["all_pairs", "anchor"])
def test_rg_unusable_h2_pairs_are_failed_without_losing_valid_traits(tmp_path, pair_kind):
    from matplotlib import pyplot as plt
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path, pair_kind=pair_kind)
    invalid = [(None, 0.1), (float("inf"), 0.1), ("bad", 0.1), (0.2, None),
               (0.2, float("inf")), (0.2, "bad"), (0.2, -0.1)]
    names = ["Anchor", "Trait B", "Trait C"] + [f"Invalid {i}" for i in range(len(invalid))]
    path = result_dir / "diagnostics/metadata.json"
    metadata = json.loads(path.read_text())
    metadata["trait_names"] = names
    path.write_text(json.dumps(metadata))
    rows = [{"trait_name": "Anchor", "total_h2_obs": 0.25, "total_h2_obs_se": 0.03}]
    rows.extend({"trait_name": name, "total_h2_obs": estimate, "total_h2_obs_se": se}
                for name, (estimate, se) in zip(names[3:], invalid))
    pd.DataFrame(rows).to_csv(result_dir / "heritabilities.tsv", sep="\t", index=False)
    artifact = plot_result(result_dir)
    try:
        if pair_kind == "all_pairs":
            labels = {text.get_position(): text.get_text() for text in artifact.axes.texts}
            assert labels[(0.5, 0.5)] == "0.25\n(0.03)"
            assert all(labels[(i + 0.5, i + 0.5)] == "failed" for i in range(1, len(names)))
        else:
            assert "0.25 (0.03)" in artifact.axes.get_title()
            # Column entries lie to the left of the scientific plotting area.
            artifact.figure.canvas.draw()
            renderer = artifact.figure.canvas.get_renderer()
            column = [text for text in artifact.axes.texts
                      if text.get_text() == "failed" and text.get_window_extent(renderer).x1
                      < artifact.axes.get_window_extent(renderer).x0]
            assert len(column) == len(names) - 1
    finally:
        plt.close(artifact.figure)


@pytest.mark.parametrize("pair_kind", ["all_pairs", "anchor"])
@pytest.mark.parametrize("problem", ["columns", "duplicate", "empty", "parser", "extra_field", "encoding", "directory"])
def test_rg_rejects_malformed_heritability_sources_before_saving(tmp_path, pair_kind, problem):
    from ldsc.errors import LDSCInputError
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path, pair_kind=pair_kind)
    path = result_dir / "heritabilities.tsv"
    header = "trait_name\ttotal_h2_obs\ttotal_h2_obs_se\n"
    if problem == "columns":
        path.write_text("trait_name\ttotal_h2_obs\nAnchor\t0.25\n")
    elif problem == "duplicate":
        path.write_text(header + "Anchor\t0.25\t0.03\nAnchor\t0.30\t0.04\n")
    elif problem == "empty":
        path.write_text("")
    elif problem == "parser":
        path.write_text(header + '"Anchor\t0.25\t0.03\n')
    elif problem == "extra_field":
        path.write_text(header + "Anchor\t0.25\t0.03\textra\n")
    elif problem == "encoding":
        path.write_bytes(b"\xff\xfe\xff")
    else:
        path.unlink()
        path.mkdir()
    with pytest.raises(LDSCInputError, match="heritabilit|h2_per_trait"):
        plot_result(result_dir)
    assert not (result_dir / "plots").exists()


@pytest.mark.parametrize("path_kind", ["absolute", "traversal", "symlink", "invalid"])
def test_rg_rejects_unsafe_h2_declarations_even_if_target_is_missing(tmp_path, path_kind):
    from ldsc.errors import LDSCInputError
    from ldsc.plotting import plot_result

    result_dir = _write_rg_with_h2(tmp_path)
    if path_kind == "absolute":
        declared = str(tmp_path / "missing.tsv")
    elif path_kind == "traversal":
        declared = "../missing.tsv"
    elif path_kind == "symlink":
        (result_dir / "linked.tsv").symlink_to(tmp_path / "missing.tsv")
        declared = "linked.tsv"
    else:
        declared = 42
    path = result_dir / "diagnostics/metadata.json"
    metadata = json.loads(path.read_text())
    metadata["files"]["h2_per_trait"] = declared
    path.write_text(json.dumps(metadata))
    with pytest.raises(LDSCInputError, match="files.h2_per_trait"):
        plot_result(result_dir)
    assert not (result_dir / "plots").exists()


def test_functional_and_quantile_enrichment_use_horizontal_bars_with_null_behind(tmp_path):
    from ldsc.plotting import plot_result

    functional = _write_result(
        tmp_path,
        name="functional",
        filename="partitioned_h2.tsv",
        metadata={
            "artifact_type": "partitioned_h2_result",
            "analysis_type": "functional_category",
            "headline_metric": "enrichment",
            "files": {"summary": "partitioned_h2.tsv"},
        },
        rows=[
            {"category": "Conserved", "enrichment": 2.5, "enrichment_se": 0.4},
            {"category": "Enhancer", "enrichment": 1.35, "enrichment_se": 0.2},
        ],
    )
    quantile = _write_result(
        tmp_path,
        name="quantile",
        filename="quantile_h2.tsv",
        metadata={
            "artifact_type": "quantile_h2_result",
            "target_annotation": "Chromatin score",
            "files": {"quantile_h2": "quantile_h2.tsv"},
        },
        rows=[
            {
                "quantile": 1,
                "target_value_lower": 0,
                "target_value_upper": 0.5,
                "enrichment": 0.6,
                "enrichment_se": 0.15,
            },
            {
                "quantile": 2,
                "target_value_lower": 0.5,
                "target_value_upper": 1,
                "enrichment": 1.75,
                "enrichment_se": 0.25,
            },
        ],
    )

    functional_artifact = plot_result(functional)
    quantile_artifact = plot_result(quantile)

    np.testing.assert_allclose([patch.get_width() for patch in functional_artifact.axes.patches], [2.5, 1.35])
    np.testing.assert_allclose([patch.get_width() for patch in quantile_artifact.axes.patches], [0.6, 1.75])
    for artifact in (functional_artifact, quantile_artifact):
        null_line = next(line for line in artifact.axes.lines if np.allclose(line.get_xdata(), [1, 1]))
        assert null_line.get_zorder() < min(patch.get_zorder() for patch in artifact.axes.patches)
    colors = np.asarray([patch.get_facecolor()[:3] for patch in quantile_artifact.axes.patches])
    luminance = colors @ np.asarray([0.2126, 0.7152, 0.0722])
    assert luminance[1] < luminance[0]
    functional_artifact.figure.clf()
    quantile_artifact.figure.clf()


def test_cell_type_summary_is_horizontal_nominal_pvalue_scatter(tmp_path):
    from ldsc.plotting import plot_result

    result_dir = _write_result(
        tmp_path,
        name="queries",
        filename="partitioned_h2.tsv",
        metadata={
            "artifact_type": "partitioned_h2_result",
            "analysis_type": "cell_type_specific",
            "headline_metric": "coefficient",
            "files": {"summary": "partitioned_h2.tsv"},
        },
        rows=[
            {"category": "Astrocytes", "coefficient_p": 0.2, "coefficient_z": 0.84},
            {"category": "Microglia", "coefficient_p": 1e-5, "coefficient_z": 4.26},
            {"category": "Extreme", "coefficient_p": 0.0, "coefficient_z": 40.0},
        ],
    )

    artifact = plot_result(result_dir)

    assert artifact.kind == "cell_type_query_pvalues"
    assert [tick.get_text() for tick in artifact.axes.get_yticklabels()] == [
        "Extreme",
        "Microglia",
        "Astrocytes",
    ]
    offsets = artifact.axes.collections[0].get_offsets()
    assert float(offsets[0, 0]) == pytest.approx(349.4370064593458)
    np.testing.assert_allclose(np.asarray(offsets[1:, 0], dtype=float), [5.0, -np.log10(0.2)])
    assert "separate fit" in artifact.axes.get_title()
    artifact.figure.clf()


def test_plot_cli_writes_default_destination_and_rejects_output_dir(tmp_path):
    from ldsc import cli

    result_dir = _write_result(
        tmp_path,
        name="functional-cli",
        filename="partitioned_h2.tsv",
        metadata={
            "artifact_type": "partitioned_h2_result",
            "analysis_type": "functional_category",
            "headline_metric": "enrichment",
            "files": {"summary": "partitioned_h2.tsv"},
        },
        rows=[{"category": "Conserved", "enrichment": 2.5, "enrichment_se": 0.4}],
    )

    artifact = cli.main(["plot", "--result-dir", str(result_dir)])
    assert artifact.path == result_dir / "plots" / "functional_h2_enrichment.png"
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["plot", "--result-dir", str(result_dir), "--output-dir", str(tmp_path / "bad")]
        )
