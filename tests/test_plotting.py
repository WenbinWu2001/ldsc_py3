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
