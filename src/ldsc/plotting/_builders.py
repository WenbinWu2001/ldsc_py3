"""Private Matplotlib builders for the LDSC plot dispatcher."""

from __future__ import annotations

import math
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from scipy.special import log_ndtr

from ..errors import LDSCInputError


_RG_CMAP = LinearSegmentedColormap.from_list(
    "ldsc_rg",
    [
        (0.138500, 0.413312, 0.740520),
        (0.609341, 0.682122, 0.794230),
        (0.980600, 0.961552, 0.958131),
        (0.827403, 0.592091, 0.581679),
        (0.660807, 0.215267, 0.230695),
    ],
)


def build_plot(kind: str, table: pd.DataFrame, metadata: dict[str, Any]):
    """Build one selected plot from its already loaded source table."""
    if kind == "ld_score_regression":
        return _plot_h2_bins(table, trait_name=str(metadata.get("trait_name") or "Trait"))
    if kind == "rg_heatmap":
        return _plot_rg_heatmap(table, _trait_order(metadata, table))
    if kind == "rg_anchor_forest":
        return _plot_rg_anchor_forest(table, _trait_order(metadata, table))
    if kind == "functional_h2_enrichment":
        _require_columns(table, {"category", "enrichment", "enrichment_se"}, label="partitioned_h2.tsv")
        return _plot_horizontal_enrichment_bars(
            table["category"].astype(str).tolist(),
            pd.to_numeric(table["enrichment"], errors="coerce"),
            pd.to_numeric(table["enrichment_se"], errors="coerce"),
            title="Functional heritability enrichment",
            y_label="Functional annotation",
        )
    if kind == "cell_type_query_pvalues":
        return _plot_query_pvalues(table)
    if kind == "continuous_annotation_quantile_enrichment":
        return _plot_quantile_profile(table, target_annotation=metadata["target_annotation"])
    raise LDSCInputError(f"Unknown internal plot kind {kind!r}.")


def close_figure(figure) -> None:
    """Close a figure owned by a completed CLI path."""
    plt.close(figure)


def _plot_h2_bins(table: pd.DataFrame, *, trait_name: str):
    required = {
        "bin",
        "mean_ld_score",
        "mean_chi_square",
        "mean_fitted_chi_square",
        "mean_regression_weight",
    }
    _require_columns(table, required, label="ld_score_regression_bins.tsv")
    if table.empty:
        raise LDSCInputError("ld_score_regression_bins.tsv contains no bins.")
    working = table.copy()
    for column in required:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    if working["bin"].isna().any() or working["bin"].duplicated().any():
        raise LDSCInputError("ld_score_regression_bins.tsv must contain unique numeric bin values.")
    working = working.sort_values("bin", kind="stable")
    numeric = working[["mean_ld_score", "mean_chi_square", "mean_fitted_chi_square"]].to_numpy(float)
    if not np.all(np.isfinite(numeric)):
        raise LDSCInputError("Every h2 regression bin requires finite observed and fitted means.")
    weights = working["mean_regression_weight"].to_numpy(float)
    if np.any(~np.isfinite(weights)) or np.any(weights <= 0):
        raise LDSCInputError("Every h2 regression bin requires a finite positive mean_regression_weight.")
    relative_weights = weights / np.max(weights)
    x = working["mean_ld_score"].to_numpy(float)
    observed = working["mean_chi_square"].to_numpy(float)
    fitted = working["mean_fitted_chi_square"].to_numpy(float)

    figure, axes = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    points = axes.scatter(
        x,
        observed,
        c=relative_weights,
        cmap="Blues",
        vmin=0,
        vmax=1,
        s=46,
        edgecolor="white",
        linewidth=0.5,
        zorder=3,
    )
    axes.plot(x, fitted, color="#C44E52", linewidth=1.8, label="Fitted LDSC expectation", zorder=2)
    axes.set_xlabel("Mean LD Score")
    axes.set_ylabel(r"Mean chi-square statistic ($\chi^2$)")
    axes.set_title(f"LD Score regression diagnostic\n{trait_name}")
    axes.legend(frameon=False, loc="upper left")
    axes.grid(color="#E5E5E5", linewidth=0.8)
    axes.spines[["top", "right"]].set_visible(False)
    weight_scale = figure.colorbar(points, ax=axes, fraction=0.05, pad=0.03)
    weight_scale.set_label("Relative regression weight")
    return figure, axes


def _trait_order(metadata: dict[str, Any], table: pd.DataFrame) -> list[str]:
    _require_columns(table, {"trait_1", "trait_2", "rg", "rg_se"}, label="rg.tsv")
    names = metadata.get("trait_names")
    if not isinstance(names, list) or not names or any(not isinstance(value, str) for value in names):
        raise LDSCInputError("rg plotting metadata must declare a nonempty trait_names list.")
    if len(names) != len(set(names)):
        raise LDSCInputError("rg plotting metadata trait_names must be unique.")
    pair_names = set(table["trait_1"].astype(str)) | set(table["trait_2"].astype(str))
    missing = sorted(pair_names.difference(names))
    if missing:
        raise LDSCInputError("rg metadata trait_names omits table traits: " + ", ".join(missing) + ".")
    return names


def _plot_rg_heatmap(table: pd.DataFrame, trait_names: list[str]):
    n_traits = len(trait_names)
    index = {name: position for position, name in enumerate(trait_names)}
    estimates = np.full((n_traits, n_traits), np.nan)
    standard_errors = np.full((n_traits, n_traits), np.nan)
    seen: set[tuple[int, int]] = set()
    for _, row in table.iterrows():
        first = index[str(row["trait_1"])]
        second = index[str(row["trait_2"])]
        if first == second:
            raise LDSCInputError("rg.tsv must not contain self-pairs for the lower-triangular plot.")
        location = (max(first, second), min(first, second))
        if location in seen:
            raise LDSCInputError(
                f"rg.tsv contains pair {row['trait_1']!r} and {row['trait_2']!r} more than once."
            )
        seen.add(location)
        estimate = _optional_float(row["rg"])
        standard_error = _optional_float(row["rg_se"])
        if math.isfinite(estimate) and (not math.isfinite(standard_error) or standard_error < 0):
            raise LDSCInputError("Every available genetic-correlation estimate requires a non-negative rg_se.")
        estimates[location] = estimate
        standard_errors[location] = standard_error

    size = max(5.2, 0.92 * n_traits + 1.8)
    figure, axes = plt.subplots(figsize=(size, size), constrained_layout=True)
    mask = np.triu(np.ones_like(estimates, dtype=bool), k=0) | ~np.isfinite(estimates)
    axes.pcolormesh(
        np.ma.masked_where(mask, estimates),
        cmap=_RG_CMAP,
        vmin=-1,
        vmax=1,
        edgecolors="white",
        linewidth=0.8,
        shading="flat",
    )
    axes.set_xlim(0, n_traits)
    axes.set_ylim(n_traits, 0)
    axes.set_aspect("equal")
    for row_index in range(1, n_traits):
        for column_index in range(row_index):
            estimate = estimates[row_index, column_index]
            standard_error = standard_errors[row_index, column_index]
            if math.isfinite(estimate):
                label = _format_estimate_se(estimate, standard_error, multiline=True)
                color = "white" if abs(estimate) >= 0.55 else "#202020"
            else:
                label = "failed"
                color = "#555555"
                axes.add_patch(
                    plt.Rectangle(
                        (column_index, row_index),
                        1,
                        1,
                        facecolor="#E0E0E0",
                        edgecolor="white",
                        linewidth=0.8,
                    )
                )
            axes.text(
                column_index + 0.5,
                row_index + 0.5,
                label,
                ha="center",
                va="center",
                fontsize=max(6.5, 10.5 - 0.15 * n_traits),
                color=color,
            )
    axes.set_xticks(np.arange(n_traits) + 0.5, labels=trait_names)
    axes.set_yticks(np.arange(n_traits) + 0.5, labels=trait_names)
    axes.xaxis.tick_top()
    axes.tick_params(axis="x", labelrotation=45, length=0)
    axes.tick_params(axis="y", length=0)
    axes.set_title(r"Genetic correlation ($r_g$)" "\nCells show estimate (jackknife SE)", pad=14)
    axes.spines[:].set_visible(False)
    return figure, axes


def _plot_rg_anchor_forest(table: pd.DataFrame, trait_names: list[str]):
    anchors = table["trait_1"].astype(str).unique().tolist()
    if len(anchors) != 1:
        raise LDSCInputError("An anchor rg result must contain one common trait_1 across every pair.")
    anchor = anchors[0]
    if anchor not in trait_names:
        raise LDSCInputError(f"Anchor trait {anchor!r} is absent from metadata trait_names.")
    partners = [name for name in trait_names if name != anchor]
    indexed = table.assign(trait_2=table["trait_2"].astype(str)).set_index("trait_2", verify_integrity=True)
    unexpected = sorted(set(indexed.index).difference(partners))
    if unexpected:
        raise LDSCInputError("Anchor rg table contains unexpected partners: " + ", ".join(unexpected) + ".")
    figure, axes = plt.subplots(
        figsize=(7.4, max(3.2, 0.48 * len(partners) + 1.7)), constrained_layout=True
    )
    finite_intervals: list[tuple[float, float]] = []
    labels: list[tuple[int, float, str, bool]] = []
    for position, partner in enumerate(partners):
        if partner not in indexed.index:
            labels.append((position, 0.0, "failed", False))
            continue
        row = indexed.loc[partner]
        estimate = _optional_float(row["rg"])
        standard_error = _optional_float(row["rg_se"])
        if not math.isfinite(estimate):
            labels.append((position, 0.0, "failed", False))
            continue
        if not math.isfinite(standard_error) or standard_error < 0:
            raise LDSCInputError("Every available anchor genetic correlation requires a non-negative rg_se.")
        axes.errorbar(
            estimate,
            position,
            xerr=standard_error,
            fmt="o",
            color="#2F5D8A",
            ecolor="#4C78A8",
            elinewidth=1.6,
            capsize=3,
            markersize=5,
        )
        lower, upper = estimate - standard_error, estimate + standard_error
        finite_intervals.append((lower, upper))
        labels.append((position, upper, _format_estimate_se(estimate, standard_error), True))
    minimum = min(0.0, *(value[0] for value in finite_intervals)) if finite_intervals else -1.0
    maximum = max(0.0, *(value[1] for value in finite_intervals)) if finite_intervals else 1.0
    span = max(maximum - minimum, 0.2)
    left, right = minimum - 0.1 * span, maximum + 0.38 * span
    axes.set_xlim(left, right)
    for position, reference, label, available in labels:
        x = reference + 0.04 * span if available else left + 0.02 * (right - left)
        axes.text(x, position, label, va="center", fontsize=9, color="#202020" if available else "#666666")
    axes.axvline(0, color="#777777", linewidth=1, linestyle="--", zorder=0)
    axes.set_yticks(np.arange(len(partners)), labels=partners)
    axes.invert_yaxis()
    axes.set_xlabel(r"Genetic correlation ($r_g$; estimate ± 1 SE)")
    axes.set_title(f"Genetic correlation with {anchor}")
    _finish_horizontal_axes(axes)
    return figure, axes


def _plot_horizontal_enrichment_bars(
    labels: Sequence[str],
    estimates: Sequence[float],
    standard_errors: Sequence[float],
    *,
    title: str,
    y_label: str,
    bar_colors: str | Sequence[Any] = "#4D4D4D",
):
    label_values = list(labels)
    estimate_values = np.asarray(estimates, dtype=float)
    se_values = np.asarray(standard_errors, dtype=float)
    if len(label_values) != len(estimate_values) or len(estimate_values) != len(se_values):
        raise LDSCInputError("Enrichment labels, estimates, and standard errors must have equal lengths.")
    available = np.isfinite(estimate_values)
    invalid_se = available & ((~np.isfinite(se_values)) | (se_values < 0))
    if np.any(invalid_se):
        raise LDSCInputError("Every available enrichment estimate requires a finite non-negative enrichment_se.")
    if isinstance(bar_colors, str):
        colors: Any = bar_colors
    else:
        supplied = np.asarray(list(bar_colors))
        if len(supplied) != len(label_values):
            raise LDSCInputError("Enrichment bar colors must have one value per category.")
        colors = supplied[available]
    positions = np.arange(len(label_values))
    figure, axes = plt.subplots(
        figsize=(7.8, max(3.2, 0.48 * len(label_values) + 1.7)), constrained_layout=True
    )
    axes.axvline(1.0, color="#737373", linewidth=1.2, linestyle="--", zorder=0)
    axes.barh(
        positions[available],
        estimate_values[available],
        xerr=se_values[available],
        height=0.62,
        color=colors,
        edgecolor="none",
        error_kw={"ecolor": "#202020", "elinewidth": 1.3, "capsize": 3},
        zorder=2,
    )
    lower_values = estimate_values[available] - se_values[available]
    upper_values = estimate_values[available] + se_values[available]
    minimum = min(0.0, 1.0, *lower_values.tolist()) if lower_values.size else 0.0
    maximum = max(0.0, 1.0, *upper_values.tolist()) if upper_values.size else 1.0
    span = max(maximum - minimum, 0.2)
    axes.set_xlim(minimum - 0.06 * span, maximum + 0.34 * span)
    for position, estimate, standard_error in zip(positions, estimate_values, se_values):
        if math.isfinite(estimate):
            axes.text(
                estimate + standard_error + 0.03 * span,
                position,
                _format_estimate_se(estimate, standard_error),
                va="center",
                fontsize=8.5,
            )
        else:
            axes.text(minimum + 0.02 * span, position, "not available", va="center", color="#666666", fontsize=8.5)
    axes.set_yticks(positions, labels=label_values)
    axes.invert_yaxis()
    axes.set_xlabel("Heritability enrichment (estimate ± 1 SE)")
    axes.set_ylabel(y_label)
    axes.set_title(title)
    _finish_horizontal_axes(axes)
    return figure, axes


def _plot_query_pvalues(table: pd.DataFrame):
    _require_columns(table, {"category", "coefficient_p"}, label="partitioned_h2.tsv")
    columns = ["category", "coefficient_p"]
    if "coefficient_z" in table.columns:
        columns.append("coefficient_z")
    working = table.loc[:, columns].copy()
    working["coefficient_p"] = pd.to_numeric(working["coefficient_p"], errors="coerce")
    finite = np.isfinite(working["coefficient_p"])
    invalid = finite & ~working["coefficient_p"].between(0, 1, inclusive="both")
    if invalid.any():
        raise LDSCInputError("Every available coefficient_p must lie in [0, 1].")
    scores = pd.Series(np.nan, index=working.index, dtype=float)
    positive = finite & (working["coefficient_p"] > 0)
    scores.loc[positive] = -np.log10(working.loc[positive, "coefficient_p"])
    underflowed = finite & (working["coefficient_p"] == 0)
    if underflowed.any():
        if "coefficient_z" not in working:
            raise LDSCInputError("coefficient_p=0 requires coefficient_z to recover the one-sided tail probability.")
        z_scores = pd.to_numeric(working["coefficient_z"], errors="coerce")
        if (underflowed & ~np.isfinite(z_scores)).any():
            raise LDSCInputError("Queries with coefficient_p=0 require a finite coefficient_z.")
        scores.loc[underflowed] = -log_ndtr(-z_scores.loc[underflowed]) / np.log(10)
    working["neg_log10_p"] = scores
    working = working.sort_values("neg_log10_p", ascending=False, kind="stable", na_position="last")
    labels = working["category"].astype(str).tolist()
    score_values = working["neg_log10_p"].to_numpy(float)
    positions = np.arange(len(labels))
    figure, axes = plt.subplots(
        figsize=(8.0, max(4.0, 0.48 * len(labels) + 1.8)), constrained_layout=True
    )
    available = np.isfinite(score_values)
    axes.scatter(score_values[available], positions[available], color="#2F5D8A", s=42, zorder=3)
    for position in positions[~available]:
        axes.text(0, position, "not available", va="center", color="#666666", fontsize=8.5)
    axes.set_yticks(positions, labels=labels)
    axes.invert_yaxis()
    axes.set_xlim(left=0)
    axes.set_xlabel(r"Nominal conditional evidence ($-\log_{10}(P)$; one-sided test of $\tau > 0$)")
    axes.set_ylabel("Query annotation")
    axes.set_title(
        "Cell-type/query annotation evidence\n"
        "Each query is evaluated in a separate fit conditional on the baseline annotations"
    )
    _finish_horizontal_axes(axes)
    return figure, axes


def _plot_quantile_profile(table: pd.DataFrame, *, target_annotation: str):
    required = {"quantile", "target_value_lower", "target_value_upper", "enrichment", "enrichment_se"}
    _require_columns(table, required, label="quantile_h2.tsv")
    if table.empty:
        raise LDSCInputError("quantile_h2.tsv contains no quantile rows.")
    working = table.copy()
    for column in required:
        working[column] = pd.to_numeric(working[column], errors="coerce")
    if working["quantile"].isna().any() or working["quantile"].duplicated().any():
        raise LDSCInputError("quantile_h2.tsv must contain unique numeric quantile labels.")
    working = working.sort_values("quantile", kind="stable")
    lower = working["target_value_lower"].to_numpy(float)
    upper = working["target_value_upper"].to_numpy(float)
    if np.any(~np.isfinite(lower)) or np.any(~np.isfinite(upper)) or np.any(lower > upper):
        raise LDSCInputError("Every quantile requires finite ordered target-value bounds.")
    tick_labels = [
        f"Q{_format_quantile(quantile)}  [{low:g}, {high:g}]"
        for quantile, low, high in zip(working["quantile"], lower, upper)
    ]
    color_map = plt.get_cmap("YlOrRd")
    colors = color_map(np.linspace(0.18, 0.82, len(working)))
    return _plot_horizontal_enrichment_bars(
        tick_labels,
        working["enrichment"].to_numpy(float),
        working["enrichment_se"].to_numpy(float),
        title=f"Heritability enrichment across {target_annotation} quantiles",
        y_label="Target annotation quantile (low to high)",
        bar_colors=colors,
    )


def _finish_horizontal_axes(axes) -> None:
    axes.set_axisbelow(True)
    axes.grid(axis="x", color="#E5E5E5", linewidth=0.8)
    axes.grid(axis="y", visible=False)
    axes.spines[["top", "right", "left"]].set_visible(False)
    axes.spines["bottom"].set_color("#CCCCCC")
    axes.spines["bottom"].set_linewidth(0.8)
    axes.tick_params(axis="x", color="#CCCCCC")
    axes.tick_params(axis="y", length=0)


def _require_columns(table: pd.DataFrame, required: set[str], *, label: str) -> None:
    missing = sorted(required.difference(table.columns))
    if missing:
        raise LDSCInputError(f"{label} is missing required plotting columns: {', '.join(missing)}.")


def _optional_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _format_estimate_se(estimate: float, standard_error: float, *, multiline: bool = False) -> str:
    separator = "\n" if multiline else " "
    value = f"{estimate:.2f}{separator}({standard_error:.2f})"
    return value.replace("-", "−")


def _format_quantile(value: object) -> str:
    numeric = float(value)
    return str(int(numeric)) if numeric.is_integer() else f"{numeric:g}"
