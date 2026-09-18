"""Utilities for exporting value-specific SHAP summaries."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict

import pandas as pd


ReuseStats = DefaultDict[tuple[int, str, str], list[float]]


def new_reuse_stats() -> ReuseStats:
    return defaultdict(lambda: [0.0, 0.0])


def merge_reuse_stats(target: ReuseStats, source: ReuseStats) -> None:
    for key, (value_sum, count) in source.items():
        target[key][0] += value_sum
        target[key][1] += count


def _add_stat(stats: ReuseStats, key: tuple[int, str, str], shap_values: pd.Series) -> None:
    values = pd.to_numeric(shap_values, errors="coerce").dropna()
    if values.empty:
        return
    stats[key][0] += float(values.sum())
    stats[key][1] += float(len(values))


def collect_site_state_reuse_stats(
    stats: ReuseStats,
    feature_values: pd.DataFrame,
    shap_values: pd.DataFrame,
    feature_cols: list[str],
    max_site: int = 550,
) -> None:
    for feature in feature_cols:
        if feature not in feature_values.columns or feature not in shap_values.columns:
            continue

        change_match = re.fullmatch(r"change_(\d+)", feature)
        if change_match:
            site = int(change_match.group(1))
            if site > max_site:
                continue
            values = pd.to_numeric(feature_values[feature], errors="coerce")
            _add_stat(stats, (site, "change", "yes"), shap_values.loc[values == 1, feature])
            _add_stat(stats, (site, "change", "no"), shap_values.loc[values == 0, feature])
            continue

        virus_match = re.fullmatch(r"([A-Z])(\d+)", feature)
        if virus_match:
            aa, site_text = virus_match.groups()
            site = int(site_text)
            if site > max_site:
                continue
            values = pd.to_numeric(feature_values[feature], errors="coerce")
            _add_stat(stats, (site, "virus aa", aa), shap_values.loc[values == 1, feature])
            continue

        serum_match = re.fullmatch(r"(\d+)([A-Z])", feature)
        if serum_match:
            site_text, aa = serum_match.groups()
            site = int(site_text)
            if site > max_site:
                continue
            values = pd.to_numeric(feature_values[feature], errors="coerce")
            _add_stat(stats, (site, "sera aa", aa), shap_values.loc[values == 1, feature])


def collect_substitution_reuse_stats(
    stats: ReuseStats,
    feature_values: pd.DataFrame,
    shap_values: pd.DataFrame,
    feature_cols: list[str],
    max_site: int = 550,
) -> None:
    for feature in feature_cols:
        match = re.fullmatch(r"sub_(\d+)_([A-Z])_([A-Z])", feature)
        if match is None or feature not in feature_values.columns or feature not in shap_values.columns:
            continue
        site_text, from_aa, to_aa = match.groups()
        site = int(site_text)
        if site > max_site:
            continue

        substitution = f"{from_aa}{site}{to_aa}"
        values = pd.to_numeric(feature_values[feature], errors="coerce")
        _add_stat(stats, (site, "substitution present", substitution), shap_values.loc[values == 1, feature])
        _add_stat(stats, (site, "substitution absent", substitution), shap_values.loc[values == 0, feature])


def reuse_stats_to_table(stats: ReuseStats) -> pd.DataFrame:
    rows = []
    for (site, covariate_name, covariate_value), (value_sum, count) in stats.items():
        if count == 0:
            continue
        model_type = "substitution" if covariate_name.startswith("substitution") else "site-state"
        rows.append(
            {
                "site": site,
                "model_type": model_type,
                "covariate_name": covariate_name,
                "covariate_value": covariate_value,
                "shap_value": value_sum / count,
                "n_observations": int(count),
            }
        )

    columns = ["site", "model_type", "covariate_name", "covariate_value", "shap_value", "n_observations"]
    return pd.DataFrame(rows, columns=columns).sort_values(columns[:4]).reset_index(drop=True)


def format_reuse_feature_label(row: pd.Series) -> str:
    site = int(row["site"])
    model_type = str(row["model_type"])
    covariate_name = str(row["covariate_name"])
    covariate_value = str(row["covariate_value"])

    if model_type == "substitution":
        return covariate_value

    if covariate_name == "change":
        change_label = "(+)change" if covariate_value.lower() in {"yes", "true", "1"} else "(-)change"
        return f"{site:>3} {change_label}"
    if covariate_name == "virus aa":
        return f"{site}{covariate_value}"
    if covariate_name == "sera aa":
        return f"{covariate_value}{site}"
    return f"{site:>3} {covariate_name}={covariate_value}"


def top_signed_reuse_rows(
    reusable_shap: pd.DataFrame,
    model_type: str,
    top_n: int = 20,
) -> pd.DataFrame:
    subset = reusable_shap[reusable_shap["model_type"].eq(model_type)].copy()
    if subset.empty:
        return subset

    subset["shap_value"] = pd.to_numeric(subset["shap_value"], errors="coerce")
    subset = subset.dropna(subset=["shap_value"])
    subset["abs_shap_value"] = subset["shap_value"].abs()
    subset = subset.sort_values(
        ["abs_shap_value", "site", "covariate_name", "covariate_value"],
        ascending=[False, True, True, True],
    ).head(top_n)
    subset["label"] = subset.apply(format_reuse_feature_label, axis=1)
    return subset


def save_signed_reusable_shap_figure(
    reusable_shap: pd.DataFrame,
    out_path: str | Path,
    top_n: int = 20,
) -> None:
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    panel_specs = [
        ("site-state", "A", "Site-state model"),
        ("substitution", "B", "Substitution model"),
    ]
    panel_rows = {
        model_type: top_signed_reuse_rows(reusable_shap, model_type=model_type, top_n=top_n)
        for model_type, _, _ in panel_specs
    }

    positive_color = "#2b8cbe"
    negative_color = "#de8f05"
    fig, axes = plt.subplots(1, 2, figsize=(16.0, 8.6))

    for ax, (model_type, panel_label, title) in zip(axes, panel_specs):
        plot_df = panel_rows[model_type].iloc[::-1].copy()
        if plot_df.empty:
            ax.text(0.5, 0.5, "No SHAP values available", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue

        values = plot_df["shap_value"].to_numpy(dtype=float)
        labels = plot_df["label"].tolist()
        y = np.arange(len(plot_df))
        colors = [positive_color if value >= 0 else negative_color for value in values]
        ax.barh(y, values, color=colors, edgecolor="white", linewidth=0.6)
        ax.axvline(0, color="#333333", linewidth=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.set_xlabel("mean signed SHAP value")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.text(
            -0.16,
            1.04,
            panel_label,
            transform=ax.transAxes,
            fontsize=14,
            fontweight="bold",
            ha="left",
            va="bottom",
        )
        axis_limit = 1.12 * (float(np.max(np.abs(values))) or 1.0)
        ax.set_xlim(-axis_limit, axis_limit)
        ax.grid(axis="x", color="#d9d9d9", linewidth=0.6)
        ax.grid(axis="y", visible=False)
        ax.set_facecolor("white")
        sns.despine(ax=ax)

    fig.tight_layout(w_pad=4.0)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
