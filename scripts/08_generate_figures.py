#!/usr/bin/env python3
"""Generate publication-oriented H3N2 comparison figures."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
import seaborn as sns

from common import append_qc_log, get_prediction_color_feature, humanize_distance_label, load_config
from paper_sites import NEHER2016_H3_SITE_ROWS, SHAH2024_H3_SITE_ROWS, WIC2023_H3_SITE_ROWS
from shap_reuse import save_signed_reusable_shap_figure


sns.set_theme(style="ticks", context="paper")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.size"] = 12
matplotlib.rcParams["axes.labelsize"] = 12
matplotlib.rcParams["xtick.labelsize"] = 12
matplotlib.rcParams["ytick.labelsize"] = 12
matplotlib.rcParams["legend.fontsize"] = 12
matplotlib.rcParams["legend.title_fontsize"] = 12
matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False
matplotlib.rcParams["axes.linewidth"] = 0.8
matplotlib.rcParams["axes.edgecolor"] = "#333333"
matplotlib.rcParams["grid.color"] = "#d9d9d9"
matplotlib.rcParams["grid.linewidth"] = 0.6
matplotlib.rcParams["grid.alpha"] = 0.8
MONO_FONT = "DejaVu Sans Mono"
KOEL_SITES = {145, 155, 156, 158, 159, 189, 193}
NEHER_SITES = {int(row["site"]) for row in NEHER2016_H3_SITE_ROWS}
HARVEY_SITES = {int(row["site"]) for row in WIC2023_H3_SITE_ROWS}
SHAH_SITES = {int(row["site"]) for row in SHAH2024_H3_SITE_ROWS}
TOP30_TEXT_SIZE = 16
STACKED_NOTE_SIZE = 14
DIAGNOSTIC_POINT_SIZE = 14
NEHER_ANALOG_POINT_SIZE = 18
RESIDUAL_POINT_SIZE = 24


def finish_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    ax.set_facecolor("white")
    if grid_axis is None:
        ax.grid(False)
    else:
        ax.grid(axis=grid_axis, color="#d9d9d9", linewidth=0.6)
        other_axis = "x" if grid_axis == "y" else "y"
        ax.grid(axis=other_axis, visible=False)
    sns.despine(ax=ax)


def maybe_set_component_ticks(ax: plt.Axes, step: float) -> None:
    xmax = ax.get_xlim()[1]
    if xmax > 0:
        ax.xaxis.set_major_locator(MultipleLocator(step))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    return parser.parse_args()


def build_reciprocal_pairs(df: pd.DataFrame) -> pd.DataFrame:
    left = df[
        [
            "virus_strain_matched",
            "serum_strain_matched",
            "standardized_titer",
            "corrected_titer",
        ]
    ].copy()
    left = left.rename(
        columns={
            "virus_strain_matched": "virus_a",
            "serum_strain_matched": "virus_b",
            "standardized_titer": "standardized_ab",
            "corrected_titer": "corrected_ab",
        }
    )
    right = left.rename(
        columns={
            "virus_a": "virus_b",
            "virus_b": "virus_a",
            "standardized_ab": "standardized_ba",
            "corrected_ab": "corrected_ba",
        }
    )
    merged = left.merge(right, on=["virus_a", "virus_b"], how="inner")
    return merged[merged["virus_a"] < merged["virus_b"]].copy()


def humanize_feature_label(feature: str) -> str:
    if feature.startswith("change_"):
        site = feature.split("_")[1]
        return f"{site.rjust(3)} "
    if feature.startswith("sub_"):
        _, site, serum_aa, virus_aa = feature.split("_", 3)
        return f"{site.rjust(3)} {serum_aa}->{virus_aa}"
    if re.fullmatch(r"[A-Z]\d+", feature):
        aa = feature[0]
        site = feature[1:]
        return f"{aa}{site.rjust(3)} "
    if re.fullmatch(r"\d+[A-Z]", feature):
        site = feature[:-1]
        aa = feature[-1]
        return f"{site.rjust(3)}{aa}"
    return {
        "temporal_distance": "temporal distance ",
        "patristic_distance": "patristic distance ",
        "lab": "lab ",
    }.get(feature, feature.replace("_", " ").lower() + " ")


def humanize_model_label(model_type: str) -> str:
    mapping = {
        "additive_only": "potency + avidity only",
        "temporal_only": "temporal distance only",
        "randomized_binary": "randomized site changes",
        "binary_site": "site changed + time",
        "site_state": "site aa state + time",
        "substitution_identity": "substitution identity + time",
    }
    return mapping.get(model_type, model_type.replace("_", " "))


def humanize_split_label(split_type: str) -> str:
    mapping = {
        "random_10pct_measurements": "random 10% measurements",
        "grouped_10pct_viruses": "all measurements for 10% viruses",
    }
    return mapping.get(split_type, split_type.replace("_", " "))


def humanize_homologous_source(label: str) -> str:
    mapping = {
        "homologous": "measured homologous",
        "max_observed": "max observed fallback",
    }
    return mapping.get(label, label.replace("_", " "))


def site_flag_suffix(row: pd.Series) -> str:
    flags = []
    if bool(row.get("is_koel_site", False)):
        flags.append("K")
    if bool(row.get("in_neher2016", row.get("paper_named", False))):
        flags.append("NB")
    if bool(row.get("in_harvey2023", False)):
        flags.append("HW")
    if bool(row.get("in_shah2024", row.get("named_in_shah2024", False))):
        flags.append("SW")
    return f" [{', '.join(flags)}]" if flags else ""


def add_reference_flags(df: pd.DataFrame) -> pd.DataFrame:
    merged = df.copy()
    merged["in_koel"] = merged["site"].astype(int).isin(KOEL_SITES)
    merged["in_neher2016"] = merged["site"].astype(int).isin(NEHER_SITES)
    merged["in_harvey2023"] = merged["site"].astype(int).isin(HARVEY_SITES)
    merged["in_shah2024"] = merged["site"].astype(int).isin(SHAH_SITES)
    return merged


def read_optional_table(path: Path, sep: str = ",") -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep=sep)
    except pd.errors.EmptyDataError:
        return None


def save_signed_feature_summary(feature_summary: pd.DataFrame, out_path: Path) -> None:
    top_features = feature_summary.head(30).iloc[::-1].copy()
    top_features["label"] = top_features["feature"].map(humanize_feature_label)
    vmax = float(np.max(np.abs(top_features["mean_signed_shap"]))) or 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm")
    colors = [cmap(norm(v)) for v in top_features["mean_signed_shap"]]

    fig, ax = plt.subplots(figsize=(10.5, 9))
    ax.barh(top_features["label"], top_features["mean_abs_shap"], color=colors)
    ax.set_xlabel("Mean |SHAP|")
    ax.set_ylabel("")
    ax.xaxis.label.set_fontfamily(MONO_FONT)
    for tick in ax.get_yticklabels() + ax.get_xticklabels():
        tick.set_fontfamily(MONO_FONT)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label("Mean signed SHAP")
    cb.ax.yaxis.label.set_fontfamily(MONO_FONT)
    for tick in cb.ax.get_yticklabels():
        tick.set_fontfamily(MONO_FONT)
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_site_importance_bar(site_importance: pd.DataFrame, out_path: Path, colors: np.ndarray, x_label: str) -> None:
    fig, ax = plt.subplots(figsize=(15, 4.8))
    ax.bar(site_importance["site"], site_importance["mean_abs_shap"], color=colors, width=1.0)
    ax.set_xlabel(x_label)
    ax.set_ylabel("mean |SHAP|")
    finish_axes(ax=ax, grid_axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_site_component_plot(site_summary: pd.DataFrame, out_path: Path) -> None:
    top_sites = site_summary.nsmallest(30, "rank").sort_values("rank", ascending=False).copy()
    component_specs = [
        ("mean_abs_shap_change", "mean_signed_shap_change"),
        ("mean_abs_shap_virus_aa", "mean_signed_shap_virus_aa"),
        ("mean_abs_shap_serum_aa", "mean_signed_shap_serum_aa"),
    ]
    signed_values = []
    for _, signed_col in component_specs:
        signed_values.extend(top_sites[signed_col].tolist())
    vmax = float(np.max(np.abs(signed_values))) or 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm")

    fig, ax = plt.subplots(figsize=(14.0, 10.8))
    y = np.arange(len(top_sites))
    left = np.zeros(len(top_sites), dtype=float)
    for abs_col, signed_col in component_specs:
        widths = top_sites[abs_col].to_numpy(dtype=float)
        colors = [cmap(norm(v)) for v in top_sites[signed_col].to_numpy(dtype=float)]
        ax.barh(y, widths, left=left, color=colors, edgecolor="white", linewidth=0.6)
        left += widths

    labels = [f"{int(row['site'])}{site_flag_suffix(row)}" for _, row in top_sites.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("mean |SHAP|")
    ax.set_ylabel("site")
    ax.xaxis.label.set_fontsize(TOP30_TEXT_SIZE)
    ax.yaxis.label.set_fontsize(TOP30_TEXT_SIZE)
    ax.tick_params(axis="both", labelsize=TOP30_TEXT_SIZE)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label("mean signed SHAP")
    cb.ax.yaxis.label.set_size(TOP30_TEXT_SIZE)
    cb.ax.tick_params(labelsize=TOP30_TEXT_SIZE)
    ax.text(
        0.985,
        0.035,
        "segment width = mean |SHAP|",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=STACKED_NOTE_SIZE,
        bbox={"facecolor": "white", "edgecolor": "#444444", "alpha": 1.0, "boxstyle": "round,pad=0.35"},
    )
    maybe_set_component_ticks(ax, 0.05)
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_substitution_site_plot(site_summary: pd.DataFrame, out_path: Path) -> None:
    top_sites = site_summary.nsmallest(30, "rank").sort_values("rank", ascending=False).copy()
    vmax = float(np.max(np.abs(top_sites["mean_signed_shap"]))) or 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm")
    colors = [cmap(norm(v)) for v in top_sites["mean_signed_shap"]]

    labels = []
    for _, row in top_sites.iterrows():
        label = f"{int(row['site'])}{site_flag_suffix(row)}"
        top_terms = row["top_substitution_terms"] if isinstance(row["top_substitution_terms"], str) else ""
        if top_terms:
            label = f"{label} | {top_terms}"
        labels.append(label)

    fig, ax = plt.subplots(figsize=(15.2, 10.2))
    ax.barh(labels, top_sites["mean_abs_shap"], color=colors)
    ax.set_xlabel("mean |SHAP|")
    ax.set_ylabel("site")
    ax.xaxis.label.set_fontsize(TOP30_TEXT_SIZE)
    ax.yaxis.label.set_fontsize(TOP30_TEXT_SIZE)
    ax.tick_params(axis="both", labelsize=TOP30_TEXT_SIZE)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label("mean signed SHAP")
    cb.ax.yaxis.label.set_size(TOP30_TEXT_SIZE)
    cb.ax.tick_params(labelsize=TOP30_TEXT_SIZE)
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_site_stability(site_stability: pd.DataFrame, out_path: Path) -> None:
    top_sites = site_stability.nsmallest(30, "rank").sort_values("mean_rank", ascending=True).copy()
    freq_norm = Normalize(0.0, 1.0)
    freq_cmap = plt.get_cmap("viridis")
    colors = [freq_cmap(freq_norm(v)) for v in top_sites["top30_frequency"].to_numpy(dtype=float)]

    fig, ax = plt.subplots(figsize=(11.3, 10.2))
    y = np.arange(len(top_sites))
    ax.hlines(
        y=y,
        xmin=top_sites["mean_rank"] - top_sites["sd_rank"],
        xmax=top_sites["mean_rank"] + top_sites["sd_rank"],
        color="#9e9e9e",
        linewidth=2,
    )
    ax.scatter(
        top_sites["mean_rank"],
        y,
        c=colors,
        s=55,
        edgecolors="black",
        linewidths=0.4,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([f"{int(row['site'])}{site_flag_suffix(row)}" for _, row in top_sites.iterrows()])
    ax.set_xlabel("mean fold rank ± SD")
    ax.set_ylabel("site")
    ax.invert_yaxis()
    sm = plt.cm.ScalarMappable(cmap=freq_cmap, norm=freq_norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label("top 30 frequency")
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_prediction_scatter(df: pd.DataFrame, pred_col: str, out_path: Path, color_col: str) -> None:
    fig, ax = plt.subplots(figsize=(8.3, 6.8))
    scatter = ax.scatter(
        df["standardized_titer"],
        df[pred_col],
        c=df[color_col],
        cmap="viridis",
        s=DIAGNOSTIC_POINT_SIZE,
        alpha=0.7,
    )
    lo = float(np.nanmin([df["standardized_titer"].min(), df[pred_col].min()]))
    hi = float(np.nanmax([df["standardized_titer"].max(), df[pred_col].max()]))
    ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("observed standardized titer")
    ax.set_ylabel("predicted standardized titer")
    cb = fig.colorbar(scatter, ax=ax, pad=0.02)
    cb.set_label(humanize_distance_label(color_col))
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_cv_performance(cv_results: pd.DataFrame, out_path: Path) -> None:
    cv_plot = cv_results[
        cv_results["model_type"].isin(
            ["additive_only", "temporal_only", "randomized_binary", "binary_site", "site_state", "substitution_identity"]
        )
    ].copy()
    order = [
        "potency + avidity only",
        "temporal distance only",
        "randomized site changes",
        "site changed + time",
        "site aa state + time",
        "substitution identity + time",
    ]
    fig, ax = plt.subplots(figsize=(4.6, 7.6))
    sns.boxplot(
        data=cv_plot,
        x="model_label",
        y="rmse",
        order=order,
        color="#9ecae1",
        width=0.55,
        showfliers=False,
        linewidth=0.9,
        ax=ax,
    )
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", va="top")
    ax.set_ylabel("RMSE")
    ax.set_xlabel("")
    finish_axes(ax=ax, grid_axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_residualization_diagnostic(titers: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 6.8))
    sns.scatterplot(
        data=titers,
        x="standardized_titer",
        y="corrected_titer",
        hue="homologous_source_label",
        alpha=0.6,
        s=RESIDUAL_POINT_SIZE,
        linewidth=0,
        ax=ax,
    )
    ax.set_xlabel("standardized titer")
    ax.set_ylabel("corrected titer")
    legend = ax.legend(title="homologous titer source", frameon=True, framealpha=1.0)
    if legend is not None:
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("#444444")
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_neher_analog(neher_pred: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.8), sharex=True, sharey=True)
    panel_map = [
        "random_10pct_measurements",
        "grouped_10pct_viruses",
    ]
    for ax, split_name in zip(axes, panel_map):
        sub = neher_pred[neher_pred["split_type"] == split_name]
        ax.scatter(
            sub["standardized_titer"],
            sub["predicted_standardized_titer_site_state"],
            s=NEHER_ANALOG_POINT_SIZE,
            alpha=0.6,
            color="#2b8cbe",
            edgecolors="none",
        )
        lo = min(sub["standardized_titer"].min(), sub["predicted_standardized_titer_site_state"].min())
        hi = max(sub["standardized_titer"].max(), sub["predicted_standardized_titer_site_state"].max())
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=1)
        ax.set_xlabel("observed standardized titer")
        finish_axes(ax=ax, grid_axis=None)
    axes[0].set_ylabel("predicted standardized titer")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_reciprocal_symmetry(titers: pd.DataFrame, out_path: Path) -> None:
    reciprocal = build_reciprocal_pairs(titers)
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    ax.scatter(
        reciprocal["standardized_ab"],
        reciprocal["standardized_ba"],
        s=DIAGNOSTIC_POINT_SIZE,
        alpha=0.5,
        color="#756bb1",
    )
    lo = min(reciprocal["standardized_ab"].min(), reciprocal["standardized_ba"].min())
    hi = max(reciprocal["standardized_ab"].max(), reciprocal["standardized_ba"].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("standardized titer A vs B")
    ax.set_ylabel("standardized titer B vs A")
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def save_heldout_comparison(neher_metrics: pd.DataFrame, out_path: Path) -> None:
    heldout_plot = neher_metrics[
        neher_metrics["model_type"].isin(["additive_only", "binary_site", "site_state", "substitution_identity"])
    ].copy()
    hue_order = [
        "potency + avidity only",
        "site changed + time",
        "site aa state + time",
        "substitution identity + time",
    ]
    split_order = ["random 10% measurements", "all measurements for 10% viruses"]
    fig, ax = plt.subplots(figsize=(10.4, 5.4))
    sns.barplot(
        data=heldout_plot,
        x="split_label",
        y="rmse",
        hue="model_label",
        hue_order=hue_order,
        order=split_order,
        palette=["#bdbdbd", "#6baed6", "#08519c", "#ef6548"],
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("RMSE")
    legend = ax.legend(
        title="model",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        borderaxespad=0.0,
        frameon=True,
        framealpha=1.0,
    )
    if legend is not None:
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("#444444")
    finish_axes(ax=ax, grid_axis="y")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def load_analysis_tables(cfg: dict) -> dict[str, pd.DataFrame]:
    out_dir = Path(cfg["output_dir"])
    subtype = cfg["subtype"]
    site_summary = add_reference_flags(pd.read_csv(out_dir / f"{subtype}_site_summary.tsv", sep="\t"))
    sub_site_summary = add_reference_flags(pd.read_csv(out_dir / f"{subtype}_substitution_site_summary.tsv", sep="\t"))
    site_importance = read_optional_table(out_dir / f"{subtype}_site_importance.csv")
    sub_site_importance = read_optional_table(out_dir / f"{subtype}_substitution_site_importance.csv")
    tables = {
        "feature_summary": read_optional_table(out_dir / f"{subtype}_feature_shap_importance.csv"),
        "reusable_shap": read_optional_table(out_dir / f"{subtype}_reusable_shap_values.tsv", sep="\t"),
        "site_importance": add_reference_flags(site_importance) if site_importance is not None else site_summary.copy(),
        "site_summary": site_summary,
        "site_stability": add_reference_flags(pd.read_csv(out_dir / f"{subtype}_site_stability.tsv", sep="\t")),
        "sub_feature_summary": read_optional_table(out_dir / f"{subtype}_substitution_feature_shap_importance.csv"),
        "sub_site_importance": (
            add_reference_flags(sub_site_importance) if sub_site_importance is not None else sub_site_summary.copy()
        ),
        "sub_site_summary": sub_site_summary,
        "sub_site_stability": add_reference_flags(
            pd.read_csv(out_dir / f"{subtype}_substitution_site_stability.tsv", sep="\t")
        ),
        "cv_results": read_optional_table(out_dir / f"{subtype}_cv_results.csv"),
        "oof": read_optional_table(out_dir / f"{subtype}_oof_predictions.csv"),
        "titers": read_optional_table(out_dir / f"{subtype}_matched_titers.csv"),
        "neher_pred": read_optional_table(out_dir / f"{subtype}_neher_analog_predictions.csv"),
        "neher_metrics": read_optional_table(out_dir / f"{subtype}_neher_analog_metrics.csv"),
    }
    if tables["cv_results"] is not None:
        tables["cv_results"]["model_label"] = tables["cv_results"]["model_type"].map(humanize_model_label)
    if tables["neher_metrics"] is not None:
        tables["neher_metrics"]["model_label"] = tables["neher_metrics"]["model_type"].map(humanize_model_label)
        tables["neher_metrics"]["split_label"] = tables["neher_metrics"]["split_type"].map(humanize_split_label)
    if tables["titers"] is not None:
        tables["titers"]["homologous_source_label"] = tables["titers"]["homologous_titer_source"].map(
            humanize_homologous_source
        )
    return tables


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    figures_dir = Path(cfg["figures_dir"])
    subtype = cfg["subtype"]
    color_col = get_prediction_color_feature(cfg)
    tables = load_analysis_tables(cfg)
    lines = []

    if tables["feature_summary"] is not None:
        save_signed_feature_summary(
            feature_summary=tables["feature_summary"],
            out_path=figures_dir / f"{subtype}_shap_summary_top30.pdf",
        )
        lines.append(f"Saved SHAP summary figure: {figures_dir / f'{subtype}_shap_summary_top30.pdf'}")
    else:
        lines.append("Skipped SHAP summary figure: missing feature_shap_importance.csv")
    site_state_colors = np.where(tables["site_importance"]["is_koel_site"].to_numpy(), "#d7301f", "#9ecae1")
    save_site_importance_bar(
        site_importance=tables["site_importance"],
        out_path=figures_dir / f"{subtype}_site_importance_bar.pdf",
        colors=site_state_colors,
        x_label="mature HA position",
    )
    lines.append(f"Saved site importance figure: {figures_dir / f'{subtype}_site_importance_bar.pdf'}")
    save_site_component_plot(
        site_summary=tables["site_summary"],
        out_path=figures_dir / f"{subtype}_site_component_top30.pdf",
    )
    lines.append(f"Saved site component figure: {figures_dir / f'{subtype}_site_component_top30.pdf'}")
    save_site_stability(
        site_stability=tables["site_stability"],
        out_path=figures_dir / f"{subtype}_site_stability_top30.pdf",
    )
    lines.append(f"Saved site stability figure: {figures_dir / f'{subtype}_site_stability_top30.pdf'}")
    if tables["oof"] is not None:
        save_prediction_scatter(
            df=tables["oof"],
            pred_col="predicted_standardized_titer_site_state",
            out_path=figures_dir / f"{subtype}_predicted_vs_actual.pdf",
            color_col=color_col,
        )
        lines.append(f"Saved predicted-vs-actual figure: {figures_dir / f'{subtype}_predicted_vs_actual.pdf'}")
    else:
        lines.append("Skipped predicted-vs-actual figure: missing oof_predictions.csv")

    if tables["sub_feature_summary"] is not None:
        save_signed_feature_summary(
            feature_summary=tables["sub_feature_summary"],
            out_path=figures_dir / f"{subtype}_substitution_shap_summary_top30.pdf",
        )
        lines.append(
            f"Saved substitution SHAP summary figure: {figures_dir / f'{subtype}_substitution_shap_summary_top30.pdf'}"
        )
    else:
        lines.append("Skipped substitution SHAP summary figure: missing substitution_feature_shap_importance.csv")
    if tables["reusable_shap"] is not None:
        save_signed_reusable_shap_figure(
            reusable_shap=tables["reusable_shap"],
            out_path=figures_dir / f"{subtype}_signed_feature_shap_top20.pdf",
        )
        lines.append(
            f"Saved signed feature SHAP figure: {figures_dir / f'{subtype}_signed_feature_shap_top20.pdf'}"
        )
    else:
        lines.append("Skipped signed feature SHAP figure: missing reusable_shap_values.tsv")
    sub_colors = np.where(tables["sub_site_importance"]["is_koel_site"].to_numpy(), "#d7301f", "#9ecae1")
    save_site_importance_bar(
        site_importance=tables["sub_site_importance"],
        out_path=figures_dir / f"{subtype}_substitution_site_importance_bar.pdf",
        colors=sub_colors,
        x_label="mature HA position",
    )
    lines.append(
        f"Saved substitution site importance figure: {figures_dir / f'{subtype}_substitution_site_importance_bar.pdf'}"
    )
    save_substitution_site_plot(
        site_summary=tables["sub_site_summary"],
        out_path=figures_dir / f"{subtype}_substitution_site_component_top30.pdf",
    )
    lines.append(
        f"Saved substitution site summary figure: {figures_dir / f'{subtype}_substitution_site_component_top30.pdf'}"
    )
    save_site_stability(
        site_stability=tables["sub_site_stability"],
        out_path=figures_dir / f"{subtype}_substitution_site_stability_top30.pdf",
    )
    lines.append(
        f"Saved substitution site stability figure: {figures_dir / f'{subtype}_substitution_site_stability_top30.pdf'}"
    )
    if tables["oof"] is not None:
        save_prediction_scatter(
            df=tables["oof"],
            pred_col="predicted_standardized_titer_substitution",
            out_path=figures_dir / f"{subtype}_substitution_predicted_vs_actual.pdf",
            color_col=color_col,
        )
        lines.append(
            f"Saved substitution predicted-vs-actual figure: {figures_dir / f'{subtype}_substitution_predicted_vs_actual.pdf'}"
        )
    else:
        lines.append("Skipped substitution predicted-vs-actual figure: missing oof_predictions.csv")

    if tables["cv_results"] is not None:
        save_cv_performance(
            cv_results=tables["cv_results"],
            out_path=figures_dir / f"{subtype}_cv_performance_summary.pdf",
        )
        lines.append(f"Saved CV summary figure: {figures_dir / f'{subtype}_cv_performance_summary.pdf'}")
    else:
        lines.append("Skipped CV summary figure: missing cv_results.csv")
    if tables["titers"] is not None:
        save_residualization_diagnostic(
            titers=tables["titers"],
            out_path=figures_dir / f"{subtype}_residualization_diagnostic.pdf",
        )
        save_reciprocal_symmetry(
            titers=tables["titers"],
            out_path=figures_dir / f"{subtype}_reciprocal_symmetry.pdf",
        )
        lines.append(f"Saved standardization diagnostic figure: {figures_dir / f'{subtype}_residualization_diagnostic.pdf'}")
        lines.append(f"Saved reciprocal symmetry figure: {figures_dir / f'{subtype}_reciprocal_symmetry.pdf'}")
    else:
        lines.append("Skipped standardization diagnostic figure: missing matched_titers.csv")
        lines.append("Skipped reciprocal symmetry figure: missing matched_titers.csv")
    if tables["neher_pred"] is not None:
        save_neher_analog(
            neher_pred=tables["neher_pred"],
            out_path=figures_dir / f"{subtype}_neher2016_fig2_analog.pdf",
        )
        lines.append(f"Saved Neher 2016 analog figure: {figures_dir / f'{subtype}_neher2016_fig2_analog.pdf'}")
    else:
        lines.append("Skipped Neher 2016 analog figure: missing neher_analog_predictions.csv")
    if tables["neher_metrics"] is not None:
        save_heldout_comparison(
            neher_metrics=tables["neher_metrics"],
            out_path=figures_dir / f"{subtype}_heldout_model_comparison.pdf",
        )
        lines.append(f"Saved held-out comparison figure: {figures_dir / f'{subtype}_heldout_model_comparison.pdf'}")
    else:
        lines.append("Skipped held-out comparison figure: missing neher_analog_metrics.csv")
    append_qc_log(cfg, "08_generate_figures", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
