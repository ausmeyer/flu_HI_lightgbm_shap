#!/usr/bin/env python3
"""Generate WIC analysis figures analogous to the original H3N2 LightGBM outputs."""

from __future__ import annotations

import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd
import seaborn as sns

from common import append_qc_log, get_prediction_color_feature, humanize_distance_label, load_config


sns.set_theme(style="ticks", context="paper")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False
matplotlib.rcParams["axes.linewidth"] = 0.8
matplotlib.rcParams["axes.edgecolor"] = "#333333"
matplotlib.rcParams["grid.color"] = "#d9d9d9"
matplotlib.rcParams["grid.linewidth"] = 0.6
matplotlib.rcParams["grid.alpha"] = 0.8
MONO_FONT = "DejaVu Sans Mono"


def finish_axes(ax: plt.Axes, grid_axis: str | None = "y") -> None:
    ax.set_facecolor("white")
    if grid_axis is None:
        ax.grid(False)
    else:
        ax.grid(axis=grid_axis, color="#d9d9d9", linewidth=0.6)
        other_axis = "x" if grid_axis == "y" else "y"
        ax.grid(axis=other_axis, visible=False)
    sns.despine(ax=ax)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    return parser.parse_args()


def build_reciprocal_pairs(df: pd.DataFrame) -> pd.DataFrame:
    left = df[
        [
            "virusStrain",
            "serumStrain",
            "standardized_titer",
            "context_corrected_titer",
        ]
    ].copy()
    left = left.rename(
        columns={
            "virusStrain": "virus_a",
            "serumStrain": "virus_b",
            "standardized_titer": "standardized_ab",
            "context_corrected_titer": "corrected_ab",
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
    merged = merged[merged["virus_a"] < merged["virus_b"]].copy()
    return merged


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
    mapping = {
        "temporal_distance": "temporal distance ",
        "patristic_distance": "patristic distance ",
        "virus_passage_class": "virus passage ",
        "reference_passage_class": "serum passage ",
        "virus_passage_known": "virus passage known ",
        "reference_passage_known": "serum passage known ",
        "antiserum_class": "antiserum class ",
        "rbc_type": "rbc type ",
        "rbc_species": "rbc species ",
        "oseltamivir_present": "oseltamivir present ",
        "oseltamivir_nm": "oseltamivir nm ",
    }
    return mapping.get(feature, feature.replace("_", " ").lower() + " ")


def humanize_model_label(model_type: str) -> str:
    mapping = {
        "additive_only": "context additive only",
        "context_only": "context only",
        "randomized_binary_context": "randomized site changes + context",
        "binary_site_context": "site changed + context",
        "site_state_context": "site aa state + context",
        "substitution_identity_context": "substitution identity + context",
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
    neher = bool(row.get("named_in_neher2016", False))
    wic = bool(row.get("named_in_wic2023", False))
    if neher and wic:
        return " [Neher, Harvey]"
    if neher:
        return " [Neher]"
    if wic:
        return " [Harvey]"
    return ""


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    figures_dir = cfg["figures_dir"]
    color_col = get_prediction_color_feature(cfg)

    feature_summary = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_feature_shap_importance.csv")
    site_importance = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_site_importance.csv")
    site_summary = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_site_summary.tsv", sep="\t")
    site_stability = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_site_stability.tsv", sep="\t")
    substitution_feature_summary = pd.read_csv(
        f"{cfg['output_dir']}/{cfg['subtype']}_substitution_feature_shap_importance.csv"
    )
    substitution_site_importance = pd.read_csv(
        f"{cfg['output_dir']}/{cfg['subtype']}_substitution_site_importance.csv"
    )
    substitution_site_summary = pd.read_csv(
        f"{cfg['output_dir']}/{cfg['subtype']}_substitution_site_summary.tsv",
        sep="\t",
    )
    substitution_site_stability = pd.read_csv(
        f"{cfg['output_dir']}/{cfg['subtype']}_substitution_site_stability.tsv",
        sep="\t",
    )
    cv_results = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_cv_results.tsv", sep="\t")
    oof = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_oof_predictions.tsv", sep="\t")
    titers = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_matched_titers.tsv", sep="\t")
    neher_pred = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_neher_analog_predictions.tsv", sep="\t")
    neher_metrics = pd.read_csv(f"{cfg['output_dir']}/{cfg['subtype']}_neher_analog_metrics.tsv", sep="\t")

    cv_results["model_label"] = cv_results["model_type"].map(humanize_model_label)
    neher_metrics["model_label"] = neher_metrics["model_type"].map(humanize_model_label)
    neher_metrics["split_label"] = neher_metrics["split_type"].map(humanize_split_label)
    titers["homologous_source_label"] = titers["homologous_titer_source"].map(humanize_homologous_source)

    top_feature_summary = feature_summary.head(30).iloc[::-1].copy()
    top_feature_summary["label"] = top_feature_summary["feature"].map(humanize_feature_label)

    vmax = float(np.max(np.abs(top_feature_summary["mean_signed_shap"]))) or 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm")
    colors = [cmap(norm(v)) for v in top_feature_summary["mean_signed_shap"]]

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.barh(top_feature_summary["label"], top_feature_summary["mean_abs_shap"], color=colors)
    ax.set_xlabel("mean |SHAP|")
    ax.set_ylabel("")
    ax.set_title("top 30 SHAP features, site-state model")
    ax.title.set_fontfamily(MONO_FONT)
    ax.xaxis.label.set_fontfamily(MONO_FONT)
    for tick in ax.get_yticklabels():
        tick.set_fontfamily(MONO_FONT)
    for tick in ax.get_xticklabels():
        tick.set_fontfamily(MONO_FONT)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, pad=0.02)
    cb.set_label("mean signed SHAP")
    cb.ax.yaxis.label.set_fontfamily(MONO_FONT)
    for tick in cb.ax.get_yticklabels():
        tick.set_fontfamily(MONO_FONT)
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    shap_summary_path = f"{figures_dir}/{cfg['subtype']}_shap_summary_top30.pdf"
    fig.savefig(shap_summary_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    site_colors = []
    for _, row in site_importance.iterrows():
        neher = bool(row["named_in_neher2016"])
        wic = bool(row["named_in_wic2023"])
        if neher and wic:
            site_colors.append("#a50f15")
        elif neher:
            site_colors.append("#ef6548")
        elif wic:
            site_colors.append("#3182bd")
        else:
            site_colors.append("#bdbdbd")

    fig, ax = plt.subplots(figsize=(14, 4.2))
    ax.bar(site_importance["site"], site_importance["mean_abs_shap"], color=site_colors, width=1.0)
    ax.set_xlabel("HA1 position")
    ax.set_ylabel("mean |SHAP|")
    ax.set_title("site-level SHAP importance, site-state model")
    finish_axes(ax=ax, grid_axis="y")
    fig.tight_layout()
    site_bar_path = f"{figures_dir}/{cfg['subtype']}_site_importance_bar.pdf"
    fig.savefig(site_bar_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    top_sites = site_summary.nsmallest(30, "rank").sort_values("rank", ascending=False).copy()
    component_specs = [
        ("mean_abs_shap_change", "mean_signed_shap_change", "change"),
        ("mean_abs_shap_virus_aa", "mean_signed_shap_virus_aa", "virus aa"),
        ("mean_abs_shap_serum_aa", "mean_signed_shap_serum_aa", "serum aa"),
    ]
    signed_vals = []
    for _, signed_col, _ in component_specs:
        signed_vals.extend(top_sites[signed_col].tolist())
    vmax_components = float(np.max(np.abs(signed_vals))) if signed_vals else 1.0
    vmax_components = vmax_components or 1.0
    comp_norm = TwoSlopeNorm(vmin=-vmax_components, vcenter=0.0, vmax=vmax_components)

    fig, ax = plt.subplots(figsize=(9.5, 9.5))
    y = np.arange(len(top_sites))
    left = np.zeros(len(top_sites), dtype=float)
    for abs_col, signed_col, _ in component_specs:
        widths = top_sites[abs_col].to_numpy(dtype=float)
        segment_colors = [cmap(comp_norm(v)) for v in top_sites[signed_col].to_numpy(dtype=float)]
        ax.barh(y, widths, left=left, color=segment_colors, edgecolor="white", linewidth=0.6)
        left += widths

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{int(row['site'])}{site_flag_suffix(row)}" for _, row in top_sites.iterrows()]
    )
    ax.set_xlabel("mean |SHAP|")
    ax.set_ylabel("site")
    ax.set_title("top 30 sites by SHAP component, site-state model")
    comp_sm = plt.cm.ScalarMappable(cmap=cmap, norm=comp_norm)
    comp_sm.set_array([])
    cb = fig.colorbar(comp_sm, ax=ax, pad=0.02)
    cb.set_label("mean signed SHAP")
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    site_component_path = f"{figures_dir}/{cfg['subtype']}_site_component_top30.pdf"
    fig.savefig(site_component_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    stability_top = site_stability.nsmallest(30, "rank").sort_values("mean_rank", ascending=True).copy()
    freq_norm = plt.Normalize(0.0, 1.0)
    freq_cmap = plt.get_cmap("viridis")
    point_colors = [freq_cmap(freq_norm(v)) for v in stability_top["top30_frequency"].to_numpy(dtype=float)]
    fig, ax = plt.subplots(figsize=(10.5, 9.5))
    y = np.arange(len(stability_top))
    ax.hlines(
        y=y,
        xmin=stability_top["mean_rank"] - stability_top["sd_rank"],
        xmax=stability_top["mean_rank"] + stability_top["sd_rank"],
        color="#9e9e9e",
        linewidth=2,
    )
    ax.scatter(
        stability_top["mean_rank"],
        y,
        c=point_colors,
        s=55,
        edgecolors="black",
        linewidths=0.4,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{int(row['site'])}{site_flag_suffix(row)}" for _, row in stability_top.iterrows()]
    )
    ax.set_xlabel("mean fold rank +/- SD")
    ax.set_ylabel("site")
    ax.set_title("top 30 site stability across CV folds, site-state model")
    ax.invert_yaxis()
    freq_sm = plt.cm.ScalarMappable(cmap=freq_cmap, norm=freq_norm)
    freq_sm.set_array([])
    cb = fig.colorbar(freq_sm, ax=ax, pad=0.02)
    cb.set_label("top 30 frequency")
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    stability_path = f"{figures_dir}/{cfg['subtype']}_site_stability_top30.pdf"
    fig.savefig(stability_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 6))
    scatter = ax.scatter(
        oof["standardized_titer"],
        oof["predicted_standardized_titer_site_state"],
        c=oof[color_col],
        cmap="viridis",
        s=8,
        alpha=0.6,
        linewidths=0,
    )
    lo = float(np.nanmin([oof["standardized_titer"].min(), oof["predicted_standardized_titer_site_state"].min()]))
    hi = float(np.nanmax([oof["standardized_titer"].max(), oof["predicted_standardized_titer_site_state"].max()]))
    ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("observed standardized titer")
    ax.set_ylabel("predicted standardized titer")
    ax.set_title("predicted vs observed out-of-fold, site-state model")
    cb = fig.colorbar(scatter, ax=ax, pad=0.02)
    cb.set_label(humanize_distance_label(color_col))
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    scatter_path = f"{figures_dir}/{cfg['subtype']}_predicted_vs_actual.pdf"
    fig.savefig(scatter_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    sub_top_feature_summary = substitution_feature_summary.head(30).iloc[::-1].copy()
    sub_top_feature_summary["label"] = sub_top_feature_summary["feature"].map(humanize_feature_label)

    sub_vmax = float(np.max(np.abs(sub_top_feature_summary["mean_signed_shap"]))) or 1.0
    sub_norm = TwoSlopeNorm(vmin=-sub_vmax, vcenter=0.0, vmax=sub_vmax)
    sub_colors = [cmap(sub_norm(v)) for v in sub_top_feature_summary["mean_signed_shap"]]

    fig, ax = plt.subplots(figsize=(10.5, 9))
    ax.barh(sub_top_feature_summary["label"], sub_top_feature_summary["mean_abs_shap"], color=sub_colors)
    ax.set_xlabel("mean |SHAP|")
    ax.set_ylabel("")
    ax.set_title("top 30 SHAP features, substitution model")
    ax.title.set_fontfamily(MONO_FONT)
    ax.xaxis.label.set_fontfamily(MONO_FONT)
    for tick in ax.get_yticklabels():
        tick.set_fontfamily(MONO_FONT)
    for tick in ax.get_xticklabels():
        tick.set_fontfamily(MONO_FONT)
    sub_sm = plt.cm.ScalarMappable(cmap=cmap, norm=sub_norm)
    sub_sm.set_array([])
    cb = fig.colorbar(sub_sm, ax=ax, pad=0.02)
    cb.set_label("mean signed SHAP")
    cb.ax.yaxis.label.set_fontfamily(MONO_FONT)
    for tick in cb.ax.get_yticklabels():
        tick.set_fontfamily(MONO_FONT)
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    sub_shap_summary_path = f"{figures_dir}/{cfg['subtype']}_substitution_shap_summary_top30.pdf"
    fig.savefig(sub_shap_summary_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    sub_site_colors = []
    for _, row in substitution_site_importance.iterrows():
        neher = bool(row["named_in_neher2016"])
        wic = bool(row["named_in_wic2023"])
        if neher and wic:
            sub_site_colors.append("#a50f15")
        elif neher:
            sub_site_colors.append("#ef6548")
        elif wic:
            sub_site_colors.append("#3182bd")
        else:
            sub_site_colors.append("#bdbdbd")

    fig, ax = plt.subplots(figsize=(14, 4.2))
    ax.bar(
        substitution_site_importance["site"],
        substitution_site_importance["mean_abs_shap"],
        color=sub_site_colors,
        width=1.0,
    )
    ax.set_xlabel("HA1 position")
    ax.set_ylabel("mean |SHAP|")
    ax.set_title("site-level SHAP importance, substitution model")
    finish_axes(ax=ax, grid_axis="y")
    fig.tight_layout()
    sub_site_bar_path = f"{figures_dir}/{cfg['subtype']}_substitution_site_importance_bar.pdf"
    fig.savefig(sub_site_bar_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    sub_top_sites = substitution_site_summary.nsmallest(30, "rank").sort_values("rank", ascending=False).copy()
    sub_site_vmax = float(np.max(np.abs(sub_top_sites["mean_signed_shap"]))) or 1.0
    sub_site_norm = TwoSlopeNorm(vmin=-sub_site_vmax, vcenter=0.0, vmax=sub_site_vmax)
    sub_bar_colors = [cmap(sub_site_norm(v)) for v in sub_top_sites["mean_signed_shap"]]
    sub_labels = []
    for _, row in sub_top_sites.iterrows():
        suffix = site_flag_suffix(row)
        top_terms = row["top_substitution_terms"] if isinstance(row["top_substitution_terms"], str) else ""
        label = f"{int(row['site'])}{suffix}"
        if top_terms:
            label = f"{label} | {top_terms}"
        sub_labels.append(label)

    fig, ax = plt.subplots(figsize=(14.5, 9.5))
    ax.barh(sub_labels, sub_top_sites["mean_abs_shap"], color=sub_bar_colors)
    ax.set_xlabel("mean |SHAP|")
    ax.set_ylabel("site")
    ax.set_title("top 30 sites, substitution model")
    sub_site_sm = plt.cm.ScalarMappable(cmap=cmap, norm=sub_site_norm)
    sub_site_sm.set_array([])
    cb = fig.colorbar(sub_site_sm, ax=ax, pad=0.02)
    cb.set_label("mean signed SHAP")
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    sub_site_component_path = f"{figures_dir}/{cfg['subtype']}_substitution_site_component_top30.pdf"
    fig.savefig(sub_site_component_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    sub_stability_top = substitution_site_stability.nsmallest(30, "rank").sort_values("mean_rank", ascending=True).copy()
    point_colors = [freq_cmap(freq_norm(v)) for v in sub_stability_top["top30_frequency"].to_numpy(dtype=float)]
    fig, ax = plt.subplots(figsize=(10.5, 9.5))
    y = np.arange(len(sub_stability_top))
    ax.hlines(
        y=y,
        xmin=sub_stability_top["mean_rank"] - sub_stability_top["sd_rank"],
        xmax=sub_stability_top["mean_rank"] + sub_stability_top["sd_rank"],
        color="#9e9e9e",
        linewidth=2,
    )
    ax.scatter(
        sub_stability_top["mean_rank"],
        y,
        c=point_colors,
        s=55,
        edgecolors="black",
        linewidths=0.4,
        zorder=3,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{int(row['site'])}{site_flag_suffix(row)}" for _, row in sub_stability_top.iterrows()]
    )
    ax.set_xlabel("mean fold rank +/- SD")
    ax.set_ylabel("site")
    ax.set_title("top 30 site stability across CV folds, substitution model")
    ax.invert_yaxis()
    sub_freq_sm = plt.cm.ScalarMappable(cmap=freq_cmap, norm=freq_norm)
    sub_freq_sm.set_array([])
    cb = fig.colorbar(sub_freq_sm, ax=ax, pad=0.02)
    cb.set_label("top 30 frequency")
    finish_axes(ax=ax, grid_axis="x")
    fig.tight_layout()
    sub_stability_path = f"{figures_dir}/{cfg['subtype']}_substitution_site_stability_top30.pdf"
    fig.savefig(sub_stability_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 6))
    scatter = ax.scatter(
        oof["standardized_titer"],
        oof["predicted_standardized_titer_substitution"],
        c=oof[color_col],
        cmap="viridis",
        s=8,
        alpha=0.6,
        linewidths=0,
    )
    lo = float(np.nanmin([oof["standardized_titer"].min(), oof["predicted_standardized_titer_substitution"].min()]))
    hi = float(np.nanmax([oof["standardized_titer"].max(), oof["predicted_standardized_titer_substitution"].max()]))
    ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("observed standardized titer")
    ax.set_ylabel("predicted standardized titer")
    ax.set_title("predicted vs observed out-of-fold, substitution model")
    cb = fig.colorbar(scatter, ax=ax, pad=0.02)
    cb.set_label(humanize_distance_label(color_col))
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    sub_scatter_path = f"{figures_dir}/{cfg['subtype']}_substitution_predicted_vs_actual.pdf"
    fig.savefig(sub_scatter_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.8, 6.8))
    cv_plot = cv_results[
        cv_results["model_type"].isin(
            [
                "additive_only",
                "context_only",
                "randomized_binary_context",
                "binary_site_context",
                "site_state_context",
                "substitution_identity_context",
            ]
        )
    ].copy()
    model_order = [
        "context additive only",
        "context only",
        "randomized site changes + context",
        "site changed + context",
        "site aa state + context",
        "substitution identity + context",
    ]
    sns.boxplot(
        data=cv_plot,
        x="model_label",
        y="rmse",
        order=model_order,
        color="#9ecae1",
        width=0.55,
        showfliers=False,
        linewidth=0.9,
        ax=ax,
    )
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", va="top")
    ax.set_ylabel("RMSE")
    ax.set_xlabel("")
    ax.set_title("cross-validation RMSE by model")
    finish_axes(ax=ax, grid_axis="y")
    fig.tight_layout()
    cv_fig_path = f"{figures_dir}/{cfg['subtype']}_cv_performance_summary.pdf"
    fig.savefig(cv_fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.scatterplot(
        data=titers,
        x="standardized_titer",
        y="context_corrected_titer",
        hue="homologous_source_label",
        alpha=0.6,
        s=10,
        linewidth=0,
        ax=ax,
    )
    ax.set_title("standardization vs corrected titer")
    ax.set_xlabel("standardized titer")
    ax.set_ylabel("context corrected titer")
    legend = ax.legend(title="homologous titer source", frameon=True, framealpha=1.0)
    if legend is not None:
        legend.get_frame().set_facecolor("white")
        legend.get_frame().set_edgecolor("#444444")
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    diag_path = f"{figures_dir}/{cfg['subtype']}_residualization_diagnostic.pdf"
    fig.savefig(diag_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
    panel_map = [
        ("random_10pct_measurements", "A) random 10% measurements"),
        ("grouped_10pct_viruses", "B) all measurements for 10% viruses"),
    ]
    for ax, (split_name, title) in zip(axes, panel_map):
        sub = neher_pred[neher_pred["split_type"] == split_name]
        ax.scatter(
            sub["standardized_titer"],
            sub["predicted_standardized_titer_site_state"],
            s=8,
            alpha=0.6,
            color="#2b8cbe",
            edgecolors="none",
        )
        lo = min(sub["standardized_titer"].min(), sub["predicted_standardized_titer_site_state"].min())
        hi = max(sub["standardized_titer"].max(), sub["predicted_standardized_titer_site_state"].max())
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=1)
        rmse = np.sqrt(
            np.mean((sub["standardized_titer"] - sub["predicted_standardized_titer_site_state"]) ** 2)
        )
        ax.set_title(f"{title}\nRMSE={rmse:.3f}, n={len(sub)}", fontsize=10)
        ax.set_xlabel("observed standardized titer")
        finish_axes(ax=ax, grid_axis=None)
    axes[0].set_ylabel("predicted standardized titer")
    fig.suptitle("Neher and Bedford 2016 figure 2 analog", fontsize=12)
    fig.tight_layout()
    neher_fig_path = f"{figures_dir}/{cfg['subtype']}_neher2016_fig2_analog.pdf"
    fig.savefig(neher_fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    reciprocal = build_reciprocal_pairs(titers)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(
        reciprocal["standardized_ab"],
        reciprocal["standardized_ba"],
        s=8,
        alpha=0.45,
        color="#756bb1",
        linewidths=0,
    )
    lo = min(reciprocal["standardized_ab"].min(), reciprocal["standardized_ba"].min())
    hi = max(reciprocal["standardized_ab"].max(), reciprocal["standardized_ba"].max())
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("standardized titer A vs B")
    ax.set_ylabel("standardized titer B vs A")
    ax.set_title("reciprocal measurement symmetry")
    finish_axes(ax=ax, grid_axis=None)
    fig.tight_layout()
    reciprocal_path = f"{figures_dir}/{cfg['subtype']}_reciprocal_symmetry.pdf"
    fig.savefig(reciprocal_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    heldout_plot = neher_metrics[
        neher_metrics["model_type"].isin(
            [
                "additive_only",
                "binary_site_context",
                "site_state_context",
                "substitution_identity_context",
            ]
        )
    ].copy()
    hue_order = [
        "context additive only",
        "site changed + context",
        "site aa state + context",
        "substitution identity + context",
    ]
    split_order = ["random 10% measurements", "all measurements for 10% viruses"]
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
    ax.set_title("held-out prediction comparison")
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
    heldout_compare_path = f"{figures_dir}/{cfg['subtype']}_heldout_model_comparison.pdf"
    fig.savefig(heldout_compare_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    lines = [
        f"Saved SHAP summary figure: {shap_summary_path}",
        f"Saved site importance figure: {site_bar_path}",
        f"Saved site component figure: {site_component_path}",
        f"Saved site stability figure: {stability_path}",
        f"Saved predicted-vs-actual figure: {scatter_path}",
        f"Saved substitution SHAP summary figure: {sub_shap_summary_path}",
        f"Saved substitution site importance figure: {sub_site_bar_path}",
        f"Saved substitution site summary figure: {sub_site_component_path}",
        f"Saved substitution site stability figure: {sub_stability_path}",
        f"Saved substitution predicted-vs-actual figure: {sub_scatter_path}",
        f"Saved CV summary figure: {cv_fig_path}",
        f"Saved residualization diagnostic figure: {diag_path}",
        f"Saved Neher 2016 analog figure: {neher_fig_path}",
        f"Saved reciprocal symmetry figure: {reciprocal_path}",
        f"Saved held-out comparison figure: {heldout_compare_path}",
    ]
    append_qc_log(cfg, "31_generate_wic_figures", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
