#!/usr/bin/env python3
"""Build the manuscript cross-study concordance figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT / "site_model_literature_comparison.tsv"
OUT_DIR = ROOT / "manuscript" / "generated_figures"

PANEL_A_MEMBERSHIP_COLS = [
    ("in_koel", "Koel\nsites"),
    ("in_neher2016", "Neher/Bedford\nsites"),
    ("in_harvey2023", "Harvey\nPIP≥0.95"),
    ("in_shah2024", "Shah/WIC\nsites"),
]

PANEL_A_RANK_COLS = [
    ("h3n2_site_state_rank", "Neher/Bedford\ndata"),
    ("wic_filtered_site_state_rank", "Harvey/WIC\nfiltered data"),
    ("lightgbm_review_mean_shap_passage_rank", "LightGBM+SHAP\npassage"),
    ("lightgbm_review_mean_shap_date_rank", "LightGBM+SHAP\ndate"),
    ("lightgbm_review_mean_entropy_rank", "Unpassaged\nentropy"),
    ("lightgbm_review_dnds_rank", "Unpassaged\ndN/dS"),
    ("virus_evolution2016_unpassaged_tips_dnds_rank", "McWhite\nunpassaged dN/dS"),
]

PANEL_B_ROWS = [
    ("Koel sites", "in_koel", "set"),
    ("Neher/Bedford sites", "in_neher2016", "set"),
    ("Harvey PIP≥0.95", "in_harvey2023", "set"),
    ("Shah/WIC sites", "in_shah2024", "set"),
    ("LightGBM+SHAP passage", "lightgbm_review_mean_shap_passage_rank", "rank"),
    ("LightGBM+SHAP date", "lightgbm_review_mean_shap_date_rank", "rank"),
    ("Unpassaged entropy", "lightgbm_review_mean_entropy_rank", "rank"),
    ("Unpassaged dN/dS", "lightgbm_review_dnds_rank", "rank"),
    ("McWhite unpassaged dN/dS", "virus_evolution2016_unpassaged_tips_dnds_rank", "rank"),
]


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "legend.fontsize": 11,
            "legend.title_fontsize": 11,
        }
    )


def rank_bin(value: float | int | None) -> int:
    if pd.isna(value):
        return 4
    if value <= 10:
        return 0
    if value <= 30:
        return 1
    if value <= 100:
        return 2
    return 3


def build_panel_a_table(df: pd.DataFrame) -> pd.DataFrame:
    mask = np.zeros(len(df), dtype=bool)
    for col, _ in PANEL_A_MEMBERSHIP_COLS:
        mask |= df[col].fillna(False).to_numpy(dtype=bool)
    sub = df.loc[mask].sort_values("site").reset_index(drop=True)
    for col, _ in PANEL_A_RANK_COLS:
        sub[f"{col}_bin"] = sub[col].map(rank_bin)
    return sub


def build_panel_b_table(df: pd.DataFrame) -> pd.DataFrame:
    h3_sites = set(df.loc[df["h3n2_site_state_rank"].notna() & (df["h3n2_site_state_rank"] <= 30), "site"])
    wic_sites = set(df.loc[df["wic_filtered_site_state_rank"].notna() & (df["wic_filtered_site_state_rank"] <= 30), "site"])

    rows: list[dict[str, object]] = []
    for label, col, kind in PANEL_B_ROWS:
        if kind == "set":
            ref_sites = set(df.loc[df[col].fillna(False).astype(bool), "site"])
        else:
            ref_sites = set(df.loc[df[col].notna() & (df[col] <= 30), "site"])
        rows.append(
            {
                "label": label,
                "h3n2_overlap": len(h3_sites & ref_sites),
                "wic_filtered_overlap": len(wic_sites & ref_sites),
                "ref_size": len(ref_sites),
            }
        )
    return pd.DataFrame(rows)


def save_supporting_tables(panel_a: pd.DataFrame, panel_b: pd.DataFrame) -> None:
    panel_a.to_csv(OUT_DIR / "fig4_panelA_sites.tsv", sep="\t", index=False)
    panel_b.to_csv(OUT_DIR / "fig4_panelB_overlaps.tsv", sep="\t", index=False)


def plot_panel_a(ax: plt.Axes, panel_a: pd.DataFrame) -> None:
    membership_cols = [c for c, _ in PANEL_A_MEMBERSHIP_COLS]
    rank_bin_cols = [f"{c}_bin" for c, _ in PANEL_A_RANK_COLS]

    membership_mat = panel_a[membership_cols].fillna(False).astype(int).to_numpy()
    rank_mat = panel_a[rank_bin_cols].to_numpy()

    membership_cmap = mpl.colors.ListedColormap(["#ffffff", "#1f1f1f"])
    rank_cmap = mpl.colors.ListedColormap(["#16324f", "#4f6d8a", "#b7c5d3", "#e4e7eb", "#ffffff"])
    membership_x = np.arange(membership_mat.shape[1] + 1) - 0.5
    rank_x = np.arange(rank_mat.shape[1] + 1) + membership_mat.shape[1] - 0.5
    y = np.arange(len(panel_a) + 1) - 0.5

    membership_norm = mpl.colors.BoundaryNorm([-0.5, 0.5, 1.5], membership_cmap.N)
    rank_norm = mpl.colors.BoundaryNorm(np.arange(-0.5, 5.5, 1), rank_cmap.N)

    ax.pcolormesh(
        membership_x,
        y,
        membership_mat,
        cmap=membership_cmap,
        norm=membership_norm,
        shading="flat",
        edgecolors="#ffffff",
        linewidth=0.8,
        antialiased=False,
    )
    ax.pcolormesh(
        rank_x,
        y,
        rank_mat,
        cmap=rank_cmap,
        norm=rank_norm,
        shading="flat",
        edgecolors="#ffffff",
        linewidth=0.8,
        antialiased=False,
    )

    total_cols = membership_mat.shape[1] + rank_mat.shape[1]
    ax.set_xlim(-0.5, total_cols - 0.5)
    ax.set_ylim(len(panel_a) - 0.5, -0.5)
    ax.set_xticks(np.arange(total_cols))
    ax.set_xticklabels(
        [label for _, label in PANEL_A_MEMBERSHIP_COLS] + [label for _, label in PANEL_A_RANK_COLS],
        rotation=90,
        va="bottom",
        ha="center",
    )
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=6)
    ax.tick_params(axis="y", length=0)
    ax.set_yticks(np.arange(len(panel_a)))
    ax.set_yticklabels(panel_a["site"].astype(int).astype(str))
    ax.set_ylabel("HA site")


def plot_panel_b(ax: plt.Axes, panel_b: pd.DataFrame) -> None:
    plot_df = panel_b.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_df))

    h3_color = "#1d3557"
    wic_color = "#c05a2b"

    for i, row in plot_df.iterrows():
        ax.plot([row["h3n2_overlap"], row["wic_filtered_overlap"]], [i, i], color="#d3d3d3", lw=1.5, zorder=1)

    ax.scatter(
        plot_df["h3n2_overlap"],
        y,
        s=42,
        color=h3_color,
        label="Neher/Bedford data top-30 overlap",
        zorder=3,
    )
    ax.scatter(
        plot_df["wic_filtered_overlap"],
        y,
        s=42,
        color=wic_color,
        label="WIC filtered data top-30 overlap",
        zorder=3,
    )

    labels = []
    for _, row in plot_df.iterrows():
        if row["label"] in {"Koel sites", "Neher/Bedford sites", "Harvey PIP≥0.95", "Shah/WIC sites"}:
            labels.append(f"{row['label']} (n={int(row['ref_size'])})")
        else:
            labels.append(f"{row['label']} (top 30)")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", labelright=False, labelleft=True, length=0, pad=4)
    ax.set_xlim(0, 30)
    ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
    ax.set_xlabel("Overlap between data top-30 sites and reference site sets")
    ax.grid(axis="x", color="#e6e6e6", lw=0.8)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=h3_color, markersize=6, label="Neher/Bedford data"),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            color=wic_color,
            markersize=6,
            label="WIC filtered data",
        ),
    ]
    legend = ax.legend(handles=handles, frameon=True, loc="lower right", fontsize=11)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#d0d0d0")
    legend.get_frame().set_linewidth(0.8)
    legend.get_frame().set_alpha(1.0)


def add_legends(fig: plt.Figure, ax_a: plt.Axes) -> None:
    bbox_a = ax_a.get_position()
    panel_a_center_x = (bbox_a.x0 + bbox_a.x1) / 2
    legend_y = bbox_a.y0 - 0.095

    rank_handles = [
        Patch(facecolor="#16324f", edgecolor="none", label="Top 10"),
        Patch(facecolor="#4f6d8a", edgecolor="none", label="11-30"),
        Patch(facecolor="#b7c5d3", edgecolor="none", label="31-100"),
        Patch(facecolor="#e4e7eb", edgecolor="none", label=">100"),
    ]
    rank_legend = fig.legend(
        handles=rank_handles,
        loc="lower center",
        bbox_to_anchor=(panel_a_center_x, legend_y),
        ncol=4,
        frameon=True,
        fontsize=11,
        title="Rank bin",
        title_fontsize=11,
        handlelength=1.3,
        handletextpad=0.5,
        borderpad=0.7,
        columnspacing=1.2,
    )
    rank_legend.get_frame().set_facecolor("white")
    rank_legend.get_frame().set_edgecolor("#d0d0d0")
    rank_legend.get_frame().set_linewidth(0.8)
    rank_legend.get_frame().set_alpha(1.0)


def add_panel_labels(fig: plt.Figure, ax_a: plt.Axes, ax_b: plt.Axes) -> None:
    bbox_a = ax_a.get_position()
    bbox_b = ax_b.get_position()
    fig.text(bbox_a.x0 - 0.06, bbox_a.y1 + 0.21, "A", fontsize=16, fontweight="bold", va="bottom")
    fig.text(bbox_b.x0 - 0.23, bbox_b.y1 + 0.21, "B", fontsize=16, fontweight="bold", va="bottom")


def main() -> None:
    configure_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(TABLE_PATH, sep="\t").sort_values("site").reset_index(drop=True)
    panel_a = build_panel_a_table(df)
    panel_b = build_panel_b_table(df)
    save_supporting_tables(panel_a, panel_b)

    fig = plt.figure(figsize=(16.3, 9.2))
    gs = fig.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[1.88, 1.14],
        left=0.08,
        right=0.988,
        top=0.88,
        bottom=0.15,
        wspace=0.68,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    plot_panel_a(ax_a, panel_a)
    plot_panel_b(ax_b, panel_b)
    add_legends(fig, ax_a)
    add_panel_labels(fig, ax_a, ax_b)

    pdf_path = OUT_DIR / "fig4_cross_study_concordance.pdf"
    png_path = OUT_DIR / "fig4_cross_study_concordance.png"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
