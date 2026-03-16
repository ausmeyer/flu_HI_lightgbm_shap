#!/usr/bin/env python3
"""Build the PNAS Figure 4 cross-study concordance figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT / "site_model_literature_comparison.tsv"
OUT_DIR = ROOT / "69b778fd3e7b181fe1c2943b" / "figures"

# Full union of Neher 2016 and Harvey 2023 H3 site sets.
# Koel adds no extra sites because all 7 Koel positions are already in Neher 2016.
PANEL_A_SITES = [53, 62, 121, 126, 131, 135, 137, 140, 144, 145, 155, 156, 157, 158, 159, 160, 173, 186, 189, 193, 196, 276]

PANEL_A_MEMBERSHIP_COLS = [
    ("in_koel", "Koel"),
    ("in_neher2016", "Neher"),
    ("in_harvey2023", "Harvey"),
]

PANEL_A_RANK_COLS = [
    ("h3n2_site_state_rank", "H3N2"),
    ("wic_filtered_site_state_rank", "WIC filtered"),
    ("lightgbm_review_mean_shap_passage_rank", "LightGBM+SHAP\npassage"),
    ("lightgbm_review_mean_shap_date_rank", "LightGBM+SHAP\ndate"),
    ("lightgbm_review_mean_entropy_rank", "Unpassaged\nentropy"),
    ("lightgbm_review_dnds_rank", "Unpassaged\ndN/dS"),
    ("virus_evolution2016_unpassaged_tips_dnds_rank", "McWhite\nunpassaged dN/dS"),
]

PANEL_B_ROWS = [
    ("Koel", "in_koel", "set"),
    ("Neher", "in_neher2016", "set"),
    ("Harvey", "in_harvey2023", "set"),
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
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
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
    sub = df.set_index("site").loc[PANEL_A_SITES].reset_index()
    for col, _ in PANEL_A_RANK_COLS:
        sub[f"{col}_bin"] = sub[col].map(rank_bin)
    return sub


def build_panel_b_table(df: pd.DataFrame) -> pd.DataFrame:
    h3_sites = set(df.loc[df["h3n2_site_state_rank"].notna() & (df["h3n2_site_state_rank"] <= 30), "site"])
    wic_sites = set(
        df.loc[
            df["wic_filtered_site_state_rank"].notna() & (df["wic_filtered_site_state_rank"] <= 30),
            "site",
        ]
    )

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
    panel_a_out = OUT_DIR / "fig4_panelA_sites.tsv"
    panel_b_out = OUT_DIR / "fig4_panelB_overlaps.tsv"
    panel_a.to_csv(panel_a_out, sep="\t", index=False)
    panel_b.to_csv(panel_b_out, sep="\t", index=False)


def plot_panel_a(ax: plt.Axes, panel_a: pd.DataFrame) -> None:
    membership_cols = [c for c, _ in PANEL_A_MEMBERSHIP_COLS]
    rank_bin_cols = [f"{c}_bin" for c, _ in PANEL_A_RANK_COLS]

    membership_mat = panel_a[membership_cols].fillna(False).astype(int).to_numpy()
    rank_mat = panel_a[rank_bin_cols].to_numpy()

    membership_cmap = mpl.colors.ListedColormap(["#ffffff", "#1f1f1f"])
    # Muted, print-friendly palette: strong signal in deep blue, weaker signal fades toward cool gray.
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

    ax.scatter(plot_df["h3n2_overlap"], y, s=42, color=h3_color, label="H3N2 top 30 overlap", zorder=3)
    ax.scatter(plot_df["wic_filtered_overlap"], y, s=42, color=wic_color, label="Filtered WIC top 30 overlap", zorder=3)

    labels = []
    for _, row in plot_df.iterrows():
        if row["label"] in {"Koel", "Neher", "Harvey"}:
            labels.append(f"{row['label']} (n={int(row['ref_size'])})")
        else:
            labels.append(f"{row['label']} (top 30)")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", labelright=False, labelleft=True, length=0, pad=4)
    ax.set_xlim(0, 30)
    ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
    ax.set_xlabel("Overlap with top-30 sites")
    ax.grid(axis="x", color="#e6e6e6", lw=0.8)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=h3_color, markersize=6, label="H3N2"),
        Line2D([0], [0], marker="o", linestyle="", color=wic_color, markersize=6, label="Filtered WIC"),
    ]
    legend = ax.legend(handles=handles, frameon=True, loc="lower right", fontsize=7)
    legend.get_frame().set_facecolor("white")
    legend.get_frame().set_edgecolor("#d0d0d0")
    legend.get_frame().set_linewidth(0.8)
    legend.get_frame().set_alpha(1.0)
def add_legends(fig: plt.Figure, ax_a: plt.Axes) -> None:
    bbox_a = ax_a.get_position()
    panel_a_center_x = (bbox_a.x0 + bbox_a.x1) / 2
    legend_y = bbox_a.y0 - 0.035

    rank_handles = [
        Patch(facecolor="#16324f", edgecolor="none", label="Top 10"),
        Patch(facecolor="#4f6d8a", edgecolor="none", label="11-30"),
        Patch(facecolor="#b7c5d3", edgecolor="none", label="31-100"),
        Patch(facecolor="#e4e7eb", edgecolor="none", label=">100"),
    ]
    rank_legend = fig.legend(
        handles=rank_handles,
        loc="lower center",
        bbox_to_anchor=(panel_a_center_x - 0.02, legend_y),
        ncol=4,
        frameon=True,
        fontsize=8,
        title="Rank bin",
        title_fontsize=8,
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
    fig.text(bbox_a.x0 - 0.025, bbox_a.y1 + 0.01, "A", fontsize=12, fontweight="bold", va="bottom")
    # Place B clearly to the left of the panel-B y tick labels, not just the plotting area.
    fig.text(bbox_b.x0 - 0.135, bbox_b.y1 + 0.01, "B", fontsize=12, fontweight="bold", va="bottom")


def main() -> None:
    configure_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(TABLE_PATH, sep="\t")
    panel_a = build_panel_a_table(df)
    panel_b = build_panel_b_table(df)
    save_supporting_tables(panel_a, panel_b)

    fig = plt.figure(figsize=(11.8, 6.0))
    gs = GridSpec(1, 2, width_ratios=[1.42, 1.0], wspace=0.56, figure=fig)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    plot_panel_a(ax_a, panel_a)
    plot_panel_b(ax_b, panel_b)
    add_legends(fig, ax_a)
    add_panel_labels(fig, ax_a, ax_b)

    fig.suptitle("Cross-study concordance of primary H3N2 antigenic-site models", x=0.49, y=0.995, fontsize=11)
    fig.subplots_adjust(bottom=0.16, top=0.76, left=0.08, right=0.97)

    pdf_path = OUT_DIR / "fig4_cross_study_concordance.pdf"
    png_path = OUT_DIR / "fig4_cross_study_concordance.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    print(f"Saved {pdf_path}")
    print(f"Saved {png_path}")


if __name__ == "__main__":
    main()
