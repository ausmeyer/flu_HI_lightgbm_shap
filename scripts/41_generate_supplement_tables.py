#!/usr/bin/env python3
"""Generate supplementary summary tables for the manuscript SI."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "manuscript" / "generated_supplement"
MASTER_TABLE = ROOT / "site_model_literature_comparison.tsv"


ANALYSES = [
    {
        "label": "Neher/Bedford data",
        "cv_path": ROOT / "H3N2" / "output" / "H3N2_cv_summary.csv",
        "sep": ",",
        "site_rank_col": "h3n2_site_state_rank",
        "sub_rank_col": "h3n2_substitution_rank",
        "site_model_type": "site_state",
        "sub_model_type": "substitution_identity",
    },
    {
        "label": "Neher/Bedford data + patristic",
        "cv_path": ROOT / "H3N2-patristic" / "output" / "H3N2_PATRISTIC_cv_summary.csv",
        "sep": ",",
        "site_rank_col": "h3n2_patristic_site_state_rank",
        "sub_rank_col": "h3n2_patristic_substitution_rank",
        "site_model_type": "site_state",
        "sub_model_type": "substitution_identity",
    },
    {
        "label": "WIC data",
        "cv_path": ROOT / "H3N2-WIC" / "output" / "modeling" / "H3N2_WIC_HA1_cv_summary.tsv",
        "sep": "\t",
        "site_rank_col": "wic_full_site_state_rank",
        "sub_rank_col": "wic_full_substitution_rank",
        "site_model_type": "site_state_context",
        "sub_model_type": "substitution_identity_context",
    },
    {
        "label": "WIC filtered data",
        "cv_path": ROOT
        / "H3N2-WIC-no-egg-no-mixed-no-unknown"
        / "output"
        / "modeling"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_cv_summary.tsv",
        "sep": "\t",
        "site_rank_col": "wic_filtered_site_state_rank",
        "sub_rank_col": "wic_filtered_substitution_rank",
        "site_model_type": "site_state_context",
        "sub_model_type": "substitution_identity_context",
    },
    {
        "label": "WIC data + patristic",
        "cv_path": ROOT / "H3N2-WIC-patristic" / "output" / "modeling" / "H3N2_WIC_HA1_PATRISTIC_cv_summary.tsv",
        "sep": "\t",
        "site_rank_col": "wic_patristic_site_state_rank",
        "sub_rank_col": "wic_patristic_substitution_rank",
        "site_model_type": "site_state_context",
        "sub_model_type": "substitution_identity_context",
    },
    {
        "label": "WIC filtered data + patristic",
        "cv_path": ROOT
        / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown"
        / "output"
        / "modeling"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_cv_summary.tsv",
        "sep": "\t",
        "site_rank_col": "wic_filtered_patristic_site_state_rank",
        "sub_rank_col": "wic_filtered_patristic_substitution_rank",
        "site_model_type": "site_state_context",
        "sub_model_type": "substitution_identity_context",
    },
]


def format_float(value: float) -> str:
    return f"{value:.3f}"


def load_cv_metrics(path: Path, sep: str, model_type: str) -> tuple[float, float]:
    df = pd.read_csv(path, sep=sep)
    df = df[df["model_type"].notna()].copy()
    row = df.loc[df["model_type"] == model_type].iloc[0]
    return float(row["rmse"]), float(row["spearman"])


def build_table_s1() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for analysis in ANALYSES:
        site_rmse, site_spearman = load_cv_metrics(
            analysis["cv_path"], analysis["sep"], analysis["site_model_type"]
        )
        sub_rmse, sub_spearman = load_cv_metrics(
            analysis["cv_path"], analysis["sep"], analysis["sub_model_type"]
        )
        rows.append(
            {
                "analysis": analysis["label"],
                "model": "site-state",
                "rmse": site_rmse,
                "spearman": site_spearman,
            }
        )
        rows.append(
            {
                "analysis": analysis["label"],
                "model": "substitution",
                "rmse": sub_rmse,
                "spearman": sub_spearman,
            }
        )
    return pd.DataFrame(rows)


def overlap_string(df: pd.DataFrame, rank_col: str, membership_col: str, denom: int) -> str:
    top_sites = set(df.loc[df[rank_col].notna() & (df[rank_col] <= 30), "site"])
    ref_sites = set(df.loc[df[membership_col].fillna(False).astype(bool), "site"])
    return f"{len(top_sites & ref_sites)}/{denom}"


def build_table_s2() -> pd.DataFrame:
    master = pd.read_csv(MASTER_TABLE, sep="\t")
    rows: list[dict[str, object]] = []
    for analysis in ANALYSES:
        for model_label, rank_col in [("site-state", analysis["site_rank_col"]), ("substitution", analysis["sub_rank_col"])]:
            rows.append(
                {
                    "analysis": analysis["label"],
                    "model": model_label,
                    "koel_top30_overlap": overlap_string(master, rank_col, "in_koel", 7),
                    "neher_top30_overlap": overlap_string(master, rank_col, "in_neher2016", 15),
                    "harvey_top30_overlap": overlap_string(master, rank_col, "in_harvey2023", 15),
                    "shah_top30_overlap": overlap_string(master, rank_col, "in_shah2024", 30),
                }
            )
    return pd.DataFrame(rows)


def write_table_s1_tex(df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        r"\begin{table}[!h]",
        r"\centering",
        r"\caption{Cross-validation performance of the primary and sensitivity-analysis models included in the manuscript and SI. Values are mean grouped-cross-validation performance across folds.}",
        r"\label{tab:si-model-performance}",
        r"\small",
        r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}llrr}",
        r"\toprule",
        r"Analysis & Model & RMSE & Spearman \\",
        r"\midrule",
    ]
    for row in df.itertuples(index=False):
        lines.append(
            f"{row.analysis} & {row.model} & {format_float(row.rmse)} & {format_float(row.spearman)} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular*}", r"\end{table}"])
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_table_s2_tex(df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        r"\begin{table}[!h]",
        r"\centering",
        r"\caption{Top-30 overlap between model-ranked sites from each data set and the Koel, Neher/Bedford, Harvey/WIC, and Shah/WIC reference site sets. The Neher/Bedford and Harvey/WIC site sets each contain 15 sites, the Shah/WIC site set contains 30 sites, and the Koel site set contains 7 sites.}",
        r"\label{tab:si-reference-overlap}",
        r"\small",
        r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}llcccc}",
        r"\toprule",
        r"Data set & Model & Koel sites & Neher/Bedford sites & Harvey/WIC sites & Shah/WIC sites \\",
        r"\midrule",
    ]
    for row in df.itertuples(index=False):
        lines.append(
            f"{row.analysis} & {row.model} & {row.koel_top30_overlap} & {row.neher_top30_overlap} & {row.harvey_top30_overlap} & {row.shah_top30_overlap} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular*}", r"\end{table}"])
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table_s1 = build_table_s1()
    table_s2 = build_table_s2()

    table_s1.to_csv(OUT_DIR / "table_s1_model_performance.tsv", sep="\t", index=False)
    table_s2.to_csv(OUT_DIR / "table_s2_reference_overlap.tsv", sep="\t", index=False)
    write_table_s1_tex(table_s1, OUT_DIR / "table_s1_model_performance.tex")
    write_table_s2_tex(table_s2, OUT_DIR / "table_s2_reference_overlap.tex")


if __name__ == "__main__":
    main()
