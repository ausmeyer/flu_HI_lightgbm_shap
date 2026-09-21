#!/usr/bin/env python3
"""Generate supplementary summary tables for the manuscript SI."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_sites import HARVEY2023_PIP95_SITE_ROWS, HARVEY2023_RESTRICTED15_SITE_ROWS


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "manuscript" / "generated_supplement"
MASTER_TABLE = ROOT / "site_model_literature_comparison.tsv"

REUSABLE_SHAP_TABLES = [
    (
        ROOT / "H3N2" / "output" / "H3N2_reusable_shap_values.tsv",
        OUT_DIR / "table_s4_neher_bedford_shap_values.tsv",
    ),
    (
        ROOT
        / "H3N2-WIC-no-egg-no-mixed-no-unknown"
        / "output"
        / "modeling"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_reusable_shap_values.tsv",
        OUT_DIR / "table_s5_wic_filtered_shap_values.tsv",
    ),
]


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


def harvey_primary_n() -> int:
    return len({int(row["site"]) for row in HARVEY2023_PIP95_SITE_ROWS})


def harvey_restricted15_n() -> int:
    return len({int(row["site"]) for row in HARVEY2023_RESTRICTED15_SITE_ROWS})


def build_table_s2() -> pd.DataFrame:
    master = pd.read_csv(MASTER_TABLE, sep="\t")
    harvey_n = harvey_primary_n()
    restricted_n = harvey_restricted15_n()
    restricted_sites = {int(row["site"]) for row in HARVEY2023_RESTRICTED15_SITE_ROWS}
    rows: list[dict[str, object]] = []
    for analysis in ANALYSES:
        for model_label, rank_col in [("site-state", analysis["site_rank_col"]), ("substitution", analysis["sub_rank_col"])]:
            top_sites = set(master.loc[master[rank_col].notna() & (master[rank_col] <= 30), "site"])
            rows.append(
                {
                    "analysis": analysis["label"],
                    "model": model_label,
                    "koel_top30_overlap": overlap_string(master, rank_col, "in_koel", 7),
                    "neher_top30_overlap": overlap_string(master, rank_col, "in_neher2016", 15),
                    "harvey_top30_overlap": overlap_string(master, rank_col, "in_harvey2023", harvey_n),
                    "harvey_restricted15_top30_overlap": f"{len(top_sites & restricted_sites)}/{restricted_n}",
                    "shah_top30_overlap": overlap_string(master, rank_col, "in_shah2024", 30),
                }
            )
    return pd.DataFrame(rows)


def write_table_s1_tex(df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        r"\begin{table}[!h]",
        r"\centering",
        r"\caption{Cross-validation performance of the primary models and the models used in sensitivity analyses. We calculated RMSE and Spearman correlation in each held-out fold of cross-validation grouped by virus name and report the mean across folds.}",
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
        r"\caption{Overlap between each model's top 30 sites and the Koel (7 sites), Neher/Bedford (15 sites), Harvey structurally aware PIP~$\geq$0.95 (14 sites), and Shah (30 sites) reference sets.}",
        r"\label{tab:si-reference-overlap}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\renewcommand{\arraystretch}{1.2}",
        r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{0.25\textwidth}>{\raggedright\arraybackslash}p{0.15\textwidth}>{\centering\arraybackslash}p{0.09\textwidth}>{\centering\arraybackslash}p{0.14\textwidth}>{\centering\arraybackslash}p{0.16\textwidth}>{\centering\arraybackslash}p{0.09\textwidth}@{}}",
        r"\toprule",
        r"Data set & Model & Koel\newline sites & Neher/Bedford\newline sites & Harvey\newline PIP~$\geq$0.95 & Shah\newline sites \\",
        r"\midrule",
    ]
    for row in df.itertuples(index=False):
        lines.append(
            f"{row.analysis} & {row.model} & {row.koel_top30_overlap} & {row.neher_top30_overlap} & {row.harvey_top30_overlap} & {row.shah_top30_overlap} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# Site notes are restricted to sources inspected for this revision. Membership
# flags are concordance with published site lists, not experimental validation.
SITE_LITERATURE_NOTES = {
    144: "Harvey structurally-aware PIP >= 0.95 (site A); Shah aggregated seasonal set.",
    145: "Koel 2013 cluster-transition site. Koel 2019: substitutions at 145 were context-independent in that reverse-genetics panel.",
    155: "Koel 2013 cluster-transition site. Koel 2019: substitutions at 155 were context-dependent in that reverse-genetics panel.",
    156: "Koel 2013 cluster-transition site.",
    158: "Koel 2013 cluster-transition site; Harvey PIP >= 0.95 (site B). Zost 2017 examined the 158--160 glycosylation motif via K160T.",
    159: "Koel 2013 cluster-transition site; Harvey PIP >= 0.95 (site B).",
    189: "Koel 2013 cluster-transition site; Harvey PIP >= 0.95 (site B).",
    193: "Koel 2013 cluster-transition site; Harvey PIP >= 0.95 (site B).",
    225: "Harvey structurally-aware PIP >= 0.95 (receptor-binding site); Shah aggregated set. Chambers 2015 lists N225D among 3C.2a/3C.3a differences relative to A/Texas/50/2012.",
    135: "Harvey structurally-aware PIP >= 0.95 (site A); Neher/Bedford named substitution.",
    173: "Harvey structurally-aware PIP >= 0.95 (site D).",
    160: "Zost 2017: K160T introduces the HA1 158--160 N-linked glycosylation motif in contemporary H3N2.",
}


def yes_no(flag: bool) -> str:
    return "yes" if flag else "no"


def build_table_s6() -> pd.DataFrame:
    master = pd.read_csv(MASTER_TABLE, sep="\t")
    neher_top = set(master.loc[master["h3n2_site_state_rank"].notna() & (master["h3n2_site_state_rank"] <= 15), "site"])
    wic_top = set(master.loc[master["wic_filtered_site_state_rank"].notna() & (master["wic_filtered_site_state_rank"] <= 15), "site"])
    sites = sorted(neher_top | wic_top)
    rows: list[dict[str, object]] = []
    for site in sites:
        rec = master.loc[master["site"] == site].iloc[0]
        rows.append(
            {
                "site": int(site),
                "neher_site_state_rank": rec["h3n2_site_state_rank"] if pd.notna(rec["h3n2_site_state_rank"]) else pd.NA,
                "wic_filtered_site_state_rank": rec["wic_filtered_site_state_rank"] if pd.notna(rec["wic_filtered_site_state_rank"]) else pd.NA,
                "in_koel": bool(rec["in_koel"]) if pd.notna(rec["in_koel"]) else False,
                "in_neher2016": bool(rec["in_neher2016"]) if pd.notna(rec["in_neher2016"]) else False,
                "in_harvey_pip95": bool(rec["in_harvey2023"]) if pd.notna(rec["in_harvey2023"]) else False,
                "in_shah2024": bool(rec["in_shah2024"]) if pd.notna(rec["in_shah2024"]) else False,
                "literature_note": SITE_LITERATURE_NOTES.get(int(site), ""),
            }
        )
    return pd.DataFrame(rows)


def write_table_s6_tex(df: pd.DataFrame, out_path: Path) -> None:
    lines = [
        r"\begin{table}[!h]",
        r"\centering",
        r"\caption{Comparison of published site lists with the top 15 sites from the Neher/Bedford and passage-filtered WIC site-state models. We included every site that appeared in either model's list. The membership columns show whether a site appears in each published list, and the notes describe the findings from those studies.}",
        r"\label{tab:si-literature-overlap}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}rrrccccp{6.6cm}}",
        r"\toprule",
        r"Site & Neher rank & WIC filt.\ rank & Koel & Neher/Bedford & Harvey PIP~$\geq$0.95 & Shah & Note \\",
        r"\midrule",
    ]
    for row in df.itertuples(index=False):
        neher = "--" if pd.isna(row.neher_site_state_rank) else str(int(row.neher_site_state_rank))
        wic = "--" if pd.isna(row.wic_filtered_site_state_rank) else str(int(row.wic_filtered_site_state_rank))
        note = str(row.literature_note).replace("&", r"\&").replace(">=", r"$\geq$")
        lines.append(
            f"{int(row.site)} & {neher} & {wic} & {yes_no(row.in_koel)} & {yes_no(row.in_neher2016)} & {yes_no(row.in_harvey_pip95)} & {yes_no(row.in_shah2024)} & {note} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular*}",
            r"\end{table}",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_reusable_shap_tables() -> None:
    for src, dst in REUSABLE_SHAP_TABLES:
        if not src.exists():
            raise FileNotFoundError(
                f"Reusable SHAP table not found: {src}. "
                "Run the corresponding SHAP analysis before generating supplement tables."
            )
        table = pd.read_csv(src, sep="\t")
        expected_cols = ["site", "model_type", "covariate_name", "covariate_value", "shap_value"]
        missing = set(expected_cols).difference(table.columns)
        if missing:
            raise ValueError(f"{src} is missing required columns: {sorted(missing)}")
        table.to_csv(dst, sep="\t", index=False)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    table_s1 = build_table_s1()
    table_s2 = build_table_s2()

    table_s1.to_csv(OUT_DIR / "table_s1_model_performance.tsv", sep="\t", index=False)
    table_s2.to_csv(OUT_DIR / "table_s2_reference_overlap.tsv", sep="\t", index=False)
    write_table_s1_tex(table_s1, OUT_DIR / "table_s1_model_performance.tex")
    write_table_s2_tex(table_s2, OUT_DIR / "table_s2_reference_overlap.tex")
    table_s6 = build_table_s6()
    table_s6.to_csv(OUT_DIR / "table_s6_literature_overlap.tsv", sep="\t", index=False)
    write_table_s6_tex(table_s6, OUT_DIR / "table_s6_literature_overlap.tex")
    copy_reusable_shap_tables()


if __name__ == "__main__":
    main()
