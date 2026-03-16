#!/usr/bin/env python3
"""Build root-level master comparison tables across local models and literature site sets."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from paper_sites import NEHER2016_H3_SITE_ROWS, WIC2023_H3_SITE_ROWS


ROOT = Path(__file__).resolve().parents[1]
LITERATURE_DIR = ROOT / "relevant_literature"

KOEL_SITES = [145, 155, 156, 158, 159, 189, 193]

PLOS2015_GEOMETRIC_DNDS_GT1 = [96, 137, 138, 143, 222, 223, 225, 226]

VIRUS_EVOLUTION_2016_KOEL_TABLE = {
    145: {
        "dn_ds": 0.672,
        "dn_ds_percentile": 0.823,
        "inv_distance_r": 0.0820,
        "inv_distance_percentile": 0.883,
    },
    155: {
        "dn_ds": 0.0,
        "dn_ds_percentile": 0.0,
        "inv_distance_r": 0.0770,
        "inv_distance_percentile": 0.867,
    },
    156: {
        "dn_ds": 0.672,
        "dn_ds_percentile": 0.832,
        "inv_distance_r": 0.1317,
        "inv_distance_percentile": 0.971,
    },
    158: {
        "dn_ds": 1.36,
        "dn_ds_percentile": 0.958,
        "inv_distance_r": 0.1797,
        "inv_distance_percentile": 0.996,
    },
    159: {
        "dn_ds": 0.49,
        "dn_ds_percentile": 0.75,
        "inv_distance_r": 0.1837,
        "inv_distance_percentile": 1.0,
    },
    189: {
        "dn_ds": 0.474,
        "dn_ds_percentile": 0.763,
        "inv_distance_r": 0.0887,
        "inv_distance_percentile": 0.905,
    },
    193: {
        "dn_ds": 0.672,
        "dn_ds_percentile": 0.845,
        "inv_distance_r": 0.0980,
        "inv_distance_percentile": 0.936,
    },
}

LIGHTGBM_REVIEW_KOEL_RANKS = {
    145: {"passage_shap_rank": 28, "date_shap_rank": 66, "dn_ds_rank": 95, "leisr_rank": 44, "entropy_rank": 30},
    155: {"passage_shap_rank": 82, "date_shap_rank": 246, "dn_ds_rank": 281, "leisr_rank": 325, "entropy_rank": 69},
    156: {"passage_shap_rank": 102, "date_shap_rank": 122, "dn_ds_rank": 266, "leisr_rank": 198, "entropy_rank": 60},
    158: {"passage_shap_rank": 4, "date_shap_rank": 28, "dn_ds_rank": 66, "leisr_rank": 79, "entropy_rank": 44},
    159: {"passage_shap_rank": 3, "date_shap_rank": 8, "dn_ds_rank": 117, "leisr_rank": 204, "entropy_rank": 3},
    189: {"passage_shap_rank": 11, "date_shap_rank": 19, "dn_ds_rank": 120, "leisr_rank": 113, "entropy_rank": 42},
    193: {"passage_shap_rank": 16, "date_shap_rank": 9, "dn_ds_rank": 18, "leisr_rank": 60, "entropy_rank": 15},
}

LIGHTGBM_REVIEW_PASSAGE_TOP_FEATURE_TERMS = {
    91: ["91N"],
    92: ["92R"],
    144: ["144N"],
    158: ["158N", "158K"],
    159: ["159Y", "159S", "159F"],
    160: ["160K"],
    186: ["186G"],
    225: ["225D", "225N", "225G"],
    326: ["326R"],
    347: ["347K"],
    489: ["489D"],
}

LIGHTGBM_REVIEW_DATE_TOP_SITES = [
    142,
    3,
    144,
    121,
    193,
    261,
    225,
    198,
    311,
    159,
    212,
    278,
    131,
    33,
    128,
    346,
    135,
    45,
    489,
    189,
]

PLOS2015_SITE_TABLE_PATH = LITERATURE_DIR / "plos2015_ha_evolution_site_table.tsv"
VIRUS_EVOLUTION2016_SITE_TABLE_PATH = LITERATURE_DIR / "virus_evolution2016_passaging_site_table.tsv"
LIGHTGBM_REVIEW_SITE_TABLE_PATH = LITERATURE_DIR / "lightgbm_review_integrated_site_table.tsv"

MODEL_SPECS = [
    {
        "analysis_id": "h3n2",
        "analysis_label": "H3N2",
        "family": "site_state",
        "site_summary_path": ROOT / "H3N2/output/H3N2_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2/output/H3N2_cv_summary.csv",
        "cv_model_type": "site_state",
    },
    {
        "analysis_id": "h3n2",
        "analysis_label": "H3N2",
        "family": "substitution",
        "site_summary_path": ROOT / "H3N2/output/H3N2_substitution_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2/output/H3N2_cv_summary.csv",
        "cv_model_type": "substitution_identity",
    },
    {
        "analysis_id": "h3n2_patristic",
        "analysis_label": "H3N2-patristic",
        "family": "site_state",
        "site_summary_path": ROOT / "H3N2-patristic/output/H3N2_PATRISTIC_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-patristic/output/H3N2_PATRISTIC_cv_summary.csv",
        "cv_model_type": "site_state",
    },
    {
        "analysis_id": "h3n2_patristic",
        "analysis_label": "H3N2-patristic",
        "family": "substitution",
        "site_summary_path": ROOT / "H3N2-patristic/output/H3N2_PATRISTIC_substitution_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-patristic/output/H3N2_PATRISTIC_cv_summary.csv",
        "cv_model_type": "substitution_identity",
    },
    {
        "analysis_id": "wic_full",
        "analysis_label": "H3N2-WIC",
        "family": "site_state",
        "site_summary_path": ROOT / "H3N2-WIC/output/modeling/H3N2_WIC_HA1_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC/output/modeling/H3N2_WIC_HA1_cv_summary.tsv",
        "cv_model_type": "site_state_context",
    },
    {
        "analysis_id": "wic_full",
        "analysis_label": "H3N2-WIC",
        "family": "substitution",
        "site_summary_path": ROOT / "H3N2-WIC/output/modeling/H3N2_WIC_HA1_substitution_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC/output/modeling/H3N2_WIC_HA1_cv_summary.tsv",
        "cv_model_type": "substitution_identity_context",
    },
    {
        "analysis_id": "wic_patristic",
        "analysis_label": "H3N2-WIC-patristic",
        "family": "site_state",
        "site_summary_path": ROOT / "H3N2-WIC-patristic/output/modeling/H3N2_WIC_HA1_PATRISTIC_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC-patristic/output/modeling/H3N2_WIC_HA1_PATRISTIC_cv_summary.tsv",
        "cv_model_type": "site_state_context",
    },
    {
        "analysis_id": "wic_patristic",
        "analysis_label": "H3N2-WIC-patristic",
        "family": "substitution",
        "site_summary_path": ROOT / "H3N2-WIC-patristic/output/modeling/H3N2_WIC_HA1_PATRISTIC_substitution_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC-patristic/output/modeling/H3N2_WIC_HA1_PATRISTIC_cv_summary.tsv",
        "cv_model_type": "substitution_identity_context",
    },
    {
        "analysis_id": "wic_filtered",
        "analysis_label": "H3N2-WIC-no-egg-no-mixed-no-unknown",
        "family": "site_state",
        "site_summary_path": ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_cv_summary.tsv",
        "cv_model_type": "site_state_context",
    },
    {
        "analysis_id": "wic_filtered",
        "analysis_label": "H3N2-WIC-no-egg-no-mixed-no-unknown",
        "family": "substitution",
        "site_summary_path": ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_substitution_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_cv_summary.tsv",
        "cv_model_type": "substitution_identity_context",
    },
    {
        "analysis_id": "wic_filtered_patristic",
        "analysis_label": "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown",
        "family": "site_state",
        "site_summary_path": ROOT / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_cv_summary.tsv",
        "cv_model_type": "site_state_context",
    },
    {
        "analysis_id": "wic_filtered_patristic",
        "analysis_label": "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown",
        "family": "substitution",
        "site_summary_path": ROOT / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_substitution_site_summary.tsv",
        "cv_summary_path": ROOT / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_cv_summary.tsv",
        "cv_model_type": "substitution_identity_context",
    },
]


def load_site_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    keep_cols = ["site", "rank", "mean_abs_shap"]
    return df[keep_cols].copy()


def load_cv_summary(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix == ".tsv" else ","
    df = pd.read_csv(path, sep=sep)
    if pd.isna(df.iloc[0, 0]):
        df = df.iloc[1:].copy()
    df = df.rename(
        columns={
            "rmse": "rmse_mean",
            "rmse.1": "rmse_std",
            "mae": "mae_mean",
            "mae.1": "mae_std",
            "r2": "r2_mean",
            "r2.1": "r2_std",
            "spearman": "spearman_mean",
            "spearman.1": "spearman_std",
        }
    )
    numeric_cols = [c for c in df.columns if c != "model_type"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def load_plos2015_site_table() -> pd.DataFrame:
    df = pd.read_csv(PLOS2015_SITE_TABLE_PATH, sep="\t")
    df["site"] = pd.to_numeric(df["site"], errors="coerce").astype("Int64")
    for col in [
        "plos2015_wiley81_class",
        "plos2015_wiley87_class",
        "plos2015_bush99_class",
        "plos2015_shih07_class",
        "plos2015_pan11_class",
        "plos2015_meyer14_group",
        "plos2015_koel13_group",
    ]:
        df[col] = df[col].replace("-", pd.NA)
    df["plos2015_pan11_site"] = df["plos2015_pan11_class"].notna()
    df["plos2015_meyer14_site"] = df["plos2015_meyer14_group"].notna()
    df["plos2015_koel13_site"] = df["plos2015_koel13_group"].notna()
    return df


def load_virus_evolution2016_site_table() -> pd.DataFrame:
    df = pd.read_csv(VIRUS_EVOLUTION2016_SITE_TABLE_PATH, sep="\t")
    df["site"] = pd.to_numeric(df["site"], errors="coerce").astype("Int64")
    return df


def load_lightgbm_review_site_table() -> pd.DataFrame:
    df = pd.read_csv(LIGHTGBM_REVIEW_SITE_TABLE_PATH, sep="\t")
    df["site"] = pd.to_numeric(df["site"], errors="coerce").astype("Int64")
    df["lightgbm_review_gene_site"] = pd.to_numeric(
        df["lightgbm_review_gene_site"], errors="coerce"
    ).astype("Int64")
    return df


def build_literature_frame() -> pd.DataFrame:
    plos2015_df = load_plos2015_site_table()
    virus_evolution2016_df = load_virus_evolution2016_site_table()
    lightgbm_review_df = load_lightgbm_review_site_table()

    neher_df = pd.DataFrame(NEHER2016_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={
            "paper_category": "neher2016_category",
            "paper_note": "neher2016_note",
        }
    )
    neher_df["in_neher2016"] = True

    harvey_df = pd.DataFrame(WIC2023_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={
            "paper_category": "harvey2023_category",
            "paper_note": "harvey2023_note",
        }
    )
    harvey_df["in_harvey2023"] = True

    all_sites = sorted(
        set(KOEL_SITES)
        | set(PLOS2015_GEOMETRIC_DNDS_GT1)
        | set(VIRUS_EVOLUTION_2016_KOEL_TABLE)
        | set(LIGHTGBM_REVIEW_KOEL_RANKS)
        | set(LIGHTGBM_REVIEW_PASSAGE_TOP_FEATURE_TERMS)
        | set(LIGHTGBM_REVIEW_DATE_TOP_SITES)
        | set(plos2015_df["site"].dropna().astype(int))
        | set(virus_evolution2016_df["site"].dropna().astype(int))
        | set(lightgbm_review_df["site"].dropna().astype(int))
        | set(neher_df["site"])
        | set(harvey_df["site"])
    )
    lit = pd.DataFrame({"site": all_sites})

    lit["in_koel"] = lit["site"].isin(KOEL_SITES)
    lit = lit.merge(neher_df, on="site", how="left")
    lit = lit.merge(harvey_df, on="site", how="left")
    lit["in_neher2016"] = lit["in_neher2016"].fillna(False)
    lit["in_harvey2023"] = lit["in_harvey2023"].fillna(False)

    lit["plos2015_geometric_predicted_dn_ds_gt1"] = lit["site"].isin(PLOS2015_GEOMETRIC_DNDS_GT1)
    lit["plos2015_note"] = lit["plos2015_geometric_predicted_dn_ds_gt1"].map(
        lambda x: "Predicted dN/dS > 1 under the RSA + inverse-distance geometric model." if x else pd.NA
    )
    lit = lit.merge(plos2015_df, on="site", how="left")

    lit["virus_evolution2016_koel_table_site"] = lit["site"].isin(VIRUS_EVOLUTION_2016_KOEL_TABLE)
    for col in ["dn_ds", "dn_ds_percentile", "inv_distance_r", "inv_distance_percentile"]:
        lit[f"virus_evolution2016_{col}"] = lit["site"].map(
            lambda s, key=col: VIRUS_EVOLUTION_2016_KOEL_TABLE.get(int(s), {}).get(key, pd.NA)
        )
    lit = lit.merge(virus_evolution2016_df, on="site", how="left")

    lit["lightgbm_review_koel_table_site"] = lit["site"].isin(LIGHTGBM_REVIEW_KOEL_RANKS)
    for col in ["passage_shap_rank", "date_shap_rank", "dn_ds_rank", "leisr_rank", "entropy_rank"]:
        lit[f"lightgbm_review_{col}"] = lit["site"].map(
            lambda s, key=col: LIGHTGBM_REVIEW_KOEL_RANKS.get(int(s), {}).get(key, pd.NA)
        )
    lit = lit.merge(lightgbm_review_df, on="site", how="left")

    lit["lightgbm_review_passage_top_feature_count"] = lit["site"].map(
        lambda s: len(LIGHTGBM_REVIEW_PASSAGE_TOP_FEATURE_TERMS.get(int(s), []))
    )
    lit["lightgbm_review_passage_top_feature_terms"] = lit["site"].map(
        lambda s: ", ".join(LIGHTGBM_REVIEW_PASSAGE_TOP_FEATURE_TERMS.get(int(s), [])) or pd.NA
    )
    date_rank_map = {site: rank for rank, site in enumerate(LIGHTGBM_REVIEW_DATE_TOP_SITES, start=1)}
    lit["lightgbm_review_date_top_site_rank"] = lit["site"].map(date_rank_map).astype("Int64")

    return lit


def build_master_site_table() -> pd.DataFrame:
    lit = build_literature_frame()

    for spec in MODEL_SPECS:
        df = load_site_summary(spec["site_summary_path"])
        prefix = f"{spec['analysis_id']}_{spec['family']}"
        df = df.rename(
            columns={
                "rank": f"{prefix}_rank",
                "mean_abs_shap": f"{prefix}_mean_abs_shap",
            }
        )
        lit = lit.merge(df, on="site", how="outer")

    bool_cols = [
        "in_koel",
        "in_neher2016",
        "in_harvey2023",
        "plos2015_pan11_site",
        "plos2015_meyer14_site",
        "plos2015_koel13_site",
        "plos2015_geometric_predicted_dn_ds_gt1",
        "virus_evolution2016_koel_table_site",
        "lightgbm_review_koel_table_site",
    ]
    for col in bool_cols:
        lit[col] = lit[col].fillna(False)
    lit["lightgbm_review_passage_top_feature_count"] = (
        lit["lightgbm_review_passage_top_feature_count"].fillna(0).astype(int)
    )
    lit["lightgbm_review_date_top_site_rank"] = lit["lightgbm_review_date_top_site_rank"].astype("Int64")
    lit = lit.sort_values("site").reset_index(drop=True)
    return lit


def build_model_performance_table() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for spec in MODEL_SPECS:
        cv = load_cv_summary(spec["cv_summary_path"])
        match = cv.loc[cv["model_type"] == spec["cv_model_type"]]
        if match.empty:
            raise ValueError(f"Model type {spec['cv_model_type']!r} not found in {spec['cv_summary_path']}")
        rec = match.iloc[0].to_dict()
        rows.append(
            {
                "analysis_id": spec["analysis_id"],
                "analysis_label": spec["analysis_label"],
                "model_family": spec["family"],
                "cv_model_type": spec["cv_model_type"],
                "rmse_mean": rec["rmse_mean"],
                "rmse_std": rec["rmse_std"],
                "mae_mean": rec["mae_mean"],
                "mae_std": rec["mae_std"],
                "r2_mean": rec["r2_mean"],
                "r2_std": rec["r2_std"],
                "spearman_mean": rec["spearman_mean"],
                "spearman_std": rec["spearman_std"],
            }
        )

    review_rows = [
        {
            "analysis_id": "lightgbm_review",
            "analysis_label": "LightGBM review manuscript",
            "model_family": "passage_classifier",
            "cv_model_type": "multiclass",
            "rmse_mean": pd.NA,
            "rmse_std": pd.NA,
            "mae_mean": pd.NA,
            "mae_std": pd.NA,
            "r2_mean": pd.NA,
            "r2_std": pd.NA,
            "spearman_mean": pd.NA,
            "spearman_std": pd.NA,
            "overall_accuracy": 0.81,
            "balanced_accuracy": 0.77,
            "notes": "Passage classification on 39,121 HA sequences.",
        },
        {
            "analysis_id": "lightgbm_review",
            "analysis_label": "LightGBM review manuscript",
            "model_family": "date_regression",
            "cv_model_type": "regression",
            "rmse_mean": pd.NA,
            "rmse_std": pd.NA,
            "mae_mean": 74.5,
            "mae_std": pd.NA,
            "r2_mean": 0.98,
            "r2_std": pd.NA,
            "spearman_mean": pd.NA,
            "spearman_std": pd.NA,
            "overall_accuracy": pd.NA,
            "balanced_accuracy": pd.NA,
            "notes": "Date regression on 20,109 unpassaged HA sequences; MAE in days.",
        },
    ]

    perf = pd.DataFrame(rows)
    perf["overall_accuracy"] = pd.NA
    perf["balanced_accuracy"] = pd.NA
    perf["notes"] = pd.NA
    perf = pd.concat([perf, pd.DataFrame(review_rows)], ignore_index=True)
    return perf


def main() -> None:
    site_table = build_master_site_table()
    site_out = ROOT / "site_model_literature_comparison.tsv"
    site_table.to_csv(site_out, sep="\t", index=False)

    perf_table = build_model_performance_table()
    perf_out = ROOT / "model_performance_literature_summary.tsv"
    perf_table.to_csv(perf_out, sep="\t", index=False)

    print(f"Saved site comparison table: {site_out}")
    print(f"Saved model summary table: {perf_out}")


if __name__ == "__main__":
    main()
