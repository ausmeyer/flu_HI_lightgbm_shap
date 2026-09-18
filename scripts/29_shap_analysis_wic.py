#!/usr/bin/env python3
"""Compute sampled cross-validated SHAP summaries and site stability for the WIC HA1 model."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
import shap
from scipy.stats import mannwhitneyu
from sklearn.model_selection import GroupKFold

from model_validation import inner_selection, fit_selected_model

from common import (
    append_qc_log,
    apply_additive_effects_multi,
    fit_additive_effects_multi,
    get_distance_feature_cols,
    get_lightgbm_params,
    load_config,
)
from paper_sites import NEHER2016_H3_SITE_ROWS, SHAH2024_H3_SITE_ROWS, WIC2023_H3_SITE_ROWS
from shap_reuse import (
    collect_site_state_reuse_stats,
    collect_substitution_reuse_stats,
    new_reuse_stats,
    reuse_stats_to_table,
)


ADDITIVE_FACTOR_COLS = [
    "serumSequenceId",
    "assay_date_label",
    "virus_passage_class",
    "reference_passage_class",
    "antiserum_class",
    "rbc_type",
]

CONTEXT_COLS = [
    "temporal_distance",
    "virus_passage_class",
    "reference_passage_class",
    "virus_passage_known",
    "reference_passage_known",
    "antiserum_class",
    "rbc_type",
    "rbc_species",
    "oseltamivir_present",
    "oseltamivir_nm",
]


def get_context_cols(cfg: dict) -> list[str]:
    context_cols = CONTEXT_COLS.copy()
    for col in get_distance_feature_cols(cfg):
        if col not in context_cols:
            context_cols.append(col)
    return context_cols


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    parser.add_argument("--n-splits", type=int, default=None)
    parser.add_argument("--max-rows-per-fold", type=int, default=2500)
    return parser.parse_args()


def prepare_frame(df: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    X = df[feature_cols].copy()
    for col in X.columns:
        if not pd.api.types.is_numeric_dtype(X[col]):
            X[col] = X[col].astype("category")
    return X


def get_site_state_columns(df: pd.DataFrame) -> tuple[list[str], list[str], list[str]]:
    change_cols = sorted(
        [c for c in df.columns if c.startswith("change_")],
        key=lambda c: int(c.split("_")[1]),
    )
    virus_aa_cols = sorted(
        [c for c in df.columns if re.fullmatch(r"[A-Z]\d+", c)],
        key=lambda c: (int(c[1:]), c[0]),
    )
    serum_aa_cols = sorted(
        [c for c in df.columns if re.fullmatch(r"\d+[A-Z]", c)],
        key=lambda c: (int(c[:-1]), c[-1]),
    )
    return change_cols, virus_aa_cols, serum_aa_cols


def fit_split_additive(
    train_df: pd.DataFrame,
    val_df: pd.DataFrame,
    target_col: str,
) -> tuple[np.ndarray, np.ndarray]:
    additive = fit_additive_effects_multi(
        df=train_df,
        target_col=target_col,
        factor_cols=ADDITIVE_FACTOR_COLS,
    )
    train_base = apply_additive_effects_multi(
        train_df,
        factor_cols=ADDITIVE_FACTOR_COLS,
        intercept=float(additive["intercept"]),
        factor_effects=additive["factor_effects"],
    )
    val_base = apply_additive_effects_multi(
        val_df,
        factor_cols=ADDITIVE_FACTOR_COLS,
        intercept=float(additive["intercept"]),
        factor_effects=additive["factor_effects"],
    )
    return train_base, val_base


def top_terms_from_feature_means(feature_means: pd.Series, cols: list[str], n: int = 3) -> str:
    if not cols:
        return ""
    ranked = feature_means.loc[cols].sort_values(ascending=False)
    ranked = ranked[ranked > 0]
    return ", ".join(ranked.head(n).index.tolist())


def parse_substitution_site(feature: str) -> int | None:
    match = re.fullmatch(r"sub_(\d+)_([A-Z])_([A-Z])", feature)
    if match is None:
        return None
    return int(match.group(1))


def build_comparison_table(
    site_importance: pd.DataFrame,
    reference_rows: list[dict[str, object]],
    top_n: int,
    named_col: str,
    category_col: str,
    note_col: str,
) -> pd.DataFrame:
    reference_df = pd.DataFrame(reference_rows).drop_duplicates(subset=["site"]).sort_values("site")
    reference_df = reference_df.rename(
        columns={
            "paper_category": category_col,
            "paper_note": note_col,
        }
    )
    reference_df[named_col] = True

    top_sites = site_importance.nsmallest(top_n, "rank")[["site", "rank", "mean_abs_shap"]]
    comparison = (
        pd.merge(
            top_sites.assign(in_model_top_n=True),
            reference_df,
            on="site",
            how="outer",
        )
        .merge(
            site_importance[["site", "rank", "mean_abs_shap"]],
            on="site",
            how="left",
            suffixes=("", "_all"),
        )
        .sort_values(["rank_all", "site"], na_position="last")
    )
    comparison["in_model_top_n"] = comparison["in_model_top_n"].fillna(False)
    comparison[named_col] = comparison[named_col].fillna(False)
    comparison["rank"] = comparison["rank"].fillna(comparison["rank_all"]).astype("Int64")
    comparison["mean_abs_shap"] = comparison["mean_abs_shap"].fillna(comparison["mean_abs_shap_all"])
    return comparison[
        ["site", "rank", "mean_abs_shap", "in_model_top_n", named_col, category_col, note_col]
    ].reset_index(drop=True)


def summarize_reference_ranks(site_importance: pd.DataFrame, site_set: set[int]) -> dict[str, float | int | None]:
    ref_ranks = site_importance.loc[site_importance["site"].isin(site_set), "rank"].to_numpy(dtype=float)
    other_ranks = site_importance.loc[~site_importance["site"].isin(site_set), "rank"].to_numpy(dtype=float)
    summary: dict[str, float | int | None] = {
        "n_reference_sites": int(len(ref_ranks)),
        "median_rank": float(np.median(ref_ranks)) if len(ref_ranks) else None,
        "mean_rank": float(np.mean(ref_ranks)) if len(ref_ranks) else None,
        "top10_overlap": int(np.sum(ref_ranks <= 10)) if len(ref_ranks) else 0,
        "top20_overlap": int(np.sum(ref_ranks <= 20)) if len(ref_ranks) else 0,
        "top30_overlap": int(np.sum(ref_ranks <= 30)) if len(ref_ranks) else 0,
    }
    if len(ref_ranks) and len(other_ranks):
        u_stat, p_value = mannwhitneyu(ref_ranks, other_ranks, alternative="less")
        summary["mannwhitney_u"] = float(u_stat)
        summary["mannwhitney_one_sided_p"] = float(p_value)
    else:
        summary["mannwhitney_u"] = None
        summary["mannwhitney_one_sided_p"] = None
    return summary


def build_substitution_outputs(
    feature_summary: pd.DataFrame,
    fold_feature_abs_rows: list[dict[str, float | int]],
    all_sites: list[int],
    neher_df: pd.DataFrame,
    wic_df: pd.DataFrame,
    shah_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sub_feature_summary = feature_summary[feature_summary["feature"].str.startswith("sub_")].copy()
    sub_feature_summary["site"] = sub_feature_summary["feature"].map(parse_substitution_site)
    sub_feature_summary = sub_feature_summary.dropna(subset=["site"]).copy()
    sub_feature_summary["site"] = sub_feature_summary["site"].astype(int)

    feature_abs_mean = sub_feature_summary.set_index("feature")["mean_abs_shap"] if not sub_feature_summary.empty else pd.Series(dtype=float)

    site_rows = []
    for site in all_sites:
        site_features = sub_feature_summary.loc[sub_feature_summary["site"] == site].copy()
        site_rows.append(
            {
                "feature": f"site_{site}",
                "site": site,
                "mean_abs_shap": float(site_features["mean_abs_shap"].sum()) if not site_features.empty else 0.0,
                "mean_signed_shap": float(site_features["mean_signed_shap"].sum()) if not site_features.empty else 0.0,
                "n_substitution_features": int(len(site_features)),
                "top_substitution_terms": (
                    top_terms_from_feature_means(feature_abs_mean, site_features["feature"].tolist())
                    if not site_features.empty
                    else ""
                ),
            }
        )

    site_importance = pd.DataFrame(site_rows).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    site_importance["rank"] = np.arange(1, len(site_importance) + 1)
    site_importance = (
        site_importance.merge(neher_df, on="site", how="left")
        .merge(wic_df, on="site", how="left")
        .merge(shah_df, on="site", how="left")
        .fillna({"named_in_neher2016": False, "named_in_wic2023": False, "named_in_shah2024": False})
    )
    site_importance["named_in_any_reference"] = (
        site_importance["named_in_neher2016"] | site_importance["named_in_wic2023"] | site_importance["named_in_shah2024"]
    )
    site_importance["in_top_30"] = site_importance["rank"] <= 30

    fold_feature_abs_df = pd.DataFrame(fold_feature_abs_rows)
    folds = sorted(set(int(row["fold"]) for row in fold_feature_abs_rows))
    if fold_feature_abs_df.empty:
        fold_site_df = pd.DataFrame(columns=["fold", "site", "fold_mean_abs_shap"])
    else:
        fold_feature_abs_df = fold_feature_abs_df[fold_feature_abs_df["feature"].str.startswith("sub_")].copy()
        fold_feature_abs_df["site"] = fold_feature_abs_df["feature"].map(parse_substitution_site)
        fold_feature_abs_df = fold_feature_abs_df.dropna(subset=["site"]).copy()
        fold_feature_abs_df["site"] = fold_feature_abs_df["site"].astype(int)
        fold_site_df = (
            fold_feature_abs_df.groupby(["fold", "site"], as_index=False)["fold_mean_abs_shap"].sum()
        )

    if folds:
        full_index = pd.MultiIndex.from_product([folds, all_sites], names=["fold", "site"])
        fold_site_df = (
            fold_site_df.set_index(["fold", "site"])
            .reindex(full_index, fill_value=0.0)
            .reset_index()
        )
        fold_site_df = fold_site_df.sort_values(
            ["fold", "fold_mean_abs_shap"],
            ascending=[True, False],
        )
        fold_site_df["fold_rank"] = fold_site_df.groupby("fold")["fold_mean_abs_shap"].rank(
            method="first",
            ascending=False,
        )
    else:
        fold_site_df["fold_rank"] = []

    site_stability = (
        fold_site_df.groupby("site")
        .agg(
            mean_rank=("fold_rank", "mean"),
            sd_rank=("fold_rank", "std"),
            min_rank=("fold_rank", "min"),
            max_rank=("fold_rank", "max"),
            mean_fold_abs_shap=("fold_mean_abs_shap", "mean"),
            sd_fold_abs_shap=("fold_mean_abs_shap", "std"),
            top10_frequency=("fold_rank", lambda s: float(np.mean(np.asarray(s) <= 10))),
            top20_frequency=("fold_rank", lambda s: float(np.mean(np.asarray(s) <= 20))),
            top30_frequency=("fold_rank", lambda s: float(np.mean(np.asarray(s) <= 30))),
        )
        .reset_index()
    )
    site_stability["sd_rank"] = site_stability["sd_rank"].fillna(0.0)
    site_stability["sd_fold_abs_shap"] = site_stability["sd_fold_abs_shap"].fillna(0.0)
    site_stability = site_stability.merge(
        site_importance[
            [
                "site",
                "rank",
                "mean_abs_shap",
                "mean_signed_shap",
                "n_substitution_features",
                "top_substitution_terms",
                "named_in_neher2016",
                "named_in_wic2023",
                "named_in_shah2024",
                "named_in_any_reference",
                "neher2016_category",
                "neher2016_note",
                "wic2023_category",
                "wic2023_note",
                "shah2024_category",
                "shah2024_note",
                "in_top_30",
            ]
        ],
        on="site",
        how="left",
    ).sort_values("rank")

    return site_importance, site_stability


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    feature_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_feature_matrix.tsv"
    substitution_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_feature_matrix.tsv"
    df = pd.read_csv(feature_path, sep="\t", low_memory=False)
    df_substitution = pd.read_csv(substitution_path, sep="\t", low_memory=False)

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(df)
    context_cols = get_context_cols(cfg)
    site_state_cols = context_cols + change_cols + virus_aa_cols + serum_aa_cols
    sub_cols = sorted([c for c in df_substitution.columns if c.startswith("sub_")])
    substitution_cols = context_cols + sub_cols
    groups = df["virusStrain"].astype(str)

    n_splits = args.n_splits or int(cfg.get("n_splits", 5))
    unique_groups = groups.nunique()
    if n_splits > unique_groups:
        n_splits = unique_groups

    gkf = GroupKFold(n_splits=n_splits)
    params = get_lightgbm_params(int(cfg.get("random_state", 42)))
    rng_seed = int(cfg.get("random_state", 42))

    feature_abs_sum = pd.Series(0.0, index=site_state_cols, dtype=float)
    feature_signed_sum = pd.Series(0.0, index=site_state_cols, dtype=float)
    fold_feature_abs_mean_rows: list[dict[str, float | int]] = []
    sample_rows = []
    total_sampled_rows = 0
    best_iterations = []

    sub_feature_abs_sum = pd.Series(0.0, index=substitution_cols, dtype=float)
    sub_feature_signed_sum = pd.Series(0.0, index=substitution_cols, dtype=float)
    sub_fold_feature_abs_mean_rows: list[dict[str, float | int]] = []
    sub_best_iterations = []
    total_sub_sampled_rows = 0
    reuse_stats = new_reuse_stats()

    for fold, (train_idx, test_idx) in enumerate(gkf.split(df, groups=groups), start=1):
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()
        train_sub = df_substitution.iloc[train_idx].copy()
        test_sub = df_substitution.iloc[test_idx].copy()

        train_base, test_base = fit_split_additive(
            train_df=train_df,
            val_df=test_df,
            target_col="standardized_titer",
        )

        X_train = prepare_frame(train_df, site_state_cols)
        X_test = prepare_frame(test_df, site_state_cols)
        y_train = train_df["standardized_titer"].to_numpy(dtype=float) - train_base
        selection = inner_selection(
            train_df, "virusStrain",
            lambda a, b: fit_split_additive(a, b, "standardized_titer"),
            random_state=params["random_state"],
        )
        categorical_cols = [c for c in X_train.columns if str(X_train[c].dtype) == "category"]

        model = fit_selected_model(X_train, y_train, params, selection, categorical_cols)
        best_iterations.append(int(model.n_estimators_))

        if args.max_rows_per_fold and len(test_idx) > args.max_rows_per_fold:
            sample_local = np.sort(
                np.random.default_rng(rng_seed + fold).choice(
                    len(test_idx),
                    size=args.max_rows_per_fold,
                    replace=False,
                )
            )
        else:
            sample_local = np.arange(len(test_idx))

        X_shap = X_test.iloc[sample_local].copy()
        shap_values = shap.TreeExplainer(model).shap_values(X_shap)
        if isinstance(shap_values, list):
            shap_values = shap_values[0]
        shap_block = pd.DataFrame(shap_values, columns=site_state_cols)
        collect_site_state_reuse_stats(
            stats=reuse_stats,
            feature_values=X_shap.reset_index(drop=True),
            shap_values=shap_block.reset_index(drop=True),
            feature_cols=site_state_cols,
        )

        feature_abs_sum += shap_block.abs().sum(axis=0)
        feature_signed_sum += shap_block.sum(axis=0)
        total_sampled_rows += len(shap_block)

        fold_abs_mean = shap_block.abs().mean(axis=0)
        fold_feature_abs_mean_rows.extend(
            {"fold": fold, "feature": feature, "fold_mean_abs_shap": float(value)}
            for feature, value in fold_abs_mean.items()
        )

        if sub_cols:
            X_sub_train = prepare_frame(train_sub, substitution_cols)
            X_sub_test = prepare_frame(test_sub, substitution_cols)
            sub_categorical_cols = [c for c in X_sub_train.columns if str(X_sub_train[c].dtype) == "category"]
            model_sub = fit_selected_model(
                X_sub_train, y_train, params, selection, sub_categorical_cols,
            )
            sub_best_iterations.append(int(model_sub.n_estimators_))

            X_sub_shap = X_sub_test.iloc[sample_local].copy()
            sub_shap_values = shap.TreeExplainer(model_sub).shap_values(X_sub_shap)
            if isinstance(sub_shap_values, list):
                sub_shap_values = sub_shap_values[0]
            sub_shap_block = pd.DataFrame(sub_shap_values, columns=substitution_cols)
            collect_substitution_reuse_stats(
                stats=reuse_stats,
                feature_values=X_sub_shap.reset_index(drop=True),
                shap_values=sub_shap_block.reset_index(drop=True),
                feature_cols=substitution_cols,
            )

            sub_feature_abs_sum += sub_shap_block.abs().sum(axis=0)
            sub_feature_signed_sum += sub_shap_block.sum(axis=0)
            total_sub_sampled_rows += len(sub_shap_block)

            sub_fold_abs_mean = sub_shap_block.abs().mean(axis=0)
            sub_fold_feature_abs_mean_rows.extend(
                {"fold": fold, "feature": feature, "fold_mean_abs_shap": float(value)}
                for feature, value in sub_fold_abs_mean.items()
            )

        sampled_meta = test_df.iloc[sample_local][
            [
                "virusStrain",
                "serumStrain",
                "virusSequenceId",
                "serumSequenceId",
                "standardized_titer",
                "context_corrected_titer",
            ]
            + [col for col in context_cols if col in test_df.columns]
        ].copy()
        sampled_meta["fold"] = fold
        sample_rows.append(sampled_meta)

    if total_sampled_rows == 0:
        raise ValueError("No sampled SHAP rows were collected.")

    feature_summary = pd.DataFrame(
        {
            "feature": site_state_cols,
            "mean_abs_shap": (feature_abs_sum / total_sampled_rows).to_numpy(dtype=float),
            "mean_signed_shap": (feature_signed_sum / total_sampled_rows).to_numpy(dtype=float),
        }
    ).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    feature_summary["rank"] = np.arange(1, len(feature_summary) + 1)

    feature_abs_mean = feature_summary.set_index("feature")["mean_abs_shap"]
    feature_signed_mean = feature_summary.set_index("feature")["mean_signed_shap"]

    site_rows = []
    for site in range(1, len(change_cols) + 1):
        change_col = f"change_{site}"
        virus_cols = [c for c in virus_aa_cols if c[1:] == str(site)]
        serum_cols = [c for c in serum_aa_cols if c[:-1] == str(site)]
        available = [change_col] + virus_cols + serum_cols
        site_rows.append(
            {
                "feature": f"site_{site}",
                "site": site,
                "mean_abs_shap": float(feature_abs_mean.loc[available].sum()),
                "mean_abs_shap_change": float(feature_abs_mean.loc[[change_col]].sum()),
                "mean_abs_shap_virus_aa": float(feature_abs_mean.loc[virus_cols].sum()) if virus_cols else 0.0,
                "mean_abs_shap_serum_aa": float(feature_abs_mean.loc[serum_cols].sum()) if serum_cols else 0.0,
                "mean_signed_shap_change": float(feature_signed_mean.loc[[change_col]].sum()),
                "mean_signed_shap_virus_aa": float(feature_signed_mean.loc[virus_cols].sum()) if virus_cols else 0.0,
                "mean_signed_shap_serum_aa": float(feature_signed_mean.loc[serum_cols].sum()) if serum_cols else 0.0,
                "top_virus_aa_terms": top_terms_from_feature_means(feature_abs_mean, virus_cols),
                "top_serum_aa_terms": top_terms_from_feature_means(feature_abs_mean, serum_cols),
            }
        )

    site_importance = pd.DataFrame(site_rows).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    site_importance["rank"] = np.arange(1, len(site_importance) + 1)

    neher_df = pd.DataFrame(NEHER2016_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "neher2016_category", "paper_note": "neher2016_note"}
    )
    neher_df["named_in_neher2016"] = True
    wic_df = pd.DataFrame(WIC2023_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "wic2023_category", "paper_note": "wic2023_note"}
    )
    wic_df["named_in_wic2023"] = True
    shah_df = pd.DataFrame(SHAH2024_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "shah2024_category", "paper_note": "shah2024_note"}
    )
    shah_df["named_in_shah2024"] = True

    site_importance = (
        site_importance.merge(neher_df, on="site", how="left")
        .merge(wic_df, on="site", how="left")
        .merge(shah_df, on="site", how="left")
        .fillna({"named_in_neher2016": False, "named_in_wic2023": False, "named_in_shah2024": False})
    )
    site_importance["named_in_any_reference"] = (
        site_importance["named_in_neher2016"] | site_importance["named_in_wic2023"] | site_importance["named_in_shah2024"]
    )
    site_importance["in_top_30"] = site_importance["rank"] <= 30

    fold_feature_abs_df = pd.DataFrame(fold_feature_abs_mean_rows)
    fold_site_rows = []
    for fold in sorted(fold_feature_abs_df["fold"].unique()):
        fold_abs_mean = fold_feature_abs_df.loc[fold_feature_abs_df["fold"] == fold].set_index("feature")[
            "fold_mean_abs_shap"
        ]
        for site in range(1, len(change_cols) + 1):
            change_col = f"change_{site}"
            virus_cols = [c for c in virus_aa_cols if c[1:] == str(site)]
            serum_cols = [c for c in serum_aa_cols if c[:-1] == str(site)]
            available = [change_col] + virus_cols + serum_cols
            fold_site_rows.append(
                {
                    "fold": int(fold),
                    "site": site,
                    "fold_mean_abs_shap": float(fold_abs_mean.loc[available].sum()),
                }
            )

    fold_site_df = pd.DataFrame(fold_site_rows).sort_values(
        ["fold", "fold_mean_abs_shap"],
        ascending=[True, False],
    )
    fold_site_df["fold_rank"] = fold_site_df.groupby("fold")["fold_mean_abs_shap"].rank(
        method="first",
        ascending=False,
    )
    site_stability = (
        fold_site_df.groupby("site")
        .agg(
            mean_rank=("fold_rank", "mean"),
            sd_rank=("fold_rank", "std"),
            min_rank=("fold_rank", "min"),
            max_rank=("fold_rank", "max"),
            mean_fold_abs_shap=("fold_mean_abs_shap", "mean"),
            sd_fold_abs_shap=("fold_mean_abs_shap", "std"),
            top10_frequency=("fold_rank", lambda s: float(np.mean(np.asarray(s) <= 10))),
            top20_frequency=("fold_rank", lambda s: float(np.mean(np.asarray(s) <= 20))),
            top30_frequency=("fold_rank", lambda s: float(np.mean(np.asarray(s) <= 30))),
        )
        .reset_index()
    )
    site_stability["sd_rank"] = site_stability["sd_rank"].fillna(0.0)
    site_stability["sd_fold_abs_shap"] = site_stability["sd_fold_abs_shap"].fillna(0.0)
    site_stability = site_stability.merge(
        site_importance[
            [
                "site",
                "rank",
                "mean_abs_shap",
                "named_in_neher2016",
                "named_in_wic2023",
                "named_in_shah2024",
                "named_in_any_reference",
                "neher2016_category",
                "neher2016_note",
                "wic2023_category",
                "wic2023_note",
                "shah2024_category",
                "shah2024_note",
                "in_top_30",
            ]
        ],
        on="site",
        how="left",
    ).sort_values("rank")

    feature_summary_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_feature_shap_importance.csv"
    feature_summary.to_csv(feature_summary_path, index=False)

    site_importance_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_site_importance.csv"
    site_importance.to_csv(site_importance_path, index=False)

    site_summary_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_site_summary.tsv"
    site_importance[
        [
            "site",
            "rank",
            "mean_abs_shap",
            "mean_abs_shap_change",
            "mean_abs_shap_virus_aa",
            "mean_abs_shap_serum_aa",
            "mean_signed_shap_change",
            "mean_signed_shap_virus_aa",
            "mean_signed_shap_serum_aa",
            "top_virus_aa_terms",
            "top_serum_aa_terms",
            "named_in_neher2016",
            "neher2016_category",
            "neher2016_note",
            "named_in_wic2023",
            "wic2023_category",
            "wic2023_note",
            "named_in_shah2024",
            "shah2024_category",
            "shah2024_note",
            "named_in_any_reference",
            "in_top_30",
        ]
    ].to_csv(site_summary_path, sep="\t", index=False)

    site_stability_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_site_stability.tsv"
    site_stability.to_csv(site_stability_path, sep="\t", index=False)

    shap_sample_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_shap_sample_rows.tsv"
    pd.concat(sample_rows, ignore_index=True).to_csv(shap_sample_path, sep="\t", index=False)

    sub_feature_summary_path = None
    sub_site_importance_path = None
    sub_site_summary_path = None
    sub_site_stability_path = None
    sub_feature_summary = None
    substitution_validation = None
    if sub_cols and total_sub_sampled_rows > 0:
        sub_feature_summary = pd.DataFrame(
            {
                "feature": substitution_cols,
                "mean_abs_shap": (sub_feature_abs_sum / total_sub_sampled_rows).to_numpy(dtype=float),
                "mean_signed_shap": (sub_feature_signed_sum / total_sub_sampled_rows).to_numpy(dtype=float),
            }
        ).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
        sub_feature_summary["rank"] = np.arange(1, len(sub_feature_summary) + 1)

        sub_site_importance, sub_site_stability = build_substitution_outputs(
            feature_summary=sub_feature_summary,
            fold_feature_abs_rows=sub_fold_feature_abs_mean_rows,
            all_sites=list(range(1, len(change_cols) + 1)),
            neher_df=neher_df,
            wic_df=wic_df,
            shah_df=shah_df,
        )

        sub_feature_summary_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_feature_shap_importance.csv"
        sub_feature_summary.to_csv(sub_feature_summary_path, index=False)
        sub_site_importance_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_site_importance.csv"
        sub_site_importance.to_csv(sub_site_importance_path, index=False)
        sub_site_summary_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_site_summary.tsv"
        sub_site_importance[
            [
                "site",
                "rank",
                "mean_abs_shap",
                "mean_signed_shap",
                "n_substitution_features",
                "top_substitution_terms",
                "named_in_neher2016",
                "neher2016_category",
                "neher2016_note",
                "named_in_wic2023",
                "wic2023_category",
                "wic2023_note",
                "named_in_shah2024",
                "shah2024_category",
                "shah2024_note",
                "named_in_any_reference",
                "in_top_30",
            ]
        ].to_csv(sub_site_summary_path, sep="\t", index=False)
        sub_site_stability_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_site_stability.tsv"
        sub_site_stability.to_csv(sub_site_stability_path, sep="\t", index=False)

        substitution_validation = {
            "n_shap_folds": int(n_splits),
            "max_rows_per_fold": int(args.max_rows_per_fold),
            "n_test_rows_for_shap_sample": int(total_sub_sampled_rows),
            "mean_best_iteration": float(np.mean(sub_best_iterations)) if sub_best_iterations else None,
            "neher2016": summarize_reference_ranks(
                site_importance=sub_site_importance,
                site_set={int(row["site"]) for row in NEHER2016_H3_SITE_ROWS},
            ),
            "wic2023": summarize_reference_ranks(
                site_importance=sub_site_importance,
                site_set={int(row["site"]) for row in WIC2023_H3_SITE_ROWS},
            ),
            "shah2024": summarize_reference_ranks(
                site_importance=sub_site_importance,
                site_set={int(row["site"]) for row in SHAH2024_H3_SITE_ROWS},
            ),
        }

    reusable_shap_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_reusable_shap_values.tsv"
    reuse_stats_to_table(reuse_stats).to_csv(reusable_shap_path, sep="\t", index=False)

    validation = {
        "n_shap_folds": int(n_splits),
        "max_rows_per_fold": int(args.max_rows_per_fold),
        "n_test_rows_for_shap_sample": int(total_sampled_rows),
        "mean_best_iteration": float(np.mean(best_iterations)) if best_iterations else None,
        "neher2016": summarize_reference_ranks(
            site_importance=site_importance,
            site_set={int(row["site"]) for row in NEHER2016_H3_SITE_ROWS},
        ),
        "wic2023": summarize_reference_ranks(
            site_importance=site_importance,
            site_set={int(row["site"]) for row in WIC2023_H3_SITE_ROWS},
        ),
        "shah2024": summarize_reference_ranks(
            site_importance=site_importance,
            site_set={int(row["site"]) for row in SHAH2024_H3_SITE_ROWS},
        ),
    }
    if substitution_validation is not None:
        validation["substitution_model"] = substitution_validation

    validation_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_paper_site_validation.json"
    with open(validation_path, "w", encoding="utf-8") as f:
        json.dump(validation, f, indent=2)

    lines = [
        f"Sampled SHAP rows across folds: {total_sampled_rows}",
        f"Max sampled rows per fold: {args.max_rows_per_fold}",
        f"Distance feature columns: {get_distance_feature_cols(cfg)}",
        f"Site-state SHAP feature count: {len(site_state_cols)}",
        f"Saved SHAP sample rows: {shap_sample_path}",
        f"Saved feature SHAP importance: {feature_summary_path}",
        f"Saved site importance: {site_importance_path}",
        f"Saved site summary: {site_summary_path}",
        f"Saved site stability: {site_stability_path}",
        f"Saved reusable SHAP value table: {reusable_shap_path}",
    ]
    if sub_feature_summary_path is not None:
        lines.extend(
            [
                f"Saved substitution feature SHAP importance: {sub_feature_summary_path}",
                f"Saved substitution site importance: {sub_site_importance_path}",
                f"Saved substitution site summary: {sub_site_summary_path}",
                f"Saved substitution site stability: {sub_site_stability_path}",
            ]
        )
    lines.append(f"Saved paper-site validation: {validation_path}")
    append_qc_log(cfg, "29_shap_analysis_wic", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
