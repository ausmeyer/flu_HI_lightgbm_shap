#!/usr/bin/env python3
"""Compute cross-validated SHAP summaries for the H3N2 direct-comparison pipeline."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
import shap
from scipy.stats import mannwhitneyu
from sklearn.model_selection import GroupKFold

from common import (
    append_qc_log,
    apply_additive_effects,
    fit_additive_effects,
    get_distance_feature_cols,
    get_lightgbm_params,
    load_config,
)
from paper_sites import NEHER2016_H3_SITE_ROWS, SHAH2024_H3_SITE_ROWS, WIC2023_H3_SITE_ROWS


PAPER_SITE_ROWS = NEHER2016_H3_SITE_ROWS


def add_reference_membership(df: pd.DataFrame) -> pd.DataFrame:
    neher_df = pd.DataFrame(NEHER2016_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "neher2016_category", "paper_note": "neher2016_note"}
    )
    neher_df["named_in_neher2016"] = True
    harvey_df = pd.DataFrame(WIC2023_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "wic2023_category", "paper_note": "wic2023_note"}
    )
    harvey_df["named_in_wic2023"] = True
    shah_df = pd.DataFrame(SHAH2024_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "shah2024_category", "paper_note": "shah2024_note"}
    )
    shah_df["named_in_shah2024"] = True

    merged = (
        df.merge(neher_df, on="site", how="left")
        .merge(harvey_df, on="site", how="left")
        .merge(shah_df, on="site", how="left")
        .fillna({"named_in_neher2016": False, "named_in_wic2023": False, "named_in_shah2024": False})
    )
    merged["named_in_any_reference"] = (
        merged["named_in_neher2016"] | merged["named_in_wic2023"] | merged["named_in_shah2024"]
    )
    merged["paper_named"] = merged["named_in_neher2016"]
    merged["paper_category"] = merged["neher2016_category"]
    merged["paper_note"] = merged["neher2016_note"]
    return merged


@dataclass
class FoldOutputs:
    shap_blocks: list[pd.DataFrame]
    feature_blocks: list[pd.DataFrame]
    best_iterations: list[int]
    sub_feature_abs_sum: pd.Series
    sub_feature_signed_sum: pd.Series
    sub_fold_feature_abs_rows: list[dict[str, float | int]]
    sub_best_iterations: list[int]
    total_sub_rows: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--n-splits", type=int, default=None)
    return parser.parse_args()


def top_terms(series: pd.Series, n: int = 3) -> str:
    ranked = series.sort_values(ascending=False)
    ranked = ranked[ranked > 0]
    return ", ".join(ranked.head(n).index.tolist())


def parse_substitution_site(feature: str) -> int | None:
    match = re.fullmatch(r"sub_(\d+)_([A-Z])_([A-Z])", feature)
    if match is None:
        return None
    return int(match.group(1))


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


def prepare_frame(df: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    X = df[feature_cols].copy()
    for col in X.columns:
        if not pd.api.types.is_numeric_dtype(X[col]):
            X[col] = X[col].astype("category")
    return X


def fit_additive_baseline(train_df: pd.DataFrame, test_df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    additive = fit_additive_effects(
        df=train_df,
        target_col="standardized_titer",
        virus_col="virus_strain_matched",
        serum_col="serumIsolate",
    )
    train_base = apply_additive_effects(
        train_df,
        serum_col="serumIsolate",
        virus_col="virus_strain_matched",
        intercept=float(additive["intercept"]),
        serum_effects=additive["serum_effects"],
        virus_effects=additive["virus_effects"],
    )
    test_base = apply_additive_effects(
        test_df,
        serum_col="serumIsolate",
        virus_col="virus_strain_matched",
        intercept=float(additive["intercept"]),
        serum_effects=additive["serum_effects"],
        virus_effects=additive["virus_effects"],
    )
    return train_base, test_base


def fit_shap_model(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    params: dict,
) -> tuple[lgb.LGBMRegressor, pd.DataFrame]:
    categorical_cols = [c for c in X_train.columns if str(X_train[c].dtype) == "category"]
    model = lgb.LGBMRegressor(**params)
    model.fit(
        X_train,
        y_train,
        eval_set=[(X_test, y_test)],
        categorical_feature=categorical_cols,
        callbacks=[lgb.early_stopping(50, verbose=False)],
    )
    shap_values = shap.TreeExplainer(model).shap_values(X_test)
    if isinstance(shap_values, list):
        shap_values = shap_values[0]
    shap_df = pd.DataFrame(shap_values, columns=X_test.columns)
    return model, shap_df


def build_feature_block(
    test_df: pd.DataFrame,
    X_test: pd.DataFrame,
    fold: int,
    test_idx: np.ndarray,
    distance_cols: list[str],
) -> pd.DataFrame:
    meta_cols = [
        "virus_strain_matched",
        "serum_strain_matched",
        "standardized_titer",
        "corrected_titer",
    ] + [col for col in distance_cols if col in test_df.columns]
    feature_block = test_df[meta_cols].reset_index(drop=True)
    feature_block = pd.concat([feature_block, X_test.reset_index(drop=True)], axis=1)
    feature_block["fold"] = fold
    feature_block["row_index"] = test_idx
    return feature_block


def run_cv_shap(
    df: pd.DataFrame,
    df_substitution: pd.DataFrame,
    main_cols: list[str],
    substitution_cols: list[str],
    params: dict,
    n_splits: int,
    distance_cols: list[str],
) -> FoldOutputs:
    groups = df["virus_strain_matched"].astype(str)
    gkf = GroupKFold(n_splits=n_splits)

    shap_blocks: list[pd.DataFrame] = []
    feature_blocks: list[pd.DataFrame] = []
    best_iterations: list[int] = []

    sub_feature_abs_sum = pd.Series(0.0, index=substitution_cols, dtype=float)
    sub_feature_signed_sum = pd.Series(0.0, index=substitution_cols, dtype=float)
    sub_fold_feature_abs_rows: list[dict[str, float | int]] = []
    sub_best_iterations: list[int] = []
    total_sub_rows = 0

    for fold, (train_idx, test_idx) in enumerate(gkf.split(df, groups=groups), start=1):
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()
        train_sub = df_substitution.iloc[train_idx].copy()
        test_sub = df_substitution.iloc[test_idx].copy()

        train_base, test_base = fit_additive_baseline(train_df=train_df, test_df=test_df)
        y_train = train_df["standardized_titer"].to_numpy(dtype=float) - train_base
        y_test = test_df["standardized_titer"].to_numpy(dtype=float) - test_base

        X_train = prepare_frame(train_df, main_cols)
        X_test = prepare_frame(test_df, main_cols)
        model, shap_df = fit_shap_model(X_train=X_train, y_train=y_train, X_test=X_test, y_test=y_test, params=params)
        best_iterations.append(int(model.best_iteration_ or params["n_estimators"]))

        shap_df["fold"] = fold
        shap_df["row_index"] = test_idx
        shap_blocks.append(shap_df)
        feature_blocks.append(
            build_feature_block(
                test_df=test_df,
                X_test=X_test,
                fold=fold,
                test_idx=test_idx,
                distance_cols=distance_cols,
            )
        )

        if not substitution_cols:
            continue

        X_sub_train = prepare_frame(train_sub, substitution_cols)
        X_sub_test = prepare_frame(test_sub, substitution_cols)
        sub_model, sub_shap_df = fit_shap_model(
            X_train=X_sub_train,
            y_train=y_train,
            X_test=X_sub_test,
            y_test=y_test,
            params=params,
        )
        sub_best_iterations.append(int(sub_model.best_iteration_ or params["n_estimators"]))
        sub_feature_abs_sum += sub_shap_df.abs().sum(axis=0)
        sub_feature_signed_sum += sub_shap_df.sum(axis=0)
        total_sub_rows += len(sub_shap_df)

        fold_abs_mean = sub_shap_df.abs().mean(axis=0)
        sub_fold_feature_abs_rows.extend(
            {"fold": fold, "feature": feature, "fold_mean_abs_shap": float(value)}
            for feature, value in fold_abs_mean.items()
        )

    return FoldOutputs(
        shap_blocks=shap_blocks,
        feature_blocks=feature_blocks,
        best_iterations=best_iterations,
        sub_feature_abs_sum=sub_feature_abs_sum,
        sub_feature_signed_sum=sub_feature_signed_sum,
        sub_fold_feature_abs_rows=sub_fold_feature_abs_rows,
        sub_best_iterations=sub_best_iterations,
        total_sub_rows=total_sub_rows,
    )


def build_feature_summary(shap_df: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    summary = pd.DataFrame(
        {
            "feature": feature_cols,
            "mean_abs_shap": shap_df[feature_cols].abs().mean(axis=0).to_numpy(dtype=float),
            "mean_signed_shap": shap_df[feature_cols].mean(axis=0).to_numpy(dtype=float),
        }
    ).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    summary["rank"] = np.arange(1, len(summary) + 1)
    return summary


def build_site_state_outputs(
    shap_df: pd.DataFrame,
    change_cols: list[str],
    virus_aa_cols: list[str],
    serum_aa_cols: list[str],
    paper_sites: pd.DataFrame,
    koel_sites: set[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    site_rows = []
    fold_site_rows = []
    for site in range(1, len(change_cols) + 1):
        change_col = f"change_{site}"
        virus_cols = [c for c in virus_aa_cols if c[1:] == str(site)]
        serum_cols = [c for c in serum_aa_cols if c[:-1] == str(site)]
        available = [change_col] + virus_cols + serum_cols
        site_rows.append(
            {
                "feature": f"site_{site}",
                "site": site,
                "mean_abs_shap": float(shap_df[available].abs().sum(axis=1).mean()),
                "mean_abs_shap_change": float(shap_df[[change_col]].abs().sum(axis=1).mean()),
                "mean_abs_shap_virus_aa": float(shap_df[virus_cols].abs().sum(axis=1).mean()) if virus_cols else 0.0,
                "mean_abs_shap_serum_aa": float(shap_df[serum_cols].abs().sum(axis=1).mean()) if serum_cols else 0.0,
                "mean_signed_shap_change": float(shap_df[[change_col]].sum(axis=1).mean()),
                "mean_signed_shap_virus_aa": float(shap_df[virus_cols].sum(axis=1).mean()) if virus_cols else 0.0,
                "mean_signed_shap_serum_aa": float(shap_df[serum_cols].sum(axis=1).mean()) if serum_cols else 0.0,
                "top_virus_aa_terms": top_terms(shap_df[virus_cols].abs().mean(axis=0), n=3) if virus_cols else "",
                "top_serum_aa_terms": top_terms(shap_df[serum_cols].abs().mean(axis=0), n=3) if serum_cols else "",
            }
        )
        for fold in sorted(shap_df["fold"].unique()):
            fold_df = shap_df.loc[shap_df["fold"] == fold, available]
            fold_site_rows.append(
                {
                    "fold": int(fold),
                    "site": site,
                    "fold_mean_abs_shap": float(fold_df.abs().sum(axis=1).mean()),
                }
            )

    site_importance = pd.DataFrame(site_rows).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    site_importance["rank"] = np.arange(1, len(site_importance) + 1)
    site_importance["is_koel_site"] = site_importance["site"].isin(koel_sites)
    site_importance["in_top_30"] = site_importance["rank"] <= 30
    site_importance = add_reference_membership(site_importance.merge(paper_sites, on="site", how="left"))

    fold_site_df = pd.DataFrame(fold_site_rows).sort_values(
        ["fold", "fold_mean_abs_shap"],
        ascending=[True, False],
    )
    fold_site_df["fold_rank"] = fold_site_df.groupby("fold")["fold_mean_abs_shap"].rank(
        method="first",
        ascending=False,
    )
    site_stability = build_site_stability(
        fold_site_df=fold_site_df,
        site_importance=site_importance,
        merge_cols=[
            "site",
            "rank",
            "mean_abs_shap",
            "paper_named",
            "paper_category",
            "paper_note",
            "named_in_neher2016",
            "named_in_wic2023",
            "named_in_shah2024",
            "named_in_any_reference",
            "is_koel_site",
            "in_top_30",
        ],
    )
    return site_importance, site_stability


def build_site_stability(
    fold_site_df: pd.DataFrame,
    site_importance: pd.DataFrame,
    merge_cols: list[str],
) -> pd.DataFrame:
    stability_df = (
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
    stability_df["sd_rank"] = stability_df["sd_rank"].fillna(0.0)
    stability_df["sd_fold_abs_shap"] = stability_df["sd_fold_abs_shap"].fillna(0.0)
    return stability_df.merge(site_importance[merge_cols], on="site", how="left").sort_values("rank")


def build_substitution_outputs(
    feature_summary: pd.DataFrame,
    fold_feature_abs_rows: list[dict[str, float | int]],
    all_sites: list[int],
    paper_sites: pd.DataFrame,
    koel_sites: set[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sub_feature_summary = feature_summary[feature_summary["feature"].str.startswith("sub_")].copy()
    sub_feature_summary["site"] = sub_feature_summary["feature"].map(parse_substitution_site)
    sub_feature_summary = sub_feature_summary.dropna(subset=["site"]).copy()
    sub_feature_summary["site"] = sub_feature_summary["site"].astype(int)

    feature_abs_mean = (
        sub_feature_summary.set_index("feature")["mean_abs_shap"] if not sub_feature_summary.empty else pd.Series(dtype=float)
    )
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
                    top_terms(feature_abs_mean.loc[site_features["feature"]], n=3) if not site_features.empty else ""
                ),
            }
        )

    site_importance = pd.DataFrame(site_rows).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    site_importance["rank"] = np.arange(1, len(site_importance) + 1)
    site_importance["is_koel_site"] = site_importance["site"].isin(koel_sites)
    site_importance["in_top_30"] = site_importance["rank"] <= 30
    site_importance = add_reference_membership(site_importance.merge(paper_sites, on="site", how="left"))

    fold_feature_abs_df = pd.DataFrame(fold_feature_abs_rows)
    if not fold_feature_abs_df.empty:
        fold_feature_abs_df = fold_feature_abs_df[fold_feature_abs_df["feature"].str.startswith("sub_")].copy()
        fold_feature_abs_df["site"] = fold_feature_abs_df["feature"].map(parse_substitution_site)
        fold_feature_abs_df = fold_feature_abs_df.dropna(subset=["site"]).copy()
        fold_feature_abs_df["site"] = fold_feature_abs_df["site"].astype(int)
        fold_site_df = fold_feature_abs_df.groupby(["fold", "site"], as_index=False)["fold_mean_abs_shap"].sum()
        folds = sorted(fold_site_df["fold"].unique())
        full_index = pd.MultiIndex.from_product([folds, all_sites], names=["fold", "site"])
        fold_site_df = fold_site_df.set_index(["fold", "site"]).reindex(full_index, fill_value=0.0).reset_index()
        fold_site_df = fold_site_df.sort_values(["fold", "fold_mean_abs_shap"], ascending=[True, False])
        fold_site_df["fold_rank"] = fold_site_df.groupby("fold")["fold_mean_abs_shap"].rank(
            method="first",
            ascending=False,
        )
    else:
        fold_site_df = pd.DataFrame(columns=["fold", "site", "fold_mean_abs_shap", "fold_rank"])

    site_stability = build_site_stability(
        fold_site_df=fold_site_df,
        site_importance=site_importance,
        merge_cols=[
            "site",
            "rank",
            "mean_abs_shap",
            "mean_signed_shap",
            "n_substitution_features",
            "top_substitution_terms",
            "paper_named",
            "paper_category",
            "paper_note",
            "named_in_neher2016",
            "named_in_wic2023",
            "named_in_shah2024",
            "named_in_any_reference",
            "is_koel_site",
            "in_top_30",
        ],
    )
    return site_importance, site_stability


def build_koel_validation(
    site_importance: pd.DataFrame,
    koel_sites: set[int],
    mean_best_iteration: float | None,
) -> dict[str, object]:
    koel_ranks = site_importance.loc[site_importance["is_koel_site"], "rank"].to_numpy(dtype=float)
    non_koel_ranks = site_importance.loc[~site_importance["is_koel_site"], "rank"].to_numpy(dtype=float)
    if len(koel_ranks) > 0 and len(non_koel_ranks) > 0:
        u_stat, p_value = mannwhitneyu(koel_ranks, non_koel_ranks, alternative="less")
        median_rank = float(np.median(koel_ranks))
    else:
        u_stat, p_value, median_rank = np.nan, np.nan, np.nan
    return {
        "koel_sites": sorted(koel_sites),
        "koel_sites_found": sorted(int(x) for x in site_importance.loc[site_importance["is_koel_site"], "site"].tolist()),
        "koel_median_rank": median_rank,
        "mannwhitney_u": float(u_stat) if np.isfinite(u_stat) else None,
        "mannwhitney_one_sided_p": float(p_value) if np.isfinite(p_value) else None,
        "mean_best_iteration": mean_best_iteration,
    }


def save_outputs(
    cfg: dict,
    shap_df: pd.DataFrame,
    shap_feature_df: pd.DataFrame,
    feature_summary: pd.DataFrame,
    site_importance: pd.DataFrame,
    site_stability: pd.DataFrame,
    substitution_feature_summary: pd.DataFrame | None,
    substitution_site_importance: pd.DataFrame | None,
    substitution_site_stability: pd.DataFrame | None,
    validation: dict[str, object],
) -> list[str]:
    out_dir = Path(cfg["output_dir"])
    subtype = cfg["subtype"]

    shap_path = out_dir / f"{subtype}_shap_values.csv"
    shap_df.to_csv(shap_path, index=False)

    shap_feature_path = out_dir / f"{subtype}_shap_test_features.csv"
    shap_feature_df.to_csv(shap_feature_path, index=False)

    feature_summary_path = out_dir / f"{subtype}_feature_shap_importance.csv"
    feature_summary.to_csv(feature_summary_path, index=False)

    site_importance_path = out_dir / f"{subtype}_site_importance.csv"
    site_importance.to_csv(site_importance_path, index=False)

    site_summary_path = out_dir / f"{subtype}_site_summary.tsv"
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
            "paper_named",
            "paper_category",
            "paper_note",
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
            "is_koel_site",
            "in_top_30",
        ]
    ].to_csv(site_summary_path, sep="\t", index=False)

    site_stability_path = out_dir / f"{subtype}_site_stability.tsv"
    site_stability.to_csv(site_stability_path, sep="\t", index=False)

    lines = [
        f"SHAP rows across folds: {len(shap_feature_df)}",
        f"SHAP feature count: {feature_summary['feature'].nunique()}",
        f"Saved SHAP values: {shap_path}",
        f"Saved SHAP test features: {shap_feature_path}",
        f"Saved feature SHAP importance: {feature_summary_path}",
        f"Saved site importance: {site_importance_path}",
        f"Saved site summary: {site_summary_path}",
        f"Saved site stability: {site_stability_path}",
    ]

    if substitution_feature_summary is not None and substitution_site_importance is not None and substitution_site_stability is not None:
        sub_feature_summary_path = out_dir / f"{subtype}_substitution_feature_shap_importance.csv"
        substitution_feature_summary.to_csv(sub_feature_summary_path, index=False)

        sub_site_importance_path = out_dir / f"{subtype}_substitution_site_importance.csv"
        substitution_site_importance.to_csv(sub_site_importance_path, index=False)

        sub_site_summary_path = out_dir / f"{subtype}_substitution_site_summary.tsv"
        substitution_site_importance[
            [
                "site",
                "rank",
                "mean_abs_shap",
                "mean_signed_shap",
                "n_substitution_features",
                "top_substitution_terms",
                "paper_named",
                "paper_category",
                "paper_note",
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
                "is_koel_site",
                "in_top_30",
            ]
        ].to_csv(sub_site_summary_path, sep="\t", index=False)

        sub_site_stability_path = out_dir / f"{subtype}_substitution_site_stability.tsv"
        substitution_site_stability.to_csv(sub_site_stability_path, sep="\t", index=False)

        lines.extend(
            [
                f"Saved substitution feature SHAP importance: {sub_feature_summary_path}",
                f"Saved substitution site importance: {sub_site_importance_path}",
                f"Saved substitution site summary: {sub_site_summary_path}",
                f"Saved substitution site stability: {sub_site_stability_path}",
            ]
        )

    validation_path = out_dir / f"{subtype}_koel_validation.json"
    with open(validation_path, "w", encoding="utf-8") as f:
        json.dump(validation, f, indent=2)
    lines.extend(
        [
            f"Saved Koel validation: {validation_path}",
            f"Koel median rank: {validation['koel_median_rank']}",
            f"Mann-Whitney one-sided p-value: {validation['mannwhitney_one_sided_p']}",
        ]
    )
    return lines


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    out_dir = Path(cfg["output_dir"])
    subtype = cfg["subtype"]
    df = pd.read_csv(out_dir / f"{subtype}_feature_matrix.csv")
    df_substitution = pd.read_csv(out_dir / f"{subtype}_substitution_feature_matrix.csv")

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(df)
    distance_cols = get_distance_feature_cols(cfg)
    main_cols = change_cols + virus_aa_cols + serum_aa_cols + distance_cols
    substitution_cols = sorted([c for c in df_substitution.columns if c.startswith("sub_")]) + distance_cols

    n_splits = min(args.n_splits or int(cfg.get("n_splits", 5)), df["virus_strain_matched"].nunique())
    params = get_lightgbm_params(int(cfg.get("random_state", 42)))
    fold_outputs = run_cv_shap(
        df=df,
        df_substitution=df_substitution,
        main_cols=main_cols,
        substitution_cols=substitution_cols,
        params=params,
        n_splits=n_splits,
        distance_cols=distance_cols,
    )

    shap_df = pd.concat(fold_outputs.shap_blocks, ignore_index=True)
    shap_feature_df = pd.concat(fold_outputs.feature_blocks, ignore_index=True)
    feature_summary = build_feature_summary(shap_df=shap_df, feature_cols=main_cols)

    paper_sites = pd.DataFrame(PAPER_SITE_ROWS).drop_duplicates(subset=["site"])
    koel_sites = set(int(x) for x in cfg.get("koel_sites", []))
    site_importance, site_stability = build_site_state_outputs(
        shap_df=shap_df,
        change_cols=change_cols,
        virus_aa_cols=virus_aa_cols,
        serum_aa_cols=serum_aa_cols,
        paper_sites=paper_sites,
        koel_sites=koel_sites,
    )

    substitution_feature_summary = None
    substitution_site_importance = None
    substitution_site_stability = None
    substitution_validation = None
    if substitution_cols and fold_outputs.total_sub_rows > 0:
        substitution_feature_summary = pd.DataFrame(
            {
                "feature": substitution_cols,
                "mean_abs_shap": (fold_outputs.sub_feature_abs_sum / fold_outputs.total_sub_rows).to_numpy(dtype=float),
                "mean_signed_shap": (fold_outputs.sub_feature_signed_sum / fold_outputs.total_sub_rows).to_numpy(dtype=float),
            }
        ).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
        substitution_feature_summary["rank"] = np.arange(1, len(substitution_feature_summary) + 1)
        substitution_site_importance, substitution_site_stability = build_substitution_outputs(
            feature_summary=substitution_feature_summary,
            fold_feature_abs_rows=fold_outputs.sub_fold_feature_abs_rows,
            all_sites=list(range(1, len(change_cols) + 1)),
            paper_sites=paper_sites,
            koel_sites=koel_sites,
        )
        substitution_validation = build_koel_validation(
            site_importance=substitution_site_importance,
            koel_sites=koel_sites,
            mean_best_iteration=float(np.mean(fold_outputs.sub_best_iterations)) if fold_outputs.sub_best_iterations else None,
        )

    validation = build_koel_validation(
        site_importance=site_importance,
        koel_sites=koel_sites,
        mean_best_iteration=float(np.mean(fold_outputs.best_iterations)) if fold_outputs.best_iterations else None,
    )
    validation["n_test_rows_for_shap"] = int(len(shap_feature_df))
    validation["n_shap_folds"] = int(n_splits)
    if substitution_validation is not None:
        validation["substitution_model"] = {
            **substitution_validation,
            "n_test_rows_for_shap": int(fold_outputs.total_sub_rows),
            "n_shap_folds": int(n_splits),
        }

    lines = save_outputs(
        cfg=cfg,
        shap_df=shap_df,
        shap_feature_df=shap_feature_df,
        feature_summary=feature_summary,
        site_importance=site_importance,
        site_stability=site_stability,
        substitution_feature_summary=substitution_feature_summary,
        substitution_site_importance=substitution_site_importance,
        substitution_site_stability=substitution_site_stability,
        validation=validation,
    )
    append_qc_log(cfg, "07_shap_analysis", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
