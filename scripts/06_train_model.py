#!/usr/bin/env python3
"""Train direct-comparison models against paper-style standardized HI titers."""

from __future__ import annotations

import argparse
import json
import re

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, train_test_split

from model_validation import inner_selection, fit_selected_model

from common import (
    append_qc_log,
    apply_additive_effects,
    fit_additive_effects,
    get_distance_feature_cols,
    get_lightgbm_params,
    load_config,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--n-splits", type=int, default=None)
    return parser.parse_args()


def compute_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    mae = float(mean_absolute_error(y_true, y_pred))
    r2 = float(r2_score(y_true, y_pred))
    rho = spearmanr(y_true, y_pred).correlation
    rho = float(rho) if np.isfinite(rho) else float("nan")
    return {"rmse": rmse, "mae": mae, "r2": r2, "spearman": rho}


def fit_predict(
    X_train: pd.DataFrame,
    y_train: np.ndarray,
    X_val: pd.DataFrame,
    selection: tuple,
    params: dict,
    categorical_cols: list[str],
) -> tuple[lgb.LGBMRegressor, np.ndarray]:
    model = fit_selected_model(X_train, y_train, params, selection, categorical_cols)
    preds = model.predict(X_val)
    return model, preds


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
    serum_col: str,
    virus_col: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    additive = fit_additive_effects(
        df=train_df,
        target_col=target_col,
        virus_col=virus_col,
        serum_col=serum_col,
    )
    train_base = apply_additive_effects(
        train_df,
        serum_col=serum_col,
        virus_col=virus_col,
        intercept=float(additive["intercept"]),
        serum_effects=additive["serum_effects"],
        virus_effects=additive["virus_effects"],
    )
    val_base = apply_additive_effects(
        val_df,
        serum_col=serum_col,
        virus_col=virus_col,
        intercept=float(additive["intercept"]),
        serum_effects=additive["serum_effects"],
        virus_effects=additive["virus_effects"],
    )
    return train_base, val_base, additive


def run_models_for_split(
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    df_binary: pd.DataFrame,
    df_substitution: pd.DataFrame,
    params: dict,
    rng: np.random.Generator,
    distance_cols: list[str],
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    target_col = "standardized_titer"
    serum_col = "serumIsolate"
    virus_col = "virus_strain_matched"

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(df_binary)
    sub_cols = sorted([c for c in df_substitution.columns if c.startswith("sub_")])

    site_state_cols = change_cols + virus_aa_cols + serum_aa_cols + distance_cols
    binary_cols = change_cols + distance_cols
    substitution_cols = sub_cols + distance_cols
    temporal_cols = ["temporal_distance"]

    train_bin = df_binary.iloc[train_idx].copy()
    val_bin = df_binary.iloc[val_idx].copy()
    train_sub = df_substitution.iloc[train_idx].copy()
    val_sub = df_substitution.iloc[val_idx].copy()

    y_train = train_bin[target_col].to_numpy(dtype=float)
    y_val = val_bin[target_col].to_numpy(dtype=float)

    train_base, val_base, additive = fit_split_additive(
        train_df=train_bin,
        val_df=val_bin,
        target_col=target_col,
        serum_col=serum_col,
        virus_col=virus_col,
    )
    y_train_corr = y_train - train_base
    selection = inner_selection(
        train_bin, "virus_strain_matched",
        lambda a, b: fit_split_additive(a, b, target_col, serum_col, virus_col),
        random_state=params["random_state"],
    )

    results = []
    pred_meta_cols = [
        "virus_strain_matched",
        "serum_strain_matched",
        "virusYear",
        "serumYear",
        "standardized_titer",
        "corrected_titer",
    ] + [col for col in distance_cols if col in val_bin.columns]
    pred_frame = val_bin[pred_meta_cols].copy()
    pred_frame["additive_only_pred"] = val_base

    additive_metrics = compute_metrics(y_val, val_base)
    results.append({"model_type": "additive_only", **additive_metrics, "best_iteration": 0})

    X_state_train = prepare_frame(train_bin, site_state_cols)
    X_state_val = prepare_frame(val_bin, site_state_cols)
    state_cats = [c for c in X_state_train.columns if str(X_state_train[c].dtype) == "category"]
    model_state, pred_state_corr = fit_predict(
        X_train=X_state_train,
        y_train=y_train_corr,
        X_val=X_state_val,
        selection=selection,
        params=params,
        categorical_cols=state_cats,
    )
    pred_state = val_base + pred_state_corr
    pred_frame["predicted_standardized_titer_site_state"] = pred_state
    pred_frame["predicted_corrected_titer_site_state"] = pred_state_corr
    results.append(
        {
            "model_type": "site_state",
            **compute_metrics(y_val, pred_state),
            "best_iteration": int(model_state.n_estimators_),
        }
    )

    X_bin_train = prepare_frame(train_bin, binary_cols)
    X_bin_val = prepare_frame(val_bin, binary_cols)
    bin_cats = [c for c in X_bin_train.columns if str(X_bin_train[c].dtype) == "category"]
    model_bin, pred_bin_corr = fit_predict(
        X_train=X_bin_train,
        y_train=y_train_corr,
        X_val=X_bin_val,
        selection=selection,
        params=params,
        categorical_cols=bin_cats,
    )
    pred_bin = val_base + pred_bin_corr
    pred_frame["predicted_standardized_titer_binary"] = pred_bin
    pred_frame["predicted_corrected_titer_binary"] = pred_bin_corr
    results.append(
        {
            "model_type": "binary_site",
            **compute_metrics(y_val, pred_bin),
            "best_iteration": int(model_bin.n_estimators_),
        }
    )

    X_temp_train = prepare_frame(train_bin, temporal_cols)
    X_temp_val = prepare_frame(val_bin, temporal_cols)
    model_temp, pred_temp_corr = fit_predict(
        X_train=X_temp_train,
        y_train=y_train_corr,
        X_val=X_temp_val,
        selection=selection,
        params=params,
        categorical_cols=[],
    )
    pred_temp = val_base + pred_temp_corr
    pred_frame["predicted_standardized_titer_temporal"] = pred_temp
    results.append(
        {
            "model_type": "temporal_only",
            **compute_metrics(y_val, pred_temp),
            "best_iteration": int(model_temp.n_estimators_),
        }
    )

    X_rand_train = X_bin_train.copy()
    X_rand_val = X_bin_val.copy()
    for col in change_cols:
        X_rand_train[col] = rng.permutation(X_rand_train[col].to_numpy())
        X_rand_val[col] = rng.permutation(X_rand_val[col].to_numpy())
    model_rand, pred_rand_corr = fit_predict(
        X_train=X_rand_train,
        y_train=y_train_corr,
        X_val=X_rand_val,
        selection=selection,
        params=params,
        categorical_cols=bin_cats,
    )
    pred_rand = val_base + pred_rand_corr
    pred_frame["predicted_standardized_titer_randomized_binary"] = pred_rand
    results.append(
        {
            "model_type": "randomized_binary",
            **compute_metrics(y_val, pred_rand),
            "best_iteration": int(model_rand.n_estimators_),
        }
    )

    if sub_cols:
        X_sub_train = prepare_frame(train_sub, substitution_cols)
        X_sub_val = prepare_frame(val_sub, substitution_cols)
        sub_cats = [c for c in X_sub_train.columns if str(X_sub_train[c].dtype) == "category"]
        model_sub, pred_sub_corr = fit_predict(
            X_train=X_sub_train,
            y_train=y_train_corr,
            X_val=X_sub_val,
            selection=selection,
            params=params,
            categorical_cols=sub_cats,
        )
        pred_sub = val_base + pred_sub_corr
        pred_frame["predicted_standardized_titer_substitution"] = pred_sub
        pred_frame["predicted_corrected_titer_substitution"] = pred_sub_corr
        results.append(
            {
                "model_type": "substitution_identity",
                **compute_metrics(y_val, pred_sub),
                "best_iteration": int(model_sub.n_estimators_),
            }
        )
    else:
        pred_frame["predicted_standardized_titer_substitution"] = np.nan
        pred_frame["predicted_corrected_titer_substitution"] = np.nan

    pred_frame["train_additive_intercept"] = float(additive["intercept"])
    return results, pred_frame


def run_neher_holdouts(
    df_binary: pd.DataFrame,
    df_substitution: pd.DataFrame,
    params: dict,
    rng_seed: int,
    distance_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    idx_all = np.arange(len(df_binary))
    groups = df_binary["virus_strain_matched"].astype(str)

    split_defs = []
    rand_train, rand_test = train_test_split(
        idx_all, test_size=0.10, random_state=rng_seed, shuffle=True
    )
    split_defs.append(("random_10pct_measurements", rand_train, rand_test))

    gss = GroupShuffleSplit(n_splits=1, test_size=0.10, random_state=rng_seed)
    grp_train, grp_test = next(gss.split(df_binary, groups=groups, y=None))
    split_defs.append(("grouped_10pct_viruses", grp_train, grp_test))

    pred_parts = []
    metric_rows = []
    for split_name, train_idx, test_idx in split_defs:
        split_results, split_pred = run_models_for_split(
            train_idx=train_idx,
            val_idx=test_idx,
            df_binary=df_binary,
            df_substitution=df_substitution,
            params=params,
            rng=np.random.default_rng(rng_seed),
            distance_cols=distance_cols,
        )
        split_pred = split_pred.reset_index(drop=True)
        split_pred["split_type"] = split_name
        pred_parts.append(split_pred)

        for row in split_results:
            metric_rows.append(
                {
                    "split_type": split_name,
                    "model_type": row["model_type"],
                    "rmse": row["rmse"],
                    "mae": row["mae"],
                    "r2": row["r2"],
                    "spearman": row["spearman"],
                    "best_iteration": row["best_iteration"],
                    "n_test": len(test_idx),
                    "n_test_unique_viruses": int(df_binary.iloc[test_idx]["virus_strain_matched"].nunique()),
                }
            )

    return pd.concat(pred_parts, ignore_index=True), pd.DataFrame(metric_rows)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    binary_path = f"{cfg['output_dir']}/{cfg['subtype']}_feature_matrix.csv"
    substitution_path = f"{cfg['output_dir']}/{cfg['subtype']}_substitution_feature_matrix.csv"
    df_binary = pd.read_csv(binary_path, low_memory=False)
    df_substitution = pd.read_csv(substitution_path, low_memory=False)

    groups = df_binary["virus_strain_matched"].astype(str)
    n_splits = args.n_splits or int(cfg.get("n_splits", 5))
    unique_groups = groups.nunique()
    if unique_groups < 2:
        raise ValueError("Need at least 2 unique virus groups for grouped CV.")
    if n_splits > unique_groups:
        n_splits = unique_groups

    gkf = GroupKFold(n_splits=n_splits)
    params = get_lightgbm_params(int(cfg.get("random_state", 42)))
    rng_seed = int(cfg.get("random_state", 42))
    distance_cols = get_distance_feature_cols(cfg)

    results = []
    oof_parts = []
    site_state_best_iterations = []
    binary_best_iterations = []
    substitution_best_iterations = []

    for fold, (train_idx, val_idx) in enumerate(gkf.split(df_binary, groups=groups), start=1):
        split_results, pred_frame = run_models_for_split(
            train_idx=train_idx,
            val_idx=val_idx,
            df_binary=df_binary,
            df_substitution=df_substitution,
            params=params,
            rng=np.random.default_rng(rng_seed + fold),
            distance_cols=distance_cols,
        )
        for row in split_results:
            row["fold"] = fold
            row["n_train"] = len(train_idx)
            row["n_val"] = len(val_idx)
            results.append(row)
            if row["model_type"] == "site_state" and row["best_iteration"] > 0:
                site_state_best_iterations.append(int(row["best_iteration"]))
            if row["model_type"] == "binary_site" and row["best_iteration"] > 0:
                binary_best_iterations.append(int(row["best_iteration"]))
            if row["model_type"] == "substitution_identity" and row["best_iteration"] > 0:
                substitution_best_iterations.append(int(row["best_iteration"]))

        pred_frame["fold"] = fold
        oof_parts.append(pred_frame)

    results_df = pd.DataFrame(results)
    cv_path = f"{cfg['output_dir']}/{cfg['subtype']}_cv_results.csv"
    results_df.to_csv(cv_path, index=False)

    oof_df = pd.concat(oof_parts, ignore_index=True)
    oof_path = f"{cfg['output_dir']}/{cfg['subtype']}_oof_predictions.csv"
    oof_df.to_csv(oof_path, index=False)

    summary = (
        results_df.groupby("model_type")[["rmse", "mae", "r2", "spearman"]]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary_path = f"{cfg['output_dir']}/{cfg['subtype']}_cv_summary.csv"
    summary.to_csv(summary_path, index=False)

    full_additive = fit_additive_effects(
        df=df_binary,
        target_col="standardized_titer",
        virus_col="virus_strain_matched",
        serum_col="serumIsolate",
    )
    full_base = apply_additive_effects(
        df_binary,
        serum_col="serumIsolate",
        virus_col="virus_strain_matched",
        intercept=float(full_additive["intercept"]),
        serum_effects=full_additive["serum_effects"],
        virus_effects=full_additive["virus_effects"],
    )
    y_full_corr = df_binary["standardized_titer"].to_numpy(dtype=float) - full_base

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(df_binary)
    sub_cols = sorted([c for c in df_substitution.columns if c.startswith("sub_")])
    site_state_cols = change_cols + virus_aa_cols + serum_aa_cols + distance_cols
    binary_cols = change_cols + distance_cols
    substitution_cols = sub_cols + distance_cols

    state_params = params.copy()
    if site_state_best_iterations:
        state_params["n_estimators"] = max(1, int(np.round(np.mean(site_state_best_iterations))))
    final_state = lgb.LGBMRegressor(**state_params)
    X_state_full = prepare_frame(df_binary, site_state_cols)
    state_cats_full = [c for c in X_state_full.columns if str(X_state_full[c].dtype) == "category"]
    final_state.fit(X_state_full, y_full_corr, categorical_feature=state_cats_full)
    state_model_path = f"{cfg['output_dir']}/{cfg['subtype']}_site_state_model.txt"
    final_state.booster_.save_model(state_model_path)

    binary_params = params.copy()
    if binary_best_iterations:
        binary_params["n_estimators"] = max(1, int(np.round(np.mean(binary_best_iterations))))
    final_binary = lgb.LGBMRegressor(**binary_params)
    X_binary_full = prepare_frame(df_binary, binary_cols)
    binary_cats_full = [c for c in X_binary_full.columns if str(X_binary_full[c].dtype) == "category"]
    final_binary.fit(X_binary_full, y_full_corr, categorical_feature=binary_cats_full)
    binary_model_path = f"{cfg['output_dir']}/{cfg['subtype']}_binary_model.txt"
    final_binary.booster_.save_model(binary_model_path)

    substitution_model_path = ""
    if substitution_cols:
        substitution_params = params.copy()
        if substitution_best_iterations:
            substitution_params["n_estimators"] = max(
                1, int(np.round(np.mean(substitution_best_iterations)))
            )
        final_substitution = lgb.LGBMRegressor(**substitution_params)
        X_sub_full = prepare_frame(df_substitution, substitution_cols)
        sub_cats_full = [c for c in X_sub_full.columns if str(X_sub_full[c].dtype) == "category"]
        final_substitution.fit(
            X_sub_full,
            y_full_corr,
            categorical_feature=sub_cats_full,
        )
        substitution_model_path = f"{cfg['output_dir']}/{cfg['subtype']}_substitution_model.txt"
        final_substitution.booster_.save_model(substitution_model_path)

    additive_path = f"{cfg['output_dir']}/{cfg['subtype']}_additive_effects.json"
    with open(additive_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "intercept": float(full_additive["intercept"]),
                "n_sera": int(full_additive["n_sera"]),
                "n_viruses": int(full_additive["n_viruses"]),
                "serum_effects": full_additive["serum_effects"],
                "virus_effects": full_additive["virus_effects"],
            },
            f,
            indent=2,
        )

    neher_pred, neher_metrics = run_neher_holdouts(
        df_binary=df_binary,
        df_substitution=df_substitution,
        params=params,
        rng_seed=rng_seed,
        distance_cols=distance_cols,
    )
    neher_pred_path = f"{cfg['output_dir']}/{cfg['subtype']}_neher_analog_predictions.csv"
    neher_pred.to_csv(neher_pred_path, index=False)
    neher_metrics_path = f"{cfg['output_dir']}/{cfg['subtype']}_neher_analog_metrics.csv"
    neher_metrics.to_csv(neher_metrics_path, index=False)

    lines = [
        f"Feature rows: {len(df_binary)}",
        f"Distance feature columns: {distance_cols}",
        f"Site-state feature columns: {len(site_state_cols)}",
        f"Binary feature columns: {len(binary_cols)}",
        f"Substitution feature columns: {len(sub_cols)}",
        f"Grouped CV folds: {n_splits}",
        f"Unique virus groups: {unique_groups}",
        f"Saved CV results: {cv_path}",
        f"Saved CV summary: {summary_path}",
        f"Saved OOF predictions: {oof_path}",
        f"Saved additive effects: {additive_path}",
        f"Saved final site-state model: {state_model_path}",
        f"Saved final binary model: {binary_model_path}",
        f"Saved final substitution model: {substitution_model_path or 'none'}",
        f"Saved Neher-analog predictions: {neher_pred_path}",
        f"Saved Neher-analog metrics: {neher_metrics_path}",
        "Site-state CV metrics (mean +/- std): "
        + "; ".join(
            f"{metric}={vals['mean']:.4f}+/-{vals['std']:.4f}"
            for metric, vals in (
                results_df[results_df["model_type"] == "site_state"][
                    ["rmse", "mae", "r2", "spearman"]
                ]
                .agg(["mean", "std"])
                .to_dict()
                .items()
            )
        ),
    ]
    append_qc_log(cfg, "06_train_model", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
