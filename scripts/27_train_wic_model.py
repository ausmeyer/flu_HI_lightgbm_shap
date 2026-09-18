#!/usr/bin/env python3
"""Train HA1-only LightGBM models for the WIC dataset."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, train_test_split

from model_validation import inner_selection, fit_selected_model

from common import (
    append_qc_log,
    apply_additive_effects_multi,
    fit_additive_effects_multi,
    get_distance_feature_cols,
    get_lightgbm_params,
    load_config,
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
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
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
    return train_base, val_base, additive


def run_models_for_split(
    train_idx: np.ndarray,
    val_idx: np.ndarray,
    df_binary: pd.DataFrame,
    df_substitution: pd.DataFrame,
    params: dict,
    rng: np.random.Generator,
    context_cols: list[str],
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    target_col = "standardized_titer"

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(df_binary)
    sub_cols = sorted([c for c in df_substitution.columns if c.startswith("sub_")])

    site_state_cols = context_cols + change_cols + virus_aa_cols + serum_aa_cols
    binary_cols = context_cols + change_cols
    substitution_cols = context_cols + sub_cols

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
    )
    y_train_corr = y_train - train_base
    selection = inner_selection(
        train_bin, "virusStrain",
        lambda a, b: fit_split_additive(a, b, target_col),
        random_state=params["random_state"],
    )

    results = []
    pred_meta_cols = [
        "virusStrain",
        "serumStrain",
        "virusSequenceId",
        "serumSequenceId",
        "standardized_titer",
        "context_corrected_titer",
    ] + [col for col in context_cols if col in val_bin.columns]
    pred_frame = val_bin[pred_meta_cols].copy()
    pred_frame["additive_only_pred"] = val_base
    results.append({"model_type": "additive_only", **compute_metrics(y_val, val_base), "best_iteration": 0})

    X_context_train = prepare_frame(train_bin, context_cols)
    X_context_val = prepare_frame(val_bin, context_cols)
    context_cats = [c for c in X_context_train.columns if str(X_context_train[c].dtype) == "category"]
    model_context, pred_context_corr = fit_predict(
        X_train=X_context_train,
        y_train=y_train_corr,
        X_val=X_context_val,
        selection=selection,
        params=params,
        categorical_cols=context_cats,
    )
    pred_context = val_base + pred_context_corr
    pred_frame["predicted_standardized_titer_context"] = pred_context
    results.append(
        {
            "model_type": "context_only",
            **compute_metrics(y_val, pred_context),
            "best_iteration": int(model_context.n_estimators_),
        }
    )

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
    results.append(
        {
            "model_type": "site_state_context",
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
    results.append(
        {
            "model_type": "binary_site_context",
            **compute_metrics(y_val, pred_bin),
            "best_iteration": int(model_bin.n_estimators_),
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
            "model_type": "randomized_binary_context",
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
        results.append(
            {
                "model_type": "substitution_identity_context",
                **compute_metrics(y_val, pred_sub),
                "best_iteration": int(model_sub.n_estimators_),
            }
        )
    else:
        pred_frame["predicted_standardized_titer_substitution"] = np.nan

    pred_frame["train_context_additive_intercept"] = float(additive["intercept"])
    return results, pred_frame


def run_neher_holdouts(
    df_binary: pd.DataFrame,
    df_substitution: pd.DataFrame,
    params: dict,
    rng_seed: int,
    context_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    idx_all = np.arange(len(df_binary))
    groups = df_binary["virusStrain"].astype(str)

    split_defs = []
    rand_train, rand_test = train_test_split(
        idx_all,
        test_size=0.10,
        random_state=rng_seed,
        shuffle=True,
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
            context_cols=context_cols,
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
                    "n_test_unique_viruses": int(df_binary.iloc[test_idx]["virusStrain"].nunique()),
                }
            )

    return pd.concat(pred_parts, ignore_index=True), pd.DataFrame(metric_rows)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    binary_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_feature_matrix.tsv"
    substitution_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_feature_matrix.tsv"
    df_binary = pd.read_csv(binary_path, sep="\t", low_memory=False)
    df_substitution = pd.read_csv(substitution_path, sep="\t", low_memory=False)

    groups = df_binary["virusStrain"].astype(str)
    n_splits = args.n_splits or int(cfg.get("n_splits", 5))
    unique_groups = groups.nunique()
    if unique_groups < 2:
        raise ValueError("Need at least 2 unique virus-strain groups for grouped CV.")
    if n_splits > unique_groups:
        n_splits = unique_groups

    gkf = GroupKFold(n_splits=n_splits)
    params = get_lightgbm_params(int(cfg.get("random_state", 42)))
    rng_seed = int(cfg.get("random_state", 42))
    context_cols = get_context_cols(cfg)

    results = []
    oof_parts = []
    best_iters: dict[str, list[int]] = {
        "context_only": [],
        "site_state_context": [],
        "binary_site_context": [],
        "substitution_identity_context": [],
    }

    for fold, (train_idx, val_idx) in enumerate(gkf.split(df_binary, groups=groups), start=1):
        split_results, pred_frame = run_models_for_split(
            train_idx=train_idx,
            val_idx=val_idx,
            df_binary=df_binary,
            df_substitution=df_substitution,
            params=params,
            rng=np.random.default_rng(rng_seed + fold),
            context_cols=context_cols,
        )
        for row in split_results:
            row["fold"] = fold
            row["n_train"] = len(train_idx)
            row["n_val"] = len(val_idx)
            results.append(row)
            if row["model_type"] in best_iters and row["best_iteration"] > 0:
                best_iters[row["model_type"]].append(int(row["best_iteration"]))

        pred_frame["fold"] = fold
        oof_parts.append(pred_frame)

    results_df = pd.DataFrame(results)
    cv_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_cv_results.tsv"
    results_df.to_csv(cv_path, sep="\t", index=False)

    oof_df = pd.concat(oof_parts, ignore_index=True)
    oof_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_oof_predictions.tsv"
    oof_df.to_csv(oof_path, sep="\t", index=False)

    summary = (
        results_df.groupby("model_type")[["rmse", "mae", "r2", "spearman"]]
        .agg(["mean", "std"])
        .reset_index()
    )
    summary_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_cv_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    full_additive = fit_additive_effects_multi(
        df=df_binary,
        target_col="standardized_titer",
        factor_cols=ADDITIVE_FACTOR_COLS,
    )
    full_base = apply_additive_effects_multi(
        df_binary,
        factor_cols=ADDITIVE_FACTOR_COLS,
        intercept=float(full_additive["intercept"]),
        factor_effects=full_additive["factor_effects"],
    )
    y_full_corr = df_binary["standardized_titer"].to_numpy(dtype=float) - full_base

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(df_binary)
    sub_cols = sorted([c for c in df_substitution.columns if c.startswith("sub_")])
    feature_sets = {
        "context_model": context_cols,
        "site_state_model": context_cols + change_cols + virus_aa_cols + serum_aa_cols,
        "binary_model": context_cols + change_cols,
    }
    if sub_cols:
        feature_sets["substitution_model"] = context_cols + sub_cols

    final_model_paths = {}
    for model_key, feature_cols in feature_sets.items():
        model_params = params.copy()
        iter_key = {
            "context_model": "context_only",
            "site_state_model": "site_state_context",
            "binary_model": "binary_site_context",
            "substitution_model": "substitution_identity_context",
        }[model_key]
        if best_iters.get(iter_key):
            model_params["n_estimators"] = max(1, int(np.round(np.mean(best_iters[iter_key]))))
        model = lgb.LGBMRegressor(**model_params)
        full_df = df_substitution if model_key == "substitution_model" else df_binary
        X_full = prepare_frame(full_df, feature_cols)
        cat_cols = [c for c in X_full.columns if str(X_full[c].dtype) == "category"]
        model.fit(X_full, y_full_corr, categorical_feature=cat_cols)
        model_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_{model_key}.txt"
        model.booster_.save_model(str(model_path))
        final_model_paths[model_key] = str(model_path)

    additive_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_context_additive_effects_refit.json"
    with open(additive_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "intercept": float(full_additive["intercept"]),
                "factor_cols": ADDITIVE_FACTOR_COLS,
                "factor_effects": full_additive["factor_effects"],
            },
            f,
            indent=2,
        )

    neher_pred, neher_metrics = run_neher_holdouts(
        df_binary=df_binary,
        df_substitution=df_substitution,
        params=params,
        rng_seed=rng_seed,
        context_cols=context_cols,
    )
    neher_pred_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_neher_analog_predictions.tsv"
    neher_pred.to_csv(neher_pred_path, sep="\t", index=False)
    neher_metrics_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_neher_analog_metrics.tsv"
    neher_metrics.to_csv(neher_metrics_path, sep="\t", index=False)

    lines = [
        f"Feature rows: {len(df_binary)}",
        f"Grouped CV folds: {n_splits}",
        f"Unique virus-strain groups: {unique_groups}",
        f"Context columns: {len(context_cols)}",
        f"Distance feature columns: {get_distance_feature_cols(cfg)}",
        f"Binary feature columns: {len(context_cols) + len(change_cols)}",
        f"Site-state feature columns: {len(context_cols) + len(change_cols) + len(virus_aa_cols) + len(serum_aa_cols)}",
        f"Substitution feature columns: {len(context_cols) + len(sub_cols)}",
        f"Saved CV results: {cv_path}",
        f"Saved CV summary: {summary_path}",
        f"Saved OOF predictions: {oof_path}",
        f"Saved refit additive effects: {additive_path}",
        f"Saved Neher-analog predictions: {neher_pred_path}",
        f"Saved Neher-analog metrics: {neher_metrics_path}",
        f"Saved final models: {final_model_paths}",
        "Primary WIC model family is sequence-plus-context only; no virus or serum identity terms were used as predictors.",
    ]
    append_qc_log(cfg, "27_train_wic_model", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
