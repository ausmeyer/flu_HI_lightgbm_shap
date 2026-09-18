"""Select boosting rounds using only the outer training viruses."""

import lightgbm as lgb
import numpy as np
from sklearn.model_selection import GroupShuffleSplit


def inner_selection(train_df, group_col, baseline_fitter, random_state=42):
    """Refit nuisance effects inside a grouped, training-only 80/20 split."""
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=random_state)
    fit_idx, stop_idx = next(splitter.split(train_df, groups=train_df[group_col].astype(str)))
    fit_df, stop_df = train_df.iloc[fit_idx], train_df.iloc[stop_idx]
    fit_base, stop_base = baseline_fitter(fit_df, stop_df)[:2]
    return (fit_idx, stop_idx,
            fit_df["standardized_titer"].to_numpy(dtype=float) - fit_base,
            stop_df["standardized_titer"].to_numpy(dtype=float) - stop_base)


def fit_selected_model(X_train, y_train, params, selection, categorical_cols=()):
    """Choose rounds internally, then refit on all outer training observations."""
    fit_idx, stop_idx, fit_target, stop_target = selection
    selector = lgb.LGBMRegressor(**params)
    selector.fit(
        X_train.iloc[fit_idx], fit_target,
        eval_set=[(X_train.iloc[stop_idx], stop_target)],
        categorical_feature=list(categorical_cols),
        callbacks=[lgb.early_stopping(50, verbose=False)],
    )
    selected_params = dict(params)
    selected_params["n_estimators"] = int(selector.best_iteration_ or params["n_estimators"])
    model = lgb.LGBMRegressor(**selected_params)
    model.fit(X_train, np.asarray(y_train), categorical_feature=list(categorical_cols))
    return model
