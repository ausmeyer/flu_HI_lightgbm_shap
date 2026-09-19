#!/usr/bin/env python3
"""Exact Koel 2019 counterfactual test for context-dependent effects.

This analysis trains the Neher/Bedford site-state LightGBM model on the full
filtered manuscript input, then evaluates exact Koel 2019 Table 1 substitutions
at HA1 sites 145 and 155 in the available historical HA backgrounds. Each
baseline/comparison virus sequence is paired with a fixed panel of observed
Neher/Bedford serum HA sequences so the predicted effect can vary with both HA
background and serum context.
"""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "flu_hi_lightgbm_shap_matplotlib"))

import lightgbm as lgb
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

from model_validation import inner_selection, fit_selected_model

from common import (
    append_qc_log,
    apply_additive_effects,
    fit_additive_effects,
    get_distance_feature_cols,
    get_lightgbm_params,
    load_config,
    read_fasta,
)


TARGET_SITES = (145, 155)
CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Exact Koel et al. 2019 substitutions for the seven Table 1 backgrounds present
# in the Neher/Bedford HA alignment. The direction is always comparison minus
# baseline in the output predicted_effect column.
KOEL_EXACT_CONTRASTS = [
    {
        "koel_cluster": "HK68",
        "background_label": "BI/16190/68",
        "background_strain": "A/Bilthoven/16190/1968",
        "site": 145,
        "wildtype_aa": "S",
        "baseline_aa": "S",
        "comparison_aa": "K",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "HK68",
        "background_label": "BI/16190/68",
        "background_strain": "A/Bilthoven/16190/1968",
        "site": 155,
        "wildtype_aa": "T",
        "baseline_aa": "T",
        "comparison_aa": "Y",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "EN72",
        "background_label": "BI/21793/72",
        "background_strain": "A/Bilthoven/21793/1972",
        "site": 145,
        "wildtype_aa": "S",
        "baseline_aa": "S",
        "comparison_aa": "K",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "EN72",
        "background_label": "BI/21793/72",
        "background_strain": "A/Bilthoven/21793/1972",
        "site": 155,
        "wildtype_aa": "Y",
        "baseline_aa": "Y",
        "comparison_aa": "T",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "VI75",
        "background_label": "BI/1761/76",
        "background_strain": "A/Bilthoven/1761/1976",
        "site": 145,
        "wildtype_aa": "N",
        "baseline_aa": "N",
        "comparison_aa": "K",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "VI75",
        "background_label": "BI/1761/76",
        "background_strain": "A/Bilthoven/1761/1976",
        "site": 155,
        "wildtype_aa": "Y",
        "baseline_aa": "Y",
        "comparison_aa": "T",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "TX77",
        "background_label": "BI/2271/76",
        "background_strain": "A/Bilthoven/2271/1976",
        "site": 145,
        "wildtype_aa": "N",
        "baseline_aa": "N",
        "comparison_aa": "K",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "TX77",
        "background_label": "BI/2271/76",
        "background_strain": "A/Bilthoven/2271/1976",
        "site": 155,
        "wildtype_aa": "Y",
        "baseline_aa": "Y",
        "comparison_aa": "T",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SI87",
        "background_label": "NL/620/89",
        "background_strain": "A/Netherlands/620/1989",
        "site": 145,
        "wildtype_aa": "N",
        "baseline_aa": "N",
        "comparison_aa": "K",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SI87",
        "background_label": "NL/620/89",
        "background_strain": "A/Netherlands/620/1989",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "H",
        "comparison_aa": "Y",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SI87",
        "background_label": "NL/620/89",
        "background_strain": "A/Netherlands/620/1989",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "H",
        "comparison_aa": "T",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "BE89",
        "background_label": "NL/823/92",
        "background_strain": "A/Netherlands/823/1992",
        "site": 145,
        "wildtype_aa": "K",
        "baseline_aa": "K",
        "comparison_aa": "N",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "BE89",
        "background_label": "NL/823/92",
        "background_strain": "A/Netherlands/823/1992",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "H",
        "comparison_aa": "Y",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "BE89",
        "background_label": "NL/823/92",
        "background_strain": "A/Netherlands/823/1992",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "H",
        "comparison_aa": "T",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SY97",
        "background_label": "NL/427/98",
        "background_strain": "A/Netherlands/427/1998",
        "site": 145,
        "wildtype_aa": "K",
        "baseline_aa": "K",
        "comparison_aa": "N",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SY97",
        "background_label": "NL/427/98",
        "background_strain": "A/Netherlands/427/1998",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "H",
        "comparison_aa": "Y",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SY97",
        "background_label": "NL/427/98",
        "background_strain": "A/Netherlands/427/1998",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "H",
        "comparison_aa": "T",
        "contrast_type": "wildtype_vs_mutant",
    },
    {
        "koel_cluster": "SI87",
        "background_label": "NL/620/89",
        "background_strain": "A/Netherlands/620/1989",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "T",
        "comparison_aa": "Y",
        "contrast_type": "mutant_vs_mutant",
    },
    {
        "koel_cluster": "BE89",
        "background_label": "NL/823/92",
        "background_strain": "A/Netherlands/823/1992",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "T",
        "comparison_aa": "Y",
        "contrast_type": "mutant_vs_mutant",
    },
    {
        "koel_cluster": "SY97",
        "background_label": "NL/427/98",
        "background_strain": "A/Netherlands/427/1998",
        "site": 155,
        "wildtype_aa": "H",
        "baseline_aa": "T",
        "comparison_aa": "Y",
        "contrast_type": "mutant_vs_mutant",
    },
]

KOEL_BACKGROUND_NOT_IN_NEHER_ALIGNMENT = [
    "A/Netherlands/233/1982",
    "A/Netherlands/179/1993",
    "A/Netherlands/178/1995",
    "A/Netherlands/213/2003",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--n-splits", type=int, default=None)
    return parser.parse_args()


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
    return df[feature_cols].copy()


def prepare_linear_frame(df: pd.DataFrame, feature_cols: list[str]) -> pd.DataFrame:
    """Encode missing indicator values as absent for the additive linear comparator."""
    return df[feature_cols].fillna(0.0).copy()


def mature_sequence_map(cfg: dict) -> dict[str, str]:
    aligned_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA_aligned.fasta"
    seq_map = read_fasta(str(aligned_path))
    reference = str(cfg["reference_strain"])
    if reference not in seq_map:
        raise ValueError(f"Reference strain not found in aligned FASTA: {reference}")

    ref_cols = [idx for idx, aa in enumerate(seq_map[reference]) if aa != "-"]
    mature_len = int(cfg.get("mature_protein_length", 550))
    if len(ref_cols) < mature_len:
        raise ValueError(f"Reference has fewer than {mature_len} mature positions.")
    pos_cols = ref_cols[:mature_len]

    mature = {}
    for strain, seq in seq_map.items():
        mature[strain] = "".join(seq[idx] for idx in pos_cols)
    return mature


def strain_year_lookup(feature_df: pd.DataFrame) -> dict[str, float]:
    rows = []
    for strain_col, year_col in [
        ("virus_strain_matched", "virusYear"),
        ("serum_strain_matched", "serumYear"),
    ]:
        tmp = feature_df[[strain_col, year_col]].rename(columns={strain_col: "strain", year_col: "year"})
        rows.append(tmp)
    year_df = pd.concat(rows, ignore_index=True)
    year_df["year"] = pd.to_numeric(year_df["year"], errors="coerce")
    return year_df.dropna().groupby("strain")["year"].median().to_dict()


def strain_year_from_name(strain: str) -> float:
    match = re.search(r"/(\d{4})$", strain)
    if match is None:
        return np.nan
    return float(match.group(1))


def lookup_year(strain: str, year_by_strain: dict[str, float]) -> float:
    return float(year_by_strain.get(strain, strain_year_from_name(strain)))


def substitution_label(site: int, baseline_aa: str, comparison_aa: str) -> str:
    return f"{baseline_aa}{site}{comparison_aa}"


def mutate_site(seq: str, site: int, aa: str) -> str:
    if site < 1 or site > len(seq):
        raise ValueError(f"Site {site} is outside sequence length {len(seq)}.")
    chars = list(seq)
    chars[site - 1] = aa
    return "".join(chars)


def build_serum_contexts(
    feature_df: pd.DataFrame,
    seq_by_strain: dict[str, str],
    year_by_strain: dict[str, float],
) -> pd.DataFrame:
    serum_strains = sorted(feature_df["serum_strain_matched"].dropna().astype(str).unique())
    rows = []
    for strain in serum_strains:
        seq = seq_by_strain.get(strain)
        if seq is None:
            continue
        if any(seq[site - 1] not in CANONICAL_AA for site in TARGET_SITES):
            continue
        rows.append(
            {
                "serum_context_strain": strain,
                "serum_context_year": lookup_year(strain, year_by_strain),
                "serum_sequence": seq,
            }
        )
    if not rows:
        raise ValueError("No unambiguous serum sequence contexts were available.")
    return pd.DataFrame(rows).sort_values("serum_context_strain").reset_index(drop=True)


def encode_pair_features(
    virus_seq: str,
    serum_seq: str,
    feature_cols: list[str],
    distance_cols: list[str],
) -> dict[str, float]:
    feature_set = set(feature_cols)
    row = {col: 0.0 for col in feature_cols}
    for site_idx, (virus_aa, serum_aa) in enumerate(zip(virus_seq, serum_seq), start=1):
        if virus_aa in CANONICAL_AA:
            virus_col = f"{virus_aa}{site_idx}"
            if virus_col in feature_set:
                row[virus_col] = 1.0
        if serum_aa in CANONICAL_AA:
            serum_col = f"{site_idx}{serum_aa}"
            if serum_col in feature_set:
                row[serum_col] = 1.0
        change_col = f"change_{site_idx}"
        if change_col in feature_set:
            if virus_aa in CANONICAL_AA and serum_aa in CANONICAL_AA:
                row[change_col] = float(virus_aa != serum_aa)
            else:
                row[change_col] = np.nan

    for distance_col in distance_cols:
        if distance_col in feature_set:
            row[distance_col] = 0.0
    return row


def target_site_encoding_complete(
    site: int,
    baseline_aa: str,
    comparison_aa: str,
    serum_aa: str,
    feature_set: set[str],
) -> bool:
    if serum_aa not in CANONICAL_AA:
        return False
    required_cols = [
        f"{baseline_aa}{site}",
        f"{comparison_aa}{site}",
        f"{site}{serum_aa}",
        f"change_{site}",
    ]
    return all(col in feature_set for col in required_cols)


def build_exact_koel_counterfactual_rows(
    feature_df: pd.DataFrame,
    seq_by_strain: dict[str, str],
    year_by_strain: dict[str, float],
    feature_cols: list[str],
    distance_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    feature_set = set(feature_cols)
    serum_contexts = build_serum_contexts(feature_df, seq_by_strain, year_by_strain)

    for contrast in KOEL_EXACT_CONTRASTS:
        background_strain = str(contrast["background_strain"])
        if background_strain not in seq_by_strain:
            raise ValueError(f"Expected Koel background missing from alignment: {background_strain}")

        background_seq = seq_by_strain[background_strain]
        site = int(contrast["site"])
        wildtype_aa = str(contrast["wildtype_aa"])
        observed_aa = background_seq[site - 1]
        if observed_aa not in CANONICAL_AA:
            raise ValueError(f"Koel background has ambiguous residue at HA1 {site}: {background_strain}")
        if observed_aa != wildtype_aa:
            raise ValueError(
                f"Unexpected residue in {background_strain} at HA1 {site}: "
                f"expected {wildtype_aa}, observed {observed_aa}"
            )

        baseline_aa = str(contrast["baseline_aa"])
        comparison_aa = str(contrast["comparison_aa"])
        baseline_seq = mutate_site(background_seq, site, baseline_aa)
        comparison_seq = mutate_site(background_seq, site, comparison_aa)
        substitution = substitution_label(site, baseline_aa, comparison_aa)
        background_year = lookup_year(background_strain, year_by_strain)

        for serum in serum_contexts.itertuples(index=False):
            serum_seq = str(serum.serum_sequence)
            serum_aa = serum_seq[site - 1]
            encoding_complete = target_site_encoding_complete(
                site=site,
                baseline_aa=baseline_aa,
                comparison_aa=comparison_aa,
                serum_aa=serum_aa,
                feature_set=feature_set,
            )
            base_meta = {
                "koel_cluster": str(contrast["koel_cluster"]),
                "background_label": str(contrast["background_label"]),
                "background_strain": background_strain,
                "background_year": background_year,
                "serum_context_strain": str(serum.serum_context_strain),
                "serum_context_year": float(serum.serum_context_year),
                "site": site,
                "contrast_type": str(contrast["contrast_type"]),
                "substitution": substitution,
                "wildtype_aa": wildtype_aa,
                "baseline_aa": baseline_aa,
                "comparison_aa": comparison_aa,
                "source_aa": baseline_aa,
                "target_aa": comparison_aa,
                "serum_context_aa": serum_aa,
                "encoding_complete": bool(encoding_complete),
            }
            rows.append(
                {
                    **base_meta,
                    "counterfactual_state": "baseline",
                    **encode_pair_features(baseline_seq, serum_seq, feature_cols, distance_cols),
                }
            )
            rows.append(
                {
                    **base_meta,
                    "counterfactual_state": "comparison",
                    **encode_pair_features(comparison_seq, serum_seq, feature_cols, distance_cols),
                }
            )

    if not rows:
        raise ValueError("No exact Koel counterfactual rows were generated.")
    return pd.DataFrame(rows), serum_contexts


def fit_lightgbm(
    train_df: pd.DataFrame,
    val_df: pd.DataFrame | None,
    feature_cols: list[str],
    params: dict,
) -> tuple[lgb.LGBMRegressor, int]:
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
    y_train = train_df["standardized_titer"].to_numpy(dtype=float) - train_base
    X_train = prepare_frame(train_df, feature_cols)

    if val_df is not None:
        def baseline(a, b):
            fitted = fit_additive_effects(a, "standardized_titer", "virus_strain_matched", "serumIsolate")
            def predict(frame):
                return apply_additive_effects(
                    frame, "serumIsolate", "virus_strain_matched", float(fitted["intercept"]),
                    fitted["serum_effects"], fitted["virus_effects"],
                )
            return predict(a), predict(b)
        selection = inner_selection(train_df, "virus_strain_matched", baseline, params["random_state"])
        model = fit_selected_model(X_train, y_train, params, selection)
    else:
        model = lgb.LGBMRegressor(**params)
        model.fit(X_train, y_train)
    return model, int(model.n_estimators_)


def fit_ridge(feature_df: pd.DataFrame, feature_cols: list[str]) -> Ridge:
    additive = fit_additive_effects(
        df=feature_df,
        target_col="standardized_titer",
        virus_col="virus_strain_matched",
        serum_col="serumIsolate",
    )
    full_base = apply_additive_effects(
        feature_df,
        serum_col="serumIsolate",
        virus_col="virus_strain_matched",
        intercept=float(additive["intercept"]),
        serum_effects=additive["serum_effects"],
        virus_effects=additive["virus_effects"],
    )
    y = feature_df["standardized_titer"].to_numpy(dtype=float) - full_base
    model = Ridge(alpha=1.0)
    model.fit(prepare_linear_frame(feature_df, feature_cols), y)
    return model


def counterfactual_effects_for_model(
    model,
    counterfactual_df: pd.DataFrame,
    feature_cols: list[str],
    model_scope: str,
    fold: int | None = None,
    best_iteration: int | None = None,
) -> pd.DataFrame:
    baseline = counterfactual_df[counterfactual_df["counterfactual_state"] == "baseline"].copy()
    comparison = counterfactual_df[counterfactual_df["counterfactual_state"] == "comparison"].copy()
    join_cols = [
        "koel_cluster",
        "background_label",
        "background_strain",
        "background_year",
        "serum_context_strain",
        "serum_context_year",
        "site",
        "contrast_type",
        "substitution",
        "wildtype_aa",
        "baseline_aa",
        "comparison_aa",
        "source_aa",
        "target_aa",
        "serum_context_aa",
        "encoding_complete",
    ]
    baseline = baseline.sort_values(join_cols).reset_index(drop=True)
    comparison = comparison.sort_values(join_cols).reset_index(drop=True)
    if not baseline[join_cols].equals(comparison[join_cols]):
        raise ValueError("Baseline and comparison counterfactual rows are not aligned.")

    X_baseline = prepare_frame(baseline, feature_cols)
    X_comparison = prepare_frame(comparison, feature_cols)
    if isinstance(model, lgb.LGBMRegressor):
        baseline_pred = model.predict(X_baseline, num_iteration=best_iteration)
        comparison_pred = model.predict(X_comparison, num_iteration=best_iteration)
    else:
        baseline_pred = model.predict(prepare_linear_frame(baseline, feature_cols))
        comparison_pred = model.predict(prepare_linear_frame(comparison, feature_cols))

    out = baseline[join_cols].copy()
    out["model_scope"] = model_scope
    out["fold"] = fold if fold is not None else 0
    out["baseline_predicted_residual_titer"] = baseline_pred
    out["comparison_predicted_residual_titer"] = comparison_pred
    out["predicted_effect"] = comparison_pred - baseline_pred
    return out


def summarize_effects(effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    def iqr(values: pd.Series) -> float:
        return float(values.quantile(0.75) - values.quantile(0.25))

    def q25(values: pd.Series) -> float:
        return float(values.quantile(0.25))

    def q75(values: pd.Series) -> float:
        return float(values.quantile(0.75))

    def value_range(values: pd.Series) -> float:
        return float(values.max() - values.min())

    contrast_cols = [
        "model_scope",
        "contrast_type",
        "koel_cluster",
        "background_label",
        "background_strain",
        "background_year",
        "site",
        "substitution",
        "wildtype_aa",
        "baseline_aa",
        "comparison_aa",
        "source_aa",
        "target_aa",
    ]
    sub_summary = (
        effects.groupby(contrast_cols, as_index=False)
        .agg(
            n_effect_rows=("predicted_effect", "size"),
            n_serum_contexts=("serum_context_strain", "nunique"),
            n_folds=("fold", "nunique"),
            mean_effect=("predicted_effect", "mean"),
            median_effect=("predicted_effect", "median"),
            sd_effect=("predicted_effect", "std"),
            q25_effect=("predicted_effect", q25),
            q75_effect=("predicted_effect", q75),
            iqr_effect=("predicted_effect", iqr),
            min_effect=("predicted_effect", "min"),
            max_effect=("predicted_effect", "max"),
            encoding_complete=("encoding_complete", "all"),
        )
        .sort_values(["model_scope", "contrast_type", "site", "background_year", "substitution"])
    )
    sub_summary["sd_effect"] = sub_summary["sd_effect"].fillna(0.0)
    sub_summary["abs_median_effect"] = sub_summary["median_effect"].abs()
    sub_summary["relative_iqr"] = sub_summary["iqr_effect"] / (sub_summary["abs_median_effect"] + 1e-9)

    site_summary = (
        sub_summary.groupby(["model_scope", "contrast_type", "site"], as_index=False)
        .agg(
            n_background_contrasts=("substitution", "size"),
            n_backgrounds=("background_strain", "nunique"),
            n_substitutions=("substitution", "nunique"),
            n_koel_clusters=("koel_cluster", "nunique"),
            mean_background_effect=("mean_effect", "mean"),
            median_background_effect=("median_effect", "median"),
            mean_abs_background_effect=("abs_median_effect", "mean"),
            sd_background_mean_effect=("mean_effect", "std"),
            sd_background_median_effect=("median_effect", "std"),
            range_background_mean_effect=("mean_effect", value_range),
            range_background_median_effect=("median_effect", value_range),
            median_iqr_serum_effect=("iqr_effect", "median"),
            max_iqr_serum_effect=("iqr_effect", "max"),
            median_relative_iqr=("relative_iqr", "median"),
            encoding_complete=("encoding_complete", "all"),
        )
        .sort_values(["model_scope", "contrast_type", "site"])
    )
    site_summary["sd_background_mean_effect"] = site_summary["sd_background_mean_effect"].fillna(0.0)
    site_summary["sd_background_median_effect"] = site_summary["sd_background_median_effect"].fillna(0.0)

    contrast_rows = []
    wildtype_site_summary = site_summary[site_summary["contrast_type"] == "wildtype_vs_mutant"]
    for model_scope, scope_df in wildtype_site_summary.groupby("model_scope"):
        by_site = scope_df.set_index("site")
        if 145 not in by_site.index or 155 not in by_site.index:
            continue
        contrast_rows.append(
            {
                "model_scope": model_scope,
                "contrast_type": "wildtype_vs_mutant",
                "site_145_n_background_contrasts": int(by_site.loc[145, "n_background_contrasts"]),
                "site_155_n_background_contrasts": int(by_site.loc[155, "n_background_contrasts"]),
                "site_145_sd_background_median_effect": float(by_site.loc[145, "sd_background_median_effect"]),
                "site_155_sd_background_median_effect": float(by_site.loc[155, "sd_background_median_effect"]),
                "site_155_minus_145_sd_background_median_effect": float(
                    by_site.loc[155, "sd_background_median_effect"]
                    - by_site.loc[145, "sd_background_median_effect"]
                ),
                "site_145_sd_background_mean_effect": float(by_site.loc[145, "sd_background_mean_effect"]),
                "site_155_sd_background_mean_effect": float(by_site.loc[155, "sd_background_mean_effect"]),
                "site_155_minus_145_sd_background_mean_effect": float(
                    by_site.loc[155, "sd_background_mean_effect"] - by_site.loc[145, "sd_background_mean_effect"]
                ),
                "site_145_range_background_median_effect": float(
                    by_site.loc[145, "range_background_median_effect"]
                ),
                "site_155_range_background_median_effect": float(
                    by_site.loc[155, "range_background_median_effect"]
                ),
                "site_155_minus_145_range_background_median_effect": float(
                    by_site.loc[155, "range_background_median_effect"]
                    - by_site.loc[145, "range_background_median_effect"]
                ),
                "site_145_range_background_mean_effect": float(
                    by_site.loc[145, "range_background_mean_effect"]
                ),
                "site_155_range_background_mean_effect": float(
                    by_site.loc[155, "range_background_mean_effect"]
                ),
                "site_155_minus_145_range_background_mean_effect": float(
                    by_site.loc[155, "range_background_mean_effect"] - by_site.loc[145, "range_background_mean_effect"]
                ),
                "site_145_median_iqr_serum_effect": float(by_site.loc[145, "median_iqr_serum_effect"]),
                "site_155_median_iqr_serum_effect": float(by_site.loc[155, "median_iqr_serum_effect"]),
                "site_155_minus_145_median_iqr_serum_effect": float(
                    by_site.loc[155, "median_iqr_serum_effect"]
                    - by_site.loc[145, "median_iqr_serum_effect"]
                ),
            }
        )
    contrast = pd.DataFrame(contrast_rows)
    return sub_summary, site_summary, contrast


def plot_effects(sub_summary: pd.DataFrame, out_path: Path) -> None:
    plot_df = sub_summary[sub_summary["model_scope"] == "full_lightgbm"].copy()
    plot_df = plot_df.dropna(subset=["background_year"])
    if plot_df.empty:
        return

    sns.set_theme(style="ticks", context="paper")
    substitutions = sorted(plot_df["substitution"].unique())
    palette = dict(zip(substitutions, sns.color_palette("tab10", n_colors=len(substitutions))))
    markers = {"wildtype_vs_mutant": "o", "mutant_vs_mutant": "X"}

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), sharey=False)
    for ax, site in zip(axes, TARGET_SITES):
        site_df = plot_df[plot_df["site"] == site].copy()
        labels_seen: set[str] = set()
        for row in site_df.itertuples(index=False):
            color = palette[str(row.substitution)]
            marker = markers.get(str(row.contrast_type), "o")
            x_value = float(row.background_year)
            y_value = float(row.median_effect)
            y_low = max(0.0, y_value - float(row.q25_effect))
            y_high = max(0.0, float(row.q75_effect) - y_value)
            label = f"{row.substitution} ({str(row.contrast_type).replace('_', ' ')})"
            plot_label = label if label not in labels_seen else None
            labels_seen.add(label)
            ax.errorbar(
                x_value,
                y_value,
                yerr=np.array([[y_low], [y_high]]),
                fmt=marker,
                color=color,
                markeredgecolor=color,
                markersize=5,
                linewidth=0.8,
                capsize=2,
                alpha=0.85,
                label=plot_label,
            )
        ax.axhline(0, color="#666666", linewidth=0.8, linestyle="--")
        ax.set_title(f"Site {site}")
        ax.set_xlabel("HA background year")
        ax.set_ylabel("Median predicted residual titer effect")
        ax.legend(title="Contrast", fontsize=7, title_fontsize=8, frameon=False)
    fig.suptitle("Exact Koel 2019 counterfactual effects across HA backgrounds", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    out_dir = Path(cfg["output_dir"])
    subtype = str(cfg["subtype"])
    feature_path = out_dir / f"{subtype}_feature_matrix.csv"
    feature_df = pd.read_csv(feature_path, low_memory=False)

    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(feature_df)
    distance_cols = get_distance_feature_cols(cfg)
    feature_cols = change_cols + virus_aa_cols + serum_aa_cols + distance_cols

    seq_by_strain = mature_sequence_map(cfg)
    year_by_strain = strain_year_lookup(feature_df)
    counterfactual_all, serum_contexts = build_exact_koel_counterfactual_rows(
        feature_df=feature_df,
        seq_by_strain=seq_by_strain,
        year_by_strain=year_by_strain,
        feature_cols=feature_cols,
        distance_cols=distance_cols,
    )
    counterfactual_df = counterfactual_all[counterfactual_all["encoding_complete"]].copy()
    if counterfactual_df.empty:
        raise ValueError("No encoding-complete exact Koel counterfactual rows remain.")

    params = get_lightgbm_params(int(cfg.get("random_state", 42)))
    params["n_jobs"] = 1
    n_splits = args.n_splits or int(cfg.get("n_splits", 5))
    groups = feature_df["virus_strain_matched"].astype(str)
    if n_splits > groups.nunique():
        n_splits = int(groups.nunique())

    effects = []
    best_iterations = []
    gkf = GroupKFold(n_splits=n_splits)
    for fold, (train_idx, test_idx) in enumerate(gkf.split(feature_df, groups=groups), start=1):
        train_df = feature_df.iloc[train_idx].copy()
        test_df = feature_df.iloc[test_idx].copy()
        model, best_iteration = fit_lightgbm(
            train_df=train_df,
            val_df=test_df,
            feature_cols=feature_cols,
            params=params,
        )
        best_iterations.append(best_iteration)
        effects.append(
            counterfactual_effects_for_model(
                model=model,
                counterfactual_df=counterfactual_df,
                feature_cols=feature_cols,
                model_scope="cv_lightgbm",
                fold=fold,
                best_iteration=best_iteration,
            )
        )

    full_params = params.copy()
    if best_iterations:
        full_params["n_estimators"] = max(1, int(np.round(np.mean(best_iterations))))
    full_model, full_best_iteration = fit_lightgbm(
        train_df=feature_df,
        val_df=None,
        feature_cols=feature_cols,
        params=full_params,
    )
    effects.append(
        counterfactual_effects_for_model(
            model=full_model,
            counterfactual_df=counterfactual_df,
            feature_cols=feature_cols,
            model_scope="full_lightgbm",
            fold=None,
            best_iteration=full_best_iteration,
        )
    )

    ridge = fit_ridge(feature_df=feature_df, feature_cols=feature_cols)
    effects.append(
        counterfactual_effects_for_model(
            model=ridge,
            counterfactual_df=counterfactual_df,
            feature_cols=feature_cols,
            model_scope="full_ridge_additive_features",
            fold=None,
            best_iteration=None,
        )
    )

    effects_df = pd.concat(effects, ignore_index=True)
    sub_summary, site_summary, contrast = summarize_effects(effects_df)

    effects_path = out_dir / f"{subtype}_koel_context_effects.tsv"
    sub_summary_path = out_dir / f"{subtype}_koel_context_substitution_summary.tsv"
    site_summary_path = out_dir / f"{subtype}_koel_context_site_summary.tsv"
    contrast_path = out_dir / f"{subtype}_koel_context_site_contrast.tsv"
    effects_df.to_csv(effects_path, sep="\t", index=False)
    sub_summary.to_csv(sub_summary_path, sep="\t", index=False)
    site_summary.to_csv(site_summary_path, sep="\t", index=False)
    contrast.to_csv(contrast_path, sep="\t", index=False)

    figure_path = Path(cfg["figures_dir"]) / f"{subtype}_koel_context_dependence_145_155.pdf"
    plot_effects(sub_summary, figure_path)

    full_contrast = contrast[contrast["model_scope"] == "full_lightgbm"]
    contrast_note = "not available"
    if not full_contrast.empty:
        row = full_contrast.iloc[0]
        contrast_note = (
            f"full_lightgbm SD background-median site155-site145="
            f"{row['site_155_minus_145_sd_background_median_effect']:.6f}; "
            f"SD background-mean site155-site145={row['site_155_minus_145_sd_background_mean_effect']:.6f}; "
            f"median range difference={row['site_155_minus_145_range_background_median_effect']:.6f}; "
            f"mean range difference={row['site_155_minus_145_range_background_mean_effect']:.6f}; "
            f"serum-context median IQR difference={row['site_155_minus_145_median_iqr_serum_effect']:.6f}"
        )

    contrast_design = pd.DataFrame(KOEL_EXACT_CONTRASTS)
    retained_pairs = len(counterfactual_df) // 2
    generated_pairs = len(counterfactual_all) // 2
    lines = [
        "Target Koel 2019 sites: 145, 155",
        f"Training rows: {len(feature_df)}",
        f"Exact Koel backgrounds evaluated: {contrast_design['background_strain'].nunique()}",
        f"Koel backgrounds not in Neher/Bedford alignment: {len(KOEL_BACKGROUND_NOT_IN_NEHER_ALIGNMENT)}",
        f"Exact Koel contrasts evaluated: {len(contrast_design)}",
        f"Wildtype-vs-mutant contrasts: {(contrast_design['contrast_type'] == 'wildtype_vs_mutant').sum()}",
        f"Mutant-vs-mutant site 155 contrasts: {(contrast_design['contrast_type'] == 'mutant_vs_mutant').sum()}",
        f"Serum sequence contexts: {len(serum_contexts)}",
        f"Encoding-complete counterfactual pairs: {retained_pairs} of {generated_pairs}",
        f"Grouped CV folds: {n_splits}",
        f"Mean LightGBM best iteration: {float(np.mean(best_iterations)):.2f}",
        f"Context contrast: {contrast_note}",
        f"Saved context effects: {effects_path}",
        f"Saved substitution summary: {sub_summary_path}",
        f"Saved site summary: {site_summary_path}",
        f"Saved site contrast: {contrast_path}",
        f"Saved context-dependence figure: {figure_path}",
    ]
    append_qc_log(cfg, "43_analyze_koel_context_dependence", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
