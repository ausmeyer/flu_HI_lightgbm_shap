#!/usr/bin/env python3
"""Counterfactual test for the HA1 158-160 glycosylation context in WIC.

This analysis uses the filtered WIC HA1 site-state model and observed WIC assay
contexts. For each observed virus-serum row with a K or T at HA1 site 160, it
mutates only site 160 (K160T or T160K), then asks whether the predicted effect is
larger when the mutation creates or removes the N158-X159-S/T160 glycosylation
motif than when the same site-160 substitution does not change motif status.
"""

from __future__ import annotations

import argparse
import json
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

from common import append_qc_log, get_distance_feature_cols, load_config, read_fasta


CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
GLYCAN_SITES = (158, 159, 160)
SITE_160 = 160

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown.json")
    return parser.parse_args()


def get_context_cols(cfg: dict) -> list[str]:
    context_cols = CONTEXT_COLS.copy()
    for col in get_distance_feature_cols(cfg):
        if col not in context_cols:
            context_cols.append(col)
    return context_cols


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
    x = df[feature_cols].copy()
    for col in x.columns:
        if not pd.api.types.is_numeric_dtype(x[col]):
            x[col] = x[col].astype("category")
    return x


def mature_ha1_sequence_map(cfg: dict) -> dict[str, str]:
    out_dir = Path(cfg["output_dir"])
    subtype = str(cfg["subtype"])
    aligned_path = out_dir / f"{subtype}_HA1_aligned.fasta"
    metadata_path = out_dir / f"{subtype}_ha1_alignment_metadata.json"

    seq_map = read_fasta(str(aligned_path))
    with open(metadata_path, "r", encoding="utf-8") as f:
        alignment_meta = json.load(f)

    reference = str(alignment_meta["reference_sequence_id"])
    if reference not in seq_map:
        raise ValueError(f"Reference sequence id not found in aligned WIC HA1 FASTA: {reference}")

    ref_cols = [idx for idx, aa in enumerate(seq_map[reference]) if aa != "-"]
    ha1_len = int(cfg.get("ha1_length", 328))
    if len(ref_cols) < ha1_len:
        raise ValueError(f"Reference HA1 has fewer than {ha1_len} non-gap positions.")
    pos_cols = ref_cols[:ha1_len]

    return {seq_id: "".join(seq[idx] for idx in pos_cols) for seq_id, seq in seq_map.items()}


def has_nlinked_motif_158_160(seq: str) -> bool:
    aa158 = seq[157]
    aa159 = seq[158]
    aa160 = seq[159]
    return aa158 == "N" and aa159 in CANONICAL_AA and aa159 != "P" and aa160 in {"S", "T"}


def motif_state_label(is_present: bool) -> str:
    return "present" if is_present else "absent"


def mutate_site(seq: str, site: int, aa: str) -> str:
    chars = list(seq)
    chars[site - 1] = aa
    return "".join(chars)


def iqr(values: pd.Series) -> float:
    return float(values.quantile(0.75) - values.quantile(0.25))


def q25(values: pd.Series) -> float:
    return float(values.quantile(0.25))


def q75(values: pd.Series) -> float:
    return float(values.quantile(0.75))


def build_counterfactual_metadata(
    feature_df: pd.DataFrame,
    seq_by_id: dict[str, str],
    feature_cols: list[str],
) -> pd.DataFrame:
    feature_set = set(feature_cols)
    rows = []
    for row in feature_df.itertuples(index=True):
        virus_id = str(row.virusSequenceId)
        serum_id = str(row.serumSequenceId)
        virus_seq = seq_by_id.get(virus_id)
        serum_seq = seq_by_id.get(serum_id)
        if virus_seq is None or serum_seq is None:
            continue

        if any(virus_seq[site - 1] not in CANONICAL_AA for site in GLYCAN_SITES):
            continue

        baseline_aa160 = virus_seq[SITE_160 - 1]
        if baseline_aa160 not in {"K", "T"}:
            continue

        comparison_aa160 = "T" if baseline_aa160 == "K" else "K"
        comparison_seq = mutate_site(virus_seq, SITE_160, comparison_aa160)
        baseline_motif = has_nlinked_motif_158_160(virus_seq)
        comparison_motif = has_nlinked_motif_158_160(comparison_seq)
        motif_transition = f"{motif_state_label(baseline_motif)}_to_{motif_state_label(comparison_motif)}"
        substitution = f"{baseline_aa160}{SITE_160}{comparison_aa160}"

        serum_aa160 = serum_seq[SITE_160 - 1]
        required_cols = [f"{baseline_aa160}{SITE_160}", f"{comparison_aa160}{SITE_160}", f"change_{SITE_160}"]
        if serum_aa160 in CANONICAL_AA:
            required_cols.append(f"{SITE_160}{serum_aa160}")
        encoding_complete = all(col in feature_set for col in required_cols)

        rows.append(
            {
                "input_row_index": int(row.Index),
                "virusStrain": str(row.virusStrain),
                "serumStrain": str(row.serumStrain),
                "virusSequenceId": virus_id,
                "serumSequenceId": serum_id,
                "aa158": virus_seq[157],
                "aa159": virus_seq[158],
                "baseline_aa160": baseline_aa160,
                "comparison_aa160": comparison_aa160,
                "serum_aa160": serum_aa160,
                "substitution": substitution,
                "baseline_motif": bool(baseline_motif),
                "comparison_motif": bool(comparison_motif),
                "motif_transition": motif_transition,
                "motif_changed": bool(baseline_motif != comparison_motif),
                "encoding_complete": bool(encoding_complete),
            }
        )

    if not rows:
        raise ValueError("No eligible WIC glycosylation counterfactual rows were generated.")
    return pd.DataFrame(rows)


def update_site160_comparison_features(
    comparison_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    virus_aa_cols: list[str],
) -> pd.DataFrame:
    out = comparison_df.copy()
    site160_virus_cols = [col for col in virus_aa_cols if int(col[1:]) == SITE_160]
    for col in site160_virus_cols:
        out[col] = 0.0

    for comparison_aa in sorted(meta_df["comparison_aa160"].unique()):
        col = f"{comparison_aa}{SITE_160}"
        if col in out.columns:
            out.loc[meta_df["comparison_aa160"].to_numpy() == comparison_aa, col] = 1.0

    change_col = f"change_{SITE_160}"
    if change_col in out.columns:
        serum_aa = meta_df["serum_aa160"].to_numpy()
        target_aa = meta_df["comparison_aa160"].to_numpy()
        valid = np.array([aa in CANONICAL_AA for aa in serum_aa], dtype=bool)
        values = np.full(len(meta_df), np.nan, dtype=float)
        values[valid] = (target_aa[valid] != serum_aa[valid]).astype(float)
        out[change_col] = values
    return out


def predict_counterfactual_effects(
    feature_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    feature_cols: list[str],
    virus_aa_cols: list[str],
    model_path: Path,
) -> pd.DataFrame:
    complete_meta = meta_df[meta_df["encoding_complete"]].reset_index(drop=True).copy()
    if complete_meta.empty:
        raise ValueError("No encoding-complete WIC glycosylation counterfactual rows remain.")

    baseline_df = feature_df.iloc[complete_meta["input_row_index"].to_numpy()].reset_index(drop=True).copy()
    comparison_df = update_site160_comparison_features(
        comparison_df=baseline_df[feature_cols].copy(),
        meta_df=complete_meta,
        virus_aa_cols=virus_aa_cols,
    )

    booster = lgb.Booster(model_file=str(model_path))
    baseline_pred = booster.predict(prepare_frame(baseline_df, feature_cols))
    comparison_pred = booster.predict(prepare_frame(comparison_df, feature_cols))

    effects = complete_meta.copy()
    effects["model_scope"] = "full_wic_filtered_site_state_context"
    effects["baseline_predicted_residual_titer"] = baseline_pred
    effects["comparison_predicted_residual_titer"] = comparison_pred
    effects["predicted_effect"] = comparison_pred - baseline_pred
    effects["glycan_present_minus_absent_effect"] = np.nan

    absent_to_present = effects["motif_transition"] == "absent_to_present"
    present_to_absent = effects["motif_transition"] == "present_to_absent"
    effects.loc[absent_to_present, "glycan_present_minus_absent_effect"] = effects.loc[
        absent_to_present, "predicted_effect"
    ]
    effects.loc[present_to_absent, "glycan_present_minus_absent_effect"] = -effects.loc[
        present_to_absent, "predicted_effect"
    ]

    context_cols = [col for col in CONTEXT_COLS if col in baseline_df.columns]
    for col in context_cols:
        effects[col] = baseline_df[col].to_numpy()
    return effects


def summarize_effects(effects: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    background_cols = [
        "virusSequenceId",
        "virusStrain",
        "substitution",
        "motif_transition",
        "motif_changed",
        "aa158",
        "aa159",
        "baseline_aa160",
        "comparison_aa160",
        "baseline_motif",
        "comparison_motif",
    ]
    background_summary = (
        effects.groupby(background_cols, as_index=False)
        .agg(
            n_assay_contexts=("predicted_effect", "size"),
            n_serum_sequence_contexts=("serumSequenceId", "nunique"),
            mean_effect=("predicted_effect", "mean"),
            median_effect=("predicted_effect", "median"),
            sd_effect=("predicted_effect", "std"),
            q25_effect=("predicted_effect", q25),
            q75_effect=("predicted_effect", q75),
            iqr_effect=("predicted_effect", iqr),
            min_effect=("predicted_effect", "min"),
            max_effect=("predicted_effect", "max"),
            mean_glycan_present_minus_absent_effect=("glycan_present_minus_absent_effect", "mean"),
            median_glycan_present_minus_absent_effect=("glycan_present_minus_absent_effect", "median"),
        )
        .sort_values(["substitution", "motif_transition", "virusSequenceId"])
    )
    background_summary["sd_effect"] = background_summary["sd_effect"].fillna(0.0)

    row_summary = (
        effects.groupby(["substitution", "motif_transition", "motif_changed"], as_index=False)
        .agg(
            n_effect_rows=("predicted_effect", "size"),
            n_backgrounds=("virusSequenceId", "nunique"),
            n_serum_sequence_contexts=("serumSequenceId", "nunique"),
            mean_effect=("predicted_effect", "mean"),
            median_effect=("predicted_effect", "median"),
            sd_effect=("predicted_effect", "std"),
            q25_effect=("predicted_effect", q25),
            q75_effect=("predicted_effect", q75),
            iqr_effect=("predicted_effect", iqr),
            min_effect=("predicted_effect", "min"),
            max_effect=("predicted_effect", "max"),
            mean_glycan_present_minus_absent_effect=("glycan_present_minus_absent_effect", "mean"),
            median_glycan_present_minus_absent_effect=("glycan_present_minus_absent_effect", "median"),
        )
        .sort_values(["substitution", "motif_transition"])
    )
    row_summary["sd_effect"] = row_summary["sd_effect"].fillna(0.0)

    bg_transition = (
        background_summary.groupby(["substitution", "motif_transition", "motif_changed"], as_index=False)
        .agg(
            n_backgrounds=("virusSequenceId", "nunique"),
            median_background_median_effect=("median_effect", "median"),
            mean_background_mean_effect=("mean_effect", "mean"),
            sd_background_median_effect=("median_effect", "std"),
            iqr_background_median_effect=("median_effect", iqr),
            median_background_iqr_effect=("iqr_effect", "median"),
            median_background_glycan_present_minus_absent_effect=(
                "median_glycan_present_minus_absent_effect",
                "median",
            ),
        )
        .sort_values(["substitution", "motif_transition"])
    )
    bg_transition["sd_background_median_effect"] = bg_transition["sd_background_median_effect"].fillna(0.0)
    transition_summary = row_summary.merge(
        bg_transition,
        on=["substitution", "motif_transition", "motif_changed"],
        suffixes=("_row", "_background"),
    )

    comparison_rows = []
    for substitution, sub_df in transition_summary.groupby("substitution"):
        by_transition = sub_df.set_index("motif_transition")
        motif_transition = "absent_to_present" if substitution == "K160T" else "present_to_absent"
        control_transition = "absent_to_absent"
        if motif_transition not in by_transition.index or control_transition not in by_transition.index:
            continue
        motif_row = by_transition.loc[motif_transition]
        control_row = by_transition.loc[control_transition]
        comparison_rows.append(
            {
                "substitution": substitution,
                "motif_changing_transition": motif_transition,
                "control_transition": control_transition,
                "motif_changing_n_backgrounds": int(motif_row["n_backgrounds_background"]),
                "control_n_backgrounds": int(control_row["n_backgrounds_background"]),
                "motif_changing_median_background_median_effect": float(
                    motif_row["median_background_median_effect"]
                ),
                "control_median_background_median_effect": float(control_row["median_background_median_effect"]),
                "motif_minus_control_median_background_median_effect": float(
                    motif_row["median_background_median_effect"] - control_row["median_background_median_effect"]
                ),
                "motif_changing_mean_background_mean_effect": float(motif_row["mean_background_mean_effect"]),
                "control_mean_background_mean_effect": float(control_row["mean_background_mean_effect"]),
                "motif_minus_control_mean_background_mean_effect": float(
                    motif_row["mean_background_mean_effect"] - control_row["mean_background_mean_effect"]
                ),
                "motif_changing_median_background_iqr_effect": float(motif_row["median_background_iqr_effect"]),
                "control_median_background_iqr_effect": float(control_row["median_background_iqr_effect"]),
            }
        )
    comparison = pd.DataFrame(comparison_rows)
    return background_summary, transition_summary, comparison


def plot_background_effects(background_summary: pd.DataFrame, out_path: Path) -> None:
    plot_df = background_summary.copy()
    if plot_df.empty:
        return

    plot_df["contrast_label"] = plot_df["substitution"] + "\n" + plot_df["motif_transition"].str.replace("_", " ")
    order = [
        "K160T\nabsent to present",
        "K160T\nabsent to absent",
        "T160K\npresent to absent",
        "T160K\nabsent to absent",
    ]
    order = [label for label in order if label in set(plot_df["contrast_label"])]

    sns.set_theme(style="ticks", context="paper")
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    sns.boxplot(
        data=plot_df,
        x="contrast_label",
        y="median_effect",
        order=order,
        color="#c8d5e6",
        linewidth=0.9,
        fliersize=0,
        ax=ax,
    )
    sns.stripplot(
        data=plot_df,
        x="contrast_label",
        y="median_effect",
        order=order,
        color="#244b6b",
        alpha=0.28,
        size=2.4,
        jitter=0.28,
        ax=ax,
    )
    ax.axhline(0, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_xlabel("")
    ax.set_ylabel("Background median predicted residual titer effect")
    ax.set_title("WIC filtered site-state model: HA1 158-160 glycosylation context")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    out_dir = Path(cfg["output_dir"])
    subtype = str(cfg["subtype"])

    feature_path = out_dir / f"{subtype}_feature_matrix.tsv"
    model_path = out_dir / f"{subtype}_site_state_model.txt"
    if not feature_path.exists():
        raise FileNotFoundError(f"Feature matrix not found; run before WIC cleanup: {feature_path}")
    if not model_path.exists():
        raise FileNotFoundError(f"Site-state model not found; run after WIC training and before cleanup: {model_path}")

    feature_df = pd.read_csv(feature_path, sep="\t", low_memory=False)
    context_cols = get_context_cols(cfg)
    change_cols, virus_aa_cols, serum_aa_cols = get_site_state_columns(feature_df)
    feature_cols = context_cols + change_cols + virus_aa_cols + serum_aa_cols

    seq_by_id = mature_ha1_sequence_map(cfg)
    meta_df = build_counterfactual_metadata(feature_df=feature_df, seq_by_id=seq_by_id, feature_cols=feature_cols)
    effects = predict_counterfactual_effects(
        feature_df=feature_df,
        meta_df=meta_df,
        feature_cols=feature_cols,
        virus_aa_cols=virus_aa_cols,
        model_path=model_path,
    )
    background_summary, transition_summary, comparison = summarize_effects(effects)

    effects_path = out_dir / f"{subtype}_glycan_context_effects.tsv"
    background_summary_path = out_dir / f"{subtype}_glycan_context_background_summary.tsv"
    transition_summary_path = out_dir / f"{subtype}_glycan_context_transition_summary.tsv"
    comparison_path = out_dir / f"{subtype}_glycan_context_comparison.tsv"
    effects.to_csv(effects_path, sep="\t", index=False)
    background_summary.to_csv(background_summary_path, sep="\t", index=False)
    transition_summary.to_csv(transition_summary_path, sep="\t", index=False)
    comparison.to_csv(comparison_path, sep="\t", index=False)

    figure_path = out_dir / "figures" / f"{subtype}_glycan_context_dependence_158_160.pdf"
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    plot_background_effects(background_summary, figure_path)

    transition_counts = (
        effects.groupby(["substitution", "motif_transition"])["virusSequenceId"]
        .agg(["size", "nunique"])
        .reset_index()
        .rename(columns={"size": "rows", "nunique": "backgrounds"})
    )
    transition_note = "; ".join(
        f"{row.substitution} {row.motif_transition}: rows={int(row.rows)}, backgrounds={int(row.backgrounds)}"
        for row in transition_counts.itertuples(index=False)
    )
    comparison_note = "not available"
    if not comparison.empty:
        comparison_note = "; ".join(
            f"{row.substitution} motif-control background-median diff="
            f"{row.motif_minus_control_median_background_median_effect:.6f}"
            for row in comparison.itertuples(index=False)
        )

    lines = [
        "Target motif: HA1 N158-X159-S/T160",
        "Counterfactual substitutions: K160T and T160K, one virus site mutated per pair",
        f"Training/model feature rows available: {len(feature_df)}",
        f"Eligible observed WIC assay rows before encoding filter: {len(meta_df)}",
        f"Encoding-complete counterfactual rows: {len(effects)}",
        f"Unique virus backgrounds: {effects['virusSequenceId'].nunique()}",
        f"Unique serum sequence contexts: {effects['serumSequenceId'].nunique()}",
        f"Transition counts: {transition_note}",
        f"Motif-change comparisons: {comparison_note}",
        f"Saved glycan context effects: {effects_path}",
        f"Saved background summary: {background_summary_path}",
        f"Saved transition summary: {transition_summary_path}",
        f"Saved comparison summary: {comparison_path}",
        f"Saved glycan context figure: {figure_path}",
    ]
    append_qc_log(cfg, "44_analyze_wic_glycan_context_dependence", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
