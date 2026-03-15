#!/usr/bin/env python3
"""Preprocess final WIC HI data for HA1 modeling."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    annotate_homologous_titers,
    append_qc_log,
    fit_additive_effects_multi,
    load_config,
    parse_censored_titer,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    return parser.parse_args()


def normalize_text(value: object, default: str = "UNKNOWN") -> str:
    if pd.isna(value):
        return default
    s = str(value).strip()
    return s if s else default


def normalize_passage(value: object) -> str:
    s = normalize_text(value)
    if s.lower() in {"na", "nan", "none", ""}:
        return "UNKNOWN"
    return s.upper()


def normalize_antiserum(value: object) -> str:
    s = normalize_text(value)
    if s.lower() in {"na", "nan", "none", ""}:
        return "UNKNOWN"
    return s.lower()


def normalize_rbc(value: object) -> str:
    s = normalize_text(value)
    if s.lower() in {"na", "nan", "none", ""}:
        return "UNKNOWN"
    return s


def normalize_assay_date_label(value: object) -> str:
    s = normalize_text(value)
    if s.lower() in {"na", "nan", "none", ""}:
        return "UNKNOWN"
    return s


def parse_rbc_species(rbc: str) -> str:
    lower = rbc.lower()
    if "guinea pig" in lower:
        return "guinea_pig"
    if "human" in lower:
        return "human"
    if "turkey" in lower:
        return "turkey"
    return "unknown"


def parse_oseltamivir_nm(rbc: str) -> float:
    lower = rbc.lower()
    if "without oseltamivir" in lower:
        return 0.0
    if "oseltamivir" not in lower:
        return np.nan
    token = lower.replace("nm", " nm").split()
    for i, part in enumerate(token[:-1]):
        try:
            value = float(part)
        except ValueError:
            continue
        if token[i + 1] == "nm":
            return value
    return np.nan


def normalize_wic_titer_string(value: object) -> str:
    s = str(value).strip()
    if s == "<":
        return "<40"
    if s == ">":
        return ">5120"
    return s


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    df = pd.read_csv(cfg["final_hi_path"], sep="\t")
    input_rows = int(len(df))
    df["titer_raw"] = df["titer"].astype(str).str.strip()
    df["titer"] = df["titer_raw"].map(normalize_wic_titer_string)
    normalized_bare_lt = int(df["titer_raw"].eq("<").sum())
    normalized_bare_gt = int(df["titer_raw"].eq(">").sum())
    df["titer_numeric"] = df["titer"].map(parse_censored_titer)
    if (df["titer_numeric"] <= 0).any():
        raise ValueError("Non-positive WIC titers encountered after censoring conversion.")

    df["log2_titer"] = np.log2(df["titer_numeric"])
    df["virus_passage_class"] = df["virus_passage"].map(normalize_passage)
    df["reference_passage_class"] = df["reference_passage"].map(normalize_passage)
    df["antiserum_class"] = df["antiserum"].map(normalize_antiserum)
    df["rbc_type"] = df["RBC"].map(normalize_rbc)
    df["rbc_species"] = df["rbc_type"].map(parse_rbc_species)
    df["oseltamivir_nm"] = df["rbc_type"].map(parse_oseltamivir_nm)
    df["oseltamivir_present"] = df["oseltamivir_nm"].fillna(0.0).gt(0).astype(int)
    df["virus_passage_known"] = df["virus_passage_class"].ne("UNKNOWN").astype(int)
    df["reference_passage_known"] = df["reference_passage_class"].ne("UNKNOWN").astype(int)
    # WIC mixes real dates with assay-batch labels like "1993-aaa"; use the raw label as the nuisance date effect.
    df["assay_date_label"] = df["date"].map(normalize_assay_date_label)
    df["temporal_distance"] = (
        pd.to_numeric(df["virusYear"], errors="coerce") - pd.to_numeric(df["serumYear"], errors="coerce")
    ).abs()
    exclude_passages = {
        normalize_passage(value)
        for value in cfg.get("exclude_passage_classes", [])
    }
    excluded_rows = 0
    if exclude_passages:
        keep_mask = (
            ~df["virus_passage_class"].isin(exclude_passages)
            & ~df["reference_passage_class"].isin(exclude_passages)
        )
        excluded_rows = int((~keep_mask).sum())
        df = df.loc[keep_mask].copy()

    standardized = annotate_homologous_titers(
        df=df,
        serum_id_col="serumSequenceId",
        serum_strain_col="serumSequenceId",
        virus_strain_col="virusSequenceId",
        log_titer_col="log2_titer",
        fallback=str(cfg.get("homologous_fallback", "max_observed")),
    )

    additive_factor_cols = [
        "serumSequenceId",
        "assay_date_label",
        "virus_passage_class",
        "reference_passage_class",
        "antiserum_class",
        "rbc_type",
    ]
    additive = fit_additive_effects_multi(
        df=standardized,
        target_col="standardized_titer",
        factor_cols=additive_factor_cols,
    )
    standardized["context_additive_intercept"] = float(additive["intercept"])
    standardized["context_additive_fitted"] = additive["fitted"]
    standardized["context_corrected_titer"] = additive["residual"]

    out_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_matched_titers.tsv"
    standardized.to_csv(out_path, sep="\t", index=False)

    additive_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_context_additive_effects.json"
    with open(additive_out, "w", encoding="utf-8") as f:
        json.dump(
            {
                "intercept": float(additive["intercept"]),
                "factor_cols": additive_factor_cols,
                "factor_effects": additive["factor_effects"],
                "n_rows": int(additive["n_rows"]),
                "n_factors": int(additive["n_factors"]),
                "n_levels_total": int(additive["n_levels_total"]),
            },
            f,
            indent=2,
        )

    proxy_counts = dict(standardized["homologous_titer_source"].value_counts())
    lines = [
        f"Input final HI rows before passage filtering: {input_rows}",
        f"Retained final HI rows after passage filtering: {len(df)}",
        f"Normalized bare '<' titers to '<40': {normalized_bare_lt}",
        f"Normalized bare '>' titers to '>5120': {normalized_bare_gt}",
        f"Excluded passage classes: {sorted(exclude_passages) if exclude_passages else 'none'}",
        f"Rows excluded by passage filter: {excluded_rows}",
        f"Homologous standardization id: serumSequenceId",
        f"Homologous proxy counts: {proxy_counts}",
        f"Standardized titer mean: {standardized['standardized_titer'].mean():.6f}",
        f"Standardized titer std: {standardized['standardized_titer'].std():.6f}",
        f"Context additive factors: {additive_factor_cols}",
        f"Context additive intercept: {additive['intercept']:.6f}",
        f"Context corrected titer mean: {standardized['context_corrected_titer'].mean():.6f}",
        f"Context corrected titer std: {standardized['context_corrected_titer'].std():.6f}",
        f"Matched titers output: {out_path}",
        f"Context additive effects: {additive_out}",
    ]
    append_qc_log(cfg, "25_preprocess_wic_model_titers", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
