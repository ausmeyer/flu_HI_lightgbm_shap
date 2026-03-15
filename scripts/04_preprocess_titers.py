#!/usr/bin/env python3
"""Build paper-style standardized titers and additive serum/virus effects."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from common import (
    annotate_homologous_titers,
    append_qc_log,
    build_lab_from_source,
    fit_additive_effects,
    load_config,
    parse_censored_titer,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument(
        "--serum-id-column",
        default="serumIsolate",
        help="Column used for serum potency estimation.",
    )
    parser.add_argument(
        "--homologous-fallback",
        default=None,
        help="Fallback for sera without a homologous measurement. Defaults to config value.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    matched_path = f"{cfg['output_dir']}/{cfg['subtype']}_matched_pairs.csv"
    df = pd.read_csv(matched_path)

    if args.serum_id_column not in df.columns:
        raise ValueError(f"serum id column not found: {args.serum_id_column}")

    df["titer_numeric"] = df["titer"].astype(str).map(parse_censored_titer)
    if (df["titer_numeric"] <= 0).any():
        raise ValueError("Non-positive titer values encountered after censoring conversion.")

    df["log2_titer"] = np.log2(df["titer_numeric"])
    df["lab"] = df["source"].astype(str).map(build_lab_from_source)
    fallback = args.homologous_fallback or str(cfg.get("homologous_fallback", "max_observed"))

    standardized = annotate_homologous_titers(
        df=df,
        serum_id_col=args.serum_id_column,
        serum_strain_col="serum_strain_matched",
        virus_strain_col="virus_strain_matched",
        log_titer_col="log2_titer",
        fallback=fallback,
    )

    additive = fit_additive_effects(
        df=standardized,
        target_col="standardized_titer",
        virus_col="virus_strain_matched",
        serum_col=args.serum_id_column,
    )
    standardized["serum_potency"] = standardized[args.serum_id_column].map(additive["serum_effects"])
    standardized["virus_avidity"] = standardized["virus_strain_matched"].map(additive["virus_effects"])
    standardized["additive_intercept"] = additive["intercept"]
    standardized["additive_fitted"] = additive["fitted"]
    standardized["corrected_titer"] = additive["residual"]

    out_path = f"{cfg['output_dir']}/{cfg['subtype']}_matched_titers.csv"
    standardized.to_csv(out_path, index=False)

    titer_vals = df["titer"].astype(str)
    n_lt = int(titer_vals.str.startswith("<").sum())
    n_gt = int(titer_vals.str.startswith(">").sum())
    proxy_counts = dict(standardized["homologous_titer_source"].value_counts())

    lines = [
        f"Input matched rows: {len(standardized)}",
        f"Censored titers '<X': {n_lt}",
        f"Censored titers '>X': {n_gt}",
        f"Serum ID column: {args.serum_id_column}",
        f"Homologous fallback mode: {fallback}",
        f"Homologous proxy counts: {proxy_counts}",
        f"Standardized titer mean: {standardized['standardized_titer'].mean():.6f}",
        f"Standardized titer std: {standardized['standardized_titer'].std():.6f}",
        f"Additive intercept: {additive['intercept']:.6f}",
        f"Additive sera: {additive['n_sera']}",
        f"Additive viruses: {additive['n_viruses']}",
        f"Corrected titer mean: {standardized['corrected_titer'].mean():.6f}",
        f"Corrected titer std: {standardized['corrected_titer'].std():.6f}",
        f"Lab counts: {dict(standardized['lab'].value_counts())}",
        f"Output: {out_path}",
    ]
    append_qc_log(cfg, "04_preprocess_titers", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
