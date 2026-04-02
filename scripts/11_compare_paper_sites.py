#!/usr/bin/env python3
"""Compare ranked sites from the local model with H3 sites named in Neher et al. 2016."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import append_qc_log, load_config
from paper_sites import NEHER2016_H3_SITE_ROWS, SHAH2024_H3_SITE_ROWS, WIC2023_H3_SITE_ROWS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--top-n", type=int, default=30)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    site_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_site_importance.csv"
    substitution_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_site_importance.csv"
    if site_path.exists():
        site_importance = pd.read_csv(site_path)
    else:
        site_importance = pd.read_csv(Path(cfg["output_dir"]) / f"{cfg['subtype']}_site_summary.tsv", sep="\t")
    if substitution_path.exists():
        substitution_importance = pd.read_csv(substitution_path)
    else:
        substitution_importance = pd.read_csv(
            Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_site_summary.tsv",
            sep="\t",
        )
    site_importance["rank"] = site_importance["rank"].astype(int)
    substitution_importance["rank"] = substitution_importance["rank"].astype(int)

    neher_sites = pd.DataFrame(NEHER2016_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "neher2016_category", "paper_note": "neher2016_note"}
    )
    neher_sites["named_in_neher2016"] = True
    harvey_sites = pd.DataFrame(WIC2023_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "wic2023_category", "paper_note": "wic2023_note"}
    )
    harvey_sites["named_in_wic2023"] = True
    shah_sites = pd.DataFrame(SHAH2024_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "shah2024_category", "paper_note": "shah2024_note"}
    )
    shah_sites["named_in_shah2024"] = True

    comparison = (
        site_importance[["site", "rank", "mean_abs_shap", "is_koel_site"]]
        .rename(
                columns={
                    "rank": "site_state_rank",
                    "mean_abs_shap": "site_state_mean_abs_shap",
                }
        )
        .merge(neher_sites, on="site", how="outer")
        .merge(harvey_sites, on="site", how="outer")
        .merge(shah_sites, on="site", how="outer")
        .merge(
            substitution_importance[["site", "rank", "mean_abs_shap"]].rename(
                columns={
                    "rank": "substitution_rank",
                    "mean_abs_shap": "substitution_mean_abs_shap",
                }
            ),
            on="site",
            how="outer",
        )
    )
    comparison["is_koel_site"] = comparison["is_koel_site"].fillna(False)
    comparison["named_in_neher2016"] = comparison["named_in_neher2016"].fillna(False)
    comparison["named_in_wic2023"] = comparison["named_in_wic2023"].fillna(False)
    comparison["named_in_shah2024"] = comparison["named_in_shah2024"].fillna(False)
    comparison["named_in_any_reference"] = (
        comparison["named_in_neher2016"] | comparison["named_in_wic2023"] | comparison["named_in_shah2024"]
    )
    comparison["paper_named"] = comparison["named_in_neher2016"]
    comparison["paper_category"] = comparison["neher2016_category"]
    comparison["paper_note"] = comparison["neher2016_note"]
    comparison["site_state_in_model_top_n"] = comparison["site_state_rank"].le(args.top_n).fillna(False)
    comparison["substitution_in_model_top_n"] = comparison["substitution_rank"].le(args.top_n).fillna(False)
    comparison["sort_rank"] = comparison[["site_state_rank", "substitution_rank"]].min(axis=1, skipna=True)
    comparison = comparison.sort_values(["sort_rank", "site"], na_position="last")
    comparison = comparison[
        [
            "site",
            "site_state_rank",
            "site_state_mean_abs_shap",
            "site_state_in_model_top_n",
            "substitution_rank",
            "substitution_mean_abs_shap",
            "substitution_in_model_top_n",
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
            "paper_category",
            "paper_note",
        ]
    ]

    out_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_paper_site_comparison.tsv"
    comparison.to_csv(out_path, sep="\t", index=False)

    site_state_overlap = comparison[
        comparison["site_state_in_model_top_n"] & comparison["named_in_neher2016"]
    ]
    substitution_overlap = comparison[
        comparison["substitution_in_model_top_n"] & comparison["named_in_neher2016"]
    ]
    harvey_site_overlap = comparison[
        comparison["site_state_in_model_top_n"] & comparison["named_in_wic2023"]
    ]
    shah_site_overlap = comparison[
        comparison["site_state_in_model_top_n"] & comparison["named_in_shah2024"]
    ]
    lines = [
        f"Neher 2016 reference sites: {len(neher_sites)}",
        f"WIC 2023 reference sites: {len(harvey_sites)}",
        f"Shah 2024 reference sites: {len(shah_sites)}",
        f"Site-state overlap with Neher/Bedford in model top-{args.top_n}: {len(site_state_overlap)} -> {sorted(site_state_overlap['site'].tolist())}",
        f"Substitution overlap with Neher/Bedford in model top-{args.top_n}: {len(substitution_overlap)} -> {sorted(substitution_overlap['site'].tolist())}",
        f"Site-state overlap with Harvey/WIC in model top-{args.top_n}: {len(harvey_site_overlap)} -> {sorted(harvey_site_overlap['site'].tolist())}",
        f"Site-state overlap with Shah/WIC in model top-{args.top_n}: {len(shah_site_overlap)} -> {sorted(shah_site_overlap['site'].tolist())}",
        f"Output: {out_path}",
    ]
    append_qc_log(cfg, "11_compare_paper_sites", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
