#!/usr/bin/env python3
"""Compare WIC site rankings to Neher 2016 and the original WIC paper."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import append_qc_log, load_config
from paper_sites import NEHER2016_H3_SITE_ROWS, SHAH2024_H3_SITE_ROWS, WIC2023_H3_SITE_ROWS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    parser.add_argument("--top-n", type=int, default=30)
    return parser.parse_args()


def build_reference_comparison(
    site_state_importance: pd.DataFrame,
    substitution_importance: pd.DataFrame,
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

    comparison = (
        reference_df.merge(
            site_state_importance[["site", "rank", "mean_abs_shap"]].rename(
                columns={
                    "rank": "site_state_rank",
                    "mean_abs_shap": "site_state_mean_abs_shap",
                }
            ),
            on="site",
            how="outer",
        )
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
    comparison[named_col] = comparison[named_col].fillna(False)
    comparison["site_state_in_model_top_n"] = comparison["site_state_rank"].le(top_n).fillna(False)
    comparison["substitution_in_model_top_n"] = comparison["substitution_rank"].le(top_n).fillna(False)
    comparison["sort_rank"] = comparison[["site_state_rank", "substitution_rank"]].min(axis=1, skipna=True)
    comparison = comparison.sort_values(["sort_rank", "site"], na_position="last")
    return comparison[
        [
            "site",
            "site_state_rank",
            "site_state_mean_abs_shap",
            "site_state_in_model_top_n",
            "substitution_rank",
            "substitution_mean_abs_shap",
            "substitution_in_model_top_n",
            named_col,
            category_col,
            note_col,
        ]
    ].reset_index(drop=True)


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

    neher_comparison = build_reference_comparison(
        site_state_importance=site_importance,
        substitution_importance=substitution_importance,
        reference_rows=NEHER2016_H3_SITE_ROWS,
        top_n=args.top_n,
        named_col="named_in_neher2016",
        category_col="neher2016_category",
        note_col="neher2016_note",
    )
    neher_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_neher_site_comparison.tsv"
    neher_comparison.to_csv(neher_path, sep="\t", index=False)

    wic_comparison = build_reference_comparison(
        site_state_importance=site_importance,
        substitution_importance=substitution_importance,
        reference_rows=WIC2023_H3_SITE_ROWS,
        top_n=args.top_n,
        named_col="named_in_wic2023",
        category_col="wic2023_category",
        note_col="wic2023_note",
    )
    wic_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_wic_paper_site_comparison.tsv"
    wic_comparison.to_csv(wic_path, sep="\t", index=False)

    shah_comparison = build_reference_comparison(
        site_state_importance=site_importance,
        substitution_importance=substitution_importance,
        reference_rows=SHAH2024_H3_SITE_ROWS,
        top_n=args.top_n,
        named_col="named_in_shah2024",
        category_col="shah2024_category",
        note_col="shah2024_note",
    )
    shah_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_shah_wic_site_comparison.tsv"
    shah_comparison.to_csv(shah_path, sep="\t", index=False)

    neher_reference = pd.DataFrame(NEHER2016_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "neher2016_category", "paper_note": "neher2016_note"}
    )
    neher_reference["named_in_neher2016"] = True
    wic_reference = pd.DataFrame(WIC2023_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "wic2023_category", "paper_note": "wic2023_note"}
    )
    wic_reference["named_in_wic2023"] = True
    shah_reference = pd.DataFrame(SHAH2024_H3_SITE_ROWS).drop_duplicates(subset=["site"]).rename(
        columns={"paper_category": "shah2024_category", "paper_note": "shah2024_note"}
    )
    shah_reference["named_in_shah2024"] = True

    combined = site_importance[["site", "rank", "mean_abs_shap"]].rename(
        columns={"rank": "site_state_rank", "mean_abs_shap": "site_state_mean_abs_shap"}
    ).merge(
        substitution_importance[["site", "rank", "mean_abs_shap"]].rename(
            columns={
                "rank": "substitution_rank",
                "mean_abs_shap": "substitution_mean_abs_shap",
            }
        ),
        on="site",
        how="left",
    ).merge(
        neher_reference,
        on="site",
        how="left",
    ).merge(
        wic_reference,
        on="site",
        how="left",
    ).merge(
        shah_reference,
        on="site",
        how="left",
    )
    combined[["named_in_neher2016", "named_in_wic2023", "named_in_shah2024"]] = combined[
        ["named_in_neher2016", "named_in_wic2023", "named_in_shah2024"]
    ].fillna(False)
    combined["named_in_any_reference"] = (
        combined["named_in_neher2016"] | combined["named_in_wic2023"] | combined["named_in_shah2024"]
    )
    combined["site_state_in_model_top_n"] = combined["site_state_rank"] <= args.top_n
    combined["substitution_in_model_top_n"] = combined["substitution_rank"].le(args.top_n).fillna(False)
    combined_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_paper_site_comparison.tsv"
    combined.to_csv(combined_path, sep="\t", index=False)

    neher_site_state_overlap = neher_comparison[
        neher_comparison["site_state_in_model_top_n"] & neher_comparison["named_in_neher2016"]
    ]
    neher_substitution_overlap = neher_comparison[
        neher_comparison["substitution_in_model_top_n"] & neher_comparison["named_in_neher2016"]
    ]
    wic_site_state_overlap = wic_comparison[
        wic_comparison["site_state_in_model_top_n"] & wic_comparison["named_in_wic2023"]
    ]
    wic_substitution_overlap = wic_comparison[
        wic_comparison["substitution_in_model_top_n"] & wic_comparison["named_in_wic2023"]
    ]
    shah_site_state_overlap = shah_comparison[
        shah_comparison["site_state_in_model_top_n"] & shah_comparison["named_in_shah2024"]
    ]
    shah_substitution_overlap = shah_comparison[
        shah_comparison["substitution_in_model_top_n"] & shah_comparison["named_in_shah2024"]
    ]
    lines = [
        f"Neher 2016 reference sites: {len(pd.DataFrame(NEHER2016_H3_SITE_ROWS).drop_duplicates(subset=['site']))}",
        f"WIC 2023 reference sites: {len(pd.DataFrame(WIC2023_H3_SITE_ROWS).drop_duplicates(subset=['site']))}",
        f"Shah 2024 reference sites: {len(pd.DataFrame(SHAH2024_H3_SITE_ROWS).drop_duplicates(subset=['site']))}",
        f"Site-state overlap with Neher in model top-{args.top_n}: {len(neher_site_state_overlap)} -> {sorted(neher_site_state_overlap['site'].tolist())}",
        f"Substitution overlap with Neher in model top-{args.top_n}: {len(neher_substitution_overlap)} -> {sorted(neher_substitution_overlap['site'].tolist())}",
        f"Site-state overlap with WIC in model top-{args.top_n}: {len(wic_site_state_overlap)} -> {sorted(wic_site_state_overlap['site'].tolist())}",
        f"Substitution overlap with WIC in model top-{args.top_n}: {len(wic_substitution_overlap)} -> {sorted(wic_substitution_overlap['site'].tolist())}",
        f"Site-state overlap with Shah/WIC in model top-{args.top_n}: {len(shah_site_state_overlap)} -> {sorted(shah_site_state_overlap['site'].tolist())}",
        f"Substitution overlap with Shah/WIC in model top-{args.top_n}: {len(shah_substitution_overlap)} -> {sorted(shah_substitution_overlap['site'].tolist())}",
        f"Output combined comparison: {combined_path}",
        f"Output Neher comparison: {neher_path}",
        f"Output WIC comparison: {wic_path}",
        f"Output Shah/WIC comparison: {shah_path}",
    ]
    append_qc_log(cfg, "30_compare_wic_paper_sites", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
