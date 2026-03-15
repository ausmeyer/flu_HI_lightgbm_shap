#!/usr/bin/env python3
"""Compare ranked sites from the local model with H3 sites named in Neher et al. 2016."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import append_qc_log, load_config


PAPER_SITE_ROWS = [
    {"site": 62, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
    {"site": 121, "paper_category": "named_substitution", "paper_note": "FU02 cluster change attributed partly to N121T"},
    {"site": 135, "paper_category": "named_substitution", "paper_note": "Repeated effects include K135E; SI87 to BE89 includes G135N"},
    {"site": 140, "paper_category": "named_substitution", "paper_note": "Repeated effects include K140E"},
    {"site": 144, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
    {"site": 145, "paper_category": "koel7_and_cluster_transition", "paper_note": "Koel 7 site; SI87 to BE89 involves N145K"},
    {"site": 155, "paper_category": "koel7", "paper_note": "Koel 7 site listed in paper"},
    {"site": 156, "paper_category": "koel7_and_cluster_transition", "paper_note": "Koel 7 site; WU95 to SY97 set includes K156Q; FU02 change includes Q156H"},
    {"site": 158, "paper_category": "koel7_and_named_substitution", "paper_note": "Koel 7 site; repeated effects include K158R; WU95 to SY97 set includes E158K"},
    {"site": 159, "paper_category": "koel7_and_named_substitution", "paper_note": "Koel 7 site; repeated effects include Y159F; later text discusses position 159"},
    {"site": 186, "paper_category": "cluster_transition", "paper_note": "SI87 to BE89 includes I186S"},
    {"site": 189, "paper_category": "koel7_and_named_substitution", "paper_note": "Koel 7 site; repeated effects include K189N; text notes S189N can be small"},
    {"site": 193, "paper_category": "koel7_and_cluster_transition", "paper_note": "Koel 7 site; SI87 to BE89 includes N193S"},
    {"site": 196, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
    {"site": 276, "paper_category": "cluster_transition", "paper_note": "WU95 to SY97 set K62E/V144I/K156Q/E158K/V196A/N276K"},
]


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
    site_importance = pd.read_csv(site_path)
    substitution_importance = pd.read_csv(substitution_path)
    site_importance["rank"] = site_importance["rank"].astype(int)
    substitution_importance["rank"] = substitution_importance["rank"].astype(int)

    paper_sites = pd.DataFrame(PAPER_SITE_ROWS).drop_duplicates(subset=["site"]).sort_values("site")

    comparison = (
        paper_sites.assign(named_in_paper=True)
        .merge(
            site_importance[["site", "rank", "mean_abs_shap", "is_koel_site"]].rename(
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
    comparison["named_in_paper"] = comparison["named_in_paper"].fillna(False)
    comparison["is_koel_site"] = comparison["is_koel_site"].fillna(False)
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
            "named_in_paper",
            "is_koel_site",
            "paper_category",
            "paper_note",
        ]
    ]

    out_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_paper_site_comparison.tsv"
    comparison.to_csv(out_path, sep="\t", index=False)

    site_state_overlap = comparison[
        comparison["site_state_in_model_top_n"] & comparison["named_in_paper"]
    ]
    substitution_overlap = comparison[
        comparison["substitution_in_model_top_n"] & comparison["named_in_paper"]
    ]
    lines = [
        f"Paper-named sites in reference table: {paper_sites['site'].nunique()}",
        f"Site-state overlap in model top-{args.top_n}: {len(site_state_overlap)} -> {sorted(site_state_overlap['site'].tolist())}",
        f"Substitution overlap in model top-{args.top_n}: {len(substitution_overlap)} -> {sorted(substitution_overlap['site'].tolist())}",
        f"Output: {out_path}",
    ]
    append_qc_log(cfg, "11_compare_paper_sites", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
