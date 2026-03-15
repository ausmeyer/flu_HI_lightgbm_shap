#!/usr/bin/env python3
"""Remove intermediate pipeline artifacts and keep only final analysis outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import append_qc_log, load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    output_dir = Path(cfg["output_dir"])
    subtype = cfg["subtype"]
    keep_names = {
        f"{subtype}_HA_aligned.fasta",
        f"{subtype}_alignment_position_map.csv",
        f"{subtype}_cv_summary.csv",
        f"{subtype}_koel_validation.json",
        f"{subtype}_match_review.tsv",
        f"{subtype}_mature_position_map.csv",
        f"{subtype}_missing_from_seq_manifest.tsv",
        f"{subtype}_neher_analog_metrics.csv",
        f"{subtype}_paper_site_comparison.tsv",
        f"{subtype}_patristic_tree.nwk",
        f"{subtype}_patristic_tree_metadata.json",
        f"{subtype}_patristic_tree_tip_depths.tsv",
        f"{subtype}_present_but_failed_qc.tsv",
        f"{subtype}_qc_log.txt",
        f"{subtype}_site_stability.tsv",
        f"{subtype}_site_summary.tsv",
        f"{subtype}_substitution_site_stability.tsv",
        f"{subtype}_substitution_site_summary.tsv",
    }
    keep_figures = {
        f"{subtype}_cv_performance_summary.pdf",
        f"{subtype}_heldout_model_comparison.pdf",
        f"{subtype}_neher2016_fig2_analog.pdf",
        f"{subtype}_predicted_vs_actual.pdf",
        f"{subtype}_residualization_diagnostic.pdf",
        f"{subtype}_shap_summary_top30.pdf",
        f"{subtype}_site_component_top30.pdf",
        f"{subtype}_site_importance_bar.pdf",
        f"{subtype}_site_stability_top30.pdf",
        f"{subtype}_substitution_predicted_vs_actual.pdf",
        f"{subtype}_substitution_site_importance_bar.pdf",
        f"{subtype}_substitution_site_stability_top30.pdf",
    }

    removed = []
    kept = []
    for path in output_dir.iterdir():
        if path.name == "figures":
            for fig in path.iterdir():
                if fig.name in keep_figures:
                    kept.append(f"figures/{fig.name}")
                elif fig.is_file():
                    fig.unlink()
                    removed.append(f"figures/{fig.name}")
            continue
        if path.name in keep_names:
            kept.append(path.name)
            continue
        if path.is_file():
            path.unlink()
            removed.append(path.name)

    lines = [
        f"Kept final outputs: {len(kept)}",
        f"Removed intermediate outputs: {len(removed)}",
    ]
    append_qc_log(cfg, "10_cleanup_outputs", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
