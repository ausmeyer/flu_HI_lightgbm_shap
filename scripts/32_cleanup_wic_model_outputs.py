#!/usr/bin/env python3
"""Remove regenerable WIC modeling intermediates and keep publication-oriented outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import append_qc_log, load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    output_dir = Path(cfg["output_dir"])
    subtype = cfg["subtype"]

    keep_names = {
        f"{subtype}_HA1_aligned.fasta",
        f"{subtype}_cv_summary.tsv",
        f"{subtype}_ha1_alignment_metadata.json",
        f"{subtype}_ha1_extraction_summary.tsv",
        f"{subtype}_ha1_position_map.csv",
        f"{subtype}_matched_titers.tsv",
        f"{subtype}_neher_analog_metrics.tsv",
        f"{subtype}_paper_site_comparison.tsv",
        f"{subtype}_paper_site_validation.json",
        f"{subtype}_patristic_tree.nwk",
        f"{subtype}_patristic_tree_metadata.json",
        f"{subtype}_patristic_tree_tip_depths.tsv",
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
                    if not args.dry_run:
                        fig.unlink()
                    removed.append(f"figures/{fig.name}")
            continue
        if path.name in keep_names:
            kept.append(path.name)
            continue
        if path.is_file():
            if not args.dry_run:
                path.unlink()
            removed.append(path.name)

    lines = [
        "Retention policy: publication-oriented WIC modeling outputs",
        f"Kept outputs: {len(kept)}",
        f"Removed intermediates: {len(removed)}",
    ]
    if removed:
        lines.append(f"Removed names: {', '.join(sorted(removed))}")
    append_qc_log(cfg, "32_cleanup_wic_model_outputs", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
