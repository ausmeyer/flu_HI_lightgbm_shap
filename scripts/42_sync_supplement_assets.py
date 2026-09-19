#!/usr/bin/env python3
"""Copy supplementary figures and generated tables into the local SI Overleaf repository."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SUPPLEMENT_DIR = ROOT / "69b7797c762f515edcff3ad6"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement-dir", default=str(DEFAULT_SUPPLEMENT_DIR))
    return parser.parse_args()


def copy_file(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Required supplement asset not found: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()
    supplement_dir = Path(args.supplement_dir)
    figures_dir = supplement_dir / "figures"
    tables_dir = supplement_dir / "tables"

    figure_mapping = {
        ROOT / "H3N2" / "output" / "figures" / "H3N2_heldout_model_comparison.pdf":
            figures_dir / "fig_s1_h3n2_heldout_model_comparison.pdf",
        ROOT / "H3N2" / "output" / "figures" / "H3N2_predicted_vs_actual.pdf":
            figures_dir / "fig_s2_h3n2_predicted_vs_actual.pdf",
        ROOT / "H3N2" / "output" / "figures" / "H3N2_site_stability_top30.pdf":
            figures_dir / "fig_s3_h3n2_site_stability_top30.pdf",
        ROOT / "H3N2" / "output" / "figures" / "H3N2_neher2016_fig2_analog.pdf":
            figures_dir / "fig_s4_h3n2_neher_analog.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_heldout_model_comparison.pdf":
            figures_dir / "fig_s5_wic_filtered_heldout_model_comparison.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_predicted_vs_actual.pdf":
            figures_dir / "fig_s6_wic_filtered_predicted_vs_actual.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_site_stability_top30.pdf":
            figures_dir / "fig_s7_wic_filtered_site_stability_top30.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_residualization_diagnostic.pdf":
            figures_dir / "fig_s8_wic_filtered_residualization.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_neher2016_fig2_analog.pdf":
            figures_dir / "fig_s9_wic_filtered_neher_analog.pdf",
        ROOT / "H3N2-WIC" / "output" / "modeling" / "figures" / "H3N2_WIC_HA1_predicted_vs_actual.pdf":
            figures_dir / "fig_s10_wic_full_predicted_vs_actual.pdf",
        ROOT / "H3N2-WIC" / "output" / "modeling" / "figures" / "H3N2_WIC_HA1_site_component_top30.pdf":
            figures_dir / "fig_s11_wic_full_site_component_top30.pdf",
        ROOT / "H3N2-patristic" / "output" / "figures" / "H3N2_PATRISTIC_site_component_top30.pdf":
            figures_dir / "fig_s12_h3n2_patristic_site_component_top30.pdf",
        ROOT / "H3N2-patristic" / "output" / "figures" / "H3N2_PATRISTIC_heldout_model_comparison.pdf":
            figures_dir / "fig_s13_h3n2_patristic_heldout_model_comparison.pdf",
        ROOT / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_site_component_top30.pdf":
            figures_dir / "fig_s14_wic_filtered_patristic_site_component_top30.pdf",
        ROOT / "H3N2-WIC-patristic-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_PATRISTIC_heldout_model_comparison.pdf":
            figures_dir / "fig_s15_wic_filtered_patristic_heldout_model_comparison.pdf",
        ROOT / "H3N2" / "output" / "figures" / "H3N2_signed_feature_shap_top20.pdf":
            figures_dir / "fig_s16_h3n2_signed_feature_shap.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_signed_feature_shap_top20.pdf":
            figures_dir / "fig_s17_wic_signed_feature_shap.pdf",
        ROOT / "H3N2" / "output" / "figures" / "H3N2_koel_context_dependence_145_155.pdf":
            figures_dir / "fig_s18_koel_context_dependence.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_glycan_context_dependence_158_160.pdf":
            figures_dir / "fig_s19_glycan_context_dependence.pdf",
    }

    table_mapping = {
        ROOT / "manuscript" / "generated_supplement" / "table_s1_model_performance.tex":
            tables_dir / "table_s1_model_performance.tex",
        ROOT / "manuscript" / "generated_supplement" / "table_s1_model_performance.tsv":
            tables_dir / "table_s1_model_performance.tsv",
        ROOT / "manuscript" / "generated_supplement" / "table_s2_reference_overlap.tex":
            tables_dir / "table_s2_reference_overlap.tex",
        ROOT / "manuscript" / "generated_supplement" / "table_s2_reference_overlap.tsv":
            tables_dir / "table_s2_reference_overlap.tsv",
        ROOT / "manuscript" / "generated_supplement" / "table_s4_neher_bedford_shap_values.tsv":
            tables_dir / "table_s4_neher_bedford_shap_values.tsv",
        ROOT / "manuscript" / "generated_supplement" / "table_s5_wic_filtered_shap_values.tsv":
            tables_dir / "table_s5_wic_filtered_shap_values.tsv",
        ROOT / "manuscript" / "generated_supplement" / "table_s6_literature_overlap.tex":
            tables_dir / "table_s6_literature_overlap.tex",
        ROOT / "manuscript" / "generated_supplement" / "table_s6_literature_overlap.tsv":
            tables_dir / "table_s6_literature_overlap.tsv",
    }

    for src, dst in {**figure_mapping, **table_mapping}.items():
        copy_file(src, dst)


if __name__ == "__main__":
    main()
