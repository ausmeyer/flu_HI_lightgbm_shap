#!/usr/bin/env python3
"""Copy manuscript figures into the local Overleaf manuscript repository."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANUSCRIPT_DIR = ROOT / "69b778fd3e7b181fe1c2943b"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manuscript-dir", default=str(DEFAULT_MANUSCRIPT_DIR))
    return parser.parse_args()


def copy_file(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Required manuscript asset not found: {src}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def main() -> None:
    args = parse_args()
    manuscript_dir = Path(args.manuscript_dir)
    figures_dir = manuscript_dir / "figures"

    mapping = {
        ROOT / "H3N2" / "output" / "figures" / "H3N2_site_component_top30.pdf":
            figures_dir / "fig1_h3n2_site_component_top30.pdf",
        ROOT / "H3N2-WIC-no-egg-no-mixed-no-unknown" / "output" / "modeling" / "figures"
        / "H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_site_component_top30.pdf":
            figures_dir / "fig2_wic_filtered_site_component_top30.pdf",
        ROOT / "manuscript" / "generated_figures" / "fig4_cross_study_concordance.pdf":
            figures_dir / "fig3_cross_study_concordance.pdf",
    }

    for src, dst in mapping.items():
        copy_file(src, dst)


if __name__ == "__main__":
    main()
