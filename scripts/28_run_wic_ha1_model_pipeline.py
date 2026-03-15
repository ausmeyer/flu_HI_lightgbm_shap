#!/usr/bin/env python3
"""Run the WIC HA1 LightGBM modeling pipeline end to end."""

from __future__ import annotations

import argparse
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    parser.add_argument("--n-splits", type=int, default=None)
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep large regenerable intermediates instead of cleaning to publication-oriented outputs.",
    )
    return parser.parse_args()


def run_step(cmd: list[str]) -> None:
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    config = args.config

    run_step([sys.executable, "scripts/24_prepare_wic_ha1_alignment.py", "--config", config])
    run_step([sys.executable, "scripts/25_preprocess_wic_model_titers.py", "--config", config])
    run_step([sys.executable, "scripts/26_build_wic_model_features.py", "--config", config])

    train_cmd = [sys.executable, "scripts/27_train_wic_model.py", "--config", config]
    if args.n_splits is not None:
        train_cmd.extend(["--n-splits", str(args.n_splits)])
    run_step(train_cmd)
    run_step([sys.executable, "scripts/29_shap_analysis_wic.py", "--config", config])
    run_step([sys.executable, "scripts/30_compare_wic_paper_sites.py", "--config", config])
    run_step([sys.executable, "scripts/31_generate_wic_figures.py", "--config", config])
    if not args.keep_intermediates:
        run_step([sys.executable, "scripts/32_cleanup_wic_model_outputs.py", "--config", config])


if __name__ == "__main__":
    main()
