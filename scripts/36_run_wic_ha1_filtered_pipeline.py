#!/usr/bin/env python3
"""Run the WIC HA1 pipeline after excluding egg, mixed, and unknown passage classes."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from common import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown.json")
    parser.add_argument("--n-splits", type=int, default=None)
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser.parse_args()


def run_step(cmd: list[str]) -> None:
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def bootstrap_alignment_artifacts(cfg: dict) -> None:
    source_dir = cfg.get("bootstrap_alignment_from_output_dir")
    if not source_dir:
        return

    source_dir = Path(str(source_dir))
    target_dir = Path(cfg["output_dir"])
    target_dir.mkdir(parents=True, exist_ok=True)

    source_subtype = str(cfg.get("bootstrap_alignment_source_subtype", "H3N2_WIC_HA1"))
    target_subtype = str(cfg["subtype"])
    artifact_suffixes = [
        "HA1_aligned.fasta",
        "ha1_alignment_metadata.json",
        "ha1_extraction_summary.tsv",
        "ha1_position_map.csv",
    ]
    for suffix in artifact_suffixes:
        source_path = source_dir / f"{source_subtype}_{suffix}"
        if not source_path.exists():
            raise FileNotFoundError(f"Bootstrap alignment artifact not found: {source_path}")
        target_path = target_dir / f"{target_subtype}_{suffix}"
        target_path.write_bytes(source_path.read_bytes())


def main() -> None:
    args = parse_args()
    config = args.config
    cfg = load_config(config)

    bootstrap_alignment_artifacts(cfg)
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
