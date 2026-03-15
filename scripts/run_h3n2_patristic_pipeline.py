#!/usr/bin/env python3
"""Run the H3N2 pipeline with patristic distance added alongside temporal distance."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from common import load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2_patristic.json")
    parser.add_argument("--email", default=None)
    parser.add_argument("--api-key", default=None)
    parser.add_argument("--skip-fetch", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--keep-intermediates", action="store_true")
    return parser.parse_args()


def run_step(cmd: list[str]) -> None:
    print("$ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def bootstrap_alignment_artifacts(cfg: dict) -> None:
    source_dir = cfg.get("bootstrap_alignment_from_output_dir")
    if not source_dir:
        return

    source_dir = Path(str(source_dir))
    target_dir = Path(cfg["output_dir"])
    target_dir.mkdir(parents=True, exist_ok=True)

    source_subtype = str(cfg.get("bootstrap_alignment_source_subtype", "H3N2"))
    target_subtype = str(cfg["subtype"])
    artifact_suffixes = [
        "HA_aligned.fasta",
        "alignment_position_map.csv",
        "mature_position_map.csv",
    ]
    for suffix in artifact_suffixes:
        source_path = source_dir / f"{source_subtype}_{suffix}"
        if not source_path.exists():
            raise FileNotFoundError(f"Bootstrap alignment artifact not found: {source_path}")
        target_path = target_dir / f"{target_subtype}_{suffix}"
        target_path.write_bytes(source_path.read_bytes())


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    if not args.skip_fetch:
        fetch_cmd = [sys.executable, "scripts/01_fetch_sequences.py", "--config", args.config]
        if args.email:
            fetch_cmd.extend(["--email", args.email])
        if args.api_key:
            fetch_cmd.extend(["--api-key", args.api_key])
        run_step(fetch_cmd)

    bootstrap_alignment_artifacts(cfg)
    run_step([sys.executable, "scripts/33_build_h3n2_patristic_tree.py", "--config", args.config])
    run_step([sys.executable, "scripts/03_match_strains.py", "--config", args.config])
    run_step([sys.executable, "scripts/04_preprocess_titers.py", "--config", args.config])
    run_step([sys.executable, "scripts/05_build_features.py", "--config", args.config])
    run_step([sys.executable, "scripts/06_train_model.py", "--config", args.config])
    run_step([sys.executable, "scripts/07_shap_analysis.py", "--config", args.config])
    if not args.skip_figures:
        run_step([sys.executable, "scripts/08_generate_figures.py", "--config", args.config])
    run_step([sys.executable, "scripts/11_compare_paper_sites.py", "--config", args.config])
    run_step([sys.executable, "scripts/09_audit_strain_coverage.py", "--config", args.config])
    if not args.keep_intermediates:
        run_step([sys.executable, "scripts/10_cleanup_outputs.py", "--config", args.config])


if __name__ == "__main__":
    main()
