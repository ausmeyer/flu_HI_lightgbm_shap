#!/usr/bin/env python3
"""Run the full direct-comparison pipeline in sequence."""

from __future__ import annotations

import argparse
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--email", default=None, help="NCBI Entrez email for sequence fetch.")
    parser.add_argument("--api-key", default=None, help="Optional NCBI API key.")
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Skip sequence fetching and start from alignment.",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Stop after SHAP analysis.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate CSV/model artifacts instead of cleaning after the run.",
    )
    return parser.parse_args()


def run_step(cmd: list[str]) -> None:
    print("$ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()

    if not args.skip_fetch:
        fetch_cmd = [sys.executable, "scripts/01_fetch_sequences.py", "--config", args.config]
        if args.email:
            fetch_cmd.extend(["--email", args.email])
        if args.api_key:
            fetch_cmd.extend(["--api-key", args.api_key])
        run_step(fetch_cmd)

    run_step([sys.executable, "scripts/02_merge_and_align.py", "--config", args.config])
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
