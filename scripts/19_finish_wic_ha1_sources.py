#!/usr/bin/env python3
"""Finish WIC HA1 source preparation after local GISAID FASTA is available."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        default="H3N2-WIC/config/h3n2_wic_config.json",
        help="Path to WIC config JSON.",
    )
    parser.add_argument("--email", required=True, help="NCBI Entrez email.")
    parser.add_argument("--api-key", default=None, help="Optional NCBI API key.")
    parser.add_argument(
        "--max-strains",
        type=int,
        default=None,
        help="Optional limit for the public search step.",
    )
    parser.add_argument(
        "--download-repo-fasta",
        action="store_true",
        help="Force fresh download of the published repo HA1 FASTA.",
    )
    parser.add_argument(
        "--refresh-public-cache",
        action="store_true",
        help="Ignore any existing cached public lookup results.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate manifests after a successful run.",
    )
    return parser.parse_args()


def run_step(cmd: list[str], env: dict[str, str]) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True, env=env)


def cleanup_intermediate_outputs() -> None:
    files_to_remove = [
        "H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_autofill.tsv",
        "H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_repofill.tsv",
        "H3N2-WIC/manifests/H3N2_WIC_seq_manifest_gisaidfill.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_unresolved_strains.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_repo_unmatched.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_public_unresolved.tsv",
    ]
    for rel_path in files_to_remove:
        path = Path(rel_path)
        if path.exists():
            path.unlink()


def main() -> None:
    args = parse_args()
    env = os.environ.copy()
    with open(args.config) as fh:
        cfg = json.load(fh)

    run_step([sys.executable, "scripts/13_autofill_wic_manifest.py"], env)

    repo_cmd = [sys.executable, "scripts/16_match_wic_repo_ha1.py"]
    if args.download_repo_fasta:
        repo_cmd.append("--download")
    run_step(repo_cmd, env)

    run_step(
        [
            sys.executable,
            "scripts/18_match_wic_local_gisaid_ha1.py",
            "--gisaid-fasta",
            cfg["gisaid_fasta_path"],
        ],
        env,
    )

    public_cmd = [
        sys.executable,
        "scripts/15_resolve_wic_public_accessions.py",
        "--config",
        args.config,
        "--manifest",
        "H3N2-WIC/manifests/H3N2_WIC_seq_manifest_gisaidfill.tsv",
        "--email",
        args.email,
    ]
    if args.api_key:
        public_cmd.extend(["--api-key", args.api_key])
    if args.max_strains is not None:
        public_cmd.extend(["--max-strains", str(args.max_strains)])
    if args.refresh_public_cache:
        public_cmd.append("--refresh-cache")
    run_step(public_cmd, env)

    if not args.keep_intermediates:
        cleanup_intermediate_outputs()


if __name__ == "__main__":
    main()
