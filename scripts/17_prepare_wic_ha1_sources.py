#!/usr/bin/env python3
"""Prepare WIC HA1 sources through the repo-match stage and write the GISAID search batches."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--email", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--api-key", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--max-strains", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--refresh-public-cache", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument(
        "--download-repo-fasta",
        action="store_true",
        help="Force fresh download of the published repo HA1 FASTA.",
    )
    parser.add_argument(
        "--restart-clean",
        action="store_true",
        help="Delete generated WIC HA1 source-prep outputs and restart from scratch.",
    )
    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate WIC source-prep files instead of pruning them after a successful run.",
    )
    return parser.parse_args()


def run_step(cmd: list[str], env: dict[str, str]) -> None:
    print("$", " ".join(cmd))
    subprocess.run(cmd, check=True, env=env)


def remove_empty_gisaid_batch_dir() -> None:
    batch_dir = Path("H3N2-WIC/metadata/gisaid_batches")
    if batch_dir.exists() and not any(batch_dir.iterdir()):
        batch_dir.rmdir()


def clean_generated_outputs() -> None:
    files_to_remove = [
        "H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_autofill.tsv",
        "H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_repofill.tsv",
        "H3N2-WIC/manifests/H3N2_WIC_seq_manifest_gisaidfill.tsv",
        "H3N2-WIC/manifests/H3N2_WIC_seq_manifest_publicfill.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_unresolved_strains.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_repo_unmatched.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_public_unresolved.tsv",
        "H3N2-WIC/metadata/H3N2_WIC_gisaid_batch_query.txt",
        "H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.txt",
        "H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.csv",
        "H3N2-WIC/metadata/H3N2_WIC_repo_match_summary.tsv",
        "H3N2-WIC/data/intermediate/H3N2_WIC_repo_matched_ha1.fasta",
        "H3N2-WIC/metadata/H3N2_WIC_gisaid_match_summary.tsv",
        "H3N2-WIC/data/intermediate/H3N2_WIC_gisaid_matched_ha1.fasta",
        "H3N2-WIC/data/raw/repo/h3n2_ha1_aa_repo.fasta",
        "H3N2-WIC/output/H3N2_WIC_public_resolution_cache.tsv",
    ]
    glob_patterns = [
        "H3N2-WIC/metadata/gisaid_batches/wic_gisaid_batch_*.txt",
    ]

    for rel_path in files_to_remove:
        path = Path(rel_path)
        if path.exists():
            path.unlink()

    for pattern in glob_patterns:
        base = Path(pattern.split("*", 1)[0]).parent
        if base.exists():
            for path in base.glob(Path(pattern).name):
                path.unlink()

    batch_manifest = Path("H3N2-WIC/metadata/gisaid_batches/wic_gisaid_batch_manifest.tsv")
    if batch_manifest.exists():
        batch_manifest.unlink()
    remove_empty_gisaid_batch_dir()


def cleanup_intermediate_outputs() -> None:
    files_to_remove = [
        "H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_autofill.tsv",
        "H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_repofill.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_unresolved_strains.tsv",
        "H3N2-WIC/metadata/intermediate/H3N2_WIC_repo_unmatched.tsv",
        "H3N2-WIC/metadata/H3N2_WIC_gisaid_batch_query.txt",
        "H3N2-WIC/output/H3N2_WIC_public_resolution_cache.tsv",
    ]
    removed = []
    for rel_path in files_to_remove:
        path = Path(rel_path)
        if path.exists():
            path.unlink()
            removed.append(rel_path)
    batch_dir = Path("H3N2-WIC/metadata/gisaid_batches")
    if batch_dir.exists():
        for path in batch_dir.glob("wic_gisaid_batch_*.txt"):
            path.unlink()
            removed.append(str(path))
        batch_manifest = batch_dir / "wic_gisaid_batch_manifest.tsv"
        if batch_manifest.exists():
            batch_manifest.unlink()
            removed.append(str(batch_manifest))
    remove_empty_gisaid_batch_dir()
    if removed:
        print(f"Pruned intermediate WIC source-prep files: {len(removed)}")
    else:
        print("No intermediate WIC source-prep files needed pruning.")


def main() -> None:
    args = parse_args()
    env = os.environ.copy()

    if args.restart_clean:
        clean_generated_outputs()

    run_step([sys.executable, "scripts/13_autofill_wic_manifest.py"], env)

    repo_cmd = [sys.executable, "scripts/16_match_wic_repo_ha1.py"]
    if args.download_repo_fasta:
        repo_cmd.append("--download")
    run_step(repo_cmd, env)

    run_step(
        [
            sys.executable,
            "scripts/14_prepare_wic_gisaid_batches.py",
            "--input",
            "H3N2-WIC/metadata/intermediate/H3N2_WIC_repo_unmatched.tsv",
        ],
        env,
    )

    if not args.keep_intermediates:
        cleanup_intermediate_outputs()


if __name__ == "__main__":
    main()
