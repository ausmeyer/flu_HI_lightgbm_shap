#!/usr/bin/env python3
"""Build a FastTree patristic-distance tree for the aligned WIC HA1 sequences."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import append_qc_log, build_patristic_tree, load_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_patristic_model.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    aligned_fasta = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA1_aligned.fasta"
    if not aligned_fasta.exists():
        raise FileNotFoundError(f"Aligned HA1 FASTA not found: {aligned_fasta}")

    tree_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_patristic_tree.nwk"
    tip_depths_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_patristic_tree_tip_depths.tsv"
    metadata_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_patristic_tree_metadata.json"

    metadata = build_patristic_tree(
        aligned_fasta_path=str(aligned_fasta),
        tree_out_path=str(tree_path),
        tip_depths_out_path=str(tip_depths_path),
        metadata_out_path=str(metadata_path),
    )

    lines = [
        f"Aligned HA1 FASTA: {aligned_fasta}",
        f"Patristic tree output: {tree_path}",
        f"Tip-depth summary: {tip_depths_path}",
        f"Tree metadata: {metadata_path}",
        f"FastTree model: {metadata['fasttree_model']}",
        f"Leaf count: {metadata['n_tips']}",
        f"Branches with lengths: {metadata['n_branches_with_lengths']}",
        f"Total branch length: {metadata['sum_branch_length']:.6f}",
    ]
    append_qc_log(cfg, "34_build_wic_ha1_patristic_tree", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
