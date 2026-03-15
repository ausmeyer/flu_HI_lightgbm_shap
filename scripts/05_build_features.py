#!/usr/bin/env python3
"""Build modeling feature matrix from matched titers and aligned HA sequences."""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from common import (
    AMBIGUOUS_AA,
    append_qc_log,
    compute_patristic_distances,
    get_distance_feature_cols,
    load_config,
    read_fasta,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    aligned_path = f"{cfg['output_dir']}/{cfg['subtype']}_HA_aligned.fasta"
    titers_path = f"{cfg['output_dir']}/{cfg['subtype']}_matched_titers.csv"
    position_map_path = f"{cfg['output_dir']}/{cfg['subtype']}_mature_position_map.csv"

    seq_map = read_fasta(aligned_path)
    df = pd.read_csv(titers_path)
    position_map = pd.read_csv(position_map_path)

    ref_name = cfg["reference_strain"]
    if ref_name not in seq_map:
        raise ValueError(f"Reference strain not found in aligned FASTA: {ref_name}")

    ref_seq = seq_map[ref_name]
    ref_cols = [i for i, aa in enumerate(ref_seq) if aa != "-"]
    mature_len = int(cfg.get("mature_protein_length", 550))
    if len(ref_cols) < mature_len:
        raise ValueError(
            f"Reference has {len(ref_cols)} non-gap residues after trimming; expected >= {mature_len}."
        )
    pos_cols = ref_cols[:mature_len]

    strains = sorted(seq_map.keys())
    strain_to_idx = {s: i for i, s in enumerate(strains)}

    missing_virus = df[~df["virus_strain_matched"].isin(strains)]
    missing_serum = df[~df["serum_strain_matched"].isin(strains)]
    if len(missing_virus) or len(missing_serum):
        raise ValueError(
            "Matched titers contain strain IDs absent from alignment. "
            f"Missing virus rows={len(missing_virus)}, missing serum rows={len(missing_serum)}"
        )

    seq_array = np.array([list(seq_map[s]) for s in strains], dtype="U1")
    seq_pos = seq_array[:, pos_cols]

    virus_idx = df["virus_strain_matched"].map(strain_to_idx).to_numpy(dtype=int)
    serum_idx = df["serum_strain_matched"].map(strain_to_idx).to_numpy(dtype=int)

    virus_chars = seq_pos[virus_idx]
    serum_chars = seq_pos[serum_idx]

    amb = np.array(sorted(AMBIGUOUS_AA), dtype="U1")
    missing_mask = (
        (virus_chars == "-")
        | (serum_chars == "-")
        | np.isin(virus_chars, amb)
        | np.isin(serum_chars, amb)
    )

    changes = (virus_chars != serum_chars).astype(float)
    changes[missing_mask] = np.nan

    change_cols = [f"change_{i}" for i in range(1, mature_len + 1)]
    change_df = pd.DataFrame(changes, columns=change_cols)
    canonical_aa = sorted(set("ACDEFGHIKLMNPQRSTVWY"))

    virus_aa_parts = []
    virus_covariate_names = []
    serum_aa_parts = []
    serum_covariate_names = []
    for site_idx in range(mature_len):
        site = site_idx + 1
        virus_site = virus_chars[:, site_idx]
        serum_site = serum_chars[:, site_idx]

        observed_virus = sorted(set(virus_site) & set(canonical_aa))
        for aa in observed_virus:
            col = f"{aa}{site}"
            virus_covariate_names.append(col)
            virus_aa_parts.append((virus_site == aa).astype(np.uint8))

        observed_serum = sorted(set(serum_site) & set(canonical_aa))
        for aa in observed_serum:
            col = f"{site}{aa}"
            serum_covariate_names.append(col)
            serum_aa_parts.append((serum_site == aa).astype(np.uint8))

    virus_aa_df = pd.DataFrame(
        np.column_stack(virus_aa_parts) if virus_aa_parts else np.empty((len(df), 0), dtype=np.uint8),
        columns=virus_covariate_names,
    )
    serum_aa_df = pd.DataFrame(
        np.column_stack(serum_aa_parts) if serum_aa_parts else np.empty((len(df), 0), dtype=np.uint8),
        columns=serum_covariate_names,
    )

    substitution_coords = np.where((changes == 1.0) & ~missing_mask)
    substitution_labels = [
        f"sub_{pos + 1}_{serum_chars[row, pos]}_{virus_chars[row, pos]}"
        for row, pos in zip(*substitution_coords)
    ]
    unique_substitutions = sorted(set(substitution_labels))
    substitution_df = pd.DataFrame(index=np.arange(len(df)))
    if unique_substitutions:
        sub_index = {label: idx for idx, label in enumerate(unique_substitutions)}
        sub_matrix = np.zeros((len(df), len(unique_substitutions)), dtype=np.uint8)
        for row, label in zip(substitution_coords[0], substitution_labels):
            sub_matrix[row, sub_index[label]] = 1
        substitution_df = pd.DataFrame(sub_matrix, columns=unique_substitutions)

    df["virusYear"] = pd.to_numeric(df["virusYear"], errors="coerce")
    df["serumYear"] = pd.to_numeric(df["serumYear"], errors="coerce")
    df["temporal_distance"] = (df["virusYear"] - df["serumYear"]).abs()
    distance_cols = get_distance_feature_cols(cfg)
    if "patristic_distance" in distance_cols:
        tree_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_patristic_tree.nwk"
        if not tree_path.exists():
            raise FileNotFoundError(f"Patristic tree not found for configured distance feature: {tree_path}")
        df["patristic_distance"] = compute_patristic_distances(
            tree_path=str(tree_path),
            left_ids=df["virus_strain_matched"],
            right_ids=df["serum_strain_matched"],
        ).to_numpy(dtype=float)

    meta_cols = [
        "virusIsolate",
        "virusStrain",
        "virus_strain_matched",
        "virusYear",
        "serumIsolate",
        "serumStrain",
        "serum_strain_matched",
        "serumYear",
        "source",
        "lab",
        "titer_numeric",
        "log2_titer",
        "homologous_log2_titer",
        "homologous_titer_source",
        "standardized_titer",
        "serum_potency",
        "virus_avidity",
        "additive_intercept",
        "additive_fitted",
        "corrected_titer",
        "temporal_distance",
    ]
    if "patristic_distance" in df.columns:
        meta_cols.append("patristic_distance")

    feature_matrix = pd.concat(
        [df[meta_cols].reset_index(drop=True), change_df, virus_aa_df, serum_aa_df],
        axis=1,
    )

    out_path = f"{cfg['output_dir']}/{cfg['subtype']}_feature_matrix.csv"
    feature_matrix.to_csv(out_path, index=False)

    substitution_feature_matrix = pd.concat(
        [df[meta_cols].reset_index(drop=True), substitution_df.reset_index(drop=True)],
        axis=1,
    )
    substitution_out_path = f"{cfg['output_dir']}/{cfg['subtype']}_substitution_feature_matrix.csv"
    substitution_feature_matrix.to_csv(substitution_out_path, index=False)

    observed_sites = (
        pd.DataFrame(
            {
                "mature_position": np.arange(1, mature_len + 1, dtype=int),
                "n_changed_measurements": np.nansum(changes == 1.0, axis=0).astype(int),
                "n_missing_measurements": np.isnan(changes).sum(axis=0).astype(int),
            }
        )
        .merge(position_map, on="mature_position", how="left")
    )
    observed_sites_path = f"{cfg['output_dir']}/{cfg['subtype']}_observed_site_changes.csv"
    observed_sites.to_csv(observed_sites_path, index=False)

    missing_ratio = float(np.isnan(changes).sum() / changes.size)

    lines = [
        f"Input rows: {len(df)}",
        f"Aligned strains available: {len(strains)}",
        f"Reference strain: {ref_name}",
        f"Mature positions used: {mature_len}",
        f"Feature matrix shape: {feature_matrix.shape}",
        f"Primary site-state covariates: {len(change_cols) + len(virus_covariate_names) + len(serum_covariate_names)}",
        f"Distance feature columns: {distance_cols}",
        f"Substitution feature matrix shape: {substitution_feature_matrix.shape}",
        f"Observed substitution features: {len(unique_substitutions)}",
        f"Change-matrix missing fraction (gaps/ambiguous): {missing_ratio:.6f}",
        f"Output: {out_path}",
        f"Substitution-feature output: {substitution_out_path}",
        f"Observed-site summary: {observed_sites_path}",
    ]
    append_qc_log(cfg, "05_build_features", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
