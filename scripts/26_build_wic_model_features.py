#!/usr/bin/env python3
"""Build HA1 feature matrices for WIC LightGBM modeling."""

from __future__ import annotations

import argparse
import json
import re
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
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    titers_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_matched_titers.tsv"
    aligned_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA1_aligned.fasta"
    metadata_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_ha1_alignment_metadata.json"

    df = pd.read_csv(titers_path, sep="\t")
    seq_map = read_fasta(str(aligned_path))
    with open(metadata_path, "r", encoding="utf-8") as f:
        alignment_meta = json.load(f)

    ref_name = str(alignment_meta["reference_sequence_id"])
    if ref_name not in seq_map:
        raise ValueError(f"Reference sequence id not found in aligned WIC HA1 FASTA: {ref_name}")

    ref_seq = seq_map[ref_name]
    ref_cols = [i for i, aa in enumerate(ref_seq) if aa != "-"]
    ha1_len = int(cfg.get("ha1_length", 328))
    if len(ref_cols) < ha1_len:
        raise ValueError(
            f"Reference HA1 has {len(ref_cols)} non-gap residues after alignment; expected >= {ha1_len}."
        )
    pos_cols = ref_cols[:ha1_len]

    strains = sorted(seq_map.keys())
    strain_to_idx = {s: i for i, s in enumerate(strains)}

    has_virus = df["virusSequenceId"].isin(strains)
    has_serum = df["serumSequenceId"].isin(strains)
    dropped_rows = int((~(has_virus & has_serum)).sum())
    df = df[has_virus & has_serum].copy()

    seq_array = np.array([list(seq_map[s]) for s in strains], dtype="U1")
    seq_pos = seq_array[:, pos_cols]

    virus_idx = df["virusSequenceId"].map(strain_to_idx).to_numpy(dtype=int)
    serum_idx = df["serumSequenceId"].map(strain_to_idx).to_numpy(dtype=int)
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

    change_cols = [f"change_{i}" for i in range(1, ha1_len + 1)]
    change_df = pd.DataFrame(changes, columns=change_cols)

    canonical_aa = sorted(set("ACDEFGHIKLMNPQRSTVWY"))
    virus_aa_parts = []
    virus_covariate_names = []
    serum_aa_parts = []
    serum_covariate_names = []
    for site_idx in range(ha1_len):
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

    meta_cols = [
        "virusIsolate",
        "virusStrain",
        "virusSequenceId",
        "serumIsolate",
        "serumStrain",
        "serumSequenceId",
        "source",
        "titer_numeric",
        "log2_titer",
        "homologous_log2_titer",
        "homologous_titer_source",
        "standardized_titer",
        "context_additive_intercept",
        "context_additive_fitted",
        "context_corrected_titer",
        "virus_passage_class",
        "reference_passage_class",
        "virus_passage_known",
        "reference_passage_known",
        "antiserum_class",
        "rbc_type",
        "rbc_species",
        "assay_date_label",
        "oseltamivir_present",
        "oseltamivir_nm",
        "temporal_distance",
    ]
    distance_cols = get_distance_feature_cols(cfg)
    if "patristic_distance" in distance_cols:
        tree_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_patristic_tree.nwk"
        if not tree_path.exists():
            raise FileNotFoundError(f"Patristic tree not found for configured distance feature: {tree_path}")
        df["patristic_distance"] = compute_patristic_distances(
            tree_path=str(tree_path),
            left_ids=df["virusSequenceId"],
            right_ids=df["serumSequenceId"],
        ).to_numpy(dtype=float)
        meta_cols.append("patristic_distance")

    feature_matrix = pd.concat(
        [df[meta_cols].reset_index(drop=True), change_df, virus_aa_df, serum_aa_df],
        axis=1,
    )
    feature_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_feature_matrix.tsv"
    feature_matrix.to_csv(feature_out, sep="\t", index=False)

    substitution_feature_matrix = pd.concat(
        [df[meta_cols].reset_index(drop=True), substitution_df.reset_index(drop=True)],
        axis=1,
    )
    substitution_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_substitution_feature_matrix.tsv"
    substitution_feature_matrix.to_csv(substitution_out, sep="\t", index=False)

    observed_sites = pd.DataFrame(
        {
            "ha1_position": np.arange(1, ha1_len + 1, dtype=int),
            "n_changed_measurements": np.nansum(changes == 1.0, axis=0).astype(int),
            "n_missing_measurements": np.isnan(changes).sum(axis=0).astype(int),
        }
    )
    observed_sites_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_observed_site_changes.tsv"
    observed_sites.to_csv(observed_sites_out, sep="\t", index=False)

    missing_ratio = float(np.isnan(changes).sum() / changes.size)
    lines = [
        f"Input processed WIC rows retained for HA1 modeling: {len(df)}",
        f"Rows dropped for missing HA1-aligned virus or serum sequence: {dropped_rows}",
        f"Aligned HA1 sequences available: {len(strains)}",
        f"Aligned HA1 reference sequence id: {ref_name}",
        f"HA1 positions used: {ha1_len}",
        f"Feature matrix shape: {feature_matrix.shape}",
        f"Context feature columns: {len(meta_cols) - 12}",
        f"Distance feature columns: {distance_cols}",
        f"Primary site-state covariates: {len(change_cols) + len(virus_covariate_names) + len(serum_covariate_names)}",
        f"Observed substitution features: {len(unique_substitutions)}",
        f"Change-matrix missing fraction: {missing_ratio:.6f}",
        f"Feature matrix output: {feature_out}",
        f"Substitution matrix output: {substitution_out}",
        f"Observed-site summary: {observed_sites_out}",
    ]
    append_qc_log(cfg, "26_build_wic_model_features", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
