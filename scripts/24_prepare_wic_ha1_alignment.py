#!/usr/bin/env python3
"""Derive HA1 sequences from final WIC proteins and align them for modeling."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from common import (
    append_qc_log,
    build_alignment_position_map,
    load_config,
    read_fasta,
    sanitize_sequence,
    write_fasta,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="H3N2-WIC/config/wic_ha1_model.json")
    return parser.parse_args()


def extract_ha1(seq: str, signal_len: int, ha1_len: int, mature_len: int) -> tuple[str | None, str]:
    clean = sanitize_sequence(seq)
    if len(clean) < ha1_len:
        return None, "too_short"
    if len(clean) >= signal_len + mature_len:
        return clean[signal_len : signal_len + ha1_len], "full_ha_with_signal"
    if len(clean) >= mature_len:
        return clean[:ha1_len], "mature_full_ha"
    if len(clean) >= ha1_len:
        return clean[:ha1_len], "ha1_like"
    return None, "too_short"


def pick_reference_sequence_id(
    seq_manifest: pd.DataFrame,
    extracted_rows: pd.DataFrame,
    configured_reference: str | None,
) -> str:
    if configured_reference and configured_reference in set(extracted_rows["sequence_id"]):
        return configured_reference

    merged = extracted_rows.merge(
        seq_manifest[["isolate_id", "strain", "sequence_source", "year"]],
        left_on="sequence_id",
        right_on="isolate_id",
        how="left",
    )
    paper = merged[
        (merged["extraction_mode"] == "ha1_like")
        & (merged["sequence_source"].astype(str) == "PAPER_REPO_HA1")
    ].copy()
    if not paper.empty:
        paper["year"] = pd.to_numeric(paper["year"], errors="coerce")
        paper = paper.sort_values(["year", "strain", "sequence_id"], na_position="last")
        return str(paper.iloc[0]["sequence_id"])

    complete = merged[merged["ha1_length"] == 328].copy()
    complete = complete.sort_values(["sequence_source", "strain", "sequence_id"], na_position="last")
    if complete.empty:
        raise ValueError("No HA1-complete sequences available for reference selection.")
    return str(complete.iloc[0]["sequence_id"])


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    seq_manifest = pd.read_csv(cfg["final_seq_path"], sep="\t")
    seq_map = read_fasta(cfg["final_fasta_path"])

    signal_len = int(cfg.get("signal_peptide_length", 16))
    ha1_len = int(cfg.get("ha1_length", 328))
    mature_len = int(cfg.get("mature_protein_length", 550))

    extracted_records: dict[str, str] = {}
    extraction_rows = []
    for sequence_id, seq in seq_map.items():
        ha1_seq, mode = extract_ha1(seq, signal_len=signal_len, ha1_len=ha1_len, mature_len=mature_len)
        if ha1_seq is not None:
            extracted_records[sequence_id] = ha1_seq
        extraction_rows.append(
            {
                "sequence_id": sequence_id,
                "input_length": len(seq),
                "ha1_length": len(ha1_seq) if ha1_seq is not None else 0,
                "extraction_mode": mode,
                "retained": ha1_seq is not None,
            }
        )

    extraction_df = pd.DataFrame(extraction_rows)
    extraction_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_ha1_extraction_summary.tsv"
    extraction_df.to_csv(extraction_out, sep="\t", index=False)

    if not extracted_records:
        raise ValueError("No HA1 sequences could be extracted from the final FASTA.")

    reference_id = pick_reference_sequence_id(
        seq_manifest=seq_manifest,
        extracted_rows=extraction_df[extraction_df["retained"]].copy(),
        configured_reference=cfg.get("reference_sequence_id"),
    )

    unaligned_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA1_unaligned.fasta"
    write_fasta(extracted_records, str(unaligned_path))

    if shutil.which("mafft") is None:
        raise RuntimeError("MAFFT not found on PATH. Install MAFFT before running WIC HA1 alignment.")

    aligned_raw_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA1_aligned_raw.fasta"
    with open(aligned_raw_path, "w", encoding="utf-8") as out_f:
        proc = subprocess.run(
            ["mafft", "--auto", str(unaligned_path)],
            check=True,
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
        )

    aligned = read_fasta(str(aligned_raw_path))
    if reference_id not in aligned:
        raise ValueError(f"Reference sequence id not found in aligned HA1 FASTA: {reference_id}")

    aligned_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA1_aligned.fasta"
    write_fasta(aligned, str(aligned_out))

    position_map = build_alignment_position_map(
        ref_seq=aligned[reference_id],
        ha1_length=ha1_len,
    )
    position_map = position_map[position_map["mature_position"].notna()].copy()
    position_map = position_map[position_map["mature_position"] <= ha1_len].copy()

    position_map_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_ha1_position_map.csv"
    position_map.to_csv(position_map_out, index=False)

    metadata = {
        "subtype": cfg["subtype"],
        "input_fasta": cfg["final_fasta_path"],
        "unaligned_ha1_fasta": str(unaligned_path),
        "aligned_ha1_fasta": str(aligned_out),
        "reference_sequence_id": reference_id,
        "signal_peptide_length": signal_len,
        "ha1_length": ha1_len,
        "n_input_sequences": len(seq_map),
        "n_retained_ha1_sequences": len(extracted_records),
        "extraction_summary_path": str(extraction_out),
        "position_map_path": str(position_map_out),
    }
    metadata_out = Path(cfg["output_dir"]) / f"{cfg['subtype']}_ha1_alignment_metadata.json"
    with open(metadata_out, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    lines = [
        f"Input final FASTA sequences: {len(seq_map)}",
        f"Retained HA1 sequences: {len(extracted_records)}",
        f"Selected HA1 reference sequence id: {reference_id}",
        f"Unaligned HA1 FASTA: {unaligned_path}",
        f"Aligned HA1 FASTA: {aligned_out}",
        f"HA1 extraction summary: {extraction_out}",
        f"HA1 position map: {position_map_out}",
        f"Alignment metadata: {metadata_out}",
    ]
    append_qc_log(cfg, "24_prepare_wic_ha1_alignment", lines)
    print("\n".join(lines))
    if proc.stderr:
        print("MAFFT log:\n" + proc.stderr.strip())


if __name__ == "__main__":
    main()
