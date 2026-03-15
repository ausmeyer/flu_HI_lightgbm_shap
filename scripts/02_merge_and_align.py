#!/usr/bin/env python3
"""Merge sequence sources, align HA proteins, trim signal peptide, and run QC."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from common import (
    append_qc_log,
    build_alignment_position_map,
    load_config,
    pick_accession_from_header,
    sanitize_sequence,
    write_fasta,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--skip-alignment", action="store_true", help="Skip MAFFT and only merge FASTA inputs.")
    return parser.parse_args()


def read_fasta_records(path: Path):
    if not path.exists():
        return []
    return list(SeqIO.parse(str(path), "fasta"))


def pick_best_existing(existing: str, new_seq: str) -> str:
    if len(new_seq) > len(existing):
        return new_seq
    return existing


def trim_by_reference(aligned: dict[str, str], reference: str, signal_len: int) -> tuple[dict[str, str], int]:
    ref_seq = aligned.get(reference)
    if ref_seq is None:
        raise ValueError(f"Reference strain {reference!r} not found in alignment.")

    residue_count = 0
    cut_col = None
    for idx, aa in enumerate(ref_seq):
        if aa != "-":
            residue_count += 1
            if residue_count == signal_len:
                cut_col = idx + 1
                break

    if cut_col is None:
        raise ValueError(
            f"Reference strain has fewer than {signal_len} non-gap residues before trimming point."
        )

    trimmed = {strain: seq[cut_col:] for strain, seq in aligned.items()}
    return trimmed, cut_col


def extract_reference_coding_region(
    trimmed: dict[str, str],
    reference: str,
    mature_len: int,
    gap_fill: str = "X",
) -> tuple[dict[str, str], list[int], dict[str, int]]:
    ref_seq = trimmed.get(reference)
    if ref_seq is None:
        raise ValueError(f"Reference strain {reference!r} not found in trimmed alignment.")

    ref_cols = [i for i, aa in enumerate(ref_seq) if aa != "-"]
    if len(ref_cols) < mature_len:
        raise ValueError(
            f"Reference has {len(ref_cols)} coding residues after trimming; expected at least {mature_len}."
        )

    coding_cols = ref_cols[:mature_len]
    extracted = {}
    gap_replacements = {}
    for strain, seq in trimmed.items():
        chars = [seq[idx] for idx in coding_cols]
        replaced = sum(1 for aa in chars if aa == "-")
        extracted[strain] = "".join(gap_fill if aa == "-" else aa for aa in chars)
        gap_replacements[strain] = replaced

    return extracted, coding_cols, gap_replacements


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    seq_df = pd.read_csv(cfg["seq_path"], sep="\t")
    seq_df["accession"] = seq_df["accession"].astype(str).str.upper().str.split(".").str[0]

    accession_to_strain = dict(zip(seq_df["accession"], seq_df["strain"]))
    known_accessions = set(accession_to_strain)

    genbank_path = Path(cfg["genbank_fasta_path"])
    gisaid_path = Path(cfg["gisaid_fasta_path"])

    source_records = {
        "genbank": read_fasta_records(genbank_path),
        "gisaid": read_fasta_records(gisaid_path),
    }

    strain_to_seq: dict[str, str] = {}
    accession_hits = 0
    unmatched_headers = []
    duplicate_strains = 0

    for source_name, records in source_records.items():
        for rec in records:
            header = f"{rec.id} {rec.description}"
            accession = pick_accession_from_header(header, known_accessions)
            if accession is None:
                unmatched_headers.append(f"{source_name}:{rec.description}")
                continue

            strain = accession_to_strain[accession]
            seq = sanitize_sequence(str(rec.seq))
            if not seq:
                continue

            accession_hits += 1
            if strain in strain_to_seq:
                duplicate_strains += 1
                strain_to_seq[strain] = pick_best_existing(strain_to_seq[strain], seq)
            else:
                strain_to_seq[strain] = seq

    merged_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA_proteins.fasta"
    write_fasta(strain_to_seq, str(merged_path))

    missing_strains = sorted(set(seq_df["strain"]) - set(strain_to_seq))
    missing_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_missing_sequence_strains.txt"
    with open(missing_path, "w", encoding="utf-8") as f:
        for s in missing_strains:
            f.write(s + "\n")

    if args.skip_alignment:
        lines = [
            f"Merged protein FASTA: {merged_path}",
            f"Sequence table strains: {seq_df['strain'].nunique()}",
            f"Strains with merged sequence: {len(strain_to_seq)}",
            f"Accession-matched FASTA records: {accession_hits}",
            f"Duplicate strain records collapsed: {duplicate_strains}",
            f"Missing strains: {len(missing_strains)} -> {missing_path}",
            f"Unmatched FASTA headers: {len(unmatched_headers)}",
            "Skipped MAFFT/QC by request.",
        ]
        append_qc_log(cfg, "02_merge_and_align", lines)
        print("\n".join(lines))
        return

    if shutil.which("mafft") is None:
        raise RuntimeError("MAFFT not found on PATH. Install MAFFT before running alignment.")

    aligned_raw_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA_aligned_raw.fasta"
    with open(aligned_raw_path, "w", encoding="utf-8") as out_f:
        proc = subprocess.run(
            ["mafft", "--auto", str(merged_path)],
            check=True,
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
        )

    aligned = {rec.id: sanitize_sequence(str(rec.seq)) for rec in SeqIO.parse(str(aligned_raw_path), "fasta")}
    trimmed, cut_col = trim_by_reference(
        aligned=aligned,
        reference=cfg["reference_strain"],
        signal_len=int(cfg["signal_peptide_length"]),
    )

    coding_region, coding_cols, gap_replacements = extract_reference_coding_region(
        trimmed=trimmed,
        reference=cfg["reference_strain"],
        mature_len=int(cfg.get("mature_protein_length", 550)),
    )

    aligned_out_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA_aligned.fasta"
    write_fasta(coding_region, str(aligned_out_path))

    if cfg["reference_strain"] not in coding_region:
        raise RuntimeError(
            "Reference strain missing from coding-region alignment."
        )

    trimmed_alignment_len = len(next(iter(coding_region.values()))) if coding_region else 0

    position_map = build_alignment_position_map(
        coding_region[cfg["reference_strain"]],
        ha1_length=int(cfg.get("ha1_length", 328)),
    )
    position_map_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_alignment_position_map.csv"
    position_map.to_csv(position_map_path, index=False)

    mature_map_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_mature_position_map.csv"
    position_map[position_map["mature_position"].notna()].to_csv(mature_map_path, index=False)

    metadata_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_alignment_metadata.json"
    metadata = {
        "subtype": cfg["subtype"],
        "merged_fasta": str(merged_path),
        "aligned_raw_fasta": str(aligned_raw_path),
        "aligned_qc_fasta": str(aligned_out_path),
        "reference_strain": cfg["reference_strain"],
        "signal_peptide_length": int(cfg["signal_peptide_length"]),
        "signal_trim_column_in_raw_alignment": cut_col,
        "reference_coding_columns_in_trimmed_alignment": [int(x) + 1 for x in coding_cols],
        "trimmed_alignment_length": trimmed_alignment_len,
        "strains_before_qc": len(trimmed),
        "strains_after_qc": len(coding_region),
        "ha1_length": int(cfg.get("ha1_length", 328)),
        "alignment_position_map": str(position_map_path),
        "mature_position_map": str(mature_map_path),
    }
    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    lines = [
        f"Merged protein FASTA: {merged_path}",
        f"Sequence table strains: {seq_df['strain'].nunique()}",
        f"Strains with merged sequence: {len(strain_to_seq)}",
        f"Accession-matched FASTA records: {accession_hits}",
        f"Duplicate strain records collapsed: {duplicate_strains}",
        f"Missing strains from sequence table: {len(missing_strains)} -> {missing_path}",
        f"Unmatched FASTA headers: {len(unmatched_headers)}",
        f"MAFFT output: {aligned_raw_path}",
        f"Signal peptide trimmed using reference {cfg['reference_strain']} at raw column {cut_col}",
        f"Reference coding columns retained: {len(coding_cols)}",
        f"Gap characters replaced with X across retained strains: {sum(gap_replacements.values())}",
        f"Final aligned sequences retained: {len(coding_region)} -> {aligned_out_path}",
        f"Coding-region alignment length: {trimmed_alignment_len}",
        f"Alignment position map: {position_map_path}",
        f"Mature position map: {mature_map_path}",
        f"Alignment metadata: {metadata_path}",
    ]
    append_qc_log(cfg, "02_merge_and_align", lines)
    print("\n".join(lines))
    if proc.stderr:
        print("MAFFT log:\n" + proc.stderr.strip())


if __name__ == "__main__":
    main()
