#!/usr/bin/env python3
"""Fetch HA protein sequences from NCBI for IRD accessions and prep GISAID list."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO

from common import append_qc_log, load_config, write_fasta


HA_KEYWORDS = ("hemagglutinin", "ha", "ha0", "ha1", "hemaglutinin")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json", help="Path to JSON config.")
    parser.add_argument("--email", default=None, help="NCBI Entrez email.")
    parser.add_argument("--api-key", default=None, help="NCBI Entrez API key (optional).")
    parser.add_argument(
        "--sleep-seconds",
        type=float,
        default=0.34,
        help="Delay between NCBI requests (ignored if --api-key and lower than 0.11).",
    )
    parser.add_argument(
        "--skip-genbank-fetch",
        action="store_true",
        help="Only write GISAID accession list; skip NCBI download.",
    )
    return parser.parse_args()


def choose_ha_record(records):
    preferred = []
    fallback = []
    for rec in records:
        desc = rec.description.lower()
        if any(key in desc for key in HA_KEYWORDS):
            preferred.append(rec)
        else:
            fallback.append(rec)
    if preferred:
        return max(preferred, key=lambda r: len(r.seq))
    if fallback:
        return max(fallback, key=lambda r: len(r.seq))
    return None


def fetch_genbank_ha(accession: str):
    handle = Entrez.efetch(
        db="nucleotide",
        id=accession,
        rettype="fasta_cds_aa",
        retmode="text",
    )
    try:
        records = list(SeqIO.parse(handle, "fasta"))
    finally:
        handle.close()
    return choose_ha_record(records)


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    seq_df = pd.read_csv(cfg["seq_path"], sep="\t")
    seq_df["database"] = seq_df["database"].astype(str)
    seq_df["accession"] = seq_df["accession"].astype(str).str.upper().str.split(".").str[0]

    gisaid_df = seq_df[seq_df["database"].str.upper() == "GISAID"].copy()
    public_df = seq_df[seq_df["database"].str.upper().isin({"IRD", "NCBI", "GENBANK"})].copy()

    output_dir = Path(cfg["output_dir"])
    data_dir = Path(cfg["data_dir"])
    data_dir.mkdir(parents=True, exist_ok=True)

    gisaid_list_path = output_dir / f"{cfg['subtype']}_gisaid_accessions_to_download.txt"
    gisaid_df["accession"].sort_values().to_csv(gisaid_list_path, index=False, header=False)

    lines = [
        f"Total accession rows: {len(seq_df)}",
        f"GISAID accession rows: {len(gisaid_df)}",
        f"Public accession rows (IRD/NCBI/GenBank): {len(public_df)}",
        f"Wrote GISAID accession list: {gisaid_list_path}",
    ]

    if args.skip_genbank_fetch:
        append_qc_log(cfg, "01_fetch_sequences", lines + ["Skipped NCBI download by request."])
        print("\n".join(lines))
        print("Skipped GenBank fetch (--skip-genbank-fetch).")
        return

    email = args.email or cfg.get("entrez_email")
    if not email:
        raise ValueError("NCBI email is required. Pass --email or set entrez_email in config.")

    Entrez.email = email
    if args.api_key:
        Entrez.api_key = args.api_key

    min_sleep = 0.11 if args.api_key else 0.34
    sleep_seconds = max(args.sleep_seconds, min_sleep)

    accession_to_strain = dict(zip(public_df["accession"], public_df["strain"]))

    genbank_records = {}
    failures = []

    for i, accession in enumerate(sorted(accession_to_strain.keys()), start=1):
        try:
            rec = fetch_genbank_ha(accession)
            if rec is None:
                failures.append((accession, "no_ha_record_found"))
                continue
            seq = str(rec.seq).upper().replace("*", "")
            if len(seq) < 300:
                failures.append((accession, f"sequence_too_short:{len(seq)}"))
                continue
            genbank_records[accession] = seq
        except Exception as exc:  # noqa: BLE001
            failures.append((accession, f"error:{type(exc).__name__}:{exc}"))
        finally:
            if i < len(accession_to_strain):
                time.sleep(sleep_seconds)

    genbank_fasta_path = Path(cfg["genbank_fasta_path"])
    write_fasta(genbank_records, str(genbank_fasta_path))

    failed_path = output_dir / f"{cfg['subtype']}_genbank_failed_accessions.tsv"
    if failures:
        pd.DataFrame(failures, columns=["accession", "reason"]).to_csv(failed_path, sep="\t", index=False)

    lines.extend(
        [
            f"Saved GenBank HA sequences: {len(genbank_records)} -> {genbank_fasta_path}",
            f"Failed GenBank accessions: {len(failures)}",
            f"Failure file: {failed_path if failures else 'none'}",
        ]
    )

    append_qc_log(cfg, "01_fetch_sequences", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
