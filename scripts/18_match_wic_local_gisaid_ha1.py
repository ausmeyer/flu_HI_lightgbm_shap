#!/usr/bin/env python3
"""Match a locally downloaded GISAID HA1 FASTA to the WIC manifest after repo matching."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import extract_accession_candidates, read_fasta, write_fasta
from wic_name_utils import canonical_strain, extract_header_strain_candidates, normalize_key


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_repofill.tsv",
        help="Input manifest after repo matching.",
    )
    parser.add_argument(
        "--gisaid-fasta",
        default="H3N2-WIC/data/raw/gisaid/gisaid_epiflu_sequence.fasta",
        help="Locally downloaded GISAID FASTA.",
    )
    parser.add_argument(
        "--output-manifest",
        default="H3N2-WIC/manifests/H3N2_WIC_seq_manifest_gisaidfill.tsv",
        help="Manifest after local GISAID FASTA matching.",
    )
    parser.add_argument(
        "--summary",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_match_summary.tsv",
        help="Summary TSV for local GISAID matching.",
    )
    parser.add_argument(
        "--canonical-table",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_canonical_isolates.tsv",
        help="Optional canonical isolate table produced by duplicate resolution.",
    )
    parser.add_argument(
        "--matched-fasta",
        default="H3N2-WIC/data/intermediate/H3N2_WIC_gisaid_matched_ha1.fasta",
        help="Matched HA1 FASTA rewritten with WIC strain headers.",
    )
    return parser.parse_args()


def build_header_lookup(records: dict[str, str]) -> dict[str, dict]:
    lookup = {}
    for header, seq in records.items():
        candidates = extract_header_strain_candidates(header)
        accession_candidates = extract_accession_candidates(header)
        for candidate in candidates:
            key = normalize_key(candidate)
            if key and key not in lookup:
                lookup[key] = {
                    "header": header,
                    "matched_strain": canonical_strain(candidate),
                    "sequence": seq,
                    "accession": accession_candidates[0] if accession_candidates else "",
                }
    return lookup


def build_isolate_lookup(records: dict[str, str]) -> dict[str, dict]:
    lookup = {}
    for header, seq in records.items():
        accession_candidates = extract_accession_candidates(header)
        isolate_id = accession_candidates[0] if accession_candidates else header.split("|", 1)[0]
        candidates = extract_header_strain_candidates(header)
        lookup[str(isolate_id)] = {
            "header": header,
            "matched_strain": canonical_strain(candidates[0]) if candidates else "",
            "sequence": seq,
            "accession": isolate_id,
        }
    return lookup


def resolve_gisaid_fasta_path(path: Path) -> Path:
    if path.exists():
        return path

    epiflu_fallback = path.with_name("gisaid_epiflu_sequence.fasta")
    if epiflu_fallback.exists():
        return epiflu_fallback

    legacy_fallback = path.with_name("gisaid_sequences.fasta")
    if legacy_fallback.exists():
        return legacy_fallback

    return path


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest)
    gisaid_fasta_path = resolve_gisaid_fasta_path(Path(args.gisaid_fasta))
    output_manifest_path = Path(args.output_manifest)
    summary_path = Path(args.summary)
    canonical_table_path = Path(args.canonical_table)
    matched_fasta_path = Path(args.matched_fasta)

    if not gisaid_fasta_path.exists():
        raise FileNotFoundError(f"GISAID FASTA not found: {gisaid_fasta_path}")

    manifest = pd.read_csv(manifest_path, sep="\t")
    records = read_fasta(str(gisaid_fasta_path))
    header_lookup = build_header_lookup(records)
    isolate_lookup = build_isolate_lookup(records)
    canonical_lookup = {}
    if canonical_table_path.exists():
        canonical_table = pd.read_csv(canonical_table_path, sep="\t")
        for row in canonical_table.to_dict("records"):
            canonical_lookup[normalize_key(row["query_strain"])] = row

    matched_rows = []
    matched_fasta = {}
    for row in manifest.to_dict("records"):
        lookup_strain = canonical_strain(row.get("lookup_strain", row.get("base_strain", row["strain"])))
        key = normalize_key(lookup_strain)
        matched = None
        match_type = ""

        canonical_row = canonical_lookup.get(key)
        if canonical_row is not None:
            isolate_id = str(canonical_row["Isolate_Id"])
            matched = isolate_lookup.get(isolate_id)
            if matched is not None:
                match_type = "canonical_isolate_id"

        if matched is None:
            matched = header_lookup.get(key)
            if matched is not None:
                match_type = "header_strain_match"

        row["gisaid_match_type"] = match_type if matched else ""
        row["gisaid_header"] = matched["header"] if matched else ""
        row["gisaid_matched_strain"] = matched["matched_strain"] if matched else ""
        row["gisaid_matched"] = matched is not None
        row["gisaid_aa_length"] = len(matched["sequence"]) if matched else pd.NA
        if canonical_row is not None:
            row["gisaid_selected_isolate_id"] = canonical_row["Isolate_Id"]
            row["gisaid_selected_passage"] = canonical_row.get("header_passage", "")
            row["gisaid_selected_lab"] = canonical_row.get("header_lab", "")
            row["gisaid_selection_reason"] = canonical_row.get("selection_reason", "")
            row["gisaid_selection_review_flag"] = canonical_row.get("review_flag", False)
        else:
            row["gisaid_selected_isolate_id"] = ""
            row["gisaid_selected_passage"] = ""
            row["gisaid_selected_lab"] = ""
            row["gisaid_selection_reason"] = ""
            row["gisaid_selection_review_flag"] = False

        if matched and str(row.get("accession", "")).strip() in {"", "nan", "<NA>"}:
            row["database"] = "GISAID"
            row["accession"] = matched["accession"] or matched["header"]

        if matched:
            matched_fasta[str(row["strain"])] = matched["sequence"]

        matched_rows.append(row)

    matched_df = pd.DataFrame(matched_rows)
    output_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    matched_df.to_csv(output_manifest_path, sep="\t", index=False)
    write_fasta(matched_fasta, str(matched_fasta_path))

    summary = pd.DataFrame(
        [
            {"metric": "gisaid_fasta_sequences", "value": len(records)},
            {"metric": "manifest_rows", "value": len(matched_df)},
            {"metric": "gisaid_matched_rows", "value": int(matched_df["gisaid_matched"].sum())},
            {
                "metric": "gisaid_matched_distinct_strains",
                "value": int(matched_df.loc[matched_df["gisaid_matched"], "strain"].nunique()),
            },
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"GISAID FASTA: {gisaid_fasta_path}")
    print(f"GISAID FASTA sequences: {len(records)}")
    print(f"Manifest rows: {len(matched_df)}")
    print(f"GISAID-matched rows: {int(matched_df['gisaid_matched'].sum())}")
    print(f"GISAID-filled manifest: {output_manifest_path}")
    print(f"Matched HA1 FASTA: {matched_fasta_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
