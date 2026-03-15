#!/usr/bin/env python3
"""Build final WIC modeling inputs restricted to HI rows with sequences on both sides."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import extract_accession_candidates, read_fasta, write_fasta
from wic_name_utils import canonical_strain, normalize_key


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hi-path",
        default="H3N2-WIC/data/raw/WIC-HI-H3N2.csv",
        help="Input WIC HI table.",
    )
    parser.add_argument(
        "--manifest",
        default="H3N2-WIC/manifests/H3N2_WIC_seq_manifest_gisaidfill.tsv",
        help="WIC sequence manifest after repo + GISAID matching.",
    )
    parser.add_argument(
        "--repo-fasta",
        default="H3N2-WIC/data/intermediate/H3N2_WIC_repo_matched_ha1.fasta",
        help="Repo-matched HA1 FASTA.",
    )
    parser.add_argument(
        "--gisaid-fasta",
        default="H3N2-WIC/data/intermediate/H3N2_WIC_gisaid_canonical_ha.fasta",
        help="Canonical GISAID FASTA keyed by WIC strain.",
    )
    parser.add_argument(
        "--public-fasta",
        default="H3N2-WIC/genbank_sequences.fasta",
        help="Optional local public-repository FASTA.",
    )
    parser.add_argument(
        "--output-hi",
        default="H3N2-WIC/data/final/H3N2_WIC_HI_data_final.tsv",
        help="Final filtered HI table.",
    )
    parser.add_argument(
        "--output-manifest",
        default="H3N2-WIC/data/final/H3N2_WIC_seq_data_final.tsv",
        help="Final filtered sequence manifest.",
    )
    parser.add_argument(
        "--output-fasta",
        default="H3N2-WIC/data/final/H3N2_WIC_sequences_final.fasta",
        help="Final FASTA with one identifier per sequence.",
    )
    parser.add_argument(
        "--summary",
        default="H3N2-WIC/output/H3N2_WIC_final_input_summary.tsv",
        help="Summary TSV for final modeling inputs.",
    )
    return parser.parse_args()


def load_optional_fasta(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    return read_fasta(str(path))


def build_accession_lookup(records: dict[str, str]) -> dict[str, dict[str, str]]:
    lookup = {}
    for header, seq in records.items():
        for accession in extract_accession_candidates(header):
            if accession not in lookup:
                lookup[accession] = {"header": header, "sequence": seq}
    return lookup


def bool_value(value) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def build_available_sequences(
    manifest: pd.DataFrame,
    repo_fasta: dict[str, str],
    gisaid_fasta: dict[str, str],
    public_lookup: dict[str, dict[str, str]],
) -> pd.DataFrame:
    rows = []
    for row in manifest.to_dict("records"):
        strain = str(row["strain"])
        sequence = None
        sequence_source = None
        source_header = ""

        if bool_value(row.get("gisaid_matched")) and strain in gisaid_fasta:
            sequence = gisaid_fasta[strain]
            sequence_source = "GISAID"
            source_header = str(row.get("gisaid_header", ""))
        elif bool_value(row.get("repo_matched")) and strain in repo_fasta:
            sequence = repo_fasta[strain]
            sequence_source = "PAPER_REPO_HA1"
            source_header = str(row.get("repo_header", ""))
        elif str(row.get("database", "")).strip() in {"NCBI", "IRD"}:
            accession = str(row.get("accession", "")).strip()
            public_match = public_lookup.get(accession)
            if public_match is not None:
                sequence = public_match["sequence"]
                sequence_source = str(row.get("database", "")).strip()
                source_header = public_match["header"]

        if sequence is None:
            continue

        isolate_id = str(row.get("accession", "")).strip()
        if not isolate_id or isolate_id in {"nan", "<NA>"}:
            continue

        out = dict(row)
        out["isolate_id"] = isolate_id
        out["sequence_source"] = sequence_source
        out["sequence_length"] = len(sequence)
        out["source_fasta_header"] = source_header
        out["sequence"] = sequence
        out["normalized_strain_key"] = normalize_key(strain)
        out["is_decorated_variant"] = bool_value(row.get("decorated_name"))
        rows.append(out)

    return pd.DataFrame(rows)


def pick_best_manifest_row(group: pd.DataFrame) -> pd.Series:
    group = group.copy()
    group["source_priority"] = group["sequence_source"].map(
        {"GISAID": 0, "PAPER_REPO_HA1": 1, "IRD": 2, "NCBI": 3}
    ).fillna(9)
    group["sequence_length_sort"] = -group["sequence_length"]
    group["measurement_sort"] = -group.get("total_measurement_count", pd.Series(0, index=group.index)).fillna(0)
    group["decorated_sort"] = group["is_decorated_variant"].astype(int)
    group = group.sort_values(
        by=["decorated_sort", "source_priority", "measurement_sort", "sequence_length_sort", "strain"],
        ascending=[True, True, True, True, True],
    )
    return group.iloc[0]


def infer_year_from_row(row: pd.Series) -> pd.Int64Dtype:
    value = row.get("year", pd.NA)
    if pd.isna(value):
        return pd.NA
    try:
        return int(float(value))
    except Exception:
        return pd.NA


def main() -> None:
    args = parse_args()

    hi_path = Path(args.hi_path)
    manifest_path = Path(args.manifest)
    repo_fasta_path = Path(args.repo_fasta)
    gisaid_fasta_path = Path(args.gisaid_fasta)
    public_fasta_path = Path(args.public_fasta)

    output_hi_path = Path(args.output_hi)
    output_manifest_path = Path(args.output_manifest)
    output_fasta_path = Path(args.output_fasta)
    summary_path = Path(args.summary)

    hi = pd.read_csv(hi_path)
    manifest = pd.read_csv(manifest_path, sep="\t", low_memory=False)
    repo_fasta = read_fasta(str(repo_fasta_path))
    gisaid_fasta = read_fasta(str(gisaid_fasta_path))
    public_fasta = load_optional_fasta(public_fasta_path)
    public_lookup = build_accession_lookup(public_fasta)

    available = build_available_sequences(
        manifest,
        repo_fasta=repo_fasta,
        gisaid_fasta=gisaid_fasta,
        public_lookup=public_lookup,
    )
    if available.empty:
        raise RuntimeError("No locally available WIC sequences were found in the manifest.")

    chosen_by_key = {
        key: pick_best_manifest_row(group)
        for key, group in available.groupby("normalized_strain_key", sort=False)
    }

    retained_rows = []
    for row in hi.to_dict("records"):
        virus_key = normalize_key(canonical_strain(row["virus"]))
        serum_key = normalize_key(canonical_strain(row["reference"]))
        virus_meta = chosen_by_key.get(virus_key)
        serum_meta = chosen_by_key.get(serum_key)
        if virus_meta is None or serum_meta is None:
            continue

        retained_rows.append(
            {
                "virusIsolate": row["virus"],
                "virusStrain": row["virus"],
                "virusYear": infer_year_from_row(virus_meta),
                "serumIsolate": row["reference"],
                "serumStrain": row["reference"],
                "serumYear": infer_year_from_row(serum_meta),
                "titer": row["titre"],
                "source": "WIC",
                "virus_passage": row["virus_passage"],
                "reference_passage": row["reference_passage"],
                "antiserum": row["antiserum"],
                "RBC": row["RBC"],
                "date": row["date"],
                "pair": row["pair"],
                "pair_strict": row["pair_strict"],
                "virusSequenceId": virus_meta["isolate_id"],
                "serumSequenceId": serum_meta["isolate_id"],
                "virusSequenceSource": virus_meta["sequence_source"],
                "serumSequenceSource": serum_meta["sequence_source"],
                "virusSequenceLength": int(virus_meta["sequence_length"]),
                "serumSequenceLength": int(serum_meta["sequence_length"]),
            }
        )

    final_hi = pd.DataFrame(retained_rows)
    if final_hi.empty:
        raise RuntimeError("No HI rows retained after requiring sequences for both virus and serum strains.")

    used_keys = set(final_hi["virusStrain"].map(lambda x: normalize_key(canonical_strain(x))))
    used_keys.update(final_hi["serumStrain"].map(lambda x: normalize_key(canonical_strain(x))))

    used_manifest_rows = []
    final_fasta_records = {}
    seen_ids = {}
    for key in sorted(used_keys):
        meta = chosen_by_key[key]
        used_manifest_rows.append(meta.to_dict())
        isolate_id = str(meta["isolate_id"])
        seq = str(meta["sequence"])
        if isolate_id in final_fasta_records and final_fasta_records[isolate_id] != seq:
            raise RuntimeError(f"Sequence ID {isolate_id} maps to conflicting sequences.")
        final_fasta_records[isolate_id] = seq
        seen_ids[isolate_id] = meta["strain"]

    final_manifest = pd.DataFrame(used_manifest_rows).drop(columns=["sequence", "normalized_strain_key"], errors="ignore")
    final_manifest = final_manifest.drop(
        columns=["source_priority", "sequence_length_sort", "measurement_sort", "decorated_sort"],
        errors="ignore",
    )
    final_manifest = final_manifest.sort_values(["sequence_source", "strain"]).reset_index(drop=True)

    output_hi_path.parent.mkdir(parents=True, exist_ok=True)
    output_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    output_fasta_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    final_hi.to_csv(output_hi_path, sep="\t", index=False)
    final_manifest.to_csv(output_manifest_path, sep="\t", index=False)
    write_fasta(final_fasta_records, str(output_fasta_path))

    summary = pd.DataFrame(
        [
            {"metric": "input_hi_rows", "value": int(len(hi))},
            {"metric": "retained_hi_rows", "value": int(len(final_hi))},
            {"metric": "retained_fraction", "value": float(len(final_hi) / len(hi))},
            {"metric": "retained_unique_virus_strains", "value": int(final_hi["virusStrain"].nunique())},
            {"metric": "retained_unique_serum_strains", "value": int(final_hi["serumStrain"].nunique())},
            {"metric": "retained_unique_sequence_manifest_rows", "value": int(len(final_manifest))},
            {"metric": "retained_unique_fasta_ids", "value": int(len(final_fasta_records))},
            {
                "metric": "manifest_rows_gisaid",
                "value": int((final_manifest["sequence_source"] == "GISAID").sum()),
            },
            {
                "metric": "manifest_rows_paper_repo",
                "value": int((final_manifest["sequence_source"] == "PAPER_REPO_HA1").sum()),
            },
            {
                "metric": "manifest_rows_ird",
                "value": int((final_manifest["sequence_source"] == "IRD").sum()),
            },
            {
                "metric": "manifest_rows_ncbi",
                "value": int((final_manifest["sequence_source"] == "NCBI").sum()),
            },
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Final HI table: {output_hi_path}")
    print(f"Final sequence manifest: {output_manifest_path}")
    print(f"Final FASTA: {output_fasta_path}")
    print(f"Summary: {summary_path}")
    print(f"Retained HI rows: {len(final_hi)} / {len(hi)}")
    print(f"Retained unique sequence IDs: {len(final_fasta_records)}")


if __name__ == "__main__":
    main()
