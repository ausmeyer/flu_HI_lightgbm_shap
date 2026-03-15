#!/usr/bin/env python3
"""Match WIC H3N2 strains to the published HA1 amino-acid FASTA from Flu_g2p_mapping."""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import pandas as pd

from common import read_fasta, write_fasta
from wic_name_utils import canonical_strain, normalize_key

REPO_FASTA_URL = (
    "https://raw.githubusercontent.com/will-harvey/Flu_g2p_mapping/main/"
    "dataset_h3n2/h3n2_ha1_aa.fasta"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_autofill.tsv",
        help="Input WIC manifest after local autofill.",
    )
    parser.add_argument(
        "--repo-fasta",
        default="H3N2-WIC/data/raw/repo/h3n2_ha1_aa_repo.fasta",
        help="Local path for the published repo HA1 FASTA.",
    )
    parser.add_argument(
        "--repo-url",
        default=REPO_FASTA_URL,
        help="Raw URL for the published repo HA1 FASTA.",
    )
    parser.add_argument(
        "--output-manifest",
        default="H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_repofill.tsv",
        help="Manifest with repo-FASTA matches filled in.",
    )
    parser.add_argument(
        "--repo-unmatched",
        default="H3N2-WIC/metadata/intermediate/H3N2_WIC_repo_unmatched.tsv",
        help="Rows still unmatched after repo-FASTA matching.",
    )
    parser.add_argument(
        "--repo-summary",
        default="H3N2-WIC/metadata/H3N2_WIC_repo_match_summary.tsv",
        help="Match summary TSV.",
    )
    parser.add_argument(
        "--matched-fasta",
        default="H3N2-WIC/data/intermediate/H3N2_WIC_repo_matched_ha1.fasta",
        help="Matched HA1 FASTA rewritten with WIC strain headers.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download the repo FASTA from --repo-url before matching.",
    )
    return parser.parse_args()

def normalize_repo_header(header: str) -> str:
    s = str(header).strip()
    if s.startswith(">"):
        s = s[1:]
    s = s.split()[0]
    s = s.replace("_", "/").upper()
    s = "/".join([p for p in s.split("/") if p])
    return s


def download_fasta(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response:  # noqa: S310
        data = response.read().decode("utf-8")
    path.write_text(data, encoding="utf-8")


def build_repo_lookup(records: dict[str, str]) -> tuple[dict[str, dict], dict[str, dict], dict[str, dict]]:
    exact = {}
    canonical = {}
    normalized = {}
    for header, seq in records.items():
        repo_strain = normalize_repo_header(header)
        row = {
            "repo_header": header,
            "repo_strain": repo_strain,
            "sequence": seq,
        }
        exact.setdefault(repo_strain, row)

        canon = canonical_strain(repo_strain)
        canonical.setdefault(canon, row)

        key = normalize_key(repo_strain)
        if key not in normalized:
            normalized[key] = row
        else:
            normalized[key] = None
    return exact, canonical, normalized


def main() -> None:
    args = parse_args()
    manifest_path = Path(args.manifest)
    repo_fasta_path = Path(args.repo_fasta)
    output_manifest_path = Path(args.output_manifest)
    unmatched_path = Path(args.repo_unmatched)
    summary_path = Path(args.repo_summary)
    matched_fasta_path = Path(args.matched_fasta)

    if args.download or not repo_fasta_path.exists():
        download_fasta(args.repo_url, repo_fasta_path)

    manifest = pd.read_csv(manifest_path, sep="\t")
    records = read_fasta(str(repo_fasta_path))

    exact_lookup, canonical_lookup, normalized_lookup = build_repo_lookup(records)

    matched_rows = []
    matched_fasta = {}
    for row in manifest.to_dict("records"):
        strain = str(row["strain"])
        base_strain = row.get("base_strain", canonical_strain(strain))
        lookup_strain = canonical_strain(row.get("lookup_strain", base_strain))
        matched = None
        match_type = ""

        if lookup_strain in exact_lookup:
            matched = exact_lookup[lookup_strain]
            match_type = "repo_exact"
        elif canonical_strain(base_strain) in canonical_lookup:
            matched = canonical_lookup[canonical_strain(base_strain)]
            match_type = "repo_canonical"
        else:
            key = normalize_key(lookup_strain)
            if key in normalized_lookup and normalized_lookup[key] is not None:
                matched = normalized_lookup[key]
                match_type = "repo_normalized"

        row["repo_match_type"] = match_type
        row["repo_header"] = matched["repo_header"] if matched else ""
        row["repo_strain"] = matched["repo_strain"] if matched else ""
        row["repo_matched"] = matched is not None
        row["repo_aa_length"] = len(matched["sequence"]) if matched else pd.NA

        if matched:
            row["database"] = "PAPER_REPO_HA1"
            row["accession"] = matched["repo_header"]
            matched_fasta[strain] = matched["sequence"]

        matched_rows.append(row)

    matched_df = pd.DataFrame(matched_rows)
    output_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    matched_df.to_csv(output_manifest_path, sep="\t", index=False)

    unmatched = matched_df[~matched_df["repo_matched"]].copy()
    unmatched.to_csv(unmatched_path, sep="\t", index=False)

    write_fasta(matched_fasta, str(matched_fasta_path))

    summary = pd.DataFrame(
        [
            {"metric": "repo_fasta_sequences", "value": len(records)},
            {"metric": "manifest_rows", "value": len(matched_df)},
            {"metric": "repo_matched_rows", "value": int(matched_df["repo_matched"].sum())},
            {
                "metric": "repo_matched_distinct_strains",
                "value": int(matched_df.loc[matched_df["repo_matched"], "strain"].nunique()),
            },
            {"metric": "repo_unmatched_rows", "value": len(unmatched)},
            {
                "metric": "repo_unmatched_distinct_lookup_strains",
                "value": int(unmatched.get("lookup_strain", unmatched["strain"]).astype(str).nunique()),
            },
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Repo FASTA: {repo_fasta_path}")
    print(f"Repo sequences: {len(records)}")
    print(f"Manifest rows: {len(matched_df)}")
    print(f"Repo-matched rows: {int(matched_df['repo_matched'].sum())}")
    print(f"Repo-unmatched rows: {len(unmatched)}")
    print(f"Repo-filled manifest: {output_manifest_path}")
    print(f"Repo-unmatched rows: {unmatched_path}")
    print(f"Matched HA1 FASTA: {matched_fasta_path}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
