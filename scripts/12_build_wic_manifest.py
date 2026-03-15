#!/usr/bin/env python3
"""Build a deduplicated sequence-manifest template for the WIC H3N2 HI dataset."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="H3N2-WIC/data/raw/WIC-HI-H3N2.csv",
        help="Path to the WIC HI CSV.",
    )
    parser.add_argument(
        "--output",
        default="H3N2-WIC/manifests/H3N2_WIC_seq_manifest_template.tsv",
        help="Output TSV manifest template.",
    )
    parser.add_argument(
        "--summary",
        default="H3N2-WIC/metadata/H3N2_WIC_manifest_summary.tsv",
        help="Output TSV with manifest summary counts.",
    )
    return parser.parse_args()


def infer_year(strain: str) -> int | None:
    match = re.search(r"/(\d{4})$", str(strain).strip())
    if not match:
        return None
    year = int(match.group(1))
    if 1900 <= year <= 2100:
        return year
    return None


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    summary_path = Path(args.summary)

    df = pd.read_csv(input_path, usecols=["virus", "reference"])

    virus_counts = df["virus"].value_counts().rename("virus_measurement_count")
    reference_counts = df["reference"].value_counts().rename("reference_measurement_count")

    all_strains = pd.Index(sorted(set(df["virus"]).union(df["reference"])), name="strain")
    manifest = pd.DataFrame({"strain": all_strains})
    manifest = manifest.merge(virus_counts, left_on="strain", right_index=True, how="left")
    manifest = manifest.merge(reference_counts, left_on="strain", right_index=True, how="left")
    manifest["virus_measurement_count"] = manifest["virus_measurement_count"].fillna(0).astype(int)
    manifest["reference_measurement_count"] = manifest["reference_measurement_count"].fillna(0).astype(int)
    manifest["in_virus_column"] = manifest["virus_measurement_count"] > 0
    manifest["in_reference_column"] = manifest["reference_measurement_count"] > 0
    manifest["wic_role"] = "virus_only"
    manifest.loc[manifest["in_reference_column"], "wic_role"] = "reference_only"
    manifest.loc[
        manifest["in_virus_column"] & manifest["in_reference_column"],
        "wic_role",
    ] = "virus_and_reference"

    manifest["accession"] = ""
    manifest["type"] = "A"
    manifest["subtype"] = "H3N2"
    manifest["year"] = manifest["strain"].map(infer_year).astype("Int64")
    manifest["database"] = ""
    manifest = manifest[
        [
            "strain",
            "accession",
            "type",
            "subtype",
            "year",
            "database",
            "wic_role",
            "virus_measurement_count",
            "reference_measurement_count",
        ]
    ]
    manifest = manifest.sort_values(["wic_role", "strain"]).reset_index(drop=True)

    summary = pd.DataFrame(
        [
            {"metric": "input_rows", "value": len(df)},
            {"metric": "unique_virus_strains", "value": df["virus"].nunique()},
            {"metric": "unique_reference_strains", "value": df["reference"].nunique()},
            {"metric": "unique_union_strains", "value": len(all_strains)},
            {"metric": "virus_only_strains", "value": int((manifest["wic_role"] == "virus_only").sum())},
            {
                "metric": "reference_only_strains",
                "value": int((manifest["wic_role"] == "reference_only").sum()),
            },
            {
                "metric": "virus_and_reference_strains",
                "value": int((manifest["wic_role"] == "virus_and_reference").sum()),
            },
            {"metric": "year_inferred_rows", "value": int(manifest["year"].notna().sum())},
            {"metric": "year_missing_rows", "value": int(manifest["year"].isna().sum())},
        ]
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(output_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    print(f"Unique virus strains: {df['virus'].nunique()}")
    print(f"Unique reference strains: {df['reference'].nunique()}")
    print(f"Unique union strains: {len(all_strains)}")
    print(f"Manifest template: {output_path}")
    print(f"Manifest summary: {summary_path}")


if __name__ == "__main__":
    main()
