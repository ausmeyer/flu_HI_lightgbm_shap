#!/usr/bin/env python3
"""Auto-fill the WIC sequence manifest from existing local H3N2 sequence metadata."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

from wic_name_utils import canonical_strain, normalize_key, preferred_query_strain, strip_decorator


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        default="H3N2-WIC/manifests/H3N2_WIC_seq_manifest_template.tsv",
        help="Input WIC manifest template.",
    )
    parser.add_argument(
        "--reference-manifest",
        default="H3N2/H3N2_seq_data.tsv",
        help="Existing curated H3N2 sequence manifest used as a local reference source.",
    )
    parser.add_argument(
        "--output",
        default="H3N2-WIC/manifests/intermediate/H3N2_WIC_seq_manifest_autofill.tsv",
        help="Output auto-filled manifest.",
    )
    parser.add_argument(
        "--unresolved",
        default="H3N2-WIC/metadata/intermediate/H3N2_WIC_unresolved_strains.tsv",
        help="Output unresolved strains table.",
    )
    parser.add_argument(
        "--gisaid-batch",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_batch_query.txt",
        help="Output batch query list for unresolved canonical strain names.",
    )
    return parser.parse_args()

def build_reference_lookup(ref: pd.DataFrame) -> tuple[dict[str, dict], dict[str, dict], dict[str, dict], dict[str, dict]]:
    exact = {}
    exact_norm = {}
    base = {}
    base_norm = {}

    for row in ref.to_dict("records"):
        strain = str(row["strain"])
        exact.setdefault(strain.upper(), row)
        norm = re.sub(r"[^A-Z0-9]", "", strain.upper())
        exact_norm.setdefault(norm, row)

        base_strain = canonical_strain(strain)
        base.setdefault(base_strain, row)
        base_key = normalize_key(strain)
        if base_key not in base_norm:
            base_norm[base_key] = row
        else:
            # Ambiguous normalized base names should not auto-resolve.
            base_norm[base_key] = None

    return exact, exact_norm, base, base_norm


def main() -> None:
    args = parse_args()
    manifest = pd.read_csv(args.manifest, sep="\t")
    ref = pd.read_csv(args.reference_manifest, sep="\t")

    exact, exact_norm, base, base_norm = build_reference_lookup(ref)

    match_rows = []
    for row in manifest.to_dict("records"):
        strain = str(row["strain"])
        base_strain = canonical_strain(strain)
        norm = re.sub(r"[^A-Z0-9]", "", strain.upper())
        base_key = normalize_key(strain)

        matched = None
        match_type = ""
        if strain.upper() in exact:
            matched = exact[strain.upper()]
            match_type = "exact_strain"
        elif norm in exact_norm:
            matched = exact_norm[norm]
            match_type = "exact_normalized"
        elif base_strain in base:
            matched = base[base_strain]
            match_type = "base_strain"
        elif base_key in base_norm and base_norm[base_key] is not None:
            matched = base_norm[base_key]
            match_type = "base_normalized"

        out = dict(row)
        out["base_strain"] = base_strain
        out["lookup_strain"] = base_strain
        out["gisaid_query_strain"] = preferred_query_strain(base_strain)
        out["decorated_name"] = strain != base_strain
        out["match_type"] = match_type
        out["matched_existing_strain"] = matched["strain"] if matched else ""
        out["autofilled"] = matched is not None

        if matched:
            out["accession"] = matched["accession"]
            out["database"] = matched["database"]
            if pd.isna(out.get("year")) or str(out.get("year")) in {"", "<NA>", "nan"}:
                out["year"] = matched["year"]

        match_rows.append(out)

    auto = pd.DataFrame(match_rows)
    auto["year"] = pd.to_numeric(auto["year"], errors="coerce").astype("Int64")

    unresolved = auto[~auto["autofilled"]].copy()
    unresolved = unresolved.sort_values(
        ["reference_measurement_count", "virus_measurement_count", "strain"],
        ascending=[False, False, True],
    )

    batch_queries = (
        unresolved["gisaid_query_strain"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
    )

    output_path = Path(args.output)
    unresolved_path = Path(args.unresolved)
    batch_path = Path(args.gisaid_batch)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    auto.to_csv(output_path, sep="\t", index=False)
    unresolved.to_csv(unresolved_path, sep="\t", index=False)
    batch_queries.to_csv(batch_path, index=False, header=False)

    print(f"Input manifest rows: {len(manifest)}")
    print(f"Auto-filled rows: {int(auto['autofilled'].sum())}")
    print(f"Unresolved rows: {len(unresolved)}")
    print(f"Distinct unresolved lookup strains: {batch_queries.shape[0]}")
    print(f"Auto-filled manifest: {output_path}")
    print(f"Unresolved strains: {unresolved_path}")
    print(f"GISAID batch query list: {batch_path}")


if __name__ == "__main__":
    main()
