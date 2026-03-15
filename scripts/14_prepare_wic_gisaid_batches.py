#!/usr/bin/env python3
"""Write unresolved WIC strain names for GISAID search, optionally chunked into batches."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from wic_name_utils import preferred_query_strain


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        default="H3N2-WIC/metadata/intermediate/H3N2_WIC_repo_unmatched.tsv",
        help="Unresolved WIC strains TSV to batch for GISAID search.",
    )
    parser.add_argument(
        "--output-dir",
        default="H3N2-WIC/metadata/gisaid_batches",
        help="Directory for chunked GISAID query files when batching is enabled.",
    )
    parser.add_argument(
        "--query-output",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.txt",
        help="Combined GISAID query list written alongside the batch files.",
    )
    parser.add_argument(
        "--csv-output",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_remaining_query.csv",
        help="CSV version of the combined GISAID query list.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=0,
        help="Number of strain names per batch file. Use 0 to skip batch files and write only the combined query list.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    unresolved = pd.read_csv(args.input, sep="\t")
    unresolved["total_count"] = (
        unresolved["virus_measurement_count"].fillna(0)
        + unresolved["reference_measurement_count"].fillna(0)
    )

    queries = (
        unresolved.assign(
            query_strain=unresolved.get("gisaid_query_strain", unresolved["lookup_strain"]).fillna(
                unresolved["lookup_strain"].map(preferred_query_strain)
            )
        )[["lookup_strain", "query_strain", "total_count"]]
        .dropna()
        .drop_duplicates(subset=["lookup_strain"])
        .sort_values(["total_count", "lookup_strain"], ascending=[False, True])
        .reset_index(drop=True)
    )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    query_output = Path(args.query_output)
    query_output.parent.mkdir(parents=True, exist_ok=True)
    csv_output = Path(args.csv_output)
    csv_output.parent.mkdir(parents=True, exist_ok=True)

    for old in output_dir.glob("wic_gisaid_batch_*.txt"):
        old.unlink()
    manifest_path = output_dir / "wic_gisaid_batch_manifest.tsv"
    if manifest_path.exists():
        manifest_path.unlink()

    n = len(queries)
    queries["query_strain"].to_csv(query_output, index=False, header=False)
    csv_output.write_text(",".join(queries["query_strain"].tolist()) + "\n")

    print(f"Distinct unresolved lookup strains: {n}")
    print(f"Combined query list: {query_output}")
    print(f"CSV query list: {csv_output}")
    if args.chunk_size and args.chunk_size > 0:
        chunk_size = int(args.chunk_size)
        batch_rows = []
        for i in range(0, n, chunk_size):
            chunk = queries.iloc[i : i + chunk_size]
            batch_index = i // chunk_size + 1
            path = output_dir / f"wic_gisaid_batch_{batch_index:03d}.txt"
            chunk["query_strain"].to_csv(path, index=False, header=False)
            batch_rows.append(
                {
                    "batch": batch_index,
                    "n_strains": len(chunk),
                    "top_strain": chunk.iloc[0]["lookup_strain"],
                    "top_query_strain": chunk.iloc[0]["query_strain"],
                    "path": str(path),
                }
            )
        pd.DataFrame(batch_rows).to_csv(manifest_path, sep="\t", index=False)
        print(f"Chunk size: {chunk_size}")
        print(f"Batches written: {len(batch_rows)} -> {output_dir}")
        print(f"Batch manifest: {manifest_path}")
    else:
        print("Batch files skipped (chunk size 0).")


if __name__ == "__main__":
    main()
