#!/usr/bin/env python3
"""Resolve multiple GISAID isolate IDs for the same WIC strain into one canonical choice."""

from __future__ import annotations

import argparse
import hashlib
import re
from itertools import combinations
from pathlib import Path

import pandas as pd

from common import AMBIGUOUS_AA, read_fasta, write_fasta

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
HA1_COMPLETE_MIN_LEN = 328


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapping",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_isolate_id_mapping.tsv",
        help="WIC strain to GISAID isolate-ID mapping table.",
    )
    parser.add_argument(
        "--gisaid-fasta",
        default="H3N2-WIC/data/raw/gisaid/gisaid_epiflu_sequence.fasta",
        help="Local GISAID FASTA with headers containing EPI_ISL IDs.",
    )
    parser.add_argument(
        "--candidate-audit",
        default="H3N2-WIC/output/H3N2_WIC_gisaid_candidate_audit.tsv",
        help="Full candidate audit table with one row per query strain / isolate ID.",
    )
    parser.add_argument(
        "--canonical-table",
        default="H3N2-WIC/metadata/H3N2_WIC_gisaid_canonical_isolates.tsv",
        help="Canonical isolate choice per WIC query strain.",
    )
    parser.add_argument(
        "--canonical-fasta",
        default="H3N2-WIC/data/intermediate/H3N2_WIC_gisaid_canonical_ha.fasta",
        help="Canonical FASTA keyed by WIC query strain.",
    )
    parser.add_argument(
        "--summary",
        default="H3N2-WIC/output/H3N2_WIC_gisaid_duplicate_resolution_summary.tsv",
        help="Summary of duplicate-resolution outcomes.",
    )
    return parser.parse_args()


def parse_gisaid_header(header: str) -> dict[str, str]:
    parts = header.split("|")
    while len(parts) < 6:
        parts.append("")
    return {
        "isolate_id": parts[0],
        "header_strain": parts[1].replace("_", "/"),
        "header_collection_date": parts[2],
        "header_passage": parts[3],
        "header_submission_date": parts[4],
        "header_lab": parts[5],
        "header_raw": header,
    }


def ambiguity_count(seq: str) -> int:
    return sum(1 for aa in seq if aa not in STANDARD_AA)


def passage_info(passage: str) -> tuple[str, int, int]:
    raw = str(passage or "").strip()
    upper = raw.upper().replace("_", " ")
    digits = [int(x) for x in re.findall(r"\d+", upper)]
    burden = sum(digits) if digits else 0

    is_original = any(
        token in upper
        for token in ["ORIGINAL", "ORIGINAL SPECIMEN", "DIRECT", "SPECIMEN", "CLINICAL"]
    ) or upper in {"CS"}
    has_egg = "EGG" in upper or bool(re.search(r"(^|[+/ ,])E\d", upper))
    has_cell = any(token in upper for token in ["SIAT", "MDCK", "CELL", "CX"]) or bool(
        re.search(r"(^|[+/ ,])S\d", upper)
    ) or bool(re.search(r"(^|[+/ ,])C\d", upper))

    if is_original:
        return "original", 0, burden
    if has_cell and not has_egg:
        return "cell", 1, burden
    if has_egg and not has_cell:
        return "egg", 2, burden
    if has_egg and has_cell:
        return "mixed", 3, burden
    if raw:
        return "unknown", 4, burden
    return "unknown", 4, burden


def seq_distance(seq_a: str, seq_b: str) -> int:
    if len(seq_a) != len(seq_b):
        return abs(len(seq_a) - len(seq_b)) + sum(a != b for a, b in zip(seq_a, seq_b))
    return sum(a != b for a, b in zip(seq_a, seq_b))


def sequence_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:16]


def choose_reason(
    chosen: pd.Series,
    group: pd.DataFrame,
    group_n_ids: int,
    group_n_unique_sequences: int,
    selected_class: str,
    group_n_selected_class_variants: int,
    top_count: int,
    second_count: int,
) -> str:
    if group_n_ids == 1:
        return "single_isolate_id"
    if group_n_unique_sequences == 1:
        return "identical_sequences"
    if selected_class == "ha1_complete":
        incomplete_exists = bool((~group["is_ha1_complete"]).any())
        if incomplete_exists:
            return "prefer_complete_ha1"
    if group_n_selected_class_variants == 1:
        return "prefer_complete_ha1" if selected_class == "ha1_complete" else "prefer_longest_available"

    selected_group = (
        group[group["is_ha1_complete"]].copy()
        if selected_class == "ha1_complete"
        else group[group["is_max_length"]].copy()
    )
    min_passage_score = int(selected_group["passage_score"].min())
    if int(chosen["passage_score"]) == min_passage_score and min_passage_score < int(
        selected_group["passage_score"].max()
    ):
        return "prefer_passage_class"
    if int(chosen["variant_count_at_max_length"]) == top_count and top_count > second_count:
        return "prefer_modal_ha1_complete_variant" if selected_class == "ha1_complete" else "prefer_modal_longest_variant"
    if int(chosen["ambiguity_count"]) == int(selected_group["ambiguity_count"].min()):
        return "prefer_low_ambiguity"
    return "tie_break"


def main() -> None:
    args = parse_args()
    mapping_path = Path(args.mapping)
    gisaid_fasta_path = Path(args.gisaid_fasta)
    candidate_audit_path = Path(args.candidate_audit)
    canonical_table_path = Path(args.canonical_table)
    canonical_fasta_path = Path(args.canonical_fasta)
    summary_path = Path(args.summary)

    mapping = pd.read_csv(mapping_path, sep="\t")
    fasta_records = read_fasta(str(gisaid_fasta_path))

    header_meta = {}
    for header, seq in fasta_records.items():
        meta = parse_gisaid_header(header)
        meta["sequence"] = seq
        meta["length"] = len(seq)
        meta["ambiguity_count"] = ambiguity_count(seq)
        meta["sequence_hash"] = sequence_hash(seq)
        meta["sequence_raw"] = seq
        header_meta[meta["isolate_id"]] = meta

    rows = []
    for row in mapping.to_dict("records"):
        isolate_id = str(row["Isolate_Id"])
        meta = header_meta.get(isolate_id)
        if meta is None:
            continue
        passage_class, passage_score, passage_burden = passage_info(meta["header_passage"])
        rows.append(
            {
                **row,
                "header_strain": meta["header_strain"],
                "header_collection_date": meta["header_collection_date"],
                "header_passage": meta["header_passage"],
                "header_submission_date": meta["header_submission_date"],
                "header_lab": meta["header_lab"],
                "header_raw": meta["header_raw"],
                "sequence_length": meta["length"],
                "ambiguity_count": meta["ambiguity_count"],
                "sequence_hash": meta["sequence_hash"],
                "sequence_raw": meta["sequence_raw"],
                "passage_class": passage_class,
                "passage_score": passage_score,
                "passage_burden": passage_burden,
            }
        )

    if not rows:
        raise RuntimeError("No GISAID isolate IDs from the mapping table were found in the local FASTA.")

    audit = pd.DataFrame(rows)
    chosen_rows = []
    chosen_fasta = {}
    audit_parts = []

    for query_strain, group in audit.groupby("query_strain", sort=True):
        group = group.copy()
        group_n_ids = int(group["Isolate_Id"].nunique())
        group_n_unique_sequences = int(group["sequence_hash"].nunique())
        max_len = int(group["sequence_length"].max())
        group["is_max_length"] = group["sequence_length"] == max_len
        group["is_ha1_complete"] = group["sequence_length"] >= HA1_COMPLETE_MIN_LEN

        if bool(group["is_ha1_complete"].any()):
            selected_class = "ha1_complete"
            selected_group = group[group["is_ha1_complete"]].copy()
        else:
            selected_class = "longest_available"
            selected_group = group[group["is_max_length"]].copy()

        max_len_group = group[group["is_max_length"]].copy()
        variant_counts = selected_group["sequence_hash"].value_counts()
        group["variant_count_at_max_length"] = group["sequence_hash"].map(variant_counts).fillna(0).astype(int)
        group_n_max_length_variants = int(max_len_group["sequence_hash"].nunique())
        group_n_selected_class_variants = int(selected_group["sequence_hash"].nunique())

        unique_selected = (
            selected_group.drop_duplicates("sequence_hash")[["sequence_hash", "sequence_raw"]].to_dict("records")
        )
        pairwise_diffs = [
            seq_distance(a["sequence_raw"], b["sequence_raw"]) for a, b in combinations(unique_selected, 2)
        ]
        min_pairwise_diff_selected = min(pairwise_diffs) if pairwise_diffs else 0
        max_pairwise_diff_selected = max(pairwise_diffs) if pairwise_diffs else 0

        # Prefer the closest representation of the named clinical strain:
        # complete HA1 if available, then passage class, then low ambiguity, then supported exact variant.
        sort_cols = [
            "is_ha1_complete",
            "passage_score",
            "ambiguity_count",
            "variant_count_at_max_length",
            "passage_burden",
            "Collection_Date",
            "Isolate_Id",
        ]
        group["is_ha1_complete_sort"] = (~group["is_ha1_complete"]).astype(int)
        group = group.sort_values(
            by=[
                "is_ha1_complete_sort",
                "passage_score",
                "ambiguity_count",
                "variant_count_at_max_length",
                "passage_burden",
                "Collection_Date",
                "Isolate_Id",
            ],
            ascending=[True, True, True, False, True, True, True],
        ).reset_index(drop=True)

        chosen = group.iloc[0].copy()
        top_counts = variant_counts.tolist()
        top_count = top_counts[0] if top_counts else 0
        second_count = top_counts[1] if len(top_counts) > 1 else 0
        top_frac = top_count / int(max_len_group.shape[0]) if len(max_len_group) else 0.0
        selection_reason = choose_reason(
            chosen,
            group,
            group_n_ids,
            group_n_unique_sequences,
            selected_class,
            group_n_selected_class_variants,
            top_count,
            second_count,
        )
        review_flag = bool(
            (group_n_selected_class_variants > 1 and top_frac <= 0.5)
            or (group_n_selected_class_variants > 1 and max_pairwise_diff_selected >= 4)
            or (
                chosen["passage_class"] in {"egg", "mixed"}
                and int(selected_group["passage_score"].min()) < int(chosen["passage_score"])
            )
        )

        group["group_n_ids"] = group_n_ids
        group["group_n_unique_sequences"] = group_n_unique_sequences
        group["group_max_length"] = max_len
        group["group_n_max_length_variants"] = group_n_max_length_variants
        group["group_selected_class"] = selected_class
        group["group_n_selected_class_variants"] = group_n_selected_class_variants
        group["group_top_variant_count"] = top_count
        group["group_second_variant_count"] = second_count
        group["group_top_variant_fraction"] = top_frac
        group["group_min_pairwise_diff_selected_class"] = min_pairwise_diff_selected
        group["group_max_pairwise_diff_selected_class"] = max_pairwise_diff_selected
        group["selection_reason"] = selection_reason
        group["review_flag"] = review_flag
        group["chosen_isolate"] = group["Isolate_Id"] == chosen["Isolate_Id"]
        audit_parts.append(group)

        chosen_out = chosen.to_dict()
        chosen_out.update(
            {
                "group_n_ids": group_n_ids,
                "group_n_unique_sequences": group_n_unique_sequences,
                "group_max_length": max_len,
                "group_n_max_length_variants": group_n_max_length_variants,
                "group_selected_class": selected_class,
                "group_n_selected_class_variants": group_n_selected_class_variants,
                "group_top_variant_count": top_count,
                "group_second_variant_count": second_count,
                "group_top_variant_fraction": top_frac,
                "group_min_pairwise_diff_selected_class": min_pairwise_diff_selected,
                "group_max_pairwise_diff_selected_class": max_pairwise_diff_selected,
                "selection_reason": selection_reason,
                "review_flag": review_flag,
            }
        )
        chosen_rows.append(chosen_out)
        chosen_fasta[query_strain] = chosen["sequence_raw"]

    audit_out = pd.concat(audit_parts, ignore_index=True)
    canonical = pd.DataFrame(chosen_rows)

    candidate_audit_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    canonical_table_path.parent.mkdir(parents=True, exist_ok=True)
    canonical_fasta_path.parent.mkdir(parents=True, exist_ok=True)

    audit_out = audit_out.drop(columns=["sequence_raw", "is_ha1_complete_sort"], errors="ignore")
    canonical_out = canonical.drop(columns=["sequence_raw", "is_ha1_complete_sort"], errors="ignore")

    audit_out.to_csv(candidate_audit_path, sep="\t", index=False)
    canonical_out.to_csv(canonical_table_path, sep="\t", index=False)
    write_fasta(chosen_fasta, str(canonical_fasta_path))

    summary = pd.DataFrame(
        [
            {"metric": "matched_query_strains", "value": int(canonical_out["query_strain"].nunique())},
            {"metric": "total_candidate_rows", "value": int(audit_out.shape[0])},
            {"metric": "query_strains_with_multiple_ids", "value": int((canonical_out["group_n_ids"] > 1).sum())},
            {
                "metric": "duplicate_groups_identical_sequences",
                "value": int(
                    ((canonical_out["group_n_ids"] > 1) & (canonical_out["group_n_unique_sequences"] == 1)).sum()
                ),
            },
            {
                "metric": "duplicate_groups_partial_only_variation",
                "value": int(
                    (
                        (canonical_out["group_n_ids"] > 1)
                        & (canonical_out["group_n_unique_sequences"] > 1)
                        & (canonical_out["group_n_max_length_variants"] == 1)
                    ).sum()
                ),
            },
            {
                "metric": "duplicate_groups_multiple_full_length_variants",
                "value": int(
                    (
                        (canonical_out["group_n_ids"] > 1)
                        & (canonical_out["group_n_selected_class_variants"] > 1)
                    ).sum()
                ),
            },
            {
                "metric": "canonical_selected_class_ha1_complete",
                "value": int((canonical_out["group_selected_class"] == "ha1_complete").sum()),
            },
            {
                "metric": "canonical_selected_class_longest_available",
                "value": int((canonical_out["group_selected_class"] == "longest_available").sum()),
            },
            {"metric": "canonical_review_flagged", "value": int(canonical_out["review_flag"].sum())},
            {"metric": "canonical_passage_original", "value": int((canonical_out["passage_class"] == "original").sum())},
            {"metric": "canonical_passage_cell", "value": int((canonical_out["passage_class"] == "cell").sum())},
            {"metric": "canonical_passage_unknown", "value": int((canonical_out["passage_class"] == "unknown").sum())},
            {"metric": "canonical_passage_egg", "value": int((canonical_out["passage_class"] == "egg").sum())},
            {"metric": "canonical_passage_mixed", "value": int((canonical_out["passage_class"] == "mixed").sum())},
        ]
    )
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Candidate audit: {candidate_audit_path}")
    print(f"Canonical isolate table: {canonical_table_path}")
    print(f"Canonical FASTA: {canonical_fasta_path}")
    print(f"Summary: {summary_path}")
    print(f"Matched query strains: {canonical_out['query_strain'].nunique()}")
    print(f"Query strains with multiple isolate IDs: {(canonical_out['group_n_ids'] > 1).sum()}")
    print(f"Review-flagged canonical choices: {int(canonical_out['review_flag'].sum())}")


if __name__ == "__main__":
    main()
