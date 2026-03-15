#!/usr/bin/env python3
"""Match HI table strain names to aligned sequence strain IDs."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
import difflib

import pandas as pd
try:
    from rapidfuzz import fuzz, process
except Exception:  # noqa: BLE001
    fuzz = None
    process = None

from common import append_qc_log, load_config, normalize_strain_name, read_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    parser.add_argument("--fuzzy-threshold", type=float, default=90.0)
    return parser.parse_args()


def build_matcher(strains: list[str]):
    norm_to_strains: dict[str, list[str]] = defaultdict(list)
    for s in strains:
        norm_to_strains[normalize_strain_name(s)].append(s)

    normalized_choices = list(norm_to_strains.keys())

    def best_fuzzy(query_norm: str):
        if process is not None and fuzz is not None:
            candidate = process.extractOne(
                query_norm,
                normalized_choices,
                scorer=fuzz.ratio,
            )
            if candidate is None:
                return None, 0.0
            best_norm, score, _ = candidate
            return best_norm, float(score)

        best = difflib.get_close_matches(query_norm, normalized_choices, n=1, cutoff=0.0)
        if not best:
            return None, 0.0
        best_norm = best[0]
        score = difflib.SequenceMatcher(a=query_norm, b=best_norm).ratio() * 100.0
        return best_norm, float(score)

    def match(query: str, threshold: float):
        query_norm = normalize_strain_name(query)

        if query_norm in norm_to_strains:
            options = norm_to_strains[query_norm]
            if len(options) == 1:
                return options[0], "exact_norm", 100.0
            if query in options:
                return query, "exact_text", 100.0
            return options[0], "exact_norm_ambiguous", 100.0

        if not normalized_choices:
            return None, "no_choices", 0.0

        best_norm, score = best_fuzzy(query_norm)
        if best_norm is None:
            return None, "no_match", 0.0
        if score < threshold:
            return None, "below_threshold", float(score)

        matched = norm_to_strains[best_norm][0]
        return matched, "fuzzy", float(score)

    return match


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    hi = pd.read_csv(cfg["hi_path"], sep="\t")
    aligned_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA_aligned.fasta"
    if not aligned_path.exists():
        raise FileNotFoundError(f"Aligned FASTA not found: {aligned_path}")

    seq_records = read_fasta(str(aligned_path))
    sequence_strains = sorted(seq_records.keys())

    matcher = build_matcher(sequence_strains)

    for prefix, col in [("virus", "virusStrain"), ("serum", "serumStrain")]:
        matches = hi[col].astype(str).map(lambda s: matcher(s, args.fuzzy_threshold))
        hi[f"{prefix}_strain_matched"] = matches.map(lambda t: t[0])
        hi[f"{prefix}_match_method"] = matches.map(lambda t: t[1])
        hi[f"{prefix}_match_score"] = matches.map(lambda t: t[2])
        hi[f"{prefix}_match_needs_review"] = hi[f"{prefix}_match_method"].isin(
            ["fuzzy", "exact_norm_ambiguous"]
        )

    matched = hi[
        hi["virus_strain_matched"].notna() & hi["serum_strain_matched"].notna()
    ].copy()
    matched["match_needs_review"] = (
        matched["virus_match_needs_review"] | matched["serum_match_needs_review"]
    )

    out_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_matched_pairs.csv"
    matched.to_csv(out_path, index=False)

    review_cols = [
        "virusStrain",
        "virus_strain_matched",
        "virus_match_method",
        "virus_match_score",
        "serumStrain",
        "serum_strain_matched",
        "serum_match_method",
        "serum_match_score",
        "source",
    ]
    review_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_match_review.tsv"
    matched.loc[matched["match_needs_review"], review_cols].to_csv(
        review_path,
        sep="\t",
        index=False,
    )

    unmatched_counter = Counter()
    for col in ["virusStrain", "serumStrain"]:
        missing = hi[hi[f"{col.split('Strain')[0]}_strain_matched"].isna()][col].astype(str)
        unmatched_counter.update(missing)

    unmatched_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_unmatched_strains.tsv"
    with open(unmatched_path, "w", encoding="utf-8") as f:
        f.write("strain\tcount\n")
        for strain, count in unmatched_counter.most_common():
            f.write(f"{strain}\t{count}\n")

    unique_virus = hi["virusStrain"].nunique()
    unique_serum = hi["serumStrain"].nunique()
    matched_virus = matched["virusStrain"].nunique()
    matched_serum = matched["serumStrain"].nunique()

    lines = [
        f"HI rows input: {len(hi)}",
        f"HI rows retained (both virus+serum matched): {len(matched)}",
        f"HI rows dropped: {len(hi) - len(matched)}",
        f"Unique virus strains in HI: {unique_virus}",
        f"Unique serum strains in HI: {unique_serum}",
        f"Unique matched virus strains: {matched_virus}",
        f"Unique matched serum strains: {matched_serum}",
        f"Virus match method counts: {dict(hi['virus_match_method'].value_counts())}",
        f"Serum match method counts: {dict(hi['serum_match_method'].value_counts())}",
        f"Matched rows flagged for review: {int(matched['match_needs_review'].sum())} -> {review_path}",
        f"Matched output: {out_path}",
        f"Unmatched strain list: {unmatched_path}",
    ]
    append_qc_log(cfg, "03_match_strains", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
