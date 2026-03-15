#!/usr/bin/env python3
"""Audit which HI strains are missing from the manifest or lost during alignment/QC."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import append_qc_log, load_config, read_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/h3n2.json")
    return parser.parse_args()


def summarize_role_counts(hi: pd.DataFrame, strains: set[str], column: str) -> pd.Series:
    return hi.loc[hi[column].isin(strains), column].value_counts().rename(f"{column}_count")


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    hi = pd.read_csv(cfg["hi_path"], sep="\t")
    seq = pd.read_csv(cfg["seq_path"], sep="\t")

    seq_strains = set(seq["strain"].astype(str))
    aligned_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_HA_aligned.fasta"
    aligned_strains = set(read_fasta(str(aligned_path)).keys())

    hi_virus = set(hi["virusStrain"].astype(str))
    hi_serum = set(hi["serumStrain"].astype(str))
    hi_union = hi_virus | hi_serum

    missing_from_manifest = sorted(hi_union - seq_strains)
    manifest_only = hi[hi["virusStrain"].isin(missing_from_manifest) | hi["serumStrain"].isin(missing_from_manifest)].copy()
    manifest_summary = pd.DataFrame({"strain": missing_from_manifest})
    if not manifest_summary.empty:
        manifest_summary = manifest_summary.merge(
            summarize_role_counts(hi, set(missing_from_manifest), "virusStrain"),
            left_on="strain",
            right_index=True,
            how="left",
        ).merge(
            summarize_role_counts(hi, set(missing_from_manifest), "serumStrain"),
            left_on="strain",
            right_index=True,
            how="left",
        )
        manifest_summary["virusStrain_count"] = manifest_summary["virusStrain_count"].fillna(0).astype(int)
        manifest_summary["serumStrain_count"] = manifest_summary["serumStrain_count"].fillna(0).astype(int)
        manifest_summary["total_hi_measurements"] = (
            manifest_summary["virusStrain_count"] + manifest_summary["serumStrain_count"]
        )
        manifest_summary = manifest_summary.sort_values("total_hi_measurements", ascending=False)

    failed_qc = sorted((hi_union & seq_strains) - aligned_strains)
    failed_qc_summary = pd.DataFrame({"strain": failed_qc})
    if not failed_qc_summary.empty:
        failed_qc_summary = failed_qc_summary.merge(
            summarize_role_counts(hi, set(failed_qc), "virusStrain"),
            left_on="strain",
            right_index=True,
            how="left",
        ).merge(
            summarize_role_counts(hi, set(failed_qc), "serumStrain"),
            left_on="strain",
            right_index=True,
            how="left",
        ).merge(
            seq[["strain", "accession", "database", "year"]],
            on="strain",
            how="left",
        )
        failed_qc_summary["virusStrain_count"] = failed_qc_summary["virusStrain_count"].fillna(0).astype(int)
        failed_qc_summary["serumStrain_count"] = failed_qc_summary["serumStrain_count"].fillna(0).astype(int)
        failed_qc_summary["total_hi_measurements"] = (
            failed_qc_summary["virusStrain_count"] + failed_qc_summary["serumStrain_count"]
        )
        failed_qc_summary = failed_qc_summary.sort_values("total_hi_measurements", ascending=False)

    manifest_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_missing_from_seq_manifest.tsv"
    manifest_summary.to_csv(manifest_path, sep="\t", index=False)

    qc_path = Path(cfg["output_dir"]) / f"{cfg['subtype']}_present_but_failed_qc.tsv"
    failed_qc_summary.to_csv(qc_path, sep="\t", index=False)

    lines = [
        f"HI unique strains (virus or serum): {len(hi_union)}",
        f"Unique strains missing from sequence manifest: {len(missing_from_manifest)} -> {manifest_path}",
        f"Unique strains present in manifest but absent after QC: {len(failed_qc)} -> {qc_path}",
    ]
    append_qc_log(cfg, "09_audit_strain_coverage", lines)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
