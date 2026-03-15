#!/usr/bin/env python3
"""Shared helpers for the influenza HI + sequence modeling pipeline."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Mapping

import numpy as np
import pandas as pd
from Bio import Phylo

AMBIGUOUS_AA = set("XBZJ")


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = json.load(f)

    output_dir = Path(cfg["output_dir"])
    figures_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    cfg["output_dir"] = str(output_dir)
    cfg["figures_dir"] = str(figures_dir)
    cfg.setdefault("qc_log_path", str(output_dir / f"{cfg['subtype']}_qc_log.txt"))
    cfg.setdefault("n_splits", 5)
    cfg.setdefault("random_state", 42)
    cfg.setdefault("ha1_length", 328)
    cfg.setdefault("homologous_fallback", "max_observed")
    cfg.setdefault("comparison_target_column", "standardized_titer")
    cfg.setdefault("distance_feature_cols", ["temporal_distance"])
    cfg.setdefault("prediction_color_feature", cfg["distance_feature_cols"][0])
    return cfg


def append_qc_log(cfg: Mapping[str, object], stage: str, lines: Iterable[str]) -> None:
    log_path = Path(str(cfg["qc_log_path"]))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(f"[{timestamp}] {stage}\n")
        for line in lines:
            f.write(f"  - {line}\n")
        f.write("\n")


def sanitize_sequence(seq: str) -> str:
    return (
        seq.upper()
        .replace("*", "")
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
    )


def read_fasta(path: str) -> Dict[str, str]:
    records: Dict[str, str] = {}
    current_id = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = sanitize_sequence("".join(chunks))
                current_id = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line)

    if current_id is not None:
        records[current_id] = sanitize_sequence("".join(chunks))
    return records


def write_fasta(records: Mapping[str, str], path: str, line_width: int = 80) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        for record_id, seq in records.items():
            clean = sanitize_sequence(seq)
            f.write(f">{record_id}\n")
            for i in range(0, len(clean), line_width):
                f.write(clean[i : i + line_width] + "\n")


def normalize_strain_name(name: str) -> str:
    if name is None:
        return ""
    s = str(name).strip().upper()
    s = s.replace(" ", "")
    return re.sub(r"[^A-Z0-9]", "", s)


def extract_accession_candidates(text: str) -> List[str]:
    upper = text.upper()
    epi = re.findall(r"EPI_ISL_\d+", upper)
    generic = re.findall(r"[A-Z]{1,4}\d{4,9}(?:\.\d+)?", upper)
    merged = []
    seen = set()
    for c in epi + generic:
        normalized = c.split(".")[0]
        if normalized not in seen:
            merged.append(normalized)
            seen.add(normalized)
    return merged


def pick_accession_from_header(header: str, known_accessions: set[str]) -> str | None:
    for candidate in extract_accession_candidates(header):
        if candidate in known_accessions:
            return candidate

    upper = header.upper()
    for acc in known_accessions:
        if acc in upper:
            return acc
    return None


def gap_fraction(seq: str) -> float:
    if not seq:
        return 1.0
    return seq.count("-") / len(seq)


def ambiguous_fraction(seq: str) -> float:
    nongap = [c for c in seq if c != "-"]
    if not nongap:
        return 1.0
    ambig = sum(1 for c in nongap if c in AMBIGUOUS_AA)
    return ambig / len(nongap)


def ungapped_length(seq: str) -> int:
    return len(seq.replace("-", ""))


def build_lab_from_source(source: str) -> str:
    s = str(source).strip().lower()
    if s.startswith("nimr"):
        return "nimr_crick"
    if s.startswith("cdc"):
        return "cdc"
    if s.startswith("vidrl"):
        return "vidrl"
    if s.startswith("niid"):
        return "niid"
    return "other"


def parse_censored_titer(value: str) -> float:
    s = str(value).strip()
    m = re.match(r"^([<>]?)(\d+(?:\.\d+)?)$", s)
    if not m:
        raise ValueError(f"Unexpected titer value: {value!r}")
    sign, num = m.groups()
    x = float(num)
    if sign == "<":
        return x / 2.0
    if sign == ">":
        return x * 2.0
    return x


def get_lightgbm_params(random_state: int = 42) -> dict:
    return {
        "objective": "regression",
        "metric": "rmse",
        "boosting_type": "gbdt",
        "num_leaves": 31,
        "max_depth": -1,
        "learning_rate": 0.1,
        "n_estimators": 1000,
        "min_child_samples": 20,
        "subsample_for_bin": 200000,
        "min_split_gain": 0.0,
        "reg_alpha": 0.0,
        "reg_lambda": 0.0,
        "colsample_bytree": 1.0,
        "subsample": 1.0,
        "subsample_freq": 0,
        "random_state": random_state,
        "verbose": -1,
    }


def build_alignment_position_map(
    ref_seq: str,
    ha1_length: int = 328,
) -> pd.DataFrame:
    rows = []
    mature_position = 0
    for trimmed_col, aa in enumerate(ref_seq, start=1):
        if aa == "-":
            rows.append(
                {
                    "trimmed_alignment_column": trimmed_col,
                    "reference_residue": aa,
                    "mature_position": pd.NA,
                    "ha_subunit": "gap",
                    "subunit_position": pd.NA,
                }
            )
            continue

        mature_position += 1
        subunit = "HA1" if mature_position <= ha1_length else "HA2"
        subunit_position = mature_position if subunit == "HA1" else mature_position - ha1_length
        rows.append(
            {
                "trimmed_alignment_column": trimmed_col,
                "reference_residue": aa,
                "mature_position": mature_position,
                "ha_subunit": subunit,
                "subunit_position": subunit_position,
            }
        )

    out = pd.DataFrame(rows)
    out["mature_position"] = out["mature_position"].astype("Int64")
    out["subunit_position"] = out["subunit_position"].astype("Int64")
    return out


def get_distance_feature_cols(cfg: Mapping[str, object]) -> list[str]:
    cols = cfg.get("distance_feature_cols", ["temporal_distance"])
    if isinstance(cols, str):
        cols = [cols]
    return [str(col) for col in cols]


def get_prediction_color_feature(cfg: Mapping[str, object]) -> str:
    return str(cfg.get("prediction_color_feature", get_distance_feature_cols(cfg)[0]))


def humanize_distance_label(feature: str) -> str:
    mapping = {
        "temporal_distance": "temporal distance",
        "patristic_distance": "patristic distance",
    }
    return mapping.get(feature, feature.replace("_", " "))


def resolve_fasttree_binary() -> str:
    for candidate in ["FastTree", "fasttree", "FastTreeMP"]:
        path = shutil.which(candidate)
        if path is not None:
            return path
    raise RuntimeError("FastTree not found on PATH. Install FastTree before running patristic analyses.")


def build_patristic_tree(
    aligned_fasta_path: str,
    tree_out_path: str,
    tip_depths_out_path: str,
    metadata_out_path: str,
    model_name: str = "WAG",
) -> dict[str, object]:
    fasttree_bin = resolve_fasttree_binary()
    tree_out = Path(tree_out_path)
    tree_out.parent.mkdir(parents=True, exist_ok=True)

    cmd = [fasttree_bin]
    if model_name.lower() == "wag":
        cmd.append("-wag")
    cmd.extend(["-gamma", aligned_fasta_path])
    proc = subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    tree_out.write_text(proc.stdout, encoding="utf-8")

    tree = Phylo.read(str(tree_out), "newick")
    depths = tree.depths()
    tip_rows = []
    for clade in tree.get_terminals():
        tip_rows.append(
            {
                "tip": str(clade.name),
                "root_to_tip_distance": float(depths.get(clade, 0.0)),
            }
        )
    pd.DataFrame(tip_rows).sort_values("tip").to_csv(tip_depths_out_path, sep="\t", index=False)

    branch_lengths = [float(clade.branch_length) for clade in tree.find_clades() if clade.branch_length is not None]
    metadata = {
        "fasttree_binary": fasttree_bin,
        "fasttree_model": model_name,
        "n_tips": int(len(tip_rows)),
        "n_branches_with_lengths": int(len(branch_lengths)),
        "sum_branch_length": float(np.sum(branch_lengths)) if branch_lengths else 0.0,
        "mean_branch_length": float(np.mean(branch_lengths)) if branch_lengths else 0.0,
        "tree_path": str(tree_out),
        "tip_depths_path": str(tip_depths_out_path),
        "fasttree_stderr": proc.stderr.strip(),
    }
    Path(metadata_out_path).write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def compute_patristic_distances(
    tree_path: str,
    left_ids: Iterable[object],
    right_ids: Iterable[object],
) -> pd.Series:
    tree = Phylo.read(tree_path, "newick")
    terminals = {str(clade.name) for clade in tree.get_terminals()}

    pairs = [(str(left), str(right)) for left, right in zip(left_ids, right_ids)]
    missing = sorted({node for pair in pairs for node in pair if node not in terminals})
    if missing:
        preview = ", ".join(missing[:10])
        raise ValueError(f"Tree is missing {len(missing)} ids needed for patristic distance lookup: {preview}")

    cache: dict[tuple[str, str], float] = {}
    values = []
    for left, right in pairs:
        key = (left, right) if left <= right else (right, left)
        if key not in cache:
            cache[key] = float(tree.distance(key[0], key[1]))
        values.append(cache[key])
    return pd.Series(values, dtype=float)


def annotate_homologous_titers(
    df: pd.DataFrame,
    serum_id_col: str,
    serum_strain_col: str,
    virus_strain_col: str,
    log_titer_col: str,
    fallback: str = "max_observed",
) -> pd.DataFrame:
    out = df.copy()
    homologous_mask = out[virus_strain_col].astype(str) == out[serum_strain_col].astype(str)

    homologous_by_serum = (
        out.loc[homologous_mask]
        .groupby(serum_id_col)[log_titer_col]
        .max()
        .rename("homologous_log2_titer")
    )
    overall_max_by_serum = (
        out.groupby(serum_id_col)[log_titer_col]
        .max()
        .rename("fallback_log2_titer")
    )

    out["homologous_log2_titer"] = out[serum_id_col].map(homologous_by_serum)
    out["homologous_titer_source"] = np.where(
        out["homologous_log2_titer"].notna(),
        "homologous",
        "missing",
    )

    if fallback == "max_observed":
        fallback_vals = out[serum_id_col].map(overall_max_by_serum)
        missing = out["homologous_log2_titer"].isna()
        out.loc[missing, "homologous_log2_titer"] = fallback_vals.loc[missing]
        out.loc[missing, "homologous_titer_source"] = "max_observed"

    if out["homologous_log2_titer"].isna().any():
        missing_count = int(out["homologous_log2_titer"].isna().sum())
        raise ValueError(
            f"Unable to determine homologous proxy titers for {missing_count} rows. "
            "Check serum IDs and fallback settings."
        )

    out["standardized_titer"] = out["homologous_log2_titer"] - out[log_titer_col]
    return out


def fit_additive_effects(
    df: pd.DataFrame,
    target_col: str,
    virus_col: str,
    serum_col: str,
) -> dict[str, object]:
    work = df[[target_col, virus_col, serum_col]].copy()
    work[target_col] = pd.to_numeric(work[target_col], errors="coerce")
    if work[target_col].isna().any():
        raise ValueError(f"NaN values found in additive-effects target column {target_col!r}")

    y = work[target_col].to_numpy(dtype=float)
    serum = work[serum_col].astype(str)
    virus = work[virus_col].astype(str)

    serum_levels = sorted(serum.unique())
    virus_levels = sorted(virus.unique())

    design = [np.ones(len(work), dtype=float)]
    serum_effect_names = serum_levels[1:]
    virus_effect_names = virus_levels[1:]

    for level in serum_effect_names:
        design.append((serum == level).to_numpy(dtype=float))
    for level in virus_effect_names:
        design.append((virus == level).to_numpy(dtype=float))

    X = np.column_stack(design)
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    intercept = float(coef[0])
    serum_coefs = coef[1 : 1 + len(serum_effect_names)]
    virus_coefs = coef[1 + len(serum_effect_names) :]

    serum_effects = {serum_levels[0]: 0.0}
    serum_effects.update(
        {level: float(value) for level, value in zip(serum_effect_names, serum_coefs)}
    )

    virus_effects = {virus_levels[0]: 0.0}
    virus_effects.update(
        {level: float(value) for level, value in zip(virus_effect_names, virus_coefs)}
    )

    fitted = (
        intercept
        + serum.map(serum_effects).to_numpy(dtype=float)
        + virus.map(virus_effects).to_numpy(dtype=float)
    )
    residual = y - fitted

    return {
        "intercept": intercept,
        "serum_effects": serum_effects,
        "virus_effects": virus_effects,
        "fitted": fitted,
        "residual": residual,
        "n_sera": len(serum_levels),
        "n_viruses": len(virus_levels),
    }


def apply_additive_effects(
    df: pd.DataFrame,
    serum_col: str,
    virus_col: str,
    intercept: float,
    serum_effects: Mapping[str, float],
    virus_effects: Mapping[str, float],
) -> np.ndarray:
    serum_vals = df[serum_col].astype(str).map(lambda x: serum_effects.get(x, 0.0)).to_numpy(dtype=float)
    virus_vals = df[virus_col].astype(str).map(lambda x: virus_effects.get(x, 0.0)).to_numpy(dtype=float)
    return intercept + serum_vals + virus_vals


def fit_additive_effects_multi(
    df: pd.DataFrame,
    target_col: str,
    factor_cols: list[str],
) -> dict[str, object]:
    if not factor_cols:
        raise ValueError("factor_cols must contain at least one column.")

    work = df[[target_col] + factor_cols].copy()
    work[target_col] = pd.to_numeric(work[target_col], errors="coerce")
    if work[target_col].isna().any():
        raise ValueError(f"NaN values found in additive-effects target column {target_col!r}")

    y = work[target_col].to_numpy(dtype=float)
    design = [np.ones(len(work), dtype=float)]
    factor_levels: dict[str, list[str]] = {}
    factor_effect_names: dict[str, list[str]] = {}

    for col in factor_cols:
        values = work[col].astype(str)
        levels = sorted(values.unique())
        factor_levels[col] = levels
        effect_names = levels[1:]
        factor_effect_names[col] = effect_names
        for level in effect_names:
            design.append((values == level).to_numpy(dtype=float))

    X = np.column_stack(design)
    coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

    intercept = float(coef[0])
    offset = 1
    factor_effects: dict[str, dict[str, float]] = {}
    fitted = np.full(len(work), intercept, dtype=float)

    for col in factor_cols:
        levels = factor_levels[col]
        effect_names = factor_effect_names[col]
        n_effects = len(effect_names)
        effect_values = coef[offset : offset + n_effects]
        offset += n_effects

        effects = {levels[0]: 0.0}
        effects.update(
            {level: float(value) for level, value in zip(effect_names, effect_values)}
        )
        factor_effects[col] = effects
        fitted += work[col].astype(str).map(effects).to_numpy(dtype=float)

    residual = y - fitted
    return {
        "intercept": intercept,
        "factor_cols": list(factor_cols),
        "factor_effects": factor_effects,
        "fitted": fitted,
        "residual": residual,
        "n_rows": len(work),
        "n_factors": len(factor_cols),
        "n_levels_total": int(sum(len(v) for v in factor_levels.values())),
    }


def apply_additive_effects_multi(
    df: pd.DataFrame,
    factor_cols: list[str],
    intercept: float,
    factor_effects: Mapping[str, Mapping[str, float]],
) -> np.ndarray:
    fitted = np.full(len(df), intercept, dtype=float)
    for col in factor_cols:
        effects = factor_effects.get(col, {})
        fitted += (
            df[col]
            .astype(str)
            .map(lambda x: effects.get(x, 0.0))
            .to_numpy(dtype=float)
        )
    return fitted
