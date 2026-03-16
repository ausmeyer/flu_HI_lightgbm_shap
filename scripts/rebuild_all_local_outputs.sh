#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/Users/austinmeyer/numpy1_env/bin/python3}"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/flu_hi_rebuild.XXXXXX")"

cleanup_tmp() {
  rm -rf "$TMP_DIR"
}
trap cleanup_tmp EXIT

log() {
  printf '\n[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

run() {
  printf '+'
  printf ' %q' "$@"
  printf '\n'
  "$@"
}

require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    printf 'Required file not found: %s\n' "$path" >&2
    exit 1
  fi
}

snapshot_h3n2_alignment() {
  local source_dir="$ROOT_DIR/H3N2/output"
  local artifacts=(
    "H3N2_HA_aligned.fasta"
    "H3N2_alignment_position_map.csv"
    "H3N2_mature_position_map.csv"
  )

  for artifact in "${artifacts[@]}"; do
    require_file "$source_dir/$artifact"
    cp "$source_dir/$artifact" "$TMP_DIR/$artifact"
  done
}

restore_h3n2_alignment() {
  local target_dir="$ROOT_DIR/H3N2/output"
  local artifacts=(
    "H3N2_HA_aligned.fasta"
    "H3N2_alignment_position_map.csv"
    "H3N2_mature_position_map.csv"
  )

  mkdir -p "$target_dir/figures"
  for artifact in "${artifacts[@]}"; do
    cp "$TMP_DIR/$artifact" "$target_dir/$artifact"
  done
}

rebuild_h3n2_from_existing_alignment() {
  log "Rebuilding H3N2 from existing aligned artifacts"
  snapshot_h3n2_alignment
  rm -rf "$ROOT_DIR/H3N2/output"
  restore_h3n2_alignment

  run "$PYTHON_BIN" "$ROOT_DIR/scripts/03_match_strains.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/04_preprocess_titers.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/05_build_features.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/06_train_model.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/07_shap_analysis.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/08_generate_figures.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/11_compare_paper_sites.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/09_audit_strain_coverage.py" --config "$ROOT_DIR/configs/h3n2.json"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/10_cleanup_outputs.py" --config "$ROOT_DIR/configs/h3n2.json"
}

main() {
  log "Using Python: $PYTHON_BIN"

  rebuild_h3n2_from_existing_alignment

  log "Rebuilding H3N2-patristic"
  rm -rf "$ROOT_DIR/H3N2-patristic/output"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/run_h3n2_patristic_pipeline.py" \
    --config "$ROOT_DIR/configs/h3n2_patristic.json" \
    --skip-fetch

  log "Rebuilding base H3N2-WIC"
  rm -rf "$ROOT_DIR/H3N2-WIC/output/modeling"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/28_run_wic_ha1_model_pipeline.py" \
    --config "$ROOT_DIR/H3N2-WIC/config/wic_ha1_model.json"

  log "Rebuilding H3N2-WIC-no-egg-no-mixed-no-unknown"
  rm -rf "$ROOT_DIR/H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/36_run_wic_ha1_filtered_pipeline.py" \
    --config "$ROOT_DIR/H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown.json"

  log "Rebuilding H3N2-WIC-patristic"
  rm -rf "$ROOT_DIR/H3N2-WIC-patristic/output/modeling"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/35_run_wic_ha1_patristic_pipeline.py" \
    --config "$ROOT_DIR/H3N2-WIC/config/wic_ha1_patristic_model.json"

  log "Rebuilding H3N2-WIC-patristic-no-egg-no-mixed-no-unknown"
  rm -rf "$ROOT_DIR/H3N2-WIC-patristic-no-egg-no-mixed-no-unknown/output/modeling"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/37_run_wic_ha1_filtered_patristic_pipeline.py" \
    --config "$ROOT_DIR/H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown_patristic.json"

  log "Rebuilding integrated site comparison tables"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/38_build_master_site_model_comparison.py"

  log "Rebuilding manuscript cross-study figure"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/39_generate_cross_study_figure.py"

  log "Syncing manuscript figures into Overleaf repo"
  run "$PYTHON_BIN" "$ROOT_DIR/scripts/40_sync_manuscript_assets.py"

  log "Rebuild complete"
}

main "$@"
