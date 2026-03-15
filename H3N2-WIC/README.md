# H3N2-WIC

WIC-specific files are organized so the raw inputs, working manifests, final modeling inputs, and reference material are separated.

There are two distinct WIC workflows in this repo:
- source assembly from raw paper/GISAID/repo inputs
- HA1 modeling from finalized WIC inputs

For publication-facing use, the main entry points are:
- [23_build_wic_final_model_inputs.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/23_build_wic_final_model_inputs.py) for finalized modeling inputs
- [28_run_wic_ha1_model_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/28_run_wic_ha1_model_pipeline.py) for HA1 model training, comparison, figures, and cleanup
- [35_run_wic_ha1_patristic_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/35_run_wic_ha1_patristic_pipeline.py) for the parallel HA1 analysis that adds patristic distance
- [36_run_wic_ha1_filtered_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/36_run_wic_ha1_filtered_pipeline.py) for the parallel HA1 analysis excluding egg, mixed, and unknown passage classes
- [37_run_wic_ha1_filtered_patristic_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/37_run_wic_ha1_filtered_patristic_pipeline.py) for the parallel HA1 analysis excluding egg, mixed, and unknown passage classes while adding patristic distance

## Layout

- `config/`
  - `h3n2_wic_config.json`: WIC config used by the source-prep scripts
- `data/raw/`
  - `WIC-HI-H3N2.csv`: raw WIC HI table
  - `gisaid/`: raw GISAID FASTA and exported metadata spreadsheets
  - `repo/`: published repo HA1 FASTA
- `data/intermediate/`
  - matched repo and GISAID FASTA files used during source assembly
- `data/final/`
  - `H3N2_WIC_HI_data_final.tsv`
  - `H3N2_WIC_seq_data_final.tsv`
  - `H3N2_WIC_sequences_final.fasta`
- `manifests/`
  - retained working manifests used for source assembly
  - `intermediate/` is used for non-final manifests during prep
- `metadata/`
  - retained lookup products, isolate-ID mappings, query files, and match summaries
  - `intermediate/` is used for non-final lookup artifacts during prep
- `references/`
  - paper PDFs and other reference material
- `output/`
  - audit tables, summaries, and other diagnostic outputs

## Main WIC Commands

Fresh repo-first source prep:

```bash
python3 scripts/17_prepare_wic_ha1_sources.py --download-repo-fasta --restart-clean
```

Post-GISAID source completion:

```bash
python3 scripts/19_finish_wic_ha1_sources.py \
  --config H3N2-WIC/config/h3n2_wic_config.json \
  --email austin.g.meyer@gmail.com \
  --max-strains 500
```

Resolve duplicate GISAID isolate IDs:

```bash
python3 scripts/22_resolve_wic_duplicate_isolates.py
```

Build final modeling inputs:

```bash
python3 scripts/23_build_wic_final_model_inputs.py
```

Run the HA1 modeling and publication-oriented output pipeline:

```bash
/Users/austinmeyer/numpy1_env/bin/python3 scripts/28_run_wic_ha1_model_pipeline.py \
  --config H3N2-WIC/config/wic_ha1_model.json
```

By default, this wrapper removes large regenerable intermediates at the end and keeps the
publication-oriented tables, metrics, figures, and methods-reporting outputs. To keep all
intermediate artifacts for debugging, add `--keep-intermediates`.

Parallel WIC analysis variants:

Patristic-distance variant:

```bash
/Users/austinmeyer/numpy1_env/bin/python3 scripts/35_run_wic_ha1_patristic_pipeline.py \
  --config H3N2-WIC/config/wic_ha1_patristic_model.json
```

Filtered-passage variant:

```bash
/Users/austinmeyer/numpy1_env/bin/python3 scripts/36_run_wic_ha1_filtered_pipeline.py \
  --config H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown.json
```

Filtered-passage + patristic-distance variant:

```bash
/Users/austinmeyer/numpy1_env/bin/python3 scripts/37_run_wic_ha1_filtered_patristic_pipeline.py \
  --config H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown_patristic.json
```

## Publication-Oriented Modeling Outputs

The default cleanup keeps only the outputs needed to support:
- the H3N2 vs WIC comparison story
- the Neher/Bedford and WIC-paper comparison tables
- methods reporting for the HA1 modeling workflow

Retained modeling outputs:
- `output/modeling/H3N2_WIC_HA1_HA1_aligned.fasta`
- `output/modeling/H3N2_WIC_HA1_ha1_alignment_metadata.json`
- `output/modeling/H3N2_WIC_HA1_ha1_extraction_summary.tsv`
- `output/modeling/H3N2_WIC_HA1_ha1_position_map.csv`
- `output/modeling/H3N2_WIC_HA1_matched_titers.tsv`
- `output/modeling/H3N2_WIC_HA1_cv_summary.tsv`
- `output/modeling/H3N2_WIC_HA1_neher_analog_metrics.tsv`
- `output/modeling/H3N2_WIC_HA1_paper_site_comparison.tsv`
- `output/modeling/H3N2_WIC_HA1_paper_site_validation.json`
- `output/modeling/H3N2_WIC_HA1_site_summary.tsv`
- `output/modeling/H3N2_WIC_HA1_site_stability.tsv`
- `output/modeling/H3N2_WIC_HA1_substitution_site_summary.tsv`
- `output/modeling/H3N2_WIC_HA1_substitution_site_stability.tsv`
- `output/modeling/figures/*.pdf`
- `output/modeling/H3N2_WIC_HA1_qc_log.txt`

Removed by default:
- large feature matrices
- prediction row tables
- redundant per-site rank tables when the same information exists in `*_site_summary.tsv`
- extra diagnostic figures not needed for the paper story
- fitted model text files and non-primary comparison tables

## Final Data

The directory `data/final/` contains only the final modeling data products:

- `H3N2_WIC_HI_data_final.tsv`
- `H3N2_WIC_seq_data_final.tsv`
- `H3N2_WIC_sequences_final.fasta`
