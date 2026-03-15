# Scripts Overview

This directory contains three script groups.

## Primary Entry Points

- `run_full_pipeline.py`
  - Main wrapper for the original `H3N2` direct-comparison pipeline.
- `run_h3n2_patristic_pipeline.py`
  - Parallel `H3N2` analysis that adds patristic distance alongside temporal distance in a separate output folder.
- `28_run_wic_ha1_model_pipeline.py`
  - Main wrapper for the `H3N2-WIC` HA1 modeling pipeline.
- `35_run_wic_ha1_patristic_pipeline.py`
  - Parallel `H3N2-WIC` HA1 analysis that adds patristic distance in a separate output folder.
- `36_run_wic_ha1_filtered_pipeline.py`
  - Parallel `H3N2-WIC` HA1 analysis that excludes egg, mixed, and unknown passage classes.
- `37_run_wic_ha1_filtered_patristic_pipeline.py`
  - Parallel `H3N2-WIC` HA1 analysis that excludes egg, mixed, and unknown passage classes and adds patristic distance.

## H3N2 Direct-Comparison Pipeline

- `01_fetch_sequences.py` to `11_compare_paper_sites.py`
  - Sequence fetch, alignment, strain matching, titer preprocessing, feature building,
    model training, SHAP analysis, figure generation, coverage audit, cleanup, and paper comparison.

## WIC Source-Assembly Pipeline

- `12_build_wic_manifest.py` to `23_build_wic_final_model_inputs.py`
  - Manifest construction, accession resolution, local GISAID/repo matching, duplicate handling,
    and construction of finalized WIC modeling inputs.

## WIC HA1 Modeling Pipeline

- `24_prepare_wic_ha1_alignment.py` to `32_cleanup_wic_model_outputs.py`
  - HA1 alignment prep, WIC titer preprocessing, feature building, model training,
    SHAP analysis, literature comparison, figure generation, and cleanup.
- `33_build_h3n2_patristic_tree.py`
  - Builds a FastTree Newick tree with branch lengths for the aligned `H3N2` HA proteins.
- `34_build_wic_ha1_patristic_tree.py`
  - Builds a FastTree Newick tree with branch lengths for the aligned `H3N2-WIC` HA1 sequences.

## Shared Utilities

- `common.py`
  - Shared config loading, additive-effect fitting, FASTA helpers, and model parameter helpers.
- `paper_sites.py`
  - Reference site lists for Neher 2016 and the WIC paper.
- `wic_name_utils.py`
  - WIC-specific strain/name normalization helpers.

The numbered scripts are stage implementations. For normal use, prefer the two wrapper scripts above.
