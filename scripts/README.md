# Scripts Overview

This directory separates analysis stages, reporting/presentation stages, and copying into the writing projects.

## Primary Entry Points

- `run_full_pipeline.py`
  - Main wrapper for the original `H3N2` direct-comparison pipeline; its default path invokes 12 stages, including stage 43 and stage 11.
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
  - Figure generation includes signed top-20 individual covariate SHAP panels for
    the site-state and substitution models.
- `43_analyze_koel_context_dependence.py`
  - Existing supplemental-analysis source for Fig. S18; included in `run_full_pipeline.py`.

## WIC Source-Assembly Pipeline

- `12_build_wic_manifest.py` to `23_build_wic_final_model_inputs.py`
  - Manifest construction, accession resolution, local GISAID/repo matching, duplicate handling,
    and construction of finalized WIC modeling inputs.

## WIC HA1 Modeling Pipeline

- `24_prepare_wic_ha1_alignment.py` to `32_cleanup_wic_model_outputs.py`
  - HA1 alignment prep, WIC titer preprocessing, feature building, model training,
    SHAP analysis, literature comparison, figure generation, and cleanup.
  - Figure generation includes signed top-20 individual covariate SHAP panels for
    the site-state and substitution models.
- `44_analyze_wic_glycan_context_dependence.py`
  - Existing supplemental-analysis source for Fig. S19; included in the filtered WIC wrapper (36).
- `33_build_h3n2_patristic_tree.py`
  - Builds a FastTree Newick tree with branch lengths for the aligned `H3N2` HA proteins.
- `34_build_wic_ha1_patristic_tree.py`
  - Builds a FastTree Newick tree with branch lengths for the aligned `H3N2-WIC` HA1 sequences.

## Reporting and writing-project assets

- `38_build_master_site_model_comparison.py`: consolidates existing summary outputs into the root comparison TSVs.
- `39_generate_cross_study_figure.py`: creates the cross-study presentation figure in `manuscript/generated_figures/`; its internal `fig4` filename is copied as manuscript Fig. 3.
- `40_sync_manuscript_assets.py`: copies three existing figure PDFs into the main Overleaf clone, default `69b778fd3e7b181fe1c2943b`.
- `41_generate_supplement_tables.py`: creates Tables S1/S2/S6 and collects the reusable Tables S4/S5 TSVs in `manuscript/generated_supplement/`.
- `42_sync_supplement_assets.py`: copies the mapped figures and generated tables into the SI Overleaf clone, default `69b7797c762f515edcff3ad6`. Table S3 and `supplementary_data/` are outside this mapping.

Scripts 40/42 overwrite mapped destination files without a dry-run mode. Copying into the authoritative Overleaf projects is a release step, separate from analysis and presentation generation; it does not compile, commit, push, or update the dated snapshot. The SI README records the hand-maintained Table S3 sources and distinguishes accompanying data from compile dependencies.

`rebuild_all_local_outputs.sh` deletes all six analysis output directories before rebuilding them and then invokes stages 38–42, including both Overleaf copy operations. It preserves only three H3N2 alignment/position-map files temporarily, not a backup of the old outputs. Its cleanup runs by default. See the root README for the affected output families and external-tool requirements.

## External dependencies and inherited provenance

MAFFT is used by stages 02/24; FastTree is used by stages 33/34. R plus `readxl` is required by the exported-metadata helper (21). GISAIDR belongs to the optional R helper (20); the author confirmed manual GISAID search for the study, so script 20 is not evidence of historical GISAIDR execution. TeX/BibTeX/latexmk are document-build dependencies and are not supplied by `environment.yml`.

The filtered wrappers (36/37) copy the base WIC alignment, position map, extraction summary, and alignment metadata. The copied metadata retains the base source subtype and paths. It records the inherited alignment, not a separate filtered alignment; see `H3N2-WIC/README.md` for the source and destination paths.

## Shared Utilities

- `common.py`
  - Shared config loading, additive-effect fitting, FASTA helpers, and model parameter helpers.
- `paper_sites.py`
  - Reference site lists for Neher 2016 and the WIC paper.
- `wic_name_utils.py`
  - WIC-specific strain/name normalization helpers.
- `shap_reuse.py`
  - Value-specific SHAP summary helpers, including plotting for signed individual
    covariate SHAP panels.

The numbered scripts are stage implementations; the wrappers above define the existing analysis entry points. Documentation and document builds can be reviewed independently of those analyses.
