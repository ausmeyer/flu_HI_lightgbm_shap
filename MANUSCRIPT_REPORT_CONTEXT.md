# Manuscript Report Context

This file is a manuscript-facing context pack assembled from pipeline inputs, retained intermediates, generated outputs, and QC logs. It is intended to provide the factual reporting context needed to write or audit the manuscript without pulling narrative language from the manuscript source.

## Provenance Rules

- This file is derived from pipeline artifacts, not from `69b778fd3e7b181fe1c2943b/` or `69b7797c762f515edcff3ad6/` manuscript prose.
- When a value appears in multiple places, prefer the retained pipeline artifact or consolidated table listed here.
- For site-level comparisons, the main master table is `site_model_literature_comparison.tsv`.
- For model metrics, the main consolidated table is `model_performance_literature_summary.tsv`.
- For manuscript-facing final tables and figure support files, prefer `manuscript/generated_supplement/` and `manuscript/generated_figures/`.

## Primary Consolidated Tables

- `site_model_literature_comparison.tsv`
  - Master site-level comparison table.
  - Contains benchmark-set membership flags, external sequence-only rankings, and all model-derived site ranks for:
    - H3N2 site-state and substitution models
    - H3N2 patristic variants
    - WIC unfiltered site-state and substitution models
    - WIC filtered site-state and substitution models
    - WIC patristic variants
    - external sequence-only rankings and dN/dS rankings
- `model_performance_literature_summary.tsv`
  - Consolidated model metrics table.
  - Contains RMSE, MAE, R2, Spearman, and external sequence-only model summary metrics.
- `manuscript/generated_supplement/table_s1_model_performance.tsv`
  - Manuscript-facing subset of primary model performance.
- `manuscript/generated_supplement/table_s2_reference_overlap.tsv`
  - Manuscript-facing benchmark overlap counts for the primary model families.
- `manuscript/generated_figures/fig4_panelB_overlaps.tsv`
  - Manuscript-facing overlap counts for the cross-study concordance figure, including external sequence-only rankings.

## Pipeline Entry Points

- Original H3N2 direct-comparison pipeline:
  - `scripts/run_full_pipeline.py --config configs/h3n2.json`
- WIC HA1 modeling pipeline:
  - `scripts/28_run_wic_ha1_model_pipeline.py --config H3N2-WIC/config/wic_ha1_model.json`
- Parallel analysis variants:
  - `scripts/run_h3n2_patristic_pipeline.py`
  - `scripts/35_run_wic_ha1_patristic_pipeline.py`
  - `scripts/36_run_wic_ha1_filtered_pipeline.py`
  - `scripts/37_run_wic_ha1_filtered_patristic_pipeline.py`

## Analysis Families Used in the Paper

### 1. H3N2 Neher/Bedford-style dataset

- Config: `configs/h3n2.json`
- Raw inputs:
  - `H3N2/H3N2_HI_data.tsv`
  - `H3N2/H3N2_seq_data.tsv`
- Reference strain: `A/HongKong/1/1968`
- Sequence trimming/QC constants:
  - signal peptide length: `16`
  - HA1 length: `328`
  - mature HA length: `550`
  - `max_gap_fraction = 0.05`
  - `max_ambiguous_fraction = 0.02`
  - `min_ungapped_length = 500`
- Modeling target:
  - `comparison_target_column = standardized_titer`
- Cross-validation:
  - `n_splits = 5`
  - `random_state = 42`

### 2. WIC H3N2 dataset, unfiltered

- Config: `H3N2-WIC/config/wic_ha1_model.json`
- Final modeling inputs:
  - `H3N2-WIC/data/final/H3N2_WIC_HI_data_final.tsv`
  - `H3N2-WIC/data/final/H3N2_WIC_seq_data_final.tsv`
  - `H3N2-WIC/data/final/H3N2_WIC_sequences_final.fasta`
- Output directory:
  - `H3N2-WIC/output/modeling`
- Cross-validation:
  - `n_splits = 5`
  - `random_state = 42`

### 3. WIC H3N2 dataset, filtered to exclude egg/mixed/unknown passage

- Config: `H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown.json`
- Output directory:
  - `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling`
- Excluded passage classes:
  - `EGG`
  - `MIXED`
  - `UNKNOWN`

### 4. Patristic-distance variants

- H3N2 patristic config adds:
  - `distance_feature_cols = ["temporal_distance", "patristic_distance"]`
- WIC patristic config adds:
  - `distance_feature_cols = ["temporal_distance", "patristic_distance"]`

## Reference Site Sets Used for Comparison

Definitions come from `scripts/paper_sites.py` and `configs/h3n2.json`.

- Koel 7 sites:
  - `145, 155, 156, 158, 159, 189, 193`
- Neher/Bedford 15-site reference set:
  - `62, 121, 135, 140, 144, 145, 155, 156, 158, 159, 186, 189, 193, 196, 276`
- Harvey/WIC 15-site reference set:
  - `53, 121, 126, 131, 135, 137, 144, 145, 157, 158, 159, 160, 173, 189, 193`
- Shah/WIC 30-site reference set:
  - `45, 62, 122, 131, 135, 138, 140, 142, 144, 145, 158, 159, 171, 173, 183, 186, 189, 193, 194, 196, 197, 208, 213, 219, 223, 225, 241, 261, 269, 311`

## H3N2 Dataset Assembly and QC Summary

Source: `H3N2/output/H3N2_qc_log.txt`

### Match and preprocessing

- HI rows input: `10,059`
- Rows retained after virus/serum matching: `9,024`
- Rows dropped at matching stage: `1,035`
- Unique virus strains in HI: `402`
- Unique serum strains in HI: `215`
- Unique matched virus strains: `402`
- Unique matched serum strains: `164`
- Serum fuzzy matches flagged for review: `425`
- Review file:
  - `H3N2/output/H3N2_match_review.tsv`

### Titer preprocessing

- Censored `<X` titers handled: `1,180`
- Homologous proxy source counts:
  - `homologous = 7,880`
  - `max_observed = 1,144`
- Standardized titer mean/std:
  - mean `2.607678`
  - std `2.839273`
- Additive correction:
  - intercept `1.103580`
  - serum effects estimated for `437` sera
  - virus effects estimated for `402` viruses
- Corrected titer std:
  - `1.756856`
- Derived lab counts:
  - `other = 7,268`
  - `nimr_crick = 1,756`

### Feature construction

- Mature positions used: `550`
- Site-state feature matrix shape: `(9,024, 2,386)`
- Primary site-state covariates: `2,365`
- Distance feature columns: `temporal_distance`
- Substitution feature matrix shape: `(9,024, 973)`
- Observed substitution features: `952`
- Change-matrix missing fraction: `0.154799`

### Primary H3N2 site-state model performance

- RMSE: `2.0586 +/- 0.2007`
- MAE: `1.6024 +/- 0.1617`
- R2: `0.4719 +/- 0.0866`
- Spearman: `0.7076 +/- 0.0395`

### H3N2 SHAP/benchmark summary

- SHAP rows analyzed: `9,024`
- SHAP feature count: `2,366`
- Koel median rank: `7.0`
- Koel one-sided Mann-Whitney p-value: `6.607471307160945e-11`

### H3N2 coverage audit

- Unique HI strains across virus or serum columns: `465`
- Unique strains missing from sequence manifest: `63`
- Unique strains present in manifest but absent after QC: `0`

## H3N2 Patristic Variant Summary

Source: `H3N2-patristic/output/H3N2_PATRISTIC_qc_log.txt`

- Tree leaf count: `402`
- Branches with lengths: `771`
- Total branch length: `1.534635`
- FastTree model: `WAG`
- Site-state performance:
  - RMSE `2.0540 +/- 0.2079`
  - Spearman `0.7082 +/- 0.0404`
- Main conclusion for manuscript reporting:
  - Adding patristic distance gives at most marginal performance changes relative to the temporal-distance-only H3N2 model.

## WIC Source-Assembly Summary

Sources:
- `H3N2-WIC/metadata/H3N2_WIC_manifest_summary.tsv`
- `H3N2-WIC/metadata/H3N2_WIC_repo_match_summary.tsv`
- `H3N2-WIC/metadata/H3N2_WIC_gisaid_match_summary.tsv`
- `H3N2-WIC/output/H3N2_WIC_gisaid_duplicate_resolution_summary.tsv`
- `H3N2-WIC/output/H3N2_WIC_final_input_summary.tsv`

### WIC manifest scale

- Input rows in combined WIC manifest: `202,605`
- Unique virus strains: `16,213`
- Unique reference strains: `381`
- Unique strains in union: `16,223`
- Virus-only strains: `15,842`
- Reference-only strains: `10`
- Virus-and-reference strains: `371`
- Rows with inferred year: `15,867`
- Rows missing year: `356`

### Sequence matching against source repositories

- Published repo FASTA sequences: `1,739`
- Repo-matched rows: `1,817`
- Repo-matched distinct strains: `1,817`
- Repo-unmatched rows: `14,406`
- Repo-unmatched distinct lookup strains: `14,353`

- GISAID FASTA sequences: `4,212`
- GISAID-matched rows: `2,862`
- GISAID-matched distinct strains: `2,862`

### GISAID duplicate resolution

- Matched query strains: `2,823`
- Total candidate rows: `4,214`
- Query strains with multiple IDs: `909`
- Duplicate groups with identical sequences: `478`
- Duplicate groups varying only in partial sequence status: `160`
- Duplicate groups with multiple full-length variants: `419`
- Canonical selected as HA1-complete: `2,785`
- Canonical selected as longest available: `38`
- Canonical review-flagged: `364`
- Canonical selected passage labels:
  - `ORIGINAL = 477`
  - `CELL = 2,102`
  - `UNKNOWN = 212`
  - `EGG = 31`
  - `MIXED = 1`

### Final WIC modeling inputs

- Input HI rows before final sequence join: `202,605`
- Retained HI rows in final WIC input tables: `104,972`
- Retained fraction of input HI rows: `0.5181115964561586`
- Retained unique virus strains: `4,346`
- Retained unique serum strains: `267`
- Retained unique sequence manifest rows: `4,239`
- Retained unique FASTA IDs: `4,239`
- Final retained sequence provenance:
  - GISAID manifest rows: `2,501`
  - paper repo rows: `1,738`
  - IRD rows: `0`
  - NCBI rows: `0`

## WIC Unfiltered HA1 Modeling Summary

Source: `H3N2-WIC/output/modeling/H3N2_WIC_HA1_qc_log.txt`

### Alignment and preprocessing

- Final FASTA sequences entering HA1 extraction: `4,239`
- Retained HA1 sequences: `4,201`
- HA1 alignment reference: `A_BANGKOK_139_1990`
- Input final HI rows before passage filtering: `104,972`
- Rows retained after passage filtering: `104,972`
- Excluded passage classes: none
- Homologous proxy counts:
  - `homologous = 103,478`
  - `max_observed = 1,494`
- Standardized titer mean/std:
  - mean `3.640561`
  - std `2.076554`
- Context additive intercept: `-0.621606`
- Context-corrected titer std: `1.528194`

### Feature construction

- Rows retained for HA1 modeling after alignment match: `104,605`
- Rows dropped for missing aligned HA1 sequence: `367`
- Site-state feature matrix shape: `(104,605, 1,960)`
- Context feature columns: `14`
- Distance feature columns: `temporal_distance`
- Primary site-state covariates: `1,934`
- Observed substitution features: `1,239`
- Change-matrix missing fraction: `0.005564`

### Model setup

- Feature rows: `104,605`
- Unique virus-strain groups in CV: `4,308`
- Context columns: `10`
- Binary feature columns: `338`
- Site-state feature columns: `1,944`
- Substitution feature columns: `1,249`
- Important modeling note:
  - The WIC primary models are sequence-plus-context models and do not use virus or serum identity terms as predictors.

## WIC Filtered HA1 Modeling Summary

Source: `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_qc_log.txt`

### Passage filtering and preprocessing

- Input final HI rows before passage filtering: `104,972`
- Rows retained after passage filtering: `28,521`
- Rows excluded by passage filter: `76,451`
- Excluded passage classes:
  - `EGG`
  - `MIXED`
  - `UNKNOWN`
- Rows retained after alignment match for final HA1 modeling: `28,448`
- Rows dropped for missing aligned HA1 sequence: `73`
- Normalized bare `<` censored titers to `<40`: `1,470`
- Homologous proxy counts:
  - `homologous = 26,981`
  - `max_observed = 1,540`
- Standardized titer mean/std:
  - mean `2.813541`
  - std `1.653171`
- Context additive intercept: `3.354401`
- Context-corrected titer std: `1.114768`

### Feature construction

- Aligned HA1 sequences available: `4,201`
- Site-state feature matrix shape: `(28,448, 1,667)`
- Context feature columns: `14`
- Distance feature columns: `temporal_distance`
- Primary site-state covariates: `1,641`
- Observed substitution features: `702`
- Change-matrix missing fraction: `0.008644`

### Model setup

- Feature rows: `28,448`
- Unique virus-strain groups in CV: `3,432`
- Context columns: `10`
- Binary feature columns: `338`
- Site-state feature columns: `1,651`
- Substitution feature columns: `712`
- Important modeling note:
  - The filtered WIC primary models are also sequence-plus-context models with no virus or serum identity terms as predictors.

## Primary Model Performance Table

Preferred source:
- `model_performance_literature_summary.tsv`
- manuscript-facing subset: `manuscript/generated_supplement/table_s1_model_performance.tsv`

| Analysis label | Model family | RMSE | Spearman |
| --- | --- | ---: | ---: |
| Neher/Bedford data | site-state | 2.058606 | 0.707594 |
| Neher/Bedford data | substitution | 2.127685 | 0.694571 |
| Neher/Bedford data + patristic | site-state | 2.053964 | 0.708240 |
| Neher/Bedford data + patristic | substitution | 2.132453 | 0.693173 |
| Harvey/WIC data | site-state | 1.213688 | 0.822040 |
| Harvey/WIC data | substitution | 1.249554 | 0.811072 |
| Harvey/WIC filtered data | site-state | 0.974737 | 0.789671 |
| Harvey/WIC filtered data | substitution | 0.996947 | 0.780296 |
| Harvey/WIC data + patristic | site-state | 1.217008 | 0.821881 |
| Harvey/WIC data + patristic | substitution | 1.240341 | 0.814406 |
| Harvey/WIC filtered data + patristic | site-state | 0.964145 | 0.794892 |
| Harvey/WIC filtered data + patristic | substitution | 0.995371 | 0.780544 |

## Primary Benchmark Overlap Table

Preferred source:
- `manuscript/generated_supplement/table_s2_reference_overlap.tsv`

Columns are overlap counts with:
- Koel 7
- Neher/Bedford 15
- Harvey/WIC 15
- Shah/WIC 30

| Analysis label | Model family | Koel | Neher/Bedford | Harvey/WIC | Shah/WIC |
| --- | --- | --- | --- | --- | --- |
| Neher/Bedford data | site-state | `6/7` | `11/15` | `11/15` | `12/30` |
| Neher/Bedford data | substitution | `7/7` | `13/15` | `13/15` | `14/30` |
| Neher/Bedford data + patristic | site-state | `6/7` | `10/15` | `10/15` | `13/30` |
| Neher/Bedford data + patristic | substitution | `7/7` | `13/15` | `10/15` | `15/30` |
| Harvey/WIC data | site-state | `6/7` | `10/15` | `11/15` | `18/30` |
| Harvey/WIC data | substitution | `6/7` | `13/15` | `12/15` | `19/30` |
| Harvey/WIC filtered data | site-state | `5/7` | `10/15` | `12/15` | `18/30` |
| Harvey/WIC filtered data | substitution | `4/7` | `9/15` | `12/15` | `17/30` |
| Harvey/WIC data + patristic | site-state | `6/7` | `11/15` | `10/15` | `18/30` |
| Harvey/WIC data + patristic | substitution | `6/7` | `12/15` | `10/15` | `20/30` |
| Harvey/WIC filtered data + patristic | site-state | `5/7` | `8/15` | `10/15` | `17/30` |
| Harvey/WIC filtered data + patristic | substitution | `5/7` | `9/15` | `11/15` | `18/30` |

## Top-Ranked Sites in the Primary Models

Preferred full sources:
- `H3N2/output/H3N2_site_summary.tsv`
- `H3N2/output/H3N2_substitution_site_summary.tsv`
- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_site_summary.tsv`
- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_substitution_site_summary.tsv`
- master site table: `site_model_literature_comparison.tsv`

### H3N2 site-state top 15

`189, 133, 135, 145, 121, 193, 156, 144, 278, 158, 225, 159, 262, 62, 226`

### H3N2 substitution top 15

`158, 189, 186, 145, 156, 225, 133, 135, 193, 157, 62, 159, 227, 190, 226`

### Filtered WIC site-state top 15

`144, 225, 159, 140, 145, 142, 62, 53, 45, 197, 326, 193, 158, 173, 189`

### Filtered WIC substitution top 15

`225, 145, 159, 144, 138, 140, 312, 142, 53, 326, 197, 193, 189, 62, 128`

## External Sequence-Only Comparisons Used in the Paper

Primary sources:
- `site_model_literature_comparison.tsv`
- `manuscript/generated_figures/fig4_panelB_overlaps.tsv`
- `scripts/39_generate_cross_study_figure.py`

The cross-study figure compares the top-30 site-state rankings from:
- H3N2 site-state model
- filtered WIC site-state model

against:
- Koel benchmark set
- Neher/Bedford benchmark set
- Harvey/WIC benchmark set
- Shah/WIC benchmark set
- LightGBM+SHAP passage classification ranking
- LightGBM+SHAP collection-date ranking
- unpassaged sequence entropy ranking
- unpassaged dN/dS ranking
- McWhite unpassaged dN/dS ranking

### Cross-study overlap counts used in Fig. 4B

Source: `manuscript/generated_figures/fig4_panelB_overlaps.tsv`

| Reference ranking or set | H3N2 top-30 overlap | Filtered WIC top-30 overlap | Reference size |
| --- | ---: | ---: | ---: |
| Koel sites | 6 | 5 | 7 |
| Neher/Bedford sites | 11 | 10 | 15 |
| Harvey/WIC sites | 11 | 12 | 15 |
| Shah/WIC sites | 12 | 18 | 30 |
| LightGBM+SHAP passage | 11 | 16 | 30 |
| LightGBM+SHAP date | 11 | 19 | 30 |
| Unpassaged entropy | 9 | 21 | 30 |
| Unpassaged dN/dS | 7 | 9 | 30 |
| McWhite unpassaged dN/dS | 4 | 6 | 30 |

### Sequence-only model metrics

Source: `model_performance_literature_summary.tsv`

- Passage classification on 39,121 HA sequences:
  - overall accuracy `0.81`
  - balanced accuracy `0.77`
- Collection-date regression on 20,109 unpassaged HA sequences:
  - MAE `74.5` days
  - R2 `0.98`

## Manuscript-Facing Generated Assets

### Figures

- `manuscript/generated_figures/fig4_cross_study_concordance.pdf`
- `manuscript/generated_figures/fig4_cross_study_concordance.png`
- `manuscript/generated_figures/fig4_panelA_sites.tsv`
- `manuscript/generated_figures/fig4_panelB_overlaps.tsv`

### Supplement tables

- `manuscript/generated_supplement/table_s1_model_performance.tex`
- `manuscript/generated_supplement/table_s1_model_performance.tsv`
- `manuscript/generated_supplement/table_s2_reference_overlap.tex`
- `manuscript/generated_supplement/table_s2_reference_overlap.tsv`

## Retained Pipeline Outputs Most Relevant for Manuscript Assembly

### H3N2 retained outputs

- `H3N2/output/H3N2_HA_aligned.fasta`
- `H3N2/output/H3N2_alignment_position_map.csv`
- `H3N2/output/H3N2_cv_summary.csv`
- `H3N2/output/H3N2_koel_validation.json`
- `H3N2/output/H3N2_match_review.tsv`
- `H3N2/output/H3N2_mature_position_map.csv`
- `H3N2/output/H3N2_missing_from_seq_manifest.tsv`
- `H3N2/output/H3N2_neher_analog_metrics.csv`
- `H3N2/output/H3N2_paper_site_comparison.tsv`
- `H3N2/output/H3N2_present_but_failed_qc.tsv`
- `H3N2/output/H3N2_site_stability.tsv`
- `H3N2/output/H3N2_site_summary.tsv`
- `H3N2/output/H3N2_substitution_site_stability.tsv`
- `H3N2/output/H3N2_substitution_site_summary.tsv`
- `H3N2/output/H3N2_qc_log.txt`

### WIC retained outputs

- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_HA1_aligned.fasta`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_ha1_alignment_metadata.json`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_ha1_extraction_summary.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_ha1_position_map.csv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_matched_titers.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_cv_summary.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_neher_analog_metrics.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_paper_site_comparison.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_paper_site_validation.json`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_site_summary.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_site_stability.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_substitution_site_summary.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_substitution_site_stability.tsv`
- `H3N2-WIC/output/modeling/H3N2_WIC_HA1_qc_log.txt`

### Filtered WIC retained outputs

- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_cv_summary.tsv`
- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_paper_site_comparison.tsv`
- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_site_summary.tsv`
- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_substitution_site_summary.tsv`
- `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_qc_log.txt`

## Source Index for Rebuilding This Report

- Repo overview:
  - `README.md`
  - `H3N2-WIC/README.md`
- Configs:
  - `configs/h3n2.json`
  - `H3N2-WIC/config/wic_ha1_model.json`
  - `H3N2-WIC/config/wic_ha1_patristic_model.json`
  - `H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown.json`
  - `H3N2-WIC/config/wic_ha1_no_egg_mixed_unknown_patristic.json`
- Reference set definitions:
  - `scripts/paper_sites.py`
- Cross-study figure logic:
  - `scripts/39_generate_cross_study_figure.py`
- H3N2 QC summary:
  - `H3N2/output/H3N2_qc_log.txt`
- H3N2 patristic QC summary:
  - `H3N2-patristic/output/H3N2_PATRISTIC_qc_log.txt`
- WIC source-assembly summaries:
  - `H3N2-WIC/metadata/H3N2_WIC_manifest_summary.tsv`
  - `H3N2-WIC/metadata/H3N2_WIC_repo_match_summary.tsv`
  - `H3N2-WIC/metadata/H3N2_WIC_gisaid_match_summary.tsv`
  - `H3N2-WIC/output/H3N2_WIC_gisaid_duplicate_resolution_summary.tsv`
  - `H3N2-WIC/output/H3N2_WIC_final_input_summary.tsv`
- WIC modeling QC summaries:
  - `H3N2-WIC/output/modeling/H3N2_WIC_HA1_qc_log.txt`
  - `H3N2-WIC-no-egg-no-mixed-no-unknown/output/modeling/H3N2_WIC_HA1_NO_EGG_MIXED_UNKNOWN_qc_log.txt`
- Consolidated model/result tables:
  - `model_performance_literature_summary.tsv`
  - `site_model_literature_comparison.tsv`
  - `manuscript/generated_supplement/table_s1_model_performance.tsv`
  - `manuscript/generated_supplement/table_s2_reference_overlap.tsv`
  - `manuscript/generated_figures/fig4_panelB_overlaps.tsv`

## Usage Note

If a future manuscript edit needs a number, site list, overlap count, or provenance statement, first look in this file and then resolve to the listed source artifact rather than the manuscript text. This keeps the reporting chain anchored to pipeline outputs.
