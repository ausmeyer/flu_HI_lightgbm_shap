# LightGBM HI Antigenic Model Pipeline

End-to-end pipeline for:
- the original `H3N2` Neher/Bedford-style direct-comparison analysis
- the `H3N2-WIC` data-assembly workflow
- the `H3N2-WIC` HA1 modeling and comparison workflow

The repo now has two publication-facing entry points:
- [run_full_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/run_full_pipeline.py) for the original `H3N2` pipeline
- [28_run_wic_ha1_model_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/28_run_wic_ha1_model_pipeline.py) for the `H3N2-WIC` HA1 modeling pipeline

Additional parallel analysis entry points:
- [run_h3n2_patristic_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/run_h3n2_patristic_pipeline.py) for `H3N2` with patristic distance added alongside temporal distance
- [35_run_wic_ha1_patristic_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/35_run_wic_ha1_patristic_pipeline.py) for `H3N2-WIC` with patristic distance added alongside temporal distance
- [36_run_wic_ha1_filtered_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/36_run_wic_ha1_filtered_pipeline.py) for `H3N2-WIC` excluding egg, mixed, and unknown passage classes
- [37_run_wic_ha1_filtered_patristic_pipeline.py](/Users/austinmeyer/My%20Drive/Research/Faculty/flu_HI_lightgbm_shap/scripts/37_run_wic_ha1_filtered_patristic_pipeline.py) for `H3N2-WIC` excluding egg, mixed, and unknown passage classes while adding patristic distance

The numbered scripts remain as implementation stages, but the wrappers above should be treated as the main human-facing commands.

## Repository layout

- `H3N2/H3N2_HI_data.tsv`
- `H3N2/H3N2_seq_data.tsv`
- `H1N1/H1N1_HI_data.tsv`
- `scripts/`
  - `01-11`: original `H3N2` pipeline stages
  - `12-23`: WIC source-assembly stages
  - `24-32`: WIC HA1 modeling stages
  - `run_full_pipeline.py`: primary `H3N2` wrapper
  - `28_run_wic_ha1_model_pipeline.py`: primary WIC HA1 wrapper
- `scripts/common.py`
- `scripts/paper_sites.py`
- `scripts/wic_name_utils.py`
- `configs/h3n2.json`
- `configs/h1n1_template.json`
- `environment.yml`

## Setup

```bash
conda env create -f environment.yml
conda activate flu_hi_lgbm
```

## Pipeline execution (H3N2)

Full run:

```bash
python scripts/run_full_pipeline.py --config configs/h3n2.json --email your_email@example.com
```

This wrapper executes the same numbered steps below.
By default it removes intermediates at the end and keeps only the publication-oriented outputs.

1) Prepare accession lists and fetch GenBank (IRD) sequences:

```bash
python scripts/01_fetch_sequences.py --config configs/h3n2.json --email your_email@example.com
```

This creates:
- `H3N2/output/H3N2_gisaid_accessions_to_download.txt`
- `H3N2/genbank_sequences.fasta`

2) Manual GISAID step:
- Use `H3N2/output/H3N2_gisaid_accessions_to_download.txt` in GISAID EpiFlu search/upload.
- Download HA protein FASTA.
- Save it as `H3N2/gisaid_sequences.fasta`.

3) Merge and align sequences + trim signal peptide + QC:

```bash
python scripts/02_merge_and_align.py --config configs/h3n2.json
```

4) Match HI strains to sequence strains:

```bash
python scripts/03_match_strains.py --config configs/h3n2.json
```

5) Titer preprocessing + paper-style standardization:

```bash
python scripts/04_preprocess_titers.py --config configs/h3n2.json
```

This creates:
- `homologous_log2_titer`
- `homologous_titer_source`
- `standardized_titer`
- `serum_potency`
- `virus_avidity`
- `corrected_titer`

6) Build feature matrices:

```bash
python scripts/05_build_features.py --config configs/h3n2.json
```

This creates:
- binary site-change features: `change_1..change_550`
- explicit virus amino-acid covariates like `K145`
- explicit serum amino-acid covariates like `145K`
- substitution-identity features: `sub_<site>_<serumAA>_<virusAA>`
- alignment numbering summaries for direct comparison

7) Train grouped LightGBM CV + direct-comparison baselines:

```bash
python scripts/06_train_model.py --config configs/h3n2.json
```

This fits:
- additive-only baseline
- temporal-only baseline
- randomized binary baseline
- site-state LightGBM using `change`, `virus_aa`, and `serum_aa` at each site
- binary site-change LightGBM
- substitution-identity LightGBM

8) Cross-validated SHAP analysis + Koel-site validation:

```bash
python scripts/07_shap_analysis.py --config configs/h3n2.json
```

9) Generate figures:

```bash
python scripts/08_generate_figures.py --config configs/h3n2.json
```

10) Audit coverage against the HI table:

```bash
python scripts/09_audit_strain_coverage.py --config configs/h3n2.json
```

This creates:
- `H3N2/output/H3N2_missing_from_seq_manifest.tsv`
- `H3N2/output/H3N2_present_but_failed_qc.tsv`

11) Remove intermediate artifacts and keep only final outputs:

```bash
python scripts/10_cleanup_outputs.py --config configs/h3n2.json
```

## Final retained outputs

The default cleanup now keeps only the files needed for:
- the comparison story
- methods reporting
- the minimum set of publication figures

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
- `H3N2/output/figures/*.pdf`
- `H3N2/output/H3N2_qc_log.txt`

The cleanup intentionally removes larger or redundant analysis products such as:
- raw SHAP matrices
- feature-level SHAP tables
- duplicate rank tables when the same information exists in `*_site_summary.tsv`
- extra diagnostic figures not needed for the paper narrative

## Notes

- MAFFT must be available on `PATH`.
- `01_fetch_sequences.py` requires an NCBI email for Entrez.
- `lab` is derived from `source` with categories: `nimr_crick`, `cdc`, `vidrl`, `niid`, `other`.
- The direct-comparison target is `standardized_titer = homologous_log2_titer - log2_titer`, with homologous fallback controlled by config.
- Cross-validation now refits additive serum/virus effects within each split to avoid leakage.
- H1N1 extension: copy `configs/h1n1_template.json`, set reference + subtype-specific fields, and run the same scripts.
