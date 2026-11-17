# Reproducibility package — CRC drug-response label audit

Manuscript: *An error budget for colorectal cancer pharmacogenomics: censoring, batch variation, and the limits of machine-learned drug resistance*

All analyses use public data. Nothing in `data/raw/` needs to be regenerated to reproduce
the tables, except where noted.

---

## 1. Environment

| | |
|---|---|
| Python | 3.13 |
| pandas / numpy / scipy | 3.1.1 / 2.5.4 / 1.18.1 |
| scikit-learn | 1.9.1 |
| statsmodels | 0.15.0 |
| matplotlib | 3.11.2 |
| pyarrow | 26.0.1 |
| requests | 2.34.2 |

```bash
pip install pandas numpy scipy scikit-learn statsmodels matplotlib pyarrow requests
```

## 2. Data acquisition

| File | Source |
|---|---|
| `GDSC2_fitted_dose_response_24Jul22.csv` | `https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/` |
| `GDSC1_fitted_dose_response_24Jul22.csv` | same |
| `GDSC2_public_raw_data_24Jul22.csv` (2.0 GB) | same |
| `Cell_Lines_Details.xlsx` | same |
| `CCLE_RNAseq_genes_rpkm_20180929.gct.gz` | `https://data.broadinstitute.org/ccle/` |
| `mutations_latest.csv.gz` | `https://cog.sanger.ac.uk/cmp/download/` |

Patient cohorts are pulled live from the cBioPortal REST API v2
(`coad_silu_2022`, `coadread_tcga`) and are not stored in this package.

**Access notes (2026-11-16).** NCBI GEO / FTP was unreachable from the execution host, so
GSE39582 and GSE119409 could not be used. DepMap portal and figshare sit behind a
verification wall; Sanger CMP RNA-seq files return 403. Broad's legacy CCLE directory and
Sanger's GDSC FTP were open.

## 3. Script order

Run from the repository root.

| # | Script | Needs raw data | Output |
|---|---|---|---|
| 01 | `scripts/01_label_forensics.py` | fitted | `table1_drug_panel.csv`, `table2_replicate_reliability.csv` |
| 02 | `scripts/02_cross_campaign.py` | fitted ×2 | `table3_cross_campaign.csv`, `coread_ic50_two_campaigns.csv` |
| 03 | `scripts/03_reliability_by_identifiability.py` | fitted ×2 | `03_reliability_stratified.csv`, `coread_pairs_with_identifiability.csv` |
| 04 | `scripts/04_extract_raw_coread.py` | **raw 2.0 GB** | `data/processed/coread_raw_wells.parquet` |
| 05 | `scripts/05_model_vs_ceiling.py` | fitted + CCLE | `05_model_vs_ceiling.csv` |
| 06 | `scripts/06_decisive_test.py` | fitted + CCLE | `06_decisive_test.csv` |
| 08 | `scripts/08_build_patient_score.py` | fitted + CCLE + network | `patient_scores.csv` |
| 09 | `scripts/09_interaction_test.py` | `patient_scores.csv` + network | `table4_interaction.csv`, `table5_stratified_event_rates.csv` |
| 10 | `scripts/10_auc_vs_ic50.py` | fitted ×2 + CCLE | `table6_auc_vs_ic50.csv` |
| 11 | `scripts/11_well_level_noise.py` | **parquet from 04** | `table7_noise_floor.csv`, `11_well_noise.json` |
| 07 | `scripts/07_figure1.py` | fitted + tables | `figure1_forensics.png` |
| 12 | `scripts/12_figures.py` | tables | `figure2/3/4*.png` |

Standalone tool:

```bash
python preflight_audit.py \
  --fitted  data/raw/GDSC2_fitted_dose_response_24Jul22.csv \
  --fitted2 data/raw/GDSC1_fitted_dose_response_24Jul22.csv \
  --disease COREAD \
  --out results/audit_COREAD
```

## 4. Determinism

All models use fixed seeds. Two known sources of run-to-run variation:

1. **Sample ordering.** R² for 5-FU varies between 0.29 and 0.36 because
   `drop_duplicates` preserves input order and therefore changes fold assignment under a
   fixed seed. This is reported in the manuscript rather than tuned away.
2. **`cKDTree`-free elastic net is deterministic**; the variation above is purely from
   fold composition.

## 5. Key numbers and where they come from

| Claim | Value | Source file |
|---|---|---|
| COREAD curves / lines | 12,538 / 46 | `01_label_forensics.json` |
| IC₅₀ censored above tested range | 73.0% | `01_label_forensics.json` |
| Median censoring across 286 compounds | 88.9%; 95 compounds at 100% | `audit_COREAD_by_drug.csv` |
| Replicate-well CV | 9.4% | `11_well_noise.json` |
| Single-dose SNR | 1.54; 74.4% below 2 | `11_well_noise.json` |
| GDSC1 vs GDSC2, all pairs | r = 0.811, 3.0× median | `02_cross_campaign.json` |
| GDSC1 vs GDSC2, 5-FU | r = 0.667, 3.5× median | `table3_cross_campaign.csv` |
| Jointly identifiable pairs | 18.0% (5-FU: 32.6%, r = 0.048) | `03_reliability_stratified.csv` |
| 5-FU model R² | 0.29–0.36 | `05_model_vs_ceiling.csv` |
| 5-FU censoring-boundary AUC | 0.778 | manuscript Results, §3 decomposition |
| Patient interaction HR | 1.08 (0.66–1.76), p = 0.76 | `table4_interaction.csv` |
| AUC target vs IC₅₀ target | 1/8 vs 3/8 beat null | `table6_auc_vs_ic50.csv` |

## 6. Known limitations of the package

- 7 of 46 cell lines could not be matched to CCLE (listed in Additional file 1).
- The patient layer rests on one cohort (AC-ICAM). TCGA-COADREAD has an adjuvant
  chemotherapy field populated in only 36/636 patients.
- No wet-lab validation.
