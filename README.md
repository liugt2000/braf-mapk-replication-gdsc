# braf-mapk-replication-gdsc — replication-first pharmacogenomics of the BRAF–MEK axis

Reproducibility package for the manuscript *"Replication-first machine learning in
pharmacogenomics: of 1,338 gene–drug pairs only the BRAF–MEK axis replicates, with a
drug-class dissociation in colorectal cancer"* (submitted to the *Journal of
Translational Medicine* special issue "Harnessing Innovative Machine Learning
Techniques to Combat Drug Resistance in Solid Tumors").

DOI badge: *(added after Zenodo archiving)*

**One-line summary.** Of 1,338 (driver gene, drug) pairs testable across two independent
GDSC screening campaigns, five replicate — all BRAF → BRAF/MEK inhibitor. What is
reproducible in colorectal cancer is not a transcriptional programme but a drug-class
dissociation: BRAF-selective inhibition is 30-fold weaker than in melanoma, MEK inhibition
is equipotent. The association and the dissociation both survive growth-rate correction
and reproduce in two third-party campaigns (CTRPv2, PRISM).

---

## Environment

- Python 3.13; numpy, pandas, scipy, scikit-learn, statsmodels, lifelines, matplotlib,
  pyarrow, requests, rdata 1.1.0 (PharmacoSet parsing for the third-party replication)
- No GPU, no cluster. Peak memory < 6 GB. Full pipeline ≈ 30 min (the third-party
  replication adds ~15 min for PSet parsing).
- Windows has the Arial font required by the submission figures; on Linux substitute
  Liberation Sans.

## Raw data (all public; downloaded 2026-11-16/17)

Files **included** in this repository:

| file | source | size |
|---|---|---|
| `data/raw/mutations_latest.csv.gz` | Sanger Cell Model Passports bulk download | 19 MB |
| `data/raw/model_list_latest.csv.gz` | same | 142 KB |
| `data/raw/Cell_Lines_Details.xlsx` | same | 115 KB |
| `data/raw/growth_rate_20220907.csv` | same (growth-rate dataset) | 41 KB |
| `data/raw/ccle2019_braf_mutations.json` | cBioPortal REST API, study `ccle_broad_2019` | 150 KB |
| `data/raw/ccle2019_samples_all.json` | same, all 1,739 sample IDs | 41 KB |
| `data/raw/ccle2019_samples_sequenced.json` | same, sequenced subset list | 4 KB |
| `data/processed/coread_raw_wells.parquet` | extracted per-well intensities (built by script 04) | 4 MB |

Files **excluded** (too large for the repository; see `data/AVAILABILITY.md` for exact
download URLs and SHA-256 checksums):

| file | source | size |
|---|---|---|
| `GDSC2_fitted_dose_response_24Jul22.csv` | Sanger FTP release directory | 38 MB |
| `GDSC1_fitted_dose_response_24Jul22.csv` | same | 52 MB |
| `GDSC2_public_raw_data_24Jul22.csv` | same | 2.0 GB |
| `CCLE_RNAseq_genes_rpkm_20180929.gct.gz` | Broad CCLE | 137 MB |
| `PSet_CTRPv2_2015.rds` | Orcestra canonical PharmacoSets (Zenodo DOI 10.5281/zenodo.7826870) | 39 MB |
| `PSet_PRISM_2020.rds` | Orcestra canonical PharmacoSets (Zenodo DOI 10.5281/zenodo.7826864) | 82 MB |

Patient data are fetched live through the cBioPortal REST API (studies `coad_silu_2022`,
`coadread_tcga`, `ccle_broad_2019`); no credentials required.

**Not accessible in this environment** (tested, do not retry): CCLE copy number, DepMap
CRISPR gene effect, PRISM secondary screen via DepMap, figshare search API, PharmacoDB API.

---

## Script order

Phase 1 — feasibility census (decides what topic is even possible)

| # | script | output | key number |
|---|---|---|---|
| 13 | `13_gene_power_census.py` | `results/13_gene_power_census.csv` | 278 CRC-only pairs; only 5 pass a naive filter, all implausible (n=3–7) |
| 14 | `14_census_verify_and_pancancer.py` | `results/14_pancancer_census.csv` | line matching verified 46/46; CRC sign concordance 47.8% vs pan-cancer 65.1% |
| 15 | `15_braf_lineage_paradox.py` | `results/15_braf_by_lineage.csv` | BRAF V600 CRC 3.3 µM vs melanoma 0.1 µM dabrafenib |
| 16 | `16_rtk_score_vs_brafi.py` | `results/16_brafi_rtk_score.csv` | RTK score r = +0.09, p = 0.62 → hypothesis rejected |
| 17 | `17_pathway_gene_x_chemo.py` | `results/17_pathway_gene_x_chemo.csv` | 0/14 genes show gene × chemo interaction |

Phase 2 — the paper's statistical backbone

| # | script | output | key number |
|---|---|---|---|
| 18 | `18_replication_filter_and_null.py` | `results/18_summary.json`, `18_surviving_pairs.csv`, `18_null_distribution.csv` | 1,338 testable → 5 replicate, all BRAF; null median 0, empirical p = 0.005 |
| 19 | `19_mutational_burden_confounder.py` | `results/19_burden_adjusted_scan.csv` | burden-adjusted BRAF β = −2.14, p = 9.4e-11 |
| 20 | `20_codon_vs_transcriptome.py` | `results/20_codon_vs_transcriptome.csv` | expression recovers BRAF at AUC 0.906 |
| 21 | `21_incremental_value.py` | `results/21_incremental_value.csv` | ΔR² beyond BRAF = +0.056 |
| 22 | `22_crc_exception_quantified.py` | `results/22_crc_exception.csv` | BRAFi/MEKi dissociation: CRC 5.28 vs 2.66, p = 1e-4; GDSC1 3.59 vs 2.19, p = 0.005; cross-campaign r = 0.73 |
| 23 | `23_patient_layer_braf.py` | `results/23_patient_braf_analysis.csv` | BRAF V600 → CMS1 28/47 vs 10/240; OS HR 0.96 |
| 24 | `24_figures_v2.py` | `results/fig1_evidence_chain.png`, `fig2_translational.png` | internal review figures (superseded by script 33) |
| 27 | `27_msi_positive_control.py` | `results/27_msi_positive_control.csv` | MSI vs MSS 5-FU: 1.71× (GDSC2) vs 0.69× (GDSC1), sign reversal; MDE 2.3× |

Phase 3 — corrections and robustness added before submission (v3→v4 of the manuscript)

| # | script | output | key number |
|---|---|---|---|
| 30 | `30_burden_enrichment_recount.py` | `results/30_burden_enrichment.csv`, `30_burden_adjusted_scan_recomputed.csv`, `30_recount_summary.json` | 26 testable genes (7 of the 33 fail thresholds); 25/26 (96%) enriched in high-burden lines, median d = +1.40; ERCC5 +5.99; BRAF +0.20; NRAS −0.04 |
| 31 | `31_gr_dissociation.py` | `results/31_gr_dissociation.csv`, `31_gr_summary.json` | GR50 identifiable 17/43 V600 lines; Δ_GR 5.6 vs 2.1, p = 0.021; dabrafenib GR50 1.33 vs 0.16 µM; trametinib equipotent |
| 32 | `32_third_party_replication.py` | `results/32_third_party_replication.csv`, `32_summary.json`, per-line tables | CTRPv2: 5/5 MAPK drugs significant pan-cancer; CRC dissociation components replicate (vemurafenib 0.09 vs 0.20, p = 0.003). PRISM: 9/9 RAF/MEK inhibitors significant; dabrafenib+trametinib dissociation index p = 0.009 |
| 32b | `32b_prism_mela_fix.py` | patch outputs of 32 | PRISM melanoma label is `skin`, not `melanoma`; n_mel 0 → 30 |
| 33 | `33_figures_jtm.py` | `results/fig1_evidence_chain_jtm.png/.pdf`, `fig2_translational_jtm.png/.pdf` | submission figures: Arial, ≥7 pt, Okabe–Ito colour-blind-safe palette, 26-gene panel C, seeded jitter, µM units |

Submission preparation

| # | script | output | purpose |
|---|---|---|---|
| 25 | `25_verify_refs.py` | `results/ref_check.csv` | verify each bibliography entry against Europe PMC |
| 26b | `26b_insert_citations.py` | — | insert citation placeholders in the manuscript text |
| 26 | `26_apply_citations.py` | — | number by first appearance (Vancouver), rebuild the bibliography |
| 28 | `28_ref_audit.py` | stdout | numbering continuity, first-appearance order, no uncited entries, DOI coverage |
| 29 | `29_add_dois.py` | — | append DOIs and cross-check titles against Crossref |

Scripts 01–12 belong to an earlier, abandoned framing (response-measurement forensics) and
are retained because their numbers are cited in the manuscript's introduction (censoring
rates, cross-campaign agreement of 5-FU). Their README is `README_v1_forensics.md`.
`_probe_pset.py` / `_probe_gr.py` are read-only structure diagnostics used while building
scripts 31/32 and are kept for transparency.

---

## Every number in the manuscript, and where it comes from

| claim | value | file |
|---|---|---|
| colorectal curves censored above highest concentration | 73.0% | `results/audit_COREAD_by_drug.csv` |
| 5-FU / oxaliplatin / irinotecan censoring | 63.0 / 81.3 / 76.1% | same |
| compounds 100% censored in every colorectal line | 95 of 286 (33.2%) | same |
| median drug censoring rate | 88.9% | same |
| testable gene–drug pairs | 1,338 | `18_summary.json` |
| pairs surviving replication | 5, all BRAF | `18_surviving_pairs.csv` |
| permutation null | median 0, 95th pct 0, max 1 | `18_null_distribution.csv` |
| empirical p | 0.005 | `18_summary.json` |
| FDR hits, common-drug set | 61 | `18_summary.json` |
| FDR hits, all compounds | 118, 83% positive, binomial p = 1.5e-13 | inline check (manuscript Methods) |
| burden-adjusted BRAF effect | β = −2.14, p = 9.4e-11 | `30_burden_adjusted_scan_recomputed.csv` |
| burden enrichment, 26 testable genes | 25/26 (96%), median d = +1.40 | `30_burden_enrichment.csv` |
| largest / smallest / only negative enrichment | ERCC5 +5.99 / BRAF +0.20 / NRAS −0.04 | same |
| burden vs mean sensitivity | r = −0.03, p = 0.37, n = 960 | `30_recount_summary.json` |
| expression → BRAF status | AUC 0.895–0.943, median 0.906 | `20_codon_vs_transcriptome.csv` |
| expression increment beyond BRAF | median ΔR² = +0.056 | `21_incremental_value.csv` |
| BRAFi/MEKi dissociation, CRC vs others | 5.28 vs 2.66, p = 1.0e-4 | `22_crc_exception.csv` |
| same, GDSC1 | 3.59 vs 2.19, p = 0.005 | same |
| dissociation cross-campaign | r = 0.732, sign 100%, n = 33 | same |
| dabrafenib IC50, BRAF V600 CRC vs melanoma | 3.31 vs 0.11 µM (29.6×) | same |
| trametinib IC50, both lineages | 0.01 µM, equipotent | same |
| GR-corrected dissociation | Δ_GR 5.6 vs 2.1, p = 0.021 | `31_gr_dissociation.csv` |
| GR50 identifiable | 17/43 V600 lines (3 CRC, 14 other; 25 censored) | `31_gr_summary.json` |
| dabrafenib / trametinib GR50, CRC vs other | 1.33 vs 0.16 µM / 0.003 vs 0.006 µM | same |
| IC50 pipeline check, raw wells vs release fit | Spearman 0.86 (dabrafenib), 0.95 (trametinib) | same |
| CTRPv2 pan-cancer, 5 MAPK drugs | all p < 1e-6 (V600 n = 26–56) | `32_third_party_replication.csv` |
| CTRPv2 CRC vs melanoma, vemurafenib / PLX4720 / selumetinib | 0.09 vs 0.20 p=0.003 / 0.11 vs 0.20 p=0.015 / 0.33 vs 0.34 p=0.34 | same |
| PRISM pan-cancer, 9 RAF/MEK inhibitors | all p < 1e-4 (V600 n = 36–45) | same |
| PRISM CRC vs melanoma, dabrafenib / trametinib | 0.12 vs 0.43 p=0.010 / 0.40 vs 0.46 p=0.19 | same |
| PRISM dissociation index | p = 0.009 (CRC n=6, melanoma n=30) | same |
| anchor check on published vs recomputed AAC | published fails, recomputed passes | `32_summary.json` |
| BRAF V600 → CMS1 | 28/47 vs 10/240 | `23_patient_layer_braf.py` stdout |
| OS, stage I–III, adjusted | HR 0.96 (0.57–1.61), p = 0.86 | same |
| MSI models among 46 colorectal lines | 15 MSI / 31 MSS | `27_msi_positive_control.py` stdout |
| 5-FU, MSI vs MSS (positive control) | 1.71× (GDSC2, p = 0.088) vs 0.69× (GDSC1, p = 0.14), sign reversal | `27_msi_positive_control.csv` |
| minimum detectable effect at colorectal n | d = 1.34 → 2.3-fold IC50 (80% power) | same |
| BRAF × chemo interaction | HR 1.88 (0.48–7.29), p = 0.36 | same |
| TCGA BRAF V600 | n = 18, HR 1.73 (0.71–4.25), p = 0.23 | same |

---

## Unit correction provenance (v2 → v4 of the manuscript)

The GDSC fitted tables store `LN_IC50` as the natural logarithm of the IC50 **in µM**
(verified against the `MIN_CONC`/`MAX_CONC` columns, which span the µM range, and the GDSC
documentation). An early draft reported absolute potencies in nM; all absolute claims in
the current manuscript are in µM (3.31 vs 0.11 µM dabrafenib; 0.01 µM trametinib; GR50
1.33 vs 0.16 µM). The ratio claims (29.6×, rounded to 30-fold) were never affected.

## Known variance between runs

- Permutation null uses seed 20261116. Re-running with a different seed moves the empirical
  p between 0.005 and 0.010; the qualitative conclusion (observed far outside the null) does
  not change.
- Elastic-net cross-validated R² varies by ±0.08 with sample ordering. This is reported in
  the manuscript as evidence of instability rather than tuned away.
- Figure 1D jitter uses seed 20261116 (script 33); the plot is otherwise deterministic.

## Rejected hypotheses (kept deliberately)

1. Baseline RTK/feedback score predicts BRAF-inhibitor response within BRAF V600 lines → r = +0.09, p = 0.62.
2. MAPK-output score predicts response within lineage → melanoma ρ = −0.12, colorectal ρ = −0.04.
3. Mutational burden explains the directional bias of single-campaign associations → burden unrelated to mean sensitivity (r = −0.03).
4. Fourteen pre-specified pathway genes × adjuvant chemotherapy interaction in AC-ICAM → 0/14 significant.
5. CTRP's own published AAC values reproduce the anchor lines' pharmacology → they fail the anchor check (no separation between A-375 and HT-29/RKO); the PharmacoGx recomputed AAC passes.

These are documented because a paper that reports only what worked is not checkable.

---

## Reference audit (2026-11-16)

All 20 references were verified against Crossref / Europe PMC by DOI, and every entry was
checked for year, volume, page range, first author and journal. Procedure:

1. `25_verify_refs.py` — title search on Europe PMC for the original 18; flagged ref 8 page
   range (`564–76` → `564–576.e16`) and one missing subtitle (REMARK).
2. `26b_insert_citations.py` + `26_apply_citations.py` — placeholders in text, renumbered by
   first appearance (Vancouver). Two references were added because claims needed support:
   NCI-60 (Shoemaker 2006) and CRC BRAF epidemiology (Yaeger 2018).
3. `29_add_dois.py` — DOI appended to each entry and re-resolved through Crossref; title
   similarity ≥ 0.85 for all 20 (two flagged at 0.77/0.79 were a missing subtitle and an
   HTML-wrapped NEJM title, both benign).
4. `28_ref_audit.py` — final check: numbering contiguous 1–20, first-appearance order strictly
   ascending, no uncited entry, no unreferenced number, 20/20 with DOI. All pass.

Figure numbering was corrected at the same time: the two figures are Figure 1 (A–D) and
Figure 2 (A–C). The submission figures are the `*_jtm` versions (script 33).
