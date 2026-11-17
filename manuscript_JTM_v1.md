# An error budget for colorectal cancer pharmacogenomics: censoring, batch variation, and the limits of machine-learned drug resistance

**Running title:** Error budget of CRC drug-response labels

**Article type:** Research (Methodology / comparative evaluation)

**Target collection:** *Harnessing Innovative Machine Learning Techniques to Combat Drug Resistance in Solid Tumors*, Journal of Translational Medicine

**Draft version:** v1, 2026-11-16

---

## Abstract

**Background.** Machine-learning models that predict chemoresistance in colorectal cancer (CRC) are almost universally trained on cell-line dose–response data, with the half-maximal inhibitory concentration (IC₅₀) treated as a continuous, error-free regression target. This assumption has not been audited. We asked how much of the reported predictive performance is attributable to resistance biology, and how much to the measurement properties of the label itself.

**Methods.** We re-analysed the complete Genomics of Drug Sensitivity in Cancer (GDSC) resource: 242,036 fitted dose–response curves across 969 cell lines and 286 compounds, of which 12,538 curves from 46 CRC (`COREAD`) lines were analysed. We (i) quantified censoring of IC₅₀ against tested concentration ranges; (ii) measured assay precision directly from 1.33 million raw plate wells; (iii) estimated label reliability by comparing two independent screening campaigns (GDSC1 vs GDSC2) and by replicate screens within GDSC2; (iv) trained best-practice elastic-net models on CCLE transcriptomes with nested resampling, permutation nulls, and proliferation-adjusted baselines; (v) decomposed predictive performance into between-group separation and within-group ranking; and (vi) transferred the resulting signature to a prospectively annotated patient cohort (AC-ICAM, n = 348) and tested for a score × adjuvant-chemotherapy interaction.

**Results.** Censoring is the rule, not the exception: 73.0% of CRC IC₅₀ values lie above the highest concentration tested and are therefore extrapolations rather than measurements (5-fluorouracil [5-FU] 63.0%, oxaliplatin 81.3%, irinotecan 76.1%), and across the full 286-compound panel the median censoring fraction is 88.9%, with 95 compounds (33.2%) censored in every cell line. Six clinically central agents — cetuximab, panitumumab, capecitabine, trifluridine, regorafenib and encorafenib — are absent from the resource entirely. Assay precision itself is good (median replicate-well CV 9.4%), but at any single concentration the spread between cell lines is only 1.5-fold the replicate error (median signal-to-noise ratio 1.54; 74.4% of drug–concentration strata < 2). Comparing two independent screening campaigns on the same 4,307 cell line–drug pairs, the median pair disagreed 3-fold, and only 18.0% of pairs had a genuinely measurable IC₅₀ in both campaigns. Elastic-net models reached cross-validated R² = 0.29–0.36 for 5-FU and beat permutation nulls for the three fluoropyrimidine/platinum/topoisomerase agents, but failed for all five targeted agents. Decomposition showed that the 5-FU signal is a censoring boundary (AUC 0.78 for discriminating measurable from unmeasurable curves) with no reproducible within-group potency ranking (Spearman +0.33 vs −0.36 across specifications). Substituting the bounded area-under-curve metric for IC₅₀ did not help (3/8 vs 1/8 agents beating null). In patients, the derived score showed no interaction with adjuvant chemotherapy (HR 1.08, 95% CI 0.66–1.76, p = 0.76; n = 232, 90 deaths), matching the null observed for a random gene-set control.

**Conclusions.** The dominant error terms in CRC pharmacogenomics are not model capacity or pipetting noise but (i) censoring of the response label and (ii) batch-to-batch variation between screening campaigns. Under these terms the achievable ceiling for cell-line-trained resistance models is a coarse sensitive/resistant classification that does not transfer to patients. We provide a computable pre-flight audit and recommend that censoring fraction and cross-campaign reliability be reported alongside any predictive metric.

**Keywords:** colorectal cancer; drug resistance; machine learning; GDSC; censoring; reproducibility; 5-fluorouracil; pharmacogenomics

---

## Background

Drug resistance is the principal reason that patients with solid tumours die of their disease, and machine learning (ML) is now routinely proposed as the route to predict it. In colorectal cancer, the dominant computational paradigm is unchanged since the first large pharmacogenomic screens: transcriptional features from cell lines are regressed onto a drug-sensitivity metric, most often IC₅₀, and the resulting model is reported with a cross-validated R² or an area under the receiver-operating characteristic curve [1–4].

That paradigm rests on an assumption that is rarely stated and, to our knowledge, never audited end-to-end: that IC₅₀ is a measured continuous quantity whose residual variation is biological. Three properties of the underlying data challenge this. First, IC₅₀ is defined as the concentration producing 50% inhibition, yet for most compounds most cell lines are not inhibited by 50% anywhere in the tested range; the value that is reported is then a property of the curve-fitting algorithm, not of the cell. Second, dose–response assays are run in batches, and batch effects in high-throughput pharmacology are known to be large [5,6]. Third, and most consequentially for a translational read-out, a cell monolayer on plastic lacks the stromal, immune and pharmacokinetic context in which clinical resistance actually arises.

Previous work has documented inconsistencies between the two flagship pharmacogenomic resources [5,6], and benchmarking studies have compared prediction algorithms [7,8]. What has not been done is a quantitative decomposition that assigns the observed (and unobserved) predictive performance to specific, measurable error terms — in other words, an error budget. Without one, the field cannot tell whether a reported R² of 0.4 represents a biological discovery, an artefact of the assay's detection limit, or both.

Here we construct such an error budget for colorectal cancer, the tumour type with the largest chemotherapeutic evidence base. We work exclusively from primary material: raw plate-level fluorescence, unfitted dose–response data, and the two independent screening campaigns that produced them. We then ask what the budget implies for ML, and we test the implication in an independent, prospectively annotated patient cohort.

---

## Methods

### Data sources

All analyses used public primary data, downloaded on 2026-11-16.

**GDSC2 and GDSC1 fitted dose–response data** (`GDSC2_fitted_dose_response_24Jul22.csv`, `GDSC1_fitted_dose_response_24Jul22.csv`; Wellcome Sanger Institute, `ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/`). GDSC2 comprised 242,036 curves (969 cell lines, 286 compounds); GDSC1 comprised an independent campaign with a partially overlapping compound panel. Colorectal lines were defined by the resource's own annotation `TCGA_DESC == "COREAD"` (n = 46 lines; note that the resource does not separate colon from rectum).

**GDSC2 raw plate data** (`GDSC2_public_raw_data_24Jul22.csv`, 2.0 GB, 15,069,449 wells). We extracted the 1,886,619 wells belonging to the 46 CRC lines (2,032 plates) and, after excluding unused positions, analysed 1,325,787 wells. Each well carries a raw fluorescence intensity, a compound identifier, a nominal concentration, and a plate barcode.

**CCLE transcriptomes** (`CCLE_RNAseq_genes_rpkm_20180929.gct.gz`; Broad Institute). Of 46 CRC lines, 39 matched to CCLE by normalised name (the seven unmatched lines are listed in Additional file 1).

**Patient cohorts** (cBioPortal REST API v2). AC-ICAM / `coad_silu_2022` (n = 348; primary colorectal tumours with annotated adjuvant treatment, overall survival and distant recurrence) served as the patient cohort. TCGA-COADREAD (`coadread_tcga`) was examined but could not be used: the adjuvant pharmacotherapy field was populated in only 36 of 636 patients (20 YES, 16 NO).

### Censoring analysis

For every curve we compared the fitted IC₅₀ (exp of the reported natural-log value) against the lowest and highest concentration tested. A curve was defined as *identifiable* when the IC₅₀ fell inside the tested interval, and as *extrapolated* when it fell outside. Because the published files contain a point estimate only — no confidence interval and no curve-quality flag — we report the censoring fraction as the primary descriptor of label validity.

### Assay precision from raw wells

For each plate we computed the median DMSO-control intensity and normalised every treated well to it, yielding a relative viability. We then computed the coefficient of variation (CV) across replicate wells sharing plate, compound, concentration and cell line (51,372 groups; median 3 wells per group). Separately, for each (compound, concentration) stratum with ≥20 cell lines, we computed the between-cell-line standard deviation of mean viability. The ratio of between-cell-line SD to replicate CV is reported as the single-concentration signal-to-noise ratio (SNR).

### Label reliability

Two orthogonal estimates were obtained.

*Within-campaign replicates.* Within GDSC2, 379 (cell line, compound) pairs were assayed in more than one screening instance (758 curves). We computed the pooled within-pair standard deviation of ln(IC₅₀), the typical fold-difference between replicates (exp[SD·√2]), and the intraclass correlation ICC(1).

*Between-campaign replicates.* For GDSC1 and GDSC2, we intersected (cell line, compound) pairs, averaged within campaign, and computed Pearson correlation, ICC, median absolute fold-difference, and the fraction of pairs agreeing within 2-fold. This comparison is the natural estimate of test–retest reliability, because neither campaign was designed for this analysis. We additionally stratified by joint identifiability, i.e. whether the IC₅₀ was measurable in *both* campaigns.

Throughout, we use the classical result that if X₁ = T + e₁ and X₂ = T + e₂, then corr(X₁,X₂) = var(T)/var(X) = ρ, and a model that perfectly recovered T would achieve corr(M,X₁) = √ρ, so that R² ≤ ρ. Test–retest correlation therefore constitutes a hard ceiling on any cross-validated R².

### Predictive modelling

Transcriptomes were log₂(RPKM+1) transformed and reduced to the 2,000 most variable genes (unsupervised, label-blind). Targets were ln(IC₅₀) per cell line, or curve AUC where stated. The primary model was elastic net with inner cross-validation for α and l1-ratio, embedded in a scikit-learn pipeline so that scaling and penalisation were fitted within folds. Outer evaluation used 5-fold cross-validation with `cross_val_predict`; the same procedure was repeated on 10–20 label permutations to construct a null distribution. A proliferation baseline (mean expression of 14 proliferation genes) and, in the patient layer, a size-matched random gene set served as negative controls. All models used a fixed random seed; we additionally report the sensitivity of R² to sample ordering.

### Decomposition of predictive performance

For the agent with the strongest apparent signal (5-FU), we partitioned performance into three components: (i) overall R²; (ii) discrimination of identifiable from extrapolated curves by the cross-validated prediction (AUC); (iii) rank correlation within the identifiable subset and within the extrapolated subset. Component (iii) was computed both from a model trained on all lines and applied within strata, and from a model retrained on the identifiable subset alone, to expose specification dependence.

### Patient transfer and interaction testing

The non-zero elastic-net coefficients from the 5-FU model defined a 37-gene signature (35 mapped to Entrez identifiers). Patient scores were computed as Σβᵢzᵢ / √Σβᵢ², with zᵢ standardised within cohort, which is platform-robust and requires no cross-platform recalibration.

The primary estimand was the interaction between standardised score and adjuvant chemotherapy in a Cox proportional-hazards model, adjusted for age, sex and stage, restricted to stage II–III disease (the adjuvant indication). A score that captures resistance — as opposed to prognosis — should show a positive interaction: worse outcome under chemotherapy in high-score patients, no difference in untreated patients. Random-gene-set and proliferation scores were analysed identically.

### Reproducibility

All steps are scripted (Python 3.13; pandas, scikit-learn, statsmodels, lifelines-free Cox via `statsmodels.duration.hazard_regression.PHReg`). Scripts are numbered in execution order and write every intermediate table to disk.

---

## Results

### The colorectal cell-line panel is small and therapeutically incomplete

GDSC2 contains 242,036 curves, but only 12,538 (5.2%) derive from colorectal lines, and these come from just **46 cell lines**. Six agents central to CRC practice are absent from the resource altogether: cetuximab and panitumumab (whose activity in vitro depends on antibody-dependent cellular cytotoxicity and is not measurable as a proliferation IC₅₀), capecitabine (a prodrug requiring hepatic activation), trifluridine, regorafenib and encorafenib. The clinically best-established resistance rule in CRC — RAS/BRAF mutation conferring anti-EGFR resistance — therefore cannot be studied in this resource at all.

### Most IC₅₀ values are extrapolations, not measurements

Across all 12,538 CRC curves, 73.0% of fitted IC₅₀ values exceeded the highest concentration tested, and 0.1% fell below the lowest; only **26.9% were identifiable** (Fig. 1A,B). Censoring was drug-specific and, critically, worst for the agents that define CRC chemotherapy: oxaliplatin 81.3%, irinotecan 76.1%, 5-FU 63.0%. Fulvestrant, an oestrogen-receptor antagonist with no CRC indication included as a negative control, was 100% censored — no CRC line responded at any tested concentration, confirming that the metric behaves as expected in the absence of biology.

At the level of the compound panel the picture is starker still. Across all 286 compounds the median censoring fraction was **88.9%** of cell lines, and **95 compounds (33.2%) were censored in every one of the 46 lines** — that is, they produced no measurable inhibition anywhere in the CRC panel, yet contribute a row to every pan-drug analysis. Only 69 compounds (24.1%) had an identifiable IC₅₀ in at least half the panel.

For 5-FU the geometry is stark. The maximum tested concentration was 32 µM. The 17 identifiable lines had IC₅₀ between 3.4 and 27.1 µM; the 29 censored lines were assigned extrapolated values from 35.4 µM to 13,924 µM (Fig. 1A). The 400-fold spread among those 29 values is generated entirely by the curve-fitting procedure.

### Assay precision is good; single-concentration biological signal is not

Direct measurement from raw wells contradicts the intuition that the problem is technical noise. Across 51,372 replicate groups the median well-to-well CV was 9.4% (90th percentile 23.9%), which is acceptable for a cell-viability assay. However, at any given compound and concentration the standard deviation between cell lines was only 1.5 times the replicate CV (median SNR 1.54), and 74.4% of the 976 drug–concentration strata had SNR < 2 (Fig. 2). In other words, at a single dose, ranking two cell lines is barely distinguishable from ranking two replicate wells of the same line. Curve fitting across 9–13 concentrations recovers precision by averaging, which is why IC₅₀ can be estimated at all — but only where the curve actually inflects within the tested window, i.e. in the 27% of cases that are identifiable.

### Two independent campaigns disagree by 3-fold

Comparing GDSC1 and GDSC2 across 4,307 shared (cell line, compound) pairs (46 lines, 121 compounds), overall agreement was r = 0.811, but the median pair differed by 3.0-fold and 67.0% of pairs differed by more than 2-fold (Fig. 1C). For the CRC chemotherapy backbone the disagreement was larger: 5-FU r = 0.667 with a median 3.5-fold difference and 76.7% of pairs beyond 2-fold; SN-38 r = 0.532, median 5.2-fold.

The most informative stratification is by joint identifiability. Only **18.0%** of all pairs, and 35.6% of chemotherapy pairs, had a genuinely measurable IC₅₀ in both campaigns. Within that jointly identifiable stratum, 5-FU agreement fell to r = 0.048 (n = 14). The apparent concordance of 5-FU IC₅₀ across campaigns is therefore carried largely by both campaigns extrapolating in the same direction, not by measuring the same quantity.

Within GDSC2 alone, replicate screening instances of the same (line, drug) pair differed by a median 6.0-fold and 43% had a within-pair SD > 1 ln unit (≈4-fold). Because these replicates derive from different compound batches and screening instances, this estimate captures total batch-plus-technical reproducibility — precisely the error a downstream model inherits.

### Machine learning learns the censoring boundary, not potency

Elastic-net models on CCLE transcriptomes beat permutation nulls for three of eight agents: 5-FU (R² = 0.29–0.36 across fold assignments), oxaliplatin (0.11–0.20) and irinotecan (0.09–0.10). All five targeted agents (SN-38, trametinib, selumetinib, dabrafenib, fulvestrant) failed to beat their nulls (Fig. 3, Table 3). A proliferation baseline was uninformative throughout (R² < 0 for every agent), excluding the simplest confounder.

Notably, reliability and predictability were dissociated. Trametinib had the most reproducible label (r = 0.751) yet was not predictable from baseline transcriptome (R² = −0.06 to −0.19). Label quality alone is therefore not a sufficient prescription.

Decomposing the 5-FU model (Fig. 1D) resolved the apparent success: the cross-validated prediction discriminated identifiable from extrapolated curves with **AUC 0.778**, but showed no reproducible ranking within the identifiable subset (Spearman +0.33 when trained on all lines; −0.36 when retrained on the subset) or within the extrapolated subset (Spearman +0.22, R² = −0.23). The overall R² of ~0.3 is thus a between-group contrast, not a potency model.

We tested the most obvious remedy — replacing the unbounded, frequently censored IC₅₀ with the bounded area under the dose–response curve, which is defined for every curve. This made matters worse: cross-campaign reliability fell from r = 0.811 to 0.664, and the number of agents beating their permutation null fell from 3/8 to 1/8.

### The derived signature does not modify chemotherapy benefit in patients

In AC-ICAM (n = 348; 156 deaths; 104 adjuvant-treated), restricted to stage II–III (n = 232; 90 deaths; 71 treated), the 37-gene 5-FU score showed no interaction with adjuvant chemotherapy (HR 1.08, 95% CI 0.66–1.76, p = 0.76), no main effect on survival (HR 1.15, p = 0.22), and behaved indistinguishably from a size-matched random gene set (interaction HR 0.92, p = 0.73) (Fig. 4A). Chemotherapy was associated with lower mortality in every score stratum (low 43.6%→26.1%, mid 36.8%→30.0%, high 51.0%→28.6%; Fig. 4B) — that is, no subgroup defined by the model failed to benefit.

---

## Discussion

We set out to determine how much of machine-learned chemoresistance in colorectal cancer is biology. The answer, in the data currently available, is: very little of it, and the reason is measurable rather than mysterious.

Three findings organise the result. First, the regression target is mostly not a measurement. In 73% of cases the reported IC₅₀ is an extrapolation beyond the tested concentration range, and for 5-FU the extrapolated values span a 400-fold interval that is generated by the fitting function. Any model trained on this target is trained partly on the behaviour of a curve fitter.

Second, the dominant error term is batch variation, not assay noise. Replicate-well precision is 9.4%, which is good. Yet two independent screening campaigns of the same cell line and drug disagree by 3-fold in the median case, and replicate screening instances within one campaign by 6-fold. This has a specific methodological consequence: the problem cannot be solved by more replicates, more genes, or better architectures, because it is situated upstream of the model, in the irreproducibility of the label across campaigns. It is also why single-cohort cross-validation is structurally over-optimistic — folds are drawn from the same campaign, so they share the very batch structure that does not transfer.

Third, when the error budget is accounted for, the residual signal has a specific and limited shape: a sensitive-versus-resistant contrast (AUC 0.78 at the cell-line level) with no reliable potency ranking inside either group, and no transfer to patients. The failure of the bounded AUC metric to improve matters indicates that the ceiling is not an artefact of how IC₅₀ is parameterised; it survives a change of target.

That the patient-level interaction is null deserves emphasis because it is often assumed to follow from a good cell-line model. It does not. A score can be prognostic without being predictive, and our negative control makes the distinction concrete: the random gene set produced an interaction estimate indistinguishable from the model's.

**Implications.** We see three. (i) *Reporting.* Any pharmacogenomic prediction study should report the censoring fraction of its target and the cross-campaign reliability of that target for the drugs studied. Both are computable in minutes from public files and both bound what can be claimed. (ii) *Design.* Because batch variation dominates, claims of generalisation require validation in an independent screening campaign or an independent technology, not merely held-out cell lines. (iii) *Target selection.* In CRC, cell-line screens cannot address the pathway that matters most clinically (anti-EGFR resistance), because the relevant antibodies have no cell-autonomous read-out. Work in this space should either move to systems that preserve immune effectors or choose endpoints that the assay can actually measure.

**What would break the ceiling.** Our results do not imply that CRC resistance is unpredictable; they imply that the currently distributed labels cannot support prediction beyond coarse classification. Three routes seem plausible. Patient-derived organoids with paired clinical response [9,10] retain some microenvironmental context and would allow the same error budget to be computed one level closer to the clinic. Functional-genomic screens (e.g. CRISPR dependency in the same lines) provide a second, mechanistically distinct axis that is less censored. And adaptive/acquired-resistance designs — measuring response after prolonged exposure rather than at baseline — target the biology that clinical resistance actually reflects, at the cost of throughput.

**Limitations.** The colorectal panel is 46 lines, and 39 matched to CCLE; with 2,000 features this is a low-power setting, and our R² estimates vary by ±0.08 with sample ordering, which we report rather than conceal. GDSC1 and GDSC2 differ in concentration design (5-FU lower bound 0.125 vs 0.0195 µM), so the cross-campaign discrepancy conflates design and biology — although from a downstream user's perspective both are inherited error. Only one patient cohort had both transcriptome and annotated adjuvant chemotherapy; TCGA-COADREAD could not be used because its adjuvant field is populated in 6% of patients. No wet-lab validation was performed; we instead state the experiments that would falsify our conclusions (Additional file 2). Finally, we analysed a single tumour type by design, so the magnitudes reported here should not be generalised to other lineages without re-running the audit.

---

## Conclusions

In colorectal cancer, the labels on which machine-learned drug-resistance models are trained are dominated by two error terms that no modelling choice can remove: censoring of the dose–response read-out, and batch-to-batch irreproducibility between screening campaigns. Under the resulting ceiling, cell-line-trained models can distinguish sensitive from resistant lines but cannot rank potency, and the derived signatures do not modify chemotherapy benefit in patients. Auditing these two terms before modelling would redirect substantial effort toward data generation rather than architecture search.

---

## Declarations

**Ethics approval and consent to participate.** Not applicable; all data are publicly available, de-identified secondary data.

**Consent for publication.** Not applicable.

**Availability of data and materials.** All primary data are public: GDSC (`ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/`), CCLE (`data.broadinstitute.org/ccle/`), AC-ICAM and TCGA-COADREAD (cBioPortal). Analysis scripts (numbered 01–12) and all derived tables are provided; the extracted CRC raw-well table (1,886,619 wells) is available on request or via the accompanying repository.

**Competing interests.** The authors declare no competing interests. *(to be confirmed by all authors)*

**Funding.** *(to be completed)*

**Authors' contributions.** *(to be completed)*

**Acknowledgements.** The authors thank the GDSC, CCLE and cBioPortal teams for open data distribution.

---

## Figure legends

**Fig. 1** Labelling forensics of colorectal drug-response data. **A** 5-FU IC₅₀ for 46 CRC lines against the maximum tested concentration (32 µM); measured values in blue, extrapolated in red. **B** Fraction of lines with IC₅₀ beyond the tested range, by agent. **C** Agreement between two independent screening campaigns (GDSC1 vs GDSC2) for the same cell line and drug. **D** Decomposition of the 5-FU model into overall R², discrimination of the censoring boundary, and within-stratum ranking; crosses mark an alternative specification.

**Fig. 2** Error budget. Fold-uncertainty contributed by each term, on a log scale. Assay precision (well-to-well) is the smallest term; inter-campaign disagreement and censoring dominate.

**Fig. 3** Reliability and predictability are dissociated. Cross-campaign label reliability (x) against cross-validated model R² (y). Blue: beats permutation null; red: fails.

**Fig. 4** Patient-layer test. **A** Hazard ratio for the score × adjuvant chemotherapy interaction in stage II–III disease, with random-gene-set and proliferation controls. **B** Death rate by score tertile and treatment status.

---

## References

1. Barretina J, et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. *Nature*. 2012;483:603–607.
2. Garnett MJ, et al. Systematic identification of genomic markers of drug sensitivity in cancer cells. *Nature*. 2012;483:570–575.
3. Yang W, et al. Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells. *Nucleic Acids Res*. 2013;41:D955–D961.
4. Iorio F, et al. A landscape of pharmacogenomic interactions in cancer. *Cell*. 2016;166:740–754.
5. Haibe-Kains B, et al. Inconsistency in large pharmacogenomic studies. *Nature*. 2013;504:389–393.
6. Safikhani Z, et al. Revisiting inconsistency in large pharmacogenomic studies. *F1000Res*. 2016;5:2333.
7. Costello JC, et al. A community effort to assess and improve drug sensitivity prediction algorithms. *Nat Biotechnol*. 2014;32:1202–1212.
8. Papillon-Cavanagh S, et al. Comparison and validation of genomic predictors for anticancer drug sensitivity. *JCO Clin Cancer Inform*. 2019;3:1–12.
9. Vlachogiannis G, et al. Patient-derived organoids model treatment response of metastatic gastrointestinal cancers. *Science*. 2018;359:920–926.
10. Ooft SN, et al. Patient-derived organoids can predict response to chemotherapy in metastatic colorectal cancer patients. *Sci Transl Med*. 2019;11:eaay2574.
11. Guinney J, et al. The consensus molecular subtypes of colorectal cancer. *Nat Med*. 2015;21:1350–1356.
12. Ribic CM, et al. Tumor microsatellite-instability status as a predictor of benefit from fluorouracil-based adjuvant chemotherapy for colon cancer. *N Engl J Med*. 2003;349:247–257.
13. Corsello SM, et al. Discovering the anti-cancer potential of non-oncology drugs by systematic viability profiling. *Nat Cancer*. 2020;1:235–248.
14. Tsherniak A, et al. Defining a cancer dependency map. *Cell*. 2017;170:564–576.
15. Collins GS, et al. TRIPOD+AI statement: updated guidance for reporting clinical prediction models that use regression or machine learning methods. *BMJ*. 2024;385:e078378.
16. Wolff RF, et al. PROBAST: a tool to assess the risk of bias and applicability of prediction model studies. *Ann Intern Med*. 2019;170:51–58.

*Reference list to be verified against publisher records before submission; the AC-ICAM cohort citation should be added from the cBioPortal study record.*

---

## Additional files

**Additional file 1.** Cell lines in the GDSC CRC panel not matched to CCLE transcriptomes (COLO-205, COLO-320-HSR, CaR-1, DiFi, GP5d, HCC2998, SNU-C2B).

**Additional file 2.** Falsification conditions. (i) If, among CRC lines whose 5-FU IC₅₀ is measurable within the tested range, a transcriptional signature predicts potency with Spearman > 0.5 in the same direction in both GDSC campaigns, the central claim fails. (ii) If knockdown of DPYD reduces 5-FU IC₅₀ by ≥2-fold and that effect is predictable from the features we identified as noise, the signal has been underestimated. (iii) If organoid-based prediction of clinical response materially exceeds the cell-line ceiling, the attribution should shift from assay to model system.
