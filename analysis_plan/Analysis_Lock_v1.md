# Analysis Lock v1 — Azithromycin treatment-failure prediction model

**Status:** locked before the final run of `R Notebook/V4.10 Modeling Full Notebook.ipynb`
**Date written:** 8 October 2026
**Applies to:** every number, table and figure that will appear in manuscript v7.8 and the Supplementary Appendix

> This file must be committed to git **before** the final run. The commit timestamp is the evidence that the decisions below were fixed before the results were seen. Nothing below may be changed after the run because of the results. If a change becomes unavoidable (e.g. a bug), it is recorded in Section 9 with the date and the reason, and reported as a deviation.

---

## 1. Data and pre-processing (unchanged from V4.00)

| Item | Decision |
|---|---|
| Analysis set | 329 children with a determined 48-hour outcome; 132 events |
| Data file | `data/azithromycin_pediatric_clinical_data.rds` |
| Candidate variables | 23 (list in Supplementary Section S1.1) |
| Respiratory rate | Values outside 15–70 breaths/min (2 values) replaced by the mean of the remaining observations |
| Microscopy bands | Mapped to numeric midpoints as in S1.1 |
| Stool colour | Four binary indicators (green, red/bloody, yellow, brown) |

## 2. Pre-admission antibiotic variable

**Decision:** use the five-category variable defined in the statistical analysis plan (`Documents/Methodology_Standard.md`, `Documents/Details Methodology.md`):

| Category | Source regimens | n |
|---|---|---|
| None (reference) | No antibiotic | 112 |
| Ciprofloxacin | Ciprofloxacin | 81 |
| Metronidazole | Metronidazole | 85 |
| Ciprofloxacin + metronidazole | Ciprofloxacin + Metronidazole | 35 |
| Other | Cefixime; Cefixime + Metronidazole; Cefixime + Amoxicillin; Ciprofloxacin + Amoxicillin / Cefuroxime / Erythromycin; Erythromycin | 16 |

**Reason:** this grouping is pre-specified in the SAP. It is **not** chosen by comparing model performance between the 11-level, 5-level and binary versions. Grouping categories without reference to the outcome follows Harrell (Regression Modeling Strategies, 2nd ed., ch. 4).
**Note:** the near-zero-variance filter inside the recipe may remove the sparse `Other` dummy within a training fold. That behaviour is part of the pre-specified pipeline and is not overridden.

## 3. Feature selection

| Item | Decision |
|---|---|
| Procedure | Unchanged: union of LASSO (α=1) and Elastic Net (α=0.5) at `lambda.1se`, plus top-quartile mutual information; then VIF > 5 removal among numeric predictors |
| Where it runs | Inside every outer training fold of every repeat (no information from the held-out fold) |
| **Final predictor set** | Variables selected in **≥ 50 %** of all outer-fold selections (10 repeats × 5 folds = 50 selections, i.e. ≥ 25/50) |
| Not used | Selection on the complete dataset (V4.00 approach); hand-picked "old 11" set |

**Reason:** the final model must come from the same procedure that was validated, and selection frequency over resamples is a recognised stability criterion (Heinze et al., Biom J 2018; Meinshausen & Bühlmann, JRSS-B 2010). The 50 % cut-off is a majority-vote convention chosen here and fixed now.

## 4. Algorithms and hyperparameter grids

Eight algorithms, unchanged. Tuning uses **fixed grids** (no randomly generated candidates):

| Model | Grid | Size |
|---|---|---|
| Elastic Net | penalty ∈ 10^{−4, −3, −2, −1.5, −1, −0.5}; mixture ∈ {0.2, 0.5, 0.8, 1} | 24 |
| Random Forest (500 trees) | mtry ∈ {2, 3, 4, 6}; min_n ∈ {2, 5, 10, 20} | 16 |
| XGBoost | trees ∈ {100, 300}; tree_depth ∈ {2, 4, 6}; learn_rate ∈ {0.01, 0.05, 0.1} | 18 |
| KNN | neighbors ∈ {5, 10, 15, 20, 30, 40}; weight_func ∈ {rectangular, triangular, optimal} | 18 |
| SVM (RBF) | cost ∈ 2^{−2, 0, 2, 4}; rbf_sigma ∈ 10^{−3, −2, −1} | 12 |
| Decision tree | cost_complexity ∈ 10^{−4, −3, −2}; tree_depth ∈ {3, 5, 8}; min_n ∈ {10, 20} | 18 |
| Naive Bayes | smoothness ∈ {0.5, 1, 1.5}; Laplace ∈ {0, 1} | 6 |
| Neural network | hidden_units ∈ {2, 5, 10}; penalty ∈ 10^{−3, −2, −1}; epochs ∈ {100, 300} | 18 |

Inner-loop tuning metric: AUROC.

## 5. Resampling design

| Item | Decision |
|---|---|
| Outer loop | 5-fold, stratified by outcome |
| Inner loop | 10-fold, stratified, within each outer training fold |
| Repeats | **10** complete repetitions of the nested procedure |
| Repeat seeds | `SEED + 100 × r`, r = 1…10, with `SEED = 123` (i.e. 223, 323, …, 1123) |
| Outer folds | Identical for all eight algorithms within a repeat, so models are compared on the same splits |

**Reason:** a single nested CV run gave AUROC between 0.69 and 0.76 depending only on the seed (`GAPS/AUROC_Variability_Diagnosis.md`). Repeated nested CV is required for stable selection and assessment (Krstajic et al., J Cheminform 2014; Vabalas et al., PLoS One 2019).

## 6. Performance estimation

| Measure | How it is computed and reported |
|---|---|
| AUROC | Pooled out-of-fold AUROC within each repeat → **mean across 10 repeats**. 95 % CI from a patient-level bootstrap (B = 2000) of the repeat-mean AUROC. Repeat SD and min–max also reported |
| Calibration | Intercept (offset model) and slope per repeat → mean; Brier per repeat → mean; Hosmer–Lemeshow p per repeat → median |
| Threshold | Youden's J per repeat → **median** across repeats = the model's fixed threshold |
| Sens, spec, PPV, NPV, accuracy, F1 | Per repeat at the fixed threshold → mean (min–max across repeats) |
| Net benefit / DCA | On all repeats stacked (equal to the mean net benefit across repeats) |
| Figures (ROC, calibration, DCA) | All repeats stacked |
| Subgroups | Final model; AUROC per repeat within subgroup → mean, with patient-level bootstrap CI |

## 7. Final-model selection rule

Unchanged from Supplementary Section S1.5, applied to the repeated-CV results:

1. Rank the eight models by mean AUROC across repeats → the leader.
2. Paired patient-level bootstrap (B = 2000; the same resamples for every model) of the difference in repeat-mean AUROC between the leader and each other model.
3. Candidate set = the leader plus every model with p ≥ 0.05 against it.
4. Within the candidate set, keep models with adequate calibration (median Hosmer–Lemeshow p > 0.05); if none qualifies, keep all candidates.
5. Select the model whose mean calibration slope is closest to 1; ties broken by higher mean net benefit over threshold probabilities 30–40 %.

The notebook applies this rule automatically. The model it selects is the final model, **whatever it is**. Figures show the final model plus the next two models by mean AUROC.

## 8. Final model, interpretation and reproducibility

| Item | Decision |
|---|---|
| Final fit | Selected algorithm refitted on all 329 children using the Section 3 predictor set, hyperparameters tuned by stratified 5-fold CV on the same fixed grid |
| SHAP | kernelSHAP on the final fit; all 329 children as explained and background data |
| Parallel workers | Fixed at 4 (`plan(multisession, workers = 4)`) |
| Threads | ranger `num.threads = 1`; xgboost `nthread = 1` |
| Checkpointing | Each repeat saved to `results/v4.10/repeats/`; a rerun skips completed repeats |
| Environment | `sessionInfo()` written to `results/v4.10/sessionInfo.txt` |
| Reproducibility check | Before the full run, repeat 1 is run twice; the two AUROC tables must match exactly |

## 9. What will not be done after the run

- No change to variables, coding, grids, algorithms, fold numbers, seeds, threshold rule or selection rule because of the results.
- No reruns with other seeds to obtain a different AUROC.
- Any deviation forced by a bug or error is logged here:

| Date | Deviation | Reason | Effect on results |
|---|---|---|---|
| — | — | — | — |

## 10. Disclosure

The rules in Sections 2, 3, 5 and 7 were fixed after preliminary single-run analyses (May and September 2026) had been inspected. This will be stated in PROBAST+AI Domain 4 and in the Limitations, together with the fact that the final run itself is fully specified by this document.
