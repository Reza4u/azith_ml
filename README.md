# azith_ml: predicting azithromycin treatment failure in children with suspected shigellosis

Analysis code and the trained model for:

> Reza N, Islam MR, Al Fidah MF, Amin R, Chisti MJ, Ahmed T, Nuzhat S. *Development of a multivariable machine-learning model to predict azithromycin treatment failure in children with suspected shigellosis: a prognostic modelling study.* BMJ Health & Care Informatics, 2026 (manuscript bmjhci-2026-102619).

The study is a secondary analysis of a phase IIb trial at icddr,b Dhaka Hospital, Bangladesh (ClinicalTrials.gov NCT05121974). The analysis included 329 children and 132 treatment failures. Eight algorithms were compared in ten repeats of nested 5×10-fold cross-validation. A Random Forest with ten predictors was selected by a pre-specified rule: AUROC 0.71 (95% CI 0.66 to 0.76), calibration slope 1.05.

**Participant data are not in this repository** (see [Data](#data)).

## Repository layout

| Path | Contents | Paper |
|---|---|---|
| `run_all.R` | Runs the full analysis in order | All results |
| `R/00_setup.R` | Packages, fixed analysis settings (seeds, folds, repeats, workers) and labels | Methods; Analysis plan |
| `R/01_data.R` | Loads the data and applies the pre-specified coding | Methods (Predictors); S1.1 |
| `R/02_descriptive_tables.R` | Baseline table, data dictionary, antibiotic regimens, outcome components | Table 1; Tables S1, S5a–b, S6a–b |
| `R/03_feature_selection.R` | LASSO + Elastic Net + mutual information, then VIF filter | S1.2 |
| `R/04_models.R` | Eight algorithms and their fixed tuning grids | S1.3 |
| `R/05_nested_cv.R` | Repeated nested cross-validation, reproducibility check, final predictor set | S1.4; Table S2 |
| `R/06_performance.R` | Performance across repeats, bootstrap CIs, pre-specified final-model rule | Table 2; Tables S4a–d, S7, S8 |
| `R/07_decision_curve.R` | Net benefit | Table S9 |
| `R/08_final_model_shap.R` | Final model refit, kernelSHAP and exact Shapley values | S1.8; Table S10 |
| `R/09_subgroups.R` | Subgroup performance of the final model | Table S3 |
| `R/10_figures.R` | ROC, decision curve and calibration plots | Figure 2; Figure S2 |
| `R/11_sample_size.R` | Events per predictor; Riley et al. minimum sample size | Methods; S6.4 |
| `R/12_key_numbers.R` | Numbers quoted in the text | Abstract; Results |
| `python/figure3_beeswarm.py` | SHAP beeswarm plot from the exact Shapley values | Figure 3 |
| `tests/check_against_paper.R` | Compares every reported value with the code's output | — |
| `figures_src/` | Mermaid sources of the two diagrams | Figure 1; Figure S1 |
| `model/final_model_fit_light.rds` | Trained final model with all training data removed | S1.8; S7 |
| `model/predict_example.R` | Predicted probability for a new child | S7 |
| `model/make_light_model.R` | How the light model was made from the full fit | — |
| `analysis_plan/Analysis_Lock_v1.md` | Analysis plan fixed before the final run | S9 |
| `R/package_versions.csv` | Versions of the key packages used for the paper | S1.7 |
| `results/` | Aggregate result tables from the run reported in the paper, and its session information | Tables 1–2; S1–S10; S1.7 |
| `figures/` | Figures 1–3 and S1–S2 as published | Figures |

## Reproducing the analysis

Requirements: R 4.5.3 and the package versions in `R/package_versions.csv` (tidymodels 1.4.1, tune 2.0.1, recipes 1.3.1, glmnet 4.1-10, ranger 0.18.0, xgboost 3.2.1.1, kernlab 0.9-33, kknn 1.4.1, rpart 4.1.24, nnet 7.3-20, klaR 1.7-4, pROC 1.19.0.1, kernelshap 0.9.1, gtsummary 2.5.0, pmsampsize). Figure 3 needs Python 3 with numpy, pandas and matplotlib. `00_setup.R` warns if a key package version differs.

With the two data files in `data/`, run from the repository root:

```sh
Rscript run_all.R                      # all tables, Figure 2, Figure S2 (about 45 min)
python3 python/figure3_beeswarm.py     # Figure 3
Rscript tests/check_against_paper.R    # compare with the published numbers
```

Everything is written to `output/` (excluded from git):

- `output/tables/` holds aggregate tables named by their number in the paper (`Table2_model_performance.csv`, `TableS4a_auroc_all_models.csv`, …) and `key_numbers.csv`.
- `output/figures/` holds Figure 2, Figure 3 and Figure S2.
- `output/private/` holds participant-level objects: cached cross-validation repeats, out-of-fold predictions, the full final model and SHAP values. **Do not share this folder.**
- `output/verification_report.md` lists every reported value next to the value the code produced.

Seeds, tuning grids and the number of parallel workers (4) are fixed, and ranger and xgboost run single-threaded. Before the main run, the first repeat is run twice and the analysis stops if the two runs differ. Completed repeats are cached, so an interrupted run continues where it stopped.

### Verification and platform

The published results were produced on macOS (Apple silicon, R 4.5.3) with the package versions in `R/package_versions.csv`. On 26 October 2026 these scripts were rerun from scratch on the same platform:

- every out-of-fold prediction in all ten repeats, the selected predictors in all 50 outer folds and the final model's predictions were identical to the published run, and `model/final_model_fit_light.rds` gives the same predictions as the refitted final model;
- `tests/check_against_paper.R` compared 630 values reported in the paper: all values in Tables 1–2 and Supplementary Tables S1–S10 matched, and 3 values in the text differ by rounding (the largest AUROC range of one model across repeats is 0.066, reported as "up to 0.06"; exact Shapley and kernelSHAP mean absolute values for blood and mucus in stool differ by 0.001). The full report is `results/verification_report.md`.

On Linux or Windows, Random Forest results differ slightly even with the same seeds, because ranger draws its random numbers with C++ standard-library routines that are implemented differently on each platform; a few near-tied tuning choices for k-nearest neighbours and the decision tree can also change. In a test run on Linux the same ten predictors and the same final model were selected (AUROC 0.71, 95% CI 0.66 to 0.76), with small differences in the Random Forest calibration slope (1.03) and threshold (0.44). The other algorithms agreed to within floating-point error.

These scripts are a cleaned version of the analysis notebook (V4.10) described in Supplementary Section S8. The notebook also loaded the `dcurves` package for a draft decision-curve plot; no reported number uses it, so it is not required here.

## Using the model

```sh
Rscript model/predict_example.R
```

The model needs ten inputs: age (months), residence, diarrhoea duration (days), straining/tenesmus, antibiotic taken before admission (five categories), blood in stool, mucus in stool, stool white-cell count, stool pus-cell count (midpoints of the recorded bands) and faecal macrophage grade (1.5 for 1–2/HPF, 3.5 for 2–5/HPF). It returns the predicted probability of treatment failure at 48 hours. The pre-specified threshold is 0.43.

**The model is internally validated only.** It has not been externally validated and must not be used for clinical decisions. Supplementary Section S7 of the paper (model card) describes the intended use and the populations in which the model should not be used.

## Data

De-identified participant data are available on reasonable request to researchers who meet icddr,b's criteria for access to confidential data, through the Research Administration of icddr,b (Ms Shiblee Sayeed, shiblee_s@icddrb.org). `data/README.md` lists the expected files and variables.

## Licence

Code: MIT licence (see `LICENSE`). Data: not included; access through icddr,b as above.
