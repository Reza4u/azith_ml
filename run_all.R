# =============================================================================
# run_all.R — reproduce every number, table and figure in the paper
#
#   Reza N, Islam MR, Al Fidah MF, Amin R, Chisti MJ, Ahmed T, Nuzhat S.
#   Development of a multivariable machine-learning model to predict
#   azithromycin treatment failure in children with suspected shigellosis:
#   a prognostic modelling study. BMJ Health & Care Informatics, 2026.
#
# Usage (from the repository root, with the data files in data/):
#   Rscript run_all.R
#   python3 python/figure3_beeswarm.py      # Figure 3
#   Rscript tests/check_against_paper.R     # compare with the published numbers
#
# Scripts share one R session and must run in this order. Runtime: about 30 min
# on a 4-core laptop (repeated nested CV) plus about 15 min for exact Shapley
# values. Completed CV repeats are cached and skipped if the run is restarted.
# =============================================================================

t_start <- Sys.time()
scripts <- c(
  "R/00_setup.R",               # packages, analysis-lock settings, labels
  "R/01_data.R",                # load data, pre-specified coding
  "R/02_descriptive_tables.R",  # Table 1; Supplementary Tables S1, S5a-b, S6a-b
  "R/03_feature_selection.R",   # feature-selection function
  "R/04_models.R",              # eight algorithms and fixed grids
  "R/05_nested_cv.R",           # repeated nested CV; Supplementary Table S2
  "R/06_performance.R",         # Table 2; Supplementary Tables S4a-d, S7, S8; final-model rule
  "R/07_decision_curve.R",      # Supplementary Table S9
  "R/08_final_model_shap.R",    # final model; Supplementary Table S10; exact Shapley values
  "R/09_subgroups.R",           # Supplementary Table S3
  "R/10_figures.R",             # Figure 2 and Supplementary Figure S2
  "R/11_sample_size.R",         # events per predictor; Riley et al. criteria
  "R/12_key_numbers.R"          # numbers quoted in the text
)
for (s in scripts) {
  cat("\n==== ", s, " ====\n", sep = "")
  source(s, echo = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
cat("\nDone in", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), "min.",
    "Tables:", TABLE_DIR, "| Figures:", FIG_DIR, "\n")
