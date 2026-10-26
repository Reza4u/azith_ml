# =============================================================================
# predict_example.R — predicted probability of azithromycin treatment failure
#                     at 48 hours for a new child
# -----------------------------------------------------------------------------
# Research use only: the model is internally validated, not externally validated.
# Run from the repository root:   Rscript model/predict_example.R
# Needs tidymodels and ranger (versions in R/package_versions.csv).
# =============================================================================
suppressMessages({library(tidymodels); library(ranger)})

fit <- readRDS("model/final_model_fit_light.rds")   # contains no participant data

# Expected inputs: names, types and allowed categories
ptype <- hardhat::extract_mold(fit)$blueprint$ptypes$predictors
for (v in names(ptype)) {
  cat(sprintf("%-14s %s\n", v, if (is.factor(ptype[[v]]))
    paste("one of:", paste(levels(ptype[[v]]), collapse = " | ")) else "numeric"))
}

# One example child: replace the values below.
new_child <- tibble(
  age           = 36,             # months (24-59)
  resident_area = "Urban Area",   # "Urban Area" or "Rural Area"
  diarrhea_dur  = 3,              # days
  straining     = "Yes",          # straining/tenesmus: "Yes" or "No"
  abx_cat       = "None",         # None | Ciprofloxacin | Metronidazole | Cipro_Metro | Other
  blood_RE      = "Small",        # Absent | Small | Moderate | Plenty
  mucus_RE      = "Moderate",     # Small | Moderate | Plenty
  wbc_numeric   = 13,             # stool white cells, band midpoint (1-5 -> 3, 6-10 -> 8, 11-15 -> 13,
                                  #   16-20 -> 18, 21-25 -> 23, 26-30 -> 28, 31+ -> 35)
  pcell_numeric = 35,             # stool pus cells, band midpoint (as above)
  macro_numeric = 1.5             # faecal macrophages: 1.5 = 1-2/HPF, 3.5 = 2-5/HPF
)
# match the column order and factor levels used in training
new_child <- new_child[, names(ptype)]
for (v in names(ptype)) if (is.factor(ptype[[v]])) {
  stopifnot(all(new_child[[v]] %in% levels(ptype[[v]])))
  new_child[[v]] <- factor(new_child[[v]], levels = levels(ptype[[v]]))
}

p <- predict(fit, new_child, type = "prob")$.pred_Yes
cat(sprintf("\nPredicted probability of treatment failure: %.2f (pre-specified threshold 0.43)\n", p))
