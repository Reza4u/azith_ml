# =============================================================================
# 00_setup.R — packages, analysis-lock settings, output folders, labels
# -----------------------------------------------------------------------------
# Every setting below is fixed in analysis_plan/Analysis_Lock_v1.md and must not
# be changed. The values reproduce the numbers reported in the paper.
# =============================================================================

pkgs <- c(
  # Core
  "tidyverse", "tidymodels", "themis", "vip", "pROC",
  # Models
  "ranger", "xgboost", "glmnet", "kknn", "kernlab", "nnet", "rpart",
  # Utilities
  "discrim", "klaR", "infotheo", "probably", "DALEX", "doFuture", "GPfit", "future.apply",
  # Tables and interpretation
  "gtsummary", "car", "readxl", "kernelshap", "shapviz", "patchwork"
)
# The packages are attached in the order used for the published run, so function
# masking is the same. (The published notebook also attached dcurves for a draft
# decision-curve plot; no reported number uses it, so it is not needed here.)
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  stop("Install these packages first (versions in R/package_versions.csv): ",
       paste(missing_pkgs, collapse = ", "))
}
invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))

# Warn (do not stop) if a key package differs from the version used for the paper
expected_versions <- read.csv("R/package_versions.csv", stringsAsFactors = FALSE)
installed_versions <- vapply(expected_versions$package, function(p)
  tryCatch(as.character(packageVersion(p)), error = function(e) "0"), character(1))
version_diff <- package_version(installed_versions) != package_version(expected_versions$version)
if (getRversion() != "4.5.3") {
  warning("The paper used R 4.5.3; this is R ", getRversion(), ".", call. = FALSE)
}
if (any(version_diff)) {
  warning("Package versions differ from those used for the paper: ",
          paste0(expected_versions$package[version_diff], " ",
                 installed_versions[version_diff], " (paper ",
                 expected_versions$version[version_diff], ")", collapse = "; "),
          ". Results may differ slightly.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# ANALYSIS LOCK v1 (analysis_plan/Analysis_Lock_v1.md)
# -----------------------------------------------------------------------------
SEED            <- 123
N_REPEATS       <- 10                               # repeated nested CV
N_OUTER         <- 5                                # outer folds (performance)
N_INNER         <- 10                               # inner folds (tuning)
REPEAT_SEEDS    <- SEED + 100 * seq_len(N_REPEATS)  # 223, 323, ..., 1123
SEL_FREQ_MIN    <- 0.50                             # final set: selected in >= 50% of outer folds
N_BOOT          <- 2000                             # patient-level bootstrap resamples
N_WORKERS       <- 4                                # fixed, not detectCores()
RUN_REPRO_CHECK <- TRUE                             # run repeat 1 twice and compare first

set.seed(SEED)
plan(multisession, workers = N_WORKERS)
options(ranger.num.threads = 1)

# -----------------------------------------------------------------------------
# Locations. Data are not distributed with the code (see data/README.md).
# -----------------------------------------------------------------------------
DATA_DIR   <- Sys.getenv("AZITH_DATA_DIR", "data")
OUT_DIR    <- Sys.getenv("AZITH_OUT_DIR", "output")
TABLE_DIR  <- file.path(OUT_DIR, "tables")    # aggregate tables (paper numbering)
FIG_DIR    <- file.path(OUT_DIR, "figures")
PRIVATE_DIR <- file.path(OUT_DIR, "private")  # participant-level objects: never share
REPEAT_DIR <- file.path(PRIVATE_DIR, "repeats")
for (d in c(TABLE_DIR, FIG_DIR, PRIVATE_DIR, REPEAT_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
tab_path  <- function(f) file.path(TABLE_DIR, f)
fig_path  <- function(f) file.path(FIG_DIR, f)
priv_path <- function(f) file.path(PRIVATE_DIR, f)

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# -----------------------------------------------------------------------------
# Variable labels.
# `var_labels` are the internal labels used during the published run. They are
# also used to ORDER the final predictor set (see 05_nested_cv.R), and the
# column order changes the Random Forest fit, so they must not be edited.
# `paper_labels` are the labels printed in the paper's tables.
# -----------------------------------------------------------------------------
var_labels <- c(
  "age" = "Age (months)", "sex" = "Sex", "resident_area" = "Residence Area",
  "diarrhea_dur" = "Diarrhea Duration (days)", "fever" = "Fever", "cough" = "Cough",
  "resp_rate" = "Respiratory Rate", "straining" = "Straining/Tenesmus",
  "abx_cat" = "Antibiotic Taken at Home",
  "abx_cat_Ciprofloxacin" = "Antibiotic at Home: Ciprofloxacin",
  "abx_cat_Metronidazole" = "Antibiotic at Home: Metronidazole",
  "abx_cat_Cipro_Metro" = "Antibiotic at Home: Ciprofloxacin + Metronidazole",
  "abx_cat_Other" = "Antibiotic at Home: Other",
  "haz" = "Height-for-Age Z-score (HAZ)", "waz" = "Weight-for-Age Z-score (WAZ)",
  "whz" = "Weight-for-Height Z-score (WHZ)",
  "consist_RE" = "Stool Consistency (RE)", "blood_RE" = "Blood in Stool (RE)",
  "mucus_RE" = "Mucus in Stool (RE)",
  "wbc_numeric" = "WBC Count in Stool", "rbc_numeric" = "RBC Count in Stool",
  "pcell_numeric" = "Pus Cells in Stool", "macro_numeric" = "Macrophages in Stool",
  "is_green" = "Stool Color: Green", "is_red" = "Stool Color: Red/Bloody",
  "is_yellow" = "Stool Color: Yellow", "is_brown" = "Stool Color: Brown",
  "azithromycin_failure" = "Azithromycin Treatment Failure"
)
lab <- function(v) if (v %in% names(var_labels)) unname(var_labels[[v]]) else v

paper_labels <- c(
  "azithromycin_failure" = "Azithromycin treatment failure (outcome)",
  "age" = "Age (months)", "sex" = "Sex", "resident_area" = "Residence",
  "diarrhea_dur" = "Diarrhoea duration (days)", "fever" = "Fever", "cough" = "Cough",
  "resp_rate" = "Respiratory rate", "straining" = "Straining/tenesmus",
  "haz" = "Height-for-age Z-score", "waz" = "Weight-for-age Z-score",
  "whz" = "Weight-for-height Z-score", "consist_RE" = "Stool consistency",
  "blood_RE" = "Blood in stool", "mucus_RE" = "Mucus in stool",
  "abx_cat" = "Pre-admission antibiotic category",
  "wbc_numeric" = "Stool white-cell count", "macro_numeric" = "Faecal macrophage grade",
  "pcell_numeric" = "Stool pus-cell count", "rbc_numeric" = "Stool red-cell count",
  "is_green" = "Stool colour: green", "is_red" = "Stool colour: red/bloody",
  "is_yellow" = "Stool colour: yellow", "is_brown" = "Stool colour: brown"
)
plab <- function(v) unname(ifelse(v %in% names(paper_labels), paper_labels[v], v))

paper_model_names <- c(
  "Elastic Net" = "Elastic Net", "Random Forest" = "Random Forest", "XGBoost" = "XGBoost",
  "KNN" = "k-nearest neighbours", "SVM" = "Support vector machine",
  "Decision Tree" = "Decision tree", "Naive Bayes" = "Naïve Bayes",
  "Neural Network" = "Neural network"
)

cat("Setup complete:", N_REPEATS, "repeats x", N_OUTER, "outer x", N_INNER,
    "inner folds;", N_WORKERS, "workers; outputs in", OUT_DIR, "\n")
