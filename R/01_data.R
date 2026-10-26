# =============================================================================
# 01_data.R — load the analysis dataset and apply the pre-specified coding
# -----------------------------------------------------------------------------
# Input : data/azithromycin_pediatric_clinical_data.rds (not distributed)
# Output: `df` (raw), `ml_df` (329 children, outcome + 23 candidate variables)
# Paper : Methods (Predictors); Supplementary Section S1.1
# =============================================================================

data_file <- file.path(DATA_DIR, "azithromycin_pediatric_clinical_data.rds")
if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file,
       "\nParticipant data are not public. See data/README.md for how to request them.")
}
df <- readRDS(data_file)
stopifnot("antibiotic_home" %in% names(df))

midpoint <- function(x) dplyr::case_when(
  x == "1-5" ~ 3, x == "6-10" ~ 8, x == "11-15" ~ 13, x == "16-20" ~ 18,
  x == "21-25" ~ 23, x == "26-30" ~ 28, x == "31+" ~ 35, .default = NA_real_)

ml_df <- df |>
  filter(!is.na(azithromycin_failure)) |>
  dplyr::select(azithromycin_failure, age, sex, resident_area, diarrhea_dur, fever, cough,
                resp_rate, straining, haz, waz, whz,
                antibiotic_home,
                color_RE, consist_RE, blood_RE, mucus_RE, rbc_RE, pcell_RE, wbc_RE, macro_RE) |>
  mutate(
    # Respiratory rate outside 15-70 breaths/min is a data-entry error and is
    # replaced with the mean of the remaining observations (overall mean, so no
    # outcome information enters a predictor).
    resp_rate = replace(
      resp_rate,
      !(resp_rate >= 15 & resp_rate <= 70),
      round(mean(resp_rate[resp_rate >= 15 & resp_rate <= 70], na.rm = TRUE))),
    azithromycin_failure = factor(azithromycin_failure, levels = c("Yes", "No")),

    # Pre-admission antibiotic: five categories pre-specified in the statistical
    # analysis plan; the grouping was fixed without reference to the outcome.
    abx_cat = factor(case_when(
      antibiotic_home == "No Antibiotic"                 ~ "None",
      antibiotic_home == "Ciprofloxacin"                 ~ "Ciprofloxacin",
      antibiotic_home == "Metronidazole"                 ~ "Metronidazole",
      antibiotic_home == "Ciprofloxacin + Metronidazole" ~ "Cipro_Metro",
      !is.na(antibiotic_home)                            ~ "Other"),
      levels = c("None", "Ciprofloxacin", "Metronidazole", "Cipro_Metro", "Other")),

    # Stool microscopy bands -> numeric midpoints
    wbc_numeric   = midpoint(wbc_RE),
    macro_numeric = case_when(macro_RE == "1-2" ~ 1.5, macro_RE == "2-5" ~ 3.5,
                              .default = NA_real_),
    pcell_numeric = midpoint(pcell_RE),
    rbc_numeric   = midpoint(rbc_RE),

    # Free-text stool colour -> four binary indicators
    color_clean = str_squish(str_replace_all(tolower(as.character(color_RE)), "[^a-z\\s]", " ")),
    is_green  = as.integer(str_detect(color_clean, "green")),
    is_red    = as.integer(str_detect(color_clean, "red|blood|stained")),
    is_yellow = as.integer(str_detect(color_clean, "yellow")),
    is_brown  = as.integer(str_detect(color_clean, "brown"))
  ) |>
  dplyr::select(-color_clean, -color_RE, -rbc_RE, -pcell_RE, -wbc_RE, -macro_RE, -antibiotic_home)

stopifnot(!anyNA(ml_df$abx_cat), nrow(ml_df) == 329, !anyNA(ml_df))
candidate_vars <- setdiff(names(ml_df), "azithromycin_failure")
n_pat    <- nrow(ml_df)
y_all    <- ml_df$azithromycin_failure == "Yes"
n_events <- sum(y_all)

# Respiratory-rate correction actually applied (Supplementary Section S1.1)
rr_raw <- df$resp_rate[!is.na(df$azithromycin_failure)]
rr_out <- rr_raw[!(rr_raw >= 15 & rr_raw <= 70)]
rr_fill <- round(mean(rr_raw[rr_raw >= 15 & rr_raw <= 70], na.rm = TRUE))

# Model columns after dummy coding and after the near-zero-variance filter,
# estimated on the complete data (Supplementary Section S1.1)
tmp_dummy <- recipe(azithromycin_failure ~ ., data = ml_df) |>
  step_zv(all_predictors()) |> step_dummy(all_nominal_predictors()) |>
  prep() |> bake(new_data = NULL)
tmp_nzv <- recipe(azithromycin_failure ~ ., data = ml_df) |>
  step_zv(all_predictors()) |> step_dummy(all_nominal_predictors()) |>
  step_nzv(all_predictors()) |> prep() |> bake(new_data = NULL)
data_facts <- list(
  n = n_pat, events = n_events, n_candidates = length(candidate_vars),
  n_missing = sum(is.na(ml_df)),
  rr_n_replaced = length(rr_out), rr_values_replaced = sort(rr_out), rr_fill = rr_fill,
  cols_after_dummy = ncol(tmp_dummy) - 1, cols_after_nzv = ncol(tmp_nzv) - 1,
  nzv_removed = setdiff(names(tmp_dummy), names(tmp_nzv)))
rm(tmp_dummy, tmp_nzv)

cat("Data:", n_pat, "children,", n_events, "failures,", length(candidate_vars),
    "candidate variables;", data_facts$rr_n_replaced, "respiratory-rate value(s) replaced with",
    rr_fill, "\n")
