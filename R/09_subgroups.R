# =============================================================================
# 09_subgroups.R — performance of the final model by subgroup (Supplementary Table S3)
# -----------------------------------------------------------------------------
# AUROC = mean across repeats of the within-subgroup AUROC; 95% CI = patient-level
# bootstrap of that mean (children resampled within the subgroup). Calibration
# slope, Brier score, sensitivity and specificity use all repeats stacked.
# =============================================================================

M_best <- pred_mat[[best_model]]                      # n x repeats, ml_df row order
chars  <- ml_df |> dplyr::select(sex, age, resident_area, haz, waz, whz)

subgroup_row <- function(rows, domain, label) {
  rows    <- which(rows)
  is_fail <- y_all[rows]
  set.seed(SEED)
  bt <- replicate(N_BOOT, {
    i <- rows[sample.int(length(rows), length(rows), replace = TRUE)]
    if (length(unique(y_all[i])) < 2) NA_real_ else mean_auc(M_best, i)
  })
  y_st <- rep(is_fail, N_REPEATS)                     # stacked over repeats (column-major)
  p_st <- as.vector(M_best[rows, , drop = FALSE])
  pf   <- p_st >= THR
  tibble(Domain = domain, Subgroup = label, n = length(rows), Events = sum(is_fail),
         AUROC = mean_auc(M_best, rows),
         CI_low = quantile(bt, 0.025, na.rm = TRUE), CI_high = quantile(bt, 0.975, na.rm = TRUE),
         `Calibration slope` = cal_slope_fn(y_st, p_st),
         `Brier score` = mean((as.integer(y_st) - p_st)^2),
         Sensitivity = sum(pf & y_st) / sum(y_st),
         Specificity = sum(!pf & !y_st) / sum(!y_st))
}

subgroups <- bind_rows(
  subgroup_row(rep(TRUE, n_pat),                    "Overall", "All participants"),
  subgroup_row(chars$sex == "Male",                 "Sex", "Male"),
  subgroup_row(chars$sex == "Female",               "Sex", "Female"),
  subgroup_row(chars$age < 36,                      "Age band", "24-35 months"),
  subgroup_row(chars$age >= 36 & chars$age < 48,    "Age band", "36-47 months"),
  subgroup_row(chars$age >= 48,                     "Age band", "48-59 months"),
  subgroup_row(chars$resident_area == "Urban Area", "Residence", "Urban"),
  subgroup_row(chars$resident_area == "Rural Area", "Residence", "Rural"),
  subgroup_row(chars$haz <  -2, "Nutritional status", "Stunted (HAZ < -2)"),
  subgroup_row(chars$haz >= -2, "Nutritional status", "Not stunted"),
  subgroup_row(chars$waz <  -2, "Nutritional status", "Underweight (WAZ < -2)"),
  subgroup_row(chars$waz >= -2, "Nutritional status", "Not underweight"),
  subgroup_row(chars$whz <  -2, "Nutritional status", "Wasted (WHZ < -2)"),
  subgroup_row(chars$whz >= -2, "Nutritional status", "Not wasted")
)

table_s3 <- subgroups |>
  transmute(Domain, Subgroup, n, Events,
            `AUROC (95% CI)` = sprintf("%.2f (%.2f-%.2f)", round(AUROC, 3), CI_low, CI_high),
            `Calibration slope` = round(`Calibration slope`, 2),
            `Brier score` = round(`Brier score`, 3),
            Sensitivity = round(Sensitivity, 2),
            Specificity = round(Specificity, 2))
write.csv(table_s3, tab_path("TableS3_subgroup_performance.csv"), row.names = FALSE)

event_rate_by_residence <- tapply(y_all, ml_df$resident_area, mean)
