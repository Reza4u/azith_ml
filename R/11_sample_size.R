# =============================================================================
# 11_sample_size.R — events per predictor and Riley et al. minimum sample size
# -----------------------------------------------------------------------------
# Paper: Methods (Sample size); Supplementary Sections S1.2, S3 (PROBAST+AI 4.1)
# and S6.4. Uses the pmsampsize package (Riley et al., BMJ 2020).
# =============================================================================

events_per_predictor  <- n_events / length(FINAL_SET)
events_per_parameter  <- n_events / n_model_columns
events_per_candidate_column <- n_events / data_facts$cols_after_dummy

if (requireNamespace("pmsampsize", quietly = TRUE)) {
  riley <- pmsampsize::pmsampsize(type = "b",
                                  cstatistic = round(point_auc[best_model], 3),
                                  parameters = n_model_columns,
                                  prevalence = round(prev, 3),
                                  shrinkage  = 0.9)
  riley_n      <- riley$sample_size
  riley_events <- round(riley$events)
  riley_r2cs   <- riley$results_table["Final", "CS_Rsq"]
} else {
  message("Package 'pmsampsize' not installed: Riley sample size skipped.")
  riley_n <- riley_events <- riley_r2cs <- NA_real_
}
cat(sprintf("Events per predictor %.1f; per model column %.1f; per candidate column %.1f; Riley n = %s\n",
            events_per_predictor, events_per_parameter, events_per_candidate_column, riley_n))
