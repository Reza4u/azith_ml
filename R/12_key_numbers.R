# =============================================================================
# 12_key_numbers.R — numbers quoted in the text that are not in a table
# -----------------------------------------------------------------------------
# Writes output/tables/key_numbers.csv (name, value). tests/check_against_paper.R
# compares these and all tables with the values printed in the paper.
# =============================================================================

f2 <- function(x) sprintf("%.2f", x)
f3 <- function(x) sprintf("%.3f", x)
fb <- auroc_summary |> filter(Model == best_model)
cb <- cal_summary   |> filter(Model == best_model)
cl <- cls_repeat    |> filter(Model == best_model)
rf_params <- tuned_params |> filter(Model == "Random Forest")
el <- vs_leader |> filter(Model == best_model)
nb_at <- function(pt) net_benefit(best_pred$y, best_pred$.pred_Yes, pt)
shap_cmp <- if (exists("shap_exact_summary")) {
  dplyr::inner_join(shap_kernel_summary |> dplyr::select(var, k = `Mean absolute SHAP`),
                    shap_exact_summary  |> dplyr::select(var, e = `Mean absolute SHAP`), by = "var")
} else NULL

key_values <- c(
  n_children =                    n_pat,
  n_failures =                    n_events,
  pct_failures =                  sprintf("%.1f", 100 * n_events / n_pat),
  n_success =                     n_pat - n_events,
  pct_success =                   sprintf("%.1f", 100 * (n_pat - n_events) / n_pat),
  n_candidate_variables =         data_facts$n_candidates,
  n_missing_values =              data_facts$n_missing,
  resp_rate_values_replaced =     paste(data_facts$rr_values_replaced, collapse = " and "),
  resp_rate_replacement =         data_facts$rr_fill,
  columns_after_dummy =           data_facts$cols_after_dummy,
  columns_after_nzv =             data_facts$cols_after_nzv,
  nzv_removed =                   paste(sort(data_facts$nzv_removed), collapse = ", "),
  cipro_containing_failures =     sprintf("%d/%d (%.1f%%)", cipro_containing["failure_n"],
                                           cipro_containing["failure_N"],
                                           100 * cipro_containing["failure_n"] / cipro_containing["failure_N"]),
  cipro_containing_successes =    sprintf("%d/%d (%.1f%%)", cipro_containing["success_n"],
                                           cipro_containing["success_N"],
                                           100 * cipro_containing["success_n"] / cipro_containing["success_N"]),
  antibiotic_fisher_p =           f3(cat_fisher$p.value),
  failures_one_component =        sprintf("%d (%.1f%%)", components_met[["1"]],
                                           100 * components_met[["1"]] / sum(components_met)),
  n_outer_selections =            n_selections,
  n_final_predictors =            length(FINAL_SET),
  final_predictors =              paste(sort(FINAL_SET), collapse = ", "),
  events_per_predictor =          sprintf("%.1f", events_per_predictor),
  model_columns_final =           n_model_columns,
  events_per_model_column =       sprintf("%.1f", events_per_parameter),
  events_per_candidate_column =   sprintf("%.1f", events_per_candidate_column),
  riley_r2cs =                    sprintf("%.2f", riley_r2cs),
  riley_sample_size =             riley_n,
  riley_events =                  riley_events,
  leader =                        leader,
  candidates =                    paste(sort(candidates), collapse = ", "),
  final_model =                   best_model,
  final_model_repeats_won =       fb$`Repeats won`,
  final_auroc =                   f2(fb$AUROC),
  final_auroc_ci =                sprintf("%.2f-%.2f", fb$CI_low, fb$CI_high),
  final_auroc_range =             sprintf("%.2f-%.2f", fb$Min, fb$Max),
  final_auroc_sd =                f3(fb$SD_across_repeats),
  leader_minus_final_dauroc =     f3(el$dAUROC),
  leader_minus_final_ci =         sprintf("%.3f to %.3f", el$CI_low, el$CI_high),
  leader_minus_final_p =          f2(el$p),
  final_cal_intercept =           f2(cb$Cal_int),
  final_cal_slope =               f2(cb$Cal_slope),
  final_brier =                   f3(cb$Brier),
  final_hl_p_median =             f2(cb$HL_p_median),
  final_threshold =               f2(THR),
  final_sensitivity =             fmt_mr(cl$Sensitivity),
  final_specificity =             fmt_mr(cl$Specificity),
  final_ppv =                     f2(mean(cl$PPV)),
  final_npv =                     f2(mean(cl$NPV)),
  nb_final_40 =                   f3(nb_at(0.40)),
  nb_treat_all_40 =               f3(treat_all_nb(0.40)),
  incremental_nb_40 =             f3(nb_at(0.40) - treat_all_nb(0.40)),
  dca_beats_both_range =          sprintf("%.0f%%-%.0f%%", 100 * dca_range["from"], 100 * dca_range["to"]),
  dca_below_treat_all =           paste(sprintf("%.0f%%", 100 * below_all), collapse = ", "),
  rf_min_n_20_outer_folds =       sprintf("%d of %d", sum(rf_params$min_n == "20"), nrow(rf_params)),
  final_hyperparameters =         paste(names(best_params)[names(best_params) != ".config"],
                                         best_params[names(best_params) != ".config"],
                                         sep = "=", collapse = ", "),
  max_auroc_spread_one_model =    f2(max(auroc_summary$Max - auroc_summary$Min)),
  brier_noninformative =          f3(prev * (1 - prev)),
  event_rate_urban =              sprintf("%.1f%%", 100 * event_rate_by_residence[["Urban Area"]]),
  event_rate_rural =              sprintf("%.1f%%", 100 * event_rate_by_residence[["Rural Area"]]),
  overall_stacked_cal_slope =     f2(subgroups$`Calibration slope`[1]),
  shap_exact_vs_kernel_max_diff = if (is.null(shap_cmp)) NA else
                                     sprintf("%.4f", max(abs(round(shap_cmp$k, 3) - round(shap_cmp$e, 3)))),
  shap_exact_blood_mucus =        if (is.null(shap_cmp)) NA else
                                     paste(f3(shap_cmp$e[shap_cmp$var %in% c("blood_RE", "mucus_RE")]),
                                           collapse = " and ")
)
key_numbers <- tibble(name = names(key_values), value = unname(as.character(key_values)))
write.csv(key_numbers, tab_path("key_numbers.csv"), row.names = FALSE)
cat("Key numbers written to", tab_path("key_numbers.csv"), "\n")
