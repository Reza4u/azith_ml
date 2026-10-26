# =============================================================================
# 07_decision_curve.R — net benefit (Supplementary Table S9) and threshold range
# -----------------------------------------------------------------------------
# Net benefit is computed on all repeats stacked; because every repeat contains
# every child once, this equals the mean net benefit across repeats.
# =============================================================================

prev <- mean(y_all)
treat_all_nb <- function(pt) prev - (1 - prev) * (pt / (1 - pt))

key_thresholds <- c(0.10, 0.20, 0.30, 0.40, 0.50)
table_s9 <- purrr::map_dfr(key_thresholds, function(pt) {
  ref <- max(treat_all_nb(pt), 0)
  bind_rows(
    tibble(`Threshold probability` = paste0(pt * 100, "%"), `Strategy / model` = "Treat all",
           `Net benefit` = sprintf("%.4f", treat_all_nb(pt)),
           `Incremental net benefit vs treat all` = "Reference"),
    tibble(`Threshold probability` = "", `Strategy / model` = "Treat none",
           `Net benefit` = "0.0000", `Incremental net benefit vs treat all` = "-"),
    purrr::map_dfr(top3_models, function(m) {
      d  <- dplyr::filter(pr, Model == m)
      nb <- net_benefit(d$y, d$.pred_Yes, pt)
      tibble(`Threshold probability` = "", `Strategy / model` = pm(m),
             `Net benefit` = sprintf("%.4f", nb),
             `Incremental net benefit vs treat all` = sprintf("%.4f", nb - ref))
    }))
})
write.csv(table_s9, tab_path("TableS9_dca_net_benefit.csv"), row.names = FALSE)

# Final model on a 1% grid: where does it beat treat-all and treat-none?
best_pred <- dplyr::filter(pr, Model == best_model)
grid_pt   <- seq(0.01, 0.95, by = 0.01)
dca_fine  <- tibble(threshold = grid_pt,
                    model     = sapply(grid_pt, function(pt) net_benefit(best_pred$y, best_pred$.pred_Yes, pt)),
                    treat_all = sapply(grid_pt, treat_all_nb)) |>
  mutate(vs_all = model - treat_all, beats_all = vs_all > 0, beats_none = model > 0,
         beats_both = beats_all & beats_none)
write.csv(dca_fine, tab_path("dca_fine_grid_final_model.csv"), row.names = FALSE)

# longest run of consecutive thresholds at which the model beats both strategies
runs <- rle(dca_fine$beats_both)
end_i <- cumsum(runs$lengths); start_i <- end_i - runs$lengths + 1
k <- which(runs$values)[which.max(runs$lengths[runs$values])]
dca_range <- c(from = dca_fine$threshold[start_i[k]], to = dca_fine$threshold[end_i[k]])
below_all <- dca_fine$threshold[dca_fine$threshold < dca_range["from"] & !dca_fine$beats_all]

cat(sprintf("Final model beats treat-all and treat-none at every threshold from %.0f%% to %.0f%%\n",
            100 * dca_range["from"], 100 * dca_range["to"]))
