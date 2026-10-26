# =============================================================================
# 10_figures.R — Figure 2 (A ROC, B decision curve, C calibration) and
#                Supplementary Figure S2 (calibration, all eight models)
# -----------------------------------------------------------------------------
# Built from the out-of-fold predictions of all ten repeats. Figure 3 (SHAP
# beeswarm) is drawn by python/figure3_beeswarm.py from the exact Shapley values.
# Figure 1 (participant flow) and Figure S1 (process diagram) are diagrams, not
# computed from data; their Mermaid sources are in figures_src/.
# =============================================================================

FONT <- if ("Liberation Sans" %in% systemfonts::system_fonts()$family) "Liberation Sans" else "sans"
p_fig <- cv_predictions |> mutate(y = as.integer(azithromycin_failure == "Yes"))
auc_fig <- auroc_summary |> mutate(across(c(AUROC, CI_low, CI_high), ~ round(.x, 3)))
cal_fig <- per_repeat |> group_by(Model) |>
  summarise(a = mean(Cal_int), b = mean(Cal_slope), brier = mean(Brier),
            hl = median(HL_p, na.rm = TRUE), .groups = "drop")

cols <- c("Random Forest" = "#2E86AB", "Elastic Net" = "#F18F01", "XGBoost" = "#A23B72",
          "SVM" = "#6B2D5C", "Neural Network" = "#95C623", "Naive Bayes" = "#3B1F2B",
          "KNN" = "#1B998B", "Decision Tree" = "#C73E1D")
show <- top3_models                                    # final model + next two by AUROC
nice <- setNames(paper_model_names[names(cols)], names(cols))
nice[best_model] <- paste0(nice[best_model], " (final model)")
nice["SVM"] <- "Support vector machine"

th <- theme_minimal(base_size = 11, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.3),
        axis.line = element_line(colour = "#404040", linewidth = 0.5),
        axis.ticks = element_line(colour = "#404040"),
        axis.title = element_text(face = "bold"), plot.tag = element_text(face = "bold", size = 15),
        legend.background = element_rect(fill = "white", colour = "#808080", linewidth = 0.3),
        strip.text = element_text(face = "bold"), plot.margin = margin(8, 10, 8, 8))

# ---- A: ROC (all repeats stacked) ------------------------------------------------
roc_df <- p_fig |> filter(Model %in% show) |> group_by(Model) |>
  group_modify(function(d, k) {
    th_ <- c(Inf, sort(unique(d$.pred_Yes), decreasing = TRUE), -Inf)
    pos <- sum(d$y); neg <- sum(1 - d$y)
    tibble(fpr = sapply(th_, function(t) sum(d$.pred_Yes >= t & d$y == 0) / neg),
           tpr = sapply(th_, function(t) sum(d$.pred_Yes >= t & d$y == 1) / pos))
  }) |> ungroup()
roc_lab <- auc_fig |> filter(Model %in% show) |>
  mutate(label = sprintf("%s: AUROC %.2f (%.2f–%.2f)", nice[Model], AUROC, CI_low, CI_high))
roc_df <- roc_df |> left_join(roc_lab[, c("Model", "label")], by = "Model") |>
  mutate(label = factor(label, levels = roc_lab$label[match(show, roc_lab$Model)]))
pA <- ggplot(roc_df, aes(fpr, tpr, colour = label)) +
  geom_abline(linetype = 2, colour = "grey55") +
  geom_path(linewidth = 0.9) +
  scale_colour_manual(values = setNames(cols[show], levels(roc_df$label))) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(x = "1 − specificity", y = "Sensitivity", colour = NULL, tag = "A") +
  th + theme(legend.position = c(0.6, 0.15), legend.text = element_text(size = 8.5),
             legend.key.height = unit(0.45, "cm"))

# ---- calibration panels ------------------------------------------------------------
cal_bins <- function(ms) p_fig |> filter(Model %in% ms) |>
  mutate(bin = cut(.pred_Yes, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) |>
  group_by(Model, bin) |>
  summarise(pred = 100 * mean(.pred_Yes), obs = 100 * mean(y), n = n(), .groups = "drop") |>
  filter(n >= 10)                                      # suppress sparse bins
cal_lab <- function(ms) cal_fig |> filter(Model %in% ms) |>
  mutate(txt = sprintf("Intercept %.2f\nSlope %.2f\nBrier %.3f\nHL p %.2f", a, b, brier, hl),
         txt = sub("HL p 0.00", "HL p <0.01", txt), txt = gsub("-", "−", txt))
mk_cal <- function(ms, ncol, tag = NULL) {
  bins  <- cal_bins(ms) |> mutate(Model = factor(nice[Model], levels = nice[ms]))
  labs_ <- cal_lab(ms)  |> mutate(Model = factor(nice[Model], levels = nice[ms]))
  ggplot(bins, aes(pred, obs, colour = Model)) +
    geom_abline(linetype = 2, colour = "grey55") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    geom_text(data = labs_, aes(x = 3, y = 97, label = txt), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 2.9, lineheight = 0.95, family = FONT) +
    facet_wrap(~Model, ncol = ncol) +
    scale_colour_manual(values = setNames(cols[ms], nice[ms]), guide = "none") +
    scale_x_continuous(breaks = seq(0, 100, 25)) + scale_y_continuous(breaks = seq(0, 100, 25)) +
    coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
    labs(x = "Mean predicted probability (%)", y = "Observed failure rate (%)", tag = tag) + th
}
pC <- mk_cal(show, 3, "C")

# ---- B: decision curve (unsmoothed) ------------------------------------------------
pts <- seq(0.01, 0.80, by = 0.01)
nb_curve <- bind_rows(
  lapply(show, function(m) {
    d <- filter(p_fig, Model == m); n <- nrow(d)
    tibble(Strategy = nice[m], pt = pts,
           nb = sapply(pts, function(t) sum(d$.pred_Yes >= t & d$y == 1) / n -
                                       sum(d$.pred_Yes >= t & d$y == 0) / n * t / (1 - t)))
  }),
  tibble(Strategy = "Treat all", pt = pts, nb = prev - (1 - prev) * pts / (1 - pts)),
  tibble(Strategy = "Treat none", pt = pts, nb = 0)) |>
  mutate(Strategy = factor(Strategy, levels = c(nice[show], "Treat all", "Treat none")))
pB <- ggplot(nb_curve, aes(pt * 100, nb, colour = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c(setNames(cols[show], nice[show]),
                                 "Treat all" = "grey25", "Treat none" = "grey65")) +
  scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dotted")) +
  coord_cartesian(xlim = c(0, 80), ylim = c(-0.05, 0.42), expand = FALSE) +
  scale_x_continuous(breaks = seq(0, 80, 20)) +
  labs(x = "Threshold probability (%)", y = "Net benefit", colour = NULL, linetype = NULL, tag = "B") +
  th + theme(legend.position = c(0.7, 0.78), legend.text = element_text(size = 8.5),
             legend.key.width = unit(0.9, "cm"), legend.key.height = unit(0.4, "cm"))

fig2 <- (pA | pB) / pC + plot_layout(heights = c(1, 0.78))
ggsave(fig_path("Figure_2.png"), fig2, width = 9, height = 8.8, dpi = 300, bg = "white")
ggsave(fig_path("Figure_2.tiff"), fig2, width = 9, height = 8.8, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")

figS2 <- mk_cal(auc_fig$Model, 4)
ggsave(fig_path("Figure_S2.png"), figS2, width = 11, height = 6.2, dpi = 300, bg = "white")
ggsave(fig_path("Figure_S2.tiff"), figS2, width = 11, height = 6.2, dpi = 300, bg = "white",
       device = "tiff", compression = "lzw")
cat("Figure 2 and Figure S2 written to", FIG_DIR, "\n")
