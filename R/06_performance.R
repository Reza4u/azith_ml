# =============================================================================
# 06_performance.R — performance across repeats and the final-model selection rule
# -----------------------------------------------------------------------------
# Every measure is computed on the pooled out-of-fold predictions of one repeat
# (329 children) and then summarised across the 10 repeats.
# Paper: Table 2; Supplementary Tables S4a-S4d, S7, S8; Section S1.5-S1.6
# =============================================================================

# ---- helpers ------------------------------------------------------------------
auc_fast <- function(y, p) {            # y logical (TRUE = failure); Mann-Whitney AUROC
  n1 <- sum(y); n0 <- length(y) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
logit_clamp   <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
cal_intercept <- function(y, p) unname(coef(glm(as.integer(y) ~ offset(logit_clamp(p)), family = binomial))[1])
cal_slope_fn  <- function(y, p) unname(coef(glm(as.integer(y) ~ logit_clamp(p), family = binomial))[2])
hl_pvalue <- function(y, p) tryCatch({
  y <- as.integer(y)
  brks <- unique(quantile(p, probs = seq(0, 1, 0.1), type = 1))
  g    <- cut(p, breaks = brks, include.lowest = TRUE)
  obs  <- tapply(y, g, sum); ex <- tapply(p, g, sum); n <- tapply(y, g, length)
  den  <- ex * (1 - ex / n); keep <- !is.na(den) & den > 0
  stat <- sum((obs[keep] - ex[keep])^2 / den[keep]); dfh <- sum(keep) - 2
  if (dfh < 1) NA_real_ else 1 - pchisq(stat, df = dfh)
}, error = function(e) NA_real_)
net_benefit <- function(y, p, pt) {
  n <- length(y)
  sum(p >= pt & y) / n - (sum(p >= pt & !y) / n) * (pt / (1 - pt))
}
mean_nb <- function(y, p, lo = 0.30, hi = 0.40) {
  mean(sapply(seq(lo, hi, by = 0.01), function(pt) net_benefit(y, p, pt)))
}
youden_threshold <- function(y, p) {
  roc_obj <- pROC::roc(factor(ifelse(y, "Yes", "No"), levels = c("No", "Yes")), p,
                       levels = c("No", "Yes"), direction = "<", quiet = TRUE)
  cc <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = "threshold")
  as.numeric(unlist(cc))[1]
}
fmt_mr <- function(x, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), mean(x), min(x), max(x))
f3s    <- function(x) sub("^-(0\\.0+)$", "\\1", sprintf("%.3f", x))   # no "-0.000"
fmt_p  <- function(p) ifelse(p < 1 / N_BOOT, sprintf("<%.4f", 1 / N_BOOT), sprintf("%.3f", p))

pr <- cv_predictions |> mutate(y = azithromycin_failure == "Yes")

# ---- per-repeat measures --------------------------------------------------------
per_repeat <- pr |>
  group_by(Model, Repeat) |>
  summarise(AUROC      = auc_fast(y, .pred_Yes),
            Cal_int    = cal_intercept(y, .pred_Yes),
            Cal_slope  = cal_slope_fn(y, .pred_Yes),
            Brier      = mean((as.integer(y) - .pred_Yes)^2),
            HL_p       = hl_pvalue(y, .pred_Yes),
            NB_30_40   = mean_nb(y, .pred_Yes),
            Youden_thr = youden_threshold(y, .pred_Yes),
            .groups = "drop")

# fixed threshold per model = median of the per-repeat Youden thresholds
model_thresholds <- per_repeat |>
  group_by(Model) |>
  summarise(Threshold = median(Youden_thr), .groups = "drop")
THR_LOOKUP <- setNames(model_thresholds$Threshold, model_thresholds$Model)

# classification measures per repeat at the fixed threshold
cls_repeat <- pr |>
  mutate(pos = .pred_Yes >= THR_LOOKUP[Model]) |>
  group_by(Model, Repeat) |>
  summarise(TP = sum(pos & y), TN = sum(!pos & !y), FP = sum(pos & !y), FN = sum(!pos & y),
            .groups = "drop") |>
  mutate(Accuracy    = (TP + TN) / (TP + TN + FP + FN),
         Sensitivity = TP / (TP + FN),
         Specificity = TN / (TN + FP),
         PPV         = ifelse(TP + FP > 0, TP / (TP + FP), NA_real_),
         NPV         = ifelse(TN + FN > 0, TN / (TN + FN), NA_real_),
         F1          = 2 * TP / (2 * TP + FP + FN))

# AUROC by repeat and the model ranked first in each repeat
repeat_winner <- per_repeat |>
  group_by(Repeat) |> slice_max(AUROC, n = 1, with_ties = FALSE) |> ungroup() |>
  dplyr::select(Repeat, Winner = Model)
win_counts <- repeat_winner |> count(Winner, name = "Repeats won")

# calibration across repeats
cal_summary <- per_repeat |>
  group_by(Model) |>
  summarise(Cal_slope_min = min(Cal_slope), Cal_slope_max = max(Cal_slope),
            Cal_int     = mean(Cal_int),
            Cal_slope   = mean(Cal_slope),
            Brier       = mean(Brier),
            HL_p_median = median(HL_p, na.rm = TRUE),
            HL_ok       = sum(HL_p > 0.05, na.rm = TRUE),
            NB_30_40    = mean(NB_30_40),
            .groups = "drop")

# ---- patient-level bootstrap of the repeat-mean AUROC (paired across models) ---
models_vec <- names(models)
pred_mat <- lapply(setNames(models_vec, models_vec), function(m) {   # n x repeats
  d <- dplyr::filter(pr, Model == m)
  M <- matrix(NA_real_, n_pat, N_REPEATS)
  M[cbind(d$.orig_row, d$Repeat)] <- d$.pred_Yes
  stopifnot(!anyNA(M))
  M
})
mean_auc <- function(M, idx) {
  mean(vapply(seq_len(ncol(M)), function(r) auc_fast(y_all[idx], M[idx, r]), numeric(1)))
}
point_auc <- vapply(models_vec, function(m) mean_auc(pred_mat[[m]], seq_len(n_pat)), numeric(1))

set.seed(SEED)
boot_idx <- replicate(N_BOOT, sample.int(n_pat, n_pat, replace = TRUE), simplify = FALSE)
cat("Bootstrapping repeat-mean AUROC (", N_BOOT, " resamples x ", length(models_vec), " models)...\n", sep = "")
boot_auc <- vapply(models_vec, function(m)
  vapply(boot_idx, function(i) mean_auc(pred_mat[[m]], i), numeric(1)),
  numeric(N_BOOT))                                   # N_BOOT x models

auroc_summary <- per_repeat |>
  group_by(Model) |>
  summarise(SD_across_repeats = sd(AUROC), Min = min(AUROC), Max = max(AUROC), .groups = "drop") |>
  mutate(AUROC   = point_auc[Model],
         CI_low  = apply(boot_auc[, Model, drop = FALSE], 2, quantile, 0.025, na.rm = TRUE),
         CI_high = apply(boot_auc[, Model, drop = FALSE], 2, quantile, 0.975, na.rm = TRUE)) |>
  left_join(win_counts |> dplyr::rename(Model = Winner), by = "Model") |>
  mutate(`Repeats won` = coalesce(`Repeats won`, 0L)) |>
  arrange(desc(AUROC))

paired_vs <- function(ref) purrr::map_dfr(setdiff(auroc_summary$Model, ref), function(m) {
  d <- boot_auc[, ref] - boot_auc[, m]
  tibble(Reference = ref, Model = m,
         dAUROC  = unname(point_auc[ref] - point_auc[m]),
         CI_low  = unname(quantile(d, 0.025, na.rm = TRUE)),
         CI_high = unname(quantile(d, 0.975, na.rm = TRUE)),
         p       = min(1, 2 * min(mean(d <= 0, na.rm = TRUE), mean(d >= 0, na.rm = TRUE))))
})

# ---- PRE-SPECIFIED final-model selection rule (Supplementary Section S1.5) ------
# 1 rank by mean AUROC; 2 paired bootstrap vs the leader; 3 candidates = leader +
# models with p >= 0.05; 4 keep median HL p > 0.05 (if none, keep all);
# 5 calibration slope closest to 1, tie-break higher mean net benefit at 30-40%.
auc_rank   <- auroc_summary$Model
leader     <- auc_rank[1]
vs_leader  <- paired_vs(leader)
candidates <- c(leader, vs_leader$Model[vs_leader$p >= 0.05])

cand_tab <- cal_summary |>
  dplyr::filter(Model %in% candidates) |>
  mutate(AUROC = point_auc[Model], slope_gap = abs(Cal_slope - 1),
         Adequate_calibration = !is.na(HL_p_median) & HL_p_median > 0.05)
pool <- if (any(cand_tab$Adequate_calibration)) dplyr::filter(cand_tab, Adequate_calibration) else cand_tab
best_model <- pool |> arrange(slope_gap, desc(NB_30_40)) |> slice(1) |> pull(Model)

top3_models <- unique(c(best_model, auc_rank))[1:3]   # figures: final model + next two by AUROC
vs_final    <- paired_vs(best_model)
THR         <- unname(THR_LOOKUP[best_model])

cat("Leader:", leader, "| candidates:", paste(candidates, collapse = ", "),
    "| FINAL MODEL:", best_model, "| threshold", round(THR, 3), "\n")

# ---- Output tables (paper numbering) -------------------------------------------
pm <- function(m) ifelse(m == best_model, paste0(paper_model_names[m], " (final model)"),
                         paper_model_names[m])
cls_summary <- cls_repeat |>
  group_by(Model) |>
  summarise(Accuracy = fmt_mr(Accuracy), Sensitivity = fmt_mr(Sensitivity),
            Specificity = fmt_mr(Specificity), PPV = fmt_mr(PPV), NPV = fmt_mr(NPV),
            F1 = fmt_mr(F1), .groups = "drop")

# Table 2
table2 <- auroc_summary |>
  left_join(cal_summary, by = "Model") |>
  left_join(cls_summary, by = "Model") |>
  mutate(Thr = unname(THR_LOOKUP[Model])) |>          # look up before Model is relabelled
  transmute(Model = pm(Model),
            `AUROC (95% CI)` = sprintf("%.2f (%.2f-%.2f)", AUROC, CI_low, CI_high),
            `Calibration slope` = sprintf("%.2f", Cal_slope),
            `Brier score` = sprintf("%.3f", Brier),
            Threshold = sprintf("%.2f", Thr),
            Sensitivity, Specificity, PPV, NPV)
write.csv(table2, tab_path("Table2_model_performance.csv"), row.names = FALSE)

# Table S4a
table_s4a <- auroc_summary |>
  transmute(Model = pm(Model),
            `Mean AUROC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", AUROC, CI_low, CI_high),
            `SD across repeats` = round(SD_across_repeats, 3),
            `Range across repeats` = sprintf("%.3f-%.3f", Min, Max),
            `Repeats ranked first` = `Repeats won`)
write.csv(table_s4a, tab_path("TableS4a_auroc_all_models.csv"), row.names = FALSE)

# Table S4b
table_s4b <- vs_final |>
  transmute(Comparison = paste(best_model, "vs", paper_model_names[Model]),
            dAUROC = f3s(dAUROC),
            `95% CI` = paste(f3s(CI_low), "to", f3s(CI_high)),
            `p value` = fmt_p(p))
write.csv(table_s4b, tab_path("TableS4b_paired_vs_final_model.csv"), row.names = FALSE)

# Table S4c: selection trace for all eight models
table_s4c <- auroc_summary |>
  dplyr::select(Model, AUROC) |>
  left_join(vs_leader |> dplyr::select(Model, dAUROC, CI_low, CI_high, p), by = "Model") |>
  left_join(cal_summary |> dplyr::select(Model, Cal_slope, HL_p_median, NB_30_40), by = "Model") |>
  mutate(in_candidates = Model %in% candidates,
         adequate = !is.na(HL_p_median) & HL_p_median > 0.05,
         Outcome = case_when(
           Model == best_model ~ "Selected (slope closest to 1)",
           !in_candidates      ~ "Excluded at step 3 (differs from leader)",
           !adequate & any(cand_tab$Adequate_calibration) ~ "Excluded at step 4 (median HL p <=0.05)",
           TRUE                ~ "Candidate; not selected at step 5")) |>
  transmute(Model = paper_model_names[Model],
            `Mean AUROC` = sprintf("%.3f", AUROC),
            `dAUROC vs leader (95% CI); p` = ifelse(is.na(dAUROC), "Leader",
              sprintf("%.3f (%.3f to %.3f); p=%s", dAUROC, CI_low, CI_high, fmt_p(p))),
            `Mean calibration slope` = ifelse(in_candidates, sprintf("%.2f", Cal_slope), "-"),
            `Median HL p` = ifelse(in_candidates, sprintf("%.3f", HL_p_median), "-"),
            `Mean net benefit 30-40%` = ifelse(in_candidates, sprintf("%.3f", NB_30_40), "-"),
            `Outcome of rule` = Outcome)
write.csv(table_s4c, tab_path("TableS4c_model_selection_trace.csv"), row.names = FALSE)

# Table S4d
table_s4d <- per_repeat |>
  dplyr::select(Repeat, Model, AUROC) |>
  pivot_wider(names_from = Model, values_from = AUROC) |>
  left_join(repeat_winner, by = "Repeat") |>
  arrange(Repeat) |>
  mutate(Repeat = sprintf("%d (%d)", Repeat, REPEAT_SEEDS[Repeat]),
         across(where(is.double), ~ sprintf("%.3f", .x)),
         Highest = paper_model_names[Winner]) |>
  dplyr::select(Repeat, all_of(auc_rank), Highest)
write.csv(table_s4d, tab_path("TableS4d_auroc_by_repeat.csv"), row.names = FALSE)

# Table S7: calibration of all eight models
table_s7 <- auroc_summary |>
  dplyr::select(Model) |>
  left_join(cal_summary, by = "Model") |>
  transmute(Model = pm(Model),
            `Mean calibration intercept` = sprintf("%.2f", Cal_int),
            `Mean calibration slope (range)` = sprintf("%.2f (%.2f-%.2f)", Cal_slope,
                                                       Cal_slope_min, Cal_slope_max),
            `Mean Brier score` = sprintf("%.3f", Brier),
            `Median HL p` = ifelse(HL_p_median < 0.001, "<0.001", sprintf("%.3f", HL_p_median)),
            `Repeats with HL p>0.05` = sprintf("%d/%d", HL_ok, N_REPEATS))
write.csv(table_s7, tab_path("TableS7_calibration_all_models.csv"), row.names = FALSE)

# Table S8: threshold and threshold-based performance
table_s8 <- auroc_summary |>
  dplyr::select(Model) |>
  left_join(cls_summary, by = "Model") |>
  mutate(Thr = unname(THR_LOOKUP[Model])) |>
  transmute(Model = pm(Model), Threshold = sprintf("%.2f", Thr),
            Accuracy, Sensitivity, Specificity, PPV, NPV, F1)
write.csv(table_s8, tab_path("TableS8_threshold_performance.csv"), row.names = FALSE)

# Machine-readable per-repeat results (aggregate, no participant rows)
write.csv(per_repeat, tab_path("per_repeat_measures.csv"), row.names = FALSE)
write.csv(cls_repeat, tab_path("per_repeat_classification.csv"), row.names = FALSE)
write.csv(auroc_summary, tab_path("auroc_summary_unrounded.csv"), row.names = FALSE)
