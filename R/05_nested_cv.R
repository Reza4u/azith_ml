# =============================================================================
# 05_nested_cv.R — repeated nested cross-validation and the final predictor set
# -----------------------------------------------------------------------------
# 10 repeats x 5 outer folds (performance) x 10 inner folds (tuning).
# Feature selection and all pre-processing are estimated inside each outer
# training fold only. All eight algorithms share the outer folds of a repeat.
# Completed repeats are cached in output/private/repeats/ and skipped on rerun.
# Paper: Methods; Supplementary Sections S1.4 and S1.7; Supplementary Table S2
# =============================================================================

create_recipe <- function(data) {
  recipe(azithromycin_failure ~ ., data = data) |>
    step_zv(all_predictors()) |>
    step_novel(all_nominal_predictors()) |>
    step_dummy(all_nominal_predictors()) |>
    step_nzv(all_predictors()) |>
    step_normalize(all_numeric_predictors())
}

# Tune one model in the inner folds, refit on the outer training fold and
# predict the held-out outer fold
fit_predict_fold <- function(train_data, test_data, inner_folds, model_spec, model_grid, fold_seed) {
  wf <- workflow() |> add_recipe(create_recipe(train_data)) |> add_model(model_spec)
  set.seed(fold_seed)
  tune_res <- tryCatch(
    tune_grid(wf, inner_folds, grid = model_grid, metrics = metric_set(roc_auc),
              control = control_grid(verbose = FALSE, save_pred = FALSE)),
    error = function(e) { message("    tuning error: ", conditionMessage(e)); NULL })
  if (is.null(tune_res)) return(NULL)
  best <- tryCatch(select_best(tune_res, metric = "roc_auc"), error = function(e) NULL)
  if (is.null(best)) return(NULL)
  set.seed(fold_seed)
  fit_wf <- fit(finalize_workflow(wf, best), train_data)
  list(prob = predict(fit_wf, test_data, type = "prob")$.pred_Yes, params = best)
}

# One complete repetition of the nested procedure
run_one_repeat <- function(r, data = ml_df, use_cache = TRUE) {
  cache <- file.path(REPEAT_DIR, sprintf("repeat_%02d.rds", r))
  if (use_cache && file.exists(cache)) {
    cat(sprintf("  repeat %2d: loaded from cache\n", r))
    return(readRDS(cache))
  }
  seed_r <- REPEAT_SEEDS[r]
  set.seed(seed_r)
  outer_folds <- vfold_cv(data, v = N_OUTER, strata = azithromycin_failure)

  preds <- list(); feats <- vector("list", N_OUTER); params <- list()
  t0 <- Sys.time()
  for (i in seq_len(N_OUTER)) {
    split      <- outer_folds$splits[[i]]
    train_full <- training(split)
    test_full  <- testing(split)
    orig_rows  <- rsample::complement(split)   # row of each held-out child in ml_df
    fold_seed  <- seed_r + i

    # feature selection inside the outer training fold only
    selected   <- perform_feature_selection(train_full, seed = fold_seed)
    feats[[i]] <- selected
    train_data <- train_full[, c("azithromycin_failure", selected)]
    test_data  <- test_full[,  c("azithromycin_failure", selected)]

    set.seed(fold_seed)
    inner_folds <- vfold_cv(train_data, v = N_INNER, strata = azithromycin_failure)

    for (m in names(models)) {
      out <- fit_predict_fold(train_data, test_data, inner_folds,
                              models[[m]], model_grids[[m]], fold_seed)
      if (is.null(out)) stop("Model '", m, "' failed in repeat ", r, ", fold ", i, ".")
      preds[[length(preds) + 1]] <- tibble(
        Repeat = r, Fold = i, Model = m, .orig_row = orig_rows,
        azithromycin_failure = test_data$azithromycin_failure,
        .pred_Yes = out$prob)
      params[[length(params) + 1]] <- out$params |>
        mutate(across(everything(), as.character)) |>
        mutate(Repeat = r, Fold = i, Model = m)
    }
    cat(sprintf("  repeat %2d  fold %d/%d  %2d predictors  (%.1f min)\n", r, i, N_OUTER,
                length(selected), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res <- list(predictions = bind_rows(preds), features = feats,
              params = bind_rows(params), seed = seed_r)
  saveRDS(res, cache)
  res
}

# ---- Reproducibility check: repeat 1 run twice must be identical ------------
if (RUN_REPRO_CHECK && !file.exists(file.path(REPEAT_DIR, "repeat_01.rds"))) {
  cat("Reproducibility check: running repeat 1 twice...\n")
  run_a <- run_one_repeat(1, use_cache = FALSE)
  run_b <- run_one_repeat(1, use_cache = FALSE)
  if (!isTRUE(all.equal(run_a$predictions, run_b$predictions, tolerance = 0))) {
    stop("Reproducibility check FAILED: repeat 1 differs between two runs.")
  }
  cat("Reproducibility check passed: repeat 1 identical on two runs\n")
  rm(run_a, run_b)
}

# ---- All repeats --------------------------------------------------------------
cat("Repeated nested CV:", N_REPEATS, "repeats x", N_OUTER, "outer x", N_INNER, "inner folds\n")
repeat_results <- lapply(seq_len(N_REPEATS), run_one_repeat)

cv_predictions <- bind_rows(lapply(repeat_results, `[[`, "predictions"))
fold_features  <- unlist(lapply(repeat_results, `[[`, "features"), recursive = FALSE)
tuned_params   <- bind_rows(lapply(repeat_results, `[[`, "params"))

# Completeness: every child has exactly one prediction per model per repeat
completeness <- cv_predictions |> count(Model, Repeat)
stopifnot(length(fold_features) == N_REPEATS * N_OUTER,
          nrow(completeness) == length(models) * N_REPEATS,
          all(completeness$n == n_pat),
          nrow(cv_predictions |> count(Model, Repeat, .orig_row) |> filter(n > 1)) == 0,
          all(cv_predictions$azithromycin_failure ==
                ml_df$azithromycin_failure[cv_predictions$.orig_row]))

saveRDS(cv_predictions, priv_path("cv_predictions_all_repeats.rds"))   # participant-level
write.csv(tuned_params, tab_path("tuned_hyperparameters_by_fold.csv"), row.names = FALSE)

# ---- Selection stability and the FINAL PREDICTOR SET (Supplementary Table S2) --
n_selections <- length(fold_features)
sel_counts   <- table(factor(unlist(lapply(fold_features, unique)), levels = candidate_vars))

stability <- tibble(var = names(sel_counts), n_sel = as.integer(sel_counts)) |>
  mutate(Variable = vapply(var, lab, character(1)),   # internal labels fix the order
         freq     = n_sel / n_selections,
         in_final = freq >= SEL_FREQ_MIN) |>
  arrange(desc(n_sel), Variable)

FINAL_SET <- stability$var[stability$in_final]       # order matters for the final fit
if (length(FINAL_SET) == 0) stop("No variable reached the selection-frequency threshold.")

table_s2 <- stability |>
  transmute(Variable = plab(var),
            `Outer-fold selections` = sprintf("%d/%d", n_sel, n_selections),
            `Selection frequency (%)` = round(100 * freq),
            `Final predictor set (>=50%)` = ifelse(in_final, "Included", "Excluded"))
write.csv(table_s2, tab_path("TableS2_feature_selection_stability.csv"), row.names = FALSE)

fold_sizes <- lengths(fold_features)
cat("Final predictor set (", length(FINAL_SET), "): ", paste(plab(FINAL_SET), collapse = ", "),
    "\nEvents per predictor: ", round(n_events / length(FINAL_SET), 1), "\n", sep = "")
