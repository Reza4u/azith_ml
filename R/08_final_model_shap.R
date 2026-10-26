# =============================================================================
# 08_final_model_shap.R — final model refit and SHAP (Supplementary Table S10)
# -----------------------------------------------------------------------------
# The model chosen by the pre-specified rule is refitted on all 329 children
# with FINAL_SET; hyperparameters are tuned by stratified 5-fold CV over the
# same fixed grid. SHAP: kernelSHAP (Table S10) and exact Shapley values over
# all 2^10 coalitions (Figure 3). Both use all 329 children as background.
# Paper: Methods (Final model); Supplementary Sections S1.6 and S1.8
# =============================================================================

train_data_final <- ml_df[, c("azithromycin_failure", FINAL_SET)]
final_wf <- workflow() |> add_recipe(create_recipe(train_data_final)) |> add_model(models[[best_model]])

set.seed(SEED)
cv_folds <- vfold_cv(train_data_final, v = 5, strata = azithromycin_failure)
set.seed(SEED)
tune_res_global <- tune_grid(final_wf, resamples = cv_folds, grid = model_grids[[best_model]],
                             metrics = metric_set(roc_auc), control = control_grid(verbose = FALSE))
best_params <- select_best(tune_res_global, metric = "roc_auc")
final_wf    <- finalize_workflow(final_wf, best_params)

set.seed(SEED)
final_fit <- fit(final_wf, train_data_final)
saveRDS(final_fit, priv_path("final_model_fit.rds"))   # contains training data: do not share
write.csv(best_params, tab_path("final_model_hyperparameters.csv"), row.names = FALSE)
cat("Final model:", best_model, "with", length(FINAL_SET), "predictors;",
    paste(names(best_params)[names(best_params) != ".config"],
          best_params[names(best_params) != ".config"], sep = "=", collapse = ", "), "\n")

# Number of model columns after dummy coding and near-zero-variance filtering
n_model_columns <- ncol(extract_mold(final_fit)$predictors)

# ---- kernelSHAP (Supplementary Table S10) --------------------------------------
X_shap   <- train_data_final |> dplyr::select(-azithromycin_failure) |> as.data.frame()
pred_fun <- function(model, newdata) predict(model, new_data = newdata, type = "prob")$.pred_Yes

set.seed(SEED)
cat("kernelSHAP for", nrow(X_shap), "children...\n")
ks <- kernelshap(final_fit, X = X_shap, bg_X = X_shap, pred_fun = pred_fun, verbose = FALSE)
S_kernel <- ks$S                                     # n x 10, columns = FINAL_SET

shap_direction <- function(S, X) purrr::map_dfr(colnames(S), function(v) {
  xj <- X[[v]]
  multi_cat <- !is.numeric(xj) && length(unique(xj)) > 2
  val <- if (is.numeric(xj)) xj else as.numeric(as.factor(xj))
  r <- if (multi_cat) NA_real_ else suppressWarnings(cor(val, as.numeric(S[, v])))
  sj <- S[, v]
  tibble(var = v, Predictor = plab(v),
         `Mean absolute SHAP` = mean(abs(sj)),
         `Mean positive SHAP` = if (any(sj > 0)) mean(sj[sj > 0]) else 0,
         `Mean negative SHAP` = if (any(sj < 0)) mean(sj[sj < 0]) else 0,
         `Correlation of value with SHAP` = r,
         Direction = case_when(multi_cat ~ "categorical: see Figure 3",
                               r >  0.15 ~ "higher values towards failure",
                               r < -0.15 ~ "higher values towards success",
                               TRUE      ~ "weak or non-monotone"))
}) |> arrange(desc(round(`Mean absolute SHAP`, 3)))   # ties keep predictor order, as in the paper

shap_kernel_summary <- shap_direction(S_kernel, X_shap)
table_s10 <- shap_kernel_summary |>
  mutate(across(c(`Mean absolute SHAP`, `Mean positive SHAP`, `Mean negative SHAP`), ~ round(.x, 3)),
         `Correlation of value with SHAP` = round(`Correlation of value with SHAP`, 2)) |>
  dplyr::select(-var)
write.csv(table_s10, tab_path("TableS10_SHAP_summary.csv"), row.names = FALSE)
saveRDS(list(S = S_kernel, X = X_shap), priv_path("shap_kernel.rds"))   # participant-level

# ---- Exact Shapley values (Figure 3) -------------------------------------------
# All 2^p coalitions; the value of a coalition is the mean prediction over the
# 329 background children with the coalition's features fixed. Takes ~15 min.
exact_shapley <- function(fit, X) {
  rec <- extract_recipe(fit)
  eng <- extract_fit_engine(fit)
  f <- function(d) {
    d$azithromycin_failure <- factor("No", levels = c("Yes", "No"))
    predict(eng, data = bake(rec, new_data = d), num.threads = 2)$predictions[, "Yes"]
  }
  key <- function(d) if (ncol(d) == 0) rep("", nrow(d)) else
    do.call(paste, c(lapply(d, as.character), sep = "\r"))
  p <- ncol(X); n <- nrow(X); coal <- 0:(2^p - 1)
  V <- matrix(NA_real_, length(coal), n)
  for (ci in seq_along(coal)) {
    S  <- which(bitwAnd(coal[ci], 2^(0:(p - 1))) > 0)
    kx <- key(X[, S, drop = FALSE]); ux <- !duplicated(kx); xi <- match(kx, kx[ux])
    kb <- key(X[, setdiff(seq_len(p), S), drop = FALSE])
    ub <- !duplicated(kb); bj <- match(kb, kb[ub]); wb <- tabulate(bj, sum(ub)) / n
    XU <- X[ux, , drop = FALSE]; BU <- X[ub, , drop = FALSE]
    nu <- nrow(XU); nb <- nrow(BU)
    D <- BU[rep(seq_len(nb), times = nu), , drop = FALSE]
    if (length(S)) D[, S] <- XU[rep(seq_len(nu), each = nb), S, drop = FALSE]
    V[ci, ] <- colSums(matrix(f(D), nrow = nb) * wb)[xi]
    if (ci %% 128 == 0) cat("  ", ci, "of", length(coal), "coalitions\n")
  }
  phi <- matrix(0, n, p, dimnames = list(NULL, names(X)))
  sz  <- sapply(coal, function(c) sum(bitwAnd(c, 2^(0:(p - 1))) > 0))
  for (k in 1:p) {
    bit <- 2^(k - 1)
    for (c in coal[bitwAnd(coal, bit) == 0]) {
      s <- sz[c + 1]; w <- factorial(s) * factorial(p - s - 1) / factorial(p)
      phi[, k] <- phi[, k] + w * (V[c + bit + 1, ] - V[c + 1, ])
    }
  }
  list(S = phi, baseline = V[1, 1], pred = V[2^p, ],
       additivity_error = max(abs(rowSums(phi) + V[1, 1] - V[2^p, ])))
}

RUN_EXACT_SHAP <- as.logical(Sys.getenv("AZITH_EXACT_SHAP", "TRUE"))
if (RUN_EXACT_SHAP) {
  X_exact <- as.data.frame(extract_recipe(final_fit, estimated = FALSE)$template)
  X_exact$azithromycin_failure <- NULL
  shap_exact <- exact_shapley(final_fit, X_exact)
  stopifnot(shap_exact$additivity_error < 1e-8)
  shap_exact_summary <- shap_direction(shap_exact$S, X_exact)
  saveRDS(c(shap_exact, list(X = X_exact)), priv_path("shap_exact.rds"))      # participant-level
  write.csv(as.data.frame(shap_exact$S), priv_path("shap_exact_S.csv"), row.names = FALSE)
  write.csv(X_exact, priv_path("shap_exact_X.csv"), row.names = FALSE)
  write.csv(shap_exact_summary |> dplyr::select(Predictor, `Mean absolute SHAP`) |>
              mutate(`Mean absolute SHAP` = round(`Mean absolute SHAP`, 4)),
            tab_path("SHAP_exact_mean_absolute.csv"), row.names = FALSE)
}
