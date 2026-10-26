# =============================================================================
# 03_feature_selection.R — ensemble feature selection (run inside each outer fold)
# -----------------------------------------------------------------------------
# Union of LASSO (alpha = 1) and Elastic Net (alpha = 0.5) at lambda.1se and the
# top quartile by mutual information, mapped back to source variables, then
# iterative VIF > 5 removal among numeric predictors.
# Paper: Methods (Model development); Supplementary Section S1.2
# =============================================================================

perform_feature_selection <- function(data, target = "azithromycin_failure", seed = 123) {
  predictors <- setdiff(names(data), target)

  # dummy coding + near-zero-variance filter + normalisation for the penalised methods
  fs_data <- recipe(as.formula(paste(target, "~ .")), data = data) |>
    step_zv(all_predictors()) |>
    step_dummy(all_nominal_predictors()) |>
    step_nzv(all_predictors()) |>
    step_normalize(all_numeric_predictors()) |>
    prep() |> bake(new_data = NULL)
  X <- as.matrix(fs_data[, -which(names(fs_data) == target)])
  y <- as.numeric(fs_data[[target]] == "Yes")

  # 1. LASSO
  set.seed(seed)
  lasso_feats <- tryCatch({
    cv <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5)
    rownames(coef(cv, s = "lambda.1se"))[which(coef(cv, s = "lambda.1se") != 0)][-1]
  }, error = \(e) character(0))

  # 2. Elastic Net
  set.seed(seed)
  enet_feats <- tryCatch({
    cv <- cv.glmnet(X, y, family = "binomial", alpha = 0.5, nfolds = 5)
    rownames(coef(cv, s = "lambda.1se"))[which(coef(cv, s = "lambda.1se") != 0)][-1]
  }, error = \(e) character(0))

  # 3. Mutual information (top 25%)
  mi_feats <- tryCatch({
    scores <- sapply(colnames(X), \(f)
      infotheo::mutinformation(infotheo::discretize(X[, f], nbins = 5), y))
    names(scores)[scores >= quantile(scores, 0.75)]
  }, error = \(e) character(0))

  # Map dummy columns back to the original variables
  map_to_original <- function(dummy, orig_vars) {
    match <- orig_vars[sapply(orig_vars, \(v) startsWith(dummy, v))]
    if (length(match) > 0) match[which.max(nchar(match))] else NA
  }
  selected_dummies <- unique(c(lasso_feats, enet_feats, mi_feats))
  selected <- unique(na.omit(sapply(selected_dummies, map_to_original, predictors)))

  # VIF filtering (threshold 5), numeric predictors only. VIF on dummy-coded
  # categorical variables is unreliable, so categorical variables are kept.
  if (length(selected) > 1) {
    numeric_sel <- selected[sapply(data[, selected, drop = FALSE], is.numeric)]
    if (length(numeric_sel) > 1) {
      vif_df <- as.data.frame(scale(data[, numeric_sel, drop = FALSE]))
      calc_vif <- function(d) sapply(names(d), \(v) tryCatch(
        1 / (1 - summary(lm(as.formula(paste(v, "~ .")), d))$r.squared),
        error = \(e) NA))
      while (ncol(vif_df) > 1) {
        vifs <- calc_vif(vif_df)
        if (all(is.na(vifs)) || max(vifs, na.rm = TRUE) <= 5) break
        vif_df <- vif_df[, -which.max(vifs), drop = FALSE]
      }
      selected <- unique(c(colnames(vif_df), setdiff(selected, numeric_sel)))
    }
  }

  if (length(selected) == 0) selected <- predictors
  selected
}
