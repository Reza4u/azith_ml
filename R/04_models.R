# =============================================================================
# 04_models.R — eight algorithms and their fixed tuning grids
# -----------------------------------------------------------------------------
# ranger and xgboost run single-threaded so that results do not depend on the
# number of processor cores. Grids are fixed (no random candidates).
# Paper: Supplementary Section S1.3
# =============================================================================

models <- list(
  "Elastic Net"    = logistic_reg(penalty = tune(), mixture = tune()) |>
                       set_engine("glmnet") |> set_mode("classification"),
  "Random Forest"  = rand_forest(mtry = tune(), trees = 500, min_n = tune()) |>
                       set_engine("ranger", importance = "impurity", num.threads = 1) |>
                       set_mode("classification"),
  "XGBoost"        = boost_tree(trees = tune(), tree_depth = tune(), learn_rate = tune()) |>
                       set_engine("xgboost", nthread = 1) |> set_mode("classification"),
  "KNN"            = nearest_neighbor(neighbors = tune(), weight_func = tune()) |>
                       set_engine("kknn") |> set_mode("classification"),
  "SVM"            = svm_rbf(cost = tune(), rbf_sigma = tune()) |>
                       set_engine("kernlab") |> set_mode("classification"),
  "Decision Tree"  = decision_tree(cost_complexity = tune(), tree_depth = tune(), min_n = tune()) |>
                       set_engine("rpart") |> set_mode("classification"),
  "Naive Bayes"    = naive_Bayes(smoothness = tune(), Laplace = tune()) |>
                       set_engine("klaR") |> set_mode("classification"),
  "Neural Network" = mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) |>
                       set_engine("nnet") |> set_mode("classification")
)

model_grids <- list(
  "Elastic Net"    = tidyr::crossing(penalty = 10^c(-4, -3, -2, -1.5, -1, -0.5),
                                     mixture = c(0.2, 0.5, 0.8, 1)),
  "Random Forest"  = tidyr::crossing(mtry = c(2L, 3L, 4L, 6L),
                                     min_n = c(2L, 5L, 10L, 20L)),
  "XGBoost"        = tidyr::crossing(trees = c(100L, 300L),
                                     tree_depth = c(2L, 4L, 6L),
                                     learn_rate = c(0.01, 0.05, 0.1)),
  "KNN"            = tidyr::crossing(neighbors = c(5L, 10L, 15L, 20L, 30L, 40L),
                                     weight_func = c("rectangular", "triangular", "optimal")),
  "SVM"            = tidyr::crossing(cost = 2^c(-2, 0, 2, 4),
                                     rbf_sigma = 10^c(-3, -2, -1)),
  "Decision Tree"  = tidyr::crossing(cost_complexity = 10^c(-4, -3, -2),
                                     tree_depth = c(3L, 5L, 8L),
                                     min_n = c(10L, 20L)),
  "Naive Bayes"    = tidyr::crossing(smoothness = c(0.5, 1, 1.5),
                                     Laplace = c(0, 1)),
  "Neural Network" = tidyr::crossing(hidden_units = c(2L, 5L, 10L),
                                     penalty = 10^c(-3, -2, -1),
                                     epochs = c(100L, 300L))
)

# Every grid column must be a tuning parameter of its model, and vice versa
stopifnot(setequal(names(model_grids), names(models)))
for (m in names(models)) {
  stopifnot(setequal(extract_parameter_set_dials(models[[m]])$id, names(model_grids[[m]])))
}
cat("Models:", paste(names(models), collapse = ", "), "| grid sizes:",
    paste(sapply(model_grids, nrow), collapse = ", "), "\n")
