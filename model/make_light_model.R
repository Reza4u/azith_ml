# =============================================================================
# make_light_model.R — shareable copy of the final model with all training data removed
# -----------------------------------------------------------------------------
# model/final_model_fit_light.rds was made this way. Run from the repository
# root after run_all.R:   Rscript model/make_light_model.R
# Reads  output/private/final_model_fit.rds (contains training data; never share)
# Writes output/final_model_fit_light.rds only if (1) no object with one row per
# child is left anywhere in it and (2) its predictions equal the full model's.
# =============================================================================
suppressMessages({library(tidymodels); library(butcher)})
OUT_DIR <- Sys.getenv("AZITH_OUT_DIR", "output")

fit <- readRDS(file.path(OUT_DIR, "private", "final_model_fit.rds"))
stopifnot(inherits(fit, "workflow"))
train_rows <- fit$pre$actions$recipe$recipe$template      # kept in memory for the check only
N <- nrow(train_rows)

light <- butcher::butcher(fit)
# 1. training data kept by workflows (the "mold")
light$pre$mold$predictors <- light$pre$mold$predictors[0, , drop = FALSE]
light$pre$mold$outcomes   <- light$pre$mold$outcomes[0, , drop = FALSE]
light$pre$mold$extras     <- NULL
# 2. data kept inside the recipes (unprepped and prepped)
if (!is.null(light$pre$actions$recipe$recipe$template))
  light$pre$actions$recipe$recipe$template <- light$pre$actions$recipe$recipe$template[0, , drop = FALSE]
bp <- light$pre$mold$blueprint
if (!is.null(bp$recipe$template)) bp$recipe$template <- bp$recipe$template[0, , drop = FALSE]
light$pre$mold$blueprint <- bp
# 3. copies of the training data in the environments of the recipe step formulas
clean_env <- function(steps) for (s in steps) for (q in s$terms) {
  e <- attr(q, ".Environment")
  if (is.environment(e) && !identical(e, globalenv()))
    for (nm in ls(e, all.names = TRUE)) {
      v <- get(nm, envir = e); if (is.data.frame(v)) assign(nm, v[0, ], envir = e)
    }
}
clean_env(light$pre$actions$recipe$recipe$steps)
clean_env(light$pre$mold$blueprint$recipe$steps)
# 4. per-child out-of-bag predictions stored by ranger
eng <- extract_fit_engine(light)
if (!is.null(eng$predictions)) { eng$predictions <- NULL; light$fit$fit$fit <- eng }

# ---- check 1: nothing with one row per child is left ----
`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a
hits <- character(0)
walk <- function(x, path = "fit", depth = 0) {
  if (depth > 12) return(invisible())
  if (is.data.frame(x) || is.matrix(x)) { if (NROW(x) == N) hits <<- c(hits, path); return(invisible()) }
  if (is.atomic(x) && length(x) == N) { hits <<- c(hits, path); return(invisible()) }
  if (is.environment(x)) {
    if (!identical(x, globalenv()) && !isNamespace(x))
      for (nm in ls(x, all.names = TRUE)) walk(get(nm, envir = x), paste0(path, "<env>$", nm), depth + 1)
    return(invisible())
  }
  if (is.list(x)) for (nm in seq_along(x)) walk(x[[nm]], paste0(path, "$", names(x)[nm] %||% nm), depth + 1)
}
walk(light)
if (length(hits)) stop("Training-size objects remain; not saved:\n", paste(hits, collapse = "\n"))

# ---- check 2: predictions unchanged ----
stopifnot(isTRUE(all.equal(predict(fit,   train_rows, type = "prob")$.pred_Yes,
                           predict(light, train_rows, type = "prob")$.pred_Yes)))
rm(train_rows)
saveRDS(light, file.path(OUT_DIR, "final_model_fit_light.rds"))
cat("Light model written to", file.path(OUT_DIR, "final_model_fit_light.rds"), "\n")
