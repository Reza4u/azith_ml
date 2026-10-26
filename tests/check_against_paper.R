# =============================================================================
# check_against_paper.R — compare every reported number with the code's output
# -----------------------------------------------------------------------------
# Expected values are copied from the published paper (main text, Tables 1-2 and
# Supplementary Tables S1-S10). Observed values are read from output/tables/.
# Writes output/verification_report.csv and .md. Run after run_all.R:
#   Rscript tests/check_against_paper.R
# =============================================================================
suppressMessages(library(dplyr))
OUT_DIR <- Sys.getenv("AZITH_OUT_DIR", "output")
tb <- function(f) read.csv(file.path(OUT_DIR, "tables", f), check.names = FALSE,
                           stringsAsFactors = FALSE, colClasses = "character")
norm <- function(x) {
  x <- gsub("−", "-", as.character(x))          # minus sign
  x <- gsub("–", "-", x)                         # en dash
  x <- gsub("%", "", x, fixed = TRUE)
  x <- gsub("Na\xc3\xafve", "Naive", x, useBytes = TRUE)   # "Naive" with i-diaeresis (UTF-8 bytes)
  x <- gsub("Na<c3><af>ve", "Naive", x, fixed = TRUE)      # same, as printed in a non-UTF-8 locale
  x <- gsub("\\s+", " ", trimws(x))
  x
}
checks <- list()
add <- function(section, item, expected, observed) {
  n <- max(length(expected), length(observed))
  checks[[length(checks) + 1]] <<- tibble(section = section, item = rep(item, length.out = n),
                                         expected = rep(as.character(expected), length.out = n),
                                         observed = if (length(observed) == 0) NA_character_ else
                                           rep(as.character(observed), length.out = n))
}
num_eq <- function(x, d) sprintf(paste0("%.", d, "f"), as.numeric(x))

# ---- key numbers quoted in the text --------------------------------------------
kn <- tb("key_numbers.csv"); k <- setNames(kn$value, kn$name)
expect_key <- c(
  n_children = "329", n_failures = "132", pct_failures = "40.1", n_success = "197",
  pct_success = "59.9", n_candidate_variables = "23", n_missing_values = "0",
  resp_rate_values_replaced = "1 and 229", resp_rate_replacement = "31",
  columns_after_dummy = "30", columns_after_nzv = "26",
  nzv_removed = "abx_cat_Other, blood_RE_Moderate, blood_RE_Plenty, is_brown",
  cipro_containing_failures = "62/132 (47.0%)", cipro_containing_successes = "54/197 (27.4%)",
  failures_one_component = "110 (83.3%)", n_outer_selections = "50", n_final_predictors = "10",
  final_predictors = paste(sort(c("diarrhea_dur", "macro_numeric", "resident_area", "mucus_RE",
                                  "abx_cat", "wbc_numeric", "blood_RE", "straining", "age",
                                  "pcell_numeric")), collapse = ", "),
  events_per_predictor = "13.2", model_columns_final = "13", events_per_model_column = "10.2",
  events_per_candidate_column = "4.4", riley_r2cs = "0.13",
  leader = "Elastic Net",
  candidates = "Elastic Net, Naive Bayes, Random Forest, SVM, XGBoost",
  final_model = "Random Forest", final_model_repeats_won = "1",
  final_auroc = "0.71", final_auroc_ci = "0.66-0.76", final_auroc_range = "0.68-0.73",
  leader_minus_final_dauroc = "0.009", leader_minus_final_ci = "-0.020 to 0.036",
  leader_minus_final_p = "0.53",
  final_cal_intercept = "-0.01", final_cal_slope = "1.05", final_brier = "0.208",
  final_hl_p_median = "0.28", final_threshold = "0.43",
  final_sensitivity = "0.61 (0.55-0.64)", final_specificity = "0.69 (0.67-0.72)",
  final_ppv = "0.57", final_npv = "0.72",
  nb_final_40 = "0.121", nb_treat_all_40 = "0.002", incremental_nb_40 = "0.119",
  dca_beats_both_range = "7%-87%", dca_below_treat_all = "4%, 5%, 6%",
  rf_min_n_20_outer_folds = "33 of 50", final_hyperparameters = "mtry=6, min_n=20",
  max_auroc_spread_one_model = "0.06", brier_noninformative = "0.240",
  event_rate_urban = "34.8%", event_rate_rural = "61.5%", overall_stacked_cal_slope = "1.02",
  shap_exact_vs_kernel_max_diff = "0.0000", shap_exact_blood_mucus = "0.024 and 0.024")
for (nm in names(expect_key)) add("Text", nm, expect_key[[nm]], k[nm])
# Riley et al. criteria: the paper says "about 840 children (about 340 events)"
add("Text", "riley_sample_size (about 840)", "TRUE",
    if (is.na(k["riley_sample_size"]) || k["riley_sample_size"] == "NA") "not computed" else
      as.character(abs(as.numeric(k["riley_sample_size"]) - 840) <= 10))
add("Text", "riley_events (about 340)", "TRUE",
    if (is.na(k["riley_events"]) || k["riley_events"] == "NA") "not computed" else
      as.character(abs(as.numeric(k["riley_events"]) - 340) <= 10))

# ---- Table 1 ----------------------------------------------------------------------
t1 <- tb("Table1_baseline.csv")
e1 <- read.csv(text = "
variable,level,overall,failure,success,p
age,,38.4 (10.2),38.7 (10.7),38.1 (9.8),0.783
sex,,176 (53.5),68 (51.5),108 (54.8),0.556
resident_area,,264 (80.2),92 (69.7),172 (87.3),<0.001
diarrhea_dur,,2.8 (2.0),3.0 (1.9),2.6 (2.0),0.001
fever,,176 (53.5),69 (52.3),107 (54.3),0.716
cough,,60 (18.2),21 (15.9),39 (19.8),0.371
resp_rate,,30.6 (3.2),30.4 (2.9),30.7 (3.4),0.587
straining,,127 (38.6),61 (46.2),66 (33.5),0.020
haz,,-0.9 (1.1),-1.0 (0.9),-0.8 (1.1),0.145
waz,,-1.4 (1.0),-1.5 (0.9),-1.3 (1.1),0.283
whz,,-1.3 (0.9),-1.3 (0.9),-1.3 (1.0),0.809
abx_cat,,,,,0.005
abx_cat,None,112 (34.0),38 (28.8),74 (37.6),
abx_cat,Ciprofloxacin,81 (24.6),40 (30.3),41 (20.8),
abx_cat,Metronidazole,85 (25.8),27 (20.5),58 (29.4),
abx_cat,Ciprofloxacin + metronidazole,35 (10.6),22 (16.7),13 (6.6),
abx_cat,Other,16 (4.9),5 (3.8),11 (5.6),
consist_RE,,,,,0.033
consist_RE,Formed,101 (30.7),30 (22.7),71 (36.0),
consist_RE,Mucoid,130 (39.5),60 (45.5),70 (35.5),
consist_RE,Liquid,98 (29.8),42 (31.8),56 (28.4),
blood_RE,,,,,<0.001
blood_RE,Absent,72 (21.9),27 (20.5),45 (22.8),
blood_RE,Small,237 (72.0),85 (64.4),152 (77.2),
blood_RE,Moderate,12 (3.6),12 (9.1),0 (0.0),
blood_RE,Plenty,8 (2.4),8 (6.1),0 (0.0),
mucus_RE,,,,,<0.001
mucus_RE,Small,22 (6.7),20 (15.2),2 (1.0),
mucus_RE,Moderate,178 (54.1),58 (43.9),120 (60.9),
mucus_RE,Plenty,129 (39.2),54 (40.9),75 (38.1),
wbc_numeric,,12.0 (5.4),11.7 (4.6),12.2 (5.9),0.856
macro_grade,,199 (60.5),56 (42.4),143 (72.6),<0.001
pcell_numeric,,31.0 (8.1),31.6 (7.1),30.6 (8.7),0.379
rbc_numeric,,28.8 (10.1),29.5 (9.5),28.3 (10.5),0.261
is_green,,173 (52.6),76 (57.6),97 (49.2),0.138
is_red,,118 (35.9),49 (37.1),69 (35.0),0.698
is_yellow,,113 (34.3),37 (28.0),76 (38.6),0.048
is_brown,,6 (1.8),3 (2.3),3 (1.5),0.687", colClasses = "character", na.strings = "")
for (i in seq_len(nrow(e1))) {
  r <- e1[i, ]
  o <- t1[t1$variable == r$variable & (if (is.na(r$level)) t1$row_type == "label"
                                       else t1$Variable == r$level), ]
  id <- paste(r$variable, ifelse(is.na(r$level), "", r$level))
  for (col in c("overall", "failure", "success", "p")) {
    ev <- r[[col]]; if (is.na(ev)) next
    oc <- c(overall = "Overall (n=329)", failure = "Failure", success = "Success", p = "p value")[[col]]
    add("Table 1", paste(id, col), ev, if (nrow(o) == 1) o[[oc]] else NA)
  }
}

# ---- Table 2 and Supplementary Table S8 ---------------------------------------------
e2 <- read.csv(text = "
model,auroc,slope,brier,thr,sens,spec,ppv,npv,acc,f1
Elastic Net,0.72 (0.67-0.77),0.78,0.206,0.41,0.62 (0.58-0.65),0.68 (0.65-0.71),0.57 (0.54-0.59),0.73 (0.70-0.75),0.66 (0.63-0.68),0.59 (0.55-0.62)
XGBoost,0.72 (0.67-0.76),0.63,0.211,0.40,0.64 (0.60-0.67),0.67 (0.62-0.73),0.56 (0.53-0.60),0.73 (0.71-0.76),0.65 (0.63-0.68),0.60 (0.56-0.63)
Support vector machine,0.71 (0.66-0.76),0.87,0.208,0.41,0.62 (0.58-0.67),0.70 (0.61-0.75),0.58 (0.53-0.62),0.73 (0.72-0.74),0.67 (0.63-0.69),0.60 (0.58-0.61)
Random Forest (final model),0.71 (0.66-0.76),1.05,0.208,0.43,0.61 (0.55-0.64),0.69 (0.67-0.72),0.57 (0.54-0.59),0.72 (0.69-0.74),0.66 (0.63-0.68),0.59 (0.54-0.61)
Neural network,0.70 (0.65-0.74),1.23,0.215,0.43,0.64 (0.54-0.71),0.64 (0.55-0.76),0.55 (0.50-0.60),0.73 (0.69-0.75),0.64 (0.60-0.67),0.59 (0.55-0.61)
Naive Bayes,0.69 (0.64-0.74),0.40,0.278,0.10,0.62 (0.58-0.69),0.65 (0.59-0.69),0.54 (0.51-0.58),0.72 (0.69-0.74),0.64 (0.61-0.67),0.58 (0.54-0.61)
k-nearest neighbours,0.69 (0.64-0.73),0.82,0.216,0.40,0.56 (0.48-0.61),0.71 (0.65-0.78),0.56 (0.52-0.61),0.71 (0.68-0.73),0.65 (0.61-0.68),0.56 (0.53-0.60)
Decision tree,0.68 (0.63-0.72),0.22,0.238,0.51,0.50 (0.44-0.55),0.75 (0.68-0.79),0.57 (0.53-0.62),0.69 (0.67-0.72),0.65 (0.62-0.68),0.53 (0.50-0.58)",
  colClasses = "character", check.names = FALSE, encoding = "UTF-8")
t2 <- tb("Table2_model_performance.csv"); t8 <- tb("TableS8_threshold_performance.csv")
add("Table 2", "model order (by mean AUROC)", paste(e2$model, collapse = " > "),
    paste(t2$Model, collapse = " > "))
for (i in seq_len(nrow(e2))) {
  m <- e2$model[i]; o2 <- t2[norm(t2$Model) == m, ]; o8 <- t8[norm(t8$Model) == m, ]
  g <- function(o, col) if (nrow(o) == 1) o[[col]] else NA
  add("Table 2", paste(m, "AUROC (95% CI)"), e2$auroc[i], g(o2, "AUROC (95% CI)"))
  add("Table 2", paste(m, "calibration slope"), e2$slope[i], g(o2, "Calibration slope"))
  add("Table 2", paste(m, "Brier"), e2$brier[i], g(o2, "Brier score"))
  add("Table 2", paste(m, "threshold"), e2$thr[i], g(o2, "Threshold"))
  for (col in c(Sensitivity = "sens", Specificity = "spec", PPV = "ppv", NPV = "npv") |> names())
    add("Table 2", paste(m, col),
        e2[[c(Sensitivity = "sens", Specificity = "spec", PPV = "ppv", NPV = "npv")[[col]]]][i],
        g(o2, col))
  add("Table S8", paste(m, "accuracy"), e2$acc[i], g(o8, "Accuracy"))
  add("Table S8", paste(m, "F1"), e2$f1[i], g(o8, "F1"))
  for (col in c("Sensitivity", "Specificity", "PPV", "NPV"))
    add("Table S8", paste(m, col), g(o2, col), g(o8, col))
}

# ---- Supplementary Table S2 ----------------------------------------------------------
e_s2 <- c("Diarrhoea duration (days)" = 50, "Faecal macrophage grade" = 50, "Residence" = 50,
          "Mucus in stool" = 47, "Pre-admission antibiotic category" = 43,
          "Stool white-cell count" = 41, "Blood in stool" = 38, "Straining/tenesmus" = 36,
          "Age (months)" = 29, "Stool pus-cell count" = 25, "Stool red-cell count" = 21,
          "Height-for-age Z-score" = 18, "Stool consistency" = 14, "Weight-for-age Z-score" = 13,
          "Weight-for-height Z-score" = 13, "Respiratory rate" = 12, "Stool colour: yellow" = 9,
          "Sex" = 7, "Stool colour: red/bloody" = 6, "Fever" = 2, "Stool colour: green" = 2,
          "Cough" = 0, "Stool colour: brown" = 0)
s2 <- tb("TableS2_feature_selection_stability.csv")
for (v in names(e_s2)) {
  o <- s2[s2$Variable == v, ]
  add("Table S2", v, sprintf("%d/50 (%s)", e_s2[[v]], ifelse(e_s2[[v]] >= 25, "Included", "Excluded")),
      if (nrow(o) == 1) sprintf("%s (%s)", o$`Outer-fold selections`, o$`Final predictor set (>=50%)`) else NA)
}

# ---- Supplementary Table S3 (subgroups) ----------------------------------------------
e_s3 <- read.csv(text = "
sub,n,events,auroc,slope,brier,sens,spec
All participants,329,132,0.71 (0.66-0.76),1.02,0.208,0.61,0.69
Male,176,68,0.75 (0.69-0.81),1.25,0.196,0.68,0.72
Female,153,64,0.67 (0.59-0.74),0.81,0.222,0.53,0.65
24-35 months,149,55,0.70 (0.62-0.77),0.90,0.209,0.61,0.68
36-47 months,109,45,0.72 (0.64-0.80),1.12,0.206,0.59,0.70
48-59 months,71,32,0.73 (0.63-0.83),1.13,0.208,0.64,0.68
Urban,264,92,0.66 (0.61-0.72),0.80,0.212,0.51,0.70
Rural,65,40,0.80 (0.69-0.89),1.53,0.188,0.83,0.60
Stunted (HAZ < -2),46,23,0.75 (0.62-0.87),1.27,0.206,0.75,0.58
Not stunted,283,109,0.70 (0.65-0.75),0.97,0.208,0.58,0.70
Underweight (WAZ < -2),96,36,0.68 (0.59-0.78),0.90,0.212,0.64,0.61
Not underweight,233,96,0.73 (0.67-0.78),1.10,0.206,0.60,0.72
Wasted (WHZ < -2),78,33,0.67 (0.57-0.77),0.88,0.222,0.64,0.63
Not wasted,251,99,0.72 (0.66-0.77),1.07,0.203,0.60,0.70", colClasses = "character")
s3 <- tb("TableS3_subgroup_performance.csv")
for (i in seq_len(nrow(e_s3))) {
  o <- s3[s3$Subgroup == e_s3$sub[i], ]
  if (nrow(o) != 1) { add("Table S3", e_s3$sub[i], "row", NA); next }
  add("Table S3", paste(e_s3$sub[i], "n / events"), paste(e_s3$n[i], e_s3$events[i]),
      paste(o$n, o$Events))
  add("Table S3", paste(e_s3$sub[i], "AUROC (95% CI)"), e_s3$auroc[i], o$`AUROC (95% CI)`)
  add("Table S3", paste(e_s3$sub[i], "calibration slope"), e_s3$slope[i], num_eq(o$`Calibration slope`, 2))
  add("Table S3", paste(e_s3$sub[i], "Brier"), e_s3$brier[i], num_eq(o$`Brier score`, 3))
  add("Table S3", paste(e_s3$sub[i], "sensitivity"), e_s3$sens[i], num_eq(o$Sensitivity, 2))
  add("Table S3", paste(e_s3$sub[i], "specificity"), e_s3$spec[i], num_eq(o$Specificity, 2))
}

# ---- Supplementary Tables S4a-S4d -----------------------------------------------------
e_s4a <- read.csv(text = "
model,auroc,sd,range,won
Elastic Net,0.720 (0.665-0.772),0.012,0.698-0.734,4
XGBoost,0.716 (0.669-0.761),0.011,0.702-0.736,2
Support vector machine,0.711 (0.659-0.761),0.020,0.688-0.745,2
Random Forest (final model),0.711 (0.660-0.759),0.014,0.682-0.731,1
Neural network,0.699 (0.653-0.743),0.019,0.667-0.733,1
Naive Bayes,0.687 (0.636-0.736),0.019,0.665-0.721,0
k-nearest neighbours,0.687 (0.635-0.734),0.015,0.666-0.721,0
Decision tree,0.676 (0.631-0.720),0.018,0.644-0.697,0", colClasses = "character", encoding = "UTF-8")
s4a <- tb("TableS4a_auroc_all_models.csv")
for (i in seq_len(nrow(e_s4a))) {
  o <- s4a[norm(s4a$Model) == e_s4a$model[i], ]
  add("Table S4a", paste(e_s4a$model[i], "mean AUROC (95% CI)"), e_s4a$auroc[i], o$`Mean AUROC (95% CI)`)
  add("Table S4a", paste(e_s4a$model[i], "SD"), e_s4a$sd[i], num_eq(o$`SD across repeats`, 3))
  add("Table S4a", paste(e_s4a$model[i], "range"), e_s4a$range[i], o$`Range across repeats`)
  add("Table S4a", paste(e_s4a$model[i], "repeats ranked first"), e_s4a$won[i], o$`Repeats ranked first`)
}

e_s4b <- read.csv(text = "
model,d,ci,p
Elastic Net,-0.009,-0.036 to 0.020,0.528
XGBoost,-0.005,-0.023 to 0.013,0.617
Support vector machine,0.000,-0.023 to 0.023,0.991
Neural network,0.012,-0.011 to 0.035,0.293
Naive Bayes,0.024,-0.005 to 0.053,0.102
k-nearest neighbours,0.025,0.003 to 0.046,0.023
Decision tree,0.035,0.009 to 0.061,0.008", colClasses = "character", encoding = "UTF-8")
s4b <- tb("TableS4b_paired_vs_final_model.csv")
for (i in seq_len(nrow(e_s4b))) {
  o <- s4b[norm(s4b$Comparison) == paste("Random Forest vs", e_s4b$model[i]), ]
  add("Table S4b", paste("RF vs", e_s4b$model[i]), paste(e_s4b$d[i], e_s4b$ci[i], e_s4b$p[i]),
      if (nrow(o) == 1) paste(sub("^-0.000$", "0.000", o$dAUROC), o$`95% CI`, o$`p value`) else NA)
}

e_s4c <- read.csv(text = "
model,d,slope,hl,nb,outcome
Elastic Net,Leader,0.78,0.242,0.154,Candidate; not selected at step 5
XGBoost,0.004 (-0.028 to 0.037); p=0.818,0.63,0.001,0.148,Excluded at step 4 (median HL p <=0.05)
Support vector machine,0.009 (-0.011 to 0.026); p=0.367,0.87,0.337,0.153,Candidate; not selected at step 5
Random Forest,0.009 (-0.020 to 0.036); p=0.528,1.05,0.279,0.150,Selected (slope closest to 1)
Neural network,0.021 (0.001 to 0.040); p=0.036,-,-,-,Excluded at step 3 (differs from leader)
Naive Bayes,0.033 (-0.004 to 0.069); p=0.079,0.40,0.000,0.071,Excluded at step 4 (median HL p <=0.05)
k-nearest neighbours,0.033 (0.006 to 0.060); p=0.012,-,-,-,Excluded at step 3 (differs from leader)
Decision tree,0.044 (0.006 to 0.080); p=0.022,-,-,-,Excluded at step 3 (differs from leader)",
  colClasses = "character", encoding = "UTF-8")
s4c <- tb("TableS4c_model_selection_trace.csv")
for (i in seq_len(nrow(e_s4c))) {
  o <- s4c[norm(s4c$Model) == e_s4c$model[i], ]
  add("Table S4c", e_s4c$model[i],
      paste(e_s4c$d[i], e_s4c$slope[i], e_s4c$hl[i], e_s4c$nb[i], e_s4c$outcome[i], sep = " | "),
      if (nrow(o) == 1) paste(o[[3]], o[[4]], o[[5]], o[[6]], o[[7]], sep = " | ") else NA)
}

e_s4d <- read.csv(text = "
r,Elastic Net,XGBoost,SVM,Random Forest,Neural Network,Naive Bayes,KNN,Decision Tree,winner
1 (223),0.710,0.710,0.695,0.701,0.700,0.665,0.673,0.644,Elastic Net
2 (323),0.710,0.713,0.698,0.682,0.701,0.670,0.685,0.661,XGBoost
3 (423),0.731,0.711,0.721,0.718,0.709,0.680,0.685,0.668,Elastic Net
4 (523),0.725,0.702,0.715,0.706,0.677,0.668,0.688,0.678,Elastic Net
5 (623),0.727,0.727,0.745,0.728,0.716,0.714,0.721,0.696,Support vector machine
6 (723),0.698,0.712,0.688,0.720,0.667,0.677,0.666,0.688,Random Forest
7 (823),0.727,0.736,0.743,0.731,0.703,0.721,0.687,0.695,Support vector machine
8 (923),0.709,0.715,0.695,0.706,0.697,0.683,0.689,0.657,XGBoost
9 (1023),0.734,0.706,0.711,0.714,0.686,0.701,0.699,0.697,Elastic Net
10 (1123),0.727,0.729,0.702,0.708,0.733,0.691,0.673,0.677,Neural network",
  colClasses = "character", check.names = FALSE)
s4d <- tb("TableS4d_auroc_by_repeat.csv")
for (i in seq_len(nrow(e_s4d))) {
  o <- s4d[s4d$Repeat == e_s4d$r[i], ]
  for (m in setdiff(names(e_s4d), c("r", "winner")))
    add("Table S4d", paste("repeat", e_s4d$r[i], m), e_s4d[[m]][i], if (nrow(o) == 1) o[[m]] else NA)
  add("Table S4d", paste("repeat", e_s4d$r[i], "highest"), e_s4d$winner[i],
      if (nrow(o) == 1) o$Highest else NA)
}

# ---- Supplementary Tables S5a, S5b, S6a, S6b ------------------------------------------
e_s5a <- read.csv(text = "
reg,f,s,t,rate
Ciprofloxacin + Cefuroxime,1,0,1,100
Ciprofloxacin + Erythromycin,1,0,1,100
Erythromycin,1,0,1,100
Ciprofloxacin + Metronidazole,22,13,35,62.9
Cefixime + Metronidazole,1,1,2,50
Ciprofloxacin,40,41,81,49.4
No Antibiotic,38,74,112,33.9
Metronidazole,27,58,85,31.8
Cefixime,1,8,9,11.1
Cefixime + Amoxicillin,0,1,1,0
Ciprofloxacin + Amoxicillin,0,1,1,0", colClasses = "character")
s5a <- tb("TableS5a_pre_admission_regimen.csv")
for (i in seq_len(nrow(e_s5a))) {
  o <- s5a[tolower(s5a$`Pre-admission regimen`) == tolower(e_s5a$reg[i]), ]
  add("Table S5a", e_s5a$reg[i], paste(e_s5a$f[i], e_s5a$s[i], e_s5a$t[i], e_s5a$rate[i]),
      if (nrow(o) == 1) paste(o$Failure, o$Success, o$Total, o$`Failure rate (%)`) else NA)
}
e_s5b <- read.csv(text = "
cat,f,s,t,rate
None,38,74,112,33.9
Ciprofloxacin,40,41,81,49.4
Metronidazole,27,58,85,31.8
Ciprofloxacin + metronidazole,22,13,35,62.9
Other,5,11,16,31.2", colClasses = "character")
s5b <- tb("TableS5b_antibiotic_categories.csv")
for (i in seq_len(nrow(e_s5b))) {
  o <- s5b[s5b$Category == e_s5b$cat[i], ]
  add("Table S5b", e_s5b$cat[i], paste(e_s5b$f[i], e_s5b$s[i], e_s5b$t[i], e_s5b$rate[i]),
      if (nrow(o) == 1) paste(o$Failure, o$Success, o$Total, o$`Failure rate (%)`) else NA)
}
e_s6a <- c("Bloody stool at 48 h" = "118 89.4", "Watery stool at 48 h" = "26 19.7",
           "Persistent fever at 48 h" = "13 9.8", "Clinical deterioration" = "0 0")
s6a <- tb("TableS6a_outcome_components.csv")
for (c6 in names(e_s6a)) {
  o <- s6a[s6a$Component == c6, ]
  add("Table S6a", c6, e_s6a[[c6]], if (nrow(o) == 1) paste(o$n, o$`% of failures`) else NA)
}
e_s6b <- c("Bloody stool at 48 h" = 97, "Bloody stool at 48 h + Watery stool at 48 h" = 12,
           "Watery stool at 48 h" = 10, "Bloody stool at 48 h + Persistent fever at 48 h" = 6,
           "Bloody stool at 48 h + Watery stool at 48 h + Persistent fever at 48 h" = 3,
           "Persistent fever at 48 h" = 3, "Watery stool at 48 h + Persistent fever at 48 h" = 1)
s6b <- tb("TableS6b_component_combinations.csv")
for (c6 in names(e_s6b)) {
  o <- s6b[s6b$`Combination of components` == c6, ]
  add("Table S6b", c6, e_s6b[[c6]], if (nrow(o) == 1) o$n else NA)
}

# ---- Supplementary Table S7 -------------------------------------------------------------
e_s7 <- read.csv(text = "
model,int,slope,brier,hl,ok
Elastic Net,-0.01,0.78 (0.69-0.85),0.206,0.242,8/10
XGBoost,0.03,0.63 (0.50-0.85),0.211,<0.001,2/10
Support vector machine,-0.04,0.87 (0.68-1.10),0.208,0.337,9/10
Random Forest (final model),-0.01,1.05 (0.82-1.26),0.208,0.279,8/10
Neural network,-0.24,1.23 (1.06-1.43),0.215,0.106,8/10
Naive Bayes,1.77,0.40 (0.29-0.49),0.278,<0.001,0/10
k-nearest neighbours,0.15,0.82 (0.28-1.07),0.216,0.130,9/10
Decision tree,0.02,0.22 (0.09-0.44),0.238,<0.001,0/10", colClasses = "character", encoding = "UTF-8")
s7 <- tb("TableS7_calibration_all_models.csv")
for (i in seq_len(nrow(e_s7))) {
  o <- s7[norm(s7$Model) == e_s7$model[i], ]
  add("Table S7", e_s7$model[i],
      paste(e_s7$int[i], e_s7$slope[i], e_s7$brier[i], e_s7$hl[i], e_s7$ok[i], sep = " | "),
      if (nrow(o) == 1) paste(o[[2]], o[[3]], o[[4]], o[[5]], o[[6]], sep = " | ") else NA)
}

# ---- Supplementary Table S9 -------------------------------------------------------------
e_s9 <- read.csv(text = "
thr,strategy,nb,inc
10%,Treat all,0.3347,Reference
10%,Random Forest (final model),0.3367,0.0020
10%,Elastic Net,0.3335,-0.0012
10%,XGBoost,0.3306,-0.0041
20%,Treat all,0.2515,Reference
20%,Random Forest (final model),0.2636,0.0121
20%,Elastic Net,0.2571,0.0056
20%,XGBoost,0.2544,0.0029
30%,Treat all,0.1446,Reference
30%,Random Forest (final model),0.1810,0.0364
30%,Elastic Net,0.1895,0.0449
30%,XGBoost,0.1782,0.0337
40%,Treat all,0.0020,Reference
40%,Random Forest (final model),0.1208,0.1187
40%,Elastic Net,0.1219,0.1199
40%,XGBoost,0.1223,0.1203
50%,Treat all,-0.1976,Reference
50%,Random Forest (final model),0.0739,0.0739
50%,Elastic Net,0.0793,0.0793
50%,XGBoost,0.0784,0.0784", colClasses = "character")
s9 <- tb("TableS9_dca_net_benefit.csv")
s9$thr <- s9$`Threshold probability`
for (j in seq_len(nrow(s9))) if (s9$thr[j] == "") s9$thr[j] <- s9$thr[j - 1]
for (i in seq_len(nrow(e_s9))) {
  o <- s9[s9$thr == e_s9$thr[i] & s9$`Strategy / model` == e_s9$strategy[i], ]
  add("Table S9", paste(e_s9$thr[i], e_s9$strategy[i]), paste(e_s9$nb[i], e_s9$inc[i]),
      if (nrow(o) == 1) paste(o$`Net benefit`, o$`Incremental net benefit vs treat all`) else NA)
}

# ---- Supplementary Table S10 (kernelSHAP) -------------------------------------------------
e_s10 <- read.csv(text = "
pred,abs,pos,neg,cor
Faecal macrophage grade,0.162,0.205,-0.134,-0.98
Stool white-cell count,0.072,0.068,-0.078,0.65
Residence,0.046,0.115,-0.028,-0.96
Pre-admission antibiotic category,0.045,0.055,-0.038,NA
Age (months),0.041,0.042,-0.040,0.16
Diarrhoea duration (days),0.032,0.034,-0.030,0.28
Straining/tenesmus,0.032,0.041,-0.026,0.89
Mucus in stool,0.024,0.031,-0.019,NA
Blood in stool,0.023,0.042,-0.016,NA
Stool pus-cell count,0.012,0.008,-0.027,0.84", colClasses = "character")
s10 <- tb("TableS10_SHAP_summary.csv")
add("Table S10", "order by mean absolute SHAP", paste(e_s10$pred, collapse = " > "),
    paste(s10$Predictor, collapse = " > "))
for (i in seq_len(nrow(e_s10))) {
  o <- s10[s10$Predictor == e_s10$pred[i], ]
  if (nrow(o) != 1) { add("Table S10", e_s10$pred[i], "row", NA); next }
  add("Table S10", e_s10$pred[i],
      paste(e_s10$abs[i], e_s10$pos[i], e_s10$neg[i], e_s10$cor[i]),
      paste(num_eq(o[[2]], 3), num_eq(o[[3]], 3), num_eq(o[[4]], 3),
            ifelse(o[[5]] %in% c("", "NA"), "NA", num_eq(o[[5]], 2))))
}

# ---- Supplementary Table S1 (data dictionary) ---------------------------------------------
s1 <- tb("TableS1_data_dictionary.csv")
raw_s1 <- strsplit(trimws(strsplit(
  "azithromycin_failure|2|No: n=197 (59.9%); Yes: n=132 (40.1%)
age|36|mean 38.35, SD 10.17, min 24.00, max 59.00
sex|2|Male: n=176 (53.5%); Female: n=153 (46.5%)
resident_area|2|Urban: n=264 (80.2%); Rural: n=65 (19.8%)
diarrhea_dur|11|mean 2.76, SD 1.97, min 1.00, max 12.00
fever|2|Yes: n=176 (53.5%); No: n=153 (46.5%)
cough|2|No: n=269 (81.8%); Yes: n=60 (18.2%)
resp_rate|20|mean 30.58, SD 3.19, min 23.00, max 48.00
straining|2|No: n=202 (61.4%); Yes: n=127 (38.6%)
haz|229|mean -0.92, SD 1.07, min -5.03, max 2.58
waz|212|mean -1.40, SD 1.03, min -4.29, max 1.52
whz|205|mean -1.31, SD 0.93, min -3.67, max 1.56
consist_RE|3|Mucoid: n=130 (39.5%); Formed: n=101 (30.7%); Liquid: n=98 (29.8%)
blood_RE|4|Small: n=237 (72.0%); Absent: n=72 (21.9%); Moderate: n=12 (3.6%); Plenty: n=8 (2.4%)
mucus_RE|3|Moderate: n=178 (54.1%); Plenty: n=129 (39.2%); Small: n=22 (6.7%)
abx_cat|5|None: n=112 (34.0%); Metronidazole: n=85 (25.8%); Ciprofloxacin: n=81 (24.6%); Ciprofloxacin + metronidazole: n=35 (10.6%); Other: n=16 (4.9%)
wbc_numeric|7|mean 12.01, SD 5.43, min 3.00, max 35.00
macro_numeric|2|mean 2.71, SD 0.98, min 1.50, max 3.50
pcell_numeric|7|mean 30.99, SD 8.09, min 3.00, max 35.00
rbc_numeric|7|mean 28.81, SD 10.11, min 3.00, max 35.00
is_green|2|mean 0.53, SD 0.50, min 0.00, max 1.00
is_red|2|mean 0.36, SD 0.48, min 0.00, max 1.00
is_yellow|2|mean 0.34, SD 0.48, min 0.00, max 1.00
is_brown|2|mean 0.02, SD 0.13, min 0.00, max 1.00", "\n")[[1]]), "|", fixed = TRUE)
for (r in raw_s1) {
  o <- s1[s1$Variable == r[1], ]
  add("Table S1", r[1], paste(r[2], r[3]),
      if (nrow(o) == 1) paste(o$`Distinct values`, o$Distribution) else NA)
}

# ---- Result ------------------------------------------------------------------------------
res <- bind_rows(checks) |>
  mutate(status = ifelse(!is.na(observed) & norm(expected) == norm(observed), "MATCH", "DIFFERENT"))
dir.create(OUT_DIR, showWarnings = FALSE)
write.csv(res, file.path(OUT_DIR, "verification_report.csv"), row.names = FALSE)

summary_tab <- res |> group_by(section) |>
  summarise(checked = n(), match = sum(status == "MATCH"), different = sum(status != "MATCH"),
            .groups = "drop")
md <- c("# Verification against the published paper", "",
        sprintf("Checked %d reported values: %d match, %d differ.", nrow(res),
                sum(res$status == "MATCH"), sum(res$status != "MATCH")), "",
        "| Section | Checked | Match | Different |", "|---|---|---|---|",
        sprintf("| %s | %d | %d | %d |", summary_tab$section, summary_tab$checked,
                summary_tab$match, summary_tab$different))
diffs <- res |> filter(status != "MATCH")
if (nrow(diffs) > 0) {
  md <- c(md, "", "## Values that differ", "", "| Section | Item | Paper | Code |", "|---|---|---|---|",
          sprintf("| %s | %s | %s | %s |", diffs$section, diffs$item, diffs$expected, diffs$observed))
}
writeLines(md, file.path(OUT_DIR, "verification_report.md"))
cat(md, sep = "\n")
