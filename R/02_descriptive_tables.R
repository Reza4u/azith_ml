# =============================================================================
# 02_descriptive_tables.R — Table 1 and Supplementary Tables S1, S5a-b, S6a-b
# =============================================================================

# ---- Supplementary Table S1: data dictionary of the analysis dataset ---------
display_level <- function(x) {
  x <- as.character(x)
  x <- sub(" Area$", "", x)                                   # "Urban Area" -> "Urban"
  x <- ifelse(x == "Cipro_Metro", "Ciprofloxacin + metronidazole", x)
  x <- ifelse(tolower(x) == "plenty", "Plenty", x)
  x
}
table_s1 <- purrr::map_dfr(names(ml_df), function(v) {
  x <- ml_df[[v]]
  if (is.numeric(x)) {
    tibble(Variable = v, Label = plab(v), Type = "numeric", Missing = sum(is.na(x)),
           `Distinct values` = dplyr::n_distinct(x),
           Distribution = sprintf("mean %.2f, SD %.2f, min %.2f, max %.2f",
                                  mean(x), sd(x), min(x), max(x)))
  } else {
    tb <- sort(table(display_level(x)), decreasing = TRUE)
    tibble(Variable = v, Label = plab(v), Type = "categorical", Missing = sum(is.na(x)),
           `Distinct values` = length(tb),
           Distribution = paste(sprintf("%s: n=%d (%.1f%%)", names(tb), as.integer(tb),
                                        100 * as.integer(tb) / sum(tb)), collapse = "; "))
  }
})
write.csv(table_s1, tab_path("TableS1_data_dictionary.csv"), row.names = FALSE)

# ---- Table 1: baseline characteristics by outcome ----------------------------
tbl_labels <- list(
  age = "Age, months", sex = "Sex, male", resident_area = "Residence, urban",
  diarrhea_dur = "Diarrhoea duration, days", fever = "Fever", cough = "Cough",
  resp_rate = "Respiratory rate", straining = "Straining/tenesmus",
  haz = "Height-for-age Z-score", waz = "Weight-for-age Z-score",
  whz = "Weight-for-height Z-score", abx_cat = "Pre-admission antibiotic use",
  consist_RE = "Stool consistency", blood_RE = "Blood in stool", mucus_RE = "Mucus in stool",
  wbc_numeric = "Stool white-cell count", macro_grade = "Faecal macrophages, 2-5 per HPF",
  pcell_numeric = "Stool pus-cell count", rbc_numeric = "Stool red-cell count",
  is_green = "Stool colour: green", is_red = "Stool colour: red/bloody",
  is_yellow = "Stool colour: yellow", is_brown = "Stool colour: brown")

set.seed(SEED)   # Monte Carlo Fisher p value for the five-category antibiotic variable
tbl1 <- ml_df |>
  # macrophages: two-level scale, shown as n (%) with the higher grade
  mutate(macro_grade = factor(ifelse(macro_numeric < 2, "1-2 per HPF", "2-5 per HPF"),
                              levels = c("1-2 per HPF", "2-5 per HPF"))) |>
  dplyr::select(-macro_numeric) |>
  mutate(
    consist_RE = factor(consist_RE, levels = c("Formed", "Mucoid", "Liquid")),
    blood_RE   = factor(tools::toTitleCase(tolower(as.character(blood_RE))),
                        levels = c("Absent", "Small", "Moderate", "Plenty")),
    mucus_RE   = factor(tools::toTitleCase(tolower(as.character(mucus_RE))),
                        levels = c("Small", "Moderate", "Plenty")),
    abx_cat    = forcats::fct_recode(abx_cat, "Ciprofloxacin + metronidazole" = "Cipro_Metro")
  ) |>
  tbl_summary(
    by = azithromycin_failure,
    type = list(all_dichotomous() ~ "dichotomous",
                c(wbc_numeric, pcell_numeric, rbc_numeric) ~ "continuous"),
    value = list(sex ~ "Male", resident_area ~ "Urban Area", macro_grade ~ "2-5 per HPF"),
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
    label = tbl_labels
  ) |>
  add_p(
    test = list(abx_cat ~ "fisher.test"),
    test.args = list(abx_cat ~ list(simulate.p.value = TRUE, B = 1e5)),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) |>
  add_overall()

table1 <- as_tibble(tbl1, col_labels = FALSE) |>
  mutate(variable = tbl1$table_body$variable, row_type = tbl1$table_body$row_type, .before = 1) |>
  dplyr::rename(Variable = label, `Overall (n=329)` = stat_0,
                Failure = stat_1, Success = stat_2, `p value` = p.value)
write.csv(table1, tab_path("Table1_baseline.csv"), row.names = FALSE)

# ---- Supplementary Tables S5a and S5b: pre-admission antibiotic regimens -----
regimen <- df |>
  filter(!is.na(azithromycin_failure)) |>
  mutate(azithromycin_failure = factor(azithromycin_failure, levels = c("Yes", "No")))
reg_tab <- table(regimen$antibiotic_home, regimen$azithromycin_failure)

table_s5a <- as.data.frame.matrix(reg_tab) |>
  tibble::rownames_to_column("Pre-admission regimen") |>
  as_tibble() |>
  mutate(Total = Yes + No, `Failure rate (%)` = round(100 * Yes / Total, 1)) |>
  arrange(desc(`Failure rate (%)`)) |>
  dplyr::select(`Pre-admission regimen`, Failure = Yes, Success = No, Total, `Failure rate (%)`)

table_s5b <- ml_df |>
  count(abx_cat, azithromycin_failure) |>
  pivot_wider(names_from = azithromycin_failure, values_from = n, values_fill = 0) |>
  mutate(abx_cat = forcats::fct_recode(abx_cat, "Ciprofloxacin + metronidazole" = "Cipro_Metro"),
         Total = Yes + No, `Failure rate (%)` = round(100 * Yes / Total, 1)) |>
  dplyr::select(Category = abx_cat, Failure = Yes, Success = No, Total, `Failure rate (%)`)

# Ciprofloxacin-containing categories (ciprofloxacin; ciprofloxacin + metronidazole)
# by outcome (Results text; footnote of Supplementary Table S5b)
cipro_rows <- ml_df$abx_cat %in% c("Ciprofloxacin", "Cipro_Metro")
cipro_containing <- c(
  failure_n = sum(cipro_rows & y_all),  failure_N = sum(y_all),
  success_n = sum(cipro_rows & !y_all), success_N = sum(!y_all))

set.seed(SEED)
reg_fisher <- fisher.test(reg_tab, simulate.p.value = TRUE, B = 1e5)
set.seed(SEED)
cat_fisher <- fisher.test(table(ml_df$abx_cat, ml_df$azithromycin_failure),
                          simulate.p.value = TRUE, B = 1e5)

write.csv(table_s5a, tab_path("TableS5a_pre_admission_regimen.csv"), row.names = FALSE)
write.csv(table_s5b, tab_path("TableS5b_antibiotic_categories.csv"), row.names = FALSE)

# ---- Supplementary Tables S6a and S6b: components of the composite outcome ----
# The 48-hour assessment fields are held in the trial dataset (modeldata.xlsx).
md_raw <- readxl::read_excel(file.path(DATA_DIR, "modeldata.xlsx"), sheet = "Sheet1")
comp_cols <- c(`Bloody stool at 48 h`     = "blood_48h",
               `Watery stool at 48 h`     = "watery_stool_48h",
               `Persistent fever at 48 h` = "fever_48h",
               `Clinical deterioration`   = "clin_deterior_48h")
stopifnot(all(comp_cols %in% names(md_raw)))
fails <- md_raw |> filter(azithromycin_failure == "Yes")
comp  <- as.data.frame(sapply(comp_cols, function(cc) !is.na(fails[[cc]]) & fails[[cc]] == "Yes"))
names(comp) <- names(comp_cols)

table_s6a <- tibble(Component = names(comp), n = colSums(comp),
                    `% of failures` = round(100 * colSums(comp) / nrow(comp), 1))
table_s6b <- as.data.frame(sort(table(
  apply(comp, 1, function(r) paste(names(comp)[r], collapse = " + "))), decreasing = TRUE))
names(table_s6b) <- c("Combination of components", "n")
components_met <- table(rowSums(comp))

write.csv(table_s6a, tab_path("TableS6a_outcome_components.csv"), row.names = FALSE)
write.csv(table_s6b, tab_path("TableS6b_component_combinations.csv"), row.names = FALSE)

cat("Table 1, S1, S5a-b and S6a-b written. Antibiotic categories, Fisher (Monte Carlo) p =",
    signif(cat_fisher$p.value, 3), "\n")
