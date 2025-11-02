################################################################################################################################
#########################################    Script generating       ##########################################################3
#########################################    Tables OF THE PAPER      ##########################################################3
################################################################################################################################

syn_data <- syn_data %>%
  group_by(trial_name) %>%
  mutate(
    arm = dense_rank(avg_rando_timing_DOAC)  # 1 = lowest value in that trial
  ) %>%
  ungroup()

syn_data$arm <- with(
  syn_data,
  ave(avg_rando_timing_DOAC, trial_name,
      FUN = function(x) match(x, sort(unique(x))))
)


## -------------------------
## 1) Arm assignment (per trial)
## -------------------------
syn_data <- syn_data %>%
  group_by(trial_name) %>%
  mutate(
    arm = dense_rank(avg_rando_timing_DOAC)  # 1 = lowest targeted timing within trial
  ) %>%
  ungroup()

## -------------------------
## 2) Helpers for formatting
## -------------------------
fmt_mean_sd <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  paste0(round(mean(x, na.rm = TRUE), 2), " (", round(sd(x, na.rm = TRUE), 2), ")")
}
fmt_median_iqr <- function(x) {
  if (all(is.na(x))) return(NA_character_)
  q <- quantile(x, probs = c(.25, .5, .75), na.rm = TRUE, names = FALSE)
  paste0(round(q[2], 2), " [", round(q[1], 2), "–", round(q[3], 2), "]")
}
fmt_count_pct <- function(x) {
  n_events <- sum(x == 1, na.rm = TRUE)
  total    <- sum(!is.na(x))
  if (total == 0) return(NA_character_)
  paste0(n_events, " (", round(100 * n_events / total, 1), "%)")
}

fmt_pct <- function(x) {
  # Works for logical or 0/1
  if (is.logical(x)) p <- mean(x, na.rm = TRUE)
  else               p <- mean(x == 1, na.rm = TRUE)
  paste0(round(100*p, 1), "%")
}

## If your dataset uses slightly different names, these lines pick what exists.
comp360_var <- intersect(c("composite_30","composite_30"), names(syn_data))[1]
comp90_var  <- intersect(c("composite_90"), names(syn_data))[1]
isch30_var  <- intersect(c("ischemic_stroke_30","ischemic_30"), names(syn_data))[1]
sich30_var  <- intersect(c("SICH_stroke_30","sich_30"), names(syn_data))[1]
mort30_var  <- intersect(c("mortality_30","death_30"), names(syn_data))[1]

## -------------------------
## 3) Summaries per Trial x Arm
## -------------------------
sum_by_arm <- syn_data %>%
  arrange(trial_name, arm, avg_rando_timing_DOAC) %>%
  group_by(trial_name, arm) %>%
  summarise(
    `Sample size`                                     = n(),
    `Targeted timing (avg_rando_timing_DOAC)`         = first(avg_rando_timing_DOAC),
    `Observed timing mean (SD) (timing_DOAC)`         = fmt_mean_sd(timing_DOAC),
    `Missing observed timing (timing_DOAC)` = {
      miss_n <- sum(is.na(timing_DOAC))
      total  <- n()
      paste0(miss_n, " (", round(100 * miss_n / total, 1), "%)")
    },
    `Composite 360 days`                              = if (!is.na(comp360_var)) fmt_count_pct(.data[[comp360_var]]) else NA_character_,
    `Composite 90 days`                               = if (!is.na(comp90_var))  fmt_count_pct(.data[[comp90_var]])  else NA_character_,
    `Recurrent ischaemic stroke 30 days`              = if (!is.na(isch30_var))  fmt_count_pct(.data[[isch30_var]])  else NA_character_,
    `Symptomatic ICH 30 days`                         = if (!is.na(sich30_var))  fmt_count_pct(.data[[sich30_var]])  else NA_character_,
    `All-cause mortality 30 days`                     = if (!is.na(mort30_var))  fmt_count_pct(.data[[mort30_var]])  else NA_character_,
    `Age mean (SD)`                                   = fmt_mean_sd(Age),
    `% Females (sex_factor=="Female")`                = fmt_count_pct(sex_factor == "Female"),
    `NIHSS at admission median [IQR]`                 = fmt_median_iqr(nihss_adm),
    .groups = "drop"
  )




report_tbl <- sum_by_arm %>%
  tidyr::pivot_longer(
    cols = -c(trial_name, arm),
    names_to = "Row",
    values_to = "Value",
    values_transform = list(Value = as.character),  # <- make everything character
    values_ptypes     = list(Value = character())
  ) %>%
  dplyr::mutate(TrialArm = paste0(trial_name, " - Arm ", arm)) %>%
  dplyr::select(Row, TrialArm, Value) %>%
  tidyr::pivot_wider(names_from = TrialArm, values_from = Value)



## Optional: order the rows like your screenshot
row_order <- c(
  "Sample size",
  "Targeted timing (avg_rando_timing_DOAC)",
  "Observed timing mean (SD) (timing_DOAC)",
  "Missing observed timing (timing_DOAC)",
  "Composite 360 days",
  "Composite 90 days",
  "Recurrent ischaemic stroke 30 days",
  "Symptomatic ICH 30 days",
  "Mortality 30 days",
  "Age mean (SD)",
  "% Females (sex_factor==\"Female\")",
  "NIHSS at admission median [IQR]"
)

report_tbl <- sum_by_arm %>%
  pivot_longer(
    cols = -c(trial_name, arm),
    names_to = "Row",
    values_to = "Value",
    values_transform = list(Value = as.character),
    values_ptypes = list(Value = character())
  ) %>%
  filter(!is.na(Row)) %>%              # <-- remove the spurious NA row
  mutate(TrialArm = paste0(trial_name, " - Arm ", arm)) %>%
  select(Row, TrialArm, Value) %>%
  pivot_wider(names_from = TrialArm, values_from = Value)


report_tbl
# Save to Excel
write_xlsx(report_tbl, "Tables/Table1.xlsx")


