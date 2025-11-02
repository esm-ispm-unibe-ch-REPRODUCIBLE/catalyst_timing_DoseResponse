
#####################################################################################################
############  Examine association of discrepancies   ###############################33 #################
############ between observed and randomized timing & NIHSS at admission  ########################3
#######################################################################################3

#### Create the ne variable that will be
# A. for early arms: time_difference = 0 if timing_DOAC<=than latest intended time of this arm
#    and time_difference=timing_DOAC-latest intended time of this arm if timing_DOAC later than intended
# B. for later arms: time_difference=0 if timing_DOAC>= than the earliest intended time for this arm
#   time_difference=timing_DOAC-earliest intended time of this arm if timing DOAC earlier than intended
syn_data$time_difference <- NA_real_

syn_data <- syn_data %>%
  mutate(
    time_difference = case_when(
      trial_name == "ELAN minor" & treatment_factor == "Early DOAC" & timing_DOAC <= 2 ~ 0,
      trial_name == "ELAN minor" & treatment_factor == "Early DOAC" & timing_DOAC >  2 ~ timing_DOAC - 2,
      trial_name == "ELAN minor" & treatment_factor == "Later DOAC" & timing_DOAC > 2 ~ 0,
      trial_name == "ELAN minor" & treatment_factor == "Later DOAC" & timing_DOAC <=  2 ~ timing_DOAC - 2,
      trial_name == "ELAN moderate" & treatment_factor == "Early DOAC" & timing_DOAC <= 2 ~ 0,
      trial_name == "ELAN moderate" & treatment_factor == "Early DOAC" & timing_DOAC >  2 ~ timing_DOAC - 2,
      trial_name == "ELAN moderate" & treatment_factor == "Later DOAC" & timing_DOAC > 5 ~ 0,
      trial_name == "ELAN moderate" & treatment_factor == "Later DOAC" & timing_DOAC <= 5 ~ timing_DOAC - 5,
      trial_name == "ELAN major" & treatment_factor == "Early DOAC" & timing_DOAC <= 7 ~ 0,
      trial_name == "ELAN major" & treatment_factor == "Early DOAC" & timing_DOAC >  7 ~ timing_DOAC - 7,
      trial_name == "ELAN major" & treatment_factor == "Later DOAC" & timing_DOAC > 11 ~ 0,
      trial_name == "ELAN major" & treatment_factor == "Later DOAC" & timing_DOAC <=  11 ~ timing_DOAC - 11,
      trial_name == "OPTIMAS" & treatment_factor == "Early DOAC" & timing_DOAC <= 4 ~ 0,
      trial_name == "OPTIMAS" & treatment_factor == "Early DOAC" & timing_DOAC >  4 ~ timing_DOAC - 4,
      trial_name == "OPTIMAS" & treatment_factor == "Later DOAC" & timing_DOAC > 6 ~ 0,
      trial_name == "OPTIMAS" & treatment_factor == "Later DOAC" & timing_DOAC <=  6 ~ timing_DOAC - 6,
      trial_name == "TIMING" & treatment_factor == "Early DOAC" & timing_DOAC <= 4 ~ 0,
      trial_name == "TIMING" & treatment_factor == "Early DOAC" & timing_DOAC >  4 ~ timing_DOAC - 4,
      trial_name == "TIMING" & treatment_factor == "Later DOAC" & timing_DOAC > 4 ~ 0,
      trial_name == "TIMING" & treatment_factor == "Later DOAC" & timing_DOAC <=  4 ~ timing_DOAC - 4,
      trial_name == "START" & avg_rando_timing_DOAC == 3 & timing_DOAC <= 4 ~ 0,
      trial_name == "START" & avg_rando_timing_DOAC == 3 & timing_DOAC >  4 ~ timing_DOAC - 4,
      trial_name == "START" & avg_rando_timing_DOAC == 5.5 & timing_DOAC <= 6 & timing_DOAC > 5 ~ 0,
      trial_name == "START" & avg_rando_timing_DOAC == 5.5 & timing_DOAC <  5 ~ timing_DOAC - 5,
      trial_name == "START" & avg_rando_timing_DOAC == 5.5 & timing_DOAC > 6  ~ timing_DOAC - 6,
      trial_name == "START" & avg_rando_timing_DOAC == 9.5 & timing_DOAC <= 10 & timing_DOAC > 9 ~ 0,
      trial_name == "START" & avg_rando_timing_DOAC == 9.5 & timing_DOAC <=  9 ~ timing_DOAC - 9,
      trial_name == "START" & avg_rando_timing_DOAC == 9.5 & timing_DOAC > 10  ~ timing_DOAC - 10,
      trial_name == "START" & avg_rando_timing_DOAC == 13.5 & timing_DOAC > 13 ~ 0,
      trial_name == "START" & avg_rando_timing_DOAC == 13.5 & timing_DOAC <=  13 ~ timing_DOAC - 13,
      TRUE ~ NA_real_
    )
  )

summary(syn_data$time_difference)

syn_data$time_difference_yn <- NA_real_
syn_data$time_difference_yn[syn_data$time_difference==0]<-0
syn_data$time_difference_yn[syn_data$time_difference!=0]<-1
table(syn_data$time_difference_yn)

##### Plots for the association between nihss and time deviation
# your palette
##### Plots for the association between nihss and time deviation
# your palette
trial_cols <- c(
  "ELAN moderate" = "green",
  "ELAN minor"    = "lightgreen",
  "ELAN major"    = "darkgreen",
  "START"         = "#e7298a",
  "TIMING"        = "darkred",
  "OPTIMAS"       = "#e6ab02"
)

desired_order <- c("ELAN minor", "ELAN moderate", "ELAN major",
                   "TIMING", "OPTIMAS", "START")

df <- syn_data %>%
  mutate(
    nihss_adm       = as.numeric(as.character(nihss_adm)),
    time_difference = as.numeric(time_difference),
    trial_name      = factor(as.character(trial_name), levels = desired_order)
  ) %>%
  filter(!is.na(time_difference), !is.na(nihss_adm))

p <- ggplot(df, aes(x = time_difference, y = nihss_adm, color = trial_name)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  facet_wrap(~ trial_name, ncol = 3) +  # respects factor order
  scale_color_manual(values = trial_cols, breaks = desired_order, name = "Trial") +
  labs(
    title = "Time deviation from the intended arm vs NIHSS at admission",
    x = "Time deviation from the intended arm (days)",
    y = "NIHSS at admission"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    strip.text  = element_text(face = "bold"),
    legend.position = "none"
  )

p
### or 

# your custom palette
trial_cols <- c(
  "ELAN moderate" = "green",
  "ELAN minor"    = "lightgreen",
  "ELAN major"    = "darkgreen",
  "START"         = "#e7298a",
  "TIMING"        = "darkred",
  "OPTIMAS"       = "#e6ab02"
)

df <- syn_data %>%
  mutate(
    nihss_adm = as.numeric(as.character(nihss_adm)),
    time_difference = as.numeric(time_difference),
    trial_name = as.character(trial_name)
  ) %>%
  filter(!is.na(time_difference), !is.na(nihss_adm))
desired_order <- c("ELAN minor", "ELAN moderate", "ELAN major",
                   "TIMING", "OPTIMAS", "START")

p2 <- ggplot(df, aes(x = time_difference, y = nihss_adm, color = trial_name)) +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(values = trial_cols, breaks = desired_order, name = "Trial") +
  labs(
    title = "Time deviation from the intended arm vs NIHSS at admission",
    x = "Time deviation from the intended arm (days)",
    y = "NIHSS at admission"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title  = element_text(face = "bold"),
    legend.position = "right"
  )

p2


### analysis

# A. multinomial with 3 categories: less than 0, 0 or more than 0

# --- Prep

dat <- syn_data %>%
  transmute(
    nihss_adm,
    time_difference,
    # tiny tolerance around 0 to avoid floating-point issues; adjust if needed
    time_cat = case_when(
      time_difference < -1e-9 ~ "early",
      abs(time_difference) <= 1e-9 ~ "on_time",
      time_difference > 1e-9 ~ "late",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(nihss_adm), !is.na(time_cat))

dat$time_cat <- relevel(factor(dat$time_cat, levels = c("on_time","early","late")),
                        ref = "on_time")

# Fit model
m_mult <- multinom(time_cat ~ nihss_adm, data = dat)

# Get summary
summ <- summary(m_mult)

# assumes: summ <- summary(m_mult)

coef_term <- summ$coefficients[, "nihss_adm"]
se_term   <- summ$standard.errors[, "nihss_adm"]

lower <- coef_term - 1.96 * se_term
upper <- coef_term + 1.96 * se_term
OR    <- exp(coef_term)
ORlo  <- exp(lower)
ORhi  <- exp(upper)

zvals <- coef_term / se_term
pvals <- 2 * (1 - pnorm(abs(zvals)))

# give names so we can subset by "early"/"late"
nm <- names(coef_term)

OR_CI_fmt <- setNames(sprintf("%.2f (%.2f–%.2f)", OR, ORlo, ORhi), nm)
p_fmt     <- setNames(ifelse(pvals < 0.001, "<0.001", sprintf("%.3f", pvals)), nm)

# desired order & labels
ord  <- c("early", "late")
labs <- c("Earlier vs On time", "Later vs On time")

tbl <- data.frame(
  Comparison   = labs,
  `OR (95% CI)`= OR_CI_fmt[ord],
  `p-value`    = p_fmt[ord],
  check.names  = FALSE,
  row.names    = NULL,
  stringsAsFactors = FALSE
)

# B. analysis

dat2 <- syn_data %>% filter(!is.na(nihss_adm), !is.na(time_difference))

m_lm <- lm(time_difference ~ nihss_adm, data = dat2)

# Robust vcov (HC3 is a good default with outliers)
V.rob <- vcovHC(m_lm, type = "HC3")

# Robust coefficient table
rob_tab <- coeftest(m_lm, vcov. = V.rob)
rob_tab

# extract slope row only
est  <- rob_tab["nihss_adm","Estimate"]
se   <- rob_tab["nihss_adm","Std. Error"]
pval <- rob_tab["nihss_adm","Pr(>|t|)"]

# CI
CI   <- est + c(-1,1)*1.96*se

# format nicely
slope_CI <- sprintf("%.4f (%.4f–%.4f)", est, CI[1], CI[2])
p_fmt    <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))

# final one-row table
tbl2 <- data.frame(
  Slope = slope_CI,
  `p-value` = p_fmt,
  check.names = FALSE,
  row.names = NULL
)

print(tbl2, row.names = FALSE)

syn_data_timeDiff<-syn_data
### linear model with robust SEs (due to not normal residuals and not Homoskedasticity)
#### Each additional NIHSS point is associated with, on average, 0.009 days later randomization.
save(tbl, tbl2, p, p2, file="Saved_Data/TimeDeviationvsNIHSS.RData")
