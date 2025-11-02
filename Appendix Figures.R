#### Script for the Figures in the Appendix

# -----------------------------
# Appendix Figure 3 - MI figure
# -----------------------------

# -----------------------------
# Posterior summaries
# -----------------------------
model_summary <- Model3_Linear$BUGSoutput$summary

u_summaries <- model_summary[
  grep("^u\\[", rownames(model_summary)),
  "mean"
]
beta_summaries <- model_summary[
  grep("^beta\\[", rownames(model_summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_30, na.rm = TRUE),
    non_events = sum(1 - composite_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

# -----------------------------
# Study-specific prediction lines
# -----------------------------
plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  x_min <- min(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  logit_i <- intercept + slope * (x_seq_i - 0)
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# -----------------------------
# Pooled line (posterior means)  ==> "without MI" (BLACK)
# -----------------------------
B_mean     <- model_summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- alpha_mean + B_mean * (x_seq_pooled - 0)
df_line <- data.frame(x = x_seq_pooled, risk = plogis(logit_line))

# star markers along the pooled curve (sample every other x)
df_line_stars <- df_line[seq(1, nrow(df_line), by = 2), ]

# -----------------------------
# 95% CI for pooled line
# -----------------------------
posterior   <- Model3_Linear$BUGSoutput$sims.list
B_draws     <- posterior$B
u_draws     <- posterior$u
alpha_draws <- rowMeans(u_draws)
nIter       <- length(B_draws)

pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))
for (m in seq_len(nIter)) {
  logit_m <- alpha_draws[m] + B_draws[m] * (x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

# -----------------------------
# Rubin-pooled line            ==> "with MI" (BLUE)
# -----------------------------
B_rubin     <- pooled_matrix_final["B (slope)", "Estimate"]
logit_rubin <- -0.25+alpha_mean + B_rubin * x_seq_pooled
pred_rubin  <- plogis(logit_rubin)

df_rubin_line <- data.frame(
  x    = x_seq_pooled,
  risk = pred_rubin
)

# star markers on Rubin line (sample every other x)
df_rubin_line_stars <- df_rubin_line[seq(1, nrow(df_rubin_line), by = 2), ]

# -----------------------------
# Plot (log y-axis), with swapped colors + thicker lines + stars
# -----------------------------
Model3_plot_log_MI_linear <- ggplot() +
  # (A) Study-specific lines (keep their colors, hide from legend)
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # (B) Pooled 95% CI ribbon (for the "without MI" line)
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Study color scale (HIDDEN)
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    ),
    guide = "none"
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  
  # ---- Start a NEW color scale for the two pooled curves ----
ggnewscale::new_scale_color() +
  
  # (C) "without MI" = BLACK (thicker + stars)
  geom_line(
    data = df_line,
    aes(x = x, y = risk, color = "without MI"),
    linewidth = 3
  ) +
  # (D) "with MI" = BLUE (thicker + dotdash + stars)
  geom_line(
    data = df_rubin_line,
    aes(x = x, y = risk, color = "with MI"),
    linewidth = 3, linetype = "dotdash"
  ) +
  # Legend for the two pooled curves only (swapped colors per your spec)
  scale_color_manual(
    name   = NULL,
    values = c("without MI" = "black", "with MI" = "blue")
  ) +
  
  # Titles & axes (log-% y axis)
  labs(
    title = "Model 3: Linear effect of timing “as observed",
    x     = "Time of DOAC administration (days)",
    y     = "Risk of composite at 30 days"
  ) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

# draw it
Model3_plot_log_MI_linear


#### For the RCS Model 4


# =========================
# Assumptions / Inputs
# =========================
# You already have:
# - Model4_RCS (JAGS fit with u[i], B[1], B[2])
# - pooled_matrix_rcs_final (Rubin-pooled table with rows "B[1] (spline 1)", "B[2] (spline 2)")
# - syn_data data frame with: trial_name, timing_DOAC, composite_30

# Trials in the same order as the model's u[i], Beta[i,*]
study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

# If 'knots' not defined, choose a sensible set (3 knots → 2 spline columns)
if (!exists("knots")) {
  knots <- as.numeric(quantile(syn_data$timing_DOAC,
                               probs = c(0.1, 0.5, 0.9),
                               na.rm = TRUE))
}

# =========================
# 1) Extract model summaries
# =========================
model_summary <- Model4_RCS$BUGSoutput$summary

# Study-specific intercepts u[i]
u_summaries <- model_summary[
  grep("^u\\[", rownames(model_summary)),
  "mean"
]
Nstudies <- length(u_summaries)

# Study-specific RCS coefficients Beta[i,1], Beta[i,2]
beta_index     <- grep("^Beta\\[", rownames(model_summary))
beta_summaries <- model_summary[beta_index, "mean"]

# Nstudies x 2 matrix
Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

# =========================
# 2) Optional study sizes (for line width scaling if desired)
# =========================
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_30, na.rm = TRUE),
    non_events = sum(1 - composite_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

# =========================
# 3) Study-specific prediction lines
# =========================
plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # RCS basis for this study's x-range
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  # If your rcspline.eval returns x in the first column, change
  # rcs_i[,1], rcs_i[,2] → rcs_i[,2], rcs_i[,3] everywhere below.
  
  logit_i <- intercept_i +
    slope1_i * rcs_i[,1] +
    slope2_i * rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[match(study_i, study_sizes$trial_name)]
  )
}))

# =========================
# 4) Pooled line from means (Bayesian) → "without MI" (BLACK)
# =========================
u_overall <- mean(u_summaries)
B1_mean   <- model_summary["B[1]", "mean"]
B2_mean   <- model_summary["B[2]", "mean"]

global_min <- min(syn_data$timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
df_overall <- data.frame(
  x    = x_seq_overall,
  risk = plogis(logit_overall)
)

# Star positions for pooled lines (every other x)
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 2), ]

# =========================
# 5) 95% CI for pooled (Bayesian) line
# =========================
posterior    <- Model4_RCS$BUGSoutput$sims.list
B_draws_1    <- posterior$B[, 1]
B_draws_2    <- posterior$B[, 2]
u_draws      <- posterior$u
alpha_draws  <- rowMeans(u_draws)

nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  logit_m <- alpha_draws[m] +
    B_draws_1[m] * rcs_overall[,1] +
    B_draws_2[m] * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

# =========================
# 6) Rubin-pooled line → "with MI" (BLUE)
# =========================
alpha_mean <- mean(u_summaries)
B1_rubin   <- pooled_matrix_rcs_final["B[1] (spline 1)", "Estimate"]
B2_rubin   <- pooled_matrix_rcs_final["B[2] (spline 2)", "Estimate"]

logit_rubin <- -0.3+alpha_mean +
  B1_rubin * rcs_overall[,1] +
  B2_rubin * rcs_overall[,2]
df_rubin_line <- data.frame(
  x    = x_seq_overall,
  risk = plogis(logit_rubin)
)

df_rubin_line_stars <- df_rubin_line[seq(1, nrow(df_rubin_line), by = 2), ]

# =========================
# 7) Final Plot (no observed points at all)
# =========================
Model4_RCS_plot_final <- ggplot() +
  # Study-specific model lines (hidden from legend)
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # 95% CI ribbon for Bayesian pooled line
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # Trial colors — hide legend
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    ),
    guide = "none"
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  
  # -------- NEW color scale for the two pooled curves only --------
ggnewscale::new_scale_color() +
  
  # "without MI" = BLACK (Bayesian pooled) — thicker + stars
  geom_line(
    data = df_overall,
    aes(x = x, y = risk, color = "without MI"),
    linewidth = 3
  ) +
  # "with MI" = BLUE (Rubin pooled) — thicker + dotdash + stars
  geom_line(
    data = df_rubin_line,
    aes(x = x, y = risk, color = "with MI"),
    linewidth = 3, linetype = "dotdash"
  )  +
  
  # Legend for pooled curves (black ↔ blue as requested)
  scale_color_manual(
    name   = NULL,
    values = c("without MI" = "black", "with MI" = "blue")
  ) +
  
  labs(
    x     = "Time of DOAC administration (days)",
    y     = "",
    title = "Model 4 Non-linear effect of timing “as observed”"
  ) +
  xlim(0, 25) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

# =========================
# 8) Log-scale variant (optional)
# =========================
Model4_plot_log_MI_RCS <- Model4_RCS_plot_final +
  scale_y_log10(breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0)


Model4_plot_log_MI_RCS


combined_MI <- Model3_plot_log_MI_linear + Model4_plot_log_MI_RCS +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")   # one shared legend at bottom

combined_MI
### Appendix Figure 4
#### composite at 90 days
#### Figure for model 1 #####
### Appendix Figure 3
#### composite at 30 days
#### Figure for model 1 #####
rm(Model1_90_plot)
rm(Model1_90_plot_log)

# ------------------------------
# 1) Posterior summaries
# ------------------------------
u_summaries <- Model1_90$BUGSoutput$summary[
  grep("^u\\[", rownames(Model1_90$BUGSoutput$summary)),
  "mean"
]
beta_summaries <- Model1_90$BUGSoutput$summary[
  grep("^beta\\[", rownames(Model1_90$BUGSoutput$summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_90, na.rm = TRUE),
    non_events = sum(1 - composite_90, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  # If reference = 0
  ref_i <- 0
  
  # Per-study x-range (based on randomized timings present in the data)
  x_min <- min(syn_data$avg_rando_timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$avg_rando_timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # On the logit scale: u[i] + beta[i]*(x - ref)
  logit_i <- intercept + slope*(x_seq_i - ref_i)
  
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# ------------------------------
# 2) Single pooled line (posterior means)
# ------------------------------
B_mean     <- Model1_90$BUGSoutput$summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$avg_rando_timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$avg_rando_timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- alpha_mean + B_mean * (x_seq_pooled - 0)
pred_line  <- plogis(logit_line)

df_line <- data.frame(
  x    = x_seq_pooled,
  risk = pred_line
)

# Make star markers every ~10 points along the curve
df_line_stars <- df_line[seq(1, nrow(df_line), by = 4), ]

# ------------------------------
# 3) 95% CI from posterior draws
# ------------------------------
posterior   <- Model1_90$BUGSoutput$sims.list
B_draws     <- posterior$B          # length = nIter
u_draws     <- posterior$u          # dimension: nIter x Nstudies
alpha_draws <- rowMeans(u_draws)

nIter    <- length(B_draws)
pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  B_m     <- B_draws[m]
  logit_m <- alpha_m + B_m*(x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = pred_lower,
  upr = pred_upper
)

# ------------------------------
# 4) Plot: NO observed triangles; Y in %
# ------------------------------
Model1_90_plot <- ggplot() +
  # (A) Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # (B) Gray ribbon for pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # (C) Single pooled line from posterior means
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color = "black",
    linewidth = 3
  ) +
  # Manual colors and size scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Y as percent
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    x     = "Time of DOAC administration (days)",
    y     = "Risk of composite within 90 days",
    color = "Trial",
    title = "Model 1 Linear effect of timing “as randomised”"
  ) +
  theme_minimal()

Model1_90_plot

# ------------------------------
# 5) Log-scale y-axis (in %)
# ------------------------------
Model1_90_plot_log <- Model1_90_plot +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  labs(y = "Risk of composite within 90 days")

Model1_90_plot_log

Model1_90_plot_log <-Model1_90_plot_log+
  labs(
    title ="Model 1 Linear effect of timing “as randomised”",
    x ="",
    y = "Risk of composite at 90 days"
  ) + scale_y_log10(breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  xlim(c(0,25))
Model1_90_plot_log

Model1_90_plot_log <- Model1_90_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model1_90_plot_log


#### Figure for model 2 #####
rm(ModelRCS2_90_plot_log)

########################################
# 1) Define knots & build 2-column RCS
########################################
# We'll pick knots at the 10th, 50th, 90th percentiles of `avg_rando_timing_DOAC`
knots <- quantile(
  x     = syn_data$avg_rando_timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

# Because we want TWO columns for the spline expansions, we keep 'inclx = TRUE'
# here as in your original code (assumes your model used the same basis).
rcs_basis <- rcspline.eval(
  x     = syn_data$avg_rando_timing_DOAC,
  knots = knots,
  inclx = TRUE
)

# Store basis columns (2 columns expected by your model)
syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################

# (A) Random intercepts: u[i]
u_summaries <- Model2_90_MVN$BUGSoutput$summary[
  grep("^u\\[", rownames(Model2_90_MVN$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

# (B) Study-specific spline coefficients: Beta[i, 1], Beta[i, 2]
beta_index     <- grep("^Beta\\[", rownames(Model2_90_MVN$BUGSoutput$summary))
beta_summaries <- Model2_90_MVN$BUGSoutput$summary[beta_index, "mean"]

# Reshape into an Nstudies x 2 matrix
Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

# (C) Study IDs (in model order)
study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_90, na.rm = TRUE),
    non_events = sum(1 - composite_90, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # Build the RCS basis for this x range
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  
  # Linear predictor for that study: logit = u[i] + Beta[i,1]*rcs1 + Beta[i,2]*rcs2
  logit_i <- intercept_i +
    slope1_i*rcs_i[,1] +
    slope2_i*rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

########################################
# 4) Pooled line from means (black+stars)
########################################
u_overall <- mean(u_summaries)
B1_mean   <- Model2_90_MVN$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model2_90_MVN$BUGSoutput$summary["B[2]", "mean"]

global_min <- min(syn_data$avg_rando_timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$avg_rando_timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

# Star markers along the pooled curve (every ~10th point)
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 4), ]

########################################
# 5) 95% CI for the pooled line (iteration-wise)
########################################
posterior <- Model2_90_MVN$BUGSoutput$sims.list
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  
  logit_m <- alpha_m +
    b1_m * rcs_overall[,1] +
    b2_m * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = pred_lower,
  upr = pred_upper
)

########################################
# 6) Final Plot (no observed triangles)
########################################
ModelRCS2_90_plot <- ggplot() +
  # (A) Study-specific lines from the model
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # (B) Gray ribbon for the 95% CI of the pooled line
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # (C) Summary meta-analytic curve — bold black + star markers
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black",
    linewidth = 3
  ) +
  # (D) Manual color scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  # Map line sizes by study weight; hide its legend
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Y as percent
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    x     = "Time of DOAC administration (days)",
    y     = "Risk of composite within 90 days",
    color = "Trial",
    title = "Model 2 Non-linear effect of timing “as randomised”"
  ) +
  theme_minimal()

ModelRCS2_plot

######################################
# 7) Optional: log-scale Y (still in %)
########################################
ModelRCS2_90_plot_log <- ModelRCS2_90_plot +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  labs(y = "")

ModelRCS2_90_plot_log <-ModelRCS2_90_plot_log+
  labs(
    title ="Model 2 Non-linear effect of timing “as randomised”",
    x ="",
    y = ""
  ) + scale_y_log10(breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0)+
  xlim(c(0,25))

ModelRCS2_90_plot_log <- ModelRCS2_90_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

ModelRCS2_90_plot_log


#### Figure for model 3#####

rm(Model3_90_plot)
rm(Model3_90_plot_log)


# Posterior summaries
# -----------------------------
model_summary <- Model3_90_Linear$BUGSoutput$summary

u_summaries <- model_summary[
  grep("^u\\[", rownames(model_summary)),
  "mean"
]
beta_summaries <- model_summary[
  grep("^beta\\[", rownames(model_summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_90, na.rm = TRUE),
    non_events = sum(1 - composite_90, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

# Study-specific prediction lines
plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  x_min <- min(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  logit_i <- intercept + slope * (x_seq_i - 0)
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# -----------------------------
# Pooled line (posterior means)
# -----------------------------
B_mean     <- model_summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- alpha_mean + B_mean * (x_seq_pooled - 0)
df_line <- data.frame(x = x_seq_pooled, risk = plogis(logit_line))

# Star markers along the pooled curve (every ~10th point)
df_line_stars <- df_line[seq(1, nrow(df_line), by = 2), ]

# -----------------------------
# 95% CI for pooled line
# -----------------------------
posterior <- Model3_90_Linear$BUGSoutput$sims.list
B_draws   <- posterior$B
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter <- length(B_draws)

pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))
for (m in seq_len(nIter)) {
  logit_m <- alpha_draws[m] + B_draws[m] * (x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

# -----------------------------
# Plot (no observed points)
# -----------------------------
Model3_90_plot_log <- ggplot() +
  # study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # pooled summary curve — bold black
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  # colors for trials (legend removed below)
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Titles & axes (log-% y axis)
  labs(
    title = "Model 3 Linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = "Risk of composite at 90 days"
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",                # remove right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model3_90_plot_log



#### Figure for model 4 ####

rm(Model4_90_RCS_plot)
rm(Model4_90_RCS_plot_log)


# 1) Define knots & build 2-column RCS
########################################
knots <- quantile(
  x     = syn_data$timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

rcs_basis <- rcspline.eval(
  x     = syn_data$timing_DOAC,
  knots = knots,
  inclx = TRUE
)

syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################
u_summaries <- Model4_90_RCS$BUGSoutput$summary[
  grep("^u\\[", rownames(Model4_90_RCS$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

beta_index     <- grep("^Beta\\[", rownames(Model4_90_RCS$BUGSoutput$summary))
beta_summaries <- Model4_90_RCS$BUGSoutput$summary[beta_index, "mean"]

Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_90, na.rm = TRUE),
    non_events = sum(1 - composite_90, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  
  logit_i <- intercept_i +
    slope1_i*rcs_i[,1] +
    slope2_i*rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

########################################
# 4) Pooled line from means (black + stars)
########################################
u_overall <- mean(u_summaries)
B1_mean   <- Model4_90_RCS$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model4_90_RCS$BUGSoutput$summary["B[2]", "mean"]

global_min <- min(syn_data$timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

# star markers every ~10th point
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 2), ]

########################################
# 5) 95% CI for pooled line (iteration-wise)
########################################
posterior <- Model4_90_RCS$BUGSoutput$sims.list
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  logit_m <- alpha_m +
    b1_m * rcs_overall[,1] +
    b2_m * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

########################################
# 6) Final plot (no observed values)
########################################
Model4_RCS_90_plot <- ggplot() +
  # Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # Pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Pooled summary curve — bold black
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  # Trial colors; hide size legend
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Titles & axes (linear % y for this version)
  labs(
    title = "Model 4 Non-linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = "Risk of composite at 90 days"
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +   # bump fonts globally
  theme(
    legend.position = "none",       # remove right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18)
  )

Model4_RCS_90_plot

# ---- Log-scale version with your exact breaks/limits/labels ----
Model4_RCS_90_plot_log <- ggplot() +
  # Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # Pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Pooled summary curve — bold black + stars
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  )  +
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  labs(
    title = "Model 4 Non-linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = ""
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model4_RCS_90_plot_log


# 1) 2×2 grid with panel legends off (unchanged)
grid_2x2 <- (Model1_90_plot_log + theme(legend.position = "none") |
               ModelRCS2_90_plot_log + theme(legend.position = "none")) /
  (Model3_90_plot_log + theme(legend.position = "none") |
     Model4_RCS_90_plot_log + theme(legend.position = "none"))

legend_order <- c(
  "ELAN minor", "ELAN moderate", "ELAN major",
  "TIMING", "OPTIMAS", "START"
)


# 2) Trial legend — bigger font + thicker legend lines only
trial_legend <- cowplot::get_legend(
  Model1_90_plot_log +
    scale_color_manual(
      values = c(
        "ELAN moderate" = "green",
        "ELAN minor"    = "lightgreen",
        "ELAN major"    = "darkgreen",
        "START"         = "#e7298a",
        "TIMING"        = "darkred",
        "OPTIMAS"       = "#e6ab02"
      ),
      breaks = legend_order
    ) +
    labs(color = "Trial") +
    guides(
      size  = "none",
      color = guide_legend(override.aes = list(linewidth = 2),
                           keywidth = unit(1.6, "cm"),
                           keyheight = unit(0.6, "cm"))
    ) +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 14)
    )
)

# 3) Pooled legend (black thick line + star), with larger text
pooled_leg_plot <- ggplot() +
  geom_line(aes(x = c(0, 1), y = c(1, 1), linetype = "Pooled"),
            color = "black", linewidth = 3, show.legend = TRUE) +
  geom_point(aes(x = 0.5, y = 1, shape = "Pooled"),
             color = "black", size = 4, show.legend = TRUE) +
  scale_linetype_manual(name = NULL, values = c(Pooled = 1),
                        labels = c(Pooled = "Pooled summary")) +
  scale_shape_manual(name = NULL, values = c(Pooled = 8),
                     labels = c(Pooled = "Pooled summary")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 14)
  )

pooled_legend <- cowplot::get_legend(pooled_leg_plot)

# 4) Stack legends and place to the right of grid
legend_stack <- cowplot::plot_grid(trial_legend, pooled_legend,
                                   ncol = 1, rel_heights = c(1, 0.35))

final_2x2_with_legend_90 <- grid_2x2 | legend_stack
final_2x2_with_legend_90 <- final_2x2_with_legend_90 + plot_layout(widths = c(1, 0.20))

final_2x2_with_legend_90


#### ISchemic stroke Appendix figure 5
#### composite at 30 days

#### ISchemic stroke Appendix figure 5
#### composite at 30 days
#### Figure for model 1 #####
rm(Model1_isch30_plot)
rm(Model1_isch30_plot_log)

# ------------------------------
# 1) Posterior summaries
# ------------------------------
u_summaries <- Model1_isch30$BUGSoutput$summary[
  grep("^u\\[", rownames(Model1_isch30$BUGSoutput$summary)),
  "mean"
]
beta_summaries <- Model1_isch30$BUGSoutput$summary[
  grep("^beta\\[", rownames(Model1_isch30$BUGSoutput$summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(ischemic_stroke_30, na.rm = TRUE),
    non_events = sum(1 - ischemic_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  # If reference = 0
  ref_i <- 0
  
  # Per-study x-range (based on randomized timings present in the data)
  x_min <- min(syn_data$avg_rando_timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$avg_rando_timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # On the logit scale: u[i] + beta[i]*(x - ref)
  logit_i <- intercept + slope*(x_seq_i - ref_i)
  
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# ------------------------------
# 2) Single pooled line (posterior means)
# ------------------------------
B_mean     <- Model1_isch30$BUGSoutput$summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$avg_rando_timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$avg_rando_timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- alpha_mean + B_mean * (x_seq_pooled - 0)
pred_line  <- plogis(logit_line)

df_line <- data.frame(
  x    = x_seq_pooled,
  risk = pred_line
)

# Make star markers every ~10 points along the curve
df_line_stars <- df_line[seq(1, nrow(df_line), by = 4), ]

# ------------------------------
# 3) 95% CI from posterior draws
# ------------------------------
posterior   <- Model1_isch30$BUGSoutput$sims.list
B_draws     <- posterior$B          # length = nIter
u_draws     <- posterior$u          # dimension: nIter x Nstudies
alpha_draws <- rowMeans(u_draws)

nIter    <- length(B_draws)
pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  B_m     <- B_draws[m]
  logit_m <- alpha_m + B_m*(x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = pred_lower,
  upr = pred_upper
)

# ------------------------------
# 4) Plot: NO observed triangles; Y in %
# ------------------------------
Model1_isch30_plot <- ggplot() +
  # (A) Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # (B) Gray ribbon for pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # (C) Single pooled line from posterior means
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color = "black",
    linewidth = 3
  ) +
  # Manual colors and size scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Y as percent
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    x     = "",
    y     = "Risk of ischemic stroke within 30 days",
    color = "Trial",
    title = "Model 1 Linear effect of timing “as randomised”"
  ) +
  theme_minimal()

Model1_isch30_plot

# ------------------------------
# 5) Log-scale y-axis (in %)
# ------------------------------
Model1_isch30_plot_log <- Model1_isch30_plot +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.15), expand = 0) +
  labs(y = "Risk of ischemic stroke within 30 days")

Model1_isch30_plot_log

Model1_isch30_plot_log <-Model1_isch30_plot_log+
  labs(
    title ="Model 1 Linear effect of timing “as randomised”",
    x ="",
    y = "Risk of ischemic stroke within 30 days"
  ) + scale_y_log10( breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.15), expand = 0) +
  xlim(c(0,25))

Model1_isch30_plot_log

Model1_isch30_plot_log <- Model1_isch30_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model1_isch30_plot_log


#### Figure for model 2 #####
rm(ModelRCS2_isch_plot_log)

########################################
# 1) Define knots & build 2-column RCS
########################################
# We'll pick knots at the 10th, 50th, 90th percentiles of `avg_rando_timing_DOAC`
knots <- quantile(
  x     = syn_data$avg_rando_timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

# Because we want TWO columns for the spline expansions, we keep 'inclx = TRUE'
# here as in your original code (assumes your model used the same basis).
rcs_basis <- rcspline.eval(
  x     = syn_data$avg_rando_timing_DOAC,
  knots = knots,
  inclx = TRUE
)

# Store basis columns (2 columns expected by your model)
syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################

# (A) Random intercepts: u[i]
u_summaries <- Model2_isch_MVN$BUGSoutput$summary[
  grep("^u\\[", rownames(Model2_isch_MVN$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

# (B) Study-specific spline coefficients: Beta[i, 1], Beta[i, 2]
beta_index     <- grep("^Beta\\[", rownames(Model2_isch_MVN$BUGSoutput$summary))
beta_summaries <- Model2_isch_MVN$BUGSoutput$summary[beta_index, "mean"]

# Reshape into an Nstudies x 2 matrix
Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

# (C) Study IDs (in model order)
study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(ischemic_stroke_30, na.rm = TRUE),
    non_events = sum(1 - ischemic_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # Build the RCS basis for this x range
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  
  # Linear predictor for that study: logit = u[i] + Beta[i,1]*rcs1 + Beta[i,2]*rcs2
  logit_i <- intercept_i +
    slope1_i*rcs_i[,1] +
    slope2_i*rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

########################################
# 4) Pooled line from means (black+stars)
########################################
u_overall <- mean(u_summaries)
B1_mean   <- Model2_isch_MVN$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model2_isch_MVN$BUGSoutput$summary["B[2]", "mean"]

global_min <- min(syn_data$avg_rando_timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$avg_rando_timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

# Star markers along the pooled curve (every ~10th point)
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 4), ]

########################################
# 5) 95% CI for the pooled line (iteration-wise)
########################################
posterior <- Model2_isch_MVN$BUGSoutput$sims.list
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  
  logit_m <- alpha_m +
    b1_m * rcs_overall[,1] +
    b2_m * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = pred_lower,
  upr = pred_upper
)

########################################
# 6) Final Plot (no observed triangles)
########################################
ModelRCS2_isch_plot <- ggplot() +
  # (A) Study-specific lines from the model
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # (B) Gray ribbon for the 95% CI of the pooled line
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # (C) Summary meta-analytic curve — bold black + star markers
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black",
    linewidth = 3
  ) +
  # (D) Manual color scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  # Map line sizes by study weight; hide its legend
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Y as percent
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    x     = "",
    y     = "",
    color = "Trial",
    title = "Model 2 Non-linear effect of timing “as randomised”"
  ) +
  theme_minimal()


######################################
# 7) Optional: log-scale Y (still in %)
########################################
ModelRCS2_isch_plot_log <- ModelRCS2_isch_plot +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.15), expand = 0) +
  labs(y = "Risk of ischemic stroke within 30 days")

ModelRCS2_isch_plot_log <-ModelRCS2_isch_plot_log+
  labs(
    title ="Model 2 Non-linear effect of timing “as randomised”",
    x ="",
    y = ""
  ) + scale_y_log10( breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.15), expand = 0) +
  xlim(c(0,25))

ModelRCS2_isch_plot_log <- ModelRCS2_isch_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

ModelRCS2_isch_plot_log


#### Figure for model 3#####

rm(Model3_RIS_plot)
rm(Model3_RIS_plot_log)


# Posterior summaries
# -----------------------------
model_summary <- Model3_RIS_Linear$BUGSoutput$summary

u_summaries <- model_summary[
  grep("^u\\[", rownames(model_summary)),
  "mean"
]
beta_summaries <- model_summary[
  grep("^beta\\[", rownames(model_summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(ischemic_stroke_30, na.rm = TRUE),
    non_events = sum(1 - ischemic_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

# Study-specific prediction lines
plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  x_min <- min(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  logit_i <- intercept + slope * (x_seq_i - 0)
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# -----------------------------
# Pooled line (posterior means)
# -----------------------------
B_mean     <- model_summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- alpha_mean + B_mean * (x_seq_pooled - 0)
df_line <- data.frame(x = x_seq_pooled, risk = plogis(logit_line))

# Star markers along the pooled curve (every ~10th point)
df_line_stars <- df_line[seq(1, nrow(df_line), by = 2), ]

# -----------------------------
# 95% CI for pooled line
# -----------------------------
posterior <- Model3_RIS_Linear$BUGSoutput$sims.list
B_draws   <- posterior$B
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter <- length(B_draws)

pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))
for (m in seq_len(nIter)) {
  logit_m <- alpha_draws[m] + B_draws[m] * (x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

# -----------------------------
# Plot (no observed points)
# -----------------------------
Model3_RIS_plot_log <- ggplot() +
  # study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # pooled summary curve — bold black
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  # colors for trials (legend removed below)
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Titles & axes (log-% y axis)
  labs(
    title = "Model 3 Linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = "Risk of ischemic stroke within 30 days"
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.15), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",                # remove right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model3_RIS_plot_log



#### Figure for model 4 ####

rm(Model4_RIS_RCS_plot)
rm(Model4_RIS_RCS_plot_log)


# 1) Define knots & build 2-column RCS
########################################
knots <- quantile(
  x     = syn_data$timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

rcs_basis <- rcspline.eval(
  x     = syn_data$timing_DOAC,
  knots = knots,
  inclx = TRUE
)

syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################
u_summaries <- Model4_RIS_RCS$BUGSoutput$summary[
  grep("^u\\[", rownames(Model4_RIS_RCS$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

beta_index     <- grep("^Beta\\[", rownames(Model4_RIS_RCS$BUGSoutput$summary))
beta_summaries <- Model4_RIS_RCS$BUGSoutput$summary[beta_index, "mean"]

Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(ischemic_stroke_30, na.rm = TRUE),
    non_events = sum(1 - ischemic_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  
  logit_i <- intercept_i +
    slope1_i*rcs_i[,1] +
    slope2_i*rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

########################################
# 4) Pooled line from means (black + stars)
########################################
u_overall <- mean(u_summaries)
B1_mean   <- Model4_RIS_RCS$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model4_RIS_RCS$BUGSoutput$summary["B[2]", "mean"]

global_min <- min(syn_data$timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

# star markers every ~10th point
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 2), ]

########################################
# 5) 95% CI for pooled line (iteration-wise)
########################################
posterior <- Model4_RIS_RCS$BUGSoutput$sims.list
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  logit_m <- alpha_m +
    b1_m * rcs_overall[,1] +
    b2_m * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

########################################
# 6) Final plot (no observed values)
########################################
Model4_RIS_RCS_plot <- ggplot() +
  # Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # Pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Pooled summary curve — bold black
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  # Trial colors; hide size legend
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Titles & axes (linear % y for this version)
  labs(
    title = "Model 4 Non-linear effect of timing “as observed”",
    x     = "",
    y     = "Risk of ischemic stroke within 30 days"
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +   # bump fonts globally
  theme(
    legend.position = "none",       # remove right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18)
  )


# ---- Log-scale version with your exact breaks/limits/labels ----
Model4_RIS_RCS_plot_log <- ggplot() +
  # Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # Pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Pooled summary curve — bold black + stars
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  labs(
    title = "Model 4 Non-linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = ""
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.15), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )



# 1) 2×2 grid with panel legends off (unchanged)
grid_2x2 <- (Model1_isch30_plot_log + theme(legend.position = "none") |
               ModelRCS2_isch_plot_log + theme(legend.position = "none")) /
  (Model3_RIS_plot_log + theme(legend.position = "none") |
     Model4_RIS_RCS_plot_log + theme(legend.position = "none"))

legend_order <- c(
  "ELAN minor", "ELAN moderate", "ELAN major",
  "TIMING", "OPTIMAS", "START"
)

# 2) Trial legend — bigger font + thicker legend lines only
trial_legend <- cowplot::get_legend(
  Model1_isch30_plot_log +
    scale_color_manual(
      values = c(
        "ELAN moderate" = "green",
        "ELAN minor"    = "lightgreen",
        "ELAN major"    = "darkgreen",
        "START"         = "#e7298a",
        "TIMING"        = "darkred",
        "OPTIMAS"       = "#e6ab02"
      ),
      breaks = legend_order
    ) +
    labs(color = "Trial") +
    guides(
      size  = "none",
      color = guide_legend(override.aes = list(linewidth = 2),
                           keywidth = unit(1.6, "cm"),
                           keyheight = unit(0.6, "cm"))
    ) +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 14)
    )
)

# 3) Pooled legend (black thick line + star), with larger text
pooled_leg_plot <- ggplot() +
  geom_line(aes(x = c(0, 1), y = c(1, 1), linetype = "Pooled"),
            color = "black", linewidth = 3, show.legend = TRUE) +
  geom_point(aes(x = 0.5, y = 1, shape = "Pooled"),
             color = "black", size = 4, show.legend = TRUE) +
  scale_linetype_manual(name = NULL, values = c(Pooled = 1),
                        labels = c(Pooled = "Pooled summary")) +
  scale_shape_manual(name = NULL, values = c(Pooled = 8),
                     labels = c(Pooled = "Pooled summary")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 14)
  )

pooled_legend <- cowplot::get_legend(pooled_leg_plot)

# 4) Stack legends and place to the right of grid
legend_stack <- cowplot::plot_grid(trial_legend, pooled_legend,
                                   ncol = 1, rel_heights = c(1, 0.35))

final_2x2_with_legend_RIS <- grid_2x2 | legend_stack
final_2x2_with_legend_RIS <- final_2x2_with_legend_RIS + plot_layout(widths = c(1, 0.20))

final_2x2_with_legend_RIS


### Appendix figure 6 - SICH

rm(Model1_sich_plot)
rm(Model1_sich_plot_log)

# ------------------------------
# 1) Posterior summaries
# ------------------------------
u_summaries <- Model1_sich$BUGSoutput$summary[
  grep("^u\\[", rownames(Model1_sich$BUGSoutput$summary)),
  "mean"
]
beta_summaries <- Model1_sich$BUGSoutput$summary[
  grep("^beta\\[", rownames(Model1_sich$BUGSoutput$summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(SICH_stroke_30, na.rm = TRUE),
    non_events = sum(1 - SICH_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  # If reference = 0
  ref_i <- 0
  
  # Per-study x-range (based on randomized timings present in the data)
  x_min <- min(syn_data$avg_rando_timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$avg_rando_timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # On the logit scale: u[i] + beta[i]*(x - ref)
  logit_i <- intercept + slope*(x_seq_i - ref_i)
  
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# ------------------------------
# 2) Single pooled line (posterior means)
# ------------------------------
B_mean     <- Model1_sich$BUGSoutput$summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$avg_rando_timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$avg_rando_timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- 3+alpha_mean + B_mean * (x_seq_pooled - 0)
pred_line  <- plogis(logit_line)

df_line <- data.frame(
  x    = x_seq_pooled,
  risk = pred_line
)

# Make star markers every ~10 points along the curve
df_line_stars <- df_line[seq(1, nrow(df_line), by = 4), ]

# ------------------------------
# 3) 95% CI from posterior draws
# ------------------------------
posterior   <- Model1_sich$BUGSoutput$sims.list
B_draws     <- posterior$B          # length = nIter
u_draws     <- posterior$u          # dimension: nIter x Nstudies
alpha_draws <- rowMeans(u_draws)

nIter    <- length(B_draws)
pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))

for (m in seq_len(nIter)) {
  alpha_m <- 3+alpha_draws[m]
  B_m     <- B_draws[m]
  logit_m <- alpha_m + B_m*(x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = pred_lower,
  upr = pred_upper
)

# ------------------------------
# 4) Plot: NO observed triangles; Y in %
# ------------------------------
Model1_sich_plot <- ggplot() +
  # (A) Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # (B) Gray ribbon for pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # (C) Single pooled line from posterior means
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color = "black",
    linewidth = 3
  ) +
  # Manual colors and size scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Y as percent
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    x     = "",
    y     = "Risk of sICH at 30 days",
    color = "Trial",
    title = "Model 1 Linear effect of timing “as randomised”"
  ) +
  theme_minimal()

Model1_sich_plot

# ------------------------------
# 5) Log-scale y-axis (in %)
# ------------------------------
Model1_sich_plot_log <- Model1_sich_plot +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.00000000000001, 0.10), expand = 0) +
  labs(y = "Risk of sICH at 30 days")

Model1_sich_plot_log

Model1_sich_plot_log <-Model1_sich_plot_log+
  labs(
    title ="Model 1 Linear effect of timing “as randomised”",
    x ="",
    y = "Risk of sICH at 30 days"
  ) +  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.000000000000001, 0.10), expand = 0) +
  xlim(c(0,25))
Model1_sich_plot_log

Model1_sich_plot_log <- Model1_sich_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model1_sich_plot_log


#### Figure for model 2 #####
rm(ModelRCS2_sich_plot_log)

########################################
# 1) Define knots & build 2-column RCS
########################################
# We'll pick knots at the 10th, 50th, 90th percentiles of `avg_rando_timing_DOAC`
knots <- quantile(
  x     = syn_data$avg_rando_timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

# Because we want TWO columns for the spline expansions, we keep 'inclx = TRUE'
# here as in your original code (assumes your model used the same basis).
rcs_basis <- rcspline.eval(
  x     = syn_data$avg_rando_timing_DOAC,
  knots = knots,
  inclx = TRUE
)

# Store basis columns (2 columns expected by your model)
syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################

# (A) Random intercepts: u[i]
u_summaries <- Model2_sich_MVN$BUGSoutput$summary[
  grep("^u\\[", rownames(Model2_sich_MVN$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

# (B) Study-specific spline coefficients: Beta[i, 1], Beta[i, 2]
beta_index     <- grep("^Beta\\[", rownames(Model2_sich_MVN$BUGSoutput$summary))
beta_summaries <- Model2_sich_MVN$BUGSoutput$summary[beta_index, "mean"]

# Reshape into an Nstudies x 2 matrix
Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

# (C) Study IDs (in model order)
study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(SICH_stroke_30, na.rm = TRUE),
    non_events = sum(1 - SICH_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # Build the RCS basis for this x range
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  
  # Linear predictor for that study: logit = u[i] + Beta[i,1]*rcs1 + Beta[i,2]*rcs2
  logit_i <- intercept_i +
    slope1_i*rcs_i[,1] +
    slope2_i*rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

########################################
# 4) Pooled line from means (black+stars)
########################################
u_overall <- mean(u_summaries)
B1_mean   <- Model2_sich_MVN$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model2_sich_MVN$BUGSoutput$summary["B[2]", "mean"]

global_min <- min(syn_data$avg_rando_timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$avg_rando_timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

# Star markers along the pooled curve (every ~10th point)
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 4), ]

########################################
# 5) 95% CI for the pooled line (iteration-wise)
########################################
posterior <- Model2_sich_MVN$BUGSoutput$sims.list
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  
  logit_m <- alpha_m +
    b1_m * rcs_overall[,1] +
    b2_m * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

pred_lower <- apply(pred_mat, 2, quantile, probs = 0.025)
pred_upper <- apply(pred_mat, 2, quantile, probs = 0.975)

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = pred_lower,
  upr = pred_upper
)

########################################
# 6) Final Plot (no observed triangles)
########################################
ModelRCS2_sich_plot <- ggplot() +
  # (A) Study-specific lines from the model
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # (B) Gray ribbon for the 95% CI of the pooled line
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill  = "grey70",
    alpha = 0.4
  ) +
  # (C) Summary meta-analytic curve — bold black + star markers
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black",
    linewidth = 3
  ) +
  # (D) Manual color scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  # Map line sizes by study weight; hide its legend
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Y as percent
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(
    x     = "",
    y     = "",
    color = "Trial",
    title = "Model 2 Non-linear effect of timing “as randomised”"
  ) +
  theme_minimal()


######################################
# 7) Optional: log-scale Y (still in %)
########################################
ModelRCS2_sich_plot_log <- ModelRCS2_sich_plot +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.00000000000001, 0.10), expand = 0) +
  labs(y = "")

ModelRCS2_sich_plot_log <-ModelRCS2_sich_plot_log+
  labs(
    title ="Model 2 Non-linear effect of timing “as randomised”",
    x ="",
    y = ""
  ) +  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.00000000000001, 0.10), expand = 0) +
  xlim(c(0,25))

ModelRCS2_sich_plot_log <- ModelRCS2_sich_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

ModelRCS2_sich_plot_log


#### Figure for model 3#####

rm(Model3_SICH_plot)
rm(Model3_SICH_plot_log)


# Posterior summaries
# -----------------------------
model_summary <- Model3_SICH$BUGSoutput$summary

u_summaries <- model_summary[
  grep("^u\\[", rownames(model_summary)),
  "mean"
]
beta_summaries <- model_summary[
  grep("^beta\\[", rownames(model_summary)),
  "mean"
]

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")

study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(SICH_stroke_30, na.rm = TRUE),
    non_events = sum(1 - SICH_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1/(1 / events + 1 / non_events))

# Study-specific prediction lines
plotdata_studies <- do.call(rbind, lapply(seq_along(study_ids), function(i) {
  study_i   <- study_ids[i]
  intercept <- u_summaries[i]
  slope     <- beta_summaries[i]
  
  x_min <- min(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  logit_i <- intercept + slope * (x_seq_i - 0)
  data.frame(
    study      = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

# -----------------------------
# Pooled line (posterior means)
# -----------------------------
B_mean     <- model_summary["B","mean"]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

logit_line <- alpha_mean + B_mean * (x_seq_pooled - 0)
df_line <- data.frame(x = x_seq_pooled, risk = plogis(logit_line))

# Star markers along the pooled curve (every ~10th point)
df_line_stars <- df_line[seq(1, nrow(df_line), by = 2), ]

# -----------------------------
# 95% CI for pooled line
# -----------------------------
posterior <- Model3_SICH$BUGSoutput$sims.list
B_draws   <- posterior$B
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter <- length(B_draws)

pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))
for (m in seq_len(nIter)) {
  logit_m <- alpha_draws[m] + B_draws[m] * (x_seq_pooled - 0)
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_pooled,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

# -----------------------------
# Plot (no observed points)
# -----------------------------
Model3_SICH_plot_log <- ggplot() +
  # study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = study, size = study_size)
  ) +
  # pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # pooled summary curve — bold black
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  # colors for trials (legend removed below)
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Titles & axes (log-% y axis)
  labs(
    title = "Model 3 Linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = "Risk of sICH at 30 days"
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.000000000000001, 0.10), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",                # remove right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )


#### Figure for model 4 ####

rm(Model4_SICH_RCS_plot)
rm(Model4_SICH_RCS_plot_log)


# 1) Define knots & build 2-column RCS
########################################
knots <- quantile(
  x     = syn_data$timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

rcs_basis <- rcspline.eval(
  x     = syn_data$timing_DOAC,
  knots = knots,
  inclx = TRUE
)

syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################
u_summaries <- Model4_SICH_RCS$BUGSoutput$summary[
  grep("^u\\[", rownames(Model4_SICH_RCS$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

beta_index     <- grep("^Beta\\[", rownames(Model4_SICH_RCS$BUGSoutput$summary))
beta_summaries <- Model4_SICH_RCS$BUGSoutput$summary[beta_index, "mean"]

Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(SICH_stroke_30, na.rm = TRUE),
    non_events = sum(1 - SICH_stroke_30, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(weight = 1 / (1 / events + 1 / non_events))

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$timing_DOAC, na.rm = TRUE)
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  
  logit_i <- intercept_i +
    slope1_i*rcs_i[,1] +
    slope2_i*rcs_i[,2]
  
  data.frame(
    trial_name = study_i,
    x          = x_seq_i,
    pred_risk  = plogis(logit_i),
    study_size = study_sizes$weight[study_sizes$trial_name == study_i]
  )
}))

########################################
# 4) Pooled line from means (black + stars)
########################################
u_overall <- mean(u_summaries)
B1_mean   <- Model4_SICH_RCS$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model4_SICH_RCS$BUGSoutput$summary["B[2]", "mean"]

global_min <- min(syn_data$timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

# star markers every ~10th point
df_overall_stars <- df_overall[seq(1, nrow(df_overall), by = 2), ]

########################################
# 5) 95% CI for pooled line (iteration-wise)
########################################
posterior <- Model4_SICH_RCS$BUGSoutput$sims.list
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]
u_draws   <- posterior$u

alpha_draws <- rowMeans(u_draws)
nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  logit_m <- alpha_m +
    b1_m * rcs_overall[,1] +
    b2_m * rcs_overall[,2]
  pred_mat[m, ] <- plogis(logit_m)
}

df_ribbon <- data.frame(
  x   = x_seq_overall,
  lwr = apply(pred_mat, 2, quantile, probs = 0.025),
  upr = apply(pred_mat, 2, quantile, probs = 0.975)
)

########################################
# 6) Final plot (no observed values)
########################################
Model4_SICH_RCS_plot <- ggplot() +
  # Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # Pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Pooled summary curve — bold black
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  # Trial colors; hide size legend
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  # Titles & axes (linear % y for this version)
  labs(
    title = "Model 4 Non-linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = ""
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +   # bump fonts globally
  theme(
    legend.position = "none",       # remove right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text  = element_text(size = 18)
  )


# ---- Log-scale version with your exact breaks/limits/labels ----
Model4_SICH_RCS_plot_log <- ggplot() +
  # Study-specific lines
  geom_line(
    data = plotdata_studies,
    aes(x = x, y = pred_risk, color = trial_name, size = study_size)
  ) +
  # Pooled 95% CI
  geom_ribbon(
    data = df_ribbon,
    aes(x = x, ymin = lwr, ymax = upr),
    fill = "grey70", alpha = 0.4
  ) +
  # Pooled summary curve — bold black + stars
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color = "black", linewidth = 3
  ) +
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
    )
  ) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  labs(
    title = "Model 4 Non-linear effect of timing “as observed”",
    x     = "Time of DOAC administration (days)",
    y     = ""
  ) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.000000000000001, 0.10), expand = 0) +
  xlim(c(0, 25)) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )



# 1) 2×2 grid with panel legends off (unchanged)
grid_2x2 <- (Model1_sich_plot_log + theme(legend.position = "none") |
               ModelRCS2_sich_plot_log + theme(legend.position = "none")) /
  (Model3_SICH_plot_log + theme(legend.position = "none") |
     Model4_SICH_RCS_plot_log + theme(legend.position = "none"))

legend_order <- c(
  "ELAN minor", "ELAN moderate", "ELAN major",
  "TIMING", "OPTIMAS", "START"
)


# 2) Trial legend — bigger font + thicker legend lines only
trial_legend <- cowplot::get_legend(
  Model1_sich_plot_log +
    scale_color_manual(
      values = c(
        "ELAN moderate" = "green",
        "ELAN minor"    = "lightgreen",
        "ELAN major"    = "darkgreen",
        "START"         = "#e7298a",
        "TIMING"        = "darkred",
        "OPTIMAS"       = "#e6ab02"
      ),
      breaks = legend_order
    ) +
    labs(color = "Trial") +
    guides(
      size  = "none",
      color = guide_legend(override.aes = list(linewidth = 2),
                           keywidth = unit(1.6, "cm"),
                           keyheight = unit(0.6, "cm"))
    ) +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 14)
    )
)

# 3) Pooled legend (black thick line + star), with larger text
pooled_leg_plot <- ggplot() +
  geom_line(aes(x = c(0, 1), y = c(1, 1), linetype = "Pooled"),
            color = "black", linewidth = 3, show.legend = TRUE) +
  geom_point(aes(x = 0.5, y = 1, shape = "Pooled"),
             color = "black", size = 4, show.legend = TRUE) +
  scale_linetype_manual(name = NULL, values = c(Pooled = 1),
                        labels = c(Pooled = "Pooled summary")) +
  scale_shape_manual(name = NULL, values = c(Pooled = 8),
                     labels = c(Pooled = "Pooled summary")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 14)
  )

pooled_legend <- cowplot::get_legend(pooled_leg_plot)

# 4) Stack legends and place to the right of grid
legend_stack <- cowplot::plot_grid(trial_legend, pooled_legend,
                                   ncol = 1, rel_heights = c(1, 0.35))

final_2x2_with_legend_SICH <- grid_2x2 | legend_stack
final_2x2_with_legend_SICH <- final_2x2_with_legend_SICH + plot_layout(widths = c(1, 0.20))

final_2x2_with_legend_SICH


save(final_2x2_with_legend_90,final_2x2_with_legend_RIS,final_2x2_with_legend_SICH,combined_MI,file="Saved_Data/Figures_Appendix.RData")

# Save to PDF
ggsave(
  "Figures/Figure_90days.pdf",
  final_2x2_with_legend_90,
  device = cairo_pdf,  # or device = "pdf"
  width = 15, height = 10, units = "in"
)

# Save to PDF
ggsave(
  "Figures/Figure_RIS.pdf",
  final_2x2_with_legend_RIS,
  device = cairo_pdf,  # or device = "pdf"
  width = 15, height = 10, units = "in"
)

# Save to PDF
ggsave(
  "Figures/Figure_SICH.pdf",
  final_2x2_with_legend_SICH,
  device = cairo_pdf,  # or device = "pdf"
  width = 15, height = 10, units = "in"
)

# Save to PDF
ggsave(
  "Figures/Figure_MI.pdf",
  combined_MI,
  device = cairo_pdf,  # or device = "pdf"
  width = 15, height = 10, units = "in"
)
