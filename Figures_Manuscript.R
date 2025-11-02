################################################################################################################################
#########################################    Script generating       ##########################################################3
#########################################    FIGURES OF THE PAPER    ##########################################################3
################################################################################################################################


## ================================
## Data prep (expects `syn_data` already in memory)
## Required columns:
##   trial, avg_rando_timing_DOAC, timing_DOAC, composite_30
## ================================
syn_data <- syn_data %>%
  mutate(
    composite_30 = as.numeric(as.character(composite_30)),
    # arm id 1..14 (this is the stable thing we will color by)
    trial_arm    = interaction(trial, avg_rando_timing_DOAC, drop = TRUE),
    arm          = as.numeric(trial_arm),
    arm_id       = factor(arm, levels = 1:14)   # fill will use this
  )

## ---------------- Labels & facet order ----------------
new_labels <- c(
  "ELAN moderate, days \u2264 2",       # 1
  "ELAN minor, days \u2264 2",          # 2
  "TIMING, days \u2264 4",              # 3
  "OPTIMAS, days \u2264 4",             # 4
  "START, 2 < days \u2264 4",           # 5
  "ELAN minor, 2 < days \u2264 4",      # 6
  "START, 5 < days \u2264 6",           # 7
  "ELAN moderate, 5 < days \u2264 7",   # 8
  "ELAN major, 5 < days \u2264 7",      # 9
  "TIMING, 4 < days \u2264 10",         # 10
  "START, 9 < days \u2264 10",          # 11
  "OPTIMAS, 6 < days \u2264 14",        # 12
  "ELAN major, 11 < days \u2264 14",    # 13
  "START, 13 < days \u2264 14"          # 14
)

# y-axis labels shown in the plot (can be reordered later)
syn_data$arm_label <- factor(syn_data$arm, levels = 1:14, labels = new_labels)

# facet order
trial_order <- c("ELAN minor","ELAN moderate","ELAN major","TIMING","OPTIMAS","START")
if (!"trial_name" %in% names(syn_data)) syn_data$trial_name <- syn_data$trial
syn_data$trial_name <- factor(syn_data$trial_name, levels = trial_order)

## ---------------- Grey windows (row-bounded tiles) ----------------
windows_df <- data.frame(
  trial_name = c(
    "ELAN moderate","ELAN minor","TIMING","OPTIMAS","START","ELAN minor","START",
    "ELAN moderate","ELAN major","TIMING","START","OPTIMAS","ELAN major","START"
  ),
  arm_label = c(
    "ELAN moderate, days \u2264 2",
    "ELAN minor, days \u2264 2",
    "TIMING, days \u2264 4",
    "OPTIMAS, days \u2264 4",
    "START, 2 < days \u2264 4",
    "ELAN minor, 2 < days \u2264 4",
    "START, 5 < days \u2264 6",
    "ELAN moderate, 5 < days \u2264 7",
    "ELAN major, 5 < days \u2264 7",
    "TIMING, 4 < days \u2264 10",
    "START, 9 < days \u2264 10",
    "OPTIMAS, 6 < days \u2264 14",
    "ELAN major, 11 < days \u2264 14",
    "START, 13 < days \u2264 14"
  ),
  xmin = c(0,0,0,0,2,2,5,5,5,4,9,6,11,13),
  xmax = c(2,2,4,4,4,4,6,7,7,10,10,14,14,14)
) %>%
  mutate(
    xmid  = (xmin + xmax) / 2,
    width = pmax(xmax - xmin, 0.001)
  )

windows_df$trial_name <- factor(windows_df$trial_name, levels = trial_order)

## ================================
## Order rows within each trial (later window on top)
## ================================
display_levels <- windows_df %>%
  dplyr::group_by(trial_name) %>%
  dplyr::arrange(dplyr::desc(xmin), dplyr::desc(xmax), .by_group = TRUE) %>%
  dplyr::pull(arm_label) %>%
  unique()

syn_data$arm_label  <- factor(syn_data$arm_label,  levels = display_levels)
windows_df$arm_label <- factor(windows_df$arm_label, levels = display_levels)

## ================================
## Fixed colors by arm ID (1..14), not by text
## ================================
palette14 <- c(
  "#2ca25f", # 1  ELAN moderate, days ≤ 2
  "#2ca25f", # 2  ELAN minor, days ≤ 2
  "#2ca25f", # 3  TIMING, days ≤ 4
  "#2ca25f", # 4  OPTIMAS, days ≤ 4
  "#2ca25f", # 5  START, 2 < days ≤ 4
  "#d7301f", # 6  ELAN minor, 2 < days ≤ 4
  "#006d2c", # 7  START, 5 < days ≤ 6
  "#d7301f", # 8  ELAN moderate, 5 < days ≤ 7
  "#2ca25f", # 9  ELAN major, 5 < days ≤ 7
  "#d7301f", # 10 TIMING, 4 < days ≤ 10
  "#fcae91", # 11 START, 9 < days ≤ 10
  "#d7301f", # 12 OPTIMAS, 6 < days ≤ 14
  "#d7301f", # 13 ELAN major, 11 < days ≤ 14
  "#d7301f"  # 14 START, 13 < days ≤ 14
)
# Name palette by arm_id levels (1..14 as character)
palette_by_id <- setNames(palette14, levels(syn_data$arm_id))

## ================================
## Plot
## ================================
violin_plot <-
  ggplot(dplyr::filter(syn_data, !is.na(timing_DOAC)),
         aes(x = timing_DOAC, y = arm_label)) +
  # row-bounded grey windows
  geom_tile(
    data = windows_df,
    aes(x = xmid, y = arm_label, width = width, height = 1),
    inherit.aes = FALSE, fill = "grey70", alpha = 0.7
  ) +
  # violins colored by arm_id (NOT by arm_label text)
  geom_violin(aes(fill = arm_id), alpha = 0.7, trim = FALSE, scale = "width") +
  scale_fill_manual(values = palette_by_id, drop = FALSE, guide = "none") +
  facet_wrap(~ trial_name, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = "Timing DOAC (days)", y = NULL) +
  theme_minimal() +
  theme(
    strip.text.x     = element_blank(),
    strip.background = element_blank(),
    text             = element_text(size = 16),
    axis.text        = element_text(size = 16),
    axis.title       = element_text(size = 16)
  )

print(violin_plot)

# Save as a vector PDF (great for print)
ggsave(
  filename = "Figures/Figure1.pdf",
  plot     = violin_plot,
  device   = cairo_pdf,   # nicer font handling; needs Cairo installed
  width    = 12,          # inches
  height   = 9,
  units    = "in"
)

############################################# Figure 2 ##################################################


load("Saved_Data/Model1_Linear.RData") 
load("Saved_Data/Model2_rcs.RData") 
load("Saved_Data/Model3_Linear.RData")
load("Saved_Data/Model4_RCS.RData")


# ------------------------------
# 1) Posterior summaries
# ------------------------------
u_summaries <- Model1$BUGSoutput$summary[
  grep("^u\\[", rownames(Model1$BUGSoutput$summary)),
  "mean"
]
beta_summaries <- Model1$BUGSoutput$summary[
  grep("^beta\\[", rownames(Model1$BUGSoutput$summary)),
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
B_mean     <- Model1$BUGSoutput$summary["B","mean"]
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
posterior   <- Model1$BUGSoutput$sims.list
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
Model1_plot <- ggplot() +
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
    y     = "Risk of composite within 30 days",
    color = "Trial",
    title = "Model 1 Linear effect of timing “as randomised”"
  ) +
  theme_minimal()

Model1_plot

# ------------------------------
# 5) Log-scale y-axis (in %)
# ------------------------------
Model1_plot_log <- Model1_plot +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  labs(y = "Risk of composite within 30 days")

Model1_plot_log

Model1_plot_log <-Model1_plot_log+
  labs(
    title ="Model 1 Linear effect of timing “as randomised”",
    x ="",
    y = "Risk of composite at 30 days"
  ) + scale_y_log10(breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  xlim(c(0,25))
Model1_plot_log

Model1_plot_log <- Model1_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

Model1_plot_log

p_top <- Model1_plot_log +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )
# (εναλλακτικά: coord_cartesian(xlim = c(0, 25), expand = FALSE))
p_hist <- ggplot(syn_data, aes(x = avg_rando_timing_DOAC)) +
  geom_histogram(
    aes(y = after_stat(density)),   # keep density scale like the bell
    binwidth = 1,                   # adjust if you want more/less bins
    boundary = 0, closed = "left",
    fill = "steelblue",
    color = NA,                     # no black outline
    alpha = 0.25
  ) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  labs(x = "Time of DOAC administration (days)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

graph_Model1<-(p_top / p_hist) + plot_layout(heights = c(3, 1))


#### Figure for model 2 #####
rm(ModelRCS2_plot_log)

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
u_summaries <- Model2_MVN$BUGSoutput$summary[
  grep("^u\\[", rownames(Model2_MVN$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

# (B) Study-specific spline coefficients: Beta[i, 1], Beta[i, 2]
beta_index     <- grep("^Beta\\[", rownames(Model2_MVN$BUGSoutput$summary))
beta_summaries <- Model2_MVN$BUGSoutput$summary[beta_index, "mean"]

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
    events     = sum(composite_30, na.rm = TRUE),
    non_events = sum(1 - composite_30, na.rm = TRUE),
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
B1_mean   <- Model2_MVN$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model2_MVN$BUGSoutput$summary["B[2]", "mean"]

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
posterior <- Model2_MVN$BUGSoutput$sims.list
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
ModelRCS2_plot <- ggplot() +
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
    y     = "Risk of composite within 30 days",
    color = "Trial",
    title = "Model 2 Non-linear effect of timing “as randomised”"
  ) +
  theme_minimal()

ModelRCS2_plot

######################################
# 7) Optional: log-scale Y (still in %)
########################################
ModelRCS2_plot_log <- ModelRCS2_plot +
  scale_y_log10(
    breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0) +
  labs(y = "")

ModelRCS2_plot_log <-ModelRCS2_plot_log+
  labs(
    title ="Model 2 Non-linear effect of timing “as randomised”",
    x ="",
    y = ""
  ) + scale_y_log10(breaks = c(0.01, 0.03, 0.05, 0.07, 0.1,0.1), labels = label_percent(accuracy = 0.1))+
  coord_cartesian(ylim = c(0.005, 0.10), expand = 0)+
  xlim(c(0,25))

ModelRCS2_plot_log <- ModelRCS2_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

ModelRCS2_plot_log

p_top <- ModelRCS2_plot_log +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )
p_hist <- ggplot(syn_data, aes(x = avg_rando_timing_DOAC)) +
  geom_histogram(
    aes(y = after_stat(density)),
    binwidth = 1,            # adjust if you want more/less granularity
    boundary = 0, closed = "left",
    fill = "steelblue",
    color = NA,
    alpha = 0.25
  ) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  labs(x = "Time of DOAC administration (days)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

graph_Model2 <- (p_top / p_hist) + plot_layout(heights = c(3, 1))

#### Figure for model 3#####

rm(Model3_plot)
rm(Model3_plot_log)


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
posterior <- Model3_Linear$BUGSoutput$sims.list
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
Model3_plot_log <- ggplot() +
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
    y     = "Risk of composite at 30 days"
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

Model3_plot_log

ModelRCS2_plot_log <- ModelRCS2_plot_log +
  theme_minimal(base_size = 16) +   # bumps all text sizes globally
  theme(
    legend.position = "none",       # removes the right-hand legend
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 16)
  )

ModelRCS2_plot_log

p_top <- Model3_plot_log +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )
# (εναλλακτικά: coord_cartesian(xlim = c(0, 25), expand = FALSE))
p_bell <- ggplot(syn_data, aes(x = timing_DOAC)) +
  geom_area(
    stat   = "density",
    adjust = 1.0,
    alpha  = 0.25,
    fill   = "steelblue",
    color  = NA     # <- no black outline
  ) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  labs(x = "Time of DOAC administration (days)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

graph_Model3<-(p_top / p_bell) + plot_layout(heights = c(3, 1))


#### Figure for model 4 ####

rm(Model4_RCS_plot)
rm(Model4_RCS_plot_log)


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
u_summaries <- Model4_RCS$BUGSoutput$summary[
  grep("^u\\[", rownames(Model4_RCS$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

beta_index     <- grep("^Beta\\[", rownames(Model4_RCS$BUGSoutput$summary))
beta_summaries <- Model4_RCS$BUGSoutput$summary[beta_index, "mean"]

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
    events     = sum(composite_30, na.rm = TRUE),
    non_events = sum(1 - composite_30, na.rm = TRUE),
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
B1_mean   <- Model4_RCS$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model4_RCS$BUGSoutput$summary["B[2]", "mean"]

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
posterior <- Model4_RCS$BUGSoutput$sims.list
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
Model4_RCS_plot <- ggplot() +
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
    y     = "Risk of composite at 30 days"
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

Model4_RCS_plot

# ---- Log-scale version with your exact breaks/limits/labels ----
Model4_RCS_plot_log <- ggplot() +
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

Model4_RCS_plot_log

p_top <- Model4_RCS_plot_log +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin  = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  )
# (εναλλακτικά: coord_cartesian(xlim = c(0, 25), expand = FALSE))
p_bell <- ggplot(syn_data, aes(x = timing_DOAC)) +
  geom_area(
    stat   = "density",
    adjust = 1.0,
    alpha  = 0.25,
    fill   = "steelblue",
    color  = NA     # <- no black outline
  ) +
  scale_x_continuous(limits = c(0, 25), expand = c(0, 0)) +
  labs(x = "Time of DOAC administration (days)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

graph_Model4<-(p_top / p_bell) + plot_layout(heights = c(3, 1))

legend_order <- c(
  "ELAN minor", "ELAN moderate", "ELAN major",
  "TIMING", "OPTIMAS", "START"
)

# 1) 2×2 grid using the *bell* versions
grid_2x2 <- (graph_Model1 | graph_Model2) /
  (graph_Model3 | graph_Model4)

# 2) Trial legend (reuse a base plot that has the color mapping)
trial_legend <- cowplot::get_legend(
  Model1_plot +
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
      color = guide_legend(
        override.aes = list(linewidth = 2),
        keywidth  = unit(1.6, "cm"),
        keyheight = unit(0.6, "cm")
      )
    ) +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 14)
    )
)

# 3) Pooled legend (black thick line + star), same styling as before
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

# 4) Stack legends and place to the right of the grid
legend_stack <- cowplot::plot_grid(trial_legend, pooled_legend,
                                   ncol = 1, rel_heights = c(1, 0.35))

final_2x2_with_legend <- grid_2x2 | legend_stack
final_2x2_with_legend <- final_2x2_with_legend + plot_layout(widths = c(1, 0.20))

final_2x2_with_legend
# Save to PDF
ggsave(
  "Figures/Figure2.pdf",
  final_2x2_with_legend,
  device = cairo_pdf,  # or device = "pdf"
  width = 18, height = 10, units = "in"
)