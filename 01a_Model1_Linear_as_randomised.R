#####################################################################################################
############  Model 1: Linear as ransomized meta-analysis #################
###################### Primary outcome: composite 30   ###############################33
######################################################################################

# Calculate reference timing for each study (T_i0)
ref_timing <- syn_data %>%
  group_by(trial) %>%
  dplyr::summarize(ref_timing = min(avg_rando_timing_DOAC[composite_30 == 0], na.rm = TRUE))

# Merge ref_timing back into the dataset
syn_data <- syn_data %>%
  left_join(ref_timing, by = "trial")


##### GLM independent before the random-effects model

# line widths to reflect sample size:
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events = sum(composite_30, na.rm = TRUE),
    non_events = sum(1 - composite_30, na.rm = TRUE)
  ) %>%
  mutate(
    weight = 1/(1 / events + 1 / non_events)
  )

# Fit a separate glm for each study and create predicted lines
studies <- unique(syn_data$trial_name)

plotdata_studies <- lapply(studies, function(study_i) {
  # Subset data for this study
  df_sub <- syn_data %>%
    filter(trial_name == study_i)
  
  # Calculate a single reference timing for this study (if you have a single ref)
  ref_i <- unique(df_sub$ref_timing)
  if (length(ref_i) > 1) {
    ref_i <- mean(ref_i, na.rm = TRUE) # fallback if multiple arms
  }
  
   mod <- glm(
    formula = composite_30 ~ I(avg_rando_timing_DOAC - ref_i),
    data    = df_sub,
    family  = binomial
  )
  
  # Build a sequence from min to max for *this study*
  x_min <- min(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  x_max <- max(df_sub$avg_rando_timing_DOAC, na.rm = TRUE)
  
  x_seq_study <- seq(x_min, x_max, length.out = 200)
  
  # Predict for each x in that sequence
  # type = "response" => probabilities
  pred_vals <- predict(
    mod,
    newdata = data.frame(avg_rando_timing_DOAC = x_seq_study),
    type    = "response"
  )
  
  # Attach the study size if you want to map line width
  this_size <- study_sizes$weight[study_sizes$trial_name == study_i]
  
  data.frame(
    study      = study_i,
    x          = x_seq_study,
    pred_risk  = pred_vals,
    study_size = this_size
  )
})

# Combine all into one data.frame
plotdata_studies <- do.call(rbind, plotdata_studies)

# Observed means: for each (study, timing), get mean outcome
observed_means <- syn_data %>%
  group_by(trial_name, avg_rando_timing_DOAC) %>%
  dplyr::summarize(
    obs_mean  = mean(composite_30, na.rm = TRUE),
    obs_count = n(),
    .groups   = "drop"
  )

#-------------------------------------------------------
# Plot: One line per study, no pooled line
#-------------------------------------------------------

Model1_separate_plot <- ggplot() +
  # (1) Fitted lines, sized by study_size
  geom_line(
    data = plotdata_studies,
    aes(
      x     = x,
      y     = pred_risk,
      color = study,
      size  = study_size
    )
  ) +
  # (2) Observed means as points
  geom_point(
    data = observed_means,
    aes(
      x     = avg_rando_timing_DOAC,
      y     = obs_mean,
      color = trial_name,
      shape = "Observed"
    ),
    size = 3
  ) +
  # (3) Pick your color scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02",
      "Pooled"        = "blue"
    )
  ) +
  # (4) Shape for observed data
  scale_shape_manual(values = c("Observed" = 17)) +
  # (5) Map line size range
  scale_size_continuous(range = c(0.5, 3)) +
  # (6) Labels
  labs(
    x     = "Timing in Days",
    y     = "Risk (Observed vs Predicted)",
    color = "Trial",
    shape = NULL,
    size  = NULL,
    title = "Separate GLMs by Study: Linear Timing-Risk"
  ) +
  theme_minimal() +scale_size_continuous(range = c(0.5, 3), guide = "none")

Model1_separate_plot


Model1_separate_plot_logy <- Model1_separate_plot +
  scale_y_log10() 

Model1_separate_plot_logy

#### Random-effects model

# Prepare the data for JAGS
jags_data <- list(
  Np = nrow(syn_data),  # Number of patients
  Nstudies = length(unique(syn_data$trial)),  # Number of studies
  trial = syn_data$trial,  # Study ID
  composite_30 = syn_data$composite_30,  # Outcome
  avg_rando_timing_DOAC = syn_data$avg_rando_timing_DOAC,  # Timing
  ref_timing = syn_data$ref_timing  # Reference timing for each study
)

modelLinearTiming <- function() {
  for (k in 1:Np) {
    composite_30[k] ~ dbern(p[k])
    logit(p[k]) <- u[trial[k]] + beta[trial[k]] * (avg_rando_timing_DOAC[k] - ref_timing[trial[k]])
    # Here we store a "predicted risk" node, if you want it separate
    p_pred[k] <- p[k]
  
    }
  
  # Priors for study-specific intercepts and slopes
  for (i in 1:Nstudies) {
    u[i] ~ dnorm(0, 0.001)
    beta[i] ~ dnorm(B, prec)
  }
  
  # Priors for pooled slope and variance
  B ~ dnorm(0, 0.001)          # Prior for the pooled slope
  tau ~ dnorm(0, 1)%_%T(0,)   # Prior for standard deviation
  tau2<-pow(tau, 2)
  prec <- 1 / pow(tau, 2)     # Precision as the inverse of variance
}

### run the model
set.seed(1294821)

Model1 <- jags.parallel(data = jags_data,inits=NULL,parameters.to.save = c("B", "tau","tau2", "u", "beta"),model.file = modelLinearTiming,
                                        n.chains=2,n.iter = 20000,n.burnin = 2000,DIC=F,n.thin = 10)
print(Model1, digits = 2)


#traceplot(as.mcmc(Model1))

#### Tables of results
mcmc_sum <- Model1$BUGSoutput$summary

# Extract the summary rows for B and tau
sum_B   <- mcmc_sum["B",   ]     # row for B
sum_tau <- mcmc_sum["tau2", ]     # row for tau

#  Pull out mean and credible intervals
B_mean   <- sum_B["mean"]
B_lower  <- sum_B["2.5%"]
B_upper  <- sum_B["97.5%"]

tau_mean  <- sum_tau["mean"]
tau_lower <- sum_tau["2.5%"]
tau_upper <- sum_tau["97.5%"]

# Convert tau (SD) to tau^2 (variance)
tau_mean  <- tau_mean
tau_lower <- tau_lower
tau_upper <- tau_upper

# 4. Convert B from log-odds to odds ratio
OR_mean  <- exp(B_mean)
OR_lower <- exp(B_lower)
OR_upper <- exp(B_upper)

# Build the table
results_table <- tibble(
  parameter = c("B", "tau2"),
  
  # Column 1: log-odds (95% CrI)
  log_odds  = c(
    sprintf("%.3f (%.3f to %.3f)", B_mean, B_lower, B_upper),
    sprintf("%.3f (%.3f to %.3f)", tau_mean, tau_lower, tau_upper)
  ),
  
  # Column 2: odds ratio (95% CrI)
  odds_ratio = c(
    sprintf("%.3f (%.3f to %.3f)", OR_mean, OR_lower, OR_upper),
    "—"  # tau^2 doesn't have an OR meaning, so we'll leave it blank or '—'
  )
)

# Print the results
results_table1<-as.data.frame(results_table)
output1<-as.data.frame(round(mcmc_sum, 3))

# Posterior means for study-specific intercepts and slopes
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
  
  ref_i <- 0
  
  # Create a per-study x-range
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

# ------------------------------------------------------------------
# Single "mean" pooled line (as in your old code) 
# ------------------------------------------------------------------

# Posterior mean of the slope "B"
B_mean <- Model1$BUGSoutput$summary["B","mean"]

# Average (over studies) of the posterior means of u[i]
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

# ------------------------------------------------------------------
#  95% CI from iteration-wise draws of "average intercept" + B
# ------------------------------------------------------------------

posterior <- Model1$BUGSoutput$sims.list
B_draws   <- posterior$B          # length = nIter
u_draws   <- posterior$u          # dimension: nIter x Nstudies

# alpha_draws[m] = mean of the random intercepts in iteration m
alpha_draws <- rowMeans(u_draws)
nIter <- length(B_draws)

pred_mat <- matrix(NA, nrow = nIter, ncol = length(x_seq_pooled))
for (m in seq_len(nIter)) {
  # alpha_m and B_m
  alpha_m <- alpha_draws[m]
  B_m     <- B_draws[m]
  
  # logit = alpha_m + B_m * (x - 0)
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

# ------------------------------------------------------------------
#  Observed means 
# ------------------------------------------------------------------

observed_means <- syn_data %>%
  group_by(trial_name, avg_rando_timing_DOAC) %>%
  dplyr::summarize(
    obs_mean  = mean(composite_30, na.rm = TRUE),
    obs_count = n(),
    .groups   = "drop"
  )

# ------------------------------------------------------------------
# Combine everything into one ggplot
# ------------------------------------------------------------------

Model1_plot <- ggplot() +
  # (A) Study-specific lines (no changes)
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
  
  # (C) Single pooled line from the posterior means
  geom_line(
    data = df_line,
    aes(x = x, y = risk),
    color    = "blue",
    linetype = "dashed",
    size     = 1.2
  ) +
  
  # (D) Observed means as triangles
  geom_point(
    data = observed_means,
    aes(
      x = avg_rando_timing_DOAC,
      y = obs_mean,
      color = trial_name,
      shape = "Observed"
    ),
    size = 3
  ) +
  
  # (E) Manual colors, shapes, etc.
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
  scale_shape_manual(values = c("Observed" = 17)) +
  scale_size_continuous(range = c(0.5, 3), guide = "none") +
  
  labs(
    x     = "Timing in Days (ref = 0)",
    y     = "Risk (Observed vs Predicted)",
    color = "Trial",
    shape = NULL,
    size  = NULL,
    title = "Model 1: Linear Timing-Risk Meta-Analysis (Ref = 0)"
  ) +
  theme_minimal() +xlim(c(0,25))

Model1_plot

Model1_plot_log <- Model1_plot + scale_y_log10(breaks = c(0.001, 0.01, 0.03, 0.05, 0.07, 0.1, 0.10), limits = c(0.001, 0.12))

Model1_plot_log
save(Model1, Model1_plot_log, Model1_plot, Model1_separate_plot, Model1_separate_plot_logy, file="Saved_Data/Model1_Linear.RData")
