#####################################################################################################
############  Model 3: Linear as observed Meta-Analysis #################
###################### Primary outcome: composite 30   ###############################33
######################################################################################

# Ensure composite_30 is numeric
syn_data <- syn_data %>%
  mutate(composite_30 = as.numeric(as.character(composite_30)))

# Calculate reference timing for each study (T_i0)
ref_timing <- syn_data %>%
  group_by(trial) %>%
  dplyr::summarize(ref_timing = min(timing_DOAC[composite_30 == 0], na.rm = TRUE))

# Merge ref_timing back into the dataset
syn_data <- syn_data %>%
  left_join(ref_timing, by = "trial")

# Calculate adjusted timing
syn_data <- syn_data %>%
  mutate(adjusted_timing = timing_DOAC - ref_timing)

# Remove rows with missing timing_DOAC
syn_data <- syn_data[!is.na(syn_data$timing_DOAC), ]


### Independent lines 

# (A) line widths to reflect sample size:
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events = sum(composite_30, na.rm = TRUE),
    non_events = sum(1 - composite_30, na.rm = TRUE)
  ) %>%
  mutate(
    weight = 1/(1 / events + 1 / non_events)
  )
# (B) Fit a separate glm for each study and create predicted lines
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
  
  # Fit the GLM (logistic) with 'composite_30' as outcome
  # We'll do: composite_30 ~ I(avg_rando_timing_DOAC - ref_i)
  # so the intercept is the risk at reference, slope is the effect of timing
  mod <- glm(
    formula = composite_30 ~ I(adjusted_timing),
    data    = df_sub,
    family  = binomial
  )
  
  # Build a sequence from min to max for *this study*
  x_min <- min(df_sub$timing_DOAC, na.rm = TRUE)
  x_max <- max(df_sub$timing_DOAC, na.rm = TRUE)
  
  x_seq_study <- seq(x_min, x_max, length.out = 200)
  
  # Predict for each x in that sequence
  # type = "response" => probabilities
  pred_vals <- predict(
    mod,
    newdata = data.frame(adjusted_timing = x_seq_study),
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

# (C) Observed means: for each (study, timing), get mean outcome
# Create binned observed data for each trial using 2-day intervals
observed_data_binned <- syn_data %>%
  group_by(trial_name) %>%
  do({
    # 1. Identify the min & max DOAC timing for this trial
    min_t <- min(.$timing_DOAC, na.rm = TRUE)
    max_t <- max(.$timing_DOAC, na.rm = TRUE)
    
    # 2. Create the sequence of breaks (2-day steps)
    breaks_vec <- seq(min_t, max_t, by = 2)
    
    # Optionally ensure the final bin covers up to max_t, in case max_t 
    # isn't a multiple of 2:
    # if (breaks_vec[length(breaks_vec)] < max_t) {
    #   breaks_vec <- c(breaks_vec, max_t)
    # }
    
    # 3. Bin the data using cut()
    .x <- mutate(.,
                 bin = cut(
                   timing_DOAC,
                   breaks = breaks_vec,
                   right  = FALSE,       # intervals are [start, end)
                   include.lowest = TRUE
                 )
    )
    
    # 4. Summarize within each bin: compute average timing & average outcome
    .x %>%
      group_by(bin) %>%
      dplyr::summarize(
        timing_bin = mean(timing_DOAC, na.rm = TRUE),
        obs_risk   = mean(composite_30, na.rm = TRUE),
        .groups    = "drop"
      )
  }) %>%
  ungroup()

#-------------------------------------------------------
# (D) Plot: One line per study, no pooled line
#-------------------------------------------------------

Model3_separate_plot <- ggplot() +
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
    data = observed_data_binned,
    aes(
      x     = timing_bin,
      y     = obs_risk,
      color = trial_name
    ),
    shape = 17,      
    size  = 3,
    alpha = 0.9
  ) +
  # (3) color scale
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
  theme_minimal() + scale_size_continuous(range = c(0.5, 3), guide = "none")

Model3_separate_plot


Model3_separate_plot_logy <- Model3_separate_plot +
  scale_y_log10() 

Model3_separate_plot_logy

### Meta-analysis model
# Prepare the data for JAGS
jags_data <- list(
  Np = nrow(syn_data),                  # Number of patients
  Nstudies = length(unique(syn_data$trial)),  # Number of studies
  trial = syn_data$trial,               # Study ID for each patient
  composite_30 = syn_data$composite_30, # Outcome variable
  adjusted_timing = syn_data$adjusted_timing  # Adjusted timing (linear predictor)
)

# Define the JAGS model for the linear relationship
modelLinearTiming_3 <- function() {
  for (k in 1:Np) {
    composite_30[k] ~ dbern(p[k])
    logit(p[k]) <- u[trial[k]] + beta[trial[k]] * (adjusted_timing[k] )
    # Here we store a "predicted risk" node, if you want it separate
    p_pred[k] <- p[k]
    
  }
  
  # Priors for study-specific intercepts and slopes
  for (i in 1:Nstudies) {
    u[i] ~ dnorm(0, 0.001)
    beta[i] ~ dnorm(B, prec)
  }
  
  # Priors for pooled slope and variance
  B ~ dnorm(0, 0.01)          # Prior for the pooled slope
  tau ~ dnorm(0, 1)%_%T(0,)   # Prior for standard deviation
  tau2<-pow(tau, 2)
  prec <- 1 / pow(tau, 2)     # Precision as the inverse of variance
}

### run the model
set.seed(1294821)

Model3_Linear <- jags.parallel(data = jags_data,inits=NULL,parameters.to.save = c("B", "tau","tau2", "u", "beta"),model.file = modelLinearTiming_3,
                        n.chains=2,n.iter = 20000,n.burnin = 2000,DIC=F,n.thin = 10)
print(Model3_Linear)
#Model3_Linear_withPred <- jags.parallel(data = jags_data,inits=NULL,parameters.to.save = c("B", "tau", "u", "beta", "p_pred"),model.file = modelLinearTiming_3,
#                                 n.chains=2,n.iter = 20000,n.burnin = 2000,DIC=F,n.thin = 10)
#### Tables of results
mcmc_sum <- Model3_Linear$BUGSoutput$summary

# 1. Extract the summary rows for B and tau
sum_B   <- mcmc_sum["B",   ]     # row for B
sum_tau <- mcmc_sum["tau", ]     # row for tau

# 2. Pull out mean and credible intervals
B_mean   <- sum_B["mean"]
B_lower  <- sum_B["2.5%"]
B_upper  <- sum_B["97.5%"]

tau_mean  <- sum_tau["mean"]
tau_lower <- sum_tau["2.5%"]
tau_upper <- sum_tau["97.5%"]

# 3. Convert tau (SD) to tau^2 (variance)
tau2_mean  <- tau_mean^2
tau2_lower <- tau_lower^2
tau2_upper <- tau_upper^2

# 4. Convert B from log-odds to odds ratio
OR_mean  <- exp(B_mean)
OR_lower <- exp(B_lower)
OR_upper <- exp(B_upper)

# 5. Build the table
results_table <- tibble(
  parameter = c("B", "tau"),
  
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

# 6. Print the results
results_table1<-as.data.frame(results_table)
output1<-as.data.frame(round(mcmc_sum, 3))

################### PLOTS ################################

# Store the summary in a convenient object
model_summary <- Model3_Linear$BUGSoutput$summary

u_summaries <- Model3_Linear$BUGSoutput$summary[
  grep("^u\\[", rownames(Model3_Linear$BUGSoutput$summary)),
  "mean"
]
beta_summaries <- Model3_Linear$BUGSoutput$summary[
  grep("^beta\\[", rownames(Model3_Linear$BUGSoutput$summary)),
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
  
  # Create a per-study x-range
  x_min <- min(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
  x_max <- max(syn_data$timing_DOAC[syn_data$trial_name == study_i], na.rm = TRUE)
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
B_mean <- Model3_Linear$BUGSoutput$summary["B","mean"]

# Average (over studies) of the posterior means of u[i]
alpha_mean <- mean(u_summaries)

x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
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

posterior <- Model3_Linear$BUGSoutput$sims.list
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

observed_data_binned <- syn_data %>%
  group_by(trial_name) %>%
  do({
    min_t <- min(.$timing_DOAC, na.rm = TRUE)
    max_t <- max(.$timing_DOAC, na.rm = TRUE)
    
    breaks_vec <- seq(min_t, max_t, by = 2)
    
    .x <- mutate(.,
                 bin = cut(
                   timing_DOAC,
                   breaks = breaks_vec,
                   right  = FALSE,       # intervals are [start, end)
                   include.lowest = TRUE
                 )
    )
    
    .x %>%
      group_by(bin) %>%
      dplyr::summarize(
        timing_bin = mean(timing_DOAC, na.rm = TRUE),
        obs_risk   = mean(composite_30, na.rm = TRUE),
        .groups    = "drop"
      )
  }) %>%
  ungroup()


# ------------------------------------------------------------------
# Combine everything into one ggplot
# ------------------------------------------------------------------

Model3_plot<- ggplot() +
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
    data = observed_data_binned,
    aes(
      x = timing_bin,
      y = obs_risk,
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
  theme_minimal()+xlim(c(0,25))

Model3_plot

Model3_plot_log <- Model3_plot + scale_y_log10(breaks = c(0.001, 0.01, 0.03, 0.05, 0.07, 0.1))

Model3_plot_log


save(Model3_Linear, Model3_plot_log, Model3_plot, Model3_separate_plot, Model3_separate_plot_logy, file="Saved_Data/Model3_Linear.RData")

