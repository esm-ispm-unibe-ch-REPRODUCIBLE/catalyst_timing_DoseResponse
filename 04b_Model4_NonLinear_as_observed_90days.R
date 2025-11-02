
#####################################################################################################
############  Model 4: Non Linear as observed Meta-Analysis #################
###################### Secondary outcome: composite 90   ###############################33
######################################################################################

# Ensure composite_90 is numeric
syn_data <- syn_data %>%
  mutate(composite_90 = as.numeric(as.character(composite_90)))

# Calculate reference timing for each study (T_i0)
ref_timing <- syn_data %>%
  group_by(trial) %>%
  dplyr::summarize(ref_timing =min(timing_DOAC[composite_90 == 0], na.rm = TRUE))

# Merge ref_timing back into the dataset
syn_data <- syn_data %>%
  left_join(ref_timing, by = "trial")

# Calculate adjusted timing
syn_data <- syn_data %>%
  mutate(adjusted_timing = timing_DOAC - ref_timing)
syn_data <- syn_data[!is.na(syn_data$timing_DOAC), ]


syn_data_ELANmode <- syn_data[syn_data$trial_name=="ELAN moderate",]
syn_data_ELANmin <- syn_data[syn_data$trial_name=="ELAN minor",]
syn_data_ELANmaj <- syn_data[syn_data$trial_name=="ELAN major",]
syn_data_TIMING<- syn_data[syn_data$trial_name=="TIMING",]
syn_data_OPTIMAS <- syn_data[syn_data$trial_name=="OPTIMAS",]
syn_data_START <- syn_data[syn_data$trial_name=="START",]


my_knots_ELANmod <- with(syn_data, 
                         quantile(adjusted_timing, probs = c(0.1, 0.5, 0.9))
)

fit_glm_rcs_ELANmod <- glm(
  composite_90 ~ rcs(adjusted_timing, my_knots_ELANmod, 3),
  data   = syn_data_ELANmode,
  family = binomial()
)

my_knots_ELANmin <- with(syn_data, 
                         quantile(adjusted_timing, probs = c(0.1, 0.5, 0.9))
)

fit_glm_rcs_ELANmin <- glm(
  composite_90 ~ rcs(adjusted_timing, my_knots_ELANmin, 3),
  data   = syn_data_ELANmin,
  family = binomial()
)

my_knots_ELANmaj <- with(syn_data, 
                         quantile(adjusted_timing, probs = c(0.1, 0.5, 0.9))
)

fit_glm_rcs_ELANmaj <- glm(
  composite_90 ~ rcs(adjusted_timing, my_knots_ELANmaj, 3),
  data   = syn_data_ELANmaj,
  family = binomial()
)
summary(fit_glm_rcs_ELANmaj)
fit_glm_ELANmaj <- glm(
  composite_90 ~ adjusted_timing,
  data   = syn_data_ELANmaj,
  family = binomial()
)
summary(fit_glm_ELANmaj)

syn_data_TIMING<-syn_data_TIMING[!is.na(syn_data_TIMING$adjusted_timing),]

my_knots_TIMING <- with(syn_data, 
                        quantile(adjusted_timing, probs = c(0.1, 0.5, 0.9))
)

fit_glm_rcs_TIMING<- glm(
  composite_90 ~ rcs(adjusted_timing, my_knots_TIMING, 3),
  data   = syn_data_TIMING,
  family = binomial()
)

syn_data_OPTIMAS<-syn_data_OPTIMAS[!is.na(syn_data_OPTIMAS$adjusted_timing),]

my_knots_OPTIMAS <- with(syn_data, 
                         quantile(adjusted_timing, probs = c(0.1, 0.5, 0.9))
)

fit_glm_rcs_OPTIMAS <- glm(
  composite_90 ~ rcs(adjusted_timing, my_knots_OPTIMAS, 3),
  data   = syn_data_OPTIMAS,
  family = binomial()
)

syn_data_START<-syn_data_START[!is.na(syn_data_START$adjusted_timing),]

my_knots_START <- with(syn_data, 
                       quantile(adjusted_timing, probs = c(0.1, 0.5, 0.9))
)

syn_data_START<-syn_data_START[!is.na(syn_data_START$adjusted_timing),]

fit_glm_rcs_START <- glm(
  composite_90 ~ rcs(adjusted_timing, my_knots_START, 3),
  data   = syn_data_START,
  family = binomial()
)

study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START", "ELAN minor", "ELAN major")
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events = sum(composite_90, na.rm = TRUE),
    non_events = sum(1 - composite_90, na.rm = TRUE)
  ) %>%
  mutate(
    weight = 1/(1 / events + 1 / non_events)
  )

predict_glm_rcs <- function(data, fit, study_name, n_points = 100) {
  # 1) Range of timing_DOAC in this study
  min_t <- min(data$timing_DOAC, na.rm = TRUE)
  max_t <- max(data$timing_DOAC, na.rm = TRUE)
  
  # 2) Build an evenly spaced sequence from min to max
  x_seq <- seq(min_t, max_t, length.out = n_points)
  
  # 3) Convert each point to adjusted timing
  #    (assuming each study has a single unique ref_timing)
  ref_t <- unique(data$ref_timing)
  x_seq_adj <- x_seq - ref_t
  
  # 4) Predicted probability
  pred_prob <- predict(
    fit, 
    newdata = data.frame(adjusted_timing = x_seq_adj),
    type = "response"
  )
  
  # 5) Store both the original timing_DOAC (x_seq) and the predicted probability
  data.frame(
    timing_DOAC          = x_seq,      # original scale for plotting
    predicted_probability = pred_prob, 
    Study                = study_name
  )
}


df_ELANmod <- predict_glm_rcs(syn_data_ELANmode, fit_glm_rcs_ELANmod, "ELAN moderate")
df_ELANmin <- predict_glm_rcs(syn_data_ELANmin, fit_glm_rcs_ELANmin, "ELAN minor")
df_ELANmaj <- predict_glm_rcs(syn_data_ELANmaj, fit_glm_rcs_ELANmaj, "ELAN major")
df_TIMING <- predict_glm_rcs(syn_data_TIMING, fit_glm_rcs_TIMING, "TIMING")
df_OPTIMAS <- predict_glm_rcs(syn_data_OPTIMAS, fit_glm_rcs_OPTIMAS, "OPTIMAS")
df_START <- predict_glm_rcs(syn_data_START, fit_glm_rcs_START, "START")

# Combine all in one data frame
df_all_glm <- bind_rows(df_ELANmod, df_ELANmin,df_ELANmaj,df_ELANmod, df_TIMING, df_OPTIMAS, df_START)
# 1) Merge the sample size into df_all_glm
df_all_glm2 <- df_all_glm %>%
  left_join(study_sizes, by = c("Study" = "trial_name"))

observed_data_binned <- syn_data %>%
  group_by(trial_name) %>%
  do({
    # 1. Identify the min & max DOAC timing for this trial
    min_t <- min(.$timing_DOAC, na.rm = TRUE)
    max_t <- max(.$timing_DOAC, na.rm = TRUE)
    
    # 2. Create the sequence of breaks (2-day steps)
    breaks_vec <- seq(min_t, max_t, by = 2)

    
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
        obs_risk   = mean(composite_90, na.rm = TRUE),
        .groups    = "drop"
      )
  }) %>%
  ungroup()



# 2) Plot with line widths mapped to sample_size
Model5_90_IPD_seperate <- ggplot(df_all_glm2, 
                              aes(x = timing_DOAC, 
                                  y = predicted_probability, 
                                  color = Study,
                                  size = weight)) +  # <--- Map size here
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Splines for Each Study (RCS) with Line Width ~ Study Size",
    x = "Adjusted Timing",
    y = "Predicted Probability"
  ) +
  # (2) Observed means as points
  geom_point(
    data = observed_data_binned,
    aes(
      x     = timing_bin,
      y     = obs_risk,
      color = trial_name
    ),
    shape = 17,      # triangle (pick what you like: 16=circle, 17=triangle, etc.)
    size  = 3,
    alpha = 0.9
  ) +
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
  scale_size_continuous(range = c(0.5, 3), guide = "none")+
  scale_y_log10(limits = c(0.001, 0.1),
                breaks = c(0.001, 0.01, 0.03, 0.05, 0.07, 0.1))
#scale_y_log10(limits = c(0.001, 0.1),
#     breaks = c(0.001, 0.01, 0.03, 0.05, 0.07, 0.1))

Model5_90_IPD_seperate
# Explicitly generate the RCS basis with 3 knots
knots <- quantile(syn_data$timing_DOAC, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)
rcs_basis <- rcspline.eval(syn_data$timing_DOAC, knots = knots, inclx = TRUE)

# Add the RCS basis to the dataset
syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )


# Prepare the data for JAGS
jags_data <- list(
  Np = nrow(syn_data),  # Number of patients
  Nstudies = length(unique(syn_data$trial)),  # Number of studies
  trial = syn_data$trial,  # Study ID
  composite_90 = syn_data$composite_90,  # Outcome
  spline1 = syn_data$spline1,  # First spline basis
  spline2 = syn_data$spline2   # Second spline basis
)

# -- JAGS model with MVN random effects -- #

modelSplineTimingMVN <- function() {
  
  # 1) Likelihood
  for (k in 1:Np) {
    composite_90[k] ~ dbern(p[k])
    logit(p[k]) <- u[trial[k]] 
    + Beta[trial[k],1] * spline1[k]
    + Beta[trial[k],2] * spline2[k]
    p_pred[k] <- p[k]  # optional
  }
  
  # 2) Random intercept & correlated random slopes
  for (i in 1:Nstudies) {
    # Intercept
    u[i] ~ dnorm(0, 0.001)
    # Two correlated random slopes
    Beta[i, 1:2] ~ dmnorm(B[], Omega[,])
  }
  
  # 3) Priors for overall (pooled) spline coefficients
  for (j in 1:2) {
    B[j] ~ dnorm(0, 0.01)
  }
  
  # 4) Build correlation structure for Beta
  s1 ~ dnorm(0, 1)%_%T(0,)
  s2 ~ dnorm(0, 1)%_%T(0,)
  rho ~ dunif(-1, 1)
  
  # Covariance
  Sigma[1,1] <- pow(s1, 2)
  Sigma[2,2] <- pow(s2, 2)
  Sigma[1,2] <- rho * s1 * s2
  Sigma[2,1] <- rho * s1 * s2
  
  # Precision
  Omega <- inverse(Sigma)
}


# Run the JAGS model
set.seed(1294821)

Model4_90_RCS <- jags.parallel(
  data = jags_data,
  inits = NULL,
  parameters.to.save = c("B", "Sigma", "Omega", "u", "Beta", "rho"),
  model.file = modelSplineTimingMVN,
  n.chains = 2,
  n.iter = 50000,
  n.burnin = 2000,
  DIC = FALSE,
  n.thin = 10
)

print(Model4_90_RCS)

# Extract the MCMC summary from the MVN model
mcmc_sum <- Model4_90_RCS$BUGSoutput$summary

# 1. Extract the summary rows for B1, B2, Sigma[1,1], Sigma[2,2], and Sigma[1,2]
sum_B1     <- mcmc_sum["B[1]", ]        # row for B1
sum_B2     <- mcmc_sum["B[2]", ]        # row for B2
sum_Sigma11 <- mcmc_sum["Sigma[1,1]", ] # row for Sigma[1,1]
sum_Sigma22 <- mcmc_sum["Sigma[2,2]", ] # row for Sigma[2,2]
sum_Sigma12 <- mcmc_sum["Sigma[1,2]", ] # row for Sigma[1,2]
sum_rho <- mcmc_sum["rho", ] # row for Sigma[1,2]

# 2. Pull out mean and credible intervals for B1, B2, Sigma[1,1], Sigma[2,2], and Sigma[1,2]
B1_mean   <- sum_B1["mean"]
B1_lower  <- sum_B1["2.5%"]
B1_upper  <- sum_B1["97.5%"]

B2_mean   <- sum_B2["mean"]
B2_lower  <- sum_B2["2.5%"]
B2_upper  <- sum_B2["97.5%"]

rho_mean   <- sum_rho["mean"]
rho_lower  <- sum_rho["2.5%"]
rho_upper  <- sum_rho["97.5%"]

Sigma11_mean  <- sum_Sigma11["mean"]
Sigma11_lower <- sum_Sigma11["2.5%"]
Sigma11_upper <- sum_Sigma11["97.5%"]

Sigma22_mean  <- sum_Sigma22["mean"]
Sigma22_lower <- sum_Sigma22["2.5%"]
Sigma22_upper <- sum_Sigma22["97.5%"]

Sigma12_mean  <- sum_Sigma12["mean"]
Sigma12_lower <- sum_Sigma12["2.5%"]
Sigma12_upper <- sum_Sigma12["97.5%"]

# 4. Convert B1 and B2 from log-odds to odds ratios
OR1_mean  <- exp(B1_mean)
OR1_lower <- exp(B1_lower)
OR1_upper <- exp(B1_upper)

OR2_mean  <- exp(B2_mean)
OR2_lower <- exp(B2_lower)
OR2_upper <- exp(B2_upper)

# 5. Calculate the correlation between beta1 and beta2
# Correlation = Sigma12 / sqrt(Sigma11 * Sigma22)


# 6. Build the results table
results_table <- tibble(
  parameter = c("B1 (spline1)", "B2 (spline2)", 
                "Sigma[1,1] (Variance beta1)", 
                "Sigma[2,2] (Variance beta2)", 
                "Sigma[1,2] (Covariance)", 
                "rho"),
  
  # Column 1: Estimates with 95% Credible Intervals
  estimate = c(
    sprintf("%.3f (%.3f to %.3f)", B1_mean, B1_lower, B1_upper),
    sprintf("%.3f (%.3f to %.3f)", B2_mean, B2_lower, B2_upper),
    sprintf("%.3f (%.3f to %.3f)", Sigma11_mean, Sigma11_lower, Sigma11_upper),
    sprintf("%.3f (%.3f to %.3f)", Sigma22_mean, Sigma22_lower, Sigma22_upper),
    sprintf("%.3f (%.3f to %.3f)", Sigma12_mean, Sigma12_lower, Sigma12_upper),
    sprintf("%.3f (%.3f to %.3f)", rho_mean, rho_lower, rho_upper)
  ),
  
  # Column 2: Odds Ratios for B1 and B2; '—' for others
  odds_ratio = c(
    sprintf("%.3f (%.3f to %.3f)", OR1_mean, OR1_lower, OR1_upper),
    sprintf("%.3f (%.3f to %.3f)", OR2_mean, OR2_lower, OR2_upper),
    "—",  # Sigma[1,1] doesn't have an OR meaning
    "—",  # Sigma[2,2] doesn't have an OR meaning
    "—",  # Sigma[1,2] doesn't have an OR meaning
    "—"   # Correlation doesn't have an OR meaning
  )
)

# 7. Convert to a data frame for better compatibility
results_table_df <- as.data.frame(results_table)
# Print the results table to the console
print(results_table_df)


########################################
# 1) Define knots & build 2-column RCS
########################################
# We'll pick knots at the 10th, 50th, 90th percentiles of `timing_DOAC`
knots <- quantile(
  x     = syn_data$timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

# Just to illustrate, let's see them:
# print(knots)
rcs_basis <- rcspline.eval(
  x     = syn_data$timing_DOAC,
  knots = knots,
  inclx = TRUE
)

# 'rcs_basis' now has exactly 2 columns, typically something like:
#   col1 = the first RCS expansion
#   col2 = the second RCS expansion
syn_data <- syn_data %>%
  mutate(
    spline1 = rcs_basis[, 1],
    spline2 = rcs_basis[, 2]
  )

########################################
# 2) Extract posterior means for each study
########################################

# (A) Random intercepts: u[i]
u_summaries <- Model4_90_RCS$BUGSoutput$summary[
  grep("^u\\[", rownames(Model4_90_RCS$BUGSoutput$summary)),
  "mean"
]
Nstudies <- length(u_summaries)

# (B) Study-specific spline coefficients: Beta[i, 1], Beta[i, 2]
beta_index     <- grep("^Beta\\[", rownames(Model4_90_RCS$BUGSoutput$summary))
beta_summaries <- Model4_90_RCS$BUGSoutput$summary[beta_index, "mean"]

# Reshape into an Nstudies x 2 matrix
Beta_matrix <- matrix(
  data  = beta_summaries,
  nrow  = Nstudies,
  ncol  = 2,
  byrow = FALSE
)

# (C) Identify your study IDs in the model's order
study_ids <- c("ELAN moderate", "TIMING", "OPTIMAS", "START",
               "ELAN minor", "ELAN major")

########################################
# 3) Lines from the model for each study
########################################

# If you'd like line width scaled by sample size, compute that:
study_sizes <- syn_data %>%
  group_by(trial_name) %>%
  summarise(
    events     = sum(composite_90, na.rm = TRUE),
    non_events = sum(1 - composite_90, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(
    weight = 1 / (1 / events + 1 / non_events)
  )

plotdata_studies <- do.call(rbind, lapply(seq_len(Nstudies), function(i) {
  study_i     <- study_ids[i]
  intercept_i <- u_summaries[i]
  slope1_i    <- Beta_matrix[i, 1]
  slope2_i    <- Beta_matrix[i, 2]
  
  # Filter data for that study
  df_sub <- syn_data %>% filter(trial_name == study_i)
  x_min  <- min(df_sub$timing_DOAC, na.rm = TRUE)
  x_max  <- max(df_sub$timing_DOAC, na.rm = TRUE)
  
  x_seq_i <- seq(x_min, x_max, length.out = 200)
  
  # Build the RCS basis for this x range, with inclx=FALSE so we get 2 columns
  rcs_i <- rcspline.eval(x_seq_i, knots = knots, inclx = TRUE)
  # col1 = first expansion, col2 = second expansion
  
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
# 4) Single "pooled" line from means
########################################
# We'll average all u[i] for alpha, and take mean B[1], B[2]

u_overall <- mean(u_summaries)
B1_mean   <- Model4_90_RCS$BUGSoutput$summary["B[1]", "mean"]
B2_mean   <- Model4_90_RCS$BUGSoutput$summary["B[2]", "mean"]

# Build a global x range
global_min <- min(syn_data$timing_DOAC, na.rm = TRUE)
global_max <- max(syn_data$timing_DOAC, na.rm = TRUE)
x_seq_overall <- seq(global_min, global_max, length.out = 200)

# RCS expansions
rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

logit_overall <- u_overall +
  B1_mean * rcs_overall[,1] +
  B2_mean * rcs_overall[,2]
risk_overall <- plogis(logit_overall)

df_overall <- data.frame(
  x    = x_seq_overall,
  risk = risk_overall
)

########################################
# 5) 95% CI for that pooled line (iteration-wise)
########################################
posterior <- Model4_90_RCS$BUGSoutput$sims.list

# Suppose B is a matrix with 2 columns in 'posterior$B'
B_draws_1 <- posterior$B[, 1]
B_draws_2 <- posterior$B[, 2]

# 'u_draws' is nIter x Nstudies
u_draws <- posterior$u
alpha_draws <- rowMeans(u_draws)  # iteration-wise average intercept

nIter   <- length(B_draws_1)
nPoints <- length(x_seq_overall)
pred_mat <- matrix(NA, nrow = nIter, ncol = nPoints)

for (m in seq_len(nIter)) {
  alpha_m <- alpha_draws[m]
  b1_m    <- B_draws_1[m]
  b2_m    <- B_draws_2[m]
  
  # logit at each x
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
# 6) Observed means as triangles
########################################

observed_data_binned <- syn_data %>%
  group_by(trial_name) %>%
  do({
    # 1. Identify the min & max DOAC timing for this trial
    min_t <- min(.$timing_DOAC, na.rm = TRUE)
    max_t <- max(.$timing_DOAC, na.rm = TRUE)
    
    # 2. Create the sequence of breaks (2-day steps)
    breaks_vec <- seq(min_t, max_t, by = 2)
  
    
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
        obs_risk   = mean(composite_90, na.rm = TRUE),
        .groups    = "drop"
      )
  }) %>%
  ungroup()

########################################
# 7) Final Plot
########################################
Model4_90_RCS_plot <- ggplot() +
  
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
  
  # (C) Dashed single best-fit pooled line
  geom_line(
    data = df_overall,
    aes(x = x, y = risk),
    color    = "blue",
    linetype = "dashed",
    size     = 1.2
  ) +
  
  # (D) Observed mean event rates as triangles
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
  
  # (E) Manual color scale
  scale_color_manual(
    values = c(
      "ELAN moderate" = "green",
      "ELAN minor"    = "lightgreen",
      "ELAN major"    = "darkgreen",
      "START"         = "#e7298a",
      "TIMING"        = "darkred",
      "OPTIMAS"       = "#e6ab02"
      # "Pooled"       = "blue" # assigned explicitly in geom_line
    )
  ) +
  
  # (F) Triangle shape for "Observed"
  scale_shape_manual(values = c("Observed" = 17)) +
  
  # (G) Line-size range
  scale_size_continuous(range = c(0.5, 3)) +
  
  # (H) Labels & theme
  labs(
    x     = "Timing in Days (ref = 0)",
    y     = "Observed & Predicted Risk",
    color = "Trial",
    shape = NULL,
    size  = NULL,
    title = "Model 2: Non-linear Timing-Risk with 95% CI"
  ) +
  theme_minimal() +scale_size_continuous(range = c(0.5, 3), guide = "none")+xlim(c(0,25))


# Print the plot
Model4_90_RCS_plot
Model4_90_RCS_plot_log <- Model4_90_RCS_plot +  scale_y_log10(breaks = c(0.001, 0.01, 0.03, 0.05, 0.07, 0.1))

Model4_90_RCS_plot_log


save(Model5_90_IPD_seperate, Model4_90_RCS, Model4_90_RCS_plot_log, Model4_90_RCS_plot, file="Saved_Data/Model4_90_RCS.RData")

