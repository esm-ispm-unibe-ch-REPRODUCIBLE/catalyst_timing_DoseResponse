#####################################################################################################
##### Model 4: Multiple imputations Non Linear as observed Meta-Analysis #################
###################### Primary outcome: composite 30   ###############################33
######################################################################################

syn_data<-syn_data[!is.na(syn_data$composite_30),]
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

# Select relevant variables
MIdata <- syn_data[, c("timing_DOAC", "nihss_adm", "composite_30", "avg_rando_timing_DOAC", "trial_name", "trial", "ref_timing")]

# Ensure composite_30 is a factor
MIdata$composite_30 <- as.factor(MIdata$composite_30)

fml<-nihss_adm+timing_DOAC~composite_30+(1+syn_data$avg_rando_timing_DOAC)|trial

set.seed(2000)
imp3<-jomoImpute(data=MIdata,formula=fml,n.burn=1000, n.iter=1000, m=10, seed=1569) #about 1 minute
summary(imp3)
#check the convergence
#plot(imp3)
#imputed datasets
imp.list<-mitmlComplete(imp3,print = "all")

for (i in 1:10) {
  # Create adjusted_timing variable
  imp.list[[i]]$adjusted_timing <- imp.list[[i]]$timing_DOAC - imp.list[[i]]$ref_timing
  
}

model_results_rcs <- list()

# Define RCS knots (fixed across imputations)
knots <- quantile(syn_data$timing_DOAC, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

#### Run the Bayesian model for all the imputed datasets

modelSplineTimingMVN <- function() {
  
  # 1) Likelihood
  for (k in 1:Np) {
    composite_30[k] ~ dbern(p[k])
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

# Loop through each imputed dataset
for (i in 1:10) {
  cat(sprintf("Running RCS JAGS on imputed dataset %d...\n", i))
  
  imp_data <- imp.list[[i]]
  
  # Make sure composite_30 is numeric binary (0/1)
  imp_data$composite_30 <- as.numeric(as.character(imp_data$composite_30))
  
  # Create spline basis
  rcs_basis <- rcspline.eval(imp_data$timing_DOAC, knots = knots, inclx = TRUE)
  imp_data$spline1 <- rcs_basis[, 1]
  imp_data$spline2 <- rcs_basis[, 2]
  
  # Prepare JAGS data
  jags_data <- list(
    Np = nrow(imp_data),
    Nstudies = length(unique(imp_data$trial)),
    trial = as.numeric(as.factor(imp_data$trial)),
    composite_30 = imp_data$composite_30,
    spline1 = imp_data$spline1,
    spline2 = imp_data$spline2
  )
  
  # Fit JAGS model
  set.seed(1294821 + i)
  fit_rcs <- jags.parallel(
    data = jags_data,
    inits = NULL,
    parameters.to.save = c("B", "Sigma", "u", "Beta", "rho"),
    model.file = modelSplineTimingMVN,
    n.chains = 2,
    n.iter = 50000,
    n.burnin = 2000,
    DIC = FALSE,
    n.thin = 10
  )
  
  model_results_rcs[[i]] <- fit_rcs
}

### RUbin's Rules
# 1. Initialize vectors to store estimates and SEs
B1_means  <- B2_means  <- rho_means  <- numeric(10)
B1_ses    <- B2_ses    <- rho_ses    <- numeric(10)

for (i in 1:10) {
  sum_i <- model_results_rcs[[i]]$BUGSoutput$summary
  
  B1_means[i]   <- sum_i["B[1]", "mean"]
  B1_ses[i]     <- sum_i["B[1]", "sd"]
  
  B2_means[i]   <- sum_i["B[2]", "mean"]
  B2_ses[i]     <- sum_i["B[2]", "sd"]
  
  rho_means[i]  <- sum_i["rho", "mean"]
  rho_ses[i]    <- sum_i["rho", "sd"]
}

# 2. Create qhat (means) and uhat (within-imputation variances)
qhat <- rbind(
  `B1 (spline 1)` = B1_means,
  `B2 (spline 2)` = B2_means,
  `rho`           = rho_means
)

uhat <- rbind(
  `B1 (spline 1)` = B1_ses^2,
  `B2 (spline 2)` = B2_ses^2,
  `rho`           = rho_ses^2
)

# 3. Run Rubin’s rules
pooled_results_rcs <- testEstimates(qhat = qhat, uhat = uhat, var.comp = TRUE)
pooled_matrix_rcs <- as.data.frame(pooled_results_rcs$estimates)

# 4. Compute 95% CI manually
t_crit <- qt(0.975, df = pooled_matrix_rcs$df)
lower  <- pooled_matrix_rcs$Estimate - t_crit * pooled_matrix_rcs$Std.Error
upper  <- pooled_matrix_rcs$Estimate + t_crit * pooled_matrix_rcs$Std.Error

# 5. Combine into final table
pooled_matrix_rcs_final <- cbind(
  Estimate  = pooled_matrix_rcs$Estimate,
  Std.Error = pooled_matrix_rcs$Std.Error,
  `2.5 %`   = lower,
  `97.5 %`  = upper
)

rownames(pooled_matrix_rcs_final) <- c("B[1] (spline 1)", "B[2] (spline 2)", "rho")

# 6. Display results
print(round(pooled_matrix_rcs_final, 4))

######### Graph

### make a graph to compare them
load("Saved_Data/Model4_RCS.RData")

# 1. Define the RCS knot locations
knots <- quantile(
  syn_data$timing_DOAC,
  probs = c(0.1, 0.5, 0.9),
  na.rm = TRUE
)

# 2. Build global x range
x_seq_overall <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

# 3. Compute RCS basis (2 columns)
rcs_overall <- rcspline.eval(x_seq_overall, knots = knots, inclx = TRUE)

# 4. Extract from complete-case JAGS model
#    a) Pooled slope estimates
B1_jags <- Model4_RCS$BUGSoutput$summary["B[1]", "mean"]
B2_jags <- Model4_RCS$BUGSoutput$summary["B[2]", "mean"]

#    b) Average intercept (mean of u[i])
u_summaries <- Model4_RCS$BUGSoutput$summary[
  grep("^u\\[", rownames(Model4_RCS$BUGSoutput$summary)), "mean"
]
alpha_mean <- mean(u_summaries)

#    c) Compute predicted probabilities (JAGS)
logit_jags <- alpha_mean +
  B1_jags * rcs_overall[, 1] +
  B2_jags * rcs_overall[, 2]
risk_jags <- plogis(logit_jags)

df_jags <- data.frame(
  x     = x_seq_overall,
  risk  = risk_jags,
  model = "Bayesian pooled spline"
)

# 5. Extract from Rubin-pooled results
B1_rubin <- pooled_matrix_rcs_final["B[1] (spline 1)", "Estimate"]
B2_rubin <- pooled_matrix_rcs_final["B[2] (spline 2)", "Estimate"]

logit_rubin <- alpha_mean +
  B1_rubin * rcs_overall[, 1] +
  B2_rubin * rcs_overall[, 2]
risk_rubin <- plogis(logit_rubin)

df_rubin <- data.frame(
  x     = x_seq_overall,
  risk  = risk_rubin,
  model = "Rubin pooled spline"
)

# 6. Combine both lines
df_both <- rbind(df_jags, df_rubin)

# 7. Plot: Only the two pooled spline lines on log y-axis
RCS_pooled_comparison_plot <- ggplot(df_both, aes(x = x, y = risk, color = model, linetype = model)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "Bayesian pooled spline" = "blue",
    "Rubin pooled spline"    = "black"
  )) +
  scale_linetype_manual(values = c(
    "Bayesian pooled spline" = "dashed",
    "Rubin pooled spline"    = "dotdash"
  )) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.03, 0.05, 0.07, 0.1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Timing in Days (ref = 0)",
    y = "Predicted Risk (log scale)",
    color = "Model",
    linetype = "Model",
    title = "Comparison of Pooled RCS Timing-Risk Slopes"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top") +
  xlim(0, 25)

# 8. Print the plot
RCS_pooled_comparison_plot

save(pooled_matrix_rcs_final, RCS_pooled_comparison_plot,  model_results_rcs, file="Saved_Data/Model3_RCS_MI.RData")

