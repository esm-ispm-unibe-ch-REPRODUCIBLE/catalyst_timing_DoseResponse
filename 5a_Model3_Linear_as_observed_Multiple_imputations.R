#####################################################################################################
##### Model 3: Multiple imputations Linear as observed Meta-Analysis #################
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

#### Run the Bayesian model for all the imputed datasets

# Store model results
model_results_list <- list()

modelLinearTiming <- function() {
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


# Loop through each imputed dataset
for (i in 1:10) {
  cat(sprintf("Running JAGS on imputed dataset %d...\n", i))
  imp_data <- imp.list[[i]]
  
  # Make sure composite_30 is binary numeric (0/1) for JAGS
  imp_data$composite_30 <- as.numeric(as.character(imp_data$composite_30))
  
  # JAGS input
  jags_data <- list(
    Np = nrow(imp_data),
    Nstudies = length(unique(imp_data$trial)),
    trial = as.numeric(as.factor(imp_data$trial)),  # ensure numeric IDs
    composite_30 = imp_data$composite_30,
    adjusted_timing = imp_data$adjusted_timing
  )
  
  # Fit JAGS model
  set.seed(1294821 + i)  # Different seed per chain
  fit <- jags.parallel(
    data = jags_data,
    inits = NULL,
    parameters.to.save = c("B", "tau", "tau2", "u", "beta"),
    model.file = modelLinearTiming,
    n.chains = 2,
    n.iter = 20000,
    n.burnin = 2000,
    DIC = FALSE,
    n.thin = 10
  )
  
  # Save results
  model_results_list[[i]] <- fit
}


# 1. Extract point estimates (means) and standard errors for B and tau
B_means <- numeric(10)
B_ses   <- numeric(10)

tau_means <- numeric(10)
tau_ses   <- numeric(10)

for (i in 1:10) {
  B_means[i]     <- model_results_list[[i]]$BUGSoutput$summary["B", "mean"]
  B_ses[i]       <- model_results_list[[i]]$BUGSoutput$summary["B", "sd"]
  
  tau_means[i]   <- model_results_list[[i]]$BUGSoutput$summary["tau2", "mean"]
  tau_ses[i]     <- model_results_list[[i]]$BUGSoutput$summary["tau2", "sd"]
}

# 2. Combine estimates
qhat <- rbind(B = B_means, tau = tau_means)
uhat <- rbind(B = B_ses^2, tau = tau_ses^2)

# 3. Run Rubin's rules pooling
pooled_results <- testEstimates(qhat = qhat, uhat = uhat, var.comp = TRUE)
pooled_matrix <- as.data.frame(pooled_results$estimates)

# 4. Compute 95% CI manually
t_crit <- qt(0.975, df = pooled_matrix$df)
lower  <- pooled_matrix$Estimate - t_crit * pooled_matrix$Std.Error
upper  <- pooled_matrix$Estimate + t_crit * pooled_matrix$Std.Error

# 5. Add to final output
pooled_matrix_final <- cbind(
  Estimate  = pooled_matrix$Estimate,
  Std.Error = pooled_matrix$Std.Error,
  `2.5 %`   = lower,
  `97.5 %`  = upper
)

pooled_matrix_final[2,3]<-0
rownames(pooled_matrix_final) <- c("B (slope)", "tau2 (SD of study slopes)")

# 6. Display
print(round(pooled_matrix_final, 4))
### make a graph to compare them
load("Saved_Data/Model3_Linear.RData")

# 1. Extract slope from JAGS model
B_jags <- Model3_Linear$BUGSoutput$summary["B", "mean"]

# 2. Average intercept from study-specific intercepts u[1], ..., u[N]
u_summaries <- Model3_Linear$BUGSoutput$summary[
  grep("^u\\[", rownames(Model3_Linear$BUGSoutput$summary)), "mean"
]
alpha_mean <- mean(u_summaries)

# 3. Rubin-pooled slope (from Rubin’s rules)
B_rubin <- pooled_matrix_final["B (slope)", "Estimate"]

# 4. Create x-axis for both lines
x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

# 5. Compute predicted probabilities for both models
#    Logistic transformation of linear predictor: logit(p) = intercept + B * x
logit_jags  <- alpha_mean + B_jags  * x_seq_pooled
logit_rubin <- alpha_mean + B_rubin * x_seq_pooled

pred_jags  <- plogis(logit_jags)
pred_rubin <- plogis(logit_rubin)

# 6. Create data frames for ggplot
df_jags <- data.frame(
  x    = x_seq_pooled,
  risk = pred_jags,
  model = "Bayesian pooled slope"
)

df_rubin <- data.frame(
  x    = x_seq_pooled,
  risk = pred_rubin,
  model = "Rubin pooled slope"
)

df_both <- rbind(df_jags, df_rubin)

# 7. Plot: Only the two pooled lines
pooled_plot <- ggplot(df_both, aes(x = x, y = risk, linetype = model, color = model)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c(
    "Bayesian pooled slope" = "blue",
    "Rubin pooled slope"    = "black"
  )) +
  scale_linetype_manual(values = c(
    "Bayesian pooled slope" = "dashed",
    "Rubin pooled slope"    = "dotdash"
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
    title = "Comparison of Pooled Timing-Risk Slopes"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top") +
  xlim(0, 25)

# 8. Show the plot
pooled_plot

###B
# 1. Extract Rubin pooled slope from Rubin's output
B_rubin <- pooled_matrix_final["B (slope)", "Estimate"]

# 2. Use same x-axis as for Bayesian line
# If you don't have this already in environment:
x_seq_pooled <- seq(
  from = min(syn_data$timing_DOAC, na.rm = TRUE),
  to   = max(syn_data$timing_DOAC, na.rm = TRUE),
  length.out = 200
)

# 3. Use alpha_mean from JAGS intercepts
u_summaries <- Model3_Linear$BUGSoutput$summary[
  grep("^u\\[", rownames(Model3_Linear$BUGSoutput$summary)), "mean"
]
alpha_mean <- mean(u_summaries)

# 4. Compute Rubin-pooled predicted probabilities
logit_rubin <- alpha_mean + B_rubin * x_seq_pooled
pred_rubin  <- plogis(logit_rubin)

df_rubin_line <- data.frame(
  x    = x_seq_pooled,
  risk = pred_rubin
)

# 5. Add Rubin-pooled line to your existing full plot
Model3_plot_log_rubin <- Model3_plot_log +
  geom_line(
    data = df_rubin_line,
    aes(x = x, y = risk),
    color = "black",
    linetype = "dotdash",
    size = 1.2
  )

# 6. Display updated plot
Model3_plot_log_rubin

save(pooled_matrix_final, Model3_plot_log_rubin, pooled_plot,  model_results_list, file="Saved_Data/Model3_Linear_MI.RData")

