#* ***********************************#
#* Study: CATALYST IPD meta-analysis - time-response 
#* Purpose: Master File
#* Author: Konstantina Chalkou
#* Date created: 15.11.2024
#* Last update: 31.10.2025
#* **********************************#
rm(list = ls())
# Load required libraries

library(dplyr)
library(R2jags)
library(coda)
library(dplyr)
library(writexl)
library(ggplot2)
library(rms)
library(Hmisc)
library(tidyr)
library(rjags)
library(readxl)
library(jomo)
library(ggpubr)
library(mitml)  # for testEstimates
library(splines)
library(lmtest)
library(sandwich)
library(nnet)
library(scales)
library(patchwork)
library(stringr)
library(ggnewscale)

#########
### A. Load the synthetic data
load("synthetic_catalyst.RData")

#########################################################################################3
################ Primary outcome: composite at 30 days ######################
##################################################################################

### Model 1: Linear as ransomized meta-analysis
source("01a_Model1_Linear_as_randomised.R")
### Tables 
# Meta-analysis results
Model1
## Figures
# independent before meta-analysis figure
Model1_separate_plot_logy
# random-effects - after the meta-analysis figure
Model1_plot_log

### Model 2: Non-linear as randomized Meta-Analysis
 
load("synthetic_catalyst.RData")
source("02a_Model2_NonLinear_as_randomised.R")
### Tables 
# Meta-analysis results
Model2_MVN
## Figures
# random-effects - after the meta-analysis figure
ModelRCS2_plot_log

### Model 3: Linear as observed meta-analysis
 
load("synthetic_catalyst.RData")
source("03a_Model3_Linear_as_observed.R")
### Tables 
# Meta-analysis results
Model3_Linear
## Figures
# independent before meta-analysis figure
Model3_separate_plot_logy
# random-effects - after the meta-analysis figure
Model3_plot_log

### Model 4: Non Linear as observed meta-analysis
 
load("synthetic_catalyst.RData")
source("04a_Model4_NonLinear_as_observed.R")
### Tables 
# Meta-analysis results
Model4_RCS
## Figures
# random-effects - after the meta-analysis figure
Model4_RCS_plot_log

#### MULTIPLE IMPUTATIONS
### Model 3: Linear as observed meta-analysis with multiple imputations
 
load("synthetic_catalyst.RData")
source("5a_Model3_Linear_as_observed_Multiple_imputations.R")

#Results
pooled_matrix_final

### Model 4: Non Linear as observed meta-analysis with multiple imputations
 
load("synthetic_catalyst.RData")
source("5b_Model4_NonLinear_as_observed_Multiple_imputations.R")

#Results
pooled_matrix_rcs_final

#########################################################################################3
################ Secondary outcome: composite at 90 days ######################
##################################################################################


### Model 1: Linear as ransomized meta-analysis
load("synthetic_catalyst.RData")
source("01b_Model1_Linear_as_randomised_90days.R")
### Tables 
# Meta-analysis results
Model1_90
## Figures
# independent before meta-analysis figure
Model1_90_separate_plot_logy
# random-effects - after the meta-analysis figure
Model1_90_plot_log

### Model 2: Non-linear as randomized Meta-Analysis

load("synthetic_catalyst.RData")
source("02b_Model2_NonLinear_as_randomised_90days.R")
### Tables 
# Meta-analysis results
Model2_90_MVN
## Figures
# random-effects - after the meta-analysis figure
ModelRCS2_90_plot_log

### Model 3: Linear as observed meta-analysis
load("synthetic_catalyst.RData")
source("03b_Model3_Linear_as_observed_90days.R")
### Tables 
# Meta-analysis results
Model3_90_Linear
## Figures
# independent before meta-analysis figure
Model3_90_separate_plot_logy
# random-effects - after the meta-analysis figure
Model3_90_plot_log

### Model 4: Non Linear as observed meta-analysis
load("synthetic_catalyst.RData")
source("04b_Model4_NonLinear_as_observed_90days.R")
### Tables 
# Meta-analysis results
Model4_90_RCS
## Figures
# random-effects - after the meta-analysis figure
Model4_90_RCS_plot_log

#########################################################################################3
################ Secondary outcome: Ischemic stroke at 30 days ######################
##################################################################################

### Model 1: Linear as ransomized meta-analysis
load("synthetic_catalyst.RData")
source("01c_Model1_Linear_as_randomised_RIS.R")
### Tables 
# Meta-analysis results
Model1_isch30
## Figures
# independent before meta-analysis figure
Model1_isch_separate_plot_logy
# random-effects - after the meta-analysis figure
Model1_isch30_plot_log

### Model 2: Non-linear as randomized Meta-Analysis

load("synthetic_catalyst.RData")
source("02c_Model2_NonLinear_as_randomised_RIS.R")
### Tables 
# Meta-analysis results
Model2_isch_MVN
## Figures
# random-effects - after the meta-analysis figure
ModelRCS2_isch_plot_log

### Model 3: Linear as observed meta-analysis
load("synthetic_catalyst.RData")
source("03c_Model3_Linear_as_observed_RIS.R")
### Tables 
# Meta-analysis results
Model3_RIS_Linear
## Figures
# independent before meta-analysis figure
Model3_RIS_separate_plot_logy
# random-effects - after the meta-analysis figure
Model3_RIS_plot_log

### Model 4: Non Linear as observed meta-analysis
load("synthetic_catalyst.RData")
source("04c_Model4_NonLinear_as_observed_RIS.R")
### Tables 
# Meta-analysis results
Model4_RIS_RCS
## Figures
# random-effects - after the meta-analysis figure
Model4_RIS_RCS_plot_log


#########################################################################################3
################ Secondary outcome: sICH at 30 days ######################
##################################################################################

### Model 1: Linear as ransomized meta-analysis
load("synthetic_catalyst.RData")
source("01d_Model1_Linear_as_randomised_sICH.R")
### Tables 
# Meta-analysis results
Model1_sich
## Figures
# independent before meta-analysis figure
Model1_sich30_separate_plot_logy
# random-effects - after the meta-analysis figure
Model1_sich_plot_log

### Model 2: Non-linear as randomized Meta-Analysis

load("synthetic_catalyst.RData")
source("02d_Model2_NonLinear_as_randomised_sICH.R")
### Tables 
# Meta-analysis results
Model2_sich_MVN
## Figures
# random-effects - after the meta-analysis figure
ModelRCS2_sich_plot_log

### Model 3: Linear as observed meta-analysis
load("synthetic_catalyst.RData")
source("03d_Model3_Linear_as_observed_sICH.R")
### Tables 
# Meta-analysis results
Model3_SICH
## Figures
# independent before meta-analysis figure
Model3_SICH_separate_plot_logy
# random-effects - after the meta-analysis figure
Model3_SICH_plot_log

### Model 4: Non Linear as observed meta-analysis
load("synthetic_catalyst.RData")
source("04d_Model4_NonLinear_as_observed_sICH.R")
### Tables 
# Meta-analysis results
Model4_SICH_RCS
## Figures
# random-effects - after the meta-analysis figure
Model4_SICH_RCS_plot_log

#### Analysis for checking the descrepancies of randomised vs observed and its association with NIHSS
 
load("synthetic_catalyst.RData")
source("06_TimeDeviation_nihss.R")

# Results
# A. Multinomial logistic regression model
tbl
# B. linear model
tbl2
# graphs
p
p2

#### Figures of the main manuscript
 
load("synthetic_catalyst.RData")
source("Figures_Manuscript.R")
## Figure 1
violin_plot
## Figure 2
final_2x2_with_legend

#### Tables of the main manuscript
load("synthetic_catalyst.RData")
source("Tablesmanuscript.R")
report_tbl

### For the Appendix
load("synthetic_catalyst.RData")
source("Appendix Figures.R")

# Then the Appendix.Rmd provides all the Appendix material

rm(list = ls())

