#06_Predictions
#############################
#Author: Pia Prkacin
#Date: 29/10/2025
#R version: 4.2.2
#
#This script loads population and individual outputs saved from 05_Run_model.R
#Data are cleaned, prepared, and saved for plotting.  
#This script also computes predictions of probability to reach control targets across different modelling scenarios.
#
# Input data:
# - population .RDS output in the folder /Output/Population
# - individual .RDS output in the folder /Output/Individual
#
# Loaded scripts:
# - 01_Handy_functions.R, list of functions used to define the transmission model
# - 02_Parameters_Smansoni.R, a full list of parameters employed for the simulations
# - 02.2_Setting_simulation_scenario.R, allows customization of specific settings for: 
#   scenarios of regulating mechanism, endemicity, age-exposure function 
#
# Output:
# - "Population data for Figure3.RData"
# - "Individual data for Figure2.RData"
# - table with predicted probabilities
##############################

#line below only needed if R library is different than default
#.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(patchwork)
library(egg)
library(rstudioapi)
library(scales)

#Install devtools
install.packages("devtools")
devtools::install_github("eliocamp/tagger")
library(tagger)
library(rempsyc)
library(rstudioapi)

#Set folder
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)

# Load functions
source("01_Handy_functions.R")

# Load parameters
source("02_Parameters_Smansoni.R")
source("02.2_Setting_simulation_scenario.R")

######### Setting #########

# Load population output
pop.output.dir <- file.path(source.dir, "Output/Population")

# Load only the "Sow" scenario
setting <- "Sow_func_SACmda_ModerateOnly"

# Moderate endemicity only
res <- readRDS(file.path(pop.output.dir, paste0(setting, ".RDS"))) %>%
  mutate(
    Exposure = factor("Based on water contacts", levels = "Based on water contacts"),
    Endemicity = factor(levels = "Moderate")
  )

gc()
getwd()
######### Check for faded-out runs and compute averages #########

# Compute number of seeds used (should be 100, since you set it in your simulation script)
seeds <- n.un(res$seed)
print(paste("Number of seeds:", seeds))


# Check for faded-out runs at pre-control (149 years is just before MDA starts)
# time is months, so (parms$mda$start - 1) * 12
faded <- res %>%
  filter(time == (parms$mda$start - 1) * 12) %>%
  group_by(Immunity, Snails, DDF) %>%
  summarise(
    faded_seed = length(which(eggs_prev_SAC == 0)) * 100 / seeds,
    .groups = "drop"
  )

n_faded <- length(which(faded$faded_seed > 0))
print(paste("Number of groups with faded-out runs:", n_faded))

#Average by seed and update results, # Remove faded-out runs if needed
if(n_faded>0){
  res <- res[- which(res$time==149*12 & res$eggs_prev_SAC==0), ]
}

# Average over seeds, keeping decay_strength and Immunity for later plotting
data_avg <- res %>%
  group_by(time, Immunity, Snails, DDF, decay_strength, Endemicity, Exposure) %>%
  summarise(
    eggs_prev_SAC = mean(eggs_prev_SAC),
    eggs_prev_tot = mean(eggs_prev),
    PHI = mean(Heggs_prev),
    miracidiae = mean(miracidiae),
    cercariae = mean(cercariae),
    snail_inf = mean(inf_snail),
    snail_exp = mean(exp_snail),
    snail_prev = mean(inf_snail / (susc_snail + inf_snail + exp_snail)),
    .groups = "drop"
  )

save(res, data_avg, file = "Population data for Figure3.RData")



##############################
# Individual data for Figure 2
##############################

library(dplyr)

# Define exposure and directories
exposure <- "Sow"
setting <- paste(exposure, "func_SACmda_ModerateOnly", sep = "_")
ind.output.dir <- file.path(source.dir, paste("Output/Individual/", setting, sep = ""))

# Load all individual simulation output files
file_paths <- list.files(path = ind.output.dir,
                         pattern = "^Ind_out_seed", full.names = TRUE)

# Keep only the pre-MDA (baseline) time point
data_all <- lapply(file_paths, function(file) {
  df <- readRDS(file)
  df[df$time == (parms$mda$start - 9), ]   # 9 months before MDA start
}) %>%
  bind_rows()

# Filter conditions: Moderate endemicity, DDF absent, Snails absent or strong
data_all <- data_all %>%
  filter(Endemicity == "Moderate",
         DDF == "Absent",
         Snails %in% c("Absent", "Mild", "Strong"))

# Compute age-structured averages (per seed)
age_out <- data_all %>%
  group_by(time, Immunity, Snails, DDF, decay_strength, Endemicity) %>%
  mutate(age_group = cut_number(age, n = 14)) %>%
  group_by(seed, time, age_group, Immunity, Snails, DDF, decay_strength, Endemicity) %>%
  summarise(
    avg_age_group = mean(age),
    epg = mean(ec),
    geom_epg = geom_mean(ec + 1),
    wp = mean(tot_wp),
    dwp = mean(cum_dwp),
    rate = mean(rate),
    .groups = "drop"
  )

# Average over seeds and compute confidence intervals
data_toplot_ind <- age_out %>%
  group_by(age_group, time, Immunity, Snails, DDF, decay_strength, Endemicity) %>%
  summarise(
    avg_age_group = mean(avg_age_group),
    epg_mean = mean(epg),
    epg_lo = ci(epg)[1],
    epg_hi = ci(epg)[2],
    geom_epg_mean = mean(geom_epg),
    geom_epg_lo = ci(geom_epg)[1],
    geom_epg_hi = ci(geom_epg)[2],
    wp_mean = mean(wp),
    wp_lo = ci(wp)[1],
    wp_hi = ci(wp)[2],
    dwp_mean = mean(dwp),
    dwp_lo = ci(dwp)[1],
    dwp_hi = ci(dwp)[2],
    rate_mean = mean(rate),
    rate_lo = ci(rate)[1],
    rate_hi = ci(rate)[2],
    .groups = "drop"
  ) %>%
  mutate(Exposure = exposure)

# Save pre-MDA data for Figure 2 plotting
save(data_toplot_ind, file = "Individual data for Figure2.RData")

# Clean up memory
gc()


###############################################################
# Predictions of probabilities to reach the control targets
###############################################################

# Load your filtered output (Sow exposure, Moderate endemicity)
setting <- paste("Sow", "func_SACmda_ModerateOnly", sep = "_")
res_Sow <- readRDS(file.path(pop.output.dir, paste(setting, ".RDS", sep = ""))) %>%
  filter(
    Endemicity == "Moderate",                      
    Snails %in% c("Absent", "Mild", "Strong"),              
    DDF == "Absent"                                 
  ) %>%
  mutate(Exposure = "Based on water contacts")

# Ensure correct factor levels
res_Sow$Endemicity <- factor(res_Sow$Endemicity,
                             levels = "Moderate")

# Compute target probabilities at 20 years post-MDA
options(digits = 1)
pred <- res_Sow %>%
  filter(time == (parms$mda$end + 20) * 12) %>%
  group_by(Snails, Immunity, DDF, decay_strength, Exposure, Endemicity) %>%
  summarise(
    ephp_seed = length(which(Heggs_prev_SAC <= 0.01)) * 100 / seeds,
    eliminated = length(which(eggs_prev == 0)) * 100 / seeds,
    .groups = "drop"
  )

nice_table <- function(df, digits = 1) {
  df[] <- lapply(df, function(x) if(is.numeric(x)) round(x, digits) else x)
  return(df)
}
prob_table <- nice_table(pred)
print(table)
