#05.run_model
#############################
#Author: Pia Prkacin
#Date: 29/10/2025
#R version: 4.2.2
#This is the CONSOLE SCRIPT. 
#It is the main script, which allows to load input data, call and customise the parameters for simulations,
#to launch the SchiSTOP model from R source (04_Model_specification.R) and store results at two levels:
# - population output in the folder /Output/Population, i.e. results aggregated at population level (e.g. prevalence)
# - individual output in the folder /Output/Individual, i.e. results tracked over time for each individuals
#
#Input data: 
# - "Equilibrium_age_distribution.RData" generated from script 00_Demography.R
# - "prob_death_Uganda_2019.csv" (age-specific death probabilities)
#Loaded scripts:
# - 01_Handy_functions.R, list of functions used to define the transmission model
# - 01.1_Initial_conditions.R, initializes the population and the conditions to start the model
# - 02_Parameters_Smansoni.R, a full list of parameters employed for the simulations
# - 02.2_Setting_simulation_scenario.R, allows customization of specific settings for: scenarios of regulating mechanism, endemicity, age-exposure function 
# - 04_Model_specification.R, the actual definition of SchiSTOP
#
#Output:
# multiple .RDS files are created in the Output folder:
# - population output in the folder /Output/Population
# - individual output in the folder /Output/Individual
#############################

#Loading packages
#line below only needed if R library is different than default
#.libPaths(c("C:/Program Files/R/R-4.1.2/library",.libPaths()))
library(tidyverse)
library(readxl)
library(ggridges)
library(foreach)
library(doParallel)
library(patchwork)
library(rstudioapi)
library(beepr)
library(dplyr)
library(ggplot2)

################
# Setting source and output directories
################
source.dir <- dirname(getActiveDocumentContext()$path)
setwd(source.dir)
output.dir <- file.path(source.dir, "/Output")
if(!file.exists(output.dir)){
  dir.create(output.dir)
}

################
#Loadings
################

#Load initial population and conditions
source("01.1_Initial_conditions.R")

#0=male, 1=female
#Check initial human cohort
hist(cohort$age, breaks = c(0, prob_death$Age_hi[-c(1, nrow(prob_death))]+1),
     main = "Initial age distribution", xlab = "Age (ys)")
table(cohort$sex)

#Load functions
source("01_Handy_functions.R")

##Load parameters
source("02_Parameters_Smansoni.R")

################
#Initializing simulations
################

#Worms are initialized in the cohort with an artificial FOI (1 worm pair per person per month)
init$humans$pop <- init$humans$pop %>%
    mutate(Ind_sus = rgamma(nrow(cohort), shape = parms$parasite$k_w, scale = 1/parms$parasite$k_w), #individual susceptibility to infection
           complier = as.numeric(rbernoulli(nrow(cohort), 1-parms$mda$fr_excluded))) #complier factor for participation in MDA

###############################
#Preparing simulation settings for the specific regulating/exposure assumptions
#Please customise this file in advance
###############################

source("02.2_Setting_simulation_scenario.R")

################
#Run the model
################
setwd(source.dir)

time.start <- Sys.time()
source("04_Model_specification.R")
time.end <- Sys.time()
time.end - time.start
#beep()

################
#Collating and saving population-level results 
################
#Population-level results
#Set endemicity setting:
pop.output.dir <- file.path(output.dir, "Population")
if(!file.exists(pop.output.dir)){
  dir.create(pop.output.dir)
}
#Collating and saving population-level output
#Individual output is automatically saved throughout the simulations (in 04_Model_specification.R), unless silenced in the simulation settings
res <- bind_rows(results)
saveRDS(bind_rows(results), file = file.path(pop.output.dir, 
                           paste(setting, ".RDS", sep = "")))

##############################
#Prevalence plot 
##############################

#Prepare summary data
data_avg <- res %>%
  group_by(time, Immunity, Snails, decay_strength, Endemicity) %>%
  summarise(
    eggs_prev_SAC = mean(eggs_prev_SAC, na.rm = TRUE),
    eggs_prev_tot = mean(eggs_prev, na.rm = TRUE),
    PHI           = mean(Heggs_prev_SAC, na.rm = TRUE),
    snail_inf     = mean(inf_snail, na.rm = TRUE),
    snail_exp     = mean(exp_snail, na.rm = TRUE),
    susc_snail    = mean(susc_snail, na.rm = TRUE),
    snail_prev    = mean(inf_snail / (susc_snail + inf_snail + exp_snail), na.rm = TRUE)
  ) %>%
  ungroup()

data_avg <- data_avg %>%
  filter(Snails %in% c("Absent", "Strong")) %>%   # Keep only these two
  mutate(
    Snails = ifelse(Snails == "Strong", "Present", "Absent"),
    Snails = factor(Snails, levels = c("Absent", "Present"))
  )

data_avg <- data_avg %>%
  mutate(
    decay_strength = case_when(
      decay_strength %in% c("ICL", "Very strong", "very strong", "icl") ~ "Very strong",
      TRUE ~ decay_strength
    ),
    decay_strength = factor(
      decay_strength,
      levels = c("Absent", "Mild", "Strong", "Very strong")
    )
  )

snail_labels <- c(
  "Absent"  = "Snail regulation: Absent",
  "Present" = "Snail regulation: Present"
)

immunity_labels <- c(
  "Absent"      = "Immunity: Absent",
  "Mild"        = "Immunity: Mild",
  "Strong"      = "Immunity: Strong",
  "Very Strong" = "Immunity: Very strong"
)

#Plotting code
data_avg %>%
  filter(Endemicity == "Moderate") %>%
  ggplot(aes(x = time / 12, group = interaction(decay_strength, Snails, Immunity))) + 
  geom_line(aes(y = eggs_prev_SAC * 100, colour = decay_strength), alpha = 0.6) + 
  geom_hline(yintercept = 30, linetype = "dashed", color = "gray30") +
  geom_line(aes(y = PHI * 100, colour = decay_strength), linetype = "longdash") +
  geom_line(aes(y = snail_prev * 100, colour = decay_strength), linetype = "dotted") +
  facet_grid(
    Snails ~ Immunity,
    labeller = labeller(
      Snails = snail_labels,
      Immunity = immunity_labels
    )
  ) +
  scale_y_continuous(
    name = "Prevalence of infection in SAC (%)\n",
    breaks = seq(0, 100, 10),
    limits = c(0, 50),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    name = "\nTime [Years]",
    expand = c(0, 0)
  ) +
  expand_limits(x = 0, y = 0) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5, 10, 0, 10, "pt"),
    panel.spacing = unit(1.5, "lines")
  ) +
  labs(colour = "Decay strength") +
  scale_color_manual(
    name = "Decay strength",
    values = c(
      "Absent"      = "red",
      "Mild"        = "darkgreen",
      "Strong"      = "blue",
      "Very strong" = "purple"
    )
  )


##############################
#Equilibrium check 
##############################

data_avg <- res %>%
  group_by(time, Immunity, Snails, decay_strength, Endemicity) %>%
  summarise(
    eggs_prev_SAC = mean(eggs_prev_SAC, na.rm = TRUE),
    eggs_prev_tot = mean(eggs_prev, na.rm = TRUE),
    PHI = mean(Heggs_prev_SAC, na.rm = TRUE),
    snail_inf = mean(inf_snail, na.rm = TRUE),
    snail_exp = mean(exp_snail, na.rm = TRUE),
    susc_snail = mean(susc_snail, na.rm = TRUE),
    snail_prev = mean(inf_snail / (susc_snail + inf_snail + exp_snail), na.rm = TRUE)
  ) %>%
  ungroup()

df <- data_avg %>%
  mutate(year = time/12)

check_equilibrium <- function(data, var, tmin, tmax) {
  d <- filter(data, year >= tmin, year <= tmax)
  fit <- lm(reformulate("year", var), data = d)
  broom::tidy(fit) %>% filter(term == "year")
}

df %>%
  filter(Endemicity == "Moderate",
         Snails == "Absent") %>%
  group_by(decay_strength, Immunity) %>%
  group_modify(~ check_equilibrium(.x, "eggs_prev_SAC*100", 100, 150)) %>%
  arrange(desc(estimate))

