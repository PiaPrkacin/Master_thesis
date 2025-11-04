##02.2_setting_simulation_scenario
#############################
#Author: Pia Prkacin
#Date: 29/10/2025
#R version: 4.1.2
#Description: Customize simulation scenario-specific parameters. 
# Species of interest: Schistosoma mansoni
#############################

# Ensure required packages are loaded
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
library(readxl)

T <- 300 # number of simulated years
seeds <- 100 # stochastic seeds
fr <- 10 # frequency for printing to file the individual output [years]
write.output <- TRUE # disable individual output if not needed

# MDA settings
parms$mda <- list(
  age.lo = 5, 
  age.hi = 15,
  start = 150,
  end = 159,
  frequency = 1,
  coverage = 0.75,
  fr_excluded = 0.05,
  efficacy = 0.86
)

# Exposure setting: choose "Sow" or "ICL"
exposure <- "Sow"

# Define stochastic scenarios
stoch_scenarios <- expand.grid(
  seed = 1:seeds,
  DDF_strength = c("Absent"),
  imm_strength = c("Absent", "Mild", "Strong", "Very Strong"),
  snails = c("Absent", "Mild", "Strong"),
  endem = c("Moderate"),
  decay_strength = c("Absent", "Mild", "Strong", "ICL")
) %>%
  filter(!(imm_strength == "Absent" & decay_strength != "Absent")) # remove invalid combinations

# Load and prepare zeta table
# Drop extra rows per scenario from zetas (keep only 1 row per combination)
zetas <- read_excel(paste0("Zetas_", exposure, "thesis.xlsx")) %>%
  filter(Endemicity == "Moderate") %>%
  rename(
    endem = Endemicity,
    imm_strength = Immunity,
    snails = Snails,
    DDF_strength = DDF,
    decay_strength = Decay,
    zeta = Zeta_grid_search,
    worms_aggr = Kw,
    tr_snails = `Transmission on snails`,
    equilibrium = Equilibrium
  ) %>%
  mutate(across(c(endem, imm_strength, snails, decay_strength), ~ trimws(as.character(.)))) %>%
  distinct(endem, imm_strength, snails, decay_strength, .keep_all = TRUE)

# Clean stoch_scenarios
stoch_scenarios <- expand.grid(
  seed = 1:seeds,
  DDF_strength = c("Absent"),
  imm_strength = c("Absent", "Mild", "Strong", "Very Strong"),
  snails = c("Absent", "Mild", "Strong"),
  endem = c("Moderate"),
  decay_strength = c("Absent", "Mild", "Strong", "ICL")
) %>%
  filter(!(imm_strength == "Absent" & decay_strength != "Absent")) %>%
  mutate(across(c(endem, imm_strength, snails, DDF_strength, decay_strength), ~ trimws(as.character(.))))

# Join zetas (now includes decay_strength)
stoch_scenarios <- stoch_scenarios %>%
  left_join(zetas, by = c("endem", "imm_strength", "snails", "DDF_strength", "decay_strength")) %>%
  filter(equilibrium == TRUE) %>%
  mutate(
    Ext_foi_value = case_when(
      endem == "Low" ~ 0.5,
      endem == "Moderate" ~ 1,
      endem == "High" ~ 5
    ),
    Ext_foi_duration = case_when(
      endem == "Low" ~ 0.5,
      endem == "Moderate" ~ 2,
      endem == "High" ~ 2
    ),
    decay_rate = case_when(
      decay_strength == "Absent" ~ 0,
      decay_strength == "Mild" ~ 0.00576
      decay_strength == "Strong" ~ 0.0561
      decay_strength == "ICL" ~ 0.467
    )
  )

stoch_scenarios %>%
  filter(imm_strength == "Mild", snails == "Strong") %>%
  distinct(zeta, worms_aggr, tr_snails, decay_strength) %>%
  arrange(decay_strength)


# Define the title of the simulation setting. This will give the name to all the output files.
setting <- paste(exposure, "func_SACmda_ModerateOnly", sep = "_") 

################
#Set output directory to save individual results
################
write.output <- TRUE

source.dir <- getwd()

if (write.output == TRUE) {
  ind.output.dir <- file.path(source.dir, "Output", "Individual", setting)
  
  if (!dir.exists(ind.output.dir)) {
    dir.create(ind.output.dir, recursive = TRUE)
    message("Created output directory: ", ind.output.dir)
  } else {
    message("Output directory already exists: ", ind.output.dir)
  }
}
