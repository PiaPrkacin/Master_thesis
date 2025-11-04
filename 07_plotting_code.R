#07_plotting_code
##############################
#Author: Pia Prkacin
#Date: 29/10/2025
#R version: 4.2.2
#Output
# - Fig.2
# - Fig.3
##############################
# Load packages
library(tidyverse)
library(scales)
library(ggtext)
library(tagger)
##############################
# Individual data for Figure 2
##############################
# Load data
load("Individual data for Figure2.RData")

# Prepare simulation data
data_toplot_ind <- data_toplot_ind %>%
  mutate(
    Immunity = factor(Immunity, levels = c("Absent", "Mild", "Strong", "Very Strong")),
    decay_strength = factor(decay_strength, levels = c("Absent", "Mild", "Strong", "ICL")),
    Exposure = "Water-contact-based",
    # Create a combined grouping variable for unique lines
    line_group = interaction(decay_strength, Snails, Immunity)
  )

# Filter for moderate endemicity only
data_toplot <- data_toplot_ind %>%
  filter(Endemicity == "Moderate", Exposure == "Water-contact-based")

# Filter simulation data
data_toplot_filtered <- data_toplot %>%
  filter(Snails %in% c("Absent", "Strong")) %>%
  # Exclude specific cases
  filter(
    !(decay_strength == "ICL" & Immunity %in% c("Mild", "Strong", "Very Strong") & Snails == "Absent"),
    !(decay_strength == "Strong" & Immunity == "Mild" & Snails == "Absent")
  ) %>%
  mutate(
    # Relabel Snails
    Snails = recode(Snails,
                    "Absent" = "Absent",
                    "Strong" = "Present"),
    # Relabel Immunity
    Immunity = recode(Immunity,
                      "Absent" = "Immunity: Absent",
                      "Mild" = "Immunity: Mild",
                      "Strong" = "Immunity: Strong",
                      "Very Strong" = "Immunity: Very strong")
  )

# Prepare Fulford field data (only Misumi village for moderate)
Fulford.data <- list.files(path = getwd(), pattern = "^Ageprofiles", full.names = TRUE) %>% 
  lapply(read_csv2, show_col_types = FALSE) %>%
  bind_rows()
colnames(Fulford.data) <- c("Age", "Eggs", "Village", "Endemicity")
Fulford.data <- Fulford.data %>%
  filter(Endemicity == "Moderate", Village == "Misuuni")

# Plot
eggs_plot <- ggplot() +
  geom_line(
    data = data_toplot_filtered,
    aes(
      x = avg_age_group,
      y = (geom_epg_mean - 1) * 24,
      group = line_group,
      colour = decay_strength,
      linetype = Snails
    ),
    linewidth = 1
  ) +
  
  # FIELD DATA
  geom_line(
    data = Fulford.data,
    aes(x = Age, y = Eggs),
    colour = "black", linewidth = 0.8, linetype = "solid"
  ) +
  
  # FACETS
  facet_grid(Immunity ~ ., scales = "free_y") +
  tag_facets(tag_levels = "A", position = "tr") +
  
  # SCALES
  scale_y_continuous(
    name = "Average eggs per gram of faeces \n(geometric mean)",
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_x_continuous(
    name = "\nAge [years]",
    breaks = seq(0, 80, 10),
    limits = c(0, 70),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    name = "Decay strength",
    values = c(
      "Absent" = "red",
      "Mild" = "purple",
      "Strong" = "blue",
      "ICL" = "darkgreen"
    ),
    labels = c(
      "Absent" = "Absent",
      "Mild" = "Mild",
      "Strong" = "Strong",
      "ICL" = "Very strong"
    )
  ) +
  scale_linetype_manual(
    name = "Snail Regulation",
    values = c("Absent" = "dotted", "Present" = "solid"),
    guide = guide_legend(override.aes = list(colour = "black"))
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5, 10, 0, 10, "pt"),
    panel.spacing = unit(1.5, "lines"),
    legend.key.width = unit(2, "lines"),
    strip.text = element_text(size = 14),
    tagger.panel.tag.text = element_text(size = 14),
    tagger.panel.tag.background = element_blank(),
    legend.box = "vertical"
  ) +
  
  # Annotation
  annotate(
    "text", x = 60, y = max(Fulford.data$Eggs) * 0.9, 
    label = "Field Data (Misuuni)", size = 5
  )

# Save
dir.create("Plots/Manuscript", recursive = TRUE, showWarnings = FALSE)
tiff("Plots/Manuscript/Fig2.png", compression = "lzw",
     width = 12, height = 12, units = "in", res = 300)
print(eggs_plot)
dev.off()


### Fig.3 (final version with facet letters + light grey strips) - uses tagger::tag_facets()
##############################
# Load parameters & functions
source("01_Handy_functions.R")
source("02_Parameters_Smansoni.R")

# Set MDA parameters
parms$mda <- list(
  age.lo = 5,  # SAC is 5–15
  age.hi = 15,
  start = 150,
  end = 159,
  frequency = 1,        # annual
  coverage = 0.75,
  fr_excluded = 0.05,   # systematic non-compliance 
  efficacy = 0.86
)

# Load collated population-level data
load("Population data for Figure3.RData")

# Packages
library(tidyverse)
library(scales)
library(ggtext)

if (!requireNamespace("tagger", quietly = TRUE)) install.packages("tagger")
library(tagger)

res %>%
  summarise(
    Endemicity_vals = paste(unique(Endemicity), collapse = ", "),
    Snails_vals     = paste(unique(Snails), collapse = ", "),
    DDF_vals        = paste(unique(DDF), collapse = ", "),
    Exposure_vals   = paste(unique(Exposure), collapse = ", "),
    Immunity_vals   = paste(unique(Immunity), collapse = ", "),
    Decay_vals      = paste(unique(decay_strength), collapse = ", ")
  ) %>% print()

# Filter dataset
res2 <- res %>%
  filter(
    Endemicity == "Moderate",
    Snails %in% c("Absent", "Strong"),
    DDF == "Absent",
    Exposure == "Based on water contacts",
    Immunity %in% c("Absent", "Mild", "Strong", "Very Strong"),
    decay_strength %in% c("Absent", "Mild", "Strong", "ICL")
  )

# Aggregate population data
data_avg2 <- res2 %>%
  group_by(time, Immunity, Snails, decay_strength, Endemicity) %>%
  summarise(
    eggs_prev_SAC = mean(eggs_prev_SAC, na.rm = TRUE),
    .groups = "drop"
  )

# Recode labels for facets/legend
data_avg2 <- data_avg2 %>%
  mutate(
    Snails = recode(Snails,
                    "Absent" = "Snail regulation: Absent",
                    "Strong" = "Snail regulation: Present"),
    Immunity = recode(Immunity,
                      "Absent"      = "Immunity: Absent",
                      "Mild"        = "Immunity: Mild",
                      "Strong"      = "Immunity: Strong",
                      "Very Strong" = "Immunity: Very strong"),
    decay_strength = recode(decay_strength,
                            "ICL"    = "Very strong",
                            "Absent" = "Absent",
                            "Mild"   = "Mild",
                            "Strong" = "Strong")
  )

# Remove combinations that did not reach equilibrium
data_avg2 <- data_avg2 %>%
  filter(!(
    (Snails == "Snail regulation: Absent" &
       decay_strength == "Very strong" &
       Immunity %in% c("Immunity: Mild", "Immunity: Strong", "Immunity: Very strong")) |
      (Snails == "Snail regulation: Absent" &
         decay_strength == "Strong" &
         Immunity == "Immunity: Mild")
  ))

# Build Fig3
Fig3 <- data_avg2 %>%
  ggplot(aes(x = time / 12, y = eggs_prev_SAC * 100, color = decay_strength)) +
  geom_line(aes(group = interaction(Immunity, Snails, decay_strength)), linewidth = 1) +
  facet_grid(
    Snails ~ Immunity,
    labeller = label_value
  ) +
  tag_facets(tag_levels = "A", position = "tr") +   # from tagger package
  scale_y_continuous(
    name = "Probability of elimination as a public health problem (%)",
    expand = expansion(mult = c(0, 0), add = c(0, 5))
  ) +
  scale_x_continuous(
    name = "\nYears since last treatment round",
    breaks = seq(parms$mda$end - 10, parms$mda$end + 20, by = 10),
    labels = seq(-10, 20, by = 10),
    expand = c(0, 0)
  ) +
  coord_cartesian(xlim = c(parms$mda$end - 10, parms$mda$end + 20)) +
  expand_limits(x = 0, y = 0) +
  scale_color_manual(
    name = "Decay strength",
    values = c(
      "Absent"      = "red",
      "Mild"        = "purple",
      "Strong"      = "blue",
      "Very strong" = "darkgreen"
    )
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(5, 10, 0, 10, "pt"),
    panel.spacing = unit(1.5, "lines"),
    legend.key.width = unit(2, "lines"),
    strip.text = element_text(size = 13),
    # light grey facet strips with black border
    strip.background = element_rect(fill = "gray90", colour = "black", size = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    tagger.panel.tag.text = element_text(size = 14, face = "bold"),
    tagger.panel.tag.background = element_blank()
  )

# Display
print(Fig3)

# Save
dir.create("Plots/Manuscript/Rebuttal", recursive = TRUE, showWarnings = FALSE)
tiff(
  filename = "Plots/Manuscript/Rebuttal/Fig3true.png",
  compression = "lzw",
  width = 13, height = 12, units = "in", res = 300
)
print(Fig3)
dev.off()



