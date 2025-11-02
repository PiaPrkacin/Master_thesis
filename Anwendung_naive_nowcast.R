
source("Funktionen/nowcast_cleantable.R", encoding = "UTF-8")
source("Funktionen/get_naive_model.R", encoding = "UTF-8")
source("Funktionen/nowcast_naive_model.R", encoding = "UTF-8")
source("Funktionen/plot_naive_models.R", encoding = "UTF-8")

# Step 1: Clean and prepare the raw SARI data
df_raw <- nowcast_cleantable()

# Step 2: Create naive delay distribution model
naive_model <- get_naive_model(
  df_raw,
  start_week = "2024-W37",
  end_week = "2025-W25",
  pathogen_group = c("Covid-19", "Influenza"),
  age_groups = c("15 - 29", "30 - 44"),
  delay_weeks = 15
)

# Step 3: Apply naive nowcasting for a specific hospitalisation week
nowcast_results <- nowcast_naiv_model(
  df_raw,
  naive_model,
  target_hospitalisation_week = "2025-W27",
  pathogen_group = c("Covid-19", "Influenza"),
  age_groups = c("15 - 29", "30 - 44"),
  save_excel_path = "N:/MED/INF/DMOD/Public/Projekte/SARI_Nowcasting/nowcast_kw27.xlsx"
)

# Step 4: Table
print(nowcast_results)



# Step 5: Plot

plot_naive_model(nowcasted_values = nowcast_results,
                 title = "Naives Nowcasting für Hospitalisierungswoche KW 27",
                 subtitle = "Gemeldete Fälle und geschätzte Gesamtfallzahl mit Konfidenzintervall",
                 save_plot_path = "N:/MED/INF/DMOD/Public/Projekte/SARI_Nowcasting/nowcast_kw27_plot.png")























####different tries


print(nowcast_results)


plot_naive_model(nowcasted_values = nowcast_results,
                 title = "Naives Nowcasting für Hospitalisierungswoche KW 27",
                 subtitle = "Gemeldete Fälle und geschätzte Gesamtfallzahl mit Konfidenzintervall",
                 save_plot_path = "N:/MED/INF/DMOD/Public/Projekte/SARI_Nowcasting/nowcast_kw27_plot.png")




unique(df_raw$hospitalisation_week[!str_detect(df_raw$hospitalisation_week, "^\\d{4}-W\\d{2}$")])
