# ==============================================================================
# Script:  calculate_ipcc_ar6_slr_rates.R
# Purpose: Extract and calculate 2040–2060 SLR rates (mm/yr) for each SSP scenario
#          matching NASA IPCC AR6 projection tool outputs
# Inputs:  - data/supporting data/ipcc_ar6_sea_level_projection_psmsl_id_523.csv
# Outputs: - data/supporting data/IPCC_AR6_SLR_Rates_2040_2060.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------
# 1. FILE PATHS & DATA INGESTION
# ------------------------------------------------------------------------------

# Access CSV via here() inside data/supporting data
input_csv <- here("data", "supporting data", "ipcc_ar6_sea_level_projection_psmsl_id_523.csv")

# Fallback check in case the file is in the root directory
if (!file.exists(input_csv)) {
  input_csv <- here("ipcc_ar6_sea_level_projection_psmsl_id_523.csv")
}

stopifnot("Input projection CSV not found!" = file.exists(input_csv))

# Output path aligned with header specification
out_dir <- here("data", "supporting data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_csv <- file.path(out_dir, "IPCC_AR6_SLR_Rates_2040_2060.csv")

df_raw <- read_csv(input_csv, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. FILTER & FORMAT 2040–2060 RATES ACROSS SSP SCENARIOS
# ------------------------------------------------------------------------------

# Filter for total rates under medium confidence (the IPCC standard reference)
# Column '2050' represents the 2040–2060 central difference rate in mm/yr
slr_rates_2050 <- df_raw %>%
  filter(
    process    == "totalrates",
    confidence == "medium"
  ) %>%
  select(scenario, quantile, rate_2050 = `2050`) %>%
  pivot_wider(
    names_from   = quantile,
    values_from  = rate_2050,
    names_prefix = "q"
  ) %>%
  mutate(
    Scenario_Label = case_match(
      scenario,
      "ssp119" ~ "SSP1-1.9",
      "ssp126" ~ "SSP1-2.6",
      "ssp245" ~ "SSP2-4.5",
      "ssp370" ~ "SSP3-7.0",
      "ssp585" ~ "SSP5-8.5",
      .default = scenario
    ),
    # Formatted as Median (17th–83rd percentile) matching the NASA web tool
    `Rate (2040–2060)` = sprintf("%.1f (%.1f–%.1f)", q50, q17, q83)
  ) %>%
  arrange(scenario) %>%
  select(
    Scenario       = Scenario_Label,
    Median_q50     = q50,
    Lower_17th_q17 = q17,
    Upper_83rd_q83 = q83,
    `Rate (2040–2060)`
  )

# ------------------------------------------------------------------------------
# 3. DISPLAY & EXPORT
# ------------------------------------------------------------------------------

print(slr_rates_2050)

write_csv(slr_rates_2050, out_csv)
cat(sprintf("\n✅ 2040–2060 SLR rates successfully exported to:\n   -> %s\n", out_csv))