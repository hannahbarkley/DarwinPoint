# ==============================================================================
# Project:  Carbonate Budgets / Darwin Point
# Script:   generate_macro_dbase_regional_pooled.R
# Purpose:  Calculate regional pooled bioerosion rates for MHI and NWHI and 
#           generate an NCRMP-formatted database table without external templates.
# Input:    data/supporting data/HA_BMU_2013-2019.csv
# Output:   data/supporting data/macro_dbase_ncrmp_hawaii.csv
# ==============================================================================

library(here)
library(dplyr)
library(tibble)

# 1. Load Data -----------------------------------------------------------------

bmu_path <- here("data", "supporting data", "HA_BMU_2013-2019.csv")
bmu      <- read.csv(bmu_path)

# 2. Dimensions & Rate Calculations --------------------------------------------

l <- 5
w <- 2
h <- 1
surface_area_cm2 <- (2 * l * w) + (2 * l * h) + (2 * w * h) # 34 cm^2

bmu <- bmu %>%
  mutate(
    # g/cm2/yr
    Rate_g_cm2_yr = (PostCTMacroboringRate_cm3_yr * PostCTDensity_g_cm3) / surface_area_cm2,
    # kg/m2/yr
    MacroboringRate_kg_m2_yr = (Rate_g_cm2_yr / 1000) / 0.0001,
    # Assign region
    REGIONCODE = case_when(
      Island %in% c("HAW", "OAH") ~ "MHI",
      Island %in% c("FFS", "PHR") ~ "NWHI"
    )
  )

# Function for 95% confidence intervals using Student's t distribution
calc_ci95 <- function(x) {
  x <- na.omit(x)
  n <- length(x)
  if (n <= 1) return(NA_real_)
  qt(0.975, df = n - 1) * (sd(x) / sqrt(n))
}

# 3. Calculate Regional Pooled Rates -------------------------------------------

regional_summary <- bmu %>%
  group_by(REGIONCODE) %>%
  summarize(
    MACROBIOEROSION_RATE_KG_M2_YR    = round(mean(MacroboringRate_kg_m2_yr, na.rm = TRUE), 6),
    MACROBIOEROSION_RATE_KG_M2_YR_CI = round(calc_ci95(MacroboringRate_kg_m2_yr), 6),
    .groups = "drop"
  )

# 4. Construct Full NCRMP Hawaiian Islands Table -------------------------------

hawaii_islands <- tribble(
  ~REGIONCODE, ~REGION,                         ~LOCATIONCODE, ~LOCATION,
  "MHI",       "Main Hawaiian Islands",         "HAW",         "Hawaii",
  "MHI",       "Main Hawaiian Islands",         "KAH",         "Kahoolawe",
  "MHI",       "Main Hawaiian Islands",         "KAU",         "Kauai",
  "MHI",       "Main Hawaiian Islands",         "LAN",         "Lanai",
  "MHI",       "Main Hawaiian Islands",         "MAI",         "Maui",
  "MHI",       "Main Hawaiian Islands",         "MOL",         "Molokai",
  "MHI",       "Main Hawaiian Islands",         "NII",         "Niihau",
  "MHI",       "Main Hawaiian Islands",         "OAH",         "Oahu",
  "NWHI",      "Northwestern Hawaiian Islands", "FFS",         "French Frigate",
  "NWHI",      "Northwestern Hawaiian Islands", "KUR",         "Kure",
  "NWHI",      "Northwestern Hawaiian Islands", "LIS",         "Lisianski",
  "NWHI",      "Northwestern Hawaiian Islands", "PHR",         "Pearl & Hermes"
)

macro_dbase_ncrmp_ha <- hawaii_islands %>%
  inner_join(regional_summary, by = "REGIONCODE")

# 5. Export CSV ----------------------------------------------------------------

out_path <- here("data", "supporting data", "macro_dbase_ncrmp_hawaii.csv")
write.csv(macro_dbase_ncrmp_ha, out_path, row.names = FALSE)