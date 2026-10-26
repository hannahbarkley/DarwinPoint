# ==============================================================================
# Script:  4b. DarwinPoint_grigg_p90_island.R
# Purpose: Hawaiian Archipelago - Grigg-Equivalent 90th Percentile Island Aggregation
#          Calculates the 90th quantile (P90) and non-parametric bootstrap standard 
#          errors (boot_q90_se) across Grigg-equivalent stations (SW, fore-reef, mid-depth).
#          Outputs a summary table structured identically to spatial_Imean.
# Inputs:  - data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv
# Outputs: - data/grigg data/HA_2024_Grigg_IslandLevel_P90.csv
# ==============================================================================

library(tidyverse)
library(here)
library(crayon)

# ==============================================================================
# 1. PARAMETERS & HELPER FUNCTIONS
# ==============================================================================

region <- "HA"
year   <- "2024"

base_dir  <- here()
data_dir  <- file.path(base_dir, "data")
grigg_dir <- file.path(data_dir, "grigg data")

grigg_csv   <- file.path(grigg_dir, sprintf("%s_%s_Grigg_Equivalent_Sites.csv", region, year))
out_p90_csv <- file.path(grigg_dir, sprintf("%s_%s_Grigg_IslandLevel_P90.csv", region, year))

stopifnot("Grigg-equivalent sites CSV not found! Run 4a first." = file.exists(grigg_csv))

# Non-parametric bootstrap SE for the 90th quantile
boot_q90_se <- function(x, n_boot = 1000) {
  x_clean <- na.omit(x)
  if (length(x_clean) < 2) return(0)
  boot_reps <- replicate(n_boot, {
    s <- sample(x_clean, size = length(x_clean), replace = TRUE)
    quantile(s, probs = 0.90, na.rm = TRUE, names = FALSE)
  })
  sd(boot_reps)
}

# ==============================================================================
# 2. INGESTION & DATA PREPARATION
# ==============================================================================

cat(cyan(sprintf("\n🚀 CALCULATING 90th QUANTILE ISLAND AGGREGATIONS: %s %s\n", region, year)))
cat(blue(sprintf("📂 Loading Grigg sites: %s\n", basename(grigg_csv))))

df_strs_sw <- read_csv(grigg_csv, show_col_types = FALSE)

# ==============================================================================
# 3. 90th PERCENTILE AGGREGATION & BOOTSTRAP SE
# ==============================================================================

cat(blue("🧮 Running 1,000 bootstrap iterations per island for 90th percentiles...\n"))

set.seed(42)

strs_sw_p90 <- df_strs_sw %>%
  group_by(LOCATION, LOCATIONCODE) %>%
  summarize(
    LATITUDE                          = mean(LATITUDE, na.rm = TRUE),
    N_SITES                           = n(),
    
    # 90th Quantile Metric Values
    GROSS_CARB_PROD_KG_M2_YR_MEAN     = quantile(GROSS_CARB_PROD_KG_M2_YR_MEAN, probs = 0.90, na.rm = TRUE, names = FALSE),
    NET_CARB_PROD_KG_M2_YR_MEAN       = quantile(NET_CARB_PROD_KG_M2_YR_MEAN, probs = 0.90, na.rm = TRUE, names = FALSE),
    FISH_EROSION_KG_M2_YR_ALL_MEAN    = quantile(FISH_EROSION_KG_M2_YR_ALL_MEAN, probs = 0.90, na.rm = TRUE, names = FALSE),
    URCHIN_EROSION_KG_M2_YR_MEAN      = quantile(URCHIN_EROSION_KG_M2_YR_MEAN, probs = 0.90, na.rm = TRUE, names = FALSE),
    GROSS_EROSION_KG_M2_YR_MEAN       = quantile(GROSS_EROSION_KG_M2_YR_MEAN, probs = 0.90, na.rm = TRUE, names = FALSE),
    HARD_CORAL_COVER_PCT_MEAN         = quantile(HARD_CORAL_COVER_PCT_MEAN, probs = 0.90, na.rm = TRUE, names = FALSE),
    
    # Bootstrap Standard Errors (Matching spatial_Imean 'se.' prefix)
    se.GROSS_CARB_PROD_KG_M2_YR_MEAN  = boot_q90_se(GROSS_CARB_PROD_KG_M2_YR_MEAN),
    se.NET_CARB_PROD_KG_M2_YR_MEAN    = boot_q90_se(NET_CARB_PROD_KG_M2_YR_MEAN),
    se.FISH_EROSION_KG_M2_YR_ALL_MEAN = boot_q90_se(FISH_EROSION_KG_M2_YR_ALL_MEAN),
    se.URCHIN_EROSION_KG_M2_YR_MEAN   = boot_q90_se(URCHIN_EROSION_KG_M2_YR_MEAN),
    se.GROSS_EROSION_KG_M2_YR_MEAN    = boot_q90_se(GROSS_EROSION_KG_M2_YR_MEAN),
    se.HARD_CORAL_COVER_PCT_MEAN      = boot_q90_se(HARD_CORAL_COVER_PCT_MEAN),
    
    .groups = "drop"
  ) %>%
  arrange(LATITUDE)

# ==============================================================================
# 4. EXPORT
# ==============================================================================

write_csv(strs_sw_p90, out_p90_csv)
cat(green(sprintf("✅ Saved spatial_Imean-style P90 dataset (%d locations) to:\n   -> %s\n", 
                  nrow(strs_sw_p90), out_p90_csv)))