# ==============================================================================
# Script:  4. DarwinPoint_grigg_comp.R
# Purpose: Master Wrapper for Hawaiian Archipelago Grigg-Equivalent P90 Pipeline
#          1. Subsets combined StRS & Fixed stations to Grigg-equivalent criteria:
#             - SW geographic quadrant relative to island centroid
#             - Exposed Fore-reef / Slope habitat
#             - Mid-depth bin (6–18 m)
#          2. Computes the 90th percentile (P90) and non-parametric bootstrap SEs
#             (B = 1,000 iterations) across all primary carbonate budget metrics.
#          3. Exports both the filtered site-level dataset and the island summary.
# Inputs:  - data/strs data/processed data/HA_2024_StRS_Net_SiteLevel_cnet.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_cnet.csv
# Outputs: - data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv
#          - data/grigg data/HA_2024_Grigg_IslandLevel_P90.csv
# ==============================================================================

library(tidyverse)
library(here)
library(crayon)

# ==============================================================================
# 1. CORE STATISTICAL & SPATIAL HELPERS
# ==============================================================================

# Non-parametric bootstrap SE for arbitrary quantiles
boot_quantile_se <- function(x, prob = 0.90, n_boot = 1000) {
  x_clean <- na.omit(x)
  if (length(x_clean) < 2) return(0)
  boot_reps <- replicate(n_boot, {
    s <- sample(x_clean, size = length(x_clean), replace = TRUE)
    quantile(s, probs = prob, na.rm = TRUE, names = FALSE)
  })
  sd(boot_reps)
}

# Master island centroids and 3-letter codes for SW quadrant filtering
get_island_reference <- function() {
  tibble(
    LOCATION = c("French Frigate", "Hawaii", "Kahoolawe", "Kauai", "Kure", 
                 "Lanai", "Lisianski", "Maui", "Molokai", "Niihau", 
                 "Oahu", "Pearl & Hermes"),
    ISL_CODE = c("FFS", "HAW", "KAH", "KAU", "KUR", 
                 "LAN", "LIS", "MAI", "MOL", "NII", 
                 "OAH", "PHR"),
    MID_LAT  = c(23.810603, 19.6059, 20.5510, 22.0673, 28.421506, 
                 20.8251, 26.008619, 20.8037, 21.1396, 21.8957, 
                 21.4883, 27.855445),
    MID_LON  = c(-166.157355, -155.5113, -156.6322, -159.5249, -178.328858, 
                 -156.9248, -173.949556, -156.3311, -157.0201, -160.1554, 
                 -157.9942, -175.848049)
  )
}

# ==============================================================================
# 2. PIPELINE STEP FUNCTIONS
# ==============================================================================

# Step 4a: Ingest, combine, and subset to Grigg criteria
subset_grigg_equivalent <- function(strs_file, fixed_file, out_site_csv = NULL) {
  cat(blue("🔍 Step 1: Ingesting datasets and applying Grigg (1982) criteria...\n"))
  
  strs_raw  <- read_csv(strs_file, show_col_types = FALSE) %>% mutate(SURVEY_TYPE = "StRS")
  fixed_raw <- read_csv(fixed_file, show_col_types = FALSE) %>% 
    rename(SITE = OCC_SITEID) %>% 
    mutate(SURVEY_TYPE = "Fixed")
  
  # De-duplicate: Fixed stations take precedence over StRS
  combined_all <- bind_rows(
    fixed_raw,
    strs_raw %>% filter(!SITE %in% fixed_raw$SITE)
  ) %>%
    mutate(
      REEF_ZONE = if ("REEF_ZONE" %in% names(.)) REEF_ZONE else "Fore-reef",
      DEPTH_BIN = case_when(
        as.numeric(SITE_DEPTH_M) < 6   ~ "Shallow",
        as.numeric(SITE_DEPTH_M) <= 18 ~ "Mid",
        as.numeric(SITE_DEPTH_M) > 18  ~ "Deep",
        TRUE                           ~ "Mid" # Preserves unmeasured depth stations
      )
    )
  
  island_ref <- get_island_reference()
  
  grigg_subset <- combined_all %>%
    left_join(island_ref, by = "LOCATION") %>%
    mutate(LOCATIONCODE = coalesce(as.character(LOCATIONCODE), ISL_CODE)) %>%
    filter(
      # Criteria 1: SW quadrant relative to island centroid
      LATITUDE < MID_LAT & LONGITUDE < MID_LON,
      # Criteria 2: Exposed fore-reef or slope
      grepl("FORE|SLOPE", toupper(REEF_ZONE)),
      # Criteria 3: Mid-depth bin (6–18 m)
      DEPTH_BIN == "Mid"
    ) %>%
    select(-MID_LAT, -MID_LON, -ISL_CODE)
  
  cat(green(sprintf("  • Retained %d Grigg-equivalent stations (%d Fixed, %d StRS)\n", 
                    nrow(grigg_subset), 
                    sum(grigg_subset$SURVEY_TYPE == "Fixed"), 
                    sum(grigg_subset$SURVEY_TYPE == "StRS"))))
  
  if (!is.null(out_site_csv)) {
    write_csv(grigg_subset, out_site_csv)
    cat(green(sprintf("  • Saved site subset to: %s\n", out_site_csv)))
  }
  
  return(grigg_subset)
}

# Step 4b: Calculate 90th percentile metrics with bootstrap SE
aggregate_grigg_p90 <- function(grigg_df, out_island_csv = NULL, prob = 0.90, n_boot = 1000, seed = 42) {
  cat(blue(sprintf("🧮 Step 2: Bootstrapping P%d metrics (B = %d iterations per island)...\n", 
                   round(prob * 100), n_boot)))
  
  set.seed(seed)
  
  island_p90 <- grigg_df %>%
    group_by(LOCATION, LOCATIONCODE) %>%
    summarize(
      LATITUDE                          = mean(LATITUDE, na.rm = TRUE),
      N_SITES                           = n(),
      
      # Quantile estimates
      GROSS_CARB_PROD_KG_M2_YR_MEAN     = quantile(GROSS_CARB_PROD_KG_M2_YR_MEAN, probs = prob, na.rm = TRUE, names = FALSE),
      NET_CARB_PROD_KG_M2_YR_MEAN       = quantile(NET_CARB_PROD_KG_M2_YR_MEAN, probs = prob, na.rm = TRUE, names = FALSE),
      FISH_EROSION_KG_M2_YR_ALL_MEAN    = quantile(FISH_EROSION_KG_M2_YR_ALL_MEAN, probs = prob, na.rm = TRUE, names = FALSE),
      URCHIN_EROSION_KG_M2_YR_MEAN      = quantile(URCHIN_EROSION_KG_M2_YR_MEAN, probs = prob, na.rm = TRUE, names = FALSE),
      GROSS_EROSION_KG_M2_YR_MEAN       = quantile(GROSS_EROSION_KG_M2_YR_MEAN, probs = prob, na.rm = TRUE, names = FALSE),
      HARD_CORAL_COVER_PCT_MEAN         = quantile(HARD_CORAL_COVER_PCT_MEAN, probs = prob, na.rm = TRUE, names = FALSE),
      
      # Bootstrap Standard Errors (Matching spatial_Imean 'se.' prefix)
      se.GROSS_CARB_PROD_KG_M2_YR_MEAN  = boot_quantile_se(GROSS_CARB_PROD_KG_M2_YR_MEAN, prob = prob, n_boot = n_boot),
      se.NET_CARB_PROD_KG_M2_YR_MEAN    = boot_quantile_se(NET_CARB_PROD_KG_M2_YR_MEAN, prob = prob, n_boot = n_boot),
      se.FISH_EROSION_KG_M2_YR_ALL_MEAN = boot_quantile_se(FISH_EROSION_KG_M2_YR_ALL_MEAN, prob = prob, n_boot = n_boot),
      se.URCHIN_EROSION_KG_M2_YR_MEAN   = boot_quantile_se(URCHIN_EROSION_KG_M2_YR_MEAN, prob = prob, n_boot = n_boot),
      se.GROSS_EROSION_KG_M2_YR_MEAN    = boot_quantile_se(GROSS_EROSION_KG_M2_YR_MEAN, prob = prob, n_boot = n_boot),
      se.HARD_CORAL_COVER_PCT_MEAN      = boot_quantile_se(HARD_CORAL_COVER_PCT_MEAN, prob = prob, n_boot = n_boot),
      
      .groups = "drop"
    ) %>%
    arrange(LATITUDE)
  
  if (!is.null(out_island_csv)) {
    write_csv(island_p90, out_island_csv)
    cat(green(sprintf("  • Saved island P%d table to: %s\n", round(prob * 100), out_island_csv)))
  }
  
  return(island_p90)
}

# ==============================================================================
# 3. MASTER WRAPPER FUNCTION
# ==============================================================================

run_grigg_p90_pipeline <- function(region         = "HA", 
                                   year           = "2024", 
                                   base_dir       = here(),
                                   prob           = 0.90, 
                                   n_boot         = 1000, 
                                   seed           = 42) {
  
  t0 <- Sys.time()
  cat(cyan(sprintf("\n=======================================================\n")))
  cat(cyan(sprintf("🚀 RUNNING GRIGG P%d PIPELINE: %s %s\n", round(prob * 100), region, year)))
  cat(cyan(sprintf("=======================================================\n")))
  
  data_dir       <- file.path(base_dir, "data")
  strs_proc_dir  <- file.path(data_dir, "strs data", "processed data")
  fixed_proc_dir <- file.path(data_dir, "fixed data", "processed data")
  out_csv_dir    <- file.path(data_dir, "grigg data")
  
  if (!dir.exists(out_csv_dir)) dir.create(out_csv_dir, recursive = TRUE)
  
  # Resolve input file paths
  strs_file <- file.path(strs_proc_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year))
  fixed_candidates <- c(
    file.path(fixed_proc_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_cnet_2.csv", region, year)),
    file.path(fixed_proc_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_cnet.csv", region, year))
  )
  fixed_file <- fixed_candidates[file.exists(fixed_candidates)][1]
  
  stopifnot("StRS dataset not found!"  = !is.na(strs_file) && file.exists(strs_file))
  stopifnot("Fixed dataset not found!" = !is.na(fixed_file) && file.exists(fixed_file))
  
  # Output file destinations
  out_site_csv   <- file.path(out_csv_dir, sprintf("%s_%s_Grigg_Equivalent_Sites.csv", region, year))
  out_island_csv <- file.path(out_csv_dir, sprintf("%s_%s_Grigg_IslandLevel_P%d.csv", region, year, round(prob * 100)))
  
  # 1. Execute Subsetting (4a)
  grigg_sites <- subset_grigg_equivalent(
    strs_file    = strs_file, 
    fixed_file   = fixed_file, 
    out_site_csv = out_site_csv
  )
  
  # 2. Execute Bootstrap Aggregation (4b)
  island_summary <- aggregate_grigg_p90(
    grigg_df       = grigg_sites, 
    out_island_csv = out_island_csv, 
    prob           = prob, 
    n_boot         = n_boot, 
    seed           = seed
  )
  
  elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
  cat(green(sprintf("\n✅ Pipeline complete in %s seconds! Processed %d sites across %d islands.\n", 
                    elapsed, nrow(grigg_sites), nrow(island_summary))))
  
  invisible(list(sites = grigg_sites, island_p90 = island_summary))
}

# ==============================================================================
# 4. DIRECT EXECUTION
# ==============================================================================

pipeline_output <- run_grigg_p90_pipeline(
  region = region, 
  year   = year, 
  prob   = 0.90, 
  n_boot = 1000
)