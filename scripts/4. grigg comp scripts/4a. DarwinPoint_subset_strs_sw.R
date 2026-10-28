# ==============================================================================
# Script:  4a. DarwinPoint_subset_strs_sw.R
# Purpose: Hawaiian Archipelago - Combined StRS & Fixed Grigg-Equivalent Analysis
#          Ingests processed StRS and Fixed net carbonate budget datasets,
#          merges REEF_ZONE, DEPTH_BIN, and HABITAT_CODE from prod_sitelevel into
#          strs_raw, pulls site-specific framework porosity (Est_Porosity) from
#          HA_2024_StRS_Porosity.csv (rather than assuming a flat 0.383),
#          calculates RAPmax (maximum reef accretion potential from gross production),
#          tags survey type, de-duplicates stations, and filters to match 
#          Grigg's (1982) environmental criteria:
#            1) Southwest geographic quadrant (relative to island centroid)
#            2) Exposed Fore-reef habitat (REEF_ZONE containing "Fore" or "Slope")
#            3) Standardized Mid depth bin (DEPTH_BIN == "Mid", 6–18 m)
# Inputs:  - data/strs data/processed data/HA_2024_StRS_Net_SiteLevel_cnet.csv
#          - data/strs data/field data/HA_2024_StRS_Prod_SiteLevel.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_cnet.csv
#          - data/HA_2024_StRS_Porosity.csv (or data/strs data/processed data/)
# Outputs: - data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv
# ==============================================================================

library(tidyverse)
library(here)
library(crayon)

# ==============================================================================
# 1. PARAMETERS & DIRECTORIES
# ==============================================================================

region <- "HA"
year   <- "2024"

# Framework Accretion Parameters 
rho_aragonite <- 2.93  # Aragonite mineral density (g/cm^3 or kg/m^2/mm)
phi_default   <- 0.383 # Fallback porosity for stations without coral / missing entries

base_dir       <- here()
data_dir       <- file.path(base_dir, "data")
strs_proc_dir  <- file.path(data_dir, "strs data", "processed data")
strs_field_dir <- file.path(data_dir, "strs data", "field data")
fixed_proc_dir <- file.path(data_dir, "fixed data", "processed data")
out_csv_dir    <- file.path(data_dir, "grigg data")

strs_file <- file.path(strs_proc_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year))

prod_candidates <- c(
  file.path(strs_field_dir, sprintf("%s_%s_StRS_Prod_SiteLevel.csv", region, year)),
  file.path(data_dir, sprintf("%s_%s_StRS_Prod_SiteLevel.csv", region, year)),
  here(sprintf("%s_%s_StRS_Prod_SiteLevel.csv", region, year)),
  here("HA_2024_StRS_Prod_SiteLevel.csv")
)
prod_file <- prod_candidates[file.exists(prod_candidates)][1]

fixed_candidates <- c(
  file.path(fixed_proc_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_cnet.csv", region, year))
)
fixed_file <- fixed_candidates[file.exists(fixed_candidates)][1]

# Porosity file candidates
porosity_candidates <- c(
  file.path(strs_proc_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year)),
  file.path(data_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year)),
  here("data", sprintf("%s_%s_StRS_Porosity.csv", region, year)),
  here(sprintf("%s_%s_StRS_Porosity.csv", region, year)),
  here("HA_2024_StRS_Porosity.csv")
)
porosity_file <- porosity_candidates[file.exists(porosity_candidates)][1]

stopifnot("StRS dataset not found!"      = !is.na(strs_file) && file.exists(strs_file))
stopifnot("StRS Prod file not found!"    = !is.na(prod_file) && file.exists(prod_file))
stopifnot("Fixed dataset not found!"     = !is.na(fixed_file) && file.exists(fixed_file))
stopifnot("Porosity dataset not found!"  = !is.na(porosity_file) && file.exists(porosity_file))

out_subset_csv <- file.path(out_csv_dir, sprintf("%s_%s_Grigg_Equivalent_Sites.csv", region, year))
if (!dir.exists(out_csv_dir)) dir.create(out_csv_dir, recursive = TRUE)

# ==============================================================================
# 2. LOAD & COMBINE DATASETS
# ==============================================================================

cat(cyan(sprintf("\n🚀 INGESTING StRS AND FIXED CARBONATE BUDGET DATA: %s %s\n", region, year)))
cat(blue(sprintf("📂 Loading StRS:      %s\n", basename(strs_file))))
cat(blue(sprintf("📂 Loading StRS Prod: %s\n", basename(prod_file))))
cat(blue(sprintf("📂 Loading Fixed:     %s\n", basename(fixed_file))))
cat(blue(sprintf("📂 Loading Porosity:  %s\n", basename(porosity_file))))

# Ingest and isolate REEF_ZONE, DEPTH_BIN, HABITAT_CODE from prod_sitelevel
prod_meta <- read_csv(prod_file, show_col_types = FALSE) %>%
  select(SITE, REEF_ZONE, DEPTH_BIN, HABITAT_CODE) %>%
  distinct(SITE, .keep_all = TRUE)

# Ingest site-specific porosity estimates derived from morphology covers
df_porosity <- read_csv(porosity_file, show_col_types = FALSE) %>%
  select(SITE, Est_Porosity) %>%
  distinct(SITE, .keep_all = TRUE)

# Join metadata into strs_raw before combined_all
strs_raw <- read_csv(strs_file, show_col_types = FALSE) %>%
  select(-any_of(c("REEF_ZONE", "DEPTH_BIN", "HABITAT_CODE"))) %>%
  left_join(prod_meta, by = "SITE") %>%
  mutate(SURVEY_TYPE = "StRS")

fixed_raw <- read_csv(fixed_file, show_col_types = FALSE) %>% 
  rename(SITE = OCC_SITEID) %>% 
  mutate(SURVEY_TYPE = "Fixed")

# Fixed sites take precedence over any StRS duplicates
combined_all <- bind_rows(
  fixed_raw,
  strs_raw %>% filter(!SITE %in% fixed_raw$SITE)
) %>%
  mutate(
    REEF_ZONE = coalesce(REEF_ZONE, "Forereef"),
    DEPTH_BIN = coalesce(
      DEPTH_BIN,
      case_when(
        as.numeric(SITE_DEPTH_M) < 6   ~ "Shallow",
        as.numeric(SITE_DEPTH_M) <= 18 ~ "Mid",
        as.numeric(SITE_DEPTH_M) > 18  ~ "Deep",
        TRUE                           ~ "Mid" # Preserves unmeasured depth stations
      )
    )
  )

cat(green(sprintf("✅ Combined master dataset: %d unique sites (%d Fixed, %d StRS)\n", 
                  nrow(combined_all), nrow(fixed_raw), nrow(combined_all) - nrow(fixed_raw))))

# ==============================================================================
# 3. MERGE SITE-SPECIFIC POROSITY & CALCULATE RAPmax
# ==============================================================================

cat(blue("📐 Joining site-level porosity and calculating RAPmax (rho = 2.93)...\n"))

combined_all <- combined_all %>%
  select(-any_of("Est_Porosity")) %>%
  left_join(df_porosity, by = "SITE") %>%
  mutate(
    GROSS_CARB_PROD_KG_M2_YR_MEAN = coalesce(!!!select(., any_of(c("GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_PROD_KG_M2_YR", "GROSS_PROD")))),
    # Fallback to default porosity (0.383) only if site is absent from porosity CSV
    Est_Porosity                  = if_else(is.na(Est_Porosity), phi_default, as.numeric(Est_Porosity)),
    RAPmax                        = GROSS_CARB_PROD_KG_M2_YR_MEAN / (rho_aragonite * (1 - Est_Porosity))
  )

# ==============================================================================
# 4. APPLY GRIGG (1982) CRITERIA & EXPORT
# ==============================================================================

cat(blue("🔍 Applying Grigg (1982) criteria (SW Quadrant, Fore-Reef, Mid Depth)...\n"))

island_centers <- tibble(
  LOCATION     = c("French Frigate", "Hawaii", "Kahoolawe", "Kauai", "Kure", 
                   "Lanai", "Lisianski", "Maui", "Molokai", "Niihau", 
                   "Oahu", "Pearl & Hermes"),
  ISL_CODE     = c("FFS", "HAW", "KAH", "KAU", "KUR", 
                   "LAN", "LIS", "MAI", "MOL", "NII", 
                   "OAH", "PHR"),
  MID_LAT      = c(23.810603, 19.6059, 20.5510, 22.0673, 28.421506, 
                   20.8251, 26.008619, 20.8037, 21.1396, 21.8957, 
                   21.4883, 27.840445),
  MID_LON      = c(-166.157355, -155.5113, -156.6322, -159.5249, -178.328858, 
                   -156.9248, -173.949556, -156.3311, -157.0201, -160.1554, 
                   -157.9942, -175.848049)
)

grigg_subset <- combined_all %>%
  left_join(island_centers, by = "LOCATION") %>%
  mutate(LOCATIONCODE = coalesce(as.character(LOCATIONCODE), ISL_CODE)) %>%
  filter(
    # Criteria 1: Southwest quadrant of island centroid
    LATITUDE < MID_LAT & LONGITUDE < MID_LON,
    
    # Criteria 2: forereef
    REEF_ZONE == "Forereef",
    
    # Criteria 3: Mid depth bin (6–18 m)
    DEPTH_BIN == "Mid"
  ) %>%
  select(-MID_LAT, -MID_LON, -ISL_CODE)

cat(green(sprintf("✅ Filtered to %d Grigg-equivalent stations (%d Fixed, %d StRS)\n", 
                  nrow(grigg_subset), 
                  sum(grigg_subset$SURVEY_TYPE == "Fixed"), 
                  sum(grigg_subset$SURVEY_TYPE == "StRS"))))

write_csv(grigg_subset, out_subset_csv)
cat(green(sprintf("💾 Saved Grigg-equivalent dataset to: %s\n", out_subset_csv)))