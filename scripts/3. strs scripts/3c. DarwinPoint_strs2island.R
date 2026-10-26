# ==============================================================================
# Script:  3c. DarwinPoint_strs2island.R
# Purpose: Hawaiian Archipelago - StRS Carbonate Budget & Island-Level Scaling
#          Merges StRS carbonate budgets with sector spatial metadata,
#          incorporates site-level framework porosity, calculates stratum 
#          sampling weights (NH / n), constructs a survey design object 
#          (svydesign), computes survey-weighted island-level metric means, 
#          and exports the survey design object for downstream modeling (3e).
# Inputs:  - data/strs data/processed data/HA_2024_StRS_Net_SiteLevel_cnet.csv
#            (or HA_2024_StRS_Net_SiteLevel_cnet_2.csv)
#          - data/strs data/processed data/HA_2024_StRS_Porosity.csv
#          - data/supporting data/SURVEY MASTER.csv (or data/)
#          - data/supporting data/Sectors-Strata-Areas.csv (or data/)
# Outputs: - data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          - data/strs data/processed data/HA_2024_StRS_SiteLevel_Weighted.csv
#          - data/strs data/processed data/HA_2024_StRS_SurveyDesign.rds
# ==============================================================================

library(tidyverse)
library(here)
library(survey)
library(crayon)

# ==============================================================================
# 1. PARAMETERS & CONSTANTS
# ==============================================================================

region <- "HA"
year   <- "2024"

# Density of solid aragonite (g cm⁻³, which equates to 2,890 kg m⁻³)
# Dividing (kg m⁻² yr⁻¹) by (2.89 * (1 - phi)) directly yields mm yr⁻¹
rho_aragonite <- 2.89  

# Baseline fallback framework porosity for unmeasured sites (Perry et al. 2025, Nature)
phi_default <- 0.383  

base_dir <- here()
data_dir <- file.path(base_dir, "data")
proc_dir <- file.path(data_dir, "strs data", "processed data")

# Pad 4-digit StRS fish/benthic SITEs to 5 digits (e.g., FFS-4723 -> FFS-04723)
pad_strs_site <- function(site_vec) {
  s <- trimws(as.character(site_vec))
  is_2part <- grepl("^[A-Za-z]+-[0-9]+$", s)
  prefix   <- sub("-.*", "", s)
  num_val  <- suppressWarnings(as.numeric(sub("^[A-Za-z]+-", "", s)))
  
  out <- s
  idx <- which(is_2part & !is.na(num_val))
  out[idx] <- sprintf("%s-%05d", prefix[idx], num_val[idx])
  return(out)
}

# Standardized NOAA PIFSC / NCRMP 3-letter location code lookup
get_location_code <- function(loc) {
  loc_clean <- trimws(as.character(loc))
  u_loc     <- toupper(loc_clean)
  
  case_when(
    is.na(loc) | loc_clean == ""       ~ NA_character_,
    grepl("FRENCH|FFS", u_loc)         ~ "FFS",
    grepl("PEARL|PHR",  u_loc)         ~ "PHR",
    grepl("LIS",        u_loc)         ~ "LIS",
    grepl("KUR",        u_loc)         ~ "KUR",
    grepl("HAW",        u_loc)         ~ "HAW",
    grepl("MAU|MAI",    u_loc)         ~ "MAI",
    grepl("MOL",        u_loc)         ~ "MOL",
    grepl("LAN",        u_loc)         ~ "LAN",
    grepl("KAH",        u_loc)         ~ "KAH",
    grepl("OAH",        u_loc)         ~ "OAH",
    grepl("KAU",        u_loc)         ~ "KAU",
    grepl("NII",        u_loc)         ~ "NII",
    grepl("MID",        u_loc)         ~ "MID",
    grepl("MAR",        u_loc)         ~ "MAR",
    grepl("LAY",        u_loc)         ~ "LAY",
    grepl("GAR",        u_loc)         ~ "GAR",
    grepl("NEC",        u_loc)         ~ "NEC",
    grepl("NIH",        u_loc)         ~ "NIH",
    TRUE                               ~ loc_clean
  )
}

# Helper function to find files across structured subfolders or flat data/
find_input_file <- function(candidates) {
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  return(candidates[1])
}

# Helper to coalesce candidate column names into numeric values
coalesce_num <- function(df, cols, default = 0) {
  present <- intersect(cols, names(df))
  if (length(present) == 0) return(rep(default, nrow(df)))
  
  vals <- df %>% 
    select(all_of(present)) %>% 
    mutate(across(everything(), as.numeric))
  
  res <- do.call(coalesce, vals)
  coalesce(res, default)
}

strs_csv <- find_input_file(c(
  file.path(proc_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet_2.csv", region, year)),
  file.path(proc_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year)),
  file.path(data_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet_2.csv", region, year)),
  file.path(data_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year)),
  here("HA_2024_StRS_Net_SiteLevel_cnet_2.csv"),
  here("HA_2024_StRS_Net_SiteLevel_cnet.csv")
))

porosity_csv <- find_input_file(c(
  file.path(proc_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year)),
  file.path(data_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year)),
  file.path(data_dir, "HA_2024_StRS_Porosity.csv"),
  here("HA_2024_StRS_Porosity.csv")
))

sm_csv <- find_input_file(c(
  file.path(data_dir, "supporting data", "SURVEY MASTER.csv"),
  file.path(data_dir, "SURVEY MASTER.csv"),
  here("SURVEY MASTER.csv")
))

sectors_csv <- find_input_file(c(
  file.path(data_dir, "supporting data", "Sectors-Strata-Areas.csv"),
  file.path(data_dir, "Sectors-Strata-Areas.csv"),
  here("Sectors-Strata-Areas.csv")
))

out_island_csv        <- file.path(proc_dir, sprintf("%s_%s_StRS_IslandLevel.csv", region, year))
out_site_weighted_csv <- file.path(proc_dir, sprintf("%s_%s_StRS_SiteLevel_Weighted.csv", region, year))
out_des_rds           <- file.path(proc_dir, sprintf("%s_%s_StRS_SurveyDesign.rds", region, year))

stopifnot("StRS Net Site-level dataset not found!" = file.exists(strs_csv))
stopifnot("Survey Master dataset not found!"       = file.exists(sm_csv))
stopifnot("Sectors/Strata area file not found!"     = file.exists(sectors_csv))
stopifnot("Porosity dataset not found!"            = file.exists(porosity_csv))

if (!dir.exists(proc_dir)) dir.create(proc_dir, recursive = TRUE)

# ==============================================================================
# 2. DATA INGESTION & HARMONIZATION
# ==============================================================================

cat(cyan(sprintf("\n🚀 STARTING StRS-TO-ISLAND SURVEY ANALYSIS: %s %s\n", region, year)))
cat(blue(sprintf("📂 Loading StRS:          %s\n", basename(strs_csv))))
cat(blue(sprintf("📂 Loading Porosity:      %s\n", basename(porosity_csv))))
cat(blue(sprintf("📂 Loading Survey Master: %s\n", basename(sm_csv))))
cat(blue(sprintf("📂 Loading Sectors Areas: %s\n", basename(sectors_csv))))

strs     <- read_csv(strs_csv, show_col_types = FALSE)
sm       <- read_csv(sm_csv, show_col_types = FALSE)
sectors  <- read_csv(sectors_csv, show_col_types = FALSE)
porosity <- read_csv(porosity_csv, show_col_types = FALSE) %>%
  mutate(SITE = pad_strs_site(SITE)) %>%
  select(SITE, Est_Porosity) %>%
  distinct(SITE, .keep_all = TRUE)

cat(blue("🔧 Standardizing site keys, joining spatial strata, and computing accretion...\n"))

# Pull spatial stratum attributes from Survey Master
sm_meta <- sm %>%
  mutate(SITE = pad_strs_site(SITE)) %>%
  select(any_of(c("SITE", "SEC_NAME", "REEF_ZONE", "DEPTH_BIN", "LATITUDE_SV", "OBS_YEAR", "ANALYSIS_SEC"))) %>%
  distinct(SITE, .keep_all = TRUE)

# Harmonize StRS metrics using the NCEI process_net output schema
strs_clean <- strs %>%
  rename(any_of(c(SITE = "OCC_SITEID"))) %>%
  mutate(
    SITE         = pad_strs_site(SITE),
    LOCATION     = if ("LOCATION" %in% names(.)) {
      if ("ISLAND" %in% names(.)) coalesce(LOCATION, ISLAND) else LOCATION
    } else if ("ISLAND" %in% names(.)) {
      ISLAND
    } else {
      NA_character_
    },
    LOCATIONCODE = coalesce(
      if ("LOCATIONCODE" %in% names(.)) as.character(LOCATIONCODE) else NULL,
      get_location_code(LOCATION)
    ),
    REGIONCODE   = if ("REGIONCODE" %in% names(.)) REGIONCODE else if ("REGION" %in% names(.)) REGION else "MHI",
    
    # Primary Budget Metrics
    GROSS_CARB_PROD_KG_M2_YR_MEAN  = coalesce_num(., c("GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_PROD_KG_M2_YR", "PROD")),
    NET_CARB_PROD_KG_M2_YR_MEAN    = coalesce_num(., c("NET_CARB_PROD_KG_M2_YR_MEAN", "NET_PROD_KG_M2_YR")),
    FISH_EROSION_KG_M2_YR_ALL_MEAN = coalesce_num(., c("FISH_EROSION_KG_M2_YR_ALL_MEAN", "FISH_EROSION_KG_M2_YR_ALL")),
    URCHIN_EROSION_KG_M2_YR_MEAN   = coalesce_num(., "URCHIN_EROSION_KG_M2_YR_MEAN"),
    GROSS_EROSION_KG_M2_YR_MEAN    = coalesce(
      if ("GROSS_EROSION_KG_M2_YR_MEAN" %in% names(.)) as.numeric(GROSS_EROSION_KG_M2_YR_MEAN) else NULL,
      FISH_EROSION_KG_M2_YR_ALL_MEAN + URCHIN_EROSION_KG_M2_YR_MEAN
    ),
    HARD_CORAL_COVER_PCT_MEAN      = coalesce_num(., c("HARD_CORAL_COVER_PCT_MEAN", "HARD_CORAL_COVER_PCT")),
    
    # Backward-compatible metric aliases for downstream modeling scripts
    GROSS_PROD_KG_M2_YR       = GROSS_CARB_PROD_KG_M2_YR_MEAN,
    NET_PROD_KG_M2_YR         = NET_CARB_PROD_KG_M2_YR_MEAN,
    FISH_EROSION_KG_M2_YR_ALL = FISH_EROSION_KG_M2_YR_ALL_MEAN,
    GROSS_EROSION_KG_M2_YR    = GROSS_EROSION_KG_M2_YR_MEAN,
    HARD_CORAL_COVER_PCT      = HARD_CORAL_COVER_PCT_MEAN
  ) %>%
  filter(!grepl("^OCC", trimws(SITE))) %>%
  left_join(porosity, by = "SITE") %>%
  mutate(
    Est_Porosity             = if_else(is.na(Est_Porosity), phi_default, Est_Porosity),
    VERTICAL_ACCRETION_MM_YR = NET_CARB_PROD_KG_M2_YR_MEAN / (rho_aragonite * (1 - Est_Porosity))
  ) %>%
  left_join(sm_meta, by = "SITE", suffix = c("", ".sm")) %>%
  mutate(
    OBS_YEAR    = coalesce(
      if ("OBS_YEAR" %in% names(.)) as.character(OBS_YEAR) else NULL,
      if ("OBS_YEAR.sm" %in% names(.)) as.character(OBS_YEAR.sm) else NULL,
      if ("YEAR" %in% names(.)) as.character(YEAR) else NULL,
      as.character(year)
    ),
    LATITUDE_SV = coalesce(
      if ("LATITUDE_SV" %in% names(.)) LATITUDE_SV else NULL,
      if ("LATITUDE" %in% names(.)) LATITUDE else NULL
    ),
    REEF_ZONE   = coalesce(
      if ("REEF_ZONE" %in% names(.)) REEF_ZONE else NULL,
      if ("REEF_ZONE.sm" %in% names(.)) REEF_ZONE.sm else NULL,
      "Fore-reef"
    ),
    DEPTH_BIN   = coalesce(
      if ("DEPTH_BIN" %in% names(.)) DEPTH_BIN else NULL,
      if ("DEPTH_BIN.sm" %in% names(.)) DEPTH_BIN.sm else NULL,
      case_when(
        as.numeric(SITE_DEPTH_M) < 6   ~ "Shallow",
        as.numeric(SITE_DEPTH_M) <= 18 ~ "Mid",
        as.numeric(SITE_DEPTH_M) > 18  ~ "Deep",
        TRUE                           ~ "Mid"
      )
    )
  ) %>%
  select(-any_of(c("OBS_YEAR.sm", "REEF_ZONE.sm", "DEPTH_BIN.sm")))

# Merge sector spatial strata areas (NH)
meta_strata <- sectors %>%
  select(any_of(c("SEC_NAME", "REEF_ZONE", "DEPTH_BIN", "NH"))) %>%
  distinct()

new_strs <- strs_clean %>%
  left_join(meta_strata, by = c("SEC_NAME", "REEF_ZONE", "DEPTH_BIN")) %>%
  mutate(
    STRATANAME  = paste(SEC_NAME, REEF_ZONE, DEPTH_BIN, sep = "_"),
    REGION_YEAR = paste(REGIONCODE, OBS_YEAR, sep = "_"),
    DB_RZ       = paste(REEF_ZONE, DEPTH_BIN, sep = "_")
  )

# ==============================================================================
# 3. STRATUM WEIGHTING & SURVEY DESIGN
# ==============================================================================

cat(blue("⚖️ Calculating stratum sample weights (NH / n) and survey design...\n"))

site_sw <- new_strs %>%
  mutate(
    n  = n_distinct(SITE), 
    sw = NH / n,
    .by = c(OBS_YEAR, SEC_NAME, REEF_ZONE, DEPTH_BIN, NH)
  ) %>%
  mutate(LOCATION_LAT = mean(LATITUDE_SV, na.rm = TRUE), .by = LOCATION) %>%
  filter(!is.na(sw) & n > 1) %>% # Retain strata with >1 site to enable design-based variance
  mutate(
    across(any_of(c("OBS_YEAR", "LOCATION", "LOCATIONCODE", "REGION", "REGIONCODE")), as.factor),
    Strat_conc = paste(OBS_YEAR, REGIONCODE, LOCATION, SEC_NAME, DB_RZ, sep = "_")
  )

des <- svydesign(id = ~1, strata = ~Strat_conc, weights = ~sw, data = site_sw)

# ==============================================================================
# 4. EXPORT WEIGHTED SITE-LEVEL DATA & SURVEY DESIGN (FOR 3e)
# ==============================================================================

write_csv(site_sw, out_site_weighted_csv)
saveRDS(des, out_des_rds)

cat(green(sprintf("✅ Saved site-level weighted data to: %s\n", out_site_weighted_csv)))
cat(green(sprintf("✅ Saved survey design object to:    %s\n", out_des_rds)))

# ==============================================================================
# 5. WEIGHTED ISLAND/LOCATION-LEVEL AGGREGATION & EXPORT
# ==============================================================================

cat(blue("📊 Computing survey-weighted location-level means...\n"))

metrics <- c(
  "GROSS_CARB_PROD_KG_M2_YR_MEAN", 
  "NET_CARB_PROD_KG_M2_YR_MEAN", 
  "FISH_EROSION_KG_M2_YR_ALL_MEAN", 
  "URCHIN_EROSION_KG_M2_YR_MEAN",
  "GROSS_EROSION_KG_M2_YR_MEAN", 
  "HARD_CORAL_COVER_PCT_MEAN", 
  "Est_Porosity", 
  "VERTICAL_ACCRETION_MM_YR"
)

spatial_Imean <- svyby(reformulate(metrics), ~LOCATION, des, svymean, na.rm = TRUE) %>%
  as_tibble() %>%
  left_join(
    site_sw %>% 
      distinct(LOCATION, LOCATIONCODE, LOCATION_LAT) %>%
      mutate(LOCATION = as.character(LOCATION), LOCATIONCODE = as.character(LOCATIONCODE)),
    by = "LOCATION"
  ) %>% 
  rename(LATITUDE = LOCATION_LAT) %>%
  relocate(LOCATION, LOCATIONCODE, LATITUDE)

write_csv(spatial_Imean, out_island_csv)
cat(green(sprintf("✅ Weighted island/location means successfully exported (%d locations):\n   -> %s\n", 
                  nrow(spatial_Imean), out_island_csv)))

cat(green("\n✅ StRS-to-Island scaling complete!\n"))