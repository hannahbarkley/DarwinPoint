# ==============================================================================
# Script:  DarwinPoint_processCB_fixed_cnet.R
# Purpose: Hawaiian Archipelago - Net Carbonate Budget Calculation & Visualization
#          Integrates pre-calculated site-level gross production from StRS data
#          (filtered for fixed sites matching meta_file) with endolithic bioerosion, fixed-site 
#          urchin bioerosion, fixed-site parrotfish bioerosion, and live coral cover to 
#          generate site-level net carbonate budgets and figures.
# Inputs:  - data/fixed data/field data/HA_2024_Fixed_SiteSummary.csv
#          - data/strs data/field data/HA_2024_StRS_Prod_SiteLevel.csv
#          - data/fixed data/field data/HA_2024_Fixed_Urchins.csv
#          - data/fixed data/field data/HA_2024_Fixed_Fish.csv
#          - data/strs data/field data/2024_HA_StRS_Cover_Tier1.csv
#          - data/dbase/prod_dbase_ha_fixed.csv
# Outputs: - data/fixed data/processed data/HA_2024_Fixed_Urchins_TransectLevel.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Urchins_SiteLevel.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Fish_SiteLevel.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_cnet.csv
#          - outputs/supporting figures/HA_2024_gross_overview.png
#          - outputs/supporting figures/HA_2024_net_overview.png
#          - outputs/supporting figures/HA_2024_netprod_fixed.png
#          - outputs/supporting figures/HA_2024_Fixed_netprod_cnet.png
#          - outputs/supporting figures/HA_2024_Fixed_netprod_vs_coralcover_cnet.png
# ==============================================================================

library(tidyverse)
library(lubridate)
library(data.table)
library(stringr)
library(abjutils)
library(reefbudgetR)
library(crayon)
library(ggrepel)
library(ggpubr)
library(RColorBrewer)
library(here)

# ==============================================================================
# 1. PARAMETERS & DATA INGESTION
# ==============================================================================

region <- "HA"
year   <- "2024"

# Base directory relative to project root
base_dir     <- here()
data_dir     <- file.path(base_dir, "data")
plot_dir     <- file.path(base_dir, "outputs")
plot_out_dir <- file.path(base_dir, "outputs", "supporting figures")
proc_out_dir <- file.path(data_dir, "fixed data", "processed data")

# Master input files
meta_file       <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_SiteSummary.csv")
benthic_file    <- file.path(data_dir, "strs data", "field data", "HA_2024_StRS_Prod_SiteLevel.csv")
urchin_file     <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_Urchins.csv")
fish_file       <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_Fish.csv")
prod_dbase_file <- file.path(data_dir, "dbase", "prod_dbase_ha_fixed.csv")

# Dynamic resolution of coral cover file path
cover_candidates <- c(
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_Cover_Tier1.csv"),
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_CoralCover_Tier1.csv"),
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_CoralCover.csv")
)
cover_file <- cover_candidates[file.exists(cover_candidates)][1]
if (is.na(cover_file)) cover_file <- cover_candidates[1]

# Processed intermediate and final paths
urch_site_file  <- file.path(proc_out_dir, sprintf("%s_%s_Fixed_Urchins_SiteLevel.csv", region, year))
fish_site_file  <- file.path(proc_out_dir, sprintf("%s_%s_Fixed_Fish_SiteLevel.csv", region, year))
out_net_file    <- file.path(proc_out_dir, sprintf("%s_%s_Fixed_Net_SiteLevel_cnet.csv", region, year))

stopifnot("Metadata file not found!"     = file.exists(meta_file))
stopifnot("Benthic StRS file not found!" = file.exists(benthic_file))

sites_to_exclude <- c("OCC-MOL-003")

# Ensure output directories exist
invisible(lapply(c(data_dir, proc_out_dir, plot_dir, plot_out_dir), function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}))

# Helper to merge site metadata without across() scoping conflicts
merge_metadata <- function(df, meta) {
  meta_cols <- setdiff(names(meta), "OCC_SITEID")
  res <- left_join(df, meta, by = "OCC_SITEID", suffix = c("", ".meta"))
  for (col in meta_cols) {
    meta_col <- paste0(col, ".meta")
    if (meta_col %in% names(res)) {
      res[[col]] <- coalesce(res[[col]], res[[meta_col]])
      res[[meta_col]] <- NULL
    }
  }
  res
}

# Dynamic plot axis limits
calc_dynamic_lims <- function(vec_min, vec_max, pad = 0.12, force_zero_lower = FALSE) {
  min_val <- min(vec_min, na.rm = TRUE)
  max_val <- max(vec_max, na.rm = TRUE)
  
  if (force_zero_lower) {
    low <- 0
  } else {
    low <- if (min_val < 0) min_val - abs(min_val) * pad else 0
  }
  high <- if (max_val > 0) max_val + abs(max_val) * pad else 1
  
  c(floor(low), ceiling(high))
}

# Format Hawaiian location names with proper ʻOkina and Kahakō
format_okinas <- function(df, levels_vec) {
  df$LOCATION <- factor(df$LOCATION, levels = levels_vec)
  df$labels   <- fct_recode(df$LOCATION,
                            "Ni\u02bbihau"    = "Niihau",
                            "Kaua\u02bbi"     = "Kauai",
                            "O\u02bbahu"      = "Oahu",
                            "Moloka\u02bbi"   = "Molokai",
                            "Kaho\u02bbolawe" = "Kahoolawe",
                            "L\u0101na\u02bbi"   = "Lanai", 
                            "Hawai\u02bbi"    = "Hawaii")
  return(df)
}

# ==============================================================================
# 2. FIXED FISH BIOEROSION PIPELINE
# ==============================================================================

get_fixed_fish_data <- function(fish_site_file, fish_file, meta_file) {
  cat(cyan(sprintf("\n🚀 LOADING FIXED FISH BIOEROSION: %s %s\n", region, year)))
  
  if (file.exists(fish_site_file)) {
    cat(blue(sprintf("📂 Reading existing processed fish data from: %s\n", basename(fish_site_file))))
    fish_df <- fread(fish_site_file, na.strings = c("", "NA"), data.table = FALSE)
  } else {
    stopifnot("Fish field file not found to process fish data!" = file.exists(fish_file))
    cat(yellow(sprintf("⚠️ Processing raw field data from: %s\n", basename(fish_file))))
    sites_metadata <- fread(meta_file, data.table = FALSE)
    spc_data       <- fread(fish_file, na.strings = c("", "NA"), data.table = FALSE)
    
    fish_df <- process_fish(
      spc_data       = spc_data, 
      method_type    = "Fixed", 
      fixed_metadata = sites_metadata, 
      dbase_type     = "Kindinger"
    )
    fish_df$YEAR <- as.character(year)
    fwrite(fish_df, fish_site_file, na = "")
  }
  
  site_col <- intersect(c("OCC_SITEID", "SITE", "REA_SITEID"), names(fish_df))[1]
  
  # Locate Erosion Column
  eros_col <- intersect(c("FISH_EROSION_KG_M2_YR_ALL_MEAN", "FISH_EROSION_KG_M2_YR_ALL"), names(fish_df))
  if (length(eros_col) == 0) eros_col <- grep("EROSION.*MEAN", names(fish_df), value = TRUE, ignore.case = TRUE)
  if (length(eros_col) == 0) eros_col <- grep("EROSION.*ALL", names(fish_df), value = TRUE, ignore.case = TRUE)
  if (length(eros_col) == 0) eros_col <- grep("EROSION", names(fish_df), value = TRUE, ignore.case = TRUE)
  eros_col <- eros_col[1]
  
  # Locate Biomass Column
  bio_col <- intersect(c("FISH_BIOMASS_KG_HA_ALL_MEAN", "FISH_BIOMASS_KG_HA_ALL"), names(fish_df))
  if (length(bio_col) == 0) bio_col <- grep("BIOMASS.*ALL", names(fish_df), value = TRUE, ignore.case = TRUE)
  bio_col <- bio_col[1]
  
  # Locate Density Column
  den_col <- intersect(c("FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN", "FISH_DENSITY_ABUNDANCE_HA_ALL"), names(fish_df))
  if (length(den_col) == 0) den_col <- grep("DENSITY.*ALL", names(fish_df), value = TRUE, ignore.case = TRUE)
  den_col <- den_col[1]
  
  stopifnot("Could not locate a fish bioerosion rate column!" = !is.na(eros_col))
  
  fish_clean <- fish_df %>%
    transmute(
      OCC_SITEID                         = trimws(as.character(.data[[site_col]])),
      FISH_EROSION_KG_M2_YR_ALL_MEAN     = as.numeric(.data[[eros_col]]),
      FISH_EROSION_KG_M2_YR_ALL_SE       = if ("FISH_EROSION_KG_M2_YR_ALL_SE" %in% names(.)) as.numeric(FISH_EROSION_KG_M2_YR_ALL_SE) else 0,
      FISH_BIOMASS_KG_HA_ALL_MEAN        = if (!is.na(bio_col)) as.numeric(.data[[bio_col]]) else 0,
      FISH_BIOMASS_KG_HA_ALL_SE          = if ("FISH_BIOMASS_KG_HA_ALL_SE" %in% names(.)) as.numeric(FISH_BIOMASS_KG_HA_ALL_SE) else 0,
      FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN = if (!is.na(den_col)) as.numeric(.data[[den_col]]) else 0,
      FISH_DENSITY_ABUNDANCE_HA_ALL_SE   = if ("FISH_DENSITY_ABUNDANCE_HA_ALL_SE" %in% names(.)) as.numeric(FISH_DENSITY_ABUNDANCE_HA_ALL_SE) else 0
    ) %>%
    distinct(OCC_SITEID, .keep_all = TRUE)
  
  cat(green(sprintf("✅ Loaded fish metrics for %d fixed sites.\n", nrow(fish_clean))))
  return(fish_clean)
}

# ==============================================================================
# 3. FIXED URCHIN BIOEROSION PIPELINE
# ==============================================================================

get_fixed_urchin_data <- function(urch_site_file, urchin_file, meta_file) {
  cat(cyan(sprintf("\n🚀 LOADING FIXED URCHIN BIOEROSION: %s %s\n", region, year)))
  
  if (file.exists(urch_site_file)) {
    cat(blue(sprintf("📂 Reading existing processed urchin data from: %s\n", basename(urch_site_file))))
    urch_df <- fread(urch_site_file, na.strings = c("", "NA"), data.table = FALSE)
  } else {
    stopifnot("Urchin field file not found to process urchin data!" = file.exists(urchin_file))
    cat(yellow(sprintf("⚠️ Processing raw field data from: %s\n", basename(urchin_file))))
    sites_metadata      <- fread(meta_file, data.table = FALSE)
    urch_data           <- fread(urchin_file, na.strings = c("", "NA"), data.table = FALSE)
    urch_data$LOCALDATE <- lubridate::ymd(urch_data$LOCALDATE)
    
    urch_data     <- prep_urchins(data = urch_data, sites_metadata = sites_metadata)
    urch_results  <- process_urchins(data = urch_data, full_summary = TRUE)
    urch_df       <- as.data.frame(urch_results$site_erosion)
    urch_transect <- as.data.frame(urch_results$transect_erosion)
    
    fwrite(urch_transect, file.path(dirname(urch_site_file), sprintf("%s_%s_Fixed_Urchins_TransectLevel.csv", region, year)), na = "")
    fwrite(urch_df, urch_site_file, na = "")
  }
  
  site_col <- intersect(c("OCC_SITEID", "SITE"), names(urch_df))[1]
  
  u_eros_candidates <- grep("EROSION|BIOEROSION", names(urch_df), value = TRUE, ignore.case = TRUE)
  u_eros_candidates <- setdiff(u_eros_candidates, grep("DENSITY", u_eros_candidates, value = TRUE, ignore.case = TRUE))
  u_eros_col <- intersect(c("URCHIN_EROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR"), u_eros_candidates)
  if (length(u_eros_col) == 0) u_eros_col <- u_eros_candidates
  u_eros_col <- u_eros_col[1]
  
  u_den_candidates <- grep("DENSITY", names(urch_df), value = TRUE, ignore.case = TRUE)
  u_den_col <- intersect(c("URCHIN_DENSITY_NO_M2_MEAN", "URCHIN_DENSITY_NO_M2"), u_den_candidates)
  if (length(u_den_col) == 0) u_den_col <- u_den_candidates
  u_den_col <- u_den_col[1]
  
  urch_clean <- urch_df %>%
    transmute(
      OCC_SITEID                   = trimws(as.character(.data[[site_col]])),
      LOCATION                     = if ("LOCATION" %in% names(.)) LOCATION else if ("ISLAND" %in% names(.)) ISLAND else NA_character_,
      URCHIN_EROSION_KG_M2_YR_MEAN = if (!is.na(u_eros_col)) as.numeric(.data[[u_eros_col]]) else 0,
      URCHIN_EROSION_KG_M2_YR_SE   = if ("URCHIN_EROSION_KG_M2_YR_SE" %in% names(.)) as.numeric(URCHIN_EROSION_KG_M2_YR_SE) else 0,
      URCHIN_DENSITY_NO_M2_MEAN    = if (!is.na(u_den_col))  as.numeric(.data[[u_den_col]])  else 0,
      URCHIN_DENSITY_NO_M2_SE      = if ("URCHIN_DENSITY_NO_M2_SE" %in% names(.)) as.numeric(URCHIN_DENSITY_NO_M2_SE) else 0
    ) %>%
    distinct(OCC_SITEID, .keep_all = TRUE)
  
  cat(green(sprintf("✅ Loaded urchin bioerosion (column '%s') for %d fixed sites.\n", u_eros_col, nrow(urch_clean))))
  return(urch_clean)
}

# ==============================================================================
# 4. CORAL COVER INGESTION & HARMONIZATION (PT / PCT FILTERED FOR TIER_1 == 'CORAL')
# ==============================================================================

get_coral_cover_data <- function(cover_file) {
  cat(cyan(sprintf("\n🚀 LOADING CORAL COVER OBSERVATIONS: %s\n", basename(cover_file))))
  
  if (!file.exists(cover_file)) {
    cat(yellow("⚠️ Coral cover file not found. Setting cover to NA.\n"))
    return(data.frame(OCC_SITEID = character(0), CORAL_COVER_PERC = numeric(0)))
  }
  
  cover_raw <- fread(cover_file, na.strings = c("", "NA"), data.table = FALSE)
  
  site_col <- intersect(c("OCC_SITEID", "SITE", "REA_SITEID"), names(cover_raw))[1]
  pct_col  <- intersect(c("PCT", "PT", "CORAL_COVER_PERC", "HARD_CORAL_COVER_PCT", "COVER"), names(cover_raw))[1]
  if (is.na(pct_col)) {
    pct_col <- grep("PCT|PT|COVER|PERC", names(cover_raw), value = TRUE, ignore.case = TRUE)[1]
  }
  tier_col <- intersect(c("TIER_1", "TIER1", "CATEGORY", "BENTHIC_CATEGORY", "TAXA"), names(cover_raw))[1]
  
  stopifnot("Could not locate a site column in coral cover file!" = !is.na(site_col))
  stopifnot("Could not locate a coral cover percentage column!" = !is.na(pct_col))
  
  # Filter strictly for hard coral tier to avoid summing or averaging across all benthic classes
  if (!is.na(tier_col)) {
    cover_raw <- cover_raw %>%
      filter(grepl("^CORAL$", trimws(toupper(.data[[tier_col]]))))
  }
  
  cover_clean <- cover_raw %>%
    transmute(
      OCC_SITEID       = trimws(as.character(.data[[site_col]])),
      CORAL_COVER_PERC = as.numeric(.data[[pct_col]])
    ) %>%
    filter(!is.na(OCC_SITEID) & !is.na(CORAL_COVER_PERC)) %>%
    group_by(OCC_SITEID) %>%
    summarize(CORAL_COVER_PERC = mean(CORAL_COVER_PERC, na.rm = TRUE), .groups = "drop")
  
  cat(green(sprintf("✅ Loaded coral cover for %d sites (Mean: %.2f%%).\n", 
                    nrow(cover_clean), mean(cover_clean$CORAL_COVER_PERC, na.rm = TRUE))))
  return(cover_clean)
}

# ==============================================================================
# 5. NET BUDGET PIPELINE (StRS PROD + BIOEROSION + CORAL COVER)
# ==============================================================================

run_fixed_net_pipeline <- function(benthic_file, meta_file, fish_clean, urch_clean, cover_clean, out_net_file, exclude_sites = NULL) {
  cat(cyan(sprintf("\n🚀 STARTING NET BUDGET PIPELINE (FIXED SITES via process_net): %s %s\n", region, year)))
  cat(blue("📂 Loading StRS production file and pulling matched fixed sites...\n"))
  
  strs_data      <- fread(benthic_file, na.strings = c("", "NA"), data.table = FALSE)
  sites_metadata <- fread(meta_file, na.strings = c("", "NA"), data.table = FALSE)
  
  # Ensure OCC_SITEID column exists in metadata even if named SITE
  meta_site_col <- intersect(c("OCC_SITEID", "SITE", "REA_SITEID"), names(sites_metadata))[1]
  if (!"OCC_SITEID" %in% names(sites_metadata) && !is.na(meta_site_col)) {
    sites_metadata$OCC_SITEID <- sites_metadata[[meta_site_col]]
  }
  target_site_ids <- unique(na.omit(trimws(sites_metadata$OCC_SITEID)))
  
  meta_to_join <- sites_metadata %>% 
    select(any_of(c("REGION", "REGIONCODE", "ISLAND", "LOCATION", "LOCATIONCODE", 
                    "YEAR", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "OCC_SITEID", "LOCALDATE"))) %>%
    mutate(across(any_of("YEAR"), as.character)) %>%
    distinct(OCC_SITEID, .keep_all = TRUE)
  
  if (!"LOCATION" %in% names(meta_to_join) && "ISLAND" %in% names(meta_to_join)) {
    meta_to_join$LOCATION <- meta_to_join$ISLAND
  }
  
  meta_cols_order <- c("REGION", "REGIONCODE", "YEAR", "CRUISE_ID", "LOCATION", "LOCATIONCODE",
                       "OCC_SITEID", "SITEVISITID", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "LOCALDATE")
  
  # 1. Format Prod Data for process_net() & Populate Bioerosion Mean & SE
  matched_fixed <- strs_data %>%
    filter(trimws(SITE) %in% target_site_ids) %>%
    mutate(
      OCC_SITEID                    = trimws(SITE),
      GROSS_CARB_PROD_KG_M2_YR_MEAN = as.numeric(PROD),
      GROSS_CARB_PROD_KG_M2_YR_SE   = if (all(c("PROD_upr", "PROD_lwr") %in% names(.))) {
        coalesce(abs(as.numeric(PROD_upr) - as.numeric(PROD_lwr)) / (2 * 1.96), 0)
      } else { 0 },
      GROSS_CARB_PROD_KG_M2_YR_N    = 1, 
      # Ingestion of MICROMACRO_TOTAL and SE calculation from MICROMACRO_upr / MICROMACRO_lwr
      BIOEROSION_KG_M2_YR_MEAN      = if ("MICROMACRO_TOTAL" %in% names(.)) as.numeric(MICROMACRO_TOTAL) else 0,
      BIOEROSION_KG_M2_YR_SE        = if (all(c("MICROMACRO_upr", "MICROMACRO_lwr") %in% names(.))) {
        coalesce(abs(as.numeric(MICROMACRO_upr) - as.numeric(MICROMACRO_lwr)) / (2 * 1.96), 0)
      } else { 0 },
      MACROBIOEROSION_KG_M2_YR_MEAN = if ("MICROMACRO_TOTAL" %in% names(.)) as.numeric(MICROMACRO_TOTAL) else 0,
      MACROBIOEROSION_KG_M2_YR_SE   = if (all(c("MICROMACRO_upr", "MICROMACRO_lwr") %in% names(.))) {
        coalesce(abs(as.numeric(MICROMACRO_upr) - as.numeric(MICROMACRO_lwr)) / (2 * 1.96), 0)
      } else { 0 },
      MICROBIOEROSION_KG_M2_YR_MEAN = 0,
      MICROBIOEROSION_KG_M2_YR_SE   = 0,
      YEAR                          = as.character(year)
    ) %>%
    merge_metadata(meta_to_join) %>%
    left_join(cover_clean, by = "OCC_SITEID") %>%
    mutate(
      HARD_CORAL_COVER_PCT_MEAN = coalesce(CORAL_COVER_PERC, 0),
      HARD_CORAL_COVER_PCT_SE   = 0
    )
  
  if (!"REGIONCODE" %in% names(matched_fixed) && "REGION" %in% names(matched_fixed)) {
    matched_fixed$REGIONCODE <- matched_fixed$REGION
  }
  
  matched_fixed <- matched_fixed %>% relocate(any_of(meta_cols_order))
  
  # 2. Format Urchin Data (incorporating Location Mean fallback) for process_net()
  urch_isl <- urch_clean %>%
    filter(!is.na(LOCATION) & !is.na(URCHIN_EROSION_KG_M2_YR_MEAN)) %>%
    group_by(LOCATION) %>%
    summarize(URCHIN_ISL_MEAN = mean(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE), .groups = "drop")
  
  urch_formatted <- urch_clean %>%
    left_join(meta_to_join %>% select(OCC_SITEID, LOCATION_meta = any_of("LOCATION")), by = "OCC_SITEID") %>%
    mutate(LOCATION = if ("LOCATION_meta" %in% names(.)) coalesce(LOCATION, LOCATION_meta) else LOCATION) %>%
    left_join(urch_isl, by = "LOCATION") %>%
    mutate(
      URCHIN_EROSION_KG_M2_YR_MEAN = coalesce(URCHIN_EROSION_KG_M2_YR_MEAN, URCHIN_ISL_MEAN, 0),
      URCHIN_EROSION_KG_M2_YR_SE   = coalesce(URCHIN_EROSION_KG_M2_YR_SE, 0),
      URCHIN_DENSITY_NO_M2_MEAN    = coalesce(URCHIN_DENSITY_NO_M2_MEAN, 0),
      URCHIN_DENSITY_NO_M2_SE      = coalesce(URCHIN_DENSITY_NO_M2_SE, 0),
      YEAR                         = as.character(year)
    ) %>%
    select(-any_of(c("LOCATION_meta", "URCHIN_ISL_MEAN"))) %>%
    relocate(any_of(meta_cols_order))
  
  # 3. Format Fish Data for process_net()
  fish_formatted <- fish_clean %>% 
    mutate(YEAR = as.character(year)) %>%
    relocate(any_of(meta_cols_order))
  
  cat(blue("🔗 Joining Fish, Urchins, and Processing Net Budgets via process_net()...\n"))
  
  # 4. Run process_net()
  comb_all <- process_net(
    prod   = matched_fixed,
    urch   = urch_formatted,
    fish   = fish_formatted,
    sum_by = "site",
    format = "wide"
  )
  
  # Defensive normalization: Ensure OCC_SITEID is present for all subsequent plotting
  if ("SITE" %in% names(comb_all) && !"OCC_SITEID" %in% names(comb_all)) {
    comb_all <- comb_all %>% rename(OCC_SITEID = SITE)
  } else if ("SITE" %in% names(comb_all) && "OCC_SITEID" %in% names(comb_all)) {
    comb_all$OCC_SITEID <- coalesce(comb_all$OCC_SITEID, comb_all$SITE)
  }
  
  # Ensure geographic coordinates are preserved post-process_net
  if (!"LATITUDE" %in% names(comb_all) || all(is.na(comb_all$LATITUDE))) {
    comb_all <- comb_all %>%
      select(-any_of(c("LATITUDE", "LONGITUDE"))) %>%
      left_join(
        matched_fixed %>% select(OCC_SITEID, LATITUDE, LONGITUDE) %>% distinct(),
        by = "OCC_SITEID"
      )
  }
  
  # Preserve explicit internal bioerosion mean and SE columns
  if (!"BIOEROSION_KG_M2_YR_SE" %in% names(comb_all)) {
    comb_all <- comb_all %>%
      select(-any_of(c("BIOEROSION_KG_M2_YR_MEAN", "BIOEROSION_KG_M2_YR_SE"))) %>%
      left_join(
        matched_fixed %>% select(OCC_SITEID, BIOEROSION_KG_M2_YR_MEAN, BIOEROSION_KG_M2_YR_SE) %>% distinct(),
        by = "OCC_SITEID"
      )
  }
  
  # Preserve live coral cover percentage if dropped during process_net
  if (!"HARD_CORAL_COVER_PCT_MEAN" %in% names(comb_all) || all(is.na(comb_all$HARD_CORAL_COVER_PCT_MEAN))) {
    comb_all <- comb_all %>%
      select(-any_of(c("HARD_CORAL_COVER_PCT_MEAN", "HARD_CORAL_COVER_PCT_SE"))) %>%
      left_join(
        matched_fixed %>% select(OCC_SITEID, HARD_CORAL_COVER_PCT_MEAN, HARD_CORAL_COVER_PCT_SE) %>% distinct(),
        by = "OCC_SITEID"
      )
  }
  
  # Calculate Net Carbonate Production Mean and propagate SE via error quadrature
  comb_all <- comb_all %>%
    mutate(
      GROSS_EROSION_KG_M2_YR_MEAN = coalesce(BIOEROSION_KG_M2_YR_MEAN, 0) + 
        coalesce(URCHIN_EROSION_KG_M2_YR_MEAN, 0) + 
        coalesce(FISH_EROSION_KG_M2_YR_ALL_MEAN, 0),
      GROSS_EROSION_KG_M2_YR_SE   = sqrt(coalesce(BIOEROSION_KG_M2_YR_SE, 0)^2 + 
                                           coalesce(URCHIN_EROSION_KG_M2_YR_SE, 0)^2 + 
                                           coalesce(FISH_EROSION_KG_M2_YR_ALL_SE, 0)^2),
      NET_CARB_PROD_KG_M2_YR_MEAN = GROSS_CARB_PROD_KG_M2_YR_MEAN - GROSS_EROSION_KG_M2_YR_MEAN,
      NET_CARB_PROD_KG_M2_YR_SE   = sqrt(coalesce(GROSS_CARB_PROD_KG_M2_YR_SE, 0)^2 + 
                                           coalesce(BIOEROSION_KG_M2_YR_SE, 0)^2 + 
                                           coalesce(URCHIN_EROSION_KG_M2_YR_SE, 0)^2 + 
                                           coalesce(FISH_EROSION_KG_M2_YR_ALL_SE, 0)^2)
    )
  
  # Safeguard: Explicitly force all metric columns back to numeric
  metric_cols <- c("GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_CARB_PROD_KG_M2_YR_SE", 
                   "BIOEROSION_KG_M2_YR_MEAN", "BIOEROSION_KG_M2_YR_SE",
                   "MACROBIOEROSION_KG_M2_YR_MEAN", "MACROBIOEROSION_KG_M2_YR_SE",
                   "MICROBIOEROSION_KG_M2_YR_MEAN", "MICROBIOEROSION_KG_M2_YR_SE",
                   "URCHIN_EROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR_SE", 
                   "URCHIN_DENSITY_NO_M2_MEAN", "URCHIN_DENSITY_NO_M2_SE",
                   "FISH_EROSION_KG_M2_YR_ALL_MEAN", "FISH_EROSION_KG_M2_YR_ALL_SE",
                   "FISH_BIOMASS_KG_HA_ALL_MEAN", "FISH_BIOMASS_KG_HA_ALL_SE",
                   "FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN", "FISH_DENSITY_ABUNDANCE_HA_ALL_SE",
                   "GROSS_EROSION_KG_M2_YR_MEAN", "GROSS_EROSION_KG_M2_YR_SE",
                   "NET_CARB_PROD_KG_M2_YR_MEAN", "NET_CARB_PROD_KG_M2_YR_SE",
                   "HARD_CORAL_COVER_PCT_MEAN")
  
  comb_all <- comb_all %>% mutate(across(any_of(metric_cols), as.numeric))
  
  # 5. Filter excluded sites
  if (!is.null(exclude_sites)) {
    comb_all <- comb_all %>% filter(!OCC_SITEID %in% exclude_sites)
  }
  
  fwrite(comb_all, out_net_file, na = "")
  cat(green(sprintf("✅ Net Budget Pipeline Complete (%d sites). Exported to: %s\n", 
                    nrow(comb_all), out_net_file)))
  return(comb_all)
}

# ==============================================================================
# 6. VISUALIZATION FUNCTIONS
# ==============================================================================

create_benthic_plots <- function(data, site_order = NULL, exclude_sites = NULL) {
  if (!is.null(exclude_sites)) {
    data <- data %>% filter(!OCC_SITEID %in% exclude_sites)
  }
  
  if (!is.null(site_order)) {
    valid_sites <- site_order[site_order %in% unique(data$OCC_SITEID)]
    data$OCC_SITEID <- factor(data$OCC_SITEID, levels = valid_sites)
  } else if ("LATITUDE" %in% names(data)) {
    lat_order <- data %>%
      filter(!is.na(LATITUDE)) %>%
      distinct(OCC_SITEID, LATITUDE) %>%
      arrange(desc(LATITUDE)) %>%
      pull(OCC_SITEID)
    other_sites     <- setdiff(unique(data$OCC_SITEID), lat_order)
    data$OCC_SITEID <- factor(data$OCC_SITEID, levels = c(lat_order, other_sites))
  }
  
  req_cols <- c("GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_CARB_PROD_KG_M2_YR_SE",
                "URCHIN_DENSITY_NO_M2_MEAN",     "URCHIN_DENSITY_NO_M2_SE",
                "URCHIN_EROSION_KG_M2_YR_MEAN",  "URCHIN_EROSION_KG_M2_YR_SE",
                "NET_CARB_PROD_KG_M2_YR_MEAN",    "NET_CARB_PROD_KG_M2_YR_SE")
  for (col in req_cols) {
    if (!col %in% names(data)) data[[col]] <- 0
  }
  
  lims_gross <- calc_dynamic_lims(
    data$GROSS_CARB_PROD_KG_M2_YR_MEAN - coalesce(data$GROSS_CARB_PROD_KG_M2_YR_SE, 0),
    data$GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(data$GROSS_CARB_PROD_KG_M2_YR_SE, 0)
  )
  lims_urch_den <- calc_dynamic_lims(
    0, data$URCHIN_DENSITY_NO_M2_MEAN + coalesce(data$URCHIN_DENSITY_NO_M2_SE, 0), force_zero_lower = TRUE
  )
  lims_urch_eros <- calc_dynamic_lims(
    0, data$URCHIN_EROSION_KG_M2_YR_MEAN + coalesce(data$URCHIN_EROSION_KG_M2_YR_SE, 0), force_zero_lower = TRUE
  )
  lims_net <- calc_dynamic_lims(
    data$NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(data$NET_CARB_PROD_KG_M2_YR_SE, 0),
    data$NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(data$NET_CARB_PROD_KG_M2_YR_SE, 0)
  )
  
  pal <- c("#003087", "#0085ca", "#1ebec7", "#76bc21", "#f0de02", "#ff8400", "#db2207")
  
  theme_benthic <- theme_bw() + 
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 8, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y      = element_text(size = 12, color = 'black'),
      axis.title.y     = element_text(size = 12),
      strip.text.x     = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position  = 'none' 
    )
  
  add_bar_layers <- function(y_var, y_se_var, fill_var = "REGIONCODE") {
    list(
      geom_hline(yintercept = 0),
      geom_col(aes(x = .data[["OCC_SITEID"]], y = .data[[y_var]], fill = .data[[fill_var]]),
               color = 'black', position = position_dodge(0.6), width = 0.6),
      geom_errorbar(aes(x = .data[["OCC_SITEID"]], 
                        ymin = .data[[y_var]] - coalesce(.data[[y_se_var]], 0), 
                        ymax = .data[[y_var]] + coalesce(.data[[y_se_var]], 0)),
                    width = 0.2, position = position_dodge(0.6)),
      facet_grid(~REGIONCODE, scales = "free", space = "free"),
      scale_fill_manual(values = pal)
    )
  }
  
  gross_plot <- ggplot(data) +
    add_bar_layers("GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_CARB_PROD_KG_M2_YR_SE") +
    scale_y_continuous(limits = lims_gross, expand = c(0, 0)) +
    ylab(expression(paste("Gross production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_benthic
  
  urch_den_plot <- ggplot(data) +
    add_bar_layers("URCHIN_DENSITY_NO_M2_MEAN", "URCHIN_DENSITY_NO_M2_SE") +
    scale_y_continuous(limits = lims_urch_den, expand = c(0, 0)) +
    ylab(expression(paste("Urchin density (no", " m" ^ -2, ")  "))) +
    theme_benthic +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 
  
  urch_eros_plot <- ggplot(data) +
    add_bar_layers("URCHIN_EROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR_SE") +
    scale_y_continuous(limits = lims_urch_eros, expand = c(0, 0)) +
    ylab(expression(paste("Urchin erosion rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_benthic +
    theme(strip.text.x = element_blank()) 
  
  net_plot <- ggplot(data) +
    add_bar_layers("NET_CARB_PROD_KG_M2_YR_MEAN", "NET_CARB_PROD_KG_M2_YR_SE") +
    scale_y_continuous(limits = lims_net, expand = c(0, 0)) +
    ylab(expression(paste("Net production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_benthic
  
  urch_combined <- ggarrange(urch_den_plot, urch_eros_plot, nrow = 2, align = "v", heights = c(1, 1.2))
  
  return(list(gross_plot = gross_plot, urch_plot = urch_combined, net_plot = net_plot))
}

format_bygroup <- function(net_data, prod_group) {
  req_eros <- c("BIOEROSION_KG_M2_YR_MEAN", "MACROBIOEROSION_KG_M2_YR_MEAN", 
                "MICROBIOEROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR_MEAN", 
                "FISH_EROSION_KG_M2_YR_ALL_MEAN")
  for (col in req_eros) {
    if (!col %in% names(net_data)) net_data[[col]] <- 0
  }
  
  prod_clean <- prod_group %>% 
    select(any_of(c("REGION", "REGIONCODE", "LOCATION", "OCC_SITEID", "LATITUDE", "GROUP", "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN")))
  
  has_combined_bioeros <- any(net_data$BIOEROSION_KG_M2_YR_MEAN > 0, na.rm = TRUE)
  
  if (has_combined_bioeros) {
    eros_clean <- net_data %>%
      select(any_of(c("REGION", "REGIONCODE", "LOCATION", "OCC_SITEID", "LATITUDE")), 
             `Internal Bioerosion` = BIOEROSION_KG_M2_YR_MEAN,
             Urchins               = URCHIN_EROSION_KG_M2_YR_MEAN,
             Parrotfish            = FISH_EROSION_KG_M2_YR_ALL_MEAN) %>%
      pivot_longer(cols = c(`Internal Bioerosion`, Urchins, Parrotfish),
                   names_to = "GROUP", 
                   values_to = "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN") %>%
      mutate(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = -1 * SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN)
  } else {
    eros_clean <- net_data %>%
      select(any_of(c("REGION", "REGIONCODE", "LOCATION", "OCC_SITEID", "LATITUDE")), 
             Macroborers = MACROBIOEROSION_KG_M2_YR_MEAN,
             Microborers = MICROBIOEROSION_KG_M2_YR_MEAN,
             Urchins     = URCHIN_EROSION_KG_M2_YR_MEAN,
             Parrotfish  = FISH_EROSION_KG_M2_YR_ALL_MEAN) %>%
      pivot_longer(cols = c(Macroborers, Microborers, Urchins, Parrotfish),
                   names_to = "GROUP", 
                   values_to = "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN") %>%
      mutate(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = -1 * SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN)
  }
  
  bind_rows(prod_clean, eros_clean)
}

create_budget_composite_plot <- function(net_site, region, year, plot_out_dir, site_order = NULL, exclude_sites = NULL) {
  cat(blue("\n🎨 Generating Net Budget Gradient and Components Composite Plot...\n"))
  
  net_site_filtered <- net_site
  if (!is.null(exclude_sites)) {
    net_site_filtered <- net_site_filtered %>% filter(!OCC_SITEID %in% exclude_sites)
  }
  
  prod_group <- net_site_filtered %>%
    select(any_of(c("REGION", "REGIONCODE", "LOCATION", "OCC_SITEID", "LATITUDE")), 
           SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = GROSS_CARB_PROD_KG_M2_YR_MEAN) %>%
    mutate(GROUP = "Gross Production")
  
  data_group <- format_bygroup(net_site_filtered, prod_group)
  all_groups <- c("Gross Production", "Internal Bioerosion", "Macroborers", "Microborers", "Urchins", "Parrotfish")
  pal_budget <- c(
    "Gross Production"    = "#2f2e94", 
    "Internal Bioerosion" = "#fd903b", 
    "Macroborers"         = "#fccf5c", 
    "Microborers"         = "#fd903b", 
    "Urchins"             = "#ea6635", 
    "Parrotfish"          = "#d8361d"
  )
  
  present_groups <- intersect(all_groups, unique(data_group$GROUP))
  data_group$GROUP <- factor(data_group$GROUP, levels = present_groups)
  
  if (!is.null(site_order)) {
    valid_order <- site_order[site_order %in% unique(data_group$OCC_SITEID)]
  } else if ("LATITUDE" %in% names(net_site_filtered)) {
    lat_order <- net_site_filtered %>%
      filter(!is.na(LATITUDE)) %>%
      distinct(OCC_SITEID, LATITUDE) %>%
      arrange(desc(LATITUDE)) %>%
      pull(OCC_SITEID)
    other_sites <- setdiff(unique(net_site_filtered$OCC_SITEID), lat_order)
    valid_order <- c(lat_order, other_sites)
  } else {
    valid_order <- unique(data_group$OCC_SITEID)
  }
  
  data_group$OCC_SITEID        <- factor(data_group$OCC_SITEID, levels = valid_order)
  net_site_filtered$OCC_SITEID <- factor(net_site_filtered$OCC_SITEID, levels = valid_order)
  
  max_prod_bar <- max(net_site_filtered$GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(net_site_filtered$GROSS_CARB_PROD_KG_M2_YR_SE, 0), na.rm = TRUE)
  min_eros_bar <- min(-1 * (net_site_filtered$GROSS_EROSION_KG_M2_YR_MEAN + coalesce(net_site_filtered$GROSS_EROSION_KG_M2_YR_SE, 0)), na.rm = TRUE)
  min_net_val  <- min(net_site_filtered$NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(net_site_filtered$NET_CARB_PROD_KG_M2_YR_SE, 0), na.rm = TRUE)
  max_net_val  <- max(net_site_filtered$NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(net_site_filtered$NET_CARB_PROD_KG_M2_YR_SE, 0), na.rm = TRUE)
  
  overall_min   <- min(min_eros_bar, min_net_val, 0, na.rm = TRUE)
  overall_max   <- max(max_prod_bar, max_net_val, 0, na.rm = TRUE)
  budget_y_lims <- c(floor(overall_min - abs(overall_min) * 0.12), ceiling(overall_max + abs(overall_max) * 0.12))
  
  theme_budget <- theme_bw() +
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 14, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y      = element_text(size = 16, color = 'black'),
      axis.title.y     = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      legend.position  = "none" 
    )
  
  prod_bycomponent_plot <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_col(data = data_group, 
             aes(x = OCC_SITEID, y = SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, fill = GROUP), width = 0.5) +
    geom_col(data = net_site_filtered, 
             aes(x = OCC_SITEID, y = GROSS_CARB_PROD_KG_M2_YR_MEAN), 
             alpha = 0, color = "black", width = 0.5) +
    geom_col(data = net_site_filtered, 
             aes(x = OCC_SITEID, y = -1 * GROSS_EROSION_KG_M2_YR_MEAN), 
             alpha = 0, color = "black", width = 0.5) +
    geom_errorbar(data = net_site_filtered,
                  aes(x = OCC_SITEID, 
                      ymin = GROSS_CARB_PROD_KG_M2_YR_MEAN, 
                      ymax = GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(GROSS_CARB_PROD_KG_M2_YR_SE, 0)), 
                  width = 0.25) +
    geom_errorbar(data = net_site_filtered,
                  aes(x = OCC_SITEID, 
                      ymin = -1 * GROSS_EROSION_KG_M2_YR_MEAN - coalesce(GROSS_EROSION_KG_M2_YR_SE, 0), 
                      ymax = -1 * GROSS_EROSION_KG_M2_YR_MEAN), 
                  width = 0.25) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = budget_y_lims, expand = c(0, 0)) +
    scale_fill_manual(values = pal_budget) +
    ylab(expression(paste("Production / bioerosion (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_budget
  
  net_plot_gradient <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_col(data = net_site_filtered,
             aes(x = OCC_SITEID, y = NET_CARB_PROD_KG_M2_YR_MEAN, fill = NET_CARB_PROD_KG_M2_YR_MEAN),
             color = 'black', width = 0.5) +
    geom_errorbar(data = net_site_filtered,
                  aes(x = OCC_SITEID, 
                      ymin = NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(NET_CARB_PROD_KG_M2_YR_SE, 0), 
                      ymax = NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(NET_CARB_PROD_KG_M2_YR_SE, 0)), 
                  width = 0.25) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = budget_y_lims, expand = c(0, 0)) +
    scale_fill_gradient2(low = "#EF2820", mid = "#CCCCCC", high = "#38A1D0", midpoint = 0, limits = budget_y_lims) +
    ylab(expression(paste("Net carbonate production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_budget
  
  legend_components <- get_legend(
    ggplot(data_group, aes(x = OCC_SITEID, y = SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, fill = GROUP)) +
      geom_col(color = "black") +
      scale_fill_manual(values = pal_budget, name = "Component") +
      guides(fill = guide_legend(byrow = TRUE)) +
      theme(legend.text = element_text(size = 14), legend.title = element_text(size = 18), 
            legend.spacing.y = unit(0.3, 'cm'), legend.key.size = unit(0.8, "cm"))
  )
  
  final_plot <- ggarrange(
    net_plot_gradient, NULL, prod_bycomponent_plot, legend_components,
    widths     = c(1, 0.05, 1, 0.25),
    nrow       = 1,
    labels     = c("a", "", "b", ""),
    font.label = list(size = 24, face = "bold")
  )
  
  out_path <- file.path(plot_out_dir, sprintf("%s_%s_netprod_fixed.png", region, year))
  ggsave(final_plot, filename = out_path, bg = "white", height = 12, width = 24, dpi = 300)
  
  cat(green(sprintf("✅ Composite plot saved to: %s\n", out_path)))
  return(final_plot)
}

create_visualizations <- function(comb_all, plot_dir, region, year) {
  cat(cyan("\n🎨 GENERATING NET BUDGET AND CORAL COVER FIGURES...\n"))
  
  location_levels <- c("Hawaii", "Kahoolawe", "Lanai", "Maui", "Molokai", "Oahu", 
                       "Kauai", "Niihau", "French Frigate", "Lisianski", 
                       "Pearl & Hermes", "Kure")
  
  comb_all <- format_okinas(comb_all, location_levels)
  
  y_min <- min(comb_all$NET_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE)
  y_max <- max(comb_all$NET_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE)
  y_limits <- c(min(-8, floor(y_min - 1)), max(16, ceiling(y_max + 1)))
  
  has_cover <- "HARD_CORAL_COVER_PCT_MEAN" %in% names(comb_all) && 
    any(!is.na(comb_all$HARD_CORAL_COVER_PCT_MEAN) & comb_all$HARD_CORAL_COVER_PCT_MEAN > 0)
  
  ha_island_plot <- ggplot(data = comb_all) +
    {
      if (has_cover) {
        geom_jitter(aes(x = labels, y = NET_CARB_PROD_KG_M2_YR_MEAN, 
                        color = NET_CARB_PROD_KG_M2_YR_MEAN, size = HARD_CORAL_COVER_PCT_MEAN), 
                    width = 0.18, height = 0, alpha = 0.85)
      } else {
        geom_jitter(aes(x = labels, y = NET_CARB_PROD_KG_M2_YR_MEAN, 
                        color = NET_CARB_PROD_KG_M2_YR_MEAN), 
                    width = 0.18, height = 0, size = 4, alpha = 0.85)
      }
    } +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    scale_color_gradient2(low = "#EF2820", mid = "#CCCCCC", high = "#38A1D0", 
                          midpoint = 0, limits = y_limits, name = "Net Prod") +
    scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 8), expand = c(0, 0)) +
    ylab(expression(atop(paste("Net carbonate production rate"), paste("(kg CaCO"[3], " m"^-2, " yr"^-1, ")  ")))) +
    theme_bw() +
    theme(
      axis.title.x      = element_blank(),
      axis.text.x       = element_text(size = 12, color = 'black', angle = 30, vjust = 1, hjust = 1),
      axis.text.y       = element_text(size = 16, color = 'black'),
      axis.title.y      = element_text(size = 16, color = 'black'),
      legend.text       = element_text(size = 12, color = 'black'),
      legend.title      = element_text(size = 14, color = 'black'),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = 'transparent'),
      plot.background   = element_rect(fill = 'transparent', color = NA),
      legend.background = element_rect(fill = 'transparent'),
      axis.line         = element_line(colour = "black"),
      axis.ticks        = element_line(colour = "black")
    )
  
  if (has_cover) {
    ha_island_plot <- ha_island_plot + 
      scale_size_continuous(range = c(2.5, 8), name = "Coral Cover (%)") +
      guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))
  }
  
  island_path <- file.path(plot_out_dir, sprintf("%s_%s_Fixed_netprod_cnet.png", region, year))
  ggsave(ha_island_plot, filename = island_path, bg = "white", height = 6, width = 15, dpi = 300)
  cat(green(sprintf("✅ Saved island overview plot to: %s\n", island_path)))
  
  if (has_cover) {
    cover_scatter <- ggplot(data = comb_all %>% filter(!is.na(HARD_CORAL_COVER_PCT_MEAN)), 
                            aes(x = HARD_CORAL_COVER_PCT_MEAN, y = NET_CARB_PROD_KG_M2_YR_MEAN)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_point(aes(color = NET_CARB_PROD_KG_M2_YR_MEAN), size = 4, alpha = 0.85) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", linewidth = 0.8) +
      scale_color_gradient2(low = "#EF2820", mid = "#CCCCCC", high = "#38A1D0", 
                            midpoint = 0, limits = y_limits) +
      xlab("Live Coral Cover (%)") +
      ylab(expression(paste("Net carbonate production rate (", "kg CaCO"[3], " m"^-2, " yr"^-1, ")"))) +
      theme_bw() +
      theme(
        axis.text        = element_text(size = 14, color = "black"),
        axis.title       = element_text(size = 16, color = "black"),
        panel.grid.minor = element_blank(),
        legend.position  = "none"
      )
    
    scatter_path <- file.path(plot_out_dir, sprintf("%s_%s_Fixed_netprod_vs_coralcover_cnet.png", region, year))
    ggsave(cover_scatter, filename = scatter_path, bg = "white", height = 7, width = 9, dpi = 300)
    cat(green(sprintf("✅ Saved Net Prod vs. Coral Cover plot to: %s\n", scatter_path)))
  }
  
  return(ha_island_plot)
}

# ==============================================================================
# 7. EXECUTE MASTER PIPELINE
# ==============================================================================

# A. LOAD FIXED FISH BIOEROSION
fish_clean <- get_fixed_fish_data(fish_site_file, fish_file, meta_file)

# B. LOAD FIXED URCHIN BIOEROSION
urch_clean <- get_fixed_urchin_data(urch_site_file, urchin_file, meta_file)

# C. LOAD CORAL COVER (TIER 1 HARD CORAL ONLY)
cover_clean <- get_coral_cover_data(cover_file)

# D. CALCULATE NET BUDGETS FOR FIXED SITES
final_net_site <- run_fixed_net_pipeline(
  benthic_file  = benthic_file, 
  meta_file     = meta_file, 
  fish_clean    = fish_clean, 
  urch_clean    = urch_clean, 
  cover_clean   = cover_clean, 
  out_net_file  = out_net_file,
  exclude_sites = sites_to_exclude
)

# E. DERIVE SITE ORDER BY DECREASING LATITUDE (North to South, preserving uncoordinated sites)
ordered_with_lat <- final_net_site %>%
  filter(!is.na(LATITUDE)) %>%
  mutate(REGIONCODE = factor(REGIONCODE, levels = c("NWHI", "MHI"))) %>%
  arrange(REGIONCODE, desc(LATITUDE)) %>%
  distinct(OCC_SITEID) %>%
  pull(OCC_SITEID)

remaining_sites     <- setdiff(unique(final_net_site$OCC_SITEID), ordered_with_lat)
latitude_site_order <- c(ordered_with_lat, remaining_sites)

# F. GENERATE & SAVE OVERVIEW VISUALIZATIONS
overview_plots <- create_benthic_plots(
  data          = final_net_site, 
  site_order    = latitude_site_order, 
  exclude_sites = sites_to_exclude
)

ggsave(file.path(plot_out_dir, sprintf("%s_%s_gross_overview.png", region, year)), 
       overview_plots$gross_plot, width = 14, height = 6, bg = "white", dpi = 300)
ggsave(file.path(plot_out_dir, sprintf("%s_%s_net_overview.png", region, year)),   
       overview_plots$net_plot,   width = 14, height = 6, bg = "white", dpi = 300)

# G. GENERATE & SAVE COMPOSITE PLOT (NET GRADIENT + BUDGET COMPONENTS)
final_composite_plot <- create_budget_composite_plot(
  net_site      = final_net_site, 
  region        = region, 
  year          = year, 
  plot_out_dir  = plot_out_dir, 
  site_order    = latitude_site_order, 
  exclude_sites = sites_to_exclude
)

# H. GENERATE & SAVE ISLAND-LEVEL JITTER & CORAL COVER FIGURES
island_plot <- create_visualizations(
  comb_all = final_net_site, 
  plot_dir = plot_dir, 
  region   = region, 
  year     = year
)

cat(green(sprintf("\n🎉 MASTER PIPELINE COMPLETE! Site-level budgets saved as: %s\n", basename(out_net_file))))