# ==============================================================================
# Script:  3a. DarwinPoint_processCB_strs_cnet.R
# Purpose: Hawaiian Archipelago - StRS Net Carbonate Budget Calculation & Figures
#          Integrates site-level gross carbonate production from StRS surveys
#          (with fixed OCC stations removed) with endolithic bioerosion, StRS parrotfish 
#          bioerosion (nSPC), island-level urchin bioerosion fallbacks, and 
#          site-level live coral cover from Tier 1 photo-quadrats.
# Inputs:  - data/strs data/field data/HA_2024_StRS_Prod_SiteLevel.csv
#          - data/strs data/field data/HA_2024_StRS_Fish.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Urchins_SiteLevel.csv
#            (or data/fixed data/field data/HA_2024_Fixed_Urchins.csv)
#          - data/strs data/field data/2024_HA_StRS_Cover_Tier1.csv
# Outputs: - data/strs data/processed data/HA_2024_StRS_Fish_SiteLevel.csv
#          - data/strs data/processed data/HA_2024_StRS_Net_SiteLevel_cnet.csv
#          - outputs/supporting figures/HA_2024_gross_overview_strs.png
#          - outputs/supporting figures/HA_2024_net_overview_strs.png
#          - outputs/supporting figures/HA_2024_netprod_strs.png
#          - outputs/supporting figures/HA_2024_StRS_netprod.png
#          - outputs/supporting figures/HA_2024_StRS_netprod_vs_coralcover.png
# ==============================================================================

library(tidyverse)
library(data.table)
library(reefbudgetR)
library(crayon)
library(ggpubr)
library(here)

# ==============================================================================
# 1. PARAMETERS & PATHS
# ==============================================================================

region <- "HA"
year   <- "2024"

base_dir     <- here()
data_dir     <- file.path(base_dir, "data")
plot_out_dir <- file.path(base_dir, "outputs", "supporting figures")
proc_out_dir <- file.path(data_dir, "strs data", "processed data")

# Master input files
benthic_file <- file.path(data_dir, "strs data", "field data", "HA_2024_StRS_Prod_SiteLevel.csv")
fish_file    <- file.path(data_dir, "strs data", "field data", "HA_2024_StRS_Fish.csv")

# Dynamic resolution of coral cover file path
cover_candidates <- c(
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_Cover_Tier1.csv"),
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_CoralCover_Tier1.csv"),
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_CoralCover.csv")
)
cover_file <- cover_candidates[file.exists(cover_candidates)][1]
if (is.na(cover_file)) cover_file <- cover_candidates[1]

# Fixed urchin inputs (checks processed first, then field raw)
fixed_proc_dir  <- file.path(data_dir, "fixed data", "processed data")
urch_site_file  <- file.path(fixed_proc_dir, sprintf("%s_%s_Fixed_Urchins_SiteLevel.csv", region, year))
urch_field_file <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_Urchins.csv")

# Processed output paths
fish_site_file <- file.path(proc_out_dir, sprintf("%s_%s_StRS_Fish_SiteLevel.csv", region, year))
out_net_file   <- file.path(proc_out_dir, sprintf("%s_%s_StRS_Net_SiteLevel_cnet.csv", region, year))

stopifnot("Benthic StRS production file not found!" = file.exists(benthic_file))
stopifnot("StRS Fish field file not found!"         = file.exists(fish_file))

invisible(lapply(c(data_dir, proc_out_dir, plot_out_dir), function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}))

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

calc_dynamic_lims <- function(vec_min, vec_max, pad = 0.12, force_zero_lower = FALSE) {
  min_val <- min(vec_min, na.rm = TRUE)
  max_val <- max(vec_max, na.rm = TRUE)
  
  low <- if (force_zero_lower || min_val >= 0) 0 else min_val - abs(min_val) * pad
  high <- if (max_val > 0) max_val + abs(max_val) * pad else 1
  
  c(floor(low), ceiling(high))
}

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
# 2. StRS FISH BIOEROSION PIPELINE (nSPC)
# ==============================================================================

get_strs_fish_data <- function(fish_site_file, fish_file) {
  cat(cyan(sprintf("\n🚀 LOADING StRS FISH BIOEROSION: %s %s\n", region, year)))
  
  if (file.exists(fish_site_file)) {
    cat(blue(sprintf("📂 Reading existing processed fish data from: %s\n", basename(fish_site_file))))
    fish_df <- fread(fish_site_file, na.strings = c("", "NA"), data.table = FALSE)
  } else {
    cat(yellow(sprintf("⚠️ Processing raw StRS fish observations from: %s\n", basename(fish_file))))
    spc_raw <- fread(fish_file, na.strings = c("", "NA"), data.table = FALSE)
    
    spc_data <- spc_raw %>%
      mutate(
        REA_SITEID   = SITE,
        LOCATION     = ISLAND,
        LOCATIONCODE = if ("LOCATIONCODE" %in% names(.)) LOCATIONCODE else ISLAND,
        REGION       = REGION_NAME,
        REGIONCODE   = if_else(grepl("Northwestern", REGION_NAME, ignore.case = TRUE), "NWHI", "MHI"),
        CRUISE_ID    = MISSIONID,
        YEAR         = as.character(OBS_YEAR),
        METHOD       = if ("METHOD" %in% names(.)) METHOD else "nSPC"
      )
    
    fish_raw <- process_fish(
      spc_data    = spc_data, 
      method_type = "StRS", 
      dbase_type  = "Kindinger"
    )
    
    fish_df <- if (is.list(fish_raw) && "calc_strs_ero" %in% names(fish_raw)) {
      fish_raw$calc_strs_ero
    } else {
      as.data.frame(fish_raw)
    }
    
    fish_df$YEAR <- as.character(year)
    fwrite(fish_df, fish_site_file, na = "")
  }
  
  site_col <- intersect(c("SITE", "REA_SITEID", "OCC_SITEID"), names(fish_df))[1]
  
  eros_mean_col <- intersect(c("FISH_EROSION_KG_M2_YR_ALL_MEAN", "FISH_EROSION_KG_M2_YR_ALL"), names(fish_df))[1]
  eros_se_col   <- intersect(c("FISH_EROSION_KG_M2_YR_ALL_SE", "FISH_EROSION_KG_M2_YR_ALL_SD"), names(fish_df))[1]
  bio_mean_col  <- intersect(c("FISH_BIOMASS_KG_HA_ALL_MEAN", "FISH_BIOMASS_KG_HA_ALL"), names(fish_df))[1]
  bio_se_col    <- intersect(c("FISH_BIOMASS_KG_HA_ALL_SE", "FISH_BIOMASS_KG_HA_ALL_SD"), names(fish_df))[1]
  den_mean_col  <- intersect(c("FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN", "FISH_DENSITY_ABUNDANCE_HA_ALL"), names(fish_df))[1]
  den_se_col    <- intersect(c("FISH_DENSITY_ABUNDANCE_HA_ALL_SE", "FISH_DENSITY_ABUNDANCE_HA_ALL_SD"), names(fish_df))[1]
  
  fish_clean <- fish_df %>%
    transmute(
      OCC_SITEID                         = pad_strs_site(.data[[site_col]]),
      FISH_EROSION_KG_M2_YR_ALL_MEAN     = if (!is.na(eros_mean_col)) as.numeric(.data[[eros_mean_col]]) else 0,
      FISH_EROSION_KG_M2_YR_ALL_SE       = if (!is.na(eros_se_col))   as.numeric(.data[[eros_se_col]])   else 0,
      FISH_BIOMASS_KG_HA_ALL_MEAN        = if (!is.na(bio_mean_col))  as.numeric(.data[[bio_mean_col]])  else 0,
      FISH_BIOMASS_KG_HA_ALL_SE          = if (!is.na(bio_se_col))    as.numeric(.data[[bio_se_col]])    else 0,
      FISH_DENSITY_ABUNDANCE_HA_ALL_MEAN = if (!is.na(den_mean_col)) as.numeric(.data[[den_mean_col]]) else 0,
      FISH_DENSITY_ABUNDANCE_HA_ALL_SE   = if (!is.na(den_se_col))   as.numeric(.data[[den_se_col]])   else 0
    ) %>%
    mutate(
      FISH_EROSION_KG_M2_YR_ALL_SE = if_else(FISH_EROSION_KG_M2_YR_ALL_MEAN == 0, 0, FISH_EROSION_KG_M2_YR_ALL_SE)
    ) %>%
    filter(!grepl("^OCC", OCC_SITEID)) %>%
    distinct(OCC_SITEID, .keep_all = TRUE)
  
  cat(green(sprintf("✅ Loaded StRS parrotfish bioerosion for %d unique sites.\n", nrow(fish_clean))))
  return(fish_clean)
}

# ==============================================================================
# 3. FIXED URCHIN BIOEROSION PIPELINE
# ==============================================================================

get_fixed_urchin_data <- function(urch_site_file, urch_field_file) {
  cat(cyan(sprintf("\n🚀 LOADING FIXED URCHIN BIOEROSION FOR ISLAND FALLBACKS: %s %s\n", region, year)))
  
  if (file.exists(urch_site_file)) {
    cat(blue(sprintf("📂 Reading existing processed urchin data from: %s\n", basename(urch_site_file))))
    urch_df <- fread(urch_site_file, na.strings = c("", "NA"), data.table = FALSE)
  } else if (file.exists(urch_field_file)) {
    cat(yellow(sprintf("⚠️ Processing raw urchin field data from: %s\n", basename(urch_field_file))))
    urch_data    <- fread(urch_field_file, na.strings = c("", "NA"), data.table = FALSE)
    urch_results <- process_urchins(data = urch_data, full_summary = TRUE)
    urch_df      <- as.data.frame(urch_results$site_erosion)
  } else {
    warning("⚠️ No urchin file found; defaulting urchin erosion to 0.")
    urch_df <- data.frame(OCC_SITEID = character(0), LOCATION = character(0), 
                          URCHIN_EROSION_KG_M2_YR_MEAN = numeric(0), URCHIN_EROSION_KG_M2_YR_SE = numeric(0))
  }
  
  site_col   <- intersect(c("SITE", "OCC_SITEID", "REA_SITEID"), names(urch_df))[1]
  u_eros_col <- intersect(c("URCHIN_EROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR"), names(urch_df))[1]
  u_den_col  <- intersect(c("URCHIN_DENSITY_NO_M2_MEAN", "URCHIN_DENSITY_NO_M2"), names(urch_df))[1]
  
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
  
  cat(green(sprintf("✅ Loaded urchin bioerosion across %d stations.\n", nrow(urch_clean))))
  return(urch_clean)
}

# ==============================================================================
# 4. CORAL COVER INGESTION & HARMONIZATION (PT / PCT FILTERED FOR TIER_1 == 'CORAL')
# ==============================================================================

get_coral_cover_data <- function(cover_file) {
  if (!file.exists(cover_file)) {
    warning(sprintf("⚠️ Coral cover file not found at: %s", cover_file))
    return(data.frame(OCC_SITEID = character(0), CORAL_COVER_PERC = numeric(0)))
  }
  
  cat(blue(sprintf("📂 Ingesting coral cover from: %s\n", basename(cover_file))))
  cover_raw <- fread(cover_file, na.strings = c("", "NA"), data.table = FALSE)
  
  site_col <- intersect(c("SITE", "REA_SITEID", "OCC_SITEID"), names(cover_raw))[1]
  pct_col  <- intersect(c("PCT", "PT", "CORAL_COVER_PERC", "HARD_CORAL_COVER_PCT", "COVER"), names(cover_raw))[1]
  if (is.na(pct_col)) {
    pct_col <- grep("PCT|PT|COVER|PERC", names(cover_raw), value = TRUE, ignore.case = TRUE)[1]
  }
  tier_col <- intersect(c("TIER_1", "TIER1", "CATEGORY", "BENTHIC_CATEGORY", "TAXA"), names(cover_raw))[1]
  
  stopifnot("No valid SITE column in coral cover file!" = !is.na(site_col))
  stopifnot("No valid PCT/PT column in coral cover file!" = !is.na(pct_col))
  
  # Filter strictly for hard coral tier to avoid summing or averaging across all benthic classes
  if (!is.na(tier_col)) {
    cover_raw <- cover_raw %>%
      filter(grepl("^CORAL$", trimws(toupper(.data[[tier_col]]))))
  }
  
  cover_clean <- cover_raw %>%
    transmute(
      OCC_SITEID       = pad_strs_site(.data[[site_col]]),
      CORAL_COVER_PERC = as.numeric(.data[[pct_col]])
    ) %>%
    filter(!is.na(OCC_SITEID) & !is.na(CORAL_COVER_PERC) & !grepl("^OCC", OCC_SITEID)) %>%
    group_by(OCC_SITEID) %>%
    summarize(CORAL_COVER_PERC = mean(CORAL_COVER_PERC, na.rm = TRUE), .groups = "drop")
  
  cat(green(sprintf("✅ Successfully extracted live coral cover across %d StRS survey stations (Mean: %.2f%%).\n",
                    nrow(cover_clean), mean(cover_clean$CORAL_COVER_PERC, na.rm = TRUE))))
  return(cover_clean)
}

# ==============================================================================
# 5. NET BUDGET PIPELINE (PURE StRS SITES ONLY - OCC EXCLUDED)
# ==============================================================================

run_strs_net_pipeline <- function(benthic_file, fish_clean, urch_clean, cover_clean, out_net_file) {
  cat(cyan(sprintf("\n🚀 PAIRING StRS PRODUCTION (EXCLUDING OCC) WITH BIOEROSION via process_net(): %s %s\n", region, year)))
  
  strs_data <- fread(benthic_file, na.strings = c("", "NA"), data.table = FALSE)
  
  if ("DATE_" %in% names(strs_data) && !"LOCALDATE" %in% names(strs_data)) {
    strs_data$LOCALDATE <- strs_data$DATE_
  }
  
  meta_cols_order <- c("REGION", "REGIONCODE", "YEAR", "CRUISE_ID", "LOCATION", "LOCATIONCODE",
                       "OCC_SITEID", "SITEVISITID", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "LOCALDATE")
  
  matched_strs <- strs_data %>%
    filter(!grepl("^OCC", trimws(SITE))) %>%
    mutate(
      OCC_SITEID                    = pad_strs_site(SITE),
      LOCATION                      = if ("LOCATION" %in% names(.) && "ISLAND" %in% names(.)) {
        coalesce(LOCATION, ISLAND)
      } else if ("ISLAND" %in% names(.)) {
        ISLAND
      } else {
        NA_character_
      },
      GROSS_CARB_PROD_KG_M2_YR_MEAN = as.numeric(PROD),
      GROSS_CARB_PROD_KG_M2_YR_SE   = if (all(c("PROD_upr", "PROD_lwr") %in% names(.))) {
        coalesce(abs(as.numeric(PROD_upr) - as.numeric(PROD_lwr)) / (2 * 1.96), 0)
      } else { 0 },
      GROSS_CARB_PROD_KG_M2_YR_N    = 1, 
      # Bioerosion Mean & SE calculated from MICROMACRO_TOTAL, MICROMACRO_upr, and MICROMACRO_lwr
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
    left_join(cover_clean, by = "OCC_SITEID") %>%
    mutate(
      HARD_CORAL_COVER_PCT_MEAN = coalesce(CORAL_COVER_PERC, 0),
      HARD_CORAL_COVER_PCT_SE   = 0
    )
  
  if (!"REGIONCODE" %in% names(matched_strs) && "REGION" %in% names(matched_strs)) {
    matched_strs$REGIONCODE <- matched_strs$REGION
  }
  
  cat(blue(sprintf("🔍 Isolated %d pure StRS randomized survey stations (excluded OCC sites).\n", nrow(matched_strs))))
  
  for (col in meta_cols_order) {
    if (!col %in% names(matched_strs)) matched_strs[[col]] <- NA
  }
  matched_strs <- matched_strs %>% relocate(all_of(meta_cols_order))
  
  # Island fallback summary for urchin bioerosion
  urch_isl <- urch_clean %>%
    filter(!is.na(LOCATION) & !is.na(URCHIN_EROSION_KG_M2_YR_MEAN)) %>%
    group_by(LOCATION) %>%
    summarize(
      URCHIN_ISL_MEAN = mean(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE),
      URCHIN_ISL_SE   = if (n() > 1) sd(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE) / sqrt(n()) else 0,
      .groups = "drop"
    )
  
  urch_formatted <- matched_strs %>%
    select(OCC_SITEID, LOCATION_meta = any_of("LOCATION")) %>%
    distinct() %>%
    left_join(urch_clean, by = "OCC_SITEID") %>%
    mutate(LOCATION = if ("LOCATION_meta" %in% names(.)) coalesce(LOCATION, LOCATION_meta) else LOCATION) %>%
    left_join(urch_isl, by = "LOCATION") %>%
    mutate(
      URCHIN_EROSION_KG_M2_YR_MEAN = coalesce(URCHIN_EROSION_KG_M2_YR_MEAN, URCHIN_ISL_MEAN, 0),
      URCHIN_EROSION_KG_M2_YR_SE   = coalesce(URCHIN_EROSION_KG_M2_YR_SE, URCHIN_ISL_SE, 0),
      URCHIN_DENSITY_NO_M2_MEAN    = coalesce(URCHIN_DENSITY_NO_M2_MEAN, 0),
      URCHIN_DENSITY_NO_M2_SE      = coalesce(URCHIN_DENSITY_NO_M2_SE, 0),
      YEAR                         = as.character(year)
    ) %>%
    select(-any_of(c("LOCATION_meta", "URCHIN_ISL_MEAN", "URCHIN_ISL_SE")))
  
  for (col in meta_cols_order) {
    if (!col %in% names(urch_formatted)) urch_formatted[[col]] <- NA
  }
  urch_formatted <- urch_formatted %>% relocate(all_of(meta_cols_order))
  
  # Retain complete survey metadata for fish observations
  fish_formatted <- matched_strs %>%
    select(all_of(meta_cols_order)) %>%
    left_join(fish_clean, by = "OCC_SITEID") %>%
    mutate(
      across(starts_with("FISH_"), ~ coalesce(.x, 0)),
      YEAR = as.character(year)
    )
  
  cat(blue("🔗 Joining Fish, Urchins, and Processing Net Budgets via process_net()...\n"))
  
  comb_all <- process_net(
    prod   = matched_strs,
    urch   = urch_formatted,
    fish   = fish_formatted,
    sum_by = "site",
    format = "wide"
  )
  
  if ("OCC_SITEID" %in% names(comb_all) && !"SITE" %in% names(comb_all)) {
    comb_all <- comb_all %>% rename(SITE = OCC_SITEID)
  } else if ("OCC_SITEID" %in% names(comb_all) && "SITE" %in% names(comb_all)) {
    comb_all$SITE <- coalesce(comb_all$SITE, comb_all$OCC_SITEID)
  }
  
  # Ensure geographic coordinates are preserved post-process_net
  if (!"LATITUDE" %in% names(comb_all) || all(is.na(comb_all$LATITUDE))) {
    comb_all <- comb_all %>%
      select(-any_of(c("LATITUDE", "LONGITUDE"))) %>%
      left_join(
        matched_strs %>% select(SITE = OCC_SITEID, LATITUDE, LONGITUDE) %>% distinct(),
        by = "SITE"
      )
  }
  
  # Preserve explicit internal bioerosion mean and SE columns
  if (!"BIOEROSION_KG_M2_YR_SE" %in% names(comb_all)) {
    comb_all <- comb_all %>%
      select(-any_of(c("BIOEROSION_KG_M2_YR_MEAN", "BIOEROSION_KG_M2_YR_SE"))) %>%
      left_join(
        matched_strs %>% select(SITE = OCC_SITEID, BIOEROSION_KG_M2_YR_MEAN, BIOEROSION_KG_M2_YR_SE) %>% distinct(),
        by = "SITE"
      )
  }
  
  # Preserve live coral cover percentage if dropped during process_net
  if (!"HARD_CORAL_COVER_PCT_MEAN" %in% names(comb_all) || all(is.na(comb_all$HARD_CORAL_COVER_PCT_MEAN))) {
    comb_all <- comb_all %>%
      select(-any_of(c("HARD_CORAL_COVER_PCT_MEAN", "HARD_CORAL_COVER_PCT_SE"))) %>%
      left_join(
        matched_strs %>% select(SITE = OCC_SITEID, HARD_CORAL_COVER_PCT_MEAN, HARD_CORAL_COVER_PCT_SE) %>% distinct(),
        by = "SITE"
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
  
  fwrite(comb_all, out_net_file, na = "")
  
  cat(green(sprintf("✅ Saved pure StRS net carbonate budgets (%d sites) to:\n  • %s\n", 
                    nrow(comb_all), out_net_file)))
  return(comb_all)
}

# ==============================================================================
# 6. VISUALIZATION FUNCTIONS (SITE-LEVEL)
# ==============================================================================

create_benthic_plots <- function(data, site_order = NULL) {
  if (!is.null(site_order)) {
    valid_sites <- site_order[site_order %in% unique(data$SITE)]
    data$SITE <- factor(data$SITE, levels = valid_sites)
  } else if ("LATITUDE" %in% names(data)) {
    lat_order <- data %>%
      filter(!is.na(LATITUDE)) %>%
      distinct(SITE, LATITUDE) %>%
      arrange(desc(LATITUDE)) %>%
      pull(SITE)
    other_sites <- setdiff(unique(data$SITE), lat_order)
    data$SITE <- factor(data$SITE, levels = c(lat_order, other_sites))
  }
  
  lims_gross <- calc_dynamic_lims(
    data$GROSS_CARB_PROD_KG_M2_YR_MEAN - coalesce(data$GROSS_CARB_PROD_KG_M2_YR_SE, 0),
    data$GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(data$GROSS_CARB_PROD_KG_M2_YR_SE, 0)
  )
  lims_net <- calc_dynamic_lims(
    data$NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(data$NET_CARB_PROD_KG_M2_YR_SE, 0),
    data$NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(data$NET_CARB_PROD_KG_M2_YR_SE, 0)
  )
  
  pal <- c("NWHI" = "#003087", "MHI" = "#0085ca")
  
  theme_benthic <- theme_bw() + 
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 4, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y      = element_text(size = 12, color = 'black'),
      axis.title.y     = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      legend.position  = 'none' 
    )
  
  gross_plot <- ggplot(data) +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = SITE, y = GROSS_CARB_PROD_KG_M2_YR_MEAN, fill = REGIONCODE),
             color = 'black', width = 0.6) +
    geom_errorbar(aes(x = SITE, 
                      ymin = GROSS_CARB_PROD_KG_M2_YR_MEAN - coalesce(GROSS_CARB_PROD_KG_M2_YR_SE, 0), 
                      ymax = GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(GROSS_CARB_PROD_KG_M2_YR_SE, 0)), 
                  width = 0.2) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = lims_gross, expand = c(0, 0)) +
    scale_fill_manual(values = pal) +
    ylab(expression(paste("Gross production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ") "))) +
    theme_benthic
  
  net_plot <- ggplot(data) +
    geom_hline(yintercept = 0) +
    geom_col(aes(x = SITE, y = NET_CARB_PROD_KG_M2_YR_MEAN, fill = REGIONCODE),
             color = 'black', width = 0.6) +
    geom_errorbar(aes(x = SITE, 
                      ymin = NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(NET_CARB_PROD_KG_M2_YR_SE, 0), 
                      ymax = NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(NET_CARB_PROD_KG_M2_YR_SE, 0)), 
                  width = 0.2) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = lims_net, expand = c(0, 0)) +
    scale_fill_manual(values = pal) +
    ylab(expression(paste("Net production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ") "))) +
    theme_benthic
  
  return(list(gross_plot = gross_plot, net_plot = net_plot))
}

format_bygroup <- function(net_data, prod_group) {
  req_eros <- c("BIOEROSION_KG_M2_YR_MEAN", "MACROBIOEROSION_KG_M2_YR_MEAN", 
                "MICROBIOEROSION_KG_M2_YR_MEAN", "URCHIN_EROSION_KG_M2_YR_MEAN", 
                "FISH_EROSION_KG_M2_YR_ALL_MEAN")
  for (col in req_eros) {
    if (!col %in% names(net_data)) net_data[[col]] <- 0
  }
  
  has_combined_bioeros <- any(net_data$BIOEROSION_KG_M2_YR_MEAN > 0, na.rm = TRUE)
  
  if (has_combined_bioeros) {
    eros_clean <- net_data %>%
      select(any_of(c("REGION", "REGIONCODE", "LOCATION", "SITE", "LATITUDE")), 
             `Internal Bioerosion` = BIOEROSION_KG_M2_YR_MEAN,
             Urchins               = URCHIN_EROSION_KG_M2_YR_MEAN,
             Parrotfish            = FISH_EROSION_KG_M2_YR_ALL_MEAN) %>%
      pivot_longer(cols = c(`Internal Bioerosion`, Urchins, Parrotfish), 
                   names_to = "GROUP", 
                   values_to = "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN") %>%
      mutate(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = -1 * SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN)
  } else {
    eros_clean <- net_data %>%
      select(any_of(c("REGION", "REGIONCODE", "LOCATION", "SITE", "LATITUDE")), 
             Macroborers = MACROBIOEROSION_KG_M2_YR_MEAN,
             Microborers = MICROBIOEROSION_KG_M2_YR_MEAN,
             Urchins     = URCHIN_EROSION_KG_M2_YR_MEAN,
             Parrotfish  = FISH_EROSION_KG_M2_YR_ALL_MEAN) %>%
      pivot_longer(cols = c(Macroborers, Microborers, Urchins, Parrotfish), 
                   names_to = "GROUP", 
                   values_to = "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN") %>%
      mutate(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = -1 * SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN)
  }
  
  bind_rows(prod_group, eros_clean)
}

create_budget_composite_plot <- function(net_site, region, year, plot_out_dir, site_order = NULL) {
  cat(blue("\n🎨 Generating StRS Composite Net Budget Figure...\n"))
  
  prod_group <- net_site %>%
    select(any_of(c("REGION", "REGIONCODE", "LOCATION", "SITE", "LATITUDE")), 
           SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = GROSS_CARB_PROD_KG_M2_YR_MEAN) %>%
    mutate(GROUP = "Gross Production")
  
  data_group <- format_bygroup(net_site, prod_group)
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
    valid_order <- site_order[site_order %in% unique(data_group$SITE)]
    data_group$SITE <- factor(data_group$SITE, levels = valid_order)
    net_site$SITE   <- factor(net_site$SITE, levels = valid_order)
  } else if ("LATITUDE" %in% names(net_site)) {
    lat_order <- net_site %>%
      filter(!is.na(LATITUDE)) %>%
      distinct(SITE, LATITUDE) %>%
      arrange(desc(LATITUDE)) %>%
      pull(SITE)
    other_sites     <- setdiff(unique(net_site$SITE), lat_order)
    valid_order     <- c(lat_order, other_sites)
    data_group$SITE <- factor(data_group$SITE, levels = valid_order)
    net_site$SITE   <- factor(net_site$SITE, levels = valid_order)
  }
  
  max_prod_bar <- max(net_site$GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(net_site$GROSS_CARB_PROD_KG_M2_YR_SE, 0), na.rm = TRUE)
  min_eros_bar <- min(-1 * (net_site$GROSS_EROSION_KG_M2_YR_MEAN + coalesce(net_site$GROSS_EROSION_KG_M2_YR_SE, 0)), na.rm = TRUE)
  min_net_val  <- min(net_site$NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(net_site$NET_CARB_PROD_KG_M2_YR_SE, 0), na.rm = TRUE)
  max_net_val  <- max(net_site$NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(net_site$NET_CARB_PROD_KG_M2_YR_SE, 0), na.rm = TRUE)
  
  overall_min   <- min(min_eros_bar, min_net_val, 0, na.rm = TRUE)
  overall_max   <- max(max_prod_bar, max_net_val, 0, na.rm = TRUE)
  budget_y_lims <- c(floor(overall_min - abs(overall_min) * 0.12), ceiling(overall_max + abs(overall_max) * 0.12))
  
  theme_budget <- theme_bw() + 
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 4, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y      = element_text(size = 14, color = 'black'),
      axis.title.y     = element_text(size = 16),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      strip.text       = element_blank(),
      legend.position  = "none" 
    )
  
  prod_bycomponent_plot <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_col(data = data_group, 
             aes(x = SITE, y = SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, fill = GROUP), width = 0.5) +
    geom_col(data = net_site, 
             aes(x = SITE, y = GROSS_CARB_PROD_KG_M2_YR_MEAN), 
             alpha = 0, color = "black", width = 0.5) +
    geom_col(data = net_site, 
             aes(x = SITE, y = -1 * GROSS_EROSION_KG_M2_YR_MEAN), 
             alpha = 0, color = "black", width = 0.5) +
    geom_errorbar(data = net_site, 
                  aes(x = SITE, 
                      ymin = GROSS_CARB_PROD_KG_M2_YR_MEAN, 
                      ymax = GROSS_CARB_PROD_KG_M2_YR_MEAN + coalesce(GROSS_CARB_PROD_KG_M2_YR_SE, 0)), 
                  width = 0.25) +
    geom_errorbar(data = net_site, 
                  aes(x = SITE, 
                      ymin = -1 * GROSS_EROSION_KG_M2_YR_MEAN - coalesce(GROSS_EROSION_KG_M2_YR_SE, 0), 
                      ymax = -1 * GROSS_EROSION_KG_M2_YR_MEAN), 
                  width = 0.25) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = budget_y_lims, expand = c(0, 0)) +
    scale_fill_manual(values = pal_budget) +
    ylab(expression(paste("Production / bioerosion (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ") "))) +
    theme_budget
  
  net_plot_gradient <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_col(data = net_site,
             aes(x = SITE, y = NET_CARB_PROD_KG_M2_YR_MEAN, fill = NET_CARB_PROD_KG_M2_YR_MEAN),
             color = 'black', width = 0.5) +
    geom_errorbar(data = net_site, 
                  aes(x = SITE, 
                      ymin = NET_CARB_PROD_KG_M2_YR_MEAN - coalesce(NET_CARB_PROD_KG_M2_YR_SE, 0), 
                      ymax = NET_CARB_PROD_KG_M2_YR_MEAN + coalesce(NET_CARB_PROD_KG_M2_YR_SE, 0)), 
                  width = 0.25) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = budget_y_lims, expand = c(0, 0)) +
    scale_fill_gradient2(low = "#EF2820", mid = "#CCCCCC", high = "#38A1D0", midpoint = 0, limits = budget_y_lims) +
    ylab(expression(paste("Net carbonate production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ") "))) +
    theme_budget
  
  legend_components <- get_legend(
    ggplot(data_group, aes(x = SITE, y = SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, fill = GROUP)) +
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
  
  out_path <- file.path(plot_out_dir, sprintf("%s_%s_netprod_strs.png", region, year))
  ggsave(final_plot, filename = out_path, bg = "white", height = 12, width = 24, dpi = 300)
  
  cat(green(sprintf("✅ Composite plot saved to: %s\n", out_path)))
  return(final_plot)
}

create_visualizations <- function(comb_all, plot_dir, region, year) {
  cat(cyan("\n🎨 GENERATING StRS LOCATION-LEVEL AND CORAL COVER FIGURES...\n"))
  
  location_levels <- c("Hawaii", "Kahoolawe", "Lanai", "Maui", "Molokai", "Oahu", 
                       "Kauai", "Niihau", "French Frigate", "Lisianski", 
                       "Pearl & Hermes", "Kure")
  
  comb_all <- format_okinas(comb_all, location_levels)
  
  y_min <- min(comb_all$NET_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE)
  y_max <- max(comb_all$NET_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE)
  y_limits <- c(min(-8, floor(y_min - 1)), max(20, ceiling(y_max + 1)))
  
  has_cover <- "HARD_CORAL_COVER_PCT_MEAN" %in% names(comb_all) && 
    any(!is.na(comb_all$HARD_CORAL_COVER_PCT_MEAN) & comb_all$HARD_CORAL_COVER_PCT_MEAN > 0)
  
  ha_island_plot <- ggplot(data = comb_all) +
    {
      if (has_cover) {
        geom_jitter(aes(x = labels, y = NET_CARB_PROD_KG_M2_YR_MEAN, 
                        color = NET_CARB_PROD_KG_M2_YR_MEAN, size = HARD_CORAL_COVER_PCT_MEAN), 
                    width = 0.22, height = 0, alpha = 0.75)
      } else {
        geom_jitter(aes(x = labels, y = NET_CARB_PROD_KG_M2_YR_MEAN, 
                        color = NET_CARB_PROD_KG_M2_YR_MEAN), 
                    width = 0.22, height = 0, size = 3, alpha = 0.75)
      }
    } +
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    scale_color_gradient2(low = "#EF2820", mid = "#CCCCCC", high = "#38A1D0", 
                          midpoint = 0, limits = y_limits, name = "Net Prod") +
    scale_y_continuous(limits = y_limits, breaks = scales::pretty_breaks(n = 8), expand = c(0, 0)) +
    ylab(expression(atop(paste("Net carbonate production rate"), paste("(kg CaCO"[3], " m"^-2, " yr"^-1, ") ")))) +
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
      axis.line         = element_line(colour = "black"),
      axis.ticks        = element_line(colour = "black")
    )
  
  if (has_cover) {
    ha_island_plot <- ha_island_plot + 
      scale_size_continuous(range = c(2, 7), name = "Coral Cover (%)") +
      guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))
  }
  
  island_path <- file.path(plot_dir, sprintf("%s_%s_StRS_netprod.png", region, year))
  ggsave(ha_island_plot, filename = island_path, bg = "white", height = 6, width = 15, dpi = 300)
  cat(green(sprintf("✅ Saved location overview plot to: %s\n", island_path)))
  
  if (has_cover) {
    cover_scatter <- ggplot(data = comb_all %>% filter(!is.na(HARD_CORAL_COVER_PCT_MEAN)), 
                            aes(x = HARD_CORAL_COVER_PCT_MEAN, y = NET_CARB_PROD_KG_M2_YR_MEAN)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_point(aes(color = NET_CARB_PROD_KG_M2_YR_MEAN), size = 3, alpha = 0.75) +
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
    
    scatter_path <- file.path(plot_dir, sprintf("%s_%s_StRS_netprod_vs_coralcover.png", region, year))
    ggsave(cover_scatter, filename = scatter_path, bg = "white", height = 7, width = 9, dpi = 300)
    cat(green(sprintf("✅ Saved Net Prod vs. Coral Cover plot to: %s\n", scatter_path)))
  }
  
  return(ha_island_plot)
}

# ==============================================================================
# 7. MASTER EXECUTION
# ==============================================================================

# A. LOAD StRS FISH BIOEROSION
fish_clean <- get_strs_fish_data(fish_site_file, fish_file)

# B. LOAD FIXED URCHIN BIOEROSION
urch_clean <- get_fixed_urchin_data(urch_site_file, urch_field_file)

# C. LOAD CORAL COVER (TIER 1 HARD CORAL ONLY)
cover_clean <- get_coral_cover_data(cover_file)

# D. CALCULATE NET BUDGETS FOR StRS SITES (OCC SITES EXCLUDED)
final_net_site <- run_strs_net_pipeline(
  benthic_file = benthic_file, 
  fish_clean   = fish_clean, 
  urch_clean   = urch_clean, 
  cover_clean  = cover_clean, 
  out_net_file = out_net_file
)

# E. DERIVE SITE ORDER BY DECREASING LATITUDE (North to South, preserving all sites)
region_order <- c("NWHI", "MHI")
final_net_site$REGIONCODE <- factor(final_net_site$REGIONCODE, levels = region_order)

ordered_with_lat <- final_net_site %>%
  filter(!is.na(LATITUDE)) %>%
  arrange(REGIONCODE, desc(LATITUDE)) %>%
  pull(SITE) %>%
  unique()

remaining_sites     <- setdiff(unique(final_net_site$SITE), ordered_with_lat)
latitude_site_order <- c(ordered_with_lat, remaining_sites)

# F. GENERATE & SAVE OVERVIEW PLOTS
overview_plots <- create_benthic_plots(
  data       = final_net_site, 
  site_order = latitude_site_order
)

ggsave(file.path(plot_out_dir, sprintf("%s_%s_gross_overview_strs.png", region, year)), 
       overview_plots$gross_plot, width = 14, height = 6, bg = "white", dpi = 300)
ggsave(file.path(plot_out_dir, sprintf("%s_%s_net_overview_strs.png", region, year)),   
       overview_plots$net_plot,   width = 14, height = 6, bg = "white", dpi = 300)

# G. GENERATE & SAVE COMPOSITE PLOT (NET GRADIENT + BUDGET COMPONENTS)
final_composite_plot <- create_budget_composite_plot(
  net_site      = final_net_site, 
  region        = region, 
  year          = year, 
  plot_out_dir  = plot_out_dir, 
  site_order    = latitude_site_order
)

# H. GENERATE & SAVE LOCATION-LEVEL JITTER & CORAL COVER FIGURES
island_plot <- create_visualizations(
  comb_all = final_net_site, 
  plot_dir = plot_out_dir, 
  region   = region, 
  year     = year
)

cat(green(sprintf("\n🎉 MASTER PIPELINE COMPLETE! %d randomized StRS stations (0 OCC) exported to: %s\n", 
                  nrow(final_net_site), basename(out_net_file))))