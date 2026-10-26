# ==============================================================================
# Script:  2. DarwinPoint_processCB_fixed_sfm.R
# Purpose: Hawaiian Archipelago - Net Carbonate Budget Calculation & Visualization
#          Integrates benthic gross production (Structure-from-Motion photogrammetry),
#          urchin bioerosion, and fish bioerosion to generate transect- and site-level
#          net carbonate budgets and figures.
# Inputs:  - data/fixed data/field data/HA_2024_Fixed_SiteSummary.csv
#          - data/fixed data/field data/HA_2024_Fixed_Benthic_sfm.csv (or HA_2024_Fixed_Benthic.csv)
#          - data/fixed data/field data/HA_2024_Fixed_Urchins.csv
#          - data/fixed data/field data/HA_2024_Fixed_Fish.csv
#          - data/dbase/prod_dbase_ha_fixed.csv
# Outputs: - data/fixed data/processed data/HA_2024_Fixed_Net_TransectLevel_sfm.csv
#          - data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_sfm.csv
#          - outputs/supporting figures/HA_2024_gross_overview_sfm.png
#          - outputs/supporting figures/HA_2024_net_overview_sfm.png
#          - outputs/supporting figures/HA_2024_netprod_fixed_sfm.png
# ==============================================================================

library(tidyverse)
library(lubridate)
library(data.table)
library(abjutils)
library(reefbudgetR)
library(crayon)
library(ggpubr) 
library(here)

# ==============================================================================
# 1. PARAMETERS & DATA INGESTION
# ==============================================================================

region <- "HA"
year   <- "2024"

# Base directory relative to project root
base_dir     <- here()
data_dir     <- file.path(base_dir, "data")
plot_out_dir <- file.path(base_dir, "outputs", "supporting figures")

# Master input files stored inside data/
meta_file       <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_SiteSummary.csv")
urchin_file     <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_Urchins.csv")
fish_file       <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_Fish.csv")
prod_dbase_file <- file.path(data_dir, "dbase", "prod_dbase_ha_fixed.csv")
benthic_file    <- file.path(data_dir, "fixed data", "field data", "HA_2024_Fixed_Benthic.csv")

sites_to_exclude <- c("OCC-MOL-003")

# Standard NCRMP Coral Group Name Lookup Table
coral_group_lookup <- c(
  "ACBR" = "Branching Acropora",      "ACTA" = "Table Acropora",
  "ASMA" = "Massive Astreopora",      "BR"   = "Other branching corals",
  "ENC"  = "Other encrusting corals",  "FOL"  = "Foliose corals",
  "FREE" = "Free-living corals",      "GOMA" = "Massive Goniastrea",
  "MASS" = "Other massive corals",    "MOBR" = "Other branching Montipora",
  "MOEN" = "Encrusting Montipora",    "MOFO" = "Foliose Montipora",
  "PAEN" = "Encrusting Pavona",       "PAMA" = "Massive Pavona",
  "POBR" = "Branching Porites",       "POCS" = "Branching Pocillopora",
  "POEN" = "Encrusting Porites",      "POLC" = "Laminar columnar Porites",
  "POMA" = "Massive Porites"
)

# Ensure output directories exist
invisible(lapply(c(data_dir, plot_out_dir, file.path(data_dir, "fixed data", "processed data")), function(d) {
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

# Repair metadata and group names in coral summaries
clean_prod_coral <- function(prod_coral, prod_site) {
  meta_cols <- setdiff(names(prod_site), c("OCC_SITEID", grep("CARB_PROD|COVER", names(prod_site), value = TRUE)))
  
  df <- prod_coral %>%
    mutate(CORAL_GROUP_NAME = coalesce(CORAL_GROUP_NAME, coral_group_lookup[CORAL_GROUP])) %>%
    left_join(prod_site %>% select(any_of(c("OCC_SITEID", meta_cols))), by = "OCC_SITEID", suffix = c("", ".meta"))
  
  for (col in meta_cols) {
    meta_col <- paste0(col, ".meta")
    if (col %in% names(df) && meta_col %in% names(df)) {
      df[[col]] <- coalesce(df[[col]], df[[meta_col]])
      df[[meta_col]] <- NULL
    }
  }
  
  df %>% select(-ends_with(".meta"))
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

# ==============================================================================
# 2. GROSS PRODUCTION PIPELINE (DIRECT BENTHIC CSV INPUT - SfM)
# ==============================================================================

run_benthic_pipeline <- function(data_dir, region, year, meta_file, prod_dbase_file, benthic_file) {
  cat(cyan(sprintf("\n🚀 STARTING BENTHIC PRODUCTION PIPELINE (SfM): %s %s\n", region, year)))
  cat(blue(sprintf("📂 Loading reference databases and benthic data from: %s\n", basename(benthic_file))))
  
  sites_metadata <- fread(meta_file, na.strings = c("", "NA"), data.table = FALSE)
  prod_dbase     <- fread(prod_dbase_file, na.strings = c("", "NA"), data.table = FALSE)
  benthic_data   <- fread(benthic_file, na.strings = c("", "NA"), data.table = FALSE)
  
  cat(blue("🧮 Calculating Gross Carbonate Production with Custom Lookup...\n"))
  
  # Strategy 1: Run whole-dataset at once (standard reefbudgetR design)
  results <- tryCatch(
    {
      process_prod(
        data              = as.data.frame(benthic_data), 
        dbase_type        = "Custom", 
        prod_dbase_custom = prod_dbase, 
        macro_rates       = "NCRMP", 
        micro_rates       = "NCRMP"
      )
    },
    error = function(e) {
      cat(yellow("⚠️ Full batch execution failed; falling back to protected per-site processing...\n"))
      NULL
    }
  )
  
  # Strategy 2: Protected per-site fallback handling zero-coral sites
  if (is.null(results)) {
    site_data_list <- split(as.data.frame(benthic_data), benthic_data$OCC_SITEID)
    
    results_list <- lapply(names(site_data_list), function(s_id) {
      df <- site_data_list[[s_id]]
      has_coral <- any(toupper(df$SUBSTRATE_CLASS) == "CORAL", na.rm = TRUE)
      
      if (!has_coral) {
        cat(yellow(sprintf("ℹ️ Site %s contains 0 hard coral observations.\n", s_id)))
      }
      
      tryCatch(
        process_prod(
          data              = df, 
          dbase_type        = "Custom", 
          prod_dbase_custom = prod_dbase, 
          macro_rates       = "NCRMP", 
          micro_rates       = "NCRMP"
        ),
        error = function(err) {
          cat(red(sprintf("⚠️ Skipped site %s due to processing error: %s\n", s_id, err$message)))
          NULL
        }
      )
    })
    
    results_list  <- purrr::compact(results_list)
    prod_site     <- rbindlist(lapply(results_list, `[[`, "summary_site"), fill = TRUE)
    prod_transect <- rbindlist(lapply(results_list, `[[`, "summary_transect"), fill = TRUE)
    prod_coral    <- rbindlist(lapply(results_list, `[[`, "summary_site_coral"), fill = TRUE)
  } else {
    prod_site     <- results$summary_site
    prod_transect <- results$summary_transect
    prod_coral    <- results$summary_site_coral
  }
  
  # Ensure site-level metadata is bound to site and transect outputs
  meta_to_join <- sites_metadata %>% 
    select(any_of(c("REGION", "REGIONCODE", "ISLAND", "LOCATION", "LOCATIONCODE", 
                    "YEAR", "LATITUDE", "LONGITUDE", "SITE_DEPTH_M", "OCC_SITEID"))) %>%
    distinct(OCC_SITEID, .keep_all = TRUE)
  
  prod_site     <- merge_metadata(as.data.frame(prod_site), meta_to_join)
  prod_transect <- merge_metadata(as.data.frame(prod_transect), meta_to_join)
  
  # Clean prod_coral to repair zero-coral site metadata and taxonomy labels
  prod_coral <- clean_prod_coral(as.data.frame(prod_coral), as.data.frame(prod_site))
  
  cat(green("✅ Benthic Production Pipeline Complete.\n"))
  return(list(
    prod_transect = as.data.frame(prod_transect), 
    prod_site     = as.data.frame(prod_site), 
    prod_coral    = as.data.frame(prod_coral)
  ))
}

# ==============================================================================
# 3. URCHIN PIPELINE
# ==============================================================================

run_urchin_pipeline <- function(data_dir, region, year, meta_file, urchin_file) {
  cat(cyan(sprintf("\n🚀 STARTING URCHIN BIOEROSION PIPELINE: %s %s\n", region, year)))
  
  sites_metadata <- fread(meta_file, data.table = FALSE)
  urch_data      <- fread(urchin_file, na.strings = c("", "NA"), data.table = FALSE)
  urch_data$LOCALDATE <- lubridate::ymd(urch_data$LOCALDATE)
  
  urch_data <- prep_urchins(data = urch_data, sites_metadata = sites_metadata) 
  
  cat(blue("🦔 Processing Urchin Bioerosion...\n"))
  urch_results  <- process_urchins(data = urch_data, full_summary = TRUE)
  urch_site     <- as.data.frame(urch_results$site_erosion)
  urch_transect <- as.data.frame(urch_results$transect_erosion)
  
  # Export processed urchin datasets
  out_proc_dir <- file.path(data_dir, "fixed data", "processed data")
  fwrite(urch_transect, file.path(out_proc_dir, sprintf("%s_%s_Fixed_Urchins_TransectLevel.csv", region, year)), na = "")
  fwrite(urch_site,     file.path(out_proc_dir, sprintf("%s_%s_Fixed_Urchins_SiteLevel.csv", region, year)), na = "")
  
  cat(green("✅ Urchin Bioerosion Pipeline Complete. Transect & Site urchin data exported.\n"))
  return(list(urch_site = urch_site, urch_transect = urch_transect))
}

# ==============================================================================
# 4. FISH PIPELINE
# ==============================================================================

run_fish_pipeline <- function(data_dir, region, year, meta_file, fish_file) {
  cat(cyan(sprintf("\n🚀 STARTING FISH BIOEROSION PIPELINE: %s %s\n", region, year)))
  
  sites_metadata <- fread(meta_file, data.table = FALSE)
  spc_data       <- fread(fish_file, na.strings = c("", "NA"), data.table = FALSE)
  
  cat(blue("🐟 Processing Parrotfish Bioerosion...\n"))
  fish_site <- process_fish(
    spc_data       = spc_data, 
    method_type    = "Fixed", 
    fixed_metadata = sites_metadata, 
    dbase_type     = "Kindinger"
  ) 
  fish_site$YEAR <- as.character(year) 
  
  # Export processed fish dataset
  out_proc_dir  <- file.path(data_dir, "fixed data", "processed data")
  fish_out_file <- file.path(out_proc_dir, sprintf("%s_%s_Fixed_Fish_SiteLevel.csv", region, year))
  fwrite(fish_site, fish_out_file, na = "")
  
  cat(green(sprintf("✅ Fish Bioerosion Pipeline Complete. Exported to: %s\n", fish_out_file)))
  return(as.data.frame(fish_site))
}

# ==============================================================================
# 5. NET BUDGET PIPELINE (TRANSECT LEVEL - SfM)
# ==============================================================================

run_net_pipeline <- function(data_dir, region, year, prod_data_list, urchin_data_list, fish_site) {
  cat(cyan(sprintf("\n🚀 STARTING NET BUDGET PIPELINE (SfM TRANSECT LEVEL): %s %s\n", region, year)))
  cat(blue("⚖️ Calculating Net Carbonate Budgets per Transect...\n"))
  
  # Extract inputs
  prod_tr <- prod_data_list$prod_transect
  urch_tr <- if (is.list(urchin_data_list) && "urch_transect" %in% names(urchin_data_list)) urchin_data_list$urch_transect else urchin_data_list
  
  # Clean fish data to ensure it aligns with process_net() expectations
  # process_net() expects columns starting with "FISH_" and matching the NCEI schema format
  fish_clean <- fish_site %>%
    mutate(
      FISH_EROSION_KG_M2_YR_ALL_MEAN = coalesce(
        if ("FISH_EROSION_KG_M2_YR_ALL_MEAN" %in% names(.)) as.numeric(FISH_EROSION_KG_M2_YR_ALL_MEAN) else NULL,
        if ("FISH_EROSION_KG_M2_YR_ALL" %in% names(.))      as.numeric(FISH_EROSION_KG_M2_YR_ALL)      else NULL,
        0
      ),
      FISH_EROSION_KG_M2_YR_ALL_SE = coalesce(
        if ("FISH_EROSION_KG_M2_YR_ALL_SE" %in% names(.))   as.numeric(FISH_EROSION_KG_M2_YR_ALL_SE)   else NULL,
        0
      )
    )
  
  # process_net joins transect data on OCC_SITEID_TRANSECT and YEAR
  # We must ensure YEAR is present and correctly formatted across all inputs
  prod_tr$YEAR    <- as.character(year)
  urch_tr$YEAR    <- as.character(year)
  fish_clean$YEAR <- as.character(year)
  
  # Run process_net 
  net_site <- process_net(
    prod   = prod_tr,
    urch   = urch_tr,
    fish   = fish_clean,
    sum_by = "transect", 
    format = "wide"       # Returns the NCEI format expected by the plotting functions
  )
  
  # Export site-level SfM net budgets
  site_out_file <- file.path(data_dir, "fixed data", "processed data", sprintf("%s_%s_Fixed_Net_SiteLevel_sfm.csv", region, year))
  fwrite(net_site, site_out_file, na = "")
  
  cat(green(sprintf("✅ PIPELINE COMPLETE! Saved SfM Net Carbonate Budgets to:\n  • %s\n", site_out_file)))
  
  # Note: net_transect is no longer returned because process_net() aggregates it internally.
  # We pass a NULL or empty dataframe for net_transect to prevent downstream errors if it's referenced elsewhere.
  return(list(net_transect = NULL, net_site = as.data.frame(net_site))) 
}

# ==============================================================================
# 6. VISUALIZATION FUNCTIONS (DYNAMIC Y-AXES & LATITUDE ORDERING)
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
    data$OCC_SITEID <- factor(data$OCC_SITEID, levels = lat_order)
  }
  
  req_cols <- c("GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_CARB_PROD_KG_M2_YR_SE",
                "URCHIN_DENSITY_NO_M2_MEAN",     "URCHIN_DENSITY_NO_M2_SE",
                "URCHIN_EROSION_KG_M2_YR_MEAN",  "URCHIN_EROSION_KG_M2_YR_SE",
                "NET_CARB_PROD_KG_M2_YR_MEAN",    "NET_CARB_PROD_KG_M2_YR_SE")
  for (col in req_cols) {
    if (!col %in% names(data)) data[[col]] <- 0
  }
  
  lims_gross <- calc_dynamic_lims(
    data$GROSS_CARB_PROD_KG_M2_YR_MEAN - data$GROSS_CARB_PROD_KG_M2_YR_SE,
    data$GROSS_CARB_PROD_KG_M2_YR_MEAN + data$GROSS_CARB_PROD_KG_M2_YR_SE
  )
  
  lims_urch_den <- calc_dynamic_lims(
    0,
    data$URCHIN_DENSITY_NO_M2_MEAN + data$URCHIN_DENSITY_NO_M2_SE,
    force_zero_lower = TRUE
  )
  
  lims_urch_eros <- calc_dynamic_lims(
    0,
    data$URCHIN_EROSION_KG_M2_YR_MEAN + data$URCHIN_EROSION_KG_M2_YR_SE,
    force_zero_lower = TRUE
  )
  
  lims_net <- calc_dynamic_lims(
    data$NET_CARB_PROD_KG_M2_YR_MEAN - data$NET_CARB_PROD_KG_M2_YR_SE,
    data$NET_CARB_PROD_KG_M2_YR_MEAN + data$NET_CARB_PROD_KG_M2_YR_SE
  )
  
  pal <- c("#003087", "#0085ca", "#1ebec7", "#76bc21", "#f0de02", "#ff8400", "#db2207")
  
  theme_benthic <- theme_bw() + 
    theme(
      axis.title.x     = element_blank(),
      axis.text.x      = element_text(size = 8, color = 'black', angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y      = element_text(size = 12, color = 'black'),
      axis.title.y     = element_text(size = 12),
      strip.text.x     = element_text(size = 12),
      legend.text      = element_text(size = 12),
      legend.title     = element_text(size = 12),
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
                        ymin = .data[[y_var]] - .data[[y_se_var]], 
                        ymax = .data[[y_var]] + .data[[y_se_var]]),
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
  req_eros <- c("MACROBIOEROSION_KG_M2_YR_MEAN", 
                "MICROBIOEROSION_KG_M2_YR_MEAN", 
                "URCHIN_EROSION_KG_M2_YR_MEAN", 
                "FISH_EROSION_KG_M2_YR_ALL_MEAN")
  for (col in req_eros) {
    if (!col %in% names(net_data)) net_data[[col]] <- 0
  }
  
  prod_clean <- prod_group %>% 
    select(any_of(c("REGION", "REGIONCODE", "ISLAND", "OCC_SITEID", "LATITUDE", "GROUP", "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN")))
  
  eros_clean <- net_data %>%
    select(any_of(c("REGION", "REGIONCODE", "ISLAND", "OCC_SITEID", "LATITUDE")), 
           Macroborers = MACROBIOEROSION_KG_M2_YR_MEAN,
           Microborers = MICROBIOEROSION_KG_M2_YR_MEAN,
           Urchins     = URCHIN_EROSION_KG_M2_YR_MEAN,
           Parrotfish  = FISH_EROSION_KG_M2_YR_ALL_MEAN) %>%
    pivot_longer(cols = c(Macroborers, Microborers, Urchins, Parrotfish),
                 names_to = "GROUP", 
                 values_to = "SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN") %>%
    mutate(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = -1 * SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN)
  
  bind_rows(prod_clean, eros_clean)
}

create_taxa_plots <- function(net_site, prod_coral, prod_site, region, year, plot_out_dir, site_order = NULL, exclude_sites = NULL) {
  cat(blue("\n🎨 Generating Taxa and Budget Plots (SfM)...\n"))
  
  prod_coral <- prod_coral %>%
    mutate(
      SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = coalesce(
        if ("SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) .data[["SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN"]] else NULL,
        if ("CORAL_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) .data[["CORAL_CARB_PROD_KG_M2_YR_MEAN"]] else NULL,
        if ("CORAL_CARB_PROD_KG_M2_YR" %in% names(.)) .data[["CORAL_CARB_PROD_KG_M2_YR"]] else NULL,
        0
      )
    )
  
  prod_coral <- prod_coral %>%
    mutate(GROUP = case_when(
      grepl("Porites", CORAL_GROUP_NAME)     ~ "Porites spp.",
      grepl("Acropora", CORAL_GROUP_NAME)    ~ "Acropora spp.",
      grepl("Montipora", CORAL_GROUP_NAME)   ~ "Montipora spp.",
      grepl("Pavona", CORAL_GROUP_NAME)      ~ "Pavona spp.",
      grepl("Pocillopora", CORAL_GROUP_NAME) ~ "Pocillopora spp.",
      TRUE                                   ~ "Other corals" 
    ))
  
  prod_coral_agg <- prod_coral %>%
    group_by(across(any_of(c("REGION", "REGIONCODE", "ISLAND", "OCC_SITEID", "LATITUDE", "GROUP")))) %>%
    summarise(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = sum(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE), .groups = "drop")
  
  cca_prod_val <- if ("CCA_CARB_PROD_KG_M2_YR_MEAN" %in% names(net_site)) {
    net_site %>% select(any_of(c("REGION", "REGIONCODE", "ISLAND", "OCC_SITEID", "LATITUDE")), SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = CCA_CARB_PROD_KG_M2_YR_MEAN)
  } else if ("CCA_CARB_PROD_KG_M2_YR" %in% names(prod_site)) {
    prod_site %>% select(any_of(c("REGION", "REGIONCODE", "ISLAND", "OCC_SITEID", "LATITUDE")), SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = CCA_CARB_PROD_KG_M2_YR)
  } else {
    net_site %>% select(any_of(c("REGION", "REGIONCODE", "ISLAND", "OCC_SITEID", "LATITUDE"))) %>% mutate(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN = 0)
  }
  cca <- cca_prod_val %>% mutate(GROUP = "CCA")
  
  prod_group <- bind_rows(prod_coral_agg, cca)
  data_group <- format_bygroup(net_site, prod_group)
  net_site_filtered <- net_site
  
  if (!is.null(exclude_sites)) {
    data_group        <- data_group %>% filter(!OCC_SITEID %in% exclude_sites)
    net_site_filtered <- net_site_filtered %>% filter(!OCC_SITEID %in% exclude_sites)
  }
  
  producer_levels <- c("CCA", "Other corals", "Acropora spp.", "Pavona spp.", 
                       "Pocillopora spp.", "Montipora spp.", "Porites spp.")
  eroder_levels   <- c("Macroborers", "Microborers", "Urchins", "Parrotfish")
  all_groups      <- c(producer_levels, eroder_levels)
  
  pal_full <- c(
    "CCA"              = "#F7B2DC", 
    "Other corals"     = "#fdffcc", 
    "Acropora spp."    = "#c5eab0", 
    "Pavona spp."      = "#7dcdbb", 
    "Pocillopora spp." = "#92bedb", 
    "Montipora spp."   = "#4989d0", 
    "Porites spp."     = "#2f2e94", 
    "Macroborers"      = "#fccf5c", 
    "Microborers"      = "#fd903b", 
    "Urchins"          = "#ea6635", 
    "Parrotfish"       = "#d8361d"
  )
  
  data_group$GROUP <- factor(data_group$GROUP, levels = all_groups)
  
  prod_data <- data_group %>% filter(GROUP %in% producer_levels)
  eros_data <- data_group %>% filter(GROUP %in% eroder_levels)
  
  if (!is.null(site_order)) {
    valid_order <- site_order[site_order %in% unique(data_group$OCC_SITEID)]
  } else if ("LATITUDE" %in% names(net_site_filtered)) {
    valid_order <- net_site_filtered %>%
      filter(!is.na(LATITUDE)) %>%
      distinct(OCC_SITEID, LATITUDE) %>%
      arrange(desc(LATITUDE)) %>%
      pull(OCC_SITEID)
  } else {
    valid_order <- unique(data_group$OCC_SITEID)
  }
  
  data_group$OCC_SITEID        <- factor(data_group$OCC_SITEID, levels = valid_order)
  net_site_filtered$OCC_SITEID <- factor(net_site_filtered$OCC_SITEID, levels = valid_order)
  prod_data$OCC_SITEID         <- factor(prod_data$OCC_SITEID, levels = valid_order)
  eros_data$OCC_SITEID         <- factor(eros_data$OCC_SITEID, levels = valid_order)
  
  site_prod_sums <- data_group %>% 
    filter(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN > 0) %>% 
    group_by(OCC_SITEID) %>% 
    summarise(total = sum(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE), .groups = "drop")
  
  site_eros_sums <- data_group %>% 
    filter(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN < 0) %>% 
    group_by(OCC_SITEID) %>% 
    summarise(total = sum(SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, na.rm = TRUE), .groups = "drop")
  
  max_prod_bar <- max(c(site_prod_sums$total, net_site_filtered$GROSS_CARB_PROD_KG_M2_YR_MEAN + net_site_filtered$GROSS_CARB_PROD_KG_M2_YR_SE), na.rm = TRUE)
  min_eros_bar <- min(c(site_eros_sums$total, -1 * (net_site_filtered$GROSS_EROSION_KG_M2_YR_MEAN + net_site_filtered$GROSS_EROSION_KG_M2_YR_SE)), na.rm = TRUE)
  
  min_net_val <- min(net_site_filtered$NET_CARB_PROD_KG_M2_YR_MEAN - net_site_filtered$NET_CARB_PROD_KG_M2_YR_SE, na.rm = TRUE)
  max_net_val <- max(net_site_filtered$NET_CARB_PROD_KG_M2_YR_MEAN + net_site_filtered$NET_CARB_PROD_KG_M2_YR_SE, na.rm = TRUE)
  
  overall_min <- min(min_eros_bar, min_net_val, 0, na.rm = TRUE)
  overall_max <- max(max_prod_bar, max_net_val, 0, na.rm = TRUE)
  
  taxa_y_lims <- c(floor(overall_min - abs(overall_min) * 0.12), ceiling(overall_max + abs(overall_max) * 0.12))
  
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
  
  prod_bytaxa_plot <- ggplot() +
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
                      ymax = GROSS_CARB_PROD_KG_M2_YR_MEAN + GROSS_CARB_PROD_KG_M2_YR_SE), 
                  width = 0.25) +
    geom_errorbar(data = net_site_filtered,
                  aes(x = OCC_SITEID, 
                      ymin = -1 * GROSS_EROSION_KG_M2_YR_MEAN - GROSS_EROSION_KG_M2_YR_SE, 
                      ymax = -1 * GROSS_EROSION_KG_M2_YR_MEAN), 
                  width = 0.25) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = taxa_y_lims, expand = c(0, 0)) +
    scale_fill_manual(values = pal_full) +
    ylab(expression(paste("Carbonate production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_budget
  
  net_plot_gradient <- ggplot() +
    geom_hline(yintercept = 0) +
    geom_col(data = net_site_filtered,
             aes(x = OCC_SITEID, y = NET_CARB_PROD_KG_M2_YR_MEAN, fill = NET_CARB_PROD_KG_M2_YR_MEAN),
             color = 'black', width = 0.5) +
    geom_errorbar(data = net_site_filtered,
                  aes(x = OCC_SITEID, 
                      ymin = NET_CARB_PROD_KG_M2_YR_MEAN - NET_CARB_PROD_KG_M2_YR_SE, 
                      ymax = NET_CARB_PROD_KG_M2_YR_MEAN + NET_CARB_PROD_KG_M2_YR_SE), 
                  width = 0.25) +
    facet_grid(~REGIONCODE, scales = "free", space = "free") +
    scale_y_continuous(limits = taxa_y_lims, expand = c(0, 0)) +
    scale_fill_gradient2(low = "#EF2820", mid = "#CCCCCC", high = "#38A1D0", midpoint = 0, limits = taxa_y_lims) +
    ylab(expression(paste("Net carbonate production rate (", "kg CaCO"[3], " m" ^ -2, " yr" ^ -1, ")  "))) +
    theme_budget
  
  legend_prod <- get_legend(
    ggplot(prod_data, aes(x = OCC_SITEID, y = SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, fill = GROUP)) +
      geom_col(color = "black") +
      scale_fill_manual(values = pal_full[producer_levels], name = "Producers") +
      guides(fill = guide_legend(byrow = TRUE)) +
      theme(legend.text = element_text(size = 15), legend.title = element_text(size = 22), 
            legend.spacing.y = unit(0.5, 'cm'), legend.key.size = unit(1, "cm"))
  )
  
  legend_eros <- get_legend(
    ggplot(eros_data, aes(x = OCC_SITEID, y = SUBSTRATE_CARB_PROD_KG_M2_YR_MEAN, fill = GROUP)) +
      geom_col(color = "black") +
      scale_fill_manual(values = pal_full[eroder_levels], name = "Eroders") +
      guides(fill = guide_legend(reverse = TRUE, byrow = TRUE)) +
      theme(legend.text = element_text(size = 15), legend.title = element_text(size = 22), 
            legend.spacing.y = unit(0.5, 'cm'), legend.key.size = unit(1, "cm"))
  )
  
  combined_legend <- ggarrange(NULL, legend_prod, legend_eros, NULL, ncol = 1, heights = c(0.1, 1, 1, 0.1))
  
  final_plot <- ggarrange(
    net_plot_gradient, NULL, prod_bytaxa_plot, combined_legend,
    widths     = c(1, 0.05, 1, 0.25),
    nrow       = 1,
    labels     = c("a", "", "b", ""),
    font.label = list(size = 24, face = "bold")
  )
  
  out_path <- file.path(plot_out_dir, sprintf("%s_%s_netprod_fixed_sfm.png", region, year))
  ggsave(final_plot, filename = out_path, bg = "white", height = 12, width = 24, dpi = 300)
  
  cat(green(sprintf("✅ Composite plot saved to: %s\n", out_path)))
  return(final_plot)
}

# ==============================================================================
# 7. EXECUTE MASTER PIPELINE
# ==============================================================================

# A. RUN BENTHIC PRODUCTION PIPELINE (SfM)
prod_data_list <- run_benthic_pipeline(data_dir, region, year, meta_file, prod_dbase_file, benthic_file)

# B. RUN URCHIN BIOEROSION PIPELINE
urchin_data_list <- run_urchin_pipeline(data_dir, region, year, meta_file, urchin_file)

# C. RUN FISH BIOEROSION PIPELINE
fish_data <- run_fish_pipeline(data_dir, region, year, meta_file, fish_file)

# D. RUN NET BUDGET PIPELINE (SfM)
net_results <- run_net_pipeline(
  data_dir         = data_dir, 
  region           = region, 
  year             = year, 
  prod_data_list   = prod_data_list, 
  urchin_data_list = urchin_data_list, 
  fish_site        = fish_data
)

final_net_transect <- net_results$net_transect
final_net_site     <- net_results$net_site

# E. DERIVE SITE ORDER BY DECREASING LATITUDE (North to South)
latitude_site_order <- final_net_site %>%
  filter(!is.na(LATITUDE)) %>%
  mutate(REGIONCODE = factor(REGIONCODE, levels = c("NWHI", "MHI"))) %>%
  arrange(REGIONCODE, desc(LATITUDE)) %>%
  distinct(OCC_SITEID) %>%
  pull(OCC_SITEID)

# F. GENERATE & SAVE OVERVIEW VISUALIZATIONS (SfM)
overview_plots <- create_benthic_plots(
  data          = final_net_site, 
  site_order    = latitude_site_order, 
  exclude_sites = sites_to_exclude
)

ggsave(file.path(plot_out_dir, sprintf("%s_%s_gross_overview_sfm.png", region, year)), 
       overview_plots$gross_plot, width = 14, height = 6, bg = "white", dpi = 300)
ggsave(file.path(plot_out_dir, sprintf("%s_%s_net_overview_sfm.png", region, year)),   
       overview_plots$net_plot,   width = 14, height = 6, bg = "white", dpi = 300)

# G. GENERATE & SAVE TAXA COMPOSITE PLOT (SfM)
final_taxa_plot <- create_taxa_plots(
  net_site      = final_net_site, 
  prod_coral    = prod_data_list$prod_coral, 
  prod_site     = prod_data_list$prod_site, 
  region        = region, 
  year          = year, 
  plot_out_dir  = plot_out_dir, 
  site_order    = latitude_site_order, 
  exclude_sites = sites_to_exclude
)