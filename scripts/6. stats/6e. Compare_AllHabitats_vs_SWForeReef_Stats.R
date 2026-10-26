# ==============================================================================
# Script:  6e. Compare_AllHabitats_vs_SWForeReef_Stats.R
# Purpose: Statistical quantification of the offset, percentage inflation, and 
#          fold-difference between representative island-scale all-habitat means 
#          (StRS) and historical shallow-water (SW) fore-reef 90th percentiles 
#          (Grigg-equivalent subset) for Coral Cover and Gross Production.
# Inputs:  - data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          - data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv 
# Outputs: - outputs/stats/HA_2024_Sensitivity_AllHabitats_vs_SW90th_Stats.csv
#          - outputs/stats/stats_sensitivity_all_habitats_vs_sw_forereef.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------
# 1. DIRECTORIES & FILE INGESTION HELPERS
# ------------------------------------------------------------------------------
out_dir <- here("outputs", "stats")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

island_lookup <- c(
  "kure"                         = "H\u014dlanik\u016b",
  "kure atoll"                   = "H\u014dlanik\u016b",
  "midway"                       = "Kuaihelani",
  "midway atoll"                 = "Kuaihelani",
  "pearl & hermes"               = "Manawai",
  "pearl and hermes"             = "Manawai",
  "pearl & hermes atoll"         = "Manawai",
  "pearl_hermes"                 = "Manawai",
  "lisianski"                    = "Kapou",
  "lisianski island"             = "Kapou",
  "laysan"                       = "Kamole",
  "laysan island"                = "Kamole",
  "maro"                         = "Kamokuokamohoali'i",
  "maro reef"                    = "Kamokuokamohoali'i",
  "gardner"                      = "P\u016bh\u0101honu",
  "gardner pinnacles"            = "P\u016bh\u0101honu",
  "gardner_pinnacles"            = "P\u016bh\u0101honu",
  "french frigate"               = "Lalo",
  "french frigate shoals"        = "Lalo",
  "french_frigate"               = "Lalo",
  "french_frigate_shoals"        = "Lalo",
  "ffs"                          = "Lalo",
  "necker"                       = "Mokumanamana",
  "necker island"                = "Mokumanamana",
  "nihoa"                        = "Nihoa",
  "nihoa island"                 = "Nihoa",
  "kauai"                        = "Kaua'i",
  "niihau"                       = "Ni'ihau",
  "oahu"                         = "O'ahu",
  "molokai"                      = "Moloka'i",
  "lanai"                        = "L\u0101na'i",
  "maui"                         = "Maui",
  "kahoolawe"                    = "Kaho'olawe",
  "hawaii"                       = "Hawai'i"
)

format_islands <- function(x) {
  clean_key <- tolower(trimws(as.character(x)))
  clean_key <- gsub("_", " ", clean_key)
  mapped    <- island_lookup[clean_key]
  coalesce(unname(mapped), as.character(x))
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "strs data", "processed data"),
    here("data", "grigg data"),
    here("data", "strs data", "field data"),
    here("data"),
    here()
  )
  for (f in candidates) {
    for (d in search_dirs) {
      target <- file.path(d, f)
      if (file.exists(target)) return(read_csv(target, show_col_types = FALSE))
    }
  }
  stop(paste("None of the candidate files found:", paste(candidates, collapse = " | ")))
}

# ------------------------------------------------------------------------------
# 2. LOAD & HARMONIZE DATASETS
# ------------------------------------------------------------------------------
df_all_raw <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>%
  mutate(
    LOCATION = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                        if ("ISLAND" %in% names(.)) ISLAND else NULL),
    ISLAND_LABEL = format_islands(LOCATION),
    CORAL_COVER  = coalesce(if ("HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) HARD_CORAL_COVER_PCT_MEAN else NULL,
                            if ("CORAL_COVER_PERC" %in% names(.)) CORAL_COVER_PERC else NULL,
                            if ("HARD_CORAL_COVER_PCT" %in% names(.)) HARD_CORAL_COVER_PCT else NULL),
    GROSS_PROD   = coalesce(if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL,
                            if ("GROSS_PROD_KG_M2_YR" %in% names(.)) GROSS_PROD_KG_M2_YR else NULL,
                            if ("PROD" %in% names(.)) PROD else NULL)
  ) %>%
  select(ISLAND_LABEL, LATITUDE, ALL_CORAL_COVER = CORAL_COVER, ALL_GROSS_PROD = GROSS_PROD) %>%
  distinct(ISLAND_LABEL, .keep_all = TRUE)

df_sw_raw <- saferead_csv("HA_2024_Grigg_Equivalent_Sites.csv") %>%
  mutate(
    LOCATION = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                        if ("ISLAND" %in% names(.)) ISLAND else NULL),
    ISLAND_LABEL = format_islands(LOCATION),
    CORAL_COVER  = coalesce(if ("HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) HARD_CORAL_COVER_PCT_MEAN else NULL,
                            if ("CORAL_COVER_PERC" %in% names(.)) CORAL_COVER_PERC else NULL,
                            if ("HARD_CORAL_COVER_PCT" %in% names(.)) HARD_CORAL_COVER_PCT else NULL),
    GROSS_PROD   = coalesce(if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL,
                            if ("GROSS_PROD_KG_M2_YR" %in% names(.)) GROSS_PROD_KG_M2_YR else NULL,
                            if ("PROD" %in% names(.)) PROD else NULL)
  ) %>%
  select(ISLAND_LABEL, SITE = any_of(c("SITE", "OCC_SITEID")), CORAL_COVER, GROSS_PROD)

# ------------------------------------------------------------------------------
# 3. COMPUTE STATISTICAL OFFSETS (ISLAND-LEVEL)
# ------------------------------------------------------------------------------
compute_metric_stats <- function(sw_df, all_df, metric_col, all_col, metric_label, unit_str) {
  sw_summary <- sw_df %>%
    group_by(ISLAND_LABEL) %>%
    summarize(
      N_SW_Sites = n(),
      SW_Mean    = mean(.data[[metric_col]], na.rm = TRUE),
      SW_Median  = median(.data[[metric_col]], na.rm = TRUE),
      SW_75th    = quantile(.data[[metric_col]], probs = 0.75, na.rm = TRUE),
      SW_90th    = quantile(.data[[metric_col]], probs = 0.90, na.rm = TRUE),
      SW_Max     = max(.data[[metric_col]], na.rm = TRUE),
      .groups    = "drop"
    )
  
  comp <- inner_join(all_df %>% select(ISLAND_LABEL, LATITUDE, All_Mean = all_of(all_col)), 
                     sw_summary, 
                     by = "ISLAND_LABEL") %>%
    mutate(
      Metric             = metric_label,
      Unit               = unit_str,
      Absolute_Offset    = SW_90th - All_Mean,
      Percent_Offset     = ((SW_90th - All_Mean) / All_Mean) * 100,
      Fold_Difference    = SW_90th / All_Mean
    ) %>%
    arrange(desc(LATITUDE))
  
  return(comp)
}

stats_coral <- compute_metric_stats(df_sw_raw, df_all_raw, "CORAL_COVER", "ALL_CORAL_COVER", "Coral Cover", "%")
stats_prod  <- compute_metric_stats(df_sw_raw, df_all_raw, "GROSS_PROD", "ALL_GROSS_PROD", "Gross Production", "kg m^-2 yr^-1")

combined_stats <- bind_rows(stats_prod, stats_coral) %>%
  relocate(Metric, ISLAND_LABEL, LATITUDE, N_SW_Sites, All_Mean, SW_Mean, SW_75th, SW_90th, SW_Max, Absolute_Offset, Percent_Offset, Fold_Difference, Unit)

# ------------------------------------------------------------------------------
# 4. ARCHIPELAGO HYPOTHESIS TESTING & AGGREGATE SUMMARY
# ------------------------------------------------------------------------------
# Paired differences across islands (testing whether SW 90th systematically inflates island estimates)
t_prod  <- t.test(stats_prod$SW_90th, stats_prod$All_Mean, paired = TRUE)
w_prod  <- wilcox.test(stats_prod$SW_90th, stats_prod$All_Mean, paired = TRUE, exact = FALSE)

t_coral <- t.test(stats_coral$SW_90th, stats_coral$All_Mean, paired = TRUE)
w_coral <- wilcox.test(stats_coral$SW_90th, stats_coral$All_Mean, paired = TRUE, exact = FALSE)

# Generate detailed text log
report_lines <- c(
  "================================================================================",
  "STATISTICAL SENSITIVITY REPORT: All-Habitat Means vs. SW Fore-Reef 90th Percentile",
  "================================================================================",
  sprintf("Number of compared islands: %d", nrow(stats_prod)),
  "",
  "1. GROSS PRODUCTION RATE (kg m^-2 yr^-1):",
  sprintf("   - Mean Archipelago All-Habitat Baseline:  %.2f kg m^-2 yr^-1", mean(stats_prod$All_Mean, na.rm = TRUE)),
  sprintf("   - Mean SW Fore-Reef 90th Percentile:      %.2f kg m^-2 yr^-1", mean(stats_prod$SW_90th, na.rm = TRUE)),
  sprintf("   - Mean Absolute Offset:                   %+.2f +/- %.2f kg m^-2 yr^-1 (mean +/- SD)", 
          mean(stats_prod$Absolute_Offset), sd(stats_prod$Absolute_Offset)),
  sprintf("   - Mean Percent Inflation:                 %+.1f%% +/- %.1f%% (Median: %+.1f%%, IQR: [%.1f%%, %.1f%%])", 
          mean(stats_prod$Percent_Offset), sd(stats_prod$Percent_Offset),
          median(stats_prod$Percent_Offset), quantile(stats_prod$Percent_Offset, 0.25), quantile(stats_prod$Percent_Offset, 0.75)),
  sprintf("   - Mean Fold-Difference:                   %.2fx (Range: %.2fx - %.2fx)", 
          mean(stats_prod$Fold_Difference), min(stats_prod$Fold_Difference), max(stats_prod$Fold_Difference)),
  sprintf("   - Paired t-test:                          t(%d) = %.3f, p = %.4e", 
          t_prod$parameter, t_prod$statistic, t_prod$p.value),
  sprintf("   - Wilcoxon Signed-Rank:                   W = %.1f, p = %.4e", 
          w_prod$statistic, w_prod$p.value),
  "",
  "2. HARD CORAL COVER (%):",
  sprintf("   - Mean Archipelago All-Habitat Baseline:  %.2f%%", mean(stats_coral$All_Mean, na.rm = TRUE)),
  sprintf("   - Mean SW Fore-Reef 90th Percentile:      %.2f%%", mean(stats_coral$SW_90th, na.rm = TRUE)),
  sprintf("   - Mean Absolute Offset:                   %+.2f%% +/- %.2f%%", 
          mean(stats_coral$Absolute_Offset), sd(stats_coral$Absolute_Offset)),
  sprintf("   - Mean Percent Inflation:                 %+.1f%% +/- %.1f%% (Median: %+.1f%%, IQR: [%.1f%%, %.1f%%])", 
          mean(stats_coral$Percent_Offset), sd(stats_coral$Percent_Offset),
          median(stats_coral$Percent_Offset), quantile(stats_coral$Percent_Offset, 0.25), quantile(stats_coral$Percent_Offset, 0.75)),
  sprintf("   - Mean Fold-Difference:                   %.2fx (Range: %.2fx - %.2fx)", 
          mean(stats_coral$Fold_Difference), min(stats_coral$Fold_Difference), max(stats_coral$Fold_Difference)),
  sprintf("   - Paired t-test:                          t(%d) = %.3f, p = %.4e", 
          t_coral$parameter, t_coral$statistic, t_coral$p.value),
  sprintf("   - Wilcoxon Signed-Rank:                   W = %.1f, p = %.4e", 
          w_coral$statistic, w_coral$p.value),
  "================================================================================"
)

cat(paste(report_lines, collapse = "\n"), "\n\n")

# ------------------------------------------------------------------------------
# 5. CONSOLE TABLE & EXPORT
# ------------------------------------------------------------------------------
cat("=== ISLAND-LEVEL SENSITIVITY SUMMARY TABLE ===\n")
print(
  combined_stats %>%
    select(Metric, ISLAND_LABEL, All_Mean, SW_90th, Absolute_Offset, Percent_Offset, Fold_Difference),
  n = Inf
)

out_txt <- file.path(out_dir, "stats_sensitivity_all_habitats_vs_sw_forereef.txt")
out_csv <- file.path(out_dir, "HA_2024_Sensitivity_AllHabitats_vs_SW90th_Stats.csv")

writeLines(report_lines, out_txt)
write_csv(combined_stats, out_csv)

cat(sprintf("\n✅ Summary report saved to: %s\n", out_txt))
cat(sprintf("✅ Island-level comparative stats exported to: %s\n", out_csv))

# ------------------------------------------------------------------------------
# SENSITIVITY OFFSETS & TEST STATISTICS (Mean ± SE & Mean ± SD)
# ------------------------------------------------------------------------------
# 1. Hard Coral Cover
mean_diff_coral <- unname(t_coral$estimate)
se_diff_coral   <- t_coral$stderr
sd_diff_coral   <- sd(stats_coral$Absolute_Offset, na.rm = TRUE)

# 2. Gross Carbonate Production
mean_diff_gprod <- unname(t_prod$estimate)
se_diff_gprod   <- t_prod$stderr
sd_diff_gprod   <- sd(stats_prod$Absolute_Offset, na.rm = TRUE)

# Console Output: Mean ± SE and Mean ± SD
cat("\n=== SENSITIVITY OFFSETS (SW 90th vs. All-Habitat Baseline) ===\n")
cat(sprintf("Coral cover:      +%4.1f ± %.1f%%             (mean ± SE) | ± %.1f%% (mean ± SD)\n", 
            mean_diff_coral, se_diff_coral, sd_diff_coral))
cat(sprintf("Gross production: +%4.2f ± %.2f kg m⁻² yr⁻¹ (mean ± SE) | ± %.2f kg m⁻² yr⁻¹ (mean ± SD)\n", 
            mean_diff_gprod, se_diff_gprod, sd_diff_gprod))
