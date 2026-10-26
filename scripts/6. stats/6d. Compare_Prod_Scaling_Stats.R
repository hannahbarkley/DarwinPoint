# ==============================================================================
# Script:  6d. Compare_Prod_Scaling_Stats.R
# Purpose: Statistically compare scaled vs. unscaled Gross Carbonate Production 
#          (PROD) estimates. Computes site-level paired tests and survey-weighted 
#          differences at the overall (archipelago) and island levels.
# Inputs:  - HA_2024_StRS_Prod_SiteLevel.csv
#          - HA_2024_StRS_Prod_SiteLevel_notscaled.csv
#          - Sectors-Strata-Areas.csv
# Outputs: - Console summary of paired t-test and archipelago statistics
#          - outputs/stats/HA_2024_StRS_PROD_Scaling_Comparison_Stats.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survey)
  library(here)
})

# ------------------------------------------------------------------------------
# 1. FILE LOCATOR FUNCTION (Matches Figure S3)
# ------------------------------------------------------------------------------
saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "strs data", "field data"),
    here("data", "strs data", "processed data"),
    here("data", "supporting data"),
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
# 2. LOAD & MERGE DATA
# ------------------------------------------------------------------------------
df_scaled   <- saferead_csv("HA_2024_StRS_Prod_SiteLevel.csv") %>%
  select(SITE, OBS_YEAR, REGION, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, PROD_Scaled = PROD)

df_unscaled <- saferead_csv("HA_2024_StRS_Prod_SiteLevel_notscaled.csv") %>%
  select(SITE, PROD_Unscaled = PROD)

sectors <- saferead_csv("Sectors-Strata-Areas.csv") %>%
  select(SEC_NAME, REEF_ZONE, DEPTH_BIN, NH) %>%
  distinct()

# Create paired site-level dataset
df_paired <- inner_join(df_scaled, df_unscaled, by = "SITE") %>%
  mutate(Diff_Site = PROD_Scaled - PROD_Unscaled)

# ------------------------------------------------------------------------------
# 3. SITE-LEVEL PAIRED STATISTICS (UNWEIGHTED)
# ------------------------------------------------------------------------------
cat("\n=== SITE-LEVEL PAIRED STATISTICAL TEST ===\n")
paired_test <- t.test(df_paired$PROD_Scaled, df_paired$PROD_Unscaled, paired = TRUE)
print(paired_test)

# ------------------------------------------------------------------------------
# 4. SURVEY-WEIGHTED STATISTICS
# ------------------------------------------------------------------------------
# Join strata areas and calculate sampling weights (NH / n)
df_weighted <- df_paired %>%
  left_join(sectors, by = c("SEC_NAME", "REEF_ZONE", "DEPTH_BIN")) %>%
  mutate(
    n = n_distinct(SITE),
    sw = NH / n,
    .by = c(OBS_YEAR, SEC_NAME, REEF_ZONE, DEPTH_BIN, NH)
  ) %>%
  filter(!is.na(sw) & n > 1) %>%
  mutate(
    Strat_conc = paste(OBS_YEAR, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, sep = "_"),
    ISLAND = as.factor(ISLAND)
  )

# Construct survey design
des <- svydesign(id = ~1, strata = ~Strat_conc, weights = ~sw, data = df_weighted)

# --- Overall Archipelago Comparison ---
archipelago_mean_scaled   <- svymean(~PROD_Scaled, des, na.rm = TRUE)
archipelago_mean_unscaled <- svymean(~PROD_Unscaled, des, na.rm = TRUE)
archipelago_diff          <- svymean(~Diff_Site, des, na.rm = TRUE)

cat("\n=== SURVEY-WEIGHTED ARCHIPELAGO MEANS ===\n")
cat(sprintf("Scaled Mean:           %.3f (SE: %.3f)\n", coef(archipelago_mean_scaled), SE(archipelago_mean_scaled)))
cat(sprintf("Unscaled Mean:         %.3f (SE: %.3f)\n", coef(archipelago_mean_unscaled), SE(archipelago_mean_unscaled)))
cat(sprintf("Archipelago Mean Diff: %.3f (SE: %.3f)\n", coef(archipelago_diff), SE(archipelago_diff)))

# --- Island-Level Comparison ---
island_scaled   <- svyby(~PROD_Scaled, ~ISLAND, des, svymean, na.rm = TRUE) %>% as_tibble()
island_unscaled <- svyby(~PROD_Unscaled, ~ISLAND, des, svymean, na.rm = TRUE) %>% as_tibble()
island_diff     <- svyby(~Diff_Site, ~ISLAND, des, svymean, na.rm = TRUE) %>% as_tibble() %>%
  rename(Mean_Difference = Diff_Site, SE_Difference = se)

island_stats <- left_join(
  island_scaled %>% rename(Mean_Scaled = PROD_Scaled, SE_Scaled = se),
  island_unscaled %>% rename(Mean_Unscaled = PROD_Unscaled, SE_Unscaled = se),
  by = "ISLAND"
) %>%
  left_join(island_diff, by = "ISLAND") %>%
  mutate(
    Percent_Change = (Mean_Difference / Mean_Unscaled) * 100
  ) %>%
  arrange(desc(Percent_Change)) %>%
  relocate(ISLAND, Mean_Scaled, SE_Scaled, Mean_Unscaled, SE_Unscaled, Mean_Difference, SE_Difference, Percent_Change)

cat("\n=== SURVEY-WEIGHTED ISLAND-LEVEL STATISTICS ===\n")
print(island_stats, n = Inf)

# ------------------------------------------------------------------------------
# 5. EXPORT RESULTS
# ------------------------------------------------------------------------------
out_dir <- here("outputs", "stats")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, "HA_2024_StRS_PROD_Scaling_Comparison_Stats.csv")
write_csv(island_stats, out_file)
cat(sprintf("\n✅ Island-level comparative statistics exported to: %s\n", out_file))