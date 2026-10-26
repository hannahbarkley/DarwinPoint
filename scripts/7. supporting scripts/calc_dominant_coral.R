# ==============================================================================
# Script:  calculate_dominant_coral_taxa_morphology.R
# Purpose: Calculate the most common coral taxa (T3_DESC) and morphologies across 
#          the 2024 Hawaiian StRS survey domain from Tier-3 benthic cover data.
# Inputs:  - 2024_HA_StRS_Cover_Tier3.csv
#          - HA_Benthic_Codes.csv
# Outputs: - Ranked tables for Morphology, T3_DESC, and T3_DESC x Morphology
# ==============================================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------------------------
# 1. LOAD INPUT DATASETS
# ------------------------------------------------------------------------------
find_input_file <- function(candidates) {
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  return(candidates[1])
}

cover_csv <- find_input_file(c(
  here("data", "strs data", "field data", "2024_HA_StRS_Cover_Tier3.csv"),
  here("data", "2024_HA_StRS_Cover_Tier3.csv"),
  here("2024_HA_StRS_Cover_Tier3.csv")
))

lookup_csv <- find_input_file(c(
  here("data", "supporting data", "HA_Benthic_Codes.csv"),
  here("data", "dbase", "HA_Benthic_Codes.csv"),
  here("data", "HA_Benthic_Codes.csv"),
  here("HA_Benthic_Codes.csv")
))

stopifnot("2024_HA_StRS_Cover_Tier3.csv not found!" = file.exists(cover_csv))
stopifnot("HA_Benthic_Codes.csv not found!" = file.exists(lookup_csv))

cover_df  <- read_csv(cover_csv, show_col_types = FALSE)
lookup_df <- read_csv(lookup_csv, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. HARMONIZE TAXONOMIC & MORPHOLOGY MAPPING VIA T3_DESC
# ------------------------------------------------------------------------------
meta_cols    <- c("SITEVISITID", "SITE", "ISLAND", "REGION", "REEF_ZONE", "DEPTH_BIN", "LATITUDE", "LONGITUDE")
meta_present <- intersect(meta_cols, names(cover_df))

ignore_cols <- c(
  meta_present, "OBS_YEAR", "HABITAT_CODE", "DATE_", "MISSIONID", 
  "OCC_SITEID", "ANALYSIS_YEAR", "ANALYSIS_SCHEME", "SEC_NAME", 
  "EXCLUDE_FLAG", "TRANSECT_PHOTOS", "Oceanography", 
  "new_MIN_DEPTH_M", "new_MAX_DEPTH_M", "METHOD", 
  "PERM_SITE", "CLIMATE_STATION_YN", "N", "new.N"
)
code_cols <- setdiff(names(cover_df), ignore_cols)

# Extract CODE, TIER_1, TIER_2, and T3_DESC from lookup table
code_map <- lookup_df %>%
  select(CODE, TIER_1, TIER_2, any_of("T3_DESC")) %>%
  distinct(CODE, .keep_all = TRUE) %>%
  bind_rows(tibble(
    CODE    = c("HAL", "HYCO"), 
    TIER_1  = c("ALGAE", "CORAL"), 
    TIER_2  = c("UPMA", "HYCO"),
    T3_DESC = c("Halimeda", "Hydrocoral")
  ))

# Pivot observations and map to T3_DESC and MORPHOLOGY
long_coral <- cover_df %>%
  pivot_longer(cols = all_of(code_cols), names_to = "CODE", values_to = "COVER_PERC") %>%
  left_join(code_map, by = "CODE") %>%
  filter(TIER_1 == "CORAL") %>%
  mutate(
    MORPHOLOGY = case_when(
      TIER_2 == "BR"   ~ "Branching",
      TIER_2 == "MASS" ~ "Massive",
      TIER_2 == "ENC"  ~ "Encrusting",
      TIER_2 == "FOL"  ~ "Foliose",
      TIER_2 == "TB"   ~ "Table",
      TIER_2 == "FREE" ~ "Free-living",
      TRUE             ~ "Other / Unclassified"
    ),
    T3_DESC = coalesce(T3_DESC, CODE)
  )

total_sites <- n_distinct(cover_df$SITE)
total_coral_cover_sum <- sum(long_coral$COVER_PERC, na.rm = TRUE)

# ------------------------------------------------------------------------------
# 3. RANKING: MOST COMMON CORAL MORPHOLOGIES
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("  1. MOST COMMON CORAL MORPHOLOGIES (ALL HAWAIIAN StRS SITES)\n")
cat("==============================================================================\n")

morph_summary <- long_coral %>%
  group_by(MORPHOLOGY) %>%
  summarize(
    Mean_Site_Cover_Pct = round(sum(COVER_PERC, na.rm = TRUE) / total_sites, 2),
    Rel_Share_Of_Coral  = round((sum(COVER_PERC, na.rm = TRUE) / total_coral_cover_sum) * 100, 1),
    Site_Prevalence_Pct = round((n_distinct(SITE[COVER_PERC > 0]) / total_sites) * 100, 1),
    .groups             = "drop"
  ) %>%
  arrange(desc(Mean_Site_Cover_Pct))

print(as.data.frame(morph_summary), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. RANKING: MOST COMMON CORAL TAXA (T3_DESC)
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("  2. MOST COMMON CORAL TAXA (T3_DESC)\n")
cat("==============================================================================\n")

t3_summary <- long_coral %>%
  group_by(T3_DESC) %>%
  summarize(
    Mean_Site_Cover_Pct = round(sum(COVER_PERC, na.rm = TRUE) / total_sites, 2),
    Rel_Share_Of_Coral  = round((sum(COVER_PERC, na.rm = TRUE) / total_coral_cover_sum) * 100, 1),
    Site_Prevalence_Pct = round((n_distinct(SITE[COVER_PERC > 0]) / total_sites) * 100, 1),
    .groups             = "drop"
  ) %>%
  arrange(desc(Mean_Site_Cover_Pct))

print(head(as.data.frame(t3_summary), 20), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. RANKING: T3_DESC x MORPHOLOGY FUNCTIONAL COMBINATIONS
# ------------------------------------------------------------------------------
cat("\n==============================================================================\n")
cat("  3. MOST COMMON T3_DESC x MORPHOLOGY COMBINATIONS\n")
cat("==============================================================================\n")

combo_summary <- long_coral %>%
  group_by(T3_DESC, MORPHOLOGY) %>%
  summarize(
    Mean_Site_Cover_Pct = round(sum(COVER_PERC, na.rm = TRUE) / total_sites, 2),
    Rel_Share_Of_Coral  = round((sum(COVER_PERC, na.rm = TRUE) / total_coral_cover_sum) * 100, 1),
    Site_Prevalence_Pct = round((n_distinct(SITE[COVER_PERC > 0]) / total_sites) * 100, 1),
    .groups             = "drop"
  ) %>%
  arrange(desc(Mean_Site_Cover_Pct))

print(head(as.data.frame(combo_summary), 20), row.names = FALSE)