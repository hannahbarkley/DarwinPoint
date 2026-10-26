# ==============================================================================
# Script: calculate_t3_colony_length_summary.R
# Project: Coral Reef Ecosystem Demography & Benthic Assessment (StRS 2025)
# Purpose: Calculate mean colony length, standard error, and sample size for 
#          key Tier-3 (T3) coral taxa/morphology groups using field demography
#          data (HA_2025_StRS_demography.csv) and merge into combo_summary.
# ==============================================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------------------------
# 1. Load Coral Demography Data
# ------------------------------------------------------------------------------
demo_file <- here("data", "strs data", "field data", "HA_2025_StRS_demography.csv")
if (!file.exists(demo_file)) {
  demo_file <- here("HA_2025_StRS_demography.csv")
}

demo_df <- read_csv(demo_file, show_col_types = FALSE)

# ------------------------------------------------------------------------------
# 2. Harmonize Taxa & Morphology to Match T3 Functional Groups
# ------------------------------------------------------------------------------
# In NOAA NCRMP field demography:
#   MD = Mounding -> maps to "Massive"
#   EM = Encrusting (mounding) -> maps to "Encrusting"
#   BR = Branching -> maps to "Branching"

target_demography <- demo_df %>%
  filter(SPCODE != "AAAA", !is.na(COLONYLENGTH), !is.na(MORPH_CODE)) %>%
  mutate(
    Genus = word(TAXONNAME, 1),
    T3_DESC = case_when(
      Genus == "Porites"     & MORPH_CODE == "MD" ~ "Porites massive",
      Genus == "Montipora"   & MORPH_CODE == "EM" ~ "Montipora encrusting",
      Genus == "Porites"     & MORPH_CODE == "BR" ~ "Porites branching",
      Genus == "Porites"     & MORPH_CODE == "EM" ~ "Porites encrusting",
      Genus == "Pocillopora" & MORPH_CODE == "BR" ~ "Pocillopora sp",
      TRUE ~ NA_character_
    ),
    MORPHOLOGY = case_when(
      MORPH_CODE == "MD" ~ "Massive",
      MORPH_CODE == "EM" ~ "Encrusting",
      MORPH_CODE == "BR" ~ "Branching",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(T3_DESC))

# ------------------------------------------------------------------------------
# 3. Compute Summary Length Statistics
# ------------------------------------------------------------------------------
taxa_length_summary <- target_demography %>%
  summarize(
    N_Colonies    = n(),
    Mean_Length   = round(mean(COLONYLENGTH), 2),
    Median_Length = round(median(COLONYLENGTH), 2),
    SD_Length     = round(sd(COLONYLENGTH), 2),
    SE_Length     = round(sd(COLONYLENGTH) / sqrt(n()), 2),
    Min_Length    = round(min(COLONYLENGTH), 2),
    Max_Length    = round(max(COLONYLENGTH), 2),
    .by           = c(T3_DESC, MORPHOLOGY)
  )

# ------------------------------------------------------------------------------
# 4. Merge Length Metrics with Benthic Cover Summary
# ------------------------------------------------------------------------------
combo_summary_updated <- combo_summary %>%
  left_join(taxa_length_summary, by = c("T3_DESC", "MORPHOLOGY"))

# ------------------------------------------------------------------------------
# 5. Export and Preview Output
# ------------------------------------------------------------------------------
output_dir <- here("outputs")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write_csv(
  combo_summary_updated,
  file.path(output_dir, "StRS_t3_taxa_cover_and_colony_length.csv")
)

print(as.data.frame(combo_summary_updated), row.names = FALSE)