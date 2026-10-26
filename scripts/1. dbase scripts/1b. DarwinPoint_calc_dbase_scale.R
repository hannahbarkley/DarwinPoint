# ==============================================================================
# Script: 2b. DarwinPoint_calc_dbase_scale.R
# Purpose: 
#   1. Load empirical reference data (prod_dbase_ha_ref.csv) and NCRMP schema.
#   2. Dynamically generate an NCRMP template crossing islands with Corals & CCA.
#   3. Apply Grigg (1982) latitudinal proportional scaling.
#   4. Integrate CCA Calcification rates by island.
#   5. Append static, non-coral benthic substrates at the end (LOCATION = NA).
#   6. Round all numeric outputs to 3 decimal places.
# Outputs: 
#   - data/prod_dbase_ha_fixed.csv
#   - data/prod_dbase_ha_strs.csv
# ==============================================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------------------------
# HELPER FUNCTION
# ------------------------------------------------------------------------------

find_file <- function(candidates) {
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(paste("File not found in paths:", paste(candidates, collapse = " | ")))
  }
  found[1]
}

# ==============================================================================
# 1. FILE LOADING 
# ==============================================================================

ref_csv     <- find_file(c(here("data", "dbase", "prod_dbase_ha_ref.csv"), here("prod_dbase_ha_ref.csv")))
ncrmp_csv   <- find_file(c(here("data", "dbase", "prod_dbase_ncrmp.csv"), here("prod_dbase_ncrmp.csv")))
islands_csv <- find_file(c(here("data",  "supporting data", "islands_midpoints_all.csv"), here("islands_midpoints_all.csv")))
cau_csv     <- find_file(c(here("data", "supporting data", "HA_CAU_data.csv"), here("HA_CAU_data.csv")))

combined_ncrmp_ref_df <- read_csv(ref_csv, show_col_types = FALSE)
ncrmp_df              <- read_csv(ncrmp_csv, show_col_types = FALSE)
islands_df            <- read_csv(islands_csv, show_col_types = FALSE)
cau_df                <- read_csv(cau_csv, show_col_types = FALSE)

# ==============================================================================
# PART B: DYNAMIC TEMPLATE GENERATION & GRIGG (1982) SCALING (Corals & CCA Only)
# ==============================================================================

# 1. Base reference metrics for Corals
base_metrics <- combined_ncrmp_ref_df %>%
  select(
    SUBSTRATE_CODE,
    MORPHOLOGYCODE,
    LATITUDE_REF, 
    BASE_EXTENSION_CM_YR    = EXTENSION_CM_YR,
    BASE_EXTENSION_CM_YR_CI = EXTENSION_CM_YR_CI,
    BASE_DENSITY_G_CM3      = DENSITY_G_CM3,
    BASE_DENSITY_G_CM3_CI   = DENSITY_G_CM3_CI,
    CONVERSION_FACTOR
  )

# Identify target islands where NCRMP == TRUE
target_lats <- islands_df %>%
  filter(NCRMP == TRUE) %>%
  select(LOCATIONCODE = ISLANDCODE, TARGET_LATITUDE = LATITUDE) %>%
  distinct(LOCATIONCODE, .keep_all = TRUE)

# 2. Base taxa descriptions for Corals
base_taxa <- combined_ncrmp_ref_df %>%
  select(SUBSTRATE_CLASS, CORAL_GROUP_NAME, CORAL_GROUP, SUBSTRATE_NAME,
         SUBSTRATE_CODE, MORPHOLOGY, MORPHOLOGYCODE, TAXA_LEVEL) %>%
  distinct()

# Explicitly add a CCA placeholder row so crossing() generates it for every island
base_taxa <- base_taxa %>%
  bind_rows(
    tibble(
      SUBSTRATE_CLASS  = "CCA",
      CORAL_GROUP_NAME = "Crustose Coralline Algae",
      CORAL_GROUP      = "CCA",
      SUBSTRATE_NAME   = "Crustose Coralline Algae",
      SUBSTRATE_CODE   = "CCA",
      MORPHOLOGY       = NA_character_,
      MORPHOLOGYCODE   = NA_character_,
      TAXA_LEVEL       = NA_character_ 
    )
  ) %>%
  distinct(SUBSTRATE_CODE, MORPHOLOGYCODE, .keep_all = TRUE)

# Dynamically cross islands with corals + CCA, then apply scaling
populated_df <- crossing(LOCATIONCODE = target_lats$LOCATIONCODE, base_taxa) %>%
  left_join(target_lats, by = "LOCATIONCODE") %>%
  left_join(base_metrics, by = c("SUBSTRATE_CODE", "MORPHOLOGYCODE")) %>%
  mutate(
    # If no empirical reference latitude exists (NA), default scale factor to 1 to preserve baseline rate
    ExtensionScaleFactor = if_else(is.na(LATITUDE_REF), 1, (3.35 - 0.108 * TARGET_LATITUDE) / (3.35 - 0.108 * LATITUDE_REF)),
    DensityScaleFactor   = if_else(is.na(LATITUDE_REF), 1, (0.67 + 0.036 * TARGET_LATITUDE) / (0.67 + 0.036 * LATITUDE_REF)),
    
    EST_EXTENSION = BASE_EXTENSION_CM_YR * ExtensionScaleFactor,
    EST_DENSITY   = BASE_DENSITY_G_CM3 * DensityScaleFactor,
    
    EXTENSION_CM_YR = if_else(EST_EXTENSION < 0, 0, EST_EXTENSION),
    DENSITY_G_CM3   = EST_DENSITY,
    
    # Pass CIs and reference latitude directly
    EXTENSION_CM_YR_CI = BASE_EXTENSION_CM_YR_CI,
    DENSITY_G_CM3_CI   = BASE_DENSITY_G_CM3_CI
  ) %>%
  # Clean up intermediate generation columns
  select(
    SUBSTRATE_CLASS, CORAL_GROUP_NAME, CORAL_GROUP, SUBSTRATE_NAME,
    SUBSTRATE_CODE, MORPHOLOGY, MORPHOLOGYCODE, TAXA_LEVEL,
    EXTENSION_CM_YR, EXTENSION_CM_YR_CI, DENSITY_G_CM3, DENSITY_G_CM3_CI,
    CONVERSION_FACTOR, LOCATIONCODE, TARGET_LATITUDE, LATITUDE_REF
  )

# ==============================================================================
# PART C: INTEGRATE CRUSTOSE CORALLINE ALGAE (CCA) DATA
# ==============================================================================

# Calculate mean calcification rate and standard error per island for CCA
cca_summary <- cau_df %>%
  filter(!is.na(CALCIFICATION_RATE_BY_UNIT)) %>%
  group_by(LOCATIONCODE) %>%
  summarize(
    CCA_EXT_MEAN = mean(CALCIFICATION_RATE_BY_UNIT, na.rm = TRUE),
    CCA_EXT_SE   = if_else(n() > 1, sd(CALCIFICATION_RATE_BY_UNIT, na.rm = TRUE) / sqrt(n()), NA_real_),
    .groups      = "drop"
  )

# Impute KAH CCA rates using the average of MOL, MAI, and LAN
kah_cca_impute <- cca_summary %>%
  filter(LOCATIONCODE %in% c("MOL", "MAI", "LAN")) %>%
  summarize(
    LOCATIONCODE = "KAH",
    CCA_EXT_MEAN = mean(CCA_EXT_MEAN, na.rm = TRUE),
    CCA_EXT_SE   = mean(CCA_EXT_SE, na.rm = TRUE)
  )

# Bind imputed KAH rates into the CCA summary
cca_summary <- bind_rows(
  cca_summary %>% filter(LOCATIONCODE != "KAH"),
  kah_cca_impute
)

# Inject CCA calculations into the populated template and set density to 1
populated_df <- populated_df %>%
  left_join(cca_summary, by = "LOCATIONCODE") %>%
  mutate(
    EXTENSION_CM_YR    = if_else(SUBSTRATE_CODE == "CCA", CCA_EXT_MEAN, EXTENSION_CM_YR),
    EXTENSION_CM_YR_CI = if_else(SUBSTRATE_CODE == "CCA", CCA_EXT_SE, EXTENSION_CM_YR_CI),
    DENSITY_G_CM3      = if_else(SUBSTRATE_CODE == "CCA", 1, DENSITY_G_CM3),
    # Set CCA reference latitude to the target island's own latitude
    LATITUDE_REF       = if_else(SUBSTRATE_CODE == "CCA", TARGET_LATITUDE, LATITUDE_REF) 
  ) %>%
  select(
    SUBSTRATE_CLASS, CORAL_GROUP_NAME, CORAL_GROUP, SUBSTRATE_NAME,
    SUBSTRATE_CODE, MORPHOLOGY, MORPHOLOGYCODE, TAXA_LEVEL,
    EXTENSION_CM_YR, EXTENSION_CM_YR_CI, DENSITY_G_CM3, DENSITY_G_CM3_CI,
    CONVERSION_FACTOR, LOCATIONCODE, LATITUDE = TARGET_LATITUDE, LATITUDE_REF
  )

# ==============================================================================
# PART D: APPEND NON-CORAL / NON-CCA SUBSTRATES (Static rows)
# ==============================================================================

other_substrates <- ncrmp_df %>%
  filter(toupper(SUBSTRATE_CLASS) != "CORAL" & SUBSTRATE_CODE != "CCA") %>%
  select(
    SUBSTRATE_CLASS, CORAL_GROUP_NAME, CORAL_GROUP, SUBSTRATE_NAME,
    SUBSTRATE_CODE, MORPHOLOGYCODE,
    DENSITY_G_CM3, DENSITY_G_CM3_CI, CONVERSION_FACTOR
  ) %>%
  distinct(SUBSTRATE_CODE, MORPHOLOGYCODE, .keep_all = TRUE) %>%
  mutate(
    MORPHOLOGY         = NA_character_,
    TAXA_LEVEL         = NA_character_,
    LOCATIONCODE       = NA_character_,
    LATITUDE           = NA_real_,
    LATITUDE_REF       = NA_real_,
    EXTENSION_CM_YR    = NA_real_,
    EXTENSION_CM_YR_CI = NA_real_
  )

# Specific columns to round (protecting coordinates from precision loss)
rate_cols <- c("EXTENSION_CM_YR", "EXTENSION_CM_YR_CI", "DENSITY_G_CM3", "DENSITY_G_CM3_CI", "CONVERSION_FACTOR")

# Bind the static rows to the bottom of the populated dataframe
populated_df <- bind_rows(populated_df, other_substrates) %>%
  mutate(across(any_of(rate_cols), ~ round(., 3)))

# ==============================================================================
# 5. EXPORT DATASETS
# ==============================================================================

out_dir <- here("data")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1. Export Original Populated Dataset ---
write_csv(populated_df, file.path(out_dir, "dbase", "prod_dbase_ha_fixed.csv"), na = "")

# --- 2. Export Populated Dataset (Coral Genus, Morphology, and CCA only) ---
populated_df_nospecies <- populated_df %>%
  filter(
    ((toupper(SUBSTRATE_CLASS) == "CORAL" & TAXA_LEVEL %in% c("GENUS", "MORPHOLOGY")) | SUBSTRATE_CODE == "CCA") &
      !(SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE %in% c("LC-CO", "LC-LA"))
  )

write_csv(populated_df_nospecies, file.path(out_dir, "dbase", "prod_dbase_ha_strs.csv"), na = "")