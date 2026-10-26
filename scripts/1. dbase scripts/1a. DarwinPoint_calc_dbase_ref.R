# ==============================================================================
# Script:  DarwinPoint_calc_dbase_ref.R
# Purpose: Clean, aggregate, and impute empirical coral growth parameters
# Inputs:  - data/dbase/prod_dbase_ha_studies.csv (Empirical studies)
#          - data/dbase/prod_dbase_ncrmp.csv      (NCRMP baseline values)
#          - data/islands_midpoints_all.csv       (Island centroid coordinates)
#          - data/HA_2024_Fixed_Benthic.csv       (Survey benthic taxa/morphologies)
# Output:  - data/dbase/prod_dbase_ha_ref.csv
#            Harmonized master reference database containing coral growth rates
#            (extension rate, skeletal density, and conversion factors), proxy 
#            imputations, with LOCATION_CODE_REF and LATITUDE_REF relocated 
#            to the final columns.
# ==============================================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------

safe_mean <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) NA_real_ else mean(x)
}

first_non_na <- function(x) {
  valid <- x[!is.na(x)]
  if (length(valid) == 0) x[NA_integer_] else valid[1]
}

calc_ci <- function(val, raw_ci) {
  valid_idx <- which(!is.na(val))
  n <- length(valid_idx)
  
  if (n > 1) {
    # Sample standard error across multiple empirical studies
    sd(val[valid_idx]) / sqrt(n)
  } else if (n == 1) {
    # Single study fallback: inherit empirical study-reported SE from ha_csv
    ci_val <- raw_ci[valid_idx]
    if (length(ci_val) > 0 && !is.na(ci_val)) ci_val else NA_real_
  } else {
    NA_real_
  }
}

get_mode <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) return(NA_character_)
  names(which.max(table(x)))
}

find_file <- function(candidates) {
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop(paste("File not found in paths:", paste(candidates, collapse = " | ")))
  }
  found[1]
}

classify_taxa_level <- function(code, name = NULL) {
  known_genus_codes <- c("POCS", "LEPT", "PAVS", "CYPS", "PSSP", "FUSP", "TUBS", "COSC", "POSP", "ACSP")
  
  is_genus_code <- coalesce(str_ends(code, "SP"), FALSE) | (code %in% known_genus_codes)
  is_genus_name <- if (!is.null(name)) {
    coalesce(str_detect(name, regex("\\b(sp|spp)\\b", ignore_case = TRUE)), FALSE)
  } else {
    FALSE
  }
  
  if_else(is_genus_code | is_genus_name, "GENUS", "SPECIES")
}

# ------------------------------------------------------------------------------
# 1. FILE LOADING & LOOKUPS
# ------------------------------------------------------------------------------

ha_csv <- find_file(c(
  here("data", "dbase", "prod_dbase_ha_studies.csv"),
  here("data", "prod_dbase_ha_studies.csv"),
  here("prod_dbase_ha_studies.csv")
))

ncrmp_csv <- find_file(c(
  here("data", "dbase", "prod_dbase_ncrmp.csv"),
  here("data", "prod_dbase_ncrmp.csv"),
  here("prod_dbase_ncrmp.csv")
))

islands_csv <- find_file(c(
  here("data", "supporting data", "islands_midpoints_all.csv"),
  here("data", "islands_midpoints_all.csv"),
  here("islands_midpoints_all.csv")
))

ha_2024_csv <- find_file(c(
  here("data", "fixed data", "field data", "HA_2024_Fixed_Benthic.csv"),
  here("data", "HA_2024_Fixed_Benthic.csv"),
  here("HA_2024_Fixed_Benthic.csv")
))

ha_df      <- read_csv(ha_csv, show_col_types = FALSE)
ncrmp_df   <- read_csv(ncrmp_csv, show_col_types = FALSE)
islands_df <- read_csv(islands_csv, show_col_types = FALSE)
ha_2024_df <- read_csv(ha_2024_csv, show_col_types = FALSE)

ncrmp_lookup <- ncrmp_df %>%
  group_by(SUBSTRATE_CODE, MORPHOLOGYCODE) %>%
  summarize(
    SUBSTRATE_CLASS          = coalesce(first_non_na(SUBSTRATE_CLASS), "CORAL"),
    CORAL_GROUP_NAME         = first_non_na(CORAL_GROUP_NAME),
    CORAL_GROUP              = first_non_na(CORAL_GROUP),
    SUBSTRATE_NAME           = first_non_na(SUBSTRATE_NAME),
    MORPHOLOGY_NCRMP         = first_non_na(MORPHOLOGY),
    EXTENSION_CM_YR_NCRMP    = first_non_na(EXTENSION_CM_YR),
    EXTENSION_CM_YR_CI_NCRMP = first_non_na(EXTENSION_CM_YR_CI),
    DENSITY_G_CM3_NCRMP      = first_non_na(DENSITY_G_CM3),
    DENSITY_G_CM3_CI_NCRMP   = first_non_na(DENSITY_G_CM3_CI),
    CONVERSION_FACTOR        = first_non_na(CONVERSION_FACTOR),
    .groups                  = "drop"
  )

islands_lookup <- islands_df %>%
  select(LOCATION_CODE_REF = ISLANDCODE, LATITUDE) %>%
  distinct(LOCATION_CODE_REF, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 2. REFERENCE STUDIES AGGREGATION
# ------------------------------------------------------------------------------

# Species Aggregation
species_summary <- ha_df %>%
  filter(!is.na(EXTENSION_CM_YR) | !is.na(DENSITY_G_CM3)) %>%
  group_by(SUBSTRATE_CODE, MORPHOLOGYCODE) %>%
  summarize(
    MORPHOLOGY         = first_non_na(MORPHOLOGY),
    LOCATION_CODE_REF  = get_mode(LOCATION_CODE_REF),
    EXTENSION_CM_YR_CI = calc_ci(EXTENSION_CM_YR, EXTENSION_CM_YR_CI),
    EXTENSION_CM_YR    = safe_mean(EXTENSION_CM_YR),
    DENSITY_G_CM3_CI   = calc_ci(DENSITY_G_CM3, DENSITY_G_CM3_CI),
    DENSITY_G_CM3      = safe_mean(DENSITY_G_CM3),
    .groups            = "drop"
  ) %>%
  mutate(TAXA_LEVEL = classify_taxa_level(SUBSTRATE_CODE))

# Partition PRUS LC into sub-components
prus_split <- species_summary %>% filter(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC")
if (nrow(prus_split) > 0) {
  prus_components <- bind_rows(
    prus_split %>% mutate(MORPHOLOGYCODE = "LC-LA", MORPHOLOGY = "Laminar Columnar - Laminar component"),
    prus_split %>% mutate(MORPHOLOGYCODE = "LC-CO", MORPHOLOGY = "Laminar Columnar - Columnar component")
  )
  species_summary <- species_summary %>%
    filter(!(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC")) %>%
    bind_rows(prus_components)
}

# Genus Aggregation
genus_summary <- ha_df %>%
  filter(!is.na(GENUS_CODE)) %>%
  filter(!is.na(EXTENSION_CM_YR) | !is.na(DENSITY_G_CM3)) %>%
  group_by(GENUS_CODE, MORPHOLOGYCODE) %>%
  summarize(
    MORPHOLOGY         = first_non_na(MORPHOLOGY),
    LOCATION_CODE_REF  = get_mode(LOCATION_CODE_REF),
    EXTENSION_CM_YR_CI = calc_ci(EXTENSION_CM_YR, EXTENSION_CM_YR_CI),
    EXTENSION_CM_YR    = safe_mean(EXTENSION_CM_YR),
    DENSITY_G_CM3_CI   = calc_ci(DENSITY_G_CM3, DENSITY_G_CM3_CI),
    DENSITY_G_CM3      = safe_mean(DENSITY_G_CM3),
    .groups            = "drop"
  ) %>%
  mutate(
    SUBSTRATE_CODE = GENUS_CODE,
    TAXA_LEVEL     = "GENUS"
  ) %>%
  select(-GENUS_CODE)

# Partition POSP LC (retaining parent LC and appending sub-components)
posp_split <- genus_summary %>% filter(SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC")
if (nrow(posp_split) > 0) {
  posp_components <- bind_rows(
    posp_split %>% mutate(MORPHOLOGYCODE = "LC-LA", MORPHOLOGY = "Laminar Columnar - Laminar component"),
    posp_split %>% mutate(MORPHOLOGYCODE = "LC-CO", MORPHOLOGY = "Laminar Columnar - Columnar component")
  )
  genus_summary <- bind_rows(genus_summary, posp_components)
}

# Morphology Aggregation (BR [including BRSL, BRFA], ENC, and MASS)
# 1. Mean conversion factor across branching morphologies from NCRMP baseline
avg_br_conv <- ncrmp_lookup %>%
  filter(MORPHOLOGYCODE %in% c("BR", "BRSL", "BRFA")) %>%
  pull(CONVERSION_FACTOR) %>%
  safe_mean()

# 2. Mean density and standard error across morphologies from prod_dbase_ncrmp
ncrmp_density_summary <- ncrmp_df %>%
  filter(coalesce(toupper(SUBSTRATE_CLASS), "CORAL") == "CORAL") %>%
  filter(!is.na(DENSITY_G_CM3)) %>%
  mutate(
    GEN_MORPH = case_when(
      MORPHOLOGYCODE %in% c("BR", "BRSL", "BRFA") ~ "BR",
      MORPHOLOGYCODE == "EM"                      ~ "EM",
      MORPHOLOGYCODE == "MD"                      ~ "MD",
      TRUE                                        ~ NA_character_
    )
  ) %>%
  filter(!is.na(GEN_MORPH)) %>%
  group_by(GEN_MORPH) %>%
  summarize(
    NCRMP_DENSITY_G_CM3    = safe_mean(DENSITY_G_CM3),
    NCRMP_DENSITY_G_CM3_CI = calc_ci(DENSITY_G_CM3, DENSITY_G_CM3_CI),
    .groups                = "drop"
  )

# 3. Aggregate empirical extension rates from ha_df and merge NCRMP density baselines
morphology_summary <- ha_df %>%
  filter(!is.na(EXTENSION_CM_YR) | !is.na(DENSITY_G_CM3)) %>%
  mutate(
    GEN_MORPH = case_when(
      MORPHOLOGYCODE %in% c("BR", "BRSL", "BRFA") ~ "BR",
      MORPHOLOGYCODE == "EM"                      ~ "EM",
      MORPHOLOGYCODE == "MD"                      ~ "MD",
      TRUE                                        ~ NA_character_
    )
  ) %>%
  filter(!is.na(GEN_MORPH)) %>%
  group_by(GEN_MORPH) %>%
  summarize(
    MORPHOLOGYCODE     = first(GEN_MORPH),
    MORPHOLOGY         = case_when(
      first(GEN_MORPH) == "BR" ~ "Branching",
      first(GEN_MORPH) == "EM" ~ "Encrusting",
      first(GEN_MORPH) == "MD" ~ "Massive",
      TRUE                     ~ first_non_na(MORPHOLOGY)
    ),
    LOCATION_CODE_REF  = get_mode(LOCATION_CODE_REF),
    EXTENSION_CM_YR_CI = calc_ci(EXTENSION_CM_YR, EXTENSION_CM_YR_CI),
    EXTENSION_CM_YR    = safe_mean(EXTENSION_CM_YR),
    .groups            = "drop"
  ) %>%
  left_join(ncrmp_density_summary, by = "GEN_MORPH") %>%
  mutate(
    SUBSTRATE_CODE   = case_when(
      MORPHOLOGYCODE == "BR" ~ "BR",
      MORPHOLOGYCODE == "EM" ~ "ENC",
      MORPHOLOGYCODE == "MD" ~ "MASS"
    ),
    DENSITY_G_CM3    = NCRMP_DENSITY_G_CM3,
    DENSITY_G_CM3_CI = NCRMP_DENSITY_G_CM3_CI,
    TAXA_LEVEL       = "MORPHOLOGY"
  ) %>%
  select(
    SUBSTRATE_CODE, MORPHOLOGYCODE, MORPHOLOGY, LOCATION_CODE_REF,
    EXTENSION_CM_YR_CI, EXTENSION_CM_YR, DENSITY_G_CM3_CI, DENSITY_G_CM3,
    TAXA_LEVEL
  )

# ------------------------------------------------------------------------------
# 3. MERGE REFERENCE DATA & NCRMP ATTRIBUTES
# ------------------------------------------------------------------------------

comb_cols <- c(
  "SUBSTRATE_CLASS", "CORAL_GROUP_NAME", "CORAL_GROUP", "SUBSTRATE_NAME",
  "SUBSTRATE_CODE", "LOCATION_CODE_REF", "LATITUDE", "MORPHOLOGY", "MORPHOLOGYCODE", "TAXA_LEVEL",
  "EXTENSION_CM_YR", "EXTENSION_CM_YR_CI", "DENSITY_G_CM3",
  "DENSITY_G_CM3_CI", "CONVERSION_FACTOR"
)

# Reference values for Porites components
prus_la_den        <- ncrmp_lookup %>% filter(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC-LA") %>% pull(DENSITY_G_CM3_NCRMP) %>% first()
prus_la_den_ci     <- ncrmp_lookup %>% filter(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC-LA") %>% pull(DENSITY_G_CM3_CI_NCRMP) %>% first()
prus_co_den        <- ncrmp_lookup %>% filter(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC-CO") %>% pull(DENSITY_G_CM3_NCRMP) %>% first()
prus_co_den_ci     <- ncrmp_lookup %>% filter(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC-CO") %>% pull(DENSITY_G_CM3_CI_NCRMP) %>% first()

prus_lc_avg_den    <- safe_mean(c(prus_la_den, prus_co_den))
prus_lc_avg_den_ci <- safe_mean(c(prus_la_den_ci, prus_co_den_ci))
prus_co_conv       <- ncrmp_lookup %>% filter(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC-CO") %>% pull(CONVERSION_FACTOR) %>% first()

posp_sub_name <- ncrmp_lookup %>% 
  filter(SUBSTRATE_CODE == "POSP" & !is.na(SUBSTRATE_NAME)) %>% 
  pull(SUBSTRATE_NAME) %>% 
  first() %>% 
  coalesce("Porites sp")

combined_ncrmp_ref_df <- bind_rows(species_summary, genus_summary, morphology_summary) %>%
  distinct(SUBSTRATE_CODE, MORPHOLOGYCODE, .keep_all = TRUE) %>% 
  left_join(ncrmp_lookup, by = c("SUBSTRATE_CODE", "MORPHOLOGYCODE")) %>%
  left_join(islands_lookup, by = "LOCATION_CODE_REF") %>%
  mutate(
    CORAL_GROUP_NAME = case_when(
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE %in% c("LC", "LC-LA", "LC-CO") ~ coalesce(CORAL_GROUP_NAME, "Laminar columnar Porites"),
      SUBSTRATE_CODE == "BR"   ~ "Branching hard coral",
      SUBSTRATE_CODE == "ENC"  ~ "Encrusting hard coral",
      SUBSTRATE_CODE == "MASS" ~ "Massive hard coral",
      TRUE                     ~ CORAL_GROUP_NAME
    ),
    CORAL_GROUP = case_when(
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE %in% c("LC", "LC-LA", "LC-CO") ~ coalesce(CORAL_GROUP, "POLC"),
      SUBSTRATE_CODE %in% c("BR", "ENC", "MASS") ~ SUBSTRATE_CODE,
      TRUE                                        ~ CORAL_GROUP
    ),
    SUBSTRATE_NAME = case_when(
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE %in% c("LC", "LC-LA", "LC-CO") ~ coalesce(SUBSTRATE_NAME, posp_sub_name),
      SUBSTRATE_CODE == "BR"   ~ "Branching hard coral",
      SUBSTRATE_CODE == "ENC"  ~ "Encrusting hard coral",
      SUBSTRATE_CODE == "MASS" ~ "Massive hard coral",
      TRUE                     ~ SUBSTRATE_NAME
    ),
    DENSITY_G_CM3 = case_when(
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC-LA" ~ coalesce(DENSITY_G_CM3, prus_la_den),
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC-CO" ~ coalesce(DENSITY_G_CM3, prus_co_den),
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC"    ~ coalesce(DENSITY_G_CM3, prus_lc_avg_den),
      TRUE ~ coalesce(DENSITY_G_CM3, DENSITY_G_CM3_NCRMP)
    ),
    DENSITY_G_CM3_CI = case_when(
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC-LA" ~ coalesce(DENSITY_G_CM3_CI, prus_la_den_ci),
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC-CO" ~ coalesce(DENSITY_G_CM3_CI, prus_co_den_ci),
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE == "LC"    ~ coalesce(DENSITY_G_CM3_CI, prus_lc_avg_den_ci),
      TRUE ~ coalesce(DENSITY_G_CM3_CI, DENSITY_G_CM3_CI_NCRMP)
    ),
    CONVERSION_FACTOR = case_when(
      SUBSTRATE_CODE == "POSP" & MORPHOLOGYCODE %in% c("LC", "LC-CO") ~ coalesce(CONVERSION_FACTOR, prus_co_conv),
      SUBSTRATE_CODE == "BR"                                          ~ coalesce(CONVERSION_FACTOR, avg_br_conv),
      TRUE ~ CONVERSION_FACTOR
    ),
    SUBSTRATE_CLASS = coalesce(SUBSTRATE_CLASS, "CORAL"),
    MORPHOLOGY      = coalesce(MORPHOLOGY, MORPHOLOGY_NCRMP)
  ) %>%
  select(all_of(comb_cols))

# ------------------------------------------------------------------------------
# 4. IDENTIFY & APPEND MISSING COMBOS FROM 2024 DATA
# ------------------------------------------------------------------------------

missing_combos <- ha_2024_df %>%
  filter(!is.na(SUBSTRATE_CODE)) %>% 
  filter(toupper(SUBSTRATE_CLASS) == "CORAL") %>% 
  filter(!(SUBSTRATE_CODE == "PRUS" & MORPHOLOGYCODE == "LC")) %>% 
  distinct(SUBSTRATE_CODE, MORPHOLOGYCODE) %>%
  bind_rows(tibble(SUBSTRATE_CODE = "LESP", MORPHOLOGYCODE = "EM")) %>%
  bind_rows(tibble(SUBSTRATE_CODE = "MOSP", MORPHOLOGYCODE = "FO")) %>%
  distinct(SUBSTRATE_CODE, MORPHOLOGYCODE) %>%
  anti_join(combined_ncrmp_ref_df, by = c("SUBSTRATE_CODE", "MORPHOLOGYCODE"))

message(paste("Appending", nrow(missing_combos), "unmatched combinations from HA 2024."))

missing_metadata <- missing_combos %>%
  left_join(ncrmp_lookup, by = c("SUBSTRATE_CODE", "MORPHOLOGYCODE")) %>%
  mutate(
    TAXA_LEVEL          = classify_taxa_level(SUBSTRATE_CODE, SUBSTRATE_NAME),
    SUBSTRATE_CLASS     = coalesce(SUBSTRATE_CLASS, "CORAL"),
    LOCATION_CODE_REF   = NA_character_,
    LATITUDE            = NA_real_,
    MORPHOLOGY          = MORPHOLOGY_NCRMP,
    EXTENSION_CM_YR     = EXTENSION_CM_YR_NCRMP,
    EXTENSION_CM_YR_CI  = EXTENSION_CM_YR_CI_NCRMP,
    DENSITY_G_CM3       = DENSITY_G_CM3_NCRMP,
    DENSITY_G_CM3_CI    = DENSITY_G_CM3_CI_NCRMP
  ) %>%
  select(all_of(comb_cols))

combined_ncrmp_ref_df <- bind_rows(combined_ncrmp_ref_df, missing_metadata)

# ------------------------------------------------------------------------------
# 5. IMPUTE VIA SISTER SPECIES PROXIES
# ------------------------------------------------------------------------------

donor_map <- tribble(
  ~REC_SUBSTRATE, ~REC_MORPH, ~DONOR_SUBSTRATE, ~DONOR_MORPH,
  "ACSP",         "BRFA",     "ACSP",           "TB",
  "ACSP",         "EM",       "ACSP",           "BRSL",
  "POCS",         "EM",       "POCS",           "BR",
  "LBEW",         "EM",       "LEPT",           "EM",
  "MPAT",         "EM",       "MOSP",           "EM",
  "LESP",         "EM",       "ENC",            "EM",
  "MOSP",         "FO",       "MOSP",            "BR"
)

# Extract donor growth and coordinate values 
donor_values <- combined_ncrmp_ref_df %>%
  inner_join(donor_map, by = c("SUBSTRATE_CODE" = "DONOR_SUBSTRATE", "MORPHOLOGYCODE" = "DONOR_MORPH")) %>%
  select(
    REC_SUBSTRATE, 
    REC_MORPH,
    DONOR_EXT        = EXTENSION_CM_YR,
    DONOR_EXT_CI     = EXTENSION_CM_YR_CI,
    DONOR_DEN        = DENSITY_G_CM3,
    DONOR_DEN_CI     = DENSITY_G_CM3_CI,
    DONOR_LOC        = LOCATION_CODE_REF,
    DONOR_LAT        = LATITUDE
  ) %>%
  distinct(REC_SUBSTRATE, REC_MORPH, .keep_all = TRUE)

rate_cols <- c("EXTENSION_CM_YR", "EXTENSION_CM_YR_CI", "DENSITY_G_CM3", "DENSITY_G_CM3_CI", "CONVERSION_FACTOR")

combined_ncrmp_ref_df <- combined_ncrmp_ref_df %>%
  left_join(donor_values, by = c("SUBSTRATE_CODE" = "REC_SUBSTRATE", "MORPHOLOGYCODE" = "REC_MORPH")) %>%
  mutate(
    EXTENSION_CM_YR    = coalesce(DONOR_EXT, EXTENSION_CM_YR),
    EXTENSION_CM_YR_CI = coalesce(DONOR_EXT_CI, EXTENSION_CM_YR_CI),
    DENSITY_G_CM3      = coalesce(DONOR_DEN, DENSITY_G_CM3),
    DENSITY_G_CM3_CI   = coalesce(DONOR_DEN_CI, DENSITY_G_CM3_CI),
    LOCATION_CODE_REF  = coalesce(DONOR_LOC, LOCATION_CODE_REF),
    LATITUDE           = coalesce(DONOR_LAT, LATITUDE)
  ) %>%
  select(-starts_with("DONOR_")) %>%
  mutate(across(all_of(rate_cols), ~ round(., 3))) %>%
  rename(LATITUDE_REF = LATITUDE) %>%
  relocate(LOCATION_CODE_REF, LATITUDE_REF, .after = last_col())

# ------------------------------------------------------------------------------
# 6. EXPORT
# ------------------------------------------------------------------------------

out_dir <- here("data", "dbase")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_csv(combined_ncrmp_ref_df, file.path(out_dir, "prod_dbase_ha_ref.csv"), na = "")