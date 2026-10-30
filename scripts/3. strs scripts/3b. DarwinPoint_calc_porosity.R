# ==============================================================================
# Script:  3b. DarwinPoint_calc_porosity.R
# Purpose: Hawaiian Archipelago - Site-Level Framework Porosity Calculation
#          Calculates structural framework porosity for StRS sites using a
#          >50% majority dominance rule (Massive = 30%, Encrusting = 25%, 
#          Branching = 47%, Mixed/Baseline = 38.3%).
# Inputs:  - data/strs data/field data/2024_HA_StRS_Cover_Tier3.csv
#          - data/dbase/HA_Benthic_Codes.csv (or data/supporting data/)
# Outputs: - data/strs data/processed data/HA_2024_StRS_Porosity.csv
# ==============================================================================

library(tidyverse)
library(here)
library(crayon)

# ==============================================================================
# 1. PARAMETERS & DIRECTORY SETUP
# ==============================================================================

region <- "HA"
year   <- "2024"

base_dir <- here()
data_dir <- file.path(base_dir, "data")
proc_dir <- file.path(data_dir, "strs data", "processed data")

# Framework Porosity Benchmarks (Perry et al. 2025, Nature)
phi_massive    <- 0.30   # Massive framework (>50% massive)
phi_encrusting <- 0.25   # Encrusting framework (>50% encrusting)
phi_branching  <- 0.47   # Branching framework (>50% branching Porites)
phi_default    <- 0.383  # Mixed assemblage baseline (no single morph >50%)

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
    is.na(loc) | loc_clean == ""        ~ NA_character_,
    grepl("FRENCH|FFS", u_loc)          ~ "FFS",
    grepl("PEARL|PHR",  u_loc)          ~ "PHR",
    grepl("LIS",        u_loc)          ~ "LIS",
    grepl("KUR",        u_loc)          ~ "KUR",
    grepl("HAW",        u_loc)          ~ "HAW",
    grepl("MAU|MAI",    u_loc)          ~ "MAI",
    grepl("MOL",        u_loc)          ~ "MOL",
    grepl("LAN",        u_loc)          ~ "LAN",
    grepl("KAH",        u_loc)          ~ "KAH",
    grepl("OAH",        u_loc)          ~ "OAH",
    grepl("KAU",        u_loc)          ~ "KAU",
    grepl("NII",        u_loc)          ~ "NII",
    grepl("MID",        u_loc)          ~ "MID",
    grepl("MAR",        u_loc)          ~ "MAR",
    grepl("LAY",        u_loc)          ~ "LAY",
    grepl("GAR",        u_loc)          ~ "GAR",
    grepl("NEC",        u_loc)          ~ "NEC",
    grepl("NIH",        u_loc)          ~ "NIH",
    TRUE                                ~ loc_clean
  )
}

# Multi-path search for input datasets
find_input_file <- function(candidates) {
  for (f in candidates) {
    if (file.exists(f)) return(f)
  }
  return(candidates[1])
}

cover_csv <- find_input_file(c(
  file.path(data_dir, "strs data", "field data", "2024_HA_StRS_Cover_Tier3.csv"),
  file.path(data_dir, "field data", "2024_HA_StRS_Cover_Tier3.csv"),
  here("2024_HA_StRS_Cover_Tier3.csv")
))

lookup_csv <- find_input_file(c(
  file.path(data_dir, "dbase", "HA_Benthic_Codes.csv"),
  file.path(data_dir, "supporting data", "HA_Benthic_Codes.csv"),
  here("HA_Benthic_Codes.csv")
))

out_porosity_csv <- file.path(proc_dir, sprintf("%s_%s_StRS_Porosity.csv", region, year))

if (!dir.exists(proc_dir)) dir.create(proc_dir, recursive = TRUE)

stopifnot("Dataset 2024_HA_StRS_Cover_Tier3.csv not found!" = file.exists(cover_csv))
stopifnot("Lookup table HA_Benthic_Codes.csv not found!"    = file.exists(lookup_csv))

# ==============================================================================
# 2. INGESTION & TAXONOMIC LOOKUP SETUP
# ==============================================================================

cat(cyan(sprintf("\n🚀 STARTING FRAMEWORK POROSITY PIPELINE: %s %s\n", region, year)))
cat(blue(sprintf("📂 Loading benthic cover observations: %s\n", basename(cover_csv))))
cat(blue(sprintf("📂 Loading benthic taxon lookup table:  %s\n", basename(lookup_csv))))

cover_df  <- read_csv(cover_csv, show_col_types = FALSE)
lookup_df <- read_csv(lookup_csv, show_col_types = FALSE)

# Build hard-coral mapping and morphology assignment table
coral_code_map <- lookup_df %>%
  select(CODE, TIER_1, TIER_2) %>%
  distinct(CODE, .keep_all = TRUE) %>%
  bind_rows(tibble(
    CODE   = c("HAL", "HYCO"), 
    TIER_1 = c("ALGAE", "CORAL"), 
    TIER_2 = c("UPMA", "HYCO")
  )) %>%
  filter(TIER_1 == "CORAL") %>%
  mutate(
    MORPH_GROUP = case_when(
      TIER_2 == "BR"   ~ "Cover_Branching",
      TIER_2 == "MASS" ~ "Cover_Massive",
      TIER_2 == "ENC"  ~ "Cover_Encrusting",
      TRUE             ~ "Cover_Other_Coral"
    )
  )

# Harmonize metadata: LOCATION, LOCATIONCODE, and padded SITE (excluding OCC stations)
candidate_site_cols <- c("SITE", "OCC_SITEID", "REA_SITEID")
site_col <- intersect(candidate_site_cols, names(cover_df))[1]

if (is.na(site_col)) {
  stop("Input dataset lacks an identifiable site column (checked: SITE, OCC_SITEID, REA_SITEID).")
}

cover_df <- cover_df %>%
  mutate(
    SITE         = pad_strs_site(.data[[site_col]]),
    LOCATION     = if ("LOCATION" %in% names(.)) {
      if ("ISLAND" %in% names(.)) coalesce(LOCATION, ISLAND) else LOCATION
    } else if ("ISLAND" %in% names(.)) {
      ISLAND
    } else {
      NA_character_
    },
    LOCATIONCODE = if ("LOCATIONCODE" %in% names(.) && !all(is.na(LOCATIONCODE))) {
      coalesce(as.character(LOCATIONCODE), get_location_code(LOCATION))
    } else {
      get_location_code(LOCATION)
    },
    REGIONCODE   = if ("REGIONCODE" %in% names(.)) {
      REGIONCODE
    } else if ("REGION" %in% names(.)) {
      if_else(grepl("Northwestern", REGION, ignore.case = TRUE) | REGION == "NWHI", "NWHI", "MHI")
    } else {
      "MHI"
    }
  ) %>%
  filter(!grepl("^OCC", trimws(SITE)))

# ==============================================================================
# 3. CALCULATE MORPHOLOGICAL COVER & MAJORITY DOMINANCE POROSITY
# ==============================================================================

cat(blue("🧮 Calculating morphological cover proportions and majority-rule framework porosity...\n"))

meta_cols <- c(
  "REGION", "REGIONCODE", "LOCATION", "LOCATIONCODE", "SITE", 
  "SITEVISITID", "REEF_ZONE", "DEPTH_BIN", "LATITUDE", "LONGITUDE"
)
meta_present <- intersect(meta_cols, names(cover_df))

# Only pivot columns corresponding to scleractinian hard coral codes
coral_cols_present <- intersect(coral_code_map$CODE, names(cover_df))

if (length(coral_cols_present) == 0) {
  stop("None of the benthic coral codes match column names in the cover table.")
}

porosity_df <- cover_df %>%
  select(all_of(meta_present), all_of(coral_cols_present)) %>%
  pivot_longer(cols = all_of(coral_cols_present), names_to = "CODE", values_to = "COVER_VAL") %>%
  left_join(coral_code_map %>% select(CODE, MORPH_GROUP), by = "CODE") %>%
  group_by(across(all_of(meta_present)), MORPH_GROUP) %>%
  summarize(Total_Cover = sum(COVER_VAL, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = MORPH_GROUP, values_from = Total_Cover, values_fill = 0)

# Guarantee all morphology columns exist even if an assemblage is absent
expected_cols <- c("Cover_Branching", "Cover_Massive", "Cover_Encrusting", "Cover_Other_Coral")
for (col in expected_cols) {
  if (!col %in% names(porosity_df)) porosity_df[[col]] <- 0
}

porosity_df <- porosity_df %>%
  mutate(
    # Total live hard coral cover
    Total_Coral               = Cover_Branching + Cover_Massive + Cover_Encrusting + Cover_Other_Coral,
    HARD_CORAL_COVER_PCT      = Total_Coral,
    HARD_CORAL_COVER_PCT_MEAN = Total_Coral,
    
    # Relative proportion of each morphology
    Rel_BR    = if_else(Total_Coral > 0, Cover_Branching   / Total_Coral, 0),
    Rel_MASS  = if_else(Total_Coral > 0, Cover_Massive     / Total_Coral, 0),
    Rel_ENC   = if_else(Total_Coral > 0, Cover_Encrusting  / Total_Coral, 0),
    Rel_OTHER = if_else(Total_Coral > 0, Cover_Other_Coral / Total_Coral, 0),
    
    # Audit column classifying the majority morphology (>50%)
    Dominant_Morphology = case_when(
      Total_Coral == 0 ~ "None",
      Rel_MASS > 0.50  ~ "Massive",
      Rel_BR   > 0.50  ~ "Branching",
      Rel_ENC  > 0.50  ~ "Encrusting",
      TRUE             ~ "Mixed"
    ),
    
    # Majority-rule framework porosity (>50% dominance)
    Est_Porosity = case_when(
      Total_Coral == 0 ~ NA_real_,
      Rel_MASS > 0.50  ~ phi_massive,
      Rel_BR   > 0.50  ~ phi_branching,
      Rel_ENC  > 0.50  ~ phi_encrusting,
      TRUE             ~ phi_default
    )
  )

# ==============================================================================
# 4. EXPORT PROCESSED DATASET
# ==============================================================================

write_csv(porosity_df, out_porosity_csv)
cat(green(sprintf("✅ Porosity dataset successfully exported (%d sites):\n    -> %s\n", nrow(porosity_df), out_porosity_csv)))