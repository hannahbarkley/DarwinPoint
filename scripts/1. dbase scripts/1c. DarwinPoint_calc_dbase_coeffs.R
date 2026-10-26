# ==============================================================================
# Script:  DarwinPoint_calc_dbase_coeffs.R
# Purpose: Pre-calculate linear regression coefficients (slopes) and intercepts 
#          (Mean, L95%, U95%) for each genus- and morphology-level substrate code 
#          in prod_dbase_ha_strs.csv.
#          Strips LATITUDE, LATITUDE_REF, LOCATION_CODE_REF, CORAL_GROUP_NAME, 
#          CORAL_GROUP, and MORPHOLOGYCODE from the final export, cleans 
#          'sp' from SUBSTRATE_NAME, and reorders remaining columns to match 
#          the clean reference table layout.
# Inputs:  - data/dbase/prod_dbase_ha_strs.csv
# Outputs: - data/dbase/prod_dbase_ha_strs_coeffs.csv
# ==============================================================================

library(tidyverse)
library(here)

# ------------------------------------------------------------------------------
# 1. CORE CALIBRATION ENGINE
# ------------------------------------------------------------------------------

compute_calibration_curve <- function(sub_class, sub_code, morph_code,
                                      g, g_ci, d, d_ci, c_conv,
                                      lookup_df, loc_code = NULL) {
  
  empty_res <- tibble(
    COEFFICIENT_MEAN = 0.0000,
    COEFFICIENT_L95  = 0.0000,
    COEFFICIENT_U95  = 0.0000,
    INTERCEPT_MEAN   = 0.0000,
    INTERCEPT_L95    = 0.0000,
    INTERCEPT_U95    = 0.0000
  )
  
  if (is.na(sub_class) || !(sub_class %in% c("CORAL", "CCA"))) {
    return(empty_res)
  }
  
  g_ci   <- coalesce(g_ci, 0.0)
  d_ci   <- coalesce(d_ci, 0.0)
  c_conv <- coalesce(c_conv, 0.0)
  
  # A. CRUSTOSE CORALLINE ALGAE (CCA)
  if (sub_class == "CCA") {
    if (is.na(g) || is.na(d)) return(empty_res)
    return(tibble(
      COEFFICIENT_MEAN = round(g * d, 4),
      COEFFICIENT_L95  = round((g - g_ci) * d, 4),
      COEFFICIENT_U95  = round((g + g_ci) * d, 4),
      INTERCEPT_MEAN   = 0.0000,
      INTERCEPT_L95    = 0.0000,
      INTERCEPT_U95    = 0.0000
    ))
  }
  
  # B. HARD CORALS
  if (sub_class == "CORAL") {
    if (is.na(g) || is.na(d)) return(empty_res)
    
    df <- data.frame(x = seq(0, 135, 1))
    
    # Case B1: Laminar-Columnar Porites (LC)
    if (!is.na(morph_code) && morph_code == "LC") {
      get_lc_sub <- function(target_morph, field) {
        filt <- lookup_df$SUBSTRATE_CODE == sub_code & lookup_df$MORPHOLOGYCODE == target_morph
        if (!is.null(loc_code) && !is.na(loc_code) && "LOCATIONCODE" %in% names(lookup_df)) {
          filt <- filt & lookup_df$LOCATIONCODE == loc_code
        }
        val <- lookup_df[[field]][which(filt)]
        if (length(val) == 0 || is.na(val[1])) 0.0 else val[1]
      }
      
      g_la    <- get_lc_sub("LC-LA", "EXTENSION_CM_YR")
      g_la_ci <- get_lc_sub("LC-LA", "EXTENSION_CM_YR_CI")
      g_co    <- get_lc_sub("LC-CO", "EXTENSION_CM_YR")
      g_co_ci <- get_lc_sub("LC-CO", "EXTENSION_CM_YR_CI")
      d_co    <- get_lc_sub("LC-CO", "DENSITY_G_CM3")
      d_co_ci <- get_lc_sub("LC-CO", "DENSITY_G_CM3_CI")
      c_co    <- get_lc_sub("LC-CO", "CONVERSION_FACTOR")
      
      if (g_la == 0.0)    g_la    <- g
      if (g_la_ci == 0.0) g_la_ci <- g_ci
      if (g_co == 0.0)    g_co    <- g
      if (g_co_ci == 0.0) g_co_ci <- g_ci
      if (d_co == 0.0)    d_co    <- d
      if (d_co_ci == 0.0) d_co_ci <- d_ci
      if (c_co == 0.0)    c_co    <- c_conv
      
      prop_la <- 0.55; prop_co <- 0.45; h <- 2
      
      df$y <- (h * g_la * d_co + 0.1 * g_la * ((prop_la * df$x) + (h * g_la)) * d_co) + 
        ((((1 - c_co) * (prop_co * df$x)) * g_co * 0.1) + (c_co * (prop_co * df$x) * g_co)) * d_co
      
      df$y_l95 <- (h * (g_la - g_la_ci) * (d_co - d_co_ci) + 0.1 * (g_la - g_la_ci) * ((prop_la * df$x) + (h * (g_la - g_la_ci))) * (d_co - d_co_ci)) + 
        ((((1 - c_co) * (prop_co * df$x)) * (g_co - g_co_ci) * 0.1) + (c_co * (prop_co * df$x) * (g_co - g_co_ci))) * (d_co - d_co_ci)
      
      df$y_u95 <- (h * (g_la + g_la_ci) * (d_co + d_co_ci) + 0.1 * (g_la + g_la_ci) * ((prop_la * df$x) + (h * (g_la + g_la_ci))) * (d_co + d_co_ci)) + 
        ((((1 - c_co) * (prop_co * df$x)) * (g_co + g_co_ci) * 0.1) + (c_co * (prop_co * df$x) * (g_co + g_co_ci))) * (d_co + d_co_ci)
      
      # Case B2: Massive, Mounding & Free-living (MD, ML, FR)
    } else if (!is.na(morph_code) && morph_code %in% c("MD", "ML", "FR")) {
      df$y     <- d * (((((df$x * 2 / pi) / 2) + g)^2 * pi / 2) - (((df$x * 2 / pi) / 2)^2 * pi / 2))
      df$y_l95 <- (d - d_ci) * (((((df$x * 2 / pi) / 2) + (g - g_ci))^2 * pi / 2) - (((df$x * 2 / pi) / 2)^2 * pi / 2))
      df$y_u95 <- (d + d_ci) * (((((df$x * 2 / pi) / 2) + (g + g_ci))^2 * pi / 2) - (((df$x * 2 / pi) / 2)^2 * pi / 2))
      
      # Case B3: Encrusting, Foliose & Table (EM, PL, FO, TB)
    } else if (!is.na(morph_code) && morph_code %in% c("EM", "PL", "FO", "TB")) {
      h <- 2
      df$y     <- h * g * d + 0.1 * g * (df$x + (h * g)) * d
      df$y_l95 <- h * (g - g_ci) * (d - d_ci) + 0.1 * (g - g_ci) * (df$x + (h * (g - g_ci))) * (d - d_ci)
      df$y_u95 <- h * (g + g_ci) * (d + d_ci) + 0.1 * (g + g_ci) * (df$x + (h * (g + g_ci))) * (d + d_ci)
      
      # Case B4: Branching (BR, BRSL, BRFA) & Columnar
    } else {
      df$y     <- ((((1 - c_conv) * df$x) * g * 0.1) + (c_conv * df$x * g)) * d
      df$y_l95 <- ((((1 - c_conv) * df$x) * (g - g_ci) * 0.1) + (c_conv * df$x * (g - g_ci))) * (d - d_ci)
      df$y_u95 <- ((((1 - c_conv) * df$x) * (g + g_ci) * 0.1) + (c_conv * df$x * (g + g_ci))) * (d + d_ci)
    }
    
    # Anchor at origin
    df[1, ] <- 0
    
    fit     <- lm(y ~ x, data = df)
    fit_l95 <- lm(y_l95 ~ x, data = df)
    fit_u95 <- lm(y_u95 ~ x, data = df)
    
    return(tibble(
      COEFFICIENT_MEAN = round(unname(coef(fit)[2]), 4),
      COEFFICIENT_L95  = round(unname(coef(fit_l95)[2]), 4),
      COEFFICIENT_U95  = round(unname(coef(fit_u95)[2]), 4),
      INTERCEPT_MEAN   = round(unname(coef(fit)[1]), 4),
      INTERCEPT_L95    = round(unname(coef(fit_l95)[1]), 4),
      INTERCEPT_U95    = round(unname(coef(fit_u95)[1]), 4)
    ))
  }
  
  return(empty_res)
}

# ------------------------------------------------------------------------------
# 2. INGESTION, CALCULATION & FORMATTING
# ------------------------------------------------------------------------------

candidates <- c(
  here("data", "dbase", "prod_dbase_ha_strs.csv"),
  here("data", "prod_dbase_ha_strs.csv"),
  here("prod_dbase_ha_strs.csv"),
  "prod_dbase_ha_strs.csv"
)

found_path <- candidates[file.exists(candidates)]
stopifnot("Input file prod_dbase_ha_strs.csv not found!" = length(found_path) > 0)

prod_dbase_ha_strs_path <- found_path[1]
cat(sprintf("\n📂 Reading input database: %s\n", prod_dbase_ha_strs_path))

df_raw <- read_csv(prod_dbase_ha_strs_path, show_col_types = FALSE)

cat("🧮 Fitting growth curves to calculate calibration coefficients...\n")
calibrated_params <- purrr::pmap_dfr(
  list(
    sub_class  = df_raw$SUBSTRATE_CLASS,
    sub_code   = df_raw$SUBSTRATE_CODE,
    morph_code = df_raw$MORPHOLOGYCODE,
    g          = df_raw$EXTENSION_CM_YR,
    g_ci       = df_raw$EXTENSION_CM_YR_CI,
    d          = df_raw$DENSITY_G_CM3,
    d_ci       = df_raw$DENSITY_G_CM3_CI,
    c_conv     = df_raw$CONVERSION_FACTOR,
    loc_code   = if ("LOCATIONCODE" %in% names(df_raw)) df_raw$LOCATIONCODE else rep(NA_character_, nrow(df_raw))
  ),
  compute_calibration_curve,
  lookup_df = df_raw
)

df_combined <- bind_cols(df_raw, calibrated_params)

# Clean names, filter Acropora Fast Growth, standardize morphology, and remove unwanted columns
df_cleaned <- df_combined %>%
  mutate(
    # Clean trailing 'sp' or 'sp.' from generic names while retaining species binomials
    SUBSTRATE_NAME = str_trim(str_remove(SUBSTRATE_NAME, regex("\\s*\\bsp\\b\\.?$", ignore_case = TRUE))),
    # Standardize 'Branching - Slow Growth' to 'Branching'
    MORPHOLOGY     = if_else(MORPHOLOGY == "Branching - Slow Growth", "Branching", MORPHOLOGY)
  ) %>%
  filter(
    # Remove all Acropora Branching - Fast Growth rows
    !((SUBSTRATE_CODE == "ACSP" | SUBSTRATE_NAME == "Acropora") & 
        (MORPHOLOGY == "Branching - Fast Growth" | MORPHOLOGYCODE == "BRFA"))
  ) %>%
  select(-any_of(c(
    "LATITUDE", "LATITUDE_REF", "LOCATION_CODE_REF", "LOCATION_REF",
    "CORAL_GROUP_NAME", "CORAL_GROUP"
  )))

# Align strictly to the reference layout order
template_ordered_cols <- c(
  "SUBSTRATE_CODE",      # CODE
  "SUBSTRATE_NAME",      # Genera/Taxon
  "MORPHOLOGY",          # Morphology
  "MORPHOLOGYCODE",      # Morphology code
  "EXTENSION_CM_YR",     # Mean extension rate (cm/yr)
  "EXTENSION_CM_YR_CI",  # CI
  "DENSITY_G_CM3",       # Mean density (g/cm³)
  "DENSITY_G_CM3_CI",    # CI
  "CONVERSION_FACTOR",   # Conversion Factor
  "COEFFICIENT_MEAN",    # Coefficient mean
  "COEFFICIENT_L95",     # Coefficient L95%
  "COEFFICIENT_U95",     # Coefficient U95%
  "INTERCEPT_MEAN",      # Intercept mean
  "INTERCEPT_L95",       # Intercept L95%
  "INTERCEPT_U95"        # Intercept U95%
)

# Retain location and taxonomy classifications at the trailing end
auxiliary_cols <- intersect(
  c("LOCATIONCODE"),
  names(df_cleaned)
)

final_column_order <- c(intersect(template_ordered_cols, names(df_cleaned)), auxiliary_cols)

df_formatted <- df_cleaned %>%
  select(all_of(final_column_order))

# ------------------------------------------------------------------------------
# 3. EXPORT
# ------------------------------------------------------------------------------

out_dir <- here("data", "dbase")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

dest_path <- file.path(out_dir, "prod_dbase_ha_strs_coeffs.csv")
write_csv(df_formatted, dest_path, na = "")
