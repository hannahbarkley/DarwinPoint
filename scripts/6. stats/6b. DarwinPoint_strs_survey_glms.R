# ==============================================================================
# Script:  6b. DarwinPoint_strs_survey_glms.R
# Purpose: Hawaiian Archipelago - StRS Carbonate Budget Latitudinal Survey GLMs
#          Fits, prints, and contrasts survey-weighted GLMs across model hierarchies:
#            1. Full Interaction:  y ~ LOCATION_LAT_c * REGIONCODE
#            2. Additive Parallel: y ~ LOCATION_LAT_c + REGIONCODE
#            3. Single Continuous: y ~ LOCATION_LAT_c
# Inputs:  - data/strs data/processed data/HA_2024_StRS_SurveyDesign.rds
#            (or HA_2024_StRS_SiteLevel_Weighted.csv)
# Outputs: - outputs/stats/HA_2024_StRS_GLM_Model_Summary.csv
#          - outputs/stats/HA_2024_StRS_GLM_Wald_Tests.csv
#          - outputs/stats/HA_2024_StRS_GLM_Manuscript_Text.txt
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(survey)
  library(crayon)
})

options(survey.lonely.psu = "adjust")
options(crayon.enabled = TRUE)

# ==============================================================================
# 1. PARAMETERS & IMPORT SURVEY DESIGN
# ==============================================================================

region <- "HA"
year   <- "2024"

proc_dir      <- here("data", "strs data", "processed data")
stats_out_dir <- here("outputs", "stats")
if (!dir.exists(stats_out_dir)) dir.create(stats_out_dir, recursive = TRUE)

des_rds_file  <- file.path(proc_dir, sprintf("%s_%s_StRS_SurveyDesign.rds", region, year))
site_csv_file <- file.path(proc_dir, sprintf("%s_%s_StRS_SiteLevel_Weighted.csv", region, year))
out_coef_csv  <- file.path(stats_out_dir, sprintf("%s_%s_StRS_GLM_Model_Summary.csv", region, year))
out_wald_csv  <- file.path(stats_out_dir, sprintf("%s_%s_StRS_GLM_Wald_Tests.csv", region, year))
out_text_txt  <- file.path(stats_out_dir, sprintf("%s_%s_StRS_GLM_Manuscript_Text.txt", region, year))

cat(cyan$bold(sprintf("\n🚀 RUNNING SURVEY GLM PIPELINE: %s %s\n", region, year)))

if (file.exists(des_rds_file)) {
  cat(blue$bold(sprintf("📂 Loading survey design: %s\n", basename(des_rds_file))))
  des_full <- readRDS(des_rds_file)
} else if (file.exists(site_csv_file)) {
  cat(yellow$bold("ℹ️ .rds not found; rebuilding design from weighted CSV...\n"))
  site_sw <- read_csv(site_csv_file, show_col_types = FALSE) %>%
    mutate(across(any_of(c("OBS_YEAR", "LOCATION", "ISLAND", "REGIONCODE")), as.factor))
  des_full <- svydesign(id = ~1, strata = ~Strat_conc, weights = ~sw, data = site_sw)
} else {
  stop("Neither SurveyDesign.rds nor SiteLevel_Weighted.csv were found. Run 3c first.")
}

# Center latitude to archipelago mean
mean_lat <- mean(des_full$variables$LOCATION_LAT, na.rm = TRUE)
des_full <- update(des_full, LOCATION_LAT_c = LOCATION_LAT - mean_lat)

responses <- c(
  "Gross Production"   = "GROSS_CARB_PROD_KG_M2_YR_MEAN",
  "Net Production"     = "NET_CARB_PROD_KG_M2_YR_MEAN",
  "Gross Bioerosion"   = "GROSS_EROSION_KG_M2_YR_MEAN",
  "Coral Cover (%)"    = "HARD_CORAL_COVER_PCT_MEAN",
  "Vertical Accretion" = "VERTICAL_ACCRETION_MM_YR"
)

# ==============================================================================
# 2. HIERARCHICAL MODEL FITTING & CONSOLE PRINTING
# ==============================================================================

coef_records <- list()
wald_records <- list()

cat(cyan$bold(sprintf("\n==================================================\nRUNNING ARCHIPELAGO DESIGN (n = %d sites)\n==================================================\n", 
                      nrow(des_full$variables))))

for (mod_name in names(responses)) {
  y_col <- responses[[mod_name]]
  
  cat(yellow$bold(sprintf("\n\n##################################################\n### RESPONSE: %s\n##################################################\n", mod_name)))
  
  # ----------------------------------------------------------------------------
  # MODEL 1: FULL INTERACTION (y ~ lat * region)
  # ----------------------------------------------------------------------------
  cat(cyan$bold(sprintf("\n--- [Model 1] Interaction Model: %s ~ LOCATION_LAT_c * REGIONCODE ---\n", mod_name)))
  
  f_inter   <- as.formula(sprintf("%s ~ LOCATION_LAT_c * REGIONCODE", y_col))
  fit_inter <- svyglm(f_inter, design = des_full, family = gaussian())
  print(summary(fit_inter))
  
  w_mhi_lat <- regTermTest(fit_inter, ~LOCATION_LAT_c)
  w_nwhi_sh <- regTermTest(fit_inter, ~REGIONCODE)
  w_inter   <- regTermTest(fit_inter, ~LOCATION_LAT_c:REGIONCODE)
  
  # Derived linear contrast for NWHI slope
  nwhi_contrast <- svycontrast(fit_inter, c("LOCATION_LAT_c" = 1, "LOCATION_LAT_c:REGIONCODENWHI" = 1))
  nwhi_est <- as.numeric(nwhi_contrast)
  nwhi_se  <- as.numeric(SE(nwhi_contrast))
  nwhi_t   <- nwhi_est / nwhi_se
  nwhi_p   <- 2 * pt(-abs(nwhi_t), df = fit_inter$df.residual)
  
  cat(magenta(sprintf("  • Wald Test Lat x Region: F = %.3f (df = %d, %d), p = %.4e\n", 
                      as.numeric(w_inter$Ftest), as.numeric(w_inter$df), as.numeric(w_inter$ddf), as.numeric(w_inter$p))))
  cat(magenta(sprintf("  • NWHI Slope Contrast:    Estimate = %.4f, SE = %.4f, t = %.3f, p = %.4e\n", 
                      nwhi_est, nwhi_se, nwhi_t, nwhi_p)))
  
  wald_records[[length(wald_records) + 1]] <- tibble(
    Response_Model  = mod_name,
    Model_Structure = "Interaction (lat * region)",
    Hypothesis_Term = "Lat x Region Interaction (Slope Difference)",
    F_stat          = as.numeric(w_inter$Ftest),
    df_num          = as.numeric(w_inter$df),
    df_den          = as.numeric(w_inter$ddf),
    p_value         = as.numeric(w_inter$p)
  )
  
  coef_raw <- as.data.frame(summary(fit_inter)$coefficients) %>%
    rownames_to_column(var = "Term") %>%
    rename(Std_Error = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    mutate(
      Response_Model  = mod_name,
      Model_Structure = "Interaction (lat * region)",
      Deg_Freedom     = fit_inter$df.residual
    )
  
  contrast_row <- tibble(
    Term            = "LOCATION_LAT_c (NWHI Slope)",
    Estimate        = nwhi_est,
    Std_Error       = nwhi_se,
    t_value         = nwhi_t,
    p_value         = nwhi_p,
    Response_Model  = mod_name,
    Model_Structure = "Interaction (lat * region)",
    Deg_Freedom     = fit_inter$df.residual
  )
  
  coef_records[[length(coef_records) + 1]] <- bind_rows(coef_raw, contrast_row)
  
  # ----------------------------------------------------------------------------
  # MODEL 2: ADDITIVE PARALLEL (y ~ lat + region)
  # ----------------------------------------------------------------------------
  cat(cyan$bold(sprintf("\n--- [Model 2] Additive Model: %s ~ LOCATION_LAT_c + REGIONCODE ---\n", mod_name)))
  
  f_add   <- as.formula(sprintf("%s ~ LOCATION_LAT_c + REGIONCODE", y_col))
  fit_add <- svyglm(f_add, design = des_full, family = gaussian())
  print(summary(fit_add))
  
  w_add_lat <- regTermTest(fit_add, ~LOCATION_LAT_c)
  w_add_reg <- regTermTest(fit_add, ~REGIONCODE)
  
  cat(magenta(sprintf("  • Wald Test Shared Slope: F = %.3f, p = %.4e\n", as.numeric(w_add_lat$Ftest), as.numeric(w_add_lat$p))))
  cat(magenta(sprintf("  • Wald Test Region Shift: F = %.3f, p = %.4e\n", as.numeric(w_add_reg$Ftest), as.numeric(w_add_reg$p))))
  
  wald_records[[length(wald_records) + 1]] <- tibble(
    Response_Model  = mod_name,
    Model_Structure = "Additive (lat + region)",
    Hypothesis_Term = "Shared Latitudinal Slope",
    F_stat          = as.numeric(w_add_lat$Ftest),
    df_num          = as.numeric(w_add_lat$df),
    df_den          = as.numeric(w_add_lat$ddf),
    p_value         = as.numeric(w_add_lat$p)
  )
  
  wald_records[[length(wald_records) + 1]] <- tibble(
    Response_Model  = mod_name,
    Model_Structure = "Additive (lat + region)",
    Hypothesis_Term = "Regional Baseline Offset",
    F_stat          = as.numeric(w_add_reg$Ftest),
    df_num          = as.numeric(w_add_reg$df),
    df_den          = as.numeric(w_add_reg$ddf),
    p_value         = as.numeric(w_add_reg$p)
  )
  
  coef_add <- as.data.frame(summary(fit_add)$coefficients) %>%
    rownames_to_column(var = "Term") %>%
    rename(Std_Error = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    mutate(
      Response_Model  = mod_name,
      Model_Structure = "Additive (lat + region)",
      Deg_Freedom     = fit_add$df.residual
    )
  
  coef_records[[length(coef_records) + 1]] <- coef_add
  
  # ----------------------------------------------------------------------------
  # MODEL 3: CONTINUOUS ARCHIPELAGO (y ~ lat)
  # ----------------------------------------------------------------------------
  cat(cyan$bold(sprintf("\n--- [Model 3] Continuous Model: %s ~ LOCATION_LAT_c ---\n", mod_name)))
  
  f_biv   <- as.formula(sprintf("%s ~ LOCATION_LAT_c", y_col))
  fit_biv <- svyglm(f_biv, design = des_full, family = gaussian())
  print(summary(fit_biv))
  
  w_biv_lat <- regTermTest(fit_biv, ~LOCATION_LAT_c)
  
  cat(magenta(sprintf("  • Wald Test Archipelago Lat: F = %.3f, p = %.4e\n", as.numeric(w_biv_lat$Ftest), as.numeric(w_biv_lat$p))))
  
  wald_records[[length(wald_records) + 1]] <- tibble(
    Response_Model  = mod_name,
    Model_Structure = "Continuous (lat only)",
    Hypothesis_Term = "Archipelago Latitudinal Slope",
    F_stat          = as.numeric(w_biv_lat$Ftest),
    df_num          = as.numeric(w_biv_lat$df),
    df_den          = as.numeric(w_biv_lat$ddf),
    p_value         = as.numeric(w_biv_lat$p)
  )
  
  coef_biv <- as.data.frame(summary(fit_biv)$coefficients) %>%
    rownames_to_column(var = "Term") %>%
    rename(Std_Error = `Std. Error`, t_value = `t value`, p_value = `Pr(>|t|)`) %>%
    mutate(
      Response_Model  = mod_name,
      Model_Structure = "Continuous (lat only)",
      Deg_Freedom     = fit_biv$df.residual
    )
  
  coef_records[[length(coef_records) + 1]] <- coef_biv
}

# ==============================================================================
# 3. TABULATE, EXPORT & PRINT SUMMARY TABLES
# ==============================================================================

cat(blue$bold("\n📑 Merging and exporting GLM summary tables...\n"))

all_coefficients <- bind_rows(coef_records) %>%
  mutate(Lat_Center_Ref = mean_lat) %>%
  select(Response_Model, Model_Structure, Term, Estimate, Std_Error, t_value, p_value, Deg_Freedom, Lat_Center_Ref)

all_wald_tests <- bind_rows(wald_records) %>%
  select(Response_Model, Model_Structure, Hypothesis_Term, F_stat, df_num, df_den, p_value)

write_csv(all_coefficients, out_coef_csv)
write_csv(all_wald_tests, out_wald_csv)

cat(green$bold(sprintf("\n✅ GLM coefficients exported to: %s\n", out_coef_csv)))
cat(green$bold(sprintf("✅ Wald tests exported to:        %s\n", out_wald_csv)))

cat(yellow$bold("\n================================================================================"))
cat(yellow$bold("\n                       ALL OMNIBUS WALD HYPOTHESIS TESTS                        "))
cat(yellow$bold("\n================================================================================\n"))
print(as.data.frame(all_wald_tests), row.names = FALSE)

# ==============================================================================
# 4. MANUSCRIPT TEXT PHRASING (CONCISE CONSOLE OUTPUT ONLY)
# ==============================================================================

cat("\n=== MANUSCRIPT TEXT PHRASING ===\n")
cat(
  "Survey GLM model fits show interactive latitude and region effects for coral cover (A), additive latitude and region effects for gross production, net production, and vertical accretion (B, D), and a region difference only for gross erosion (C).\n\n"
)