# ==============================================================================
# Script:  TableS2_GLMSummary.R
# Purpose: Publication Summary Table of Carbonate Budget Latitudinal Survey GLMs
#          (Full Archipelago vs. Lisianski Excluded Sensitivity Analysis)
# Outputs: outputs/main figures/TableS2.docx
#          outputs/main figures/TableS2.rtf
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(gt)
})

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS: METRIC ORDERING, FORMATTING & SAFE INGESTION
# ------------------------------------------------------------------------------

# Define strict publication metric order:
# 1. Coral Cover, 2. Gross Production, 3. Gross Erosion, 4. Net Production, 5. Vertical Accretion
metric_order <- c(
  "Coral Cover (%)",
  "Gross Production (kg m⁻² yr⁻¹)",
  "Gross Erosion (kg m⁻² yr⁻¹)",
  "Net Production (kg m⁻² yr⁻¹)",
  "Vertical Accretion (mm yr⁻¹)"
)

format_metric <- function(x) {
  case_match(
    as.character(x),
    c("Coral Cover (%)", "HARD_CORAL_COVER_PCT_MEAN", "Coral Cover")                     ~ "Coral Cover (%)",
    c("Gross Production", "GROSS_CARB_PROD_KG_M2_YR_MEAN", "Gross Production (Interaction)") ~ "Gross Production (kg m⁻² yr⁻¹)",
    c("Gross Bioerosion", "Gross Erosion", "GROSS_EROSION_KG_M2_YR_MEAN", 
      "Gross Bioerosion (Interaction)")                                                  ~ "Gross Erosion (kg m⁻² yr⁻¹)",
    c("Net Production", "NET_CARB_PROD_KG_M2_YR_MEAN", "Net Production (Interaction)")     ~ "Net Production (kg m⁻² yr⁻¹)",
    c("Vertical Accretion", "VERTICAL_ACCRETION_MM_YR", 
      "Vertical Accretion (Interaction)")                                                ~ "Vertical Accretion (mm yr⁻¹)",
    .default = as.character(x)
  )
}

fmt_val_se <- function(val, se = NULL, digits = 3, add_stars = FALSE, pval = NULL) {
  map2_chr(val, if (is.null(se)) list(NULL) else se, function(v, s) {
    if (is.na(v)) return("—")
    v_str <- sprintf(paste0("%.", digits, "f"), v)
    if (is.null(s) || is.na(s)) {
      return(v_str)
    } else {
      s_str <- sprintf(paste0("%.", digits, "f"), s)
      return(paste0(v_str, " (", s_str, ")"))
    }
  })
}

fmt_pval <- function(p) {
  map_chr(p, function(x) {
    if (is.na(x)) return("—")
    if (x < 0.001) return("< 0.001***")
    if (x < 0.01)  return(sprintf("%.3f**", x))
    if (x < 0.05)  return(sprintf("%.3f*", x))
    if (x < 0.10)  return(sprintf("%.3f†", x))
    sprintf("%.3f", x)
  })
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("outputs", "stats"),
    here("data", "strs data", "processed data"),
    here("outputs", "main figures"),
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
# 2. DATA INGESTION & RESHAPING
# ------------------------------------------------------------------------------

# Ingest output produced by Script 6b
df_glm_raw <- saferead_csv(c(
  "HA_2024_StRS_GLM_Model_Summary.csv", 
  "Table1_GLM_Publication_Summary.csv"
))

# Clean, pivot, and enforce factor levels for ordered grouping
table_summary_df <- df_glm_raw %>%
  mutate(
    Metric_Clean = format_metric(Response_Model),
    Scope_Clean  = case_match(
      as.character(Dataset_Scope),
      "All_Islands"    ~ "All Islands",
      "Excl_Lisianski" ~ "Excl. Lisianski (Kapou)",
      .default         = as.character(Dataset_Scope)
    )
  ) %>%
  # Filter to the requested metrics and pivot terms wide
  filter(Metric_Clean %in% metric_order) %>%
  pivot_wider(
    id_cols     = c(Metric_Clean, Scope_Clean, Deg_Freedom),
    names_from  = Term,
    values_from = c(Estimate, Std_Error, p_value)
  ) %>%
  mutate(
    # Net slope in NWHI = MHI slope (beta1) + interaction slope (beta3)
    NWHI_Slope = `Estimate_LOCATION_LAT_c` + `Estimate_LOCATION_LAT_c:REGIONCODENWHI`,
    # Enforce exact response variable sequence
    Metric_Clean = factor(Metric_Clean, levels = metric_order),
    Scope_Clean  = factor(Scope_Clean, levels = c("All Islands", "Excl. Lisianski (Kapou)"))
  ) %>%
  arrange(Metric_Clean, Scope_Clean) %>%
  transmute(
    Metric           = as.character(Metric_Clean),
    Scope            = as.character(Scope_Clean),
    `Intercept`      = fmt_val_se(`Estimate_(Intercept)`, `Std_Error_(Intercept)`),
    `Lat_MHI`        = fmt_val_se(`Estimate_LOCATION_LAT_c`, `Std_Error_LOCATION_LAT_c`),
    `Region_NWHI`    = fmt_val_se(`Estimate_REGIONCODENWHI`, `Std_Error_REGIONCODENWHI`),
    `Lat_x_Region`   = fmt_val_se(`Estimate_LOCATION_LAT_c:REGIONCODENWHI`, `Std_Error_LOCATION_LAT_c:REGIONCODENWHI`),
    `Int_P`          = fmt_pval(`p_value_LOCATION_LAT_c:REGIONCODENWHI`),
    `NWHI_Net_Slope` = sprintf("%.3f", NWHI_Slope),
    df               = as.character(Deg_Freedom)
  )

# ------------------------------------------------------------------------------
# 3. GT TABLE CONSTRUCTION
# ------------------------------------------------------------------------------

gt_journal_table <- gt(table_summary_df, groupname_col = "Metric") %>%
  tab_header(
    title    = md("**Table S2. Hawaiian Archipelago Latitudinal Carbonate Budget Survey GLMs**"),
    subtitle = md("*Survey-weighted linear models (svyglm) evaluating regional slope interactions across all islands and sensitivity subset excluding Lisianski (Kapou)*")
  ) %>%
  cols_label(
    Scope            = "Dataset Scope",
    `Intercept`      = html("Intercept (&beta;₀)<br><small>MHI Ref at 23.1°N (SE)</small>"),
    `Lat_MHI`        = html("Latitude (&beta;₁)<br><small>MHI Slope /°N (SE)</small>"),
    `Region_NWHI`    = html("Region (&beta;₂)<br><small>NWHI Offset (SE)</small>"),
    `Lat_x_Region`   = html("Lat &times; Region (&beta;₃)<br><small>Interaction (SE)</small>"),
    `Int_P`          = html("Interaction<br>p-value"),
    `NWHI_Net_Slope` = html("NWHI Net Slope<br><small>(&beta;₁ + &beta;₃)</small>"),
    df               = "df"
  ) %>%
  cols_align(align = "left", columns = c(Scope)) %>%
  cols_align(align = "center", columns = c(Int_P, NWHI_Net_Slope, df)) %>%
  cols_align(align = "right", columns = c(Intercept, Lat_MHI, Region_NWHI, Lat_x_Region)) %>%
  tab_footnote(
    footnote = md("Models fit using stratified random survey (StRS) weights with solitary sampling units handled via lonely PSU variance adjustments. Latitude is centered at the archipelago mean (23.07°N); MHI is the reference region. Significance codes: *p* < 0.001 (***), *p* < 0.01 (**), *p* < 0.05 (*), *p* < 0.10 (†)."),
    locations = cells_title(groups = "title")
  ) %>%
  tab_options(
    table.font.size           = px(11),
    column_labels.font.weight = "bold",
    row_group.font.weight     = "bold",
    data_row.padding          = px(3),
    table.border.top.color    = "black",
    table.border.bottom.color = "black",
    row_group.border.top.color = "gray80"
  )

# ------------------------------------------------------------------------------
# 4. EXPORT TO DOCX (WITH RTF FALLBACK)
# ------------------------------------------------------------------------------

target_docx <- file.path(main_dir, "TableS2.docx")
if (file.exists(target_docx)) unlink(target_docx, force = TRUE)

table_saved <- FALSE
tryCatch(
  {
    gtsave(gt_journal_table, filename = "TableS2.docx", path = main_dir)
    table_saved <- file.exists(target_docx) && file.info(target_docx)$size > 0
  },
  error = function(e) {
    message("Notice: Pandoc DOCX export encountered an issue; switching to native Word RTF.")
  }
)

if (!table_saved) {
  target_rtf <- file.path(main_dir, "TableS2.rtf")
  if (file.exists(target_rtf)) unlink(target_rtf, force = TRUE)
  gtsave(gt_journal_table, filename = "TableS2.rtf", path = main_dir)
  cat(sprintf("✅ Saved Table S2 (RTF) to: %s\n", target_rtf))
} else {
  cat(sprintf("✅ Saved Table S2 (DOCX) to: %s\n", target_docx))
}