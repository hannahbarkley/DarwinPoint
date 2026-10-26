# ==============================================================================
# Script:  Table1_IslandSummary.R
# Purpose: Publication Summary Table of Carbonate Budget Metrics by Island
# Outputs: outputs/main figures/Table1_IslandSummary.docx
#          outputs/main figures/Table1_IslandSummary.rtf
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(gt)
})

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

format_islands <- function(x) {
  case_match(
    as.character(x),
    "Oahu"                  ~ "O'ahu",
    "Hawaii"                ~ "Hawai'i",
    "Lanai"                 ~ "L\u0101na'i",
    "Molokai"               ~ "Moloka'i",
    "Kauai"                 ~ "Kaua'i",
    "Niihau"                ~ "Ni'ihau",
    "Kahoolawe"             ~ "Kaho'olawe",
    "FFS"                   ~ "French Frigate Shoals (Lalo)",
    "French Frigate Shoals" ~ "French Frigate Shoals (Lalo)",
    "French Frigate"        ~ "French Frigate Shoals (Lalo)",
    "Lisianski"             ~ "Lisianski (Kapou)",
    "Pearl & Hermes"        ~ "Pearl & Hermes (Manawai)",
    "Pearl and Hermes"      ~ "Pearl & Hermes (Manawai)",
    "Kure"                  ~ "Kure (H\u014dlanik\u016b)",
    .default = as.character(x)
  )
}

fmt_val_se <- function(val, se = NULL, digits = 2) {
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

saferead_csv <- function(candidates) {
  search_dirs <- c(here("data", "strs data", "processed data"), here("data"), here())
  for (f in candidates) {
    for (d in search_dirs) {
      target <- file.path(d, f)
      if (file.exists(target)) return(read_csv(target, show_col_types = FALSE))
    }
  }
  stop(paste("None of the candidate files found:", paste(candidates, collapse = " | ")))
}

df_strs_darwin <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>%
  mutate(
    LOCATION                  = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                         if ("ISLAND" %in% names(.)) ISLAND else NULL),
    PLOT_LABEL                = format_islands(LOCATION),
    CORAL_COVER_PERC          = coalesce(if ("CORAL_COVER_PERC" %in% names(.)) CORAL_COVER_PERC else NULL,
                                         if ("HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) HARD_CORAL_COVER_PCT_MEAN else NULL),
    se.CORAL_COVER_PERC       = coalesce(if ("se.CORAL_COVER_PERC" %in% names(.)) se.CORAL_COVER_PERC else NULL,
                                         if ("se.HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) se.HARD_CORAL_COVER_PCT_MEAN else NULL),
    GROSS_PROD_KG_M2_YR       = coalesce(if ("GROSS_PROD_KG_M2_YR" %in% names(.)) GROSS_PROD_KG_M2_YR else NULL,
                                         if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL),
    se.GROSS_PROD_KG_M2_YR    = coalesce(if ("se.GROSS_PROD_KG_M2_YR" %in% names(.)) se.GROSS_PROD_KG_M2_YR else NULL,
                                         if ("se.GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) se.GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL),
    FISH_EROSION_KG_M2_YR_ALL = coalesce(if ("FISH_EROSION_KG_M2_YR_ALL" %in% names(.)) FISH_EROSION_KG_M2_YR_ALL else NULL,
                                         if ("FISH_EROSION_KG_M2_YR_ALL_MEAN" %in% names(.)) FISH_EROSION_KG_M2_YR_ALL_MEAN else NULL),
    se.FISH_EROSION_KG_M2_YR_ALL = coalesce(if ("se.FISH_EROSION_KG_M2_YR_ALL" %in% names(.)) se.FISH_EROSION_KG_M2_YR_ALL else NULL,
                                            if ("se.FISH_EROSION_KG_M2_YR_ALL_MEAN" %in% names(.)) se.FISH_EROSION_KG_M2_YR_ALL_MEAN else NULL),
    GROSS_EROSION_KG_M2_YR    = coalesce(if ("GROSS_EROSION_KG_M2_YR" %in% names(.)) GROSS_EROSION_KG_M2_YR else NULL,
                                         if ("GROSS_EROSION_KG_M2_YR_MEAN" %in% names(.)) GROSS_EROSION_KG_M2_YR_MEAN else NULL),
    se.GROSS_EROSION_KG_M2_YR = coalesce(if ("se.GROSS_EROSION_KG_M2_YR" %in% names(.)) se.GROSS_EROSION_KG_M2_YR else NULL,
                                         if ("se.GROSS_EROSION_KG_M2_YR_MEAN" %in% names(.)) se.GROSS_EROSION_KG_M2_YR_MEAN else NULL),
    NET_PROD_KG_M2_YR         = coalesce(if ("NET_PROD_KG_M2_YR" %in% names(.)) NET_PROD_KG_M2_YR else NULL,
                                         if ("NET_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) NET_CARB_PROD_KG_M2_YR_MEAN else NULL),
    se.NET_PROD_KG_M2_YR      = coalesce(if ("se.NET_PROD_KG_M2_YR" %in% names(.)) se.NET_PROD_KG_M2_YR else NULL,
                                         if ("se.NET_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) se.NET_CARB_PROD_KG_M2_YR_MEAN else NULL),
    VERTICAL_ACCRETION_MM_YR  = coalesce(if ("VERTICAL_ACCRETION_MM_YR" %in% names(.)) VERTICAL_ACCRETION_MM_YR else NULL, 0),
    se.VERTICAL_ACCRETION_MM_YR = coalesce(if ("se.VERTICAL_ACCRETION_MM_YR" %in% names(.)) se.VERTICAL_ACCRETION_MM_YR else NULL, 0)
  )

table_summary_df <- df_strs_darwin %>%
  arrange(LATITUDE) %>%
  transmute(
    Island             = PLOT_LABEL,
    `Latitude (°N)`    = sprintf("%.2f", LATITUDE),
    `Hard Coral (%)`   = fmt_val_se(CORAL_COVER_PERC, se.CORAL_COVER_PERC),
    `Gross Prod`       = fmt_val_se(GROSS_PROD_KG_M2_YR, se.GROSS_PROD_KG_M2_YR),
    `Fish Bioerosion`  = fmt_val_se(FISH_EROSION_KG_M2_YR_ALL, se.FISH_EROSION_KG_M2_YR_ALL),
    `Gross Bioerosion` = fmt_val_se(GROSS_EROSION_KG_M2_YR, se.GROSS_EROSION_KG_M2_YR),
    `Net Prod`         = fmt_val_se(NET_PROD_KG_M2_YR, se.NET_PROD_KG_M2_YR),
    `Accretion (mm/yr)`= fmt_val_se(VERTICAL_ACCRETION_MM_YR, se.VERTICAL_ACCRETION_MM_YR)
  )

gt_journal_table <- gt(table_summary_df) %>%
  tab_header(
    title    = md("**Table 1. Hawaiian Archipelago Carbonate Budget Summary by Island**"),
    subtitle = md("*Survey-weighted means and standard errors (SE) across stratified random survey (StRS) stations (2024)*")
  ) %>%
  cols_label(
    Island             = "Island / Atoll",
    `Latitude (°N)`    = "Latitude (°N)",
    `Hard Coral (%)`   = "Coral Cover (%)",
    `Gross Prod`       = html("Gross Production<br>(kg m⁻² yr⁻¹)"),
    `Fish Bioerosion`  = html("Fish Bioerosion<br>(kg m⁻² yr⁻¹)"),
    `Gross Bioerosion` = html("Gross Bioerosion<br>(kg m⁻² yr⁻¹)"),
    `Net Prod`         = html("Net Production<br>(kg m⁻² yr⁻¹)"),
    `Accretion (mm/yr)`= html("Vertical Accretion<br>(mm yr⁻¹)")
  ) %>%
  tab_options(
    table.font.size           = px(11),
    column_labels.font.weight = "bold",
    data_row.padding          = px(4),
    table.border.top.color    = "black",
    table.border.bottom.color = "black"
  )

target_docx <- file.path(main_dir, "Table1_IslandSummary.docx")
if (file.exists(target_docx)) unlink(target_docx, force = TRUE)

table_saved <- FALSE
tryCatch(
  {
    gtsave(gt_journal_table, filename = "Table1_IslandSummary.docx", path = main_dir)
    table_saved <- file.exists(target_docx) && file.info(target_docx)$size > 0
  },
  error = function(e) {
    message("Notice: Pandoc DOCX export encountered an issue; switching to native Word RTF.")
  }
)

if (!table_saved) {
  target_rtf <- file.path(main_dir, "Table1_IslandSummary.rtf")
  if (file.exists(target_rtf)) unlink(target_rtf, force = TRUE)
  gtsave(gt_journal_table, filename = "Table1_IslandSummary.rtf", path = main_dir)
  cat(sprintf("✅ Saved Table 1 (RTF) to: %s\n", target_rtf))
} else {
  cat(sprintf("✅ Saved Table 1 (DOCX) to: %s\n", target_docx))
}