# ==============================================================================
# Script:  TableS1_IslandSummary.R
# Purpose: Publication Summary Table of Carbonate Budget Metrics by Island
# Outputs: outputs/supplementary figures/TableS1_IslandSummary.docx
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(gt)
})

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# Hawaiian island and atoll names across the entire archipelago
format_islands <- function(x) {
  clean_key <- tolower(trimws(as.character(x)))
  clean_key <- gsub("_", " ", clean_key)
  
  lookup <- c(
    # Northwestern Hawaiian Islands
    "kure"                   = "H\u014dlanik\u016b",
    "kure atoll"             = "H\u014dlanik\u016b",
    "midway"                 = "Kuaihelani",
    "midway atoll"           = "Kuaihelani",
    "pearl & hermes"         = "Manawai",
    "pearl and hermes"       = "Manawai",
    "pearl & hermes atoll"   = "Manawai",
    "lisianski"              = "Kapou",
    "lisianski island"       = "Kapou",
    "laysan"                 = "Kamole",
    "laysan island"          = "Kamole",
    "maro"                   = "Kamokuokamohoali'i",
    "maro reef"              = "Kamokuokamohoali'i",
    "gardner"                = "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "gardner pinnacles"      = "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "french frigate"         = "Lalo",
    "french frigate shoals"  = "Lalo",
    "ffs"                    = "Lalo",
    "necker"                 = "Mokumanamana",
    "necker island"          = "Mokumanamana",
    "nihoa"                  = "Nihoa",
    "nihoa island"           = "Nihoa",
    
    # Main Hawaiian Islands
    "kauai"                  = "Kaua'i",
    "niihau"                 = "Ni'ihau",
    "oahu"                   = "O'ahu",
    "molokai"                = "Moloka'i",
    "lanai"                  = "L\u0101na'i",
    "maui"                   = "Maui",
    "kahoolawe"              = "Kaho'olawe",
    "hawaii"                 = "Hawai'i"
  )
  
  mapped <- lookup[clean_key]
  coalesce(unname(mapped), as.character(x))
}

# Formats mean and standard error on the same line: "mean (SE)"
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
    GROSS_EROSION_KG_M2_YR    = coalesce(if ("GROSS_EROSION_KG_M2_YR" %in% names(.)) GROSS_EROSION_KG_M2_YR else NULL,
                                         if ("GROSS_EROSION_KG_M2_YR_MEAN" %in% names(.)) GROSS_EROSION_KG_M2_YR_MEAN else NULL),
    se.GROSS_EROSION_KG_M2_YR = coalesce(if ("se.GROSS_EROSION_KG_M2_YR" %in% names(.)) se.GROSS_EROSION_KG_M2_YR else NULL,
                                         if ("se.GROSS_EROSION_KG_M2_YR_MEAN" %in% names(.)) se.GROSS_EROSION_KG_M2_YR_MEAN else NULL),
    NET_PROD_KG_M2_YR         = coalesce(if ("NET_PROD_KG_M2_YR" %in% names(.)) NET_PROD_KG_M2_YR else NULL,
                                         if ("NET_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) NET_CARB_PROD_KG_M2_YR_MEAN else NULL),
    se.NET_PROD_KG_M2_YR      = coalesce(if ("se.NET_PROD_KG_M2_YR" %in% names(.)) se.NET_PROD_KG_M2_YR else NULL,
                                         if ("se.NET_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) se.NET_CARB_PROD_KG_M2_YR_MEAN else NULL),
    RAP_NET                   = coalesce(if ("RAP_NET" %in% names(.)) RAP_NET else NULL,
                                         if ("RAP_NET_MM_YR" %in% names(.)) RAP_NET_MM_YR else NULL,
                                         if ("RAPnet" %in% names(.)) RAPnet else NULL,
                                         if ("RAPNET" %in% names(.)) RAPNET else NULL,
                                         if ("NET_RAP" %in% names(.)) NET_RAP else NULL,
                                         if ("VERTICAL_ACCRETION_MM_YR" %in% names(.)) VERTICAL_ACCRETION_MM_YR else NULL,
                                         0),
    se.RAP_NET                = coalesce(if ("se.RAP_NET" %in% names(.)) se.RAP_NET else NULL,
                                         if ("se.RAP_NET_MM_YR" %in% names(.)) se.RAP_NET_MM_YR else NULL,
                                         if ("se.RAPnet" %in% names(.)) se.RAPnet else NULL,
                                         if ("se.RAPNET" %in% names(.)) se.RAPNET else NULL,
                                         if ("se.NET_RAP" %in% names(.)) se.NET_RAP else NULL,
                                         if ("se.VERTICAL_ACCRETION_MM_YR" %in% names(.)) se.VERTICAL_ACCRETION_MM_YR else NULL,
                                         0)
  )

# Sort from highest to lowest latitude (North to South)
table_summary_df <- df_strs_darwin %>%
  arrange(desc(LATITUDE)) %>%
  transmute(
    Island             = PLOT_LABEL,
    `Latitude (°N)`    = sprintf("%.2f", LATITUDE),
    `Hard Coral (%)`   = fmt_val_se(CORAL_COVER_PERC, se.CORAL_COVER_PERC),
    `Gross Prod`       = fmt_val_se(GROSS_PROD_KG_M2_YR, se.GROSS_PROD_KG_M2_YR),
    `Gross Bioerosion` = fmt_val_se(GROSS_EROSION_KG_M2_YR, se.GROSS_EROSION_KG_M2_YR),
    `Net Prod`         = fmt_val_se(NET_PROD_KG_M2_YR, se.NET_PROD_KG_M2_YR),
    `RAPnet (mm/yr)`   = fmt_val_se(RAP_NET, se.RAP_NET)
  )

# Assemble table with RAP_net subscripted
gt_journal_table <- gt(table_summary_df) %>%
  tab_header(
    title = md("**Table S1.** Island survey-weighted means (standard error) for carbonate budget metrics calculated across stratified random survey stations in 2024.")
  ) %>%
  cols_label(
    Island             = "Island / Atoll",
    `Latitude (°N)`    = md("Midpoint Latitude<br>(°N)"),
    `Hard Coral (%)`   = md("Coral Cover<br>(%)"),
    `Gross Prod`       = md("Gross Production<br>(kg m⁻² yr⁻¹)"),
    `Gross Bioerosion` = md("Gross Bioerosion<br>(kg m⁻² yr⁻¹)"),
    `Net Prod`         = md("Net Production<br>(kg m⁻² yr⁻¹)"),
    `RAPnet (mm/yr)`   = md("RAP~net~<br>(mm yr⁻¹)")
  ) %>%
  fmt_markdown(columns = everything()) %>%
  cols_align(
    align   = "left",
    columns = Island
  ) %>%
  cols_align(
    align   = "center",
    columns = -Island
  ) %>%
  opt_table_font(
    font = list(
      "Times New Roman",
      default_fonts()
    )
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),
    locations = list(
      cells_title(),
      cells_column_labels(),
      cells_body()
    )
  ) %>%
  tab_options(
    table.font.names            = "Times New Roman",
    table.font.size             = px(10.5),
    heading.align               = "left",
    heading.border.bottom.color = "transparent",
    column_labels.font.weight   = "normal",
    data_row.padding            = px(3),
    table.border.top.color      = "black",
    table.border.top.width      = px(1.5),
    table_body.border.top.color = "black",
    table_body.border.top.width = px(1),
    table.border.bottom.color   = "black",
    table.border.bottom.width   = px(1.5),
    table_body.hlines.color     = "transparent"
  )

target_docx <- file.path(supp_dir, "TableS1_IslandSummary.docx")
if (file.exists(target_docx)) unlink(target_docx, force = TRUE)

# Export directly to DOCX
gtsave(gt_journal_table, filename = "TableS1_IslandSummary.docx", path = supp_dir, autonum = FALSE)

# Post-process Word XML to ensure Times New Roman and explicit subscript formatting for RAP_net
tryCatch({
  exdir <- tempfile("docx_unzip_")
  dir.create(exdir)
  utils::unzip(target_docx, exdir = exdir)
  
  # 1. Enforce document-level Times New Roman font defaults
  styles_xml <- file.path(exdir, "word", "styles.xml")
  if (file.exists(styles_xml)) {
    raw_xml <- paste(readLines(styles_xml, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    raw_xml <- gsub(
      "<w:rFonts\\b[^>]*/>",
      "<w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\" w:cs=\"Times New Roman\" w:eastAsia=\"Times New Roman\"/>",
      raw_xml
    )
    writeLines(raw_xml, styles_xml, useBytes = TRUE)
  }
  
  # 2. Enforce OpenXML subscript markup for RAP_net in document.xml
  doc_xml <- file.path(exdir, "word", "document.xml")
  if (file.exists(doc_xml)) {
    raw_doc <- paste(readLines(doc_xml, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    sub_replacement <- paste0(
      "<w:r><w:rPr><w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\"/></w:rPr><w:t xml:space=\"preserve\">RAP</w:t></w:r>",
      "<w:r><w:rPr><w:rFonts w:ascii=\"Times New Roman\" w:hAnsi=\"Times New Roman\"/><w:vertAlign w:val=\"subscript\"/></w:rPr><w:t>net</w:t></w:r>"
    )
    raw_doc <- gsub("<w:t[^>]*>RAPnet</w:t>", sub_replacement, raw_doc)
    raw_doc <- gsub("<w:t[^>]*>RAP~net~</w:t>", sub_replacement, raw_doc)
    writeLines(raw_doc, doc_xml, useBytes = TRUE)
  }
  
  # Re-pack the modified DOCX
  old_wd <- setwd(exdir)
  zip_files <- list.files(".", recursive = TRUE, all.files = TRUE, full.names = FALSE)
  zip_files <- zip_files[!zip_files %in% c(".", "..")]
  utils::zip(zipfile = target_docx, files = zip_files, flags = "-q")
  setwd(old_wd)
  unlink(exdir, recursive = TRUE)
}, error = function(e) NULL)

cat(sprintf("✅ Saved Table S1 (DOCX, Times New Roman, Hawaiian Names Only, RAPnet subscripted) to: %s\n", target_docx))