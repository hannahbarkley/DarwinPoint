# ==============================================================================
# Script:  TableS2_GLMSummary.R
# Purpose: Publication Summary Table of Carbonate Budget Latitudinal Survey GLMs
#          (Full Archipelago Models Only)
#          Rows: Latitude, Region, and Latitude × Region (3 rows per metric)
#          Columns: Model Effect | Estimate | SE | df | p-value
#          Erases the OpenXML SEQ Table field (no "Table 1:" prefix)
#          Subscripted RAP_net and Times New Roman font styling.
# Outputs: outputs/supplementary figures/TableS2_GLMSummary.docx
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(gt)
})

# Install zip package if not available for cross-platform archive repacking
if (!requireNamespace("zip", quietly = TRUE)) {
  install.packages("zip", repos = "https://cloud.r-project.org")
}

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS: METRIC & TERM ORDERING, VECTORIZED FORMATTING
# ------------------------------------------------------------------------------

metric_order <- c(
  "Coral Cover (%)",
  "Gross Production (kg m⁻² yr⁻¹)",
  "Gross Bioerosion (kg m⁻² yr⁻¹)",
  "Net Production (kg m⁻² yr⁻¹)",
  "RAPnet (mm yr⁻¹)"
)

term_order <- c(
  "Latitude",
  "Region",
  "Latitude \u00D7 Region"
)

format_metric <- function(x) {
  case_match(
    as.character(x),
    c("Coral Cover (%)", "HARD_CORAL_COVER_PCT_MEAN", "Coral Cover", 
      "CORAL_COVER_PERC")                                                   ~ "Coral Cover (%)",
    c("Gross Production", "GROSS_CARB_PROD_KG_M2_YR_MEAN", 
      "Gross Production (Interaction)", "GROSS_PROD_KG_M2_YR")              ~ "Gross Production (kg m⁻² yr⁻¹)",
    c("Gross Bioerosion", "Gross Erosion", "GROSS_EROSION_KG_M2_YR_MEAN", 
      "Gross Bioerosion (Interaction)", "GROSS_EROSION_KG_M2_YR")           ~ "Gross Bioerosion (kg m⁻² yr⁻¹)",
    c("Net Production", "NET_CARB_PROD_KG_M2_YR_MEAN", 
      "Net Production (Interaction)", "NET_PROD_KG_M2_YR")                  ~ "Net Production (kg m⁻² yr⁻¹)",
    c("Vertical Accretion", "VERTICAL_ACCRETION_MM_YR", 
      "Vertical Accretion (Interaction)", "RAP_NET", "RAPnet", 
      "RAPNET", "RAP_NET_MM_YR", "RAPnet (mm yr⁻¹)",
      "RAP<sub>net</sub> (mm yr⁻¹)")                                        ~ "RAPnet (mm yr⁻¹)",
    .default = as.character(x)
  )
}

format_term <- function(t) {
  t_low <- tolower(as.character(t))
  case_when(
    grepl(":", t_low) | grepl("interaction", t_low)          ~ "Latitude \u00D7 Region",
    grepl("intercept", t_low) | grepl("slope", t_low)        ~ NA_character_,
    grepl("lat", t_low)                                      ~ "Latitude",
    grepl("region|nwhi", t_low)                              ~ "Region",
    TRUE                                                     ~ NA_character_
  )
}

fmt_stat <- function(val, digits = 3) {
  ifelse(is.na(val), "—", sprintf(paste0("%.", digits, "f"), val))
}

fmt_pval <- function(p) {
  case_when(
    is.na(p)  ~ "—",
    p < 0.001 ~ "< 0.001",
    TRUE      ~ sprintf("%.3f", p)
  )
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("outputs", "stats"),
    here("data", "strs data", "processed data"),
    here("outputs", "supplementary figures"),
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
# 2. DATA INGESTION & ISOLATING FULL INTERACTION MODELS ONLY
# ------------------------------------------------------------------------------

df_glm_raw <- saferead_csv(c(
  "HA_2024_StRS_GLM_Model_Summary.csv", 
  "Table1_GLM_Publication_Summary.csv"
))

col_names <- names(df_glm_raw)

resp_col  <- col_names[grepl("^(response|model|metric)", col_names, ignore.case = TRUE)][1]
if (is.na(resp_col)) resp_col <- col_names[1]

scope_col <- col_names[grepl("(scope|subset|dataset)", col_names, ignore.case = TRUE)][1]

term_col  <- col_names[grepl("^term$", col_names, ignore.case = TRUE)][1]
if (is.na(term_col)) term_col <- col_names[grepl("term|parameter|coef", col_names, ignore.case = TRUE)][1]

est_col   <- col_names[grepl("^(est|estimate)", col_names, ignore.case = TRUE)][1]
se_col    <- col_names[grepl("^(std|se|std_error|std\\.error)", col_names, ignore.case = TRUE)][1]
p_col     <- col_names[grepl("^(p|pval|p_val|p_value|p\\.value|pr)", col_names, ignore.case = TRUE)][1]
df_col    <- col_names[grepl("^(df|deg_freedom|deg\\.freedom)", col_names, ignore.case = TRUE)][1]

table_summary_df <- df_glm_raw %>%
  # Exclude sensitivity subsets (full archipelago models only)
  filter(
    if (!is.na(scope_col)) {
      !grepl("lisianski|kapou", as.character(.data[[scope_col]]), ignore.case = TRUE)
    } else {
      TRUE
    }
  ) %>%
  mutate(
    Metric_Clean = format_metric(.data[[resp_col]]),
    Term_Raw     = as.character(.data[[term_col]]),
    Term_Clean   = format_term(Term_Raw),
    Estimate_Val = as.numeric(.data[[est_col]]),
    SE_Val       = as.numeric(.data[[se_col]]),
    P_Val        = as.numeric(.data[[p_col]]),
    DF_Val       = if (!is.na(df_col)) as.character(.data[[df_col]]) else "—"
  ) %>%
  filter(Metric_Clean %in% metric_order) %>%
  # Isolate degrees of freedom for the interaction model
  group_by(Metric_Clean) %>%
  filter(
    if (any(grepl(":", Term_Raw) | grepl("interaction", Term_Raw, ignore.case = TRUE))) {
      DF_Val == DF_Val[grepl(":", Term_Raw) | grepl("interaction", Term_Raw, ignore.case = TRUE)][1]
    } else {
      DF_Val == min(DF_Val, na.rm = TRUE)
    }
  ) %>%
  ungroup() %>%
  # Retain only Latitude, Region, and Latitude × Region
  filter(!is.na(Term_Clean)) %>%
  mutate(
    Metric_Clean   = factor(Metric_Clean, levels = metric_order),
    `Model Effect` = factor(Term_Clean, levels = term_order)
  ) %>%
  arrange(Metric_Clean, `Model Effect`) %>%
  transmute(
    Metric         = as.character(Metric_Clean),
    `Model Effect` = as.character(`Model Effect`),
    Estimate       = fmt_stat(Estimate_Val, digits = 3),
    SE             = fmt_stat(SE_Val, digits = 3),
    df             = DF_Val,
    `p-value`      = fmt_pval(P_Val)
  )

# ------------------------------------------------------------------------------
# 3. GT TABLE CONSTRUCTION
# ------------------------------------------------------------------------------

table_title <- paste(
  "**Table S2.** Survey-weighted generalized linear models (svyglm) evaluating",
  "latitudinal, regional (MHI = Main Hawaiian Islands, NWHI = Northwestern",
  "Hawaiian Islands), and interaction effects on carbonate budget metrics in 2024.",
  "Latitude is centered at the archipelago mean (23.07°N); MHI is the reference region."
)

gt_journal_table <- gt(table_summary_df, groupname_col = "Metric") %>%
  tab_header(
    title = md(table_title)
  ) %>%
  cols_label(
    `Model Effect` = "Model Effect",
    Estimate       = "Estimate",
    SE             = "SE",
    df             = "df",
    `p-value`      = "p-value"
  ) %>%
  fmt_markdown(columns = everything()) %>%
  cols_align(align = "left", columns = `Model Effect`) %>%
  cols_align(align = "right", columns = c(Estimate, SE)) %>%
  cols_align(align = "center", columns = c(df, `p-value`)) %>%
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
      cells_row_groups(),
      cells_body()
    )
  ) %>%
  tab_options(
    table.font.names            = "Times New Roman",
    table.font.size             = px(10.5),
    heading.align               = "left",
    heading.border.bottom.color = "transparent",
    column_labels.font.weight   = "bold",
    row_group.font.weight       = "bold",
    data_row.padding            = px(3),
    table.border.top.color      = "black",
    table.border.top.width      = px(1.5),
    table_body.border.top.color = "black",
    table_body.border.top.width = px(1),
    table.border.bottom.color   = "black",
    table.border.bottom.width   = px(1.5),
    table_body.hlines.color     = "transparent",
    row_group.border.top.color  = "gray85"
  )

# ------------------------------------------------------------------------------
# 4. EXPORT TO DOCX & OPENXML POST-PROCESSING
# ------------------------------------------------------------------------------

target_docx <- file.path(supp_dir, "TableS2_GLMSummary.docx")
if (file.exists(target_docx)) unlink(target_docx, force = TRUE)

gtsave(gt_journal_table, filename = "TableS2_GLMSummary.docx", path = supp_dir)

# Unzip, modify XML, and repack
exdir <- tempfile("docx_tnr_")
dir.create(exdir)
utils::unzip(target_docx, exdir = exdir)

# 1. Update font defaults to Times New Roman in styles.xml
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

# 2. Modify document.xml
doc_xml <- file.path(exdir, "word", "document.xml")
if (file.exists(doc_xml)) {
  raw_doc <- paste(readLines(doc_xml, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  
  # Change caption paragraph style to Normal so Word will not auto-number it
  raw_doc <- gsub("<w:pStyle\\s+w:val=\"[^\"]*caption[^\"]*\"[^>]*/>", "<w:pStyle w:val=\"Normal\"/>", raw_doc, ignore.case = TRUE)
  
  # Eradicate the entire OpenXML sequence field block ("Table " + SEQ Table + ": ")
  raw_doc <- gsub(
    "<w:r>\\s*<w:rPr>.*?<\\/w:rPr>\\s*<w:t[^>]*>Table\\s*<\\/w:t>\\s*<\\/w:r>.*?<w:r>\\s*<w:rPr>.*?<\\/w:rPr>\\s*<w:t[^>]*>:\\s*<\\/w:t>\\s*<\\/w:r>",
    "",
    raw_doc,
    perl = TRUE
  )
  
  # Strip any residual literal table prefixes
  raw_doc <- gsub("(<w:t[^>]*>)Table\\s*1:\\s*", "\\1", raw_doc, perl = TRUE)
  raw_doc <- gsub("(<w:t[^>]*>)Table1:\\s*", "\\1", raw_doc, perl = TRUE)
  
  # Format RAPnet with explicit OpenXML subscript tags
  sub_replacement <- paste0(
    "<w:t xml:space=\"preserve\">RAP</w:t></w:r>",
    "<w:r><w:rPr><w:rFonts w:ascii=\"'Times New Roman'\" w:hAnsi=\"'Times New Roman'\"/><w:sz w:val=\"20\"/><w:vertAlign w:val=\"subscript\"/></w:rPr><w:t>net</w:t></w:r>",
    "<w:r><w:rPr><w:rFonts w:ascii=\"'Times New Roman'\" w:hAnsi=\"'Times New Roman'\"/><w:sz w:val=\"20\"/></w:rPr><w:t xml:space=\"preserve\"> (mm yr⁻¹)</w:t>"
  )
  raw_doc <- gsub(
    "<w:t[^>]*>RAP(?:&lt;sub&gt;|<sub>)?net(?:&lt;/sub&gt;|</sub>)?\\s*\\(mm yr⁻¹\\)</w:t>",
    sub_replacement,
    raw_doc,
    perl = TRUE
  )
  
  writeLines(raw_doc, doc_xml, useBytes = TRUE)
}

# 3. Clean repacking with zipr to guarantee preservation of file structure
unlink(target_docx, force = TRUE)
orig_wd <- getwd()
setwd(exdir)
all_files <- list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE)

if (requireNamespace("zip", quietly = TRUE)) {
  zip::zipr(zipfile = target_docx, files = all_files)
} else {
  utils::zip(zipfile = target_docx, files = all_files, flags = "-q -r")
}
setwd(orig_wd)
unlink(exdir, recursive = TRUE)

cat(sprintf("✅ Saved Table S2 (DOCX, Separate Estimate & SE, No 'Table 1:' Prefix) to: %s\n", target_docx))