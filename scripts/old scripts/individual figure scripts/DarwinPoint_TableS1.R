# ==============================================================================
# Script: Generate_Table1_Hawaiian_Darwin_Point.R
# Standalone Script: Table 1 - 1982 vs 2024 Coral Cover & Production Metrics
# Exports: outputs/Table1_Hawaiian_Darwin_Point.docx
# ==============================================================================

library(tidyverse)
library(here)
library(gt)

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1. SHARED HELPER FUNCTIONS ---
format_islands <- function(x) {
  case_match(
    x,
    "Oahu"                  ~ "O'ahu", 
    "Hawaii"                ~ "Hawai'i", 
    "Lanai"                 ~ "L\u0101na'i",
    "Molokai"               ~ "Moloka'i", 
    "Kauai"                 ~ "Kaua'i", 
    "Niihau"                ~ "Ni'ihau",
    "Kahoolawe"             ~ "Kaho'olawe", 
    "FFS"                   ~ "Lalo",
    "French Frigate Shoals" ~ "Lalo",
    "French Frigate"        ~ "Lalo",
    "Lisianski"             ~ "Kapou",
    "Pearl & Hermes"        ~ "Manawai",
    "Pearl and Hermes"      ~ "Manawai",
    "Kure"                  ~ "H\u014dlanik\u016b",
    "Midway"                ~ "Kuaihelani",
    "Laysan"                ~ "Kamole",
    "Maro Reef"             ~ "Kamokuokamohoali'i",
    "Maro"                  ~ "Kamokuokamohoali'i",
    "Gardner Pinnacles"     ~ "P\u016bh\u0101honu",
    "Gardner"               ~ "P\u016bh\u0101honu",
    "Necker"                ~ "Mokumanamana",
    .default = x
  )
}

format_islands_with_english <- function(x, html = TRUE) {
  sep <- if (html) "<br>" else " "
  case_match(
    x,
    "H\u014dlanik\u016b"        ~ paste0("H\u014dlanik\u016b", sep, "(Kure Atoll)"),
    "Kuaihelani"               ~ paste0("Kuaihelani", sep, "(Midway Atoll)"),
    "Manawai"                  ~ paste0("Manawai", sep, "(Pearl & Hermes)"),
    "Kapou"                    ~ paste0("Kapou", sep, "(Lisianski Island)"),
    "Kamole"                   ~ paste0("Kamole", sep, "(Laysan Island)"),
    "Kamokuokamohoali'i"       ~ paste0("Kamokuokamohoali'i", sep, "(Maro Reef)"),
    "P\u016bh\u0101honu"       ~ paste0("P\u016bh\u0101honu", sep, "(Gardner Pinnacles)"),
    "Lalo"                     ~ paste0("Lalo", sep, "(French Frigate Shoals)"),
    "Mokumanamana"             ~ paste0("Mokumanamana", sep, "(Necker Island)"),
    .default = x
  )
}

boot_q90_se <- function(x, B = 1000) {
  x <- na.omit(x)
  if (length(x) <= 1) return(0)
  boot_vals <- replicate(B, quantile(sample(x, replace = TRUE), probs = 0.90, na.rm = TRUE))
  sd(boot_vals)
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

saferead_csv <- function(file_name) {
  p1 <- here("outputs", file_name); p2 <- here("data", file_name); p3 <- here(file_name)
  if (file.exists(p1)) return(read_csv(p1, show_col_types = FALSE))
  if (file.exists(p2)) return(read_csv(p2, show_col_types = FALSE))
  if (file.exists(p3)) return(read_csv(p3, show_col_types = FALSE))
  stop(paste("Dataset", file_name, "not found."))
}

# --- 2. DATA INGESTION ---
df_grigg       <- saferead_csv("Grigg_1982_data.csv") %>% mutate(PLOT_LABEL = format_islands(LOCATION))
df_strs_darwin <- saferead_csv("2024_HI_StRS_Islandlevel_DarwinPoint.csv") %>% mutate(PLOT_LABEL = format_islands(as.character(ISLAND)))
df_strs_sw     <- saferead_csv("HA_2024_StRS_net_sw.csv")

if (!"CORAL_COVER_PERC" %in% names(df_strs_sw) && "TOTAL" %in% names(df_strs_sw)) {
  df_strs_sw$CORAL_COVER_PERC <- df_strs_sw$TOTAL
}

se_var_col <- if ("se.VERTICAL_ACCRETION_MM_YR" %in% names(df_strs_darwin)) {
  "se.VERTICAL_ACCRETION_MM_YR"
} else if ("VERTICAL_ACCRETION_MM_YR_SE" %in% names(df_strs_darwin)) {
  "VERTICAL_ACCRETION_MM_YR_SE"
} else {
  NULL
}

# --- 3. DATA AGGREGATION & SE COMPUTATION ---
set.seed(42)
strs_sw_p90 <- df_strs_sw %>% 
  mutate(PLOT_LABEL = format_islands(ISLAND)) %>% 
  group_by(PLOT_LABEL) %>% 
  summarize(
    LATITUDE_SW                     = mean(LATITUDE, na.rm = TRUE),
    CORAL_COVER_P90                 = quantile(CORAL_COVER_PERC, probs = 0.90, na.rm = TRUE),
    CORAL_COVER_P90_SE              = boot_q90_se(CORAL_COVER_PERC),
    GROSS_CARB_PROD_KG_M2_YR_P90     = quantile(GROSS_PROD_KG_M2_YR, probs = 0.90, na.rm = TRUE),
    GROSS_CARB_PROD_KG_M2_YR_P90_SE    = boot_q90_se(GROSS_PROD_KG_M2_YR),
    .groups = "drop"
  )

df_1982_summary <- df_grigg %>%
  select(
    PLOT_LABEL,
    LATITUDE_1982   = LATITUDE,
    Cover_1982      = CORAL_COVER_PERC,
    Gross_Prod_1982 = GROSS_PROD
  )

df_2024_means <- df_strs_darwin %>%
  mutate(
    Vert_Acc_SE = if (!is.null(se_var_col)) .data[[se_var_col]] else NA_real_
  ) %>%
  select(
    PLOT_LABEL,
    LATITUDE_2024        = LATITUDE,
    Cover_2024_All       = CORAL_COVER_PERC,
    Cover_2024_All_SE    = se.CORAL_COVER_PERC,
    Gross_Prod_2024_All  = GROSS_PROD_KG_M2_YR,
    Gross_Prod_2024_SE   = se.GROSS_PROD_KG_M2_YR,
    Net_Prod_2024_All    = NET_PROD_KG_M2_YR,
    Net_Prod_2024_SE     = se.NET_PROD_KG_M2_YR,
    Vert_Acc_2024_All    = VERTICAL_ACCRETION_MM_YR,
    Vert_Acc_2024_SE     = Vert_Acc_SE
  )

# --- 4. COMBINE, SORT, FORMAT & OMIT LATITUDE ---
summary_table_gt <- df_1982_summary %>%
  full_join(df_2024_means, by = "PLOT_LABEL") %>%
  full_join(strs_sw_p90, by = "PLOT_LABEL") %>%
  mutate(LATITUDE = coalesce(LATITUDE_2024, LATITUDE_SW, LATITUDE_1982)) %>%
  arrange(desc(LATITUDE)) %>%
  mutate(
    PLOT_LABEL_HTML        = format_islands_with_english(PLOT_LABEL, html = TRUE),
    Cover_1982_fmt          = fmt_val_se(Cover_1982),
    Gross_Prod_1982_fmt     = fmt_val_se(Gross_Prod_1982),
    Cover_2024_All_fmt      = fmt_val_se(Cover_2024_All, Cover_2024_All_SE),
    Cover_2024_SW_P90_fmt   = fmt_val_se(CORAL_COVER_P90, CORAL_COVER_P90_SE),
    Gross_Prod_2024_P90_fmt = fmt_val_se(GROSS_CARB_PROD_KG_M2_YR_P90, GROSS_CARB_PROD_KG_M2_YR_P90_SE),
    Gross_Prod_2024_All_fmt = fmt_val_se(Gross_Prod_2024_All, Gross_Prod_2024_SE),
    Net_Prod_2024_All_fmt   = fmt_val_se(Net_Prod_2024_All, Net_Prod_2024_SE),
    Vert_Acc_2024_All_fmt   = fmt_val_se(Vert_Acc_2024_All, Vert_Acc_2024_SE)
  ) %>%
  select(
    PLOT_LABEL = PLOT_LABEL_HTML,
    Cover_1982_fmt,
    Gross_Prod_1982_fmt,
    Cover_2024_All_fmt,
    Cover_2024_SW_P90_fmt,
    Gross_Prod_2024_P90_fmt,
    Gross_Prod_2024_All_fmt,
    Net_Prod_2024_All_fmt,
    Vert_Acc_2024_All_fmt
  )

# --- 5. JOURNAL TABLE FORMATTING & EXPORT ---
gt_journal_table <- summary_table_gt %>%
  gt() %>%
  fmt_markdown(columns = PLOT_LABEL) %>%
  tab_header(
    title = md("**Table 1. Historical (1982) and modern (2024) island-level coral cover, carbonate production, and vertical accretion rates across the Hawaiian Archipelago.**"),
    subtitle = "Islands are ordered north to south by latitude. Values represent Mean (SE) or P90 (SE)."
  ) %>%
  tab_spanner(
    label = "1982 Baseline (Grigg)",
    columns = c(Cover_1982_fmt, Gross_Prod_1982_fmt)
  ) %>%
  tab_spanner(
    label = "2024 Modern Survey (NCRMP StRS)",
    columns = c(Cover_2024_All_fmt, Cover_2024_SW_P90_fmt, Gross_Prod_2024_P90_fmt, Gross_Prod_2024_All_fmt, Net_Prod_2024_All_fmt, Vert_Acc_2024_All_fmt)
  ) %>%
  cols_label(
    PLOT_LABEL              = "Island / Atoll",
    Cover_1982_fmt          = "Coral Cover (%)",
    Gross_Prod_1982_fmt     = "Gross Prod. (kg m⁻² yr⁻¹)",
    Cover_2024_All_fmt      = "Mean Cover (%)",
    Cover_2024_SW_P90_fmt   = "SW P90 Cover (%)",
    Gross_Prod_2024_P90_fmt = "SW P90 Gross Prod. (kg m⁻² yr⁻¹)",
    Gross_Prod_2024_All_fmt = "Mean Gross Prod. (kg m⁻² yr⁻¹)",
    Net_Prod_2024_All_fmt   = "Mean Net Prod. (kg m⁻² yr⁻¹)",
    Vert_Acc_2024_All_fmt   = "Mean Vert. Accretion (mm yr⁻¹)"
  ) %>%
  cols_align(align = "center", columns = -PLOT_LABEL) %>%
  cols_align(align = "left", columns = PLOT_LABEL) %>%
  tab_footnote(
    footnote = "SW = Shallow Water (<10 m); P90 = 90th percentile estimate; SE = Standard Error (bootstrapped with B = 1000 for P90 metrics).",
    locations = cells_title(groups = "subtitle")
  ) %>%
  opt_table_font(font = list(google_font("Times New Roman"), default_fonts())) %>%
  tab_options(
    table.border.top.color = "black",
    table.border.top.width = px(2),
    table.border.bottom.color = "black",
    table.border.bottom.width = px(2),
    column_labels.border.top.color = "black",
    column_labels.border.bottom.color = "black",
    column_labels.border.bottom.width = px(1.5),
    table_body.hlines.color = "transparent",
    data_row.padding = px(3)
  )

# Save strictly as .docx
out_file <- file.path(out_dir, "Table1_Hawaiian_Darwin_Point.docx")
gtsave(gt_journal_table, out_file)
cat("Table 1 generated successfully and saved to:", out_file, "\n")