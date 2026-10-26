# ==============================================================================
# Script:  FigureS5_Sensitivity.R
# Purpose: Sensitivity Analysis Comparing All StRS Habitats vs. SW Fore-Reef Subset
# Inputs:  data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv 
#            (or HA_2024_StRS_Net_SiteLevel_GriggSubset.csv)
# Outputs: outputs/supplementary figures/FigureS5_Sensitivity.png
#          outputs/supplementary figures/FigureS5_Sensitivity.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggpubr)
})

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

color_acc_red    <- "#D73027"
color_acc_orange <- "#F46D43"
color_acc_lblue  <- "#74ADD1"
color_acc_blue   <- "#4575B4"

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
    "FFS"                   ~ "Lalo",
    "French Frigate Shoals" ~ "Lalo",
    "French Frigate"        ~ "Lalo",
    "Lisianski"             ~ "Kapou",
    "Pearl & Hermes"        ~ "Manawai",
    "Pearl and Hermes"      ~ "Manawai",
    "Kure"                  ~ "H\u014dlanik\u016b",
    .default = as.character(x)
  )
}

no_x_title <- function(p) p + labs(x = NULL) + theme(axis.title.x = element_blank())

saferead_csv <- function(candidates) {
  search_dirs <- c(here("data", "strs data", "processed data"), here("data", "grigg data"), here("data"), here())
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
    LOCATION            = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                   if ("ISLAND" %in% names(.)) ISLAND else NULL),
    CORAL_COVER_PERC    = coalesce(if ("CORAL_COVER_PERC" %in% names(.)) CORAL_COVER_PERC else NULL,
                                   if ("HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) HARD_CORAL_COVER_PCT_MEAN else NULL),
    GROSS_PROD_KG_M2_YR = coalesce(if ("GROSS_PROD_KG_M2_YR" %in% names(.)) GROSS_PROD_KG_M2_YR else NULL,
                                   if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL)
  )

df_strs_sw <- saferead_csv(c("HA_2024_Grigg_Equivalent_Sites.csv", "HA_2024_StRS_Net_SiteLevel_GriggSubset.csv")) %>%
  mutate(
    LOCATION            = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                   if ("ISLAND" %in% names(.)) ISLAND else NULL),
    CORAL_COVER_PERC    = coalesce(if ("CORAL_COVER_PERC" %in% names(.)) CORAL_COVER_PERC else NULL,
                                   if ("HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) HARD_CORAL_COVER_PCT_MEAN else NULL),
    GROSS_PROD_KG_M2_YR = coalesce(if ("GROSS_PROD_KG_M2_YR" %in% names(.)) GROSS_PROD_KG_M2_YR else NULL,
                                   if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL)
  )

prep_sens <- function(sw_data, full_data, var_name) {
  full_domain <- full_data %>%
    mutate(PLOT_LABEL = format_islands(LOCATION)) %>%
    select(PLOT_LABEL, `All Sites` = .data[[var_name]])
  
  sw_data %>%
    mutate(PLOT_LABEL = format_islands(LOCATION)) %>%
    group_by(PLOT_LABEL) %>%
    summarize(
      LATITUDE  = mean(LATITUDE, na.rm = TRUE),
      `SW Mean` = mean(.data[[var_name]], na.rm = TRUE),
      `SW 75th` = quantile(.data[[var_name]], probs = 0.75, na.rm = TRUE),
      `SW 90th` = quantile(.data[[var_name]], probs = 0.90, na.rm = TRUE),
      `SW Max`  = max(.data[[var_name]], na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    left_join(full_domain, by = "PLOT_LABEL") %>%
    arrange(desc(LATITUDE)) %>%
    mutate(PLOT_LABEL = factor(PLOT_LABEL, levels = PLOT_LABEL))
}

coral_stats <- prep_sens(df_strs_sw, df_strs_darwin, "CORAL_COVER_PERC")
prod_stats  <- prep_sens(df_strs_sw, df_strs_darwin, "GROSS_PROD_KG_M2_YR")

build_sens_plot <- function(df_wide, y_lab, y_lim, y_brk) {
  df_long <- df_wide %>%
    pivot_longer(
      cols      = c("All Sites", "SW Mean", "SW 75th", "SW 90th", "SW Max"),
      names_to  = "Statistic",
      values_to = "Value"
    ) %>%
    mutate(Statistic = factor(Statistic, levels = c("All Sites", "SW Mean", "SW 75th", "SW 90th", "SW Max")))
  
  ggplot() +
    geom_segment(
      data = df_wide, 
      aes(x = PLOT_LABEL, xend = PLOT_LABEL, y = pmin(`All Sites`, `SW Mean`, na.rm = TRUE), yend = `SW Max`), 
      color = "gray60", linetype = "dashed", alpha = 0.7, linewidth = 0.6
    ) +
    geom_point(
      data = df_long, 
      aes(x = PLOT_LABEL, y = Value, fill = Statistic), 
      shape = 21, color = "black", size = 3.0, stroke = 0.5, alpha = 0.85
    ) +
    scale_y_continuous(limits = y_lim, breaks = y_brk) +
    scale_fill_manual(values = c(
      "All Sites" = "#1E293B", 
      "SW Mean"   = color_acc_blue, 
      "SW 75th"   = color_acc_lblue, 
      "SW 90th"   = color_acc_orange, 
      "SW Max"    = color_acc_red
    )) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
      axis.text.y        = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x        = element_text(color = "black", size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title         = element_text(color = "black", size = 11, face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.6),
      axis.ticks         = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length  = unit(0.15, "cm"),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 10, face = "bold"),
      plot.margin        = margin(10, 14, 10, 14, "pt")
    )
}

build_dumbbell_plot <- function(df_wide, y_lab, y_lim, y_brk) {
  ggplot(df_wide) +
    geom_segment(
      aes(x = PLOT_LABEL, xend = PLOT_LABEL, y = `All Sites`, yend = `SW 90th`), 
      color = "gray50", linewidth = 0.8
    ) +
    geom_point(
      aes(x = PLOT_LABEL, y = `All Sites`, fill = "All Sites"), 
      shape = 21, color = "black", size = 3.2, stroke = 0.5, alpha = 0.90
    ) +
    geom_point(
      aes(x = PLOT_LABEL, y = `SW 90th`, fill = "SW 90th"), 
      shape = 21, color = "black", size = 3.2, stroke = 0.5, alpha = 0.90
    ) +
    scale_y_continuous(limits = y_lim, breaks = y_brk) +
    scale_fill_manual(values = c("All Sites" = "#1E293B", "SW 90th" = color_acc_orange)) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
      axis.text.y        = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x        = element_text(color = "black", size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title         = element_text(color = "black", size = 11, face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.6),
      axis.ticks         = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length  = unit(0.15, "cm"),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 10, face = "bold"),
      plot.margin        = margin(10, 14, 10, 14, "pt")
    )
}

p_sens_cov  <- build_sens_plot(coral_stats, "Coral Cover (%)", c(0, 100), seq(0, 100, 20))
p_sens_prod <- build_sens_plot(prod_stats, expression(bold("Gross Production Rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), c(0, 20), seq(0, 20, 5))

p_dumb_cov  <- build_dumbbell_plot(coral_stats, "Coral Cover (%)", c(0, 100), seq(0, 100, 20))
p_dumb_prod <- build_dumbbell_plot(prod_stats, expression(bold("Gross Production Rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), c(0, 20), seq(0, 20, 5))

# ggarrange keeps a single 5-item legend at the top and places tags without overlap
figS5_combined <- ggarrange(
  no_x_title(p_sens_cov),
  no_x_title(p_dumb_cov),
  p_sens_prod,
  p_dumb_prod,
  ncol          = 2,
  nrow          = 2,
  align         = "hv",
  labels        = "AUTO",
  font.label    = list(size = 12, face = "bold"),
  common.legend = TRUE,
  legend        = "top"
)

ggsave(file.path(supp_dir, "FigureS5_Sensitivity.png"), figS5_combined, width = 28.0, height = 20.0, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS5_Sensitivity.tiff"), figS5_combined, width = 28.0, height = 20.0, units = "cm", dpi = 300, compression = "lzw", bg = "white")

cat("✅ Figure S5 saved successfully.\n")