# ==============================================================================
# Script:  FigureS4_Sensitivity.R
# Purpose: Sensitivity Analysis Comparing SW Fore-Reef 90th Percentile vs. 
#          All StRS Habitats Across the Hawaiian Archipelago
# Inputs:  data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv 
#            (or HA_2024_StRS_Net_SiteLevel_GriggSubset.csv)
# Outputs: outputs/supplementary figures/FigureS4_Sensitivity.png
#          outputs/supplementary figures/FigureS4_Sensitivity.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggpubr)
})

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# Color Palette matching Figure S1 
# ------------------------------------------------------------------------------
color_s1_orange <- "#FF8400"  # SW mid-depth forereef (Grigg-equivalent) subset
color_s1_blue   <- "#1D70B8"  # All other survey stations in 2024

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
    "Gardner"               ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "Gardner Pinnacles"     ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "Gardner_Pinnacles"     ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "gardner"               ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "gardner pinnacles"     ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "gardner_pinnacles"     ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "Pūhāhonu"              ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "P\u016bh\u0101honu"     ~ "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
    "FFS"                   ~ "Lalo",
    "French Frigate Shoals" ~ "Lalo",
    "French Frigate"        ~ "Lalo",
    "Lisianski"             ~ "Kapou",
    "Pearl & Hermes"        ~ "Manawai",
    "Pearl and Hermes"      ~ "Manawai",
    "Kure"                  ~ "H\u014dlanik\u016b",
    "Midway"                ~ "Kuaihelani",
    "Midway Atoll"          ~ "Kuaihelani",
    "Laysan"                ~ "Kamole",
    "Laysan Island"         ~ "Kamole",
    "Maro"                  ~ "Kamokuokamohoali'i",
    "Maro Reef"             ~ "Kamokuokamohoali'i",
    "Necker"                ~ "Mokumanamana",
    "Necker Island"         ~ "Mokumanamana",
    "Nihoa"                 ~ "Nihoa",
    .default = as.character(x)
  )
}

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

# Non-parametric bootstrap standard error for the 90th percentile
boot_q90_se <- function(x, B = 1000) {
  x <- na.omit(x)
  n <- length(x)
  if (n <= 1) return(0)
  samples <- matrix(sample(x, size = n * B, replace = TRUE), nrow = n, ncol = B)
  sd(apply(samples, 2, quantile, probs = 0.90, names = FALSE))
}

# ------------------------------------------------------------------------------
# 1. LOAD & HARMONIZE DATASETS
# ------------------------------------------------------------------------------
df_strs_darwin <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>%
  mutate(
    LOCATION            = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                   if ("ISLAND" %in% names(.)) ISLAND else NULL),
    CORAL_COVER_PERC    = coalesce(if ("CORAL_COVER_PERC" %in% names(.)) CORAL_COVER_PERC else NULL,
                                   if ("HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) HARD_CORAL_COVER_PCT_MEAN else NULL),
    CORAL_COVER_SE      = coalesce(if ("CORAL_COVER_SE" %in% names(.)) CORAL_COVER_SE else NULL,
                                   if ("se.HARD_CORAL_COVER_PCT_MEAN" %in% names(.)) se.HARD_CORAL_COVER_PCT_MEAN else NULL,
                                   if ("se.CORAL_COVER_PERC" %in% names(.)) se.CORAL_COVER_PERC else NULL,
                                   0),
    GROSS_PROD_KG_M2_YR = coalesce(if ("GROSS_PROD_KG_M2_YR" %in% names(.)) GROSS_PROD_KG_M2_YR else NULL,
                                   if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL),
    GROSS_PROD_SE       = coalesce(if ("GROSS_PROD_SE" %in% names(.)) GROSS_PROD_SE else NULL,
                                   if ("se.GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) se.GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL,
                                   if ("se.GROSS_PROD_KG_M2_YR" %in% names(.)) se.GROSS_PROD_KG_M2_YR else NULL,
                                   0)
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

prep_sens_long <- function(sw_data, full_data, var_name, se_name) {
  full_domain <- full_data %>%
    mutate(PLOT_LABEL = format_islands(LOCATION)) %>%
    select(PLOT_LABEL, `All Sites` = .data[[var_name]], `All Sites SE` = .data[[se_name]])
  
  df_wide <- sw_data %>%
    mutate(PLOT_LABEL = format_islands(LOCATION)) %>%
    group_by(PLOT_LABEL) %>%
    summarize(
      LATITUDE     = mean(LATITUDE, na.rm = TRUE),
      `SW 90th`    = quantile(.data[[var_name]], probs = 0.90, na.rm = TRUE),
      `SW 90th SE` = boot_q90_se(.data[[var_name]]),
      .groups      = "drop"
    ) %>%
    left_join(full_domain, by = "PLOT_LABEL") %>%
    arrange(desc(LATITUDE)) %>%  
    mutate(PLOT_LABEL = factor(PLOT_LABEL, levels = PLOT_LABEL))
  bind_rows(
    df_wide %>% transmute(PLOT_LABEL, LATITUDE, Group = "SW 90th", Value = `SW 90th`, SE = `SW 90th SE`),
    df_wide %>% transmute(PLOT_LABEL, LATITUDE, Group = "All Sites", Value = `All Sites`, SE = `All Sites SE`)
  ) %>%
    mutate(
      Group      = factor(Group, levels = c("SW 90th", "All Sites")),
      PLOT_LABEL = factor(PLOT_LABEL, levels = levels(df_wide$PLOT_LABEL))
    )
}

set.seed(42)
coral_stats_long <- prep_sens_long(df_strs_sw, df_strs_darwin, "CORAL_COVER_PERC", "CORAL_COVER_SE")
prod_stats_long  <- prep_sens_long(df_strs_sw, df_strs_darwin, "GROSS_PROD_KG_M2_YR", "GROSS_PROD_SE")

# ------------------------------------------------------------------------------
# 2. GROUPED BAR PLOT BUILDER 
# ------------------------------------------------------------------------------
build_bar_plot <- function(df_long, y_lab, y_lim, y_brk) {
  ggplot(df_long, aes(x = PLOT_LABEL, y = Value, fill = Group)) +
    geom_col(
      position  = position_dodge(width = 0.75), 
      width     = 0.65, 
      color     = "black", 
      linewidth = 0.35, 
      alpha     = 0.90
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(0, Value - SE), 
        ymax = pmin(y_lim[2], Value + SE)
      ),
      position  = position_dodge(width = 0.75),
      width     = 0.25,
      linewidth = 0.45,
      color     = "black"
    ) +
    scale_y_continuous(
      limits = y_lim, 
      breaks = y_brk, 
      expand = expansion(mult = c(0, 0.05))
    ) +
    scale_fill_manual(
      values = c("SW 90th" = color_s1_orange, "All Sites" = color_s1_blue),
      labels = c("SW 90th" = "SW Forereef 90th Percentile", "All Sites" = "All Sites")
    ) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray75"),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.text.y        = element_text(color = "black", size = 8.5),
      axis.text.x        = element_text(color = "black", size = 8.5, angle = 45, hjust = 1, vjust = 1),
      axis.title.y       = element_text(color = "black", size = 9, face = "bold", margin = margin(r = 6)),
      axis.line          = element_blank(),
      axis.ticks         = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length  = unit(0.12, "cm"),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 8.5, face = "bold"),
      legend.position    = "top",
      plot.margin        = margin(4, 12, 10, 12, "pt")
    )
}

lbl_cov   <- "Coral cover (%)"
lbl_gprod <- expression(bold(G[gross] ~ "(kg" ~ m^{-2} ~ yr^{-1} * ")"))

p_bar_cov  <- build_bar_plot(coral_stats_long, lbl_cov, c(0, 100), seq(0, 100, 20))
p_bar_prod <- build_bar_plot(prod_stats_long, lbl_gprod, c(0, 20), seq(0, 20, 5))

# ------------------------------------------------------------------------------
# 3. FIGURE ASSEMBLY & EXPORT 
# ------------------------------------------------------------------------------
figS4_final <- ggarrange(
  p_bar_cov,
  p_bar_prod,
  ncol          = 1,
  nrow          = 2,
  align         = "v",
  labels        = "AUTO",
  font.label    = list(size = 11, face = "bold"),
  common.legend = TRUE,
  legend        = "top"
)

ggsave(file.path(supp_dir, "FigureS4_Sensitivity.png"), figS4_final, width = 17.8, height = 20.0, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(supp_dir, "FigureS4_Sensitivity.tiff"), figS4_final, width = 17.8, height = 20.0, units = "cm", dpi = 600, compression = "lzw", bg = "white")