# ==============================================================================
# Script: FigureS3.R
# Hawaiian Archipelago - Latitudinal Sensitivity Analysis
# Aesthetic: Clean Vertical Cleveland dot plots, no error bars, top legend
# Outputs: outputs/FigureS3_Sensitivity.png & .tiff
# ==============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(ggpubr)

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. HELPER FUNCTIONS
# ==============================================================================

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

theme_publication <- function(base_size = 11) {
  theme_minimal(base_size = base_size) + 
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
      panel.background   = element_rect(fill = "transparent", color = NA),
      plot.background    = element_rect(fill = "transparent", color = NA),
      axis.text.y        = element_text(color = "black", size = base_size - 1, face = "bold"),
      axis.text.x        = element_text(color = "black", size = base_size - 1, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title         = element_text(color = "black", size = base_size, face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.6),
      axis.ticks         = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length  = unit(0.15, "cm"),
      legend.position    = "top",
      legend.title       = element_text(face = "bold"),
      plot.margin        = margin(t = 10, r = 14, b = 10, l = 14, unit = "pt")
    )
}

prep_sensitivity_data <- function(data, var_name) {
  data %>%
    mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
    group_by(PLOT_LABEL) %>%
    summarize(
      LATITUDE = mean(LATITUDE, na.rm = TRUE), # Capture latitude for sorting
      Mean     = mean(.data[[var_name]], na.rm = TRUE),
      `75th`   = quantile(.data[[var_name]], probs = 0.75, na.rm = TRUE),
      `90th`   = quantile(.data[[var_name]], probs = 0.90, na.rm = TRUE),
      Max      = max(.data[[var_name]], na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    # Order by Latitude descending (Northwestern Islands on left, Main Islands on right)
    arrange(desc(LATITUDE)) %>% 
    mutate(PLOT_LABEL = factor(PLOT_LABEL, levels = PLOT_LABEL)) 
}

# ==============================================================================
# 2. LOAD & PREPARE DATASETS
# ==============================================================================

strs_sw_csv <- here("data", "HA_2024_StRS_net_sw.csv")
if (!file.exists(strs_sw_csv)) strs_sw_csv <- here("HA_2024_StRS_net_sw.csv")

if (file.exists(strs_sw_csv)) {
  strs_sw <- read_csv(strs_sw_csv, show_col_types = FALSE)
  
  if (!"CORAL_COVER_PERC" %in% names(strs_sw) && "TOTAL" %in% names(strs_sw)) {
    strs_sw$CORAL_COVER_PERC <- strs_sw$TOTAL
  }
} else {
  stop("Dataset HA_2024_StRS_net_sw.csv not found.")
}

coral_stats <- prep_sensitivity_data(strs_sw, "CORAL_COVER_PERC")
prod_stats  <- prep_sensitivity_data(strs_sw, "GROSS_PROD_KG_M2_YR")

# ==============================================================================
# 3. BASE THEMES & AESTHETICS
# ==============================================================================

color_mean <- "#4575b4"
color_75th <- "#74add1"
color_90th <- "#fdae61"
color_max  <- "#d73027"

palette_sens <- c(
  "Mean" = color_mean, 
  "75th" = color_75th, 
  "90th" = color_90th, 
  "Max"  = color_max
)

# ==============================================================================
# 4. BUILD SENSITIVITY PLOTS
# ==============================================================================

build_sensitivity_plot <- function(data_wide, y_label) {
  
  data_long <- data_wide %>%
    pivot_longer(
      cols = c("Mean", "75th", "90th", "Max"), 
      names_to = "Statistic", 
      values_to = "Value"
    ) %>%
    mutate(Statistic = factor(Statistic, levels = c("Mean", "75th", "90th", "Max")))
  
  p <- ggplot() +
    # Draw dashed line perfectly centered vertically connecting Mean to Max
    geom_segment(
      data = data_wide, 
      aes(x = PLOT_LABEL, xend = PLOT_LABEL, y = Mean, yend = Max), 
      color = "gray60", linetype = "dashed", alpha = 0.7, linewidth = 0.6
    ) +
    # Draw the points stacked exactly on the dashed line
    geom_point(
      data = data_long, 
      aes(x = PLOT_LABEL, y = Value, fill = Statistic), 
      shape = 21, color = "black", size = 2.8, stroke = 0.5, alpha = 0.85
    ) +
    scale_fill_manual(values = palette_sens) +
    labs(x = NULL, y = y_label) +
    theme_publication()
  
  return(p)
}

# ==============================================================================
# 5. BUILD & ASSEMBLE MAIN PLOTS
# ==============================================================================

lbl_cov_sens  <- "Coral Cover (%)"
lbl_prod_sens <- expression(bold("Gross Production Rate (kg"~m^{-2}~yr^{-1}*")"))

p_coral <- build_sensitivity_plot(coral_stats, y_label = lbl_cov_sens)
p_prod  <- build_sensitivity_plot(prod_stats, y_label = lbl_prod_sens)

final_figure <- ggarrange(
  p_coral, p_prod, 
  ncol = 1, nrow = 2, 
  align = "v", labels = "auto",
  font.label = list(size = 12, face = "bold"),
  common.legend = TRUE, legend = "top"
)

# ==============================================================================
# 6. EXPORT FORMATS
# ==============================================================================

plot_width  <- 18.0
plot_height <- 22.0

# PNG Export (300 DPI)
ggsave(
  filename = file.path(out_dir, "FigureS3_Sensitivity.png"), 
  plot     = final_figure, 
  width    = plot_width, 
  height   = plot_height, 
  units    = "cm", 
  dpi      = 300, 
  bg       = "white"
)

# TIFF Export (300 DPI, LZW)
ggsave(
  filename    = file.path(out_dir, "FigureS3_Sensitivity.tiff"), 
  plot        = final_figure, 
  width       = plot_width, 
  height      = plot_height, 
  units       = "cm", 
  dpi         = 300, 
  compression = "lzw", 
  bg          = "white"
)

final_figure