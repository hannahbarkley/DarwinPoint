# ==============================================================================
# Script: FigureS3_Porosity_Sensitivity.R
# Hawaiian Archipelago - Porosity Sensitivity Analysis on Vertical Accretion
# Aesthetic: Matching FigureS3 Cleveland dot plots with vertical segments
# Outputs: outputs/FigureS3_Porosity_Sensitivity.png & .tiff
# ==============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(ggpubr)

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. PARAMETERS & HELPER FUNCTIONS
# ==============================================================================

porosity_levels <- c(0.30, 0.40, 0.50, 0.60, 0.70)
rho_skeletal    <- 1.66 # Coral skeletal bulk density baseline (g cm^-3)

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

# Calculates VAR (mm/yr) for each island across porosity levels
prep_porosity_data <- function(data, prod_col, porosities = porosity_levels, rho = rho_skeletal) {
  island_summary <- data %>%
    mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
    group_by(PLOT_LABEL) %>%
    summarize(
      LATITUDE  = mean(LATITUDE, na.rm = TRUE),
      Mean_Prod = mean(.data[[prod_col]], na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    arrange(desc(LATITUDE)) %>%
    mutate(PLOT_LABEL = factor(PLOT_LABEL, levels = PLOT_LABEL))
  
  # Calculate VAR = P_net / (rho * (1 - phi)) for each porosity value
  porosity_df <- expand_grid(island_summary, Porosity = porosities) %>%
    mutate(
      VAR_mmyr = Mean_Prod / (rho * (1 - Porosity)),
      Porosity_Label = factor(
        paste0(Porosity * 100, "%"),
        levels = paste0(porosities * 100, "%")
      )
    )
  
  return(porosity_df)
}

# ==============================================================================
# 2. LOAD & PREPARE DATASETS
# ==============================================================================

strs_sw_csv <- here("data", "HA_2024_StRS_net_sw.csv")
if (!file.exists(strs_sw_csv)) strs_sw_csv <- here("HA_2024_StRS_net_sw.csv")

if (file.exists(strs_sw_csv)) {
  strs_sw <- read_csv(strs_sw_csv, show_col_types = FALSE)
  
  # Determine production variable column
  prod_var <- if ("NET_PROD_KG_M2_YR" %in% names(strs_sw)) {
    "NET_PROD_KG_M2_YR"
  } else if ("GROSS_PROD_KG_M2_YR" %in% names(strs_sw)) {
    "GROSS_PROD_KG_M2_YR"
  } else {
    stop("Neither NET_PROD_KG_M2_YR nor GROSS_PROD_KG_M2_YR found in dataset.")
  }
} else {
  stop("Dataset HA_2024_StRS_net_sw.csv not found.")
}

var_porosity_stats <- prep_porosity_data(strs_sw, prod_col = prod_var)

# ==============================================================================
# 3. COLOR PALETTE & PLOT BUILDER
# ==============================================================================

palette_porosity <- c(
  "30%" = "#2c7bb6",
  "40%" = "#abd9e9",
  "50%" = "#ffffbf",
  "60%" = "#fdae61",
  "70%" = "#d7191c"
)

# Wide format bounds for vertical dashed line (30% to 70% porosity span)
var_bounds <- var_porosity_stats %>%
  group_by(PLOT_LABEL) %>%
  summarize(
    Y_Min = min(VAR_mmyr, na.rm = TRUE),
    Y_Max = max(VAR_mmyr, na.rm = TRUE),
    .groups = "drop"
  )

p_var_porosity <- ggplot() +
  # Vertical dashed segment connecting 30% to 70% porosity bounds
  geom_segment(
    data = var_bounds,
    aes(x = PLOT_LABEL, xend = PLOT_LABEL, y = Y_Min, yend = Y_Max),
    color = "gray60", linetype = "dashed", alpha = 0.7, linewidth = 0.6
  ) +
  # Points stacked along segment for each porosity level
  geom_point(
    data = var_porosity_stats,
    aes(x = PLOT_LABEL, y = VAR_mmyr, fill = Porosity_Label),
    shape = 21, color = "black", size = 2.8, stroke = 0.5, alpha = 0.85
  ) +
  scale_fill_manual(values = palette_porosity, name = "Framework Porosity (ϕ)") +
  labs(
    x = NULL,
    y = expression(bold("Vertical Accretion Rate (mm "*yr^{-1}*")"))
  ) +
  theme_publication()

# ==============================================================================
# 4. EXPORT FIGURES
# ==============================================================================

plot_width  <- 18.0
plot_height <- 12.0

# PNG Export (300 DPI)
ggsave(
  filename = file.path(out_dir, "Porosity_Sensitivity.png"), 
  plot     = p_var_porosity, 
  width    = plot_width, 
  height   = plot_height, 
  units    = "cm", 
  dpi      = 300, 
  bg       = "white"
)

# TIFF Export (300 DPI, LZW)
ggsave(
  filename    = file.path(out_dir, "Porosity_Sensitivity.tiff"), 
  plot        = p_var_porosity, 
  width       = plot_width, 
  height      = plot_height, 
  units       = plot_units <- "cm", 
  dpi         = 300, 
  compression = "lzw", 
  bg          = "white"
)

p_var_porosity