# ==============================================================================
# Script: Figure4_VARvsSLR.R
# Hawaiian Archipelago - Island-Level Reef Vertical Accretion Rates vs. SLR
# Aesthetic: Red-Gray-Blue Color Scheme, Seamless SLR Band, 100% Manual Labels
# Inputs: outputs/2024_HI_StRS_Islandlevel_DarwinPoint.csv
# Outputs: outputs/Figure4_VARvsSLR.png & .tiff
# ==============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)

# --- SLR BENCHMARK CONSTANTS ---
hi_slr_min <- 6   # Predicted min SLR rate by 2050 (mm/yr)
hi_slr_max <- 10  # Predicted max SLR rate by 2050 (mm/yr)

# --- HELPER FUNCTIONS ---
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

theme_publication <- function(base_size = 9) {
  theme_minimal(base_size = base_size) + 
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA),
      axis.text         = element_text(color = "gray30", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size, face = "bold"),
      axis.title.y      = element_text(margin = margin(r = 4)),
      axis.title.x      = element_text(margin = margin(t = 6)),
      axis.line         = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks        = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      plot.margin       = margin(t = 15, r = 20, b = 10, l = 15, unit = "pt")
    )
}

# --- 1. LOAD DATA ---
csv_path <- here("outputs", "2024_HI_StRS_Islandlevel_DarwinPoint.csv")
if (!file.exists(csv_path)) csv_path <- here("data", "2024_HI_StRS_Islandlevel_DarwinPoint.csv")
if (!file.exists(csv_path)) csv_path <- here("2024_HI_StRS_Islandlevel_DarwinPoint.csv")

if (!file.exists(csv_path)) stop("Dataset 2024_HI_StRS_Islandlevel_DarwinPoint.csv not found.")

spatial_Imean <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(PLOT_LABEL = format_islands(as.character(ISLAND)))

se_col <- if ("se.VERTICAL_ACCRETION_MM_YR" %in% names(spatial_Imean)) {
  "se.VERTICAL_ACCRETION_MM_YR"
} else if ("VERTICAL_ACCRETION_MM_YR_SE" %in% names(spatial_Imean)) {
  "VERTICAL_ACCRETION_MM_YR_SE"
} else {
  "se.VERTICAL_ACCRETION_MM_YR"
}

# --- 2. MANUAL LABEL NUDGES ---
spatial_Imean <- spatial_Imean %>%
  mutate(
    nudge_x = case_when(
      PLOT_LABEL == "Hawai'i"            ~  0.8,
      PLOT_LABEL == "Kaho'olawe"         ~  0.5,
      PLOT_LABEL == "Maui"               ~  1.0,
      PLOT_LABEL == "L\u0101na'i"        ~ -0.2,
      PLOT_LABEL == "Moloka'i"           ~ -1.0,
      PLOT_LABEL == "O'ahu"              ~  0.3,
      PLOT_LABEL == "Kaua'i"             ~ -0.8,
      PLOT_LABEL == "Ni'ihau"            ~ -0.6,
      PLOT_LABEL == "Lalo"               ~ -0.8,
      PLOT_LABEL == "Kapou"              ~ -0.5,
      PLOT_LABEL == "Manawai"            ~  0.2,
      PLOT_LABEL == "H\u014dlanik\u016b" ~  0.8,
      TRUE ~ 0
    ),
    nudge_y = case_when(
      PLOT_LABEL == "Hawai'i"            ~  0.5,
      PLOT_LABEL == "Kaho'olawe"         ~  2.5,
      PLOT_LABEL == "Maui"               ~ -0.2,
      PLOT_LABEL == "L\u0101na'i"        ~  1.5,
      PLOT_LABEL == "Moloka'i"           ~  0.8,
      PLOT_LABEL == "O'ahu"              ~ -1.5,
      PLOT_LABEL == "Kaua'i"             ~  0.8,
      PLOT_LABEL == "Ni'ihau"            ~ -1.2,
      PLOT_LABEL == "Lalo"               ~  1.0,
      PLOT_LABEL == "Kapou"              ~  1.5,
      PLOT_LABEL == "Manawai"            ~  1.5,
      PLOT_LABEL == "H\u014dlanik\u016b" ~ -1.5,
      TRUE ~ 0
    )
  )

# --- 3. BUILD FIGURE 4 ---
lbl_vacc <- expression(bold("Vertical accretion rate (mm"~yr^{-1}*")"))

p_fig4 <- ggplot(spatial_Imean, aes(x = LATITUDE, y = VERTICAL_ACCRETION_MM_YR)) +
  # 1. Seamless Full-Panel Width SLR Shaded Band (6 - 10 mm/yr)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = hi_slr_min, ymax = hi_slr_max,
           fill = "#F97316", alpha = 0.09) +
  
  # 2. Crisp Dotted Boundary Lines at 6 and 10 mm/yr
  geom_hline(yintercept = c(hi_slr_min, hi_slr_max), 
             linetype = "dotted", color = "#EA580C", linewidth = 0.45) +
  
  # 3. Direct Clean Text Annotation (No Background Box)
  annotate("text", x = 29.8, y = (hi_slr_min + hi_slr_max) / 2, 
           label = paste0("Predicted sea level rise by 2050"),
           color = "#B33F00", size = 2.8, fontface = "bold", hjust = 0) +
  
  # Net Accretion / Net Erosion Equilibrium Baseline (0 mm/yr)
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.55) +
  
  # Error Bars (SE)
  geom_errorbar(
    aes(ymin = VERTICAL_ACCRETION_MM_YR - .data[[se_col]], 
        ymax = VERTICAL_ACCRETION_MM_YR + .data[[se_col]]), 
    color = "gray55", width = 0, linewidth = 0.45
  ) +
  
  # Island Points with Red (Negative VAR) -> Gray (Zero VAR) -> Blue (Positive VAR) Scheme
  geom_point(
    aes(fill = VERTICAL_ACCRETION_MM_YR), 
    shape = 21, color = "gray20", stroke = 0.5, size = 3.4, alpha = 0.95
  ) +
  scale_fill_gradient2(
    low      = "#D7191C",  # Red for Net Framework Erosion (< 0 mm/yr)
    mid      = "#94A3B8",  # Gray for Near-Equilibrium (~ 0 mm/yr)
    high     = "#2C7BB6",  # Blue for Net Accretion (> 0 mm/yr)
    midpoint = 0,
    guide    = "none"
  ) +
  
  # Explicitly Mapped Manual Labels
  geom_text_repel(
    aes(label = PLOT_LABEL), 
    nudge_x           = spatial_Imean$nudge_x,
    nudge_y           = spatial_Imean$nudge_y,
    size              = 2.8, 
    fontface          = "bold", 
    color             = "black", 
    box.padding       = 0,         # Set to 0 to respect explicit nudges
    point.padding     = 0.2,       # Keeps leader line slightly off the point
    force             = 0.5,       # Low force allows manual coordinates to dominate
    max.overlaps      = Inf, 
    seed              = 42,
    segment.color     = "gray40", 
    segment.alpha     = 0.75, 
    segment.linewidth = 0.35, 
    min.segment.length= 0        
  ) +
  
  scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 2)) +
  coord_cartesian(ylim = c(-4.5, 12.5), clip = "off") +
  theme_publication() +
  labs(x = "Latitude (°N)", y = lbl_vacc)

# --- 4. EXPORT FIGURE ---
out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Figure4_VARvsSLR.png"),  p_fig4, width = 16.0, height = 11.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure4_VARvsSLR.tiff"), p_fig4, width = 16.0, height = 11.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")

p_fig4