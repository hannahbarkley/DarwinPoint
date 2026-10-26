# ==============================================================================
# Script: Figure3_FullStRS.R
# Hawaiian Archipelago - Full StRS Carbonate Budget Metrics Across Latitude (2x2)
# Aesthetic: Publication-Ready 2x2 Grid with Corrected Panel D Y-Axis
# Outputs: 
#   - outputs/Figure3_FullStRS.png & .tiff
# ==============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(ggpubr)

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
      axis.title.y      = element_text(margin = margin(r = 2)),
      axis.title.x      = element_text(margin = margin(t = 4)),
      axis.line         = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks        = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      plot.margin       = margin(t = 10, r = 4, b = 2, l = 22, unit = "pt")
    )
}

no_x_title <- function(p) {
  p + labs(x = NULL) + theme(axis.title.x = element_blank())
}

# --- 1. LOAD REQUIRED DATA ---
csv_path <- here("outputs", "2024_HI_StRS_Islandlevel_DarwinPoint.csv")
if (!file.exists(csv_path)) csv_path <- here("data", "2024_HI_StRS_Islandlevel_DarwinPoint.csv")
if (!file.exists(csv_path)) csv_path <- here("2024_HI_StRS_Islandlevel_DarwinPoint.csv")

if (!file.exists(csv_path)) stop("Dataset 2024_HI_StRS_Islandlevel_DarwinPoint.csv not found.")

spatial_Imean <- read_csv(csv_path, show_col_types = FALSE) %>%
  mutate(PLOT_LABEL = format_islands(as.character(ISLAND)))

# --- 2. EXTRACT UNIQUE ISLAND LATITUDES FOR LEGEND BAR ---
island_locs <- spatial_Imean %>%
  select(PLOT_LABEL, LATITUDE) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

# --- 3. PLOT BUILDER ---
build_plot <- function(data, y_var, y_se = NULL, y_label, y_lim, y_brk, 
                       add_zero_line = FALSE) {
  
  p <- ggplot(data, aes(x = LATITUDE, y = .data[[y_var]])) +
    scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
    theme_publication()
  
  if (!is.null(y_se)) {
    p <- p + geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_se]], ymax = .data[[y_var]] + .data[[y_se]]), 
      color = "gray55", width = 0, linewidth = 0.5
    )
  }
  
  p <- p + scale_y_continuous(breaks = y_brk) + coord_cartesian(ylim = y_lim, clip = "off")
  
  if (add_zero_line) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6)
  }
  
  p <- p + geom_point(aes(fill = LATITUDE), shape = 21, color = "gray30", stroke = 0.5, size = 2.8, alpha = 0.85) +
    scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none", limits = c(18, 30)) +
    labs(x = "Latitude (°N)", y = y_label)
  
  return(p)
}

# --- 4. LABELS ---
lbl_cov   <- "Coral cover (%)"
lbl_prod  <- expression(bold("Gross production rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_geros <- expression(bold("Gross erosion rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_nprod <- expression(bold("Net production rate (kg"~m^{-2}~yr^{-1}*")"))

# --- 5. BUILD SUBPLOTS ---

# Panel A: Coral Cover (0 to 30%)
p_strs_cov <- build_plot(
  spatial_Imean, "CORAL_COVER_PERC", "se.CORAL_COVER_PERC", 
  y_label = lbl_cov, y_lim = c(0, 30), y_brk = seq(0, 30, 5)
)

# Panel B: Gross Production Rate (0 to 6)
p_strs_gprod <- build_plot(
  spatial_Imean, "GROSS_PROD_KG_M2_YR", "se.GROSS_PROD_KG_M2_YR", 
  y_label = lbl_prod, y_lim = c(0, 6), y_brk = seq(0, 6, 1), 
  add_zero_line = FALSE
)

# Panel C: Gross Erosion Rate (0 to 4)
p_strs_geros <- build_plot(
  spatial_Imean, "GROSS_EROSION_KG_M2_YR", "se.GROSS_EROSION_KG_M2_YR", 
  y_label = lbl_geros, y_lim = c(0, 4), y_brk = seq(0, 4, 1)
)

# Panel D: Net Production Rate (Corrected Y-Lim: -4.5 to 4.2, Ticks: -4 to 4 by 2)
p_strs_nprod <- build_plot(
  spatial_Imean, "NET_PROD_KG_M2_YR", "se.NET_PROD_KG_M2_YR", 
  y_label = lbl_nprod, y_lim = c(-4.5, 4.2), y_brk = seq(-4, 4, by = 2), 
  add_zero_line = TRUE
)

# --- 6. ASSEMBLE 2x2 MAIN PLOTS ---
main_plots <- ggarrange(
  no_x_title(p_strs_cov),  no_x_title(p_strs_gprod),
  p_strs_geros,            p_strs_nprod,
  ncol = 2, nrow = 2, align = "hv",
  labels = "auto",
  font.label = list(size = 12, face = "bold"),
  label.x = 0.01, label.y = 1.02
)

# --- 7. BUILD SIDE LATITUDE LEGEND BAR ---
legend_bar <- data.frame(lat = seq(18, 30, length.out = 500), x = 0)

p_legend <- ggplot() +
  geom_tile(data = legend_bar, aes(x = x, y = lat, fill = lat), width = 0.5) +
  annotate("rect", xmin = -0.25, xmax = 0.25, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) +
  scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none", limits = c(18, 30)) +
  scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.02, 0.02)) +
  scale_x_continuous(limits = c(-0.25, 3.8), expand = c(0,0)) + 
  geom_point(data = island_locs, aes(x = 0.25, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2, color = "gray30", stroke = 0.5) +
  geom_text_repel(
    data = island_locs, aes(x = 0.25, y = LATITUDE, label = PLOT_LABEL),
    direction = "y", hjust = 0, nudge_x = 0.4, segment.color = "gray50", segment.size = 0.4,
    size = 2.6, fontface = "bold", color = "gray20", force = 3, box.padding = unit(0.3, "lines")
  ) +
  labs(title = "Latitude (°N)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 4)),
    panel.grid = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 5)),
    axis.ticks.y = element_line(color = "gray20"), axis.ticks.length = unit(0.15, "cm"),
    plot.margin = margin(t = 15, r = 0, b = 15, l = 10)
  )

# --- 8. FINAL ASSEMBLY & EXPORTS ---
fig3_full_strs <- ggarrange(main_plots, p_legend, ncol = 2, widths = c(1, 0.22))

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "Figure3_FullStRS.png"),  fig3_full_strs, width = 24.0, height = 14.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure3_FullStRS.tiff"), fig3_full_strs, width = 24.0, height = 14.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")