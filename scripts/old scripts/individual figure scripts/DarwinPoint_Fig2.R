# ==============================================================================
# Script: Figure2_Grigg_vs_NCRMP.R
# Hawaiian Archipelago - Historical (Grigg 1982) vs Modern (NCRMP) P90 Comparison
# Outputs: 
#   - outputs/Figure2_Grigg_vs_NCRMP.png & .tiff
# ==============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# --- PHYSICAL CONSTANTS ---
rho_aragonite <- 2.89  
reef_porosity <- 0.50  

# --- HELPER FUNCTIONS ---
format_islands <- function(x) {
  case_match(
    x,
    "Oahu"                 ~ "O'ahu", 
    "Hawaii"               ~ "Hawai'i", 
    "Lanai"                ~ "L\u0101na'i",
    "Molokai"              ~ "Moloka'i", 
    "Kauai"                ~ "Kaua'i", 
    "Niihau"               ~ "Ni'ihau",
    "Kahoolawe"            ~ "Kaho'olawe", 
    "FFS"                  ~ "Lalo",
    "French Frigate Shoals" ~ "Lalo",
    "French Frigate"       ~ "Lalo",
    "Lisianski"            ~ "Kapou",
    "Pearl & Hermes"       ~ "Manawai",
    "Pearl and Hermes"     ~ "Manawai",
    "Kure"                 ~ "H\u014dlanik\u016b",
    "Midway"               ~ "Kuaihelani",
    "Laysan"               ~ "Kamole",
    "Maro Reef"            ~ "Kamokuokamohoali'i",
    "Maro"                 ~ "Kamokuokamohoali'i",
    "Gardner Pinnacles"    ~ "P\u016bh\u0101honu",
    "Gardner"              ~ "P\u016bh\u0101honu",
    "Necker"               ~ "Mokumanamana",
    .default = x
  )
}

boot_q90_se <- function(x, B = 1000) {
  x <- na.omit(x)
  if (length(x) <= 1) return(0)
  boot_vals <- replicate(B, quantile(sample(x, replace = TRUE), probs = 0.90, na.rm = TRUE))
  sd(boot_vals)
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
      plot.margin       = margin(t = 14, r = 4, b = 2, l = 22, unit = "pt")
    )
}

no_x_title <- function(p) {
  p + labs(x = NULL) + theme(axis.title.x = element_blank())
}

# --- 1. LOAD REQUIRED DATA ---
strs_sw_csv <- here("data", "HA_2024_StRS_net_sw.csv")
if (!file.exists(strs_sw_csv)) strs_sw_csv <- here("HA_2024_StRS_net_sw.csv")

grigg_csv <- here("data", "Grigg_1982_data.csv")
if (!file.exists(grigg_csv)) grigg_csv <- here("Grigg_1982_data.csv")

strs_sw <- read_csv(strs_sw_csv, show_col_types = FALSE)
grigg   <- read_csv(grigg_csv, show_col_types = FALSE)

# Fallback column mapping if coral cover is named TOTAL in exported subset
if (!"CORAL_COVER_PERC" %in% names(strs_sw) && "TOTAL" %in% names(strs_sw)) {
  strs_sw$CORAL_COVER_PERC <- strs_sw$TOTAL
}

# --- 2. PREPARE MAIN DATA ---
grigg <- grigg %>% mutate(
  PLOT_LABEL = format_islands(LOCATION),
  VERTICAL_ACCRETION_MM_YR = GROSS_PROD / (rho_aragonite * (1 - reef_porosity))
)

set.seed(42)
strs_sw_p90 <- strs_sw %>%
  mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
  group_by(PLOT_LABEL) %>%
  summarize(
    LATITUDE                        = mean(LATITUDE, na.rm = TRUE),
    CORAL_COVER_P90                 = quantile(CORAL_COVER_PERC, probs = 0.90, na.rm = TRUE),
    CORAL_COVER_P90_SE              = boot_q90_se(CORAL_COVER_PERC),
    GROSS_CARB_PROD_KG_M2_YR_P90   = quantile(GROSS_PROD_KG_M2_YR, probs = 0.90, na.rm = TRUE),
    GROSS_CARB_PROD_KG_M2_YR_P90_SE = boot_q90_se(GROSS_PROD_KG_M2_YR),
    .groups                         = "drop"
  )

# --- 3. CALCULATE BOOTSTRAPPED P90 DELTAS ---
grigg_ref <- grigg %>% select(PLOT_LABEL, GRIGG_LAT = LATITUDE, GRIGG_COV = CORAL_COVER_PERC, GRIGG_PROD = GROSS_PROD)

shared_p90_deltas <- strs_sw_p90 %>%
  inner_join(grigg_ref, by = "PLOT_LABEL") %>%
  mutate(
    LATITUDE       = GRIGG_LAT,
    DELTA_COV_P90  = CORAL_COVER_P90 - GRIGG_COV,
    DELTA_COV_SE   = CORAL_COVER_P90_SE,
    DELTA_PROD_P90 = GROSS_CARB_PROD_KG_M2_YR_P90 - GRIGG_PROD,
    DELTA_PROD_SE  = GROSS_CARB_PROD_KG_M2_YR_P90_SE
  )

# --- 4. EXTRACT UNIQUE ISLAND LATITUDES ---
island_locs <- bind_rows(
  select(grigg, PLOT_LABEL, LATITUDE),
  select(strs_sw_p90, PLOT_LABEL, LATITUDE)
) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

# --- 5. EXACT NON-OVERLAPPING COORDINATE POSITIONS FOR ALL LABELS ---
# Note: On scale_x_reverse(), positive nudge_x shifts label LEFTward on page, negative shifts RIGHTward.

# Panel C: Delta Coral Cover Nudges
delta_cov_data <- shared_p90_deltas %>%
  mutate(
    nudge_x = case_when(
      PLOT_LABEL == "H\u014dlanik\u016b" ~  2.5, # Left towards 29.6°N (below y=0 line)
      PLOT_LABEL == "Manawai"          ~ -1.0, # Right towards 26.8°N
      PLOT_LABEL == "Kapou"            ~ -2.2, # Right (clear of tall vertical error bar)
      PLOT_LABEL == "Lalo"             ~ -1.4, # Left towards 25.0°N
      PLOT_LABEL == "Kaua'i"           ~  1.0, # Left towards 23.1°N
      PLOT_LABEL == "Ni'ihau"          ~  1.2, # Left towards 23.1°N
      PLOT_LABEL == "O'ahu"            ~ -2.0, # Left towards 22.6°N
      PLOT_LABEL == "Moloka'i"         ~ -1.0, # Right towards 20.1°N (above point)
      PLOT_LABEL == "L\u0101na'i"       ~ -1.2, # Right towards 19.6°N
      PLOT_LABEL == "Maui"             ~  1.6, # Right towards 19.6°N
      PLOT_LABEL == "Hawai'i"          ~  1.8, # Right towards 18.1°N
      TRUE                             ~  1.0
    ),
    nudge_y = case_when(
      PLOT_LABEL == "H\u014dlanik\u016b" ~ -2.0, # Down to y = -19.7% (well below y=0 baseline)
      PLOT_LABEL == "Manawai"          ~ -12.0, # Down to y = -28.4%
      PLOT_LABEL == "Kapou"            ~  8.0, # Up to y = +28.4% (right of error bar)
      PLOT_LABEL == "Lalo"             ~ -12.0, # Down to y = -70.7%
      PLOT_LABEL == "Kaua'i"           ~  8.0, # Up to y = -5.3%
      PLOT_LABEL == "Ni'ihau"          ~ -8.0, # Down to y = -38.5%
      PLOT_LABEL == "O'ahu"            ~  6.0, # Down to y = -65.9%
      PLOT_LABEL == "Moloka'i"         ~  8.0, # Up to y = +7.1%
      PLOT_LABEL == "L\u0101na'i"       ~ -8.0, # Down to y = -62.9%
      PLOT_LABEL == "Maui"             ~  8.0, # Up to y = -14.8%
      PLOT_LABEL == "Hawai'i"          ~  14.0, # Down to y = -74.1%
      TRUE                             ~  0.0
    )
  )

# Panel F: Delta Gross Production Rate Nudges
delta_prod_data <- shared_p90_deltas %>%
  mutate(
    nudge_x = case_when(
      PLOT_LABEL == "H\u014dlanik\u016b" ~  1.2, # Left towards 29.6°N
      PLOT_LABEL == "Manawai"          ~ -1.0, # Right towards 26.8°N
      PLOT_LABEL == "Kapou"            ~ -2.2, # Right
      PLOT_LABEL == "Lalo"             ~ -1.6, # Left towards 25.0°N
      PLOT_LABEL == "Kaua'i"           ~  1.0, # Left towards 23.1°N
      PLOT_LABEL == "Ni'ihau"          ~  1.2, # Left towards 23.1°N
      PLOT_LABEL == "O'ahu"            ~ -1.2, # Left towards 22.6°N
      PLOT_LABEL == "Moloka'i"         ~ -1.0, # Right towards 20.1°N
      PLOT_LABEL == "L\u0101na'i"       ~ -1.2, # Right towards 19.6°N
      PLOT_LABEL == "Maui"             ~  1.6, # Right towards 19.6°N
      PLOT_LABEL == "Hawai'i"          ~  1.6, # Right towards 18.1°N
      TRUE                             ~  1.0
    ),
    nudge_y = case_when(
      PLOT_LABEL == "H\u014dlanik\u016b" ~ -2.5, # Down to y = -2.58 (well below y=0 baseline)
      PLOT_LABEL == "Manawai"          ~  2.5, # Up to y = +2.80 (well above y=0 baseline)
      PLOT_LABEL == "Kapou"            ~  1.6, # Up to y = +9.58 (right of error bar)
      PLOT_LABEL == "Lalo"             ~ -1.3, # Down to y = -10.26
      PLOT_LABEL == "Kaua'i"           ~  2.0, # Up to y = +0.05
      PLOT_LABEL == "Ni'ihau"          ~ -2.5, # Down to y = -4.84
      PLOT_LABEL == "O'ahu"            ~ -2.8, # Down to y = -12.62
      PLOT_LABEL == "Moloka'i"         ~  2.5, # Up to y = +4.62
      PLOT_LABEL == "L\u0101na'i"       ~ -2.5, # Down to y = -7.35
      PLOT_LABEL == "Maui"             ~  2.4, # Up to y = -1.11
      PLOT_LABEL == "Hawai'i"          ~  2.5, # Down to y = -11.15
      TRUE                             ~  0.0
    )
  )

# --- 6. PLOT BUILDER ---
build_plot <- function(data, y_var, y_se = NULL, y_label, y_lim, y_brk, 
                       add_zero_line = FALSE, add_labels = FALSE) {
  
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
    # REVERSED VIRIDIS PALETTE: direction = -1 sets Yellow = 18°N (Lower Lat) & Dark Purple = 30°N
    scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none", limits = c(18, 30)) +
    labs(x = "Latitude (°N)", y = y_label)  
  
  if (add_labels) {
    p <- p + geom_text_repel(
      aes(label = PLOT_LABEL), 
      size = 2.4, 
      fontface = "bold", 
      color = "gray20", 
      box.padding = unit(0.25, "lines"),     
      point.padding = unit(0.35, "lines"),   
      nudge_x = data$nudge_x,        
      nudge_y = data$nudge_y,
      force = 1,              
      seed = 42,
      segment.color = "gray60",      # Lighter gray leader lines
      segment.alpha = 0.8, 
      segment.linewidth = 0.35, 
      min.segment.length = 0         # Always draw leader line
    )
  }
  return(p)
}

# --- 7. LABELS ---
lbl_cov        <- "Coral cover (%)"
lbl_prod       <- expression(bold("Gross production rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_delta_cov  <- "\u0394 Coral cover (%)"
lbl_delta_prod <- expression(bold(Delta~"Gross production rate (kg"~m^{-2}~yr^{-1}*")"))

# --- 8. BUILD SUBPLOTS ---
grigg_cover    <- build_plot(grigg, "CORAL_COVER_PERC", y_label = lbl_cov, y_lim = c(0, 100), y_brk = seq(0, 100, 20))
ncrmp_cover_sw <- build_plot(strs_sw_p90, "CORAL_COVER_P90", y_se = "CORAL_COVER_P90_SE", y_label = lbl_cov, y_lim = c(0, 100), y_brk = seq(0, 100, 20))

grigg_prod     <- build_plot(grigg, "GROSS_PROD", y_label = lbl_prod, y_lim = c(0, 20), y_brk = seq(0, 20, 5))
ncrmp_prod_sw  <- build_plot(strs_sw_p90, "GROSS_CARB_PROD_KG_M2_YR_P90", y_se = "GROSS_CARB_PROD_KG_M2_YR_P90_SE", y_label = lbl_prod, y_lim = c(0, 20), y_brk = seq(0, 20, 5))

# Subplots c & f with custom non-overlapping label offsets
delta_cover <- build_plot(
  delta_cov_data, "DELTA_COV_P90", y_se = "DELTA_COV_SE", 
  y_label = lbl_delta_cov, y_lim = c(-80, 40), y_brk = seq(-80, 40, 20), 
  add_zero_line = TRUE, add_labels = TRUE
)

delta_prod  <- build_plot(
  delta_prod_data, "DELTA_PROD_P90", y_se = "DELTA_PROD_SE", 
  y_label = lbl_delta_prod, y_lim = c(-15, 10), y_brk = seq(-15, 10, 5), 
  add_zero_line = TRUE, add_labels = TRUE
)

# --- 9. ASSEMBLE MAIN PLOTS ---
main_plots <- ggarrange(
  no_x_title(grigg_cover), no_x_title(ncrmp_cover_sw), no_x_title(delta_cover),
  grigg_prod,              ncrmp_prod_sw,              delta_prod,
  ncol = 3, nrow = 2, align = "hv", labels = "auto",
  font.label = list(size = 12, face = "bold"), 
  label.x = 0.01, label.y = 1.02, 
  legend = "none" 
)

# --- 10. BUILD CUSTOM LEGEND WITH REVERSED PALETTE ---
legend_bar <- data.frame(lat = seq(18, 30, length.out = 500), x = 0)

p_legend <- ggplot() +
  geom_tile(data = legend_bar, aes(x = x, y = lat, fill = lat), width = 0.5) +
  annotate("rect", xmin = -0.25, xmax = 0.25, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) +
  # REVERSED VIRIDIS PALETTE: direction = -1 sets Yellow = 18°N & Dark Purple = 30°N
  scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none", limits = c(18, 30)) +
  scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.02, 0.02)) +
  scale_x_continuous(limits = c(-0.25, 3.8), expand = c(0,0)) + 
  geom_point(data = island_locs, aes(x = 0.25, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2, color = "gray30", stroke = 0.5) +
  geom_text_repel(
    data = island_locs, aes(x = 0.25, y = LATITUDE, label = PLOT_LABEL),
    direction = "y", hjust = 0, nudge_x = 0.4, segment.color = "gray50", segment.size = 0.4,
    size = 2.6, fontface = "bold", color = "gray20", force = 2
  ) +
  labs(title = "Latitude (°N)") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 4)),
    panel.grid = element_blank(), axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 5)),
    axis.ticks.y = element_line(color = "gray20"), axis.ticks.length = unit(0.15, "cm"),
    plot.margin = margin(t = 20, r = 0, b = 20, l = 10)
  )

# --- 11. FINAL ASSEMBLY & EXPORTS ---
grigg_ncrmp_plot <- ggarrange(main_plots, p_legend, ncol = 2, widths = c(1, 0.22))

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Export grigg_ncrmp_plot (PNG & TIFF)
ggsave(file.path(out_dir, "Figure2_Grigg_vs_NCRMP.png"), grigg_ncrmp_plot, width = 28, height = 12.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "Figure2_Grigg_vs_NCRMP.tiff"), grigg_ncrmp_plot, width = 28, height = 12.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")