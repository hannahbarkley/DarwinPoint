# ==============================================================================
# Script: FigureS2_SfM_vs_CoralNet.R
# Hawaiian Archipelago - SfM vs CoralNet Gross Carbonate Production Comparison
# Compares GROSS_CARB_PROD_KG_M2_YR_MEAN (SfM) vs PROD (CoralNet)
# Outputs: 
#   - outputs/FigureS2_SfM_vs_CoralNet.png & .tiff
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

theme_publication <- function(base_size = 9.5) {
  theme_minimal(base_size = base_size) + 
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA),
      axis.text         = element_text(color = "gray25", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size, face = "bold"),
      axis.title.y      = element_text(margin = margin(r = 6)),
      axis.title.x      = element_text(margin = margin(t = 6)),
      axis.line         = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks        = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      plot.margin       = margin(t = 10, r = 4, b = 10, l = 10, unit = "pt")
    )
}

# --- 1. LOAD REQUIRED DATA ---
sfm_csv      <- here("data", "HA_2024_fixed_net.csv")
if (!file.exists(sfm_csv)) sfm_csv <- here("HA_2024_fixed_net.csv")

coralnet_csv <- here("data", "HA_2024_StRS_prod.csv")
if (!file.exists(coralnet_csv)) coralnet_csv <- here("HA_2024_StRS_prod.csv")

sfm_df      <- read_csv(sfm_csv, show_col_types = FALSE)
coralnet_df <- read_csv(coralnet_csv, show_col_types = FALSE)

# --- 2. PREPARE & MERGE MATCHED SITES ---
matched_data <- sfm_df %>%
  select(
    OCC_SITEID, LOCATION, LATITUDE, LONGITUDE,
    SFM_PROD = GROSS_CARB_PROD_KG_M2_YR_MEAN,
    SFM_SE   = GROSS_CARB_PROD_KG_M2_YR_SE
  ) %>%
  inner_join(
    coralnet_df %>% select(SITE, CORALNET_PROD = PROD, CORALNET_LWR = PROD_lwr, CORALNET_UPR = PROD_upr),
    by = c("OCC_SITEID" = "SITE")
  ) %>%
  mutate(
    ISLAND_LABEL = format_islands(LOCATION),
    # Standardize 95% CI bounds to approximate 1 SE for direct comparability with SFM_SE
    CORALNET_SE  = (CORALNET_UPR - CORALNET_LWR) / (2 * 1.96)
  )

# Island location coordinates for custom legend
island_locs <- matched_data %>%
  group_by(ISLAND_LABEL) %>%
  summarize(LATITUDE = mean(LATITUDE), .groups = "drop") %>%
  rename(PLOT_LABEL = ISLAND_LABEL) %>%
  arrange(LATITUDE)

# Maximum observed x-value for unified line cutoff
max_sfm_val <- max(matched_data$SFM_PROD)

# --- 3. BUILD MAIN SCATTER PLOT ---
x_label <- expression(bold("SfM-based Gross Production (kg"~m^{-2}~yr^{-1}*")"))
y_label <- expression(bold("CoralNet-based Gross Production (kg"~m^{-2}~yr^{-1}*")"))

p_main <- ggplot(matched_data, aes(x = SFM_PROD, y = CORALNET_PROD)) +
  # 1:1 Identity Line truncated at maximum observed data point
  geom_segment(
    x = 0, y = 0, xend = max_sfm_val, yend = max_sfm_val, 
    linetype = "dashed", color = "gray60", linewidth = 0.55
  ) +
  # Linear Regression Fit Line with 95% Confidence Band
  geom_smooth(
    method = "lm", formula = y ~ x, 
    color = "gray15", linewidth = 0.75, se = TRUE, fill = "gray80", alpha = 0.35
  ) +
  # Vertical SE (CoralNet)
  geom_errorbar(
    aes(ymin = CORALNET_PROD - CORALNET_SE, ymax = CORALNET_PROD + CORALNET_SE), 
    color = "gray70", width = 0, linewidth = 0.35
  ) +
  # Horizontal SE (SfM)
  geom_errorbarh(
    aes(xmin = SFM_PROD - SFM_SE, xmax = SFM_PROD + SFM_SE), 
    color = "gray70", height = 0, linewidth = 0.35
  ) +
  # Point markers colored by Latitude (Reversed Viridis Plasma: Yellow = 18°N, Dark Purple = 30°N)
  geom_point(
    aes(fill = LATITUDE), 
    shape = 21, color = "gray20", stroke = 0.55, size = 3.3, alpha = 0.90
  ) +
  scale_fill_viridis_c(
    option = "plasma", direction = -1, limits = c(18, 30), guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 14.0), breaks = seq(0, 14, by = 2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 14.0), breaks = seq(0, 14, by = 2), expand = c(0.01, 0.01)) +
  coord_fixed(ratio = 1, xlim = c(0, 14.0), ylim = c(0, 14.0), clip = "off") +
  labs(x = x_label, y = y_label) +
  theme_publication()

# --- 4. BUILD CUSTOM LEGEND ---
legend_bar <- data.frame(lat = seq(18, 30, length.out = 500), x = 0)

p_legend <- ggplot() +
  geom_tile(data = legend_bar, aes(x = x, y = lat, fill = lat), width = 0.45) +
  annotate("rect", xmin = -0.225, xmax = 0.225, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) +
  scale_fill_viridis_c(option = "plasma", direction = -1, guide = "none", limits = c(18, 30)) +
  scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(-0.25, 3.2), expand = c(0,0)) + 
  geom_point(data = island_locs, aes(x = 0.225, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2.1, color = "gray30", stroke = 0.5) +
  geom_text_repel(
    data = island_locs, aes(x = 0.225, y = LATITUDE, label = PLOT_LABEL),
    direction = "y", hjust = 0, nudge_x = 0.35, 
    segment.color = "gray50", segment.alpha = 0.8, segment.linewidth = 0.35,
    size = 2.7, fontface = "bold", color = "gray20", force = 2,
    min.segment.length = 0
  ) +
  labs(title = "Latitude (°N)") +
  theme_minimal() +
  theme(
    plot.title        = element_text(size = 8.5, face = "bold", hjust = 0, margin = margin(b = 4)),
    panel.grid        = element_blank(), 
    axis.title        = element_blank(), 
    axis.text.x       = element_blank(), 
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 4)),
    axis.ticks.y      = element_line(color = "gray20"), 
    axis.ticks.length = unit(0.12, "cm"),
    plot.margin       = margin(t = 12, r = 0, b = 12, l = 0)
  )

# --- 5. FINAL ASSEMBLY & EXPORT ---
final_s2_plot <- ggarrange(p_main, p_legend, ncol = 2, widths = c(1, 0.28), align = "h")

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "FigureS2_SfM_vs_CoralNet.png"), final_s2_plot, width = 17.5, height = 13.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(out_dir, "FigureS2_SfM_vs_CoralNet.tiff"), final_s2_plot, width = 17.5, height = 13.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")

final_s2_plot