# ==============================================================================
# Script:  FigureS4_SfM_vs_CoralNet.R
# Purpose: Methodological Validation of SfM vs. CoralNet Gross Production
# Inputs:  data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_sfm.csv
#            (or HA_2024_Fixed_Net_SiteLevel_cnet.csv / HA_2024_Fixed_Net_SiteLevel.csv)
#          data/strs data/field data/HA_2024_StRS_Prod_SiteLevel.csv
# Outputs: outputs/supplementary figures/FigureS4_SfM_vs_CoralNet.png
#          outputs/supplementary figures/FigureS4_SfM_vs_CoralNet.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(ggpubr)
})

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# Color Palette preserved from Figures 2 and 4 (Tol Sunset 19N to 30N)
pal_lat_grad <- c("#A50026", "#DD3D2D", "#F67E4B", "#FDB366", "#FEDA8B", 
                  "#EAECCC", "#C2E4EF", "#98CAE1", "#6EA6CD", "#4A7BB7", "#364B9A")

# Calculate relative position for the Lat 22.75 breakpoint 
# For limits c(19, 30): (22.75 - 19) / (30 - 19) = ~0.3409
# Map the yellow center (#EAECCC) to 0.3409 and distribute the rest evenly
pal_lat_vals <- c(
  seq(0, 0.3409, length.out = 6)[1:5], 
  0.3409, 
  seq(0.3409, 1, length.out = 6)[2:6]
)

# ------------------------------------------------------------------------------
# Helpers: Hawaiian Archipelago Names & Publication Theme
# ------------------------------------------------------------------------------
island_lookup <- c(
  "kure"                  = "H\u014dlanik\u016b",
  "kure atoll"            = "H\u014dlanik\u016b",
  "midway"                = "Kuaihelani",
  "midway atoll"          = "Kuaihelani",
  "pearl & hermes"        = "Manawai",
  "pearl and hermes"      = "Manawai",
  "pearl & hermes atoll"  = "Manawai",
  "pearl_hermes"          = "Manawai",
  "lisianski"             = "Kapou",
  "lisianski island"      = "Kapou",
  "laysan"                = "Kamole",
  "laysan island"         = "Kamole",
  "maro"                  = "Kamokuokamohoali'i",
  "maro reef"             = "Kamokuokamohoali'i",
  "gardner"               = "P\u016bh\u0101honu",
  "gardner pinnacles"     = "P\u016bh\u0101honu",
  "gardner_pinnacles"     = "P\u016bh\u0101honu",
  "french frigate"        = "Lalo",
  "french frigate shoals" = "Lalo",
  "french_frigate"        = "Lalo",
  "french_frigate_shoals" = "Lalo",
  "ffs"                   = "Lalo",
  "necker"                = "Mokumanamana",
  "necker island"         = "Mokumanamana",
  "nihoa"                 = "Nihoa",
  "nihoa island"          = "Nihoa",
  "kauai"                 = "Kaua'i",
  "niihau"                = "Ni'ihau",
  "oahu"                  = "O'ahu",
  "molokai"               = "Moloka'i",
  "lanai"                 = "L\u0101na'i",
  "maui"                  = "Maui",
  "kahoolawe"             = "Kaho'olawe",
  "hawaii"                = "Hawai'i"
)

format_islands <- function(x) {
  clean_key <- tolower(trimws(as.character(x)))
  clean_key <- gsub("_", " ", clean_key)
  mapped    <- island_lookup[clean_key]
  coalesce(unname(mapped), as.character(x))
}

theme_publication <- function(base_size = 9, y_margin_r = 4, x_margin_t = 6) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid        = element_blank(),
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA),
      axis.text         = element_text(color = "gray30", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size, face = "bold"),
      axis.title.y      = element_text(margin = margin(r = y_margin_r)),
      axis.title.x      = element_text(margin = margin(t = x_margin_t)),
      axis.line         = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks        = element_line(color = "gray20", linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm"),
      plot.margin       = margin(t = 12, r = 12, b = 4, l = 16, unit = "pt")
    )
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "fixed data", "processed data"),
    here("data", "strs data", "field data"),
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
# 1. LOAD & HARMONIZE DATASETS
# ------------------------------------------------------------------------------
df_fixed_net <- saferead_csv(c(
  "HA_2024_Fixed_Net_SiteLevel_sfm.csv",
  "HA_2024_Fixed_Net_SiteLevel_cnet.csv", 
  "HA_2024_Fixed_Net_SiteLevel.csv"
)) %>%
  mutate(
    LONGITUDE                     = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE),
    LOCATION                      = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                             if ("ISLAND" %in% names(.)) ISLAND else NULL),
    GROSS_CARB_PROD_KG_M2_YR_MEAN = coalesce(if ("GROSS_CARB_PROD_KG_M2_YR_MEAN" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_MEAN else NULL,
                                             if ("PROD" %in% names(.)) PROD else NULL,
                                             0),
    GROSS_CARB_PROD_KG_M2_YR_SE   = coalesce(if ("GROSS_CARB_PROD_KG_M2_YR_SE" %in% names(.)) GROSS_CARB_PROD_KG_M2_YR_SE else NULL,
                                             0)
  )

df_coralnet_prod <- saferead_csv("HA_2024_StRS_Prod_SiteLevel.csv")

matched_s3 <- df_fixed_net %>% 
  transmute(
    OCC_SITEID = coalesce(if ("OCC_SITEID" %in% names(.)) OCC_SITEID else NULL,
                          if ("SITE" %in% names(.)) SITE else NULL),
    LOCATION, 
    LATITUDE, 
    LONGITUDE, 
    SFM_PROD   = as.numeric(GROSS_CARB_PROD_KG_M2_YR_MEAN), 
    SFM_SE     = as.numeric(GROSS_CARB_PROD_KG_M2_YR_SE)
  ) %>%
  inner_join(
    df_coralnet_prod %>% select(SITE, CORALNET_PROD = PROD, CORALNET_LWR = PROD_lwr, CORALNET_UPR = PROD_upr),
    by = c("OCC_SITEID" = "SITE")
  ) %>%
  mutate(
    ISLAND_LABEL = format_islands(LOCATION),
    CORALNET_SE  = (CORALNET_UPR - CORALNET_LWR) / (2 * 1.96)
  )

island_locs_s3 <- matched_s3 %>% 
  group_by(ISLAND_LABEL) %>% 
  summarize(LATITUDE = mean(LATITUDE), .groups = "drop") %>% 
  rename(PLOT_LABEL = ISLAND_LABEL) %>% 
  arrange(LATITUDE)

# ------------------------------------------------------------------------------
# 2. PLOT CONSTRUCTION
# ------------------------------------------------------------------------------

# Main Scatter Plot
p_main_s3 <- ggplot(matched_s3, aes(x = SFM_PROD, y = CORALNET_PROD)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.55) +
  geom_smooth(
    method    = "lm", 
    formula   = y ~ x, 
    color     = "#334155", 
    linewidth = 0.75, 
    se        = TRUE, 
    fill      = "#64748B", 
    alpha     = 0.15,
    fullrange = TRUE
  ) +
  geom_errorbar(aes(ymin = CORALNET_PROD - CORALNET_SE, ymax = CORALNET_PROD + CORALNET_SE), color = "gray60", width = 0, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = SFM_PROD - SFM_SE, xmax = SFM_PROD + SFM_SE), color = "gray60", width = 0, linewidth = 0.35) +
  geom_point(aes(fill = LATITUDE), shape = 21, color = "gray20", stroke = 0.6, size = 3.0, alpha = 0.85) +
  scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
  scale_x_continuous(breaks = seq(0, 16, by = 2), expand = c(0.01, 0.01)) +
  scale_y_continuous(breaks = seq(0, 16, by = 2), expand = c(0.01, 0.01)) +
  coord_fixed(ratio = 1, xlim = c(0, 16), ylim = c(0, 16), clip = "off") +
  theme_publication(base_size = 9, y_margin_r = 4, x_margin_t = 6) +
  theme(plot.margin = margin(12, 16, 10, 14, "pt")) +
  labs(
    x = expression(bold("SfM-based gross production (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")),
    y = expression(bold("CoralNet-based gross production (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
  )

# Figure 2-Style Vertical Sidebar Legend
bar_w  <- 0.16
half_w <- bar_w / 2

df_brackets <- tibble(
  region = c("Northwestern Hawaiian Islands", "Main Hawaiian Islands"),
  ymin   = c(23.0, 19.0),
  ymax   = c(29.8, 22.7),
  ymid   = (ymin + ymax) / 2
)

p_legend_s3 <- ggplot() +
  # 1. Vertical continuous gradient bar centered at x = 0
  geom_tile(
    data = data.frame(lat = seq(19, 30, length.out = 500), x = 0), 
    aes(x = x, y = lat, fill = lat), 
    width = bar_w
  ) +
  # Bounding frame matching point border styling
  annotate(
    "rect", 
    xmin = -half_w, xmax = half_w, 
    ymin = 19, ymax = 30, 
    color = "gray20", fill = NA, linewidth = 0.55
  ) +
  scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
  
  # 2. Y-axis with 1° ticks running from 19 to 30 °N
  scale_y_continuous(
    limits = c(18.8, 30.2), 
    breaks = seq(19, 30, by = 1), 
    expand = c(0, 0)
  ) +
  # Expanded upper X limit from 2.05 to 2.45 to accommodate rightward bracket shift
  scale_x_continuous(limits = c(-half_w, 2.45), expand = c(0, 0)) +
  
  # 3. Circular island points centered directly in the middle of the color bar
  geom_point(
    data   = island_locs_s3, 
    aes(x = 0, y = LATITUDE, fill = LATITUDE), 
    shape  = 21, 
    size   = 2.4, 
    color  = "gray20", 
    stroke = 0.55
  ) +
  
  # 4. Angled/repelled island labels with leader lines
  geom_text_repel(
    data               = island_locs_s3, 
    aes(x = 0, y = LATITUDE, label = PLOT_LABEL), 
    direction          = "y", 
    hjust              = 0, 
    nudge_x            = 0.45, 
    size               = 2.6, 
    fontface           = "bold", 
    color              = "gray20", 
    box.padding        = 0.22, 
    point.padding      = 0.25, 
    force              = 2.2, 
    seed               = 42,
    segment.color      = "gray50", 
    segment.size       = 0.35, 
    segment.alpha      = 0.8, 
    min.segment.length = 0
  ) +
  
  # 5. Regional Group Brackets (Shifted rightward from x = 1.60 to x = 1.95)
  # Vertical bracket spine
  geom_segment(
    data = df_brackets, 
    aes(x = 1.95, xend = 1.95, y = ymin, yend = ymax), 
    color = "gray20", linewidth = 0.5
  ) +
  # Top tick pointing inward from 1.95 to 1.83
  geom_segment(
    data = df_brackets, 
    aes(x = 1.95, xend = 1.83, y = ymax, yend = ymax), 
    color = "gray20", linewidth = 0.5
  ) +
  # Bottom tick pointing inward from 1.95 to 1.83
  geom_segment(
    data = df_brackets, 
    aes(x = 1.95, xend = 1.83, y = ymin, yend = ymin), 
    color = "gray20", linewidth = 0.5
  ) +
  # Regional text placed alongside bracket spine at x = 2.18
  geom_text(
    data = df_brackets, 
    aes(x = 2.18, y = ymid, label = region), 
    angle    = -90, 
    size     = 2.6, 
    fontface = "bold", 
    color    = "gray20"
  ) +
  
  # 6. Minimal styling with adequate right and left margins
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.title        = element_text(size = 8.5, face = "bold", color = "black", hjust = 0.08, margin = margin(b = 6)),
    axis.text.y       = element_text(size = 7.5, color = "gray30", margin = margin(r = 3)),
    axis.ticks.y      = element_line(color = "gray20", linewidth = 0.55),
    axis.ticks.length = unit(0.12, "cm"),
    plot.margin       = margin(12, 10, 10, 18, "pt")
  ) +
  labs(title = "Latitude (°N)")

# ------------------------------------------------------------------------------
# 3. ASSEMBLY & EXPORT
# ------------------------------------------------------------------------------
# Increased legend width proportion to 0.36 to cleanly house the expanded bracket clearance
figS4_final <- ggarrange(p_main_s3, p_legend_s3, ncol = 2, widths = c(1, 0.36), align = "h")

ggsave(file.path(supp_dir, "FigureS4_SfM_vs_CoralNet.png"), figS4_final, width = 17.8, height = 13.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS4_SfM_vs_CoralNet.tiff"), figS4_final, width = 17.8, height = 13.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")