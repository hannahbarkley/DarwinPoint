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
pal_lat_vals <- c(
  seq(0, 0.3409, length.out = 6)[1:5], 
  0.3409, 
  seq(0.3409, 1, length.out = 6)[2:6]
)

# ------------------------------------------------------------------------------
# Helpers: Hawaiian Archipelago Names & Publication Theme
# ------------------------------------------------------------------------------
island_lookup <- c(
  "kure"                   = "H\u014dlanik\u016b",
  "kure atoll"             = "H\u014dlanik\u016b",
  "midway"                 = "Kuaihelani",
  "midway atoll"           = "Kuaihelani",
  "pearl & hermes"         = "Manawai",
  "pearl and hermes"       = "Manawai",
  "pearl & hermes atoll"   = "Manawai",
  "pearl_hermes"           = "Manawai",
  "lisianski"              = "Kapou",
  "lisianski island"       = "Kapou",
  "laysan"                 = "Kamole",
  "laysan island"          = "Kamole",
  "maro"                   = "Kamokuokamohoali'i",
  "maro reef"              = "Kamokuokamohoali'i",
  "gardner"                = "P\u016bh\u0101honu",
  "gardner pinnacles"      = "P\u016bh\u0101honu",
  "gardner_pinnacles"      = "P\u016bh\u0101honu",
  "french frigate"         = "Lalo",
  "french frigate shoals"  = "Lalo",
  "french_frigate"         = "Lalo",
  "french_frigate_shoals"  = "Lalo",
  "ffs"                    = "Lalo",
  "necker"                 = "Mokumanamana",
  "necker island"          = "Mokumanamana",
  "nihoa"                  = "Nihoa",
  "nihoa island"           = "Nihoa",
  "kauai"                  = "Kaua'i",
  "niihau"                 = "Ni'ihau",
  "oahu"                   = "O'ahu",
  "molokai"                = "Moloka'i",
  "lanai"                  = "L\u0101na'i",
  "maui"                   = "Maui",
  "kahoolawe"              = "Kaho'olawe",
  "hawaii"                 = "Hawai'i"
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
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.4),
      axis.text         = element_text(color = "gray30", size = base_size - 1),
      axis.title        = element_text(color = "black", size = base_size, face = "bold"),
      axis.title.y      = element_text(margin = margin(r = y_margin_r)),
      axis.title.x      = element_text(margin = margin(t = x_margin_t)),
      axis.line         = element_blank(),
      axis.ticks        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      plot.margin       = margin(t = 10, r = 10, b = 6, l = 14, unit = "pt")
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
                                             if ("PROD_SE" %in% names(.)) PROD_SE else NULL,
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
# 2. MODEL PREDICTION FOR FULL-RANGE FIT (0 to 18)
# ------------------------------------------------------------------------------
lm_sfm_cnet <- lm(CORALNET_PROD ~ SFM_PROD, data = matched_s3)
pred_grid   <- tibble(SFM_PROD = seq(0, 18, length.out = 300))
pred_fit    <- predict(lm_sfm_cnet, newdata = pred_grid, se.fit = TRUE, interval = "confidence", level = 0.95)

df_lm_fit <- pred_grid %>%
  mutate(
    fit  = pred_fit$fit[, "fit"],
    ymin = pred_fit$fit[, "lwr"],
    ymax = pred_fit$fit[, "upr"]
  )

# ------------------------------------------------------------------------------
# 3. PLOT CONSTRUCTION
# ------------------------------------------------------------------------------

# Main Scatter Plot
p_main_s3 <- ggplot(matched_s3, aes(x = SFM_PROD, y = CORALNET_PROD)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", linewidth = 0.55) +
  geom_ribbon(
    data = df_lm_fit,
    aes(x = SFM_PROD, ymin = ymin, ymax = ymax),
    fill = "#64748B", alpha = 0.15, inherit.aes = FALSE
  ) +
  geom_line(
    data = df_lm_fit,
    aes(x = SFM_PROD, y = fit),
    color = "#334155", linewidth = 0.75, inherit.aes = FALSE
  ) +
  geom_errorbar(aes(ymin = CORALNET_PROD - CORALNET_SE, ymax = CORALNET_PROD + CORALNET_SE), color = "gray60", width = 0, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = SFM_PROD - SFM_SE, xmax = SFM_PROD + SFM_SE), color = "gray60", width = 0, linewidth = 0.35) +
  geom_point(aes(fill = LATITUDE), shape = 21, color = "gray20", stroke = 0.6, size = 3.0, alpha = 0.85) +
  scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
  scale_x_continuous(breaks = seq(0, 18, by = 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 18, by = 2), expand = c(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(0, 18), ylim = c(0, 18), clip = "on") +
  theme_publication(base_size = 9, y_margin_r = 4, x_margin_t = 6) +
  labs(
    x = expression(bold("SfM" ~ G[gross] ~ "(kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")),
    y = expression(bold("CoralNet" ~ G[gross] ~ "(kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
  )

# Shorter Vertical Legend aligned to Scatter Plot Box
bar_w  <- 0.16
half_w <- bar_w / 2

df_brackets <- tibble(
  region = c("Northwestern Hawaiian Islands", "Main Hawaiian Islands"),
  ymin   = c(23.0, 19.0),
  ymax   = c(29.8, 22.7),
  ymid   = (ymin + ymax) / 2
)

p_legend_s3 <- ggplot() +
  geom_tile(
    data = data.frame(lat = seq(19, 30, length.out = 500), x = 0), 
    aes(x = x, y = lat, fill = lat), 
    width = bar_w
  ) +
  annotate(
    "rect", 
    xmin = -half_w, xmax = half_w, 
    ymin = 19, ymax = 30, 
    color = "gray20", fill = NA, linewidth = 0.55
  ) +
  scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
  scale_y_continuous(
    limits = c(18.8, 30.2), 
    breaks = seq(19, 30, by = 1), 
    expand = c(0, 0)
  ) +
  scale_x_continuous(limits = c(-half_w, 2.45), expand = c(0, 0)) +
  geom_point(
    data   = island_locs_s3, 
    aes(x = 0, y = LATITUDE, fill = LATITUDE), 
    shape  = 21, size = 2.4, color = "gray20", stroke = 0.55
  ) +
  geom_text_repel(
    data               = island_locs_s3, 
    aes(x = 0, y = LATITUDE, label = PLOT_LABEL), 
    direction          = "y", hjust = 0, nudge_x = 0.45, size = 2.6, 
    fontface           = "bold", color = "gray20", 
    box.padding        = 0.22, point.padding = 0.25, 
    force              = 2.2, seed = 42, 
    segment.color      = "gray50", segment.size = 0.35, segment.alpha = 0.8, 
    min.segment.length = 0
  ) +
  geom_segment(data = df_brackets, aes(x = 1.95, xend = 1.95, y = ymin, yend = ymax), color = "gray20", linewidth = 0.5) +
  geom_segment(data = df_brackets, aes(x = 1.95, xend = 1.83, y = ymax, yend = ymax), color = "gray20", linewidth = 0.5) +
  geom_segment(data = df_brackets, aes(x = 1.95, xend = 1.83, y = ymin, yend = ymin), color = "gray20", linewidth = 0.5) +
  geom_text(
    data = df_brackets, 
    aes(x = 2.18, y = ymid, label = region), 
    angle = -90, size = 2.6, fontface = "bold", color = "gray20"
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.title        = element_text(size = 8.5, face = "bold", color = "black", hjust = 0.08, margin = margin(b = 6)),
    axis.text.y       = element_text(size = 7.5, color = "gray30", margin = margin(r = 3)),
    axis.ticks.y      = element_line(color = "gray20", linewidth = 0.55),
    axis.ticks.length = unit(0.12, "cm"),
    
    # Kept t = 28 to preserve top flush alignment, increased b to 50 to lift the bottom
    plot.margin       = margin(t = 28, r = 10, b = 50, l = 16, unit = "pt")
  ) +
  labs(title = "Latitude (°N)")

# ------------------------------------------------------------------------------
# 4. ASSEMBLY & EXPORT
# ------------------------------------------------------------------------------
# Removed align = "h" to prevent the square aspect-ratio plot from collapsing
figS4_final <- ggarrange(p_main_s3, p_legend_s3, ncol = 2, widths = c(1, 0.36))

ggsave(file.path(supp_dir, "FigureS4_SfM_vs_CoralNet.png"), figS4_final, width = 17.8, height = 13.5, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(supp_dir, "FigureS4_SfM_vs_CoralNet.tiff"), figS4_final, width = 17.8, height = 13.5, units = "cm", dpi = 600, compression = "lzw", bg = "white")