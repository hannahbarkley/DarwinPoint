# ==============================================================================
# Script:  FigureS5_Bioerosion.R
# Purpose: Archipelago-Wide Urchin and Parrotfish Bioerosion 
# Inputs:  - data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_cnet.csv
#          - data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
# Outputs: - outputs/supplementary figures/FigureS5_Bioerosion.png
#          - outputs/supplementary figures/FigureS5_Bioerosion.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(patchwork)
})

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# Reversed Paul Tol Sunset (19N is Red, 30N is Blue)
pal_lat_grad <- c("#A50026", "#DD3D2D", "#F67E4B", "#FDB366", "#FEDA8B", 
                  "#EAECCC", "#C2E4EF", "#98CAE1", "#6EA6CD", "#4A7BB7", "#364B9A")

pal_lat_vals <- c(
  seq(0, 0.3409, length.out = 6)[1:5], 
  0.3409, 
  seq(0.3409, 1, length.out = 6)[2:6]
)

# Standardized Archipelago X-Scale Limits across all components
x_limits <- c(30.2, 18.6)
x_expand <- c(0.01, 0.01)

# ------------------------------------------------------------------------------
# 1. HELPERS: HAWAIIAN ARCHIPELAGO LOOKUP & PUBLICATION THEME
# ------------------------------------------------------------------------------
island_lookup <- c(
  "kure"                         = "H\u014dlanik\u016b",
  "kure atoll"                   = "H\u014dlanik\u016b",
  "midway"                       = "Kuaihelani",
  "midway atoll"                 = "Kuaihelani",
  "pearl & hermes"               = "Manawai",
  "pearl and hermes"             = "Manawai",
  "pearl & hermes atoll"         = "Manawai",
  "pearl_hermes"                 = "Manawai",
  "lisianski"                    = "Kapou",
  "lisianski island"             = "Kapou",
  "laysan"                       = "Kamole",
  "laysan island"                = "Kamole",
  "maro"                         = "Kamokuokamohoali'i",
  "maro reef"                    = "Kamokuokamohoali'i",
  "gardner"                      = "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
  "gardner pinnacles"            = "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
  "gardner_pinnacles"            = "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
  "p\u016bh\u0101honu"           = "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki",
  "french frigate"               = "Lalo",
  "french frigate shoals"        = "Lalo",
  "french_frigate"               = "Lalo",
  "french_frigate_shoals"        = "Lalo",
  "ffs"                          = "Lalo",
  "necker"                       = "Mokumanamana",
  "necker island"                = "Mokumanamana",
  "nihoa"                        = "Nihoa",
  "nihoa island"                 = "Nihoa",
  "kauai"                        = "Kaua'i",
  "niihau"                       = "Ni'ihau",
  "oahu"                         = "O'ahu",
  "molokai"                      = "Moloka'i",
  "lanai"                        = "L\u0101na'i",
  "maui"                         = "Maui",
  "kahoolawe"                    = "Kaho'olawe",
  "hawaii"                       = "Hawai'i"
)

format_islands <- function(x) {
  clean_key <- tolower(trimws(as.character(x)))
  clean_key <- gsub("_", " ", clean_key)
  mapped    <- island_lookup[clean_key]
  coalesce(unname(mapped), as.character(x))
}

plot_theme <- function(base_size = 9, y_margin_r = 4, x_margin_t = 5,
                              t_margin = 6, r_margin = 12, b_margin = 6, l_margin = 12) {
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
      axis.ticks.length = unit(0.12, "cm"),
      plot.tag          = element_text(color = "black", size = 11, face = "bold"),
      plot.tag.position = "topleft",
      plot.margin       = margin(t = t_margin, r = r_margin, b = b_margin, l = l_margin, unit = "pt")
    )
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "strs data", "processed data"),
    here("data", "fixed data", "processed data"),
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
# 2. DATA INGESTION & HARMONIZATION
# ------------------------------------------------------------------------------
df_fixed_net <- saferead_csv(c("HA_2024_Fixed_Net_SiteLevel_cnet.csv", "HA_2024_Fixed_Net_SiteLevel.csv")) %>%
  mutate(
    LOCATION                      = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                             if ("ISLAND" %in% names(.)) ISLAND else NULL),
    URCHIN_EROSION_KG_M2_YR_MEAN = coalesce(if ("URCHIN_EROSION_KG_M2_YR_MEAN" %in% names(.)) URCHIN_EROSION_KG_M2_YR_MEAN else NULL, 0)
  )

df_strs_darwin <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>%
  mutate(
    LOCATION                      = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                             if ("ISLAND" %in% names(.)) ISLAND else NULL),
    PLOT_LABEL                    = format_islands(LOCATION),
    FISH_EROSION_KG_M2_YR_ALL     = coalesce(if ("FISH_EROSION_KG_M2_YR_ALL" %in% names(.)) FISH_EROSION_KG_M2_YR_ALL else NULL,
                                             if ("FISH_EROSION_KG_M2_YR_ALL_MEAN" %in% names(.)) FISH_EROSION_KG_M2_YR_ALL_MEAN else NULL),
    se.FISH_EROSION_KG_M2_YR_ALL  = coalesce(if ("se.FISH_EROSION_KG_M2_YR_ALL" %in% names(.)) se.FISH_EROSION_KG_M2_YR_ALL else NULL,
                                             if ("se.FISH_EROSION_KG_M2_YR_ALL_MEAN" %in% names(.)) se.FISH_EROSION_KG_M2_YR_ALL_MEAN else NULL,
                                             0)
  )

df_urchin_island <- df_fixed_net %>%
  mutate(PLOT_LABEL = format_islands(LOCATION)) %>%
  group_by(PLOT_LABEL) %>%
  summarize(
    LATITUDE           = mean(LATITUDE, na.rm = TRUE),
    URCHIN_EROSION_AVG = mean(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE),
    se.URCHIN_EROSION  = if (n() > 1) sd(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE) / sqrt(n()) else 0,
    .groups            = "drop"
  )

# ------------------------------------------------------------------------------
# 3. PLOT CONSTRUCTION (GUARANTEED EQUAL SIZING)
# ------------------------------------------------------------------------------
build_plot_s5 <- function(data, y_var, y_se, y_label, y_lim, y_brk, panel_tag, show_x_title = TRUE) {
  p <- ggplot(data, aes(x = LATITUDE, y = .data[[y_var]])) +
    scale_x_reverse(limits = x_limits, breaks = seq(19, 30, by = 1), expand = x_expand) +
    scale_y_continuous(breaks = y_brk) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    geom_errorbar(
      aes(ymin = pmax(0, .data[[y_var]] - .data[[y_se]]), ymax = .data[[y_var]] + .data[[y_se]]),
      color = "gray55", width = 0, linewidth = 0.45
    ) +
    geom_point(
      aes(fill = LATITUDE),
      shape = 21, size = 2.6, color = "gray20", stroke = 0.5, alpha = 0.90
    ) +
    scale_fill_gradientn(
      colors   = pal_lat_grad, 
      values   = pal_lat_vals, 
      limits   = c(19, 30), 
      na.value = "transparent", 
      guide    = "none"
    ) +
    labs(
      x   = if (show_x_title) "Latitude (°N)" else NULL, 
      y   = y_label, 
      tag = panel_tag
    ) +
    plot_theme(base_size = 9)
  
  if (!show_x_title) {
    p <- p + theme(
      axis.title.x = element_blank(),
      plot.margin  = margin(t = 8, r = 12, b = 4, l = 12, unit = "pt")
    )
  } else {
    p <- p + theme(
      plot.margin  = margin(t = 8, r = 12, b = 6, l = 12, unit = "pt")
    )
  }
  return(p)
}

lbl_urchin <- expression(bold(E[urchin] ~ "(kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
lbl_fish   <- expression(bold(E[parrotfish] ~ "(kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))

# Panel A: Has tick marks and numbers, omits redundant "Latitude (°N)" title
p_urchin <- build_plot_s5(
  df_urchin_island, "URCHIN_EROSION_AVG", "se.URCHIN_EROSION", 
  lbl_urchin, c(0, 4), seq(0, 4, 1), 
  panel_tag = "A", show_x_title = FALSE
)

# Panel B: Full x-axis including "Latitude (°N)" title
p_fish <- build_plot_s5(
  df_strs_darwin, "FISH_EROSION_KG_M2_YR_ALL", "se.FISH_EROSION_KG_M2_YR_ALL", 
  lbl_fish, c(0, 2), seq(0, 2, 0.5), 
  panel_tag = "B", show_x_title = TRUE
)

# ------------------------------------------------------------------------------
# 4. TOP HORIZONTAL LEGEND
# ------------------------------------------------------------------------------
island_locs_s5 <- df_strs_darwin %>%
  select(PLOT_LABEL, LATITUDE) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

legend_nudges <- tribble(
  ~PLOT_LABEL,                             ~leg_nx, ~leg_ny,
  "H\u014dlanik\u016b",                     -0.25,   0.55,
  "Kuaihelani",                            -0.10,   0.55,
  "Manawai",                               -0.10,   0.55,
  "Kapou",                                 -0.25,   0.55,
  "Kamole",                                -0.15,   0.55,
  "Kamokuokamohoali'i",                    -0.15,   0.55,
  "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki", -0.10, 0.55,
  "Lalo",                                  -0.25,   0.55,
  "Mokumanamana",                          -0.15,   0.55,
  "Nihoa",                                 -0.15,   0.55,
  "Kaua'i",                                -0.45,   0.55,
  "Ni'ihau",                               -0.20,   0.55,
  "O'ahu",                                 -0.30,   0.55,
  "Moloka'i",                              -0.15,   0.55,
  "L\u0101na'i",                            0.00,   0.55,
  "Maui",                                  -0.25,   0.55,
  "Kaho'olawe",                             0.00,   0.55,
  "Hawai'i",                                0.00,   0.55
)

island_locs_s5 <- island_locs_s5 %>%
  left_join(legend_nudges, by = "PLOT_LABEL") %>%
  mutate(leg_nx = replace_na(leg_nx, 0), leg_ny = replace_na(leg_ny, 0.55))

p_legend_s5 <- ggplot() +
  # Color bar bounded between 19 and 30
  geom_tile(
    data = data.frame(lat = seq(19, 30, length.out = 500), y = 0),
    aes(x = lat, y = y, fill = lat),
    height = 0.45
  ) +
  annotate("rect", xmin = 19, xmax = 30, ymin = -0.225, ymax = 0.225, 
           color = "gray20", fill = NA, linewidth = 0.45) +
  
  # Latitude label placed to the left of the 30°N boundary
  annotate("text", x = 30.12, y = 0, label = "Latitude (°N)", 
           fontface = "bold", size = 2.8, hjust = 1, color = "black") +
  
  # NWHI Bracket Bar
  annotate("segment", x = 29.8, xend = 23.0, y = 2.15, yend = 2.15, color = "gray30", linewidth = 0.5) +
  annotate("segment", x = 29.8, xend = 29.8, y = 1.90, yend = 2.15, color = "gray30", linewidth = 0.5) +
  annotate("segment", x = 23.0, xend = 23.0, y = 1.90, yend = 2.15, color = "gray30", linewidth = 0.5) +
  annotate("text", x = 26.4, y = 2.55, label = "Northwestern Hawaiian Islands", fontface = "bold", size = 2.8, color = "gray20") +
  
  # MHI Bracket Bar
  annotate("segment", x = 22.5, xend = 19.0, y = 2.15, yend = 2.15, color = "gray30", linewidth = 0.5) +
  annotate("segment", x = 22.5, xend = 22.5, y = 1.90, yend = 2.15, color = "gray30", linewidth = 0.5) +
  annotate("segment", x = 19.0, xend = 19.0, y = 1.90, yend = 2.15, color = "gray30", linewidth = 0.5) +
  annotate("text", x = 20.75, y = 2.55, label = "Main Hawaiian Islands", fontface = "bold", size = 2.8, color = "gray20") +
  
  scale_fill_gradientn(
    colors   = pal_lat_grad, 
    values   = pal_lat_vals, 
    limits   = c(19, 30), 
    na.value = "transparent", 
    guide    = "none"
  ) +
  scale_x_reverse(limits = x_limits, breaks = seq(19, 30, by = 1), expand = x_expand) +
  scale_y_continuous(limits = c(-0.5, 2.95), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  
  geom_point(data = island_locs_s5, aes(x = LATITUDE, y = 0, fill = LATITUDE), shape = 21, size = 1.8, color = "gray20", stroke = 0.5) +
  
  geom_text_repel(
    data = island_locs_s5, aes(x = LATITUDE, y = 0, label = PLOT_LABEL),
    nudge_x            = island_locs_s5$leg_nx,
    nudge_y            = island_locs_s5$leg_ny,
    angle              = 45,
    hjust              = 0,
    vjust              = 0.5,
    size               = 2.4,
    fontface           = "bold",
    color              = "gray20",
    segment.color      = "gray50",
    segment.size       = 0.35,
    min.segment.length = 0,
    force              = 0.5,
    box.padding        = 0.08,
    max.overlaps       = Inf
  ) +
  theme_void() +
  theme(
    axis.text.x        = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(t = 3)),
    axis.ticks.x       = element_line(color = "gray20", linewidth = 0.4),
    axis.ticks.length  = unit(0.12, "cm"),
    plot.margin        = margin(t = 4, r = 12, b = 6, l = 12, unit = "pt")
  )

# ------------------------------------------------------------------------------
# 5. FIGURE ASSEMBLY & EXPORT (Strict Viewport & Panel Height Locking)
# ------------------------------------------------------------------------------
# heights = c(1, 1) locks data panels A and B to identical heights
main_plots_s5 <- (p_urchin / p_fish) +
  plot_layout(heights = c(1, 1))

figS5_final <- (p_legend_s5 / main_plots_s5) +
  plot_layout(heights = c(0.20, 1))

out_png  <- file.path(supp_dir, "FigureS5_Bioerosion.png")
out_tiff <- file.path(supp_dir, "FigureS5_Bioerosion.tiff")

ggsave(out_png,  figS5_final, width = 17.8, height = 18.5, units = "cm", dpi = 600, bg = "white")
ggsave(out_tiff, figS5_final, width = 17.8, height = 18.5, units = "cm", dpi = 600, compression = "lzw", bg = "white")