# ==============================================================================
# Script:  Figure2_Grigg_vs_NCRMP.R
# Purpose: Historical (Grigg 1982) vs. Modern (NCRMP) P90 Carbonate Budgets
# Inputs:  data/grigg data/HA_1982_Grigg.csv
#          data/grigg data/HA_2024_StRS_Net_SiteLevel_GriggSubset.csv
# Outputs: outputs/main figures/Figure2_Grigg_vs_NCRMP.png
#          outputs/main figures/Figure2_Grigg_vs_NCRMP.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(ggpubr)
})

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

# Reversed Paul Tol Sunset (19N is Red, 30N is Blue)
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
# Helpers: Naming & IO
# ------------------------------------------------------------------------------
island_lookup <- c(
  "kure"                         = "H\u014dlanik\u016b",
  "kure atoll"                   = "H\u014dlanik\u016b",
  "midway"                       = "Kuaihelani",
  "midway atoll"                 = "Kuaihelani",
  "pearl & hermes"               = "Manawai",
  "pearl and hermes"             = "Manawai",
  "pearl & hermes atoll"         = "Manawai",
  "pearl and hermes atoll"       = "Manawai",
  "pearl_hermes"                 = "Manawai",
  "lisianski"                    = "Kapou",
  "lisianski island"             = "Kapou",
  "laysan"                       = "Kamole",
  "laysan island"                = "Kamole",
  "maro"                         = "Kamokuokamohoali'i",
  "maro reef"                    = "Kamokuokamohoali'i",
  "gardner"                      = "P\u016bh\u0101honu",
  "gardner pinnacles"            = "P\u016bh\u0101honu",
  "gardner_pinnacles"            = "P\u016bh\u0101honu",
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

saferead_csv <- function(candidates, dir = here("data", "grigg data")) {
  targets <- file.path(dir, candidates)
  found   <- targets[file.exists(targets)][1]
  if (is.na(found)) stop("None of candidate files found: ", paste(candidates, collapse = " | "))
  read_csv(found, show_col_types = FALSE)
}

boot_q90_se <- function(x, B = 1000) {
  x <- na.omit(x)
  n <- length(x)
  if (n <= 1) return(0)
  samples <- matrix(sample(x, size = n * B, replace = TRUE), nrow = n, ncol = B)
  sd(apply(samples, 2, quantile, probs = 0.90, names = FALSE))
}

theme_publication <- function(base_size = 9, y_margin_r = 2, x_margin_t = 4) {
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
      plot.margin       = margin(t = 5, r = 15, b = 5, l = 22, unit = "pt")
    )
}

# ------------------------------------------------------------------------------
# Data Import & Processing
# ------------------------------------------------------------------------------
df_grigg <- saferead_csv("HA_1982_Grigg.csv") %>%
  mutate(PLOT_LABEL = format_islands(LOCATION))

df_strs_sw <- saferead_csv(c("HA_2024_Grigg_Equivalent_Sites.csv", "HA_2024_StRS_Net_SiteLevel_GriggSubset.csv")) %>%
  mutate(
    LOCATION            = coalesce(!!!select(., any_of(c("LOCATION", "ISLAND")))),
    PLOT_LABEL          = format_islands(LOCATION),
    CORAL_COVER_PERC    = coalesce(!!!select(., any_of(c("CORAL_COVER_PERC", "HARD_CORAL_COVER_PCT_MEAN")))),
    GROSS_PROD_KG_M2_YR = coalesce(!!!select(., any_of(c("GROSS_PROD_KG_M2_YR", "GROSS_CARB_PROD_KG_M2_YR_MEAN"))))
  )

set.seed(42)
strs_sw_p90 <- df_strs_sw %>% 
  group_by(PLOT_LABEL) %>%
  summarize(
    LATITUDE                        = mean(LATITUDE, na.rm = TRUE),
    CORAL_COVER_P90                 = quantile(CORAL_COVER_PERC, probs = 0.90, na.rm = TRUE),
    CORAL_COVER_P90_SE              = boot_q90_se(CORAL_COVER_PERC),
    GROSS_CARB_PROD_KG_M2_YR_P90    = quantile(GROSS_PROD_KG_M2_YR, probs = 0.90, na.rm = TRUE),
    GROSS_CARB_PROD_KG_M2_YR_P90_SE = boot_q90_se(GROSS_PROD_KG_M2_YR),
    .groups                         = "drop"
  )

# Determine DATA_STATUS: 1982 Only, 2024 Only, or Shared
status_df <- full_join(
  df_grigg %>% distinct(PLOT_LABEL) %>% mutate(in_82 = TRUE),
  strs_sw_p90 %>% distinct(PLOT_LABEL) %>% mutate(in_24 = TRUE),
  by = "PLOT_LABEL"
) %>%
  mutate(DATA_STATUS = case_when(
    in_82 & in_24 ~ "Shared",
    in_82 & is.na(in_24) ~ "1982 Only",
    is.na(in_82) & in_24 ~ "2024 Only"
  )) %>%
  select(PLOT_LABEL, DATA_STATUS)

# Attach status to dataframes
df_grigg <- df_grigg %>% left_join(status_df, by = "PLOT_LABEL")
strs_sw_p90 <- strs_sw_p90 %>% left_join(status_df, by = "PLOT_LABEL")

shared_p90_deltas <- strs_sw_p90 %>% 
  inner_join(
    df_grigg %>% select(PLOT_LABEL, GRIGG_LAT = LATITUDE, GRIGG_COV = CORAL_COVER_PERC, GRIGG_PROD = GROSS_PROD),
    by = "PLOT_LABEL"
  ) %>%
  mutate(
    LATITUDE       = GRIGG_LAT,
    DELTA_COV_P90  = CORAL_COVER_P90 - GRIGG_COV,
    DELTA_COV_SE   = CORAL_COVER_P90_SE,
    DELTA_PROD_P90 = GROSS_CARB_PROD_KG_M2_YR_P90 - GRIGG_PROD,
    DELTA_PROD_SE  = GROSS_CARB_PROD_KG_M2_YR_P90_SE
  )

island_locs_f2 <- bind_rows(
  select(df_grigg, PLOT_LABEL, LATITUDE, DATA_STATUS),
  select(strs_sw_p90, PLOT_LABEL, LATITUDE, DATA_STATUS)
) %>%
  filter(!is.na(LATITUDE) & !is.na(PLOT_LABEL)) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

# ------------------------------------------------------------------------------
# Label Position Offsets for Panels
# ------------------------------------------------------------------------------
nudge_lookup <- tribble(
  ~PLOT_LABEL,        ~cov_nx, ~cov_ny, ~prod_nx, ~prod_ny,
  "H\u014dlanik\u016b",   2.5,    -2.0,      1.2,     -2.5,
  "Manawai",             -1.0,   -12.0,     -1.0,      2.8,
  "Kapou",               -1.2,     8.0,     -1.2,      1.6,
  "Lalo",                -1.4,   -12.0,     -1.6,     -1.3,
  "Kaua'i",               1.0,     8.0,      1.0,      2.0,
  "Ni'ihau",              1.2,    -8.0,      1.2,     -2.5,
  "O'ahu",               -1.2,     6.0,     -1.2,     -2.8,
  "Moloka'i",            -1.0,     8.0,     -1.0,      2.5,
  "L\u0101na'i",         -1.2,    -8.0,     -1.2,     -2.5,
  "Maui",                 1.6,     8.0,      1.6,      1.0,
  "Hawai'i",             -1.6,   -13.5,     -1.2,     -3.0
)

delta_cov_data <- shared_p90_deltas %>%
  left_join(nudge_lookup %>% select(PLOT_LABEL, nudge_x = cov_nx, nudge_y = cov_ny), by = "PLOT_LABEL") %>%
  mutate(nudge_x = replace_na(nudge_x, 1.0), nudge_y = replace_na(nudge_y, 0.0))

delta_prod_data <- shared_p90_deltas %>%
  left_join(nudge_lookup %>% select(PLOT_LABEL, nudge_x = prod_nx, nudge_y = prod_ny), by = "PLOT_LABEL") %>%
  mutate(nudge_x = replace_na(nudge_x, 1.0), nudge_y = replace_na(nudge_y, 0.0))

# ------------------------------------------------------------------------------
# Plot Construction
# ------------------------------------------------------------------------------
build_plot_f2 <- function(data, y_var, y_se = NULL, y_label, y_lim, y_brk, 
                          add_zero = FALSE, add_lbl = FALSE, show_x = TRUE) {
  p <- ggplot(data, aes(x = LATITUDE, y = .data[[y_var]])) +
    scale_x_reverse(limits = c(30.2, 18.6), breaks = seq(19, 30, by = 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = y_brk) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme_publication()
  
  if (!show_x) {
    p <- p + labs(x = NULL) + theme(axis.title.x = element_blank())
  } else {
    p <- p + labs(x = "Latitude (°N)")
  }
  
  if (add_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6)
  }
  
  if (!is.null(y_se)) {
    p <- p + geom_errorbar(aes(ymin = .data[[y_var]] - .data[[y_se]], 
                               ymax = .data[[y_var]] + .data[[y_se]]), 
                           color = "gray55", width = 0, linewidth = 0.5)
  }
  
  p <- p + geom_point(aes(fill = LATITUDE), shape = 21, size = 2.4, color = "gray20", stroke = 0.5, alpha = 0.85) +
    scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
    labs(y = y_label)
  
  if (add_lbl) {
    p <- p + geom_text_repel(
      aes(label = PLOT_LABEL), size = 2.4, fontface = "bold", color = "gray20",
      box.padding = 0.25, point.padding = 0.35, 
      nudge_x = data$nudge_x, nudge_y = data$nudge_y,
      force = 1, seed = 42, segment.color = "gray60", segment.alpha = 0.8,
      segment.size = 0.35, min.segment.length = 0,
      max.overlaps = Inf
    )
  }
  p
}

lbl_cov  <- "Coral cover (%)"
lbl_prod <- expression(bold("Gross production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))

grigg_cover <- build_plot_f2(df_grigg, "CORAL_COVER_PERC", y_label = lbl_cov, y_lim = c(0, 100), y_brk = seq(0, 100, 20), show_x = FALSE)
ncrmp_cover <- build_plot_f2(strs_sw_p90, "CORAL_COVER_P90", "CORAL_COVER_P90_SE", y_label = lbl_cov, y_lim = c(0, 100), y_brk = seq(0, 100, 20), show_x = FALSE)
delta_cover <- build_plot_f2(delta_cov_data, "DELTA_COV_P90", "DELTA_COV_SE", y_label = "\u0394 Coral cover (%)", y_lim = c(-80, 40), y_brk = seq(-80, 40, 20), add_zero = TRUE, add_lbl = TRUE, show_x = FALSE)

grigg_prod  <- build_plot_f2(df_grigg, "GROSS_PROD", y_label = lbl_prod, y_lim = c(0, 20), y_brk = seq(0, 20, 5))
ncrmp_prod  <- build_plot_f2(strs_sw_p90, "GROSS_CARB_PROD_KG_M2_YR_P90", "GROSS_CARB_PROD_KG_M2_YR_P90_SE", y_label = lbl_prod, y_lim = c(0, 20), y_brk = seq(0, 20, 5))
delta_prod  <- build_plot_f2(delta_prod_data, "DELTA_PROD_P90", "DELTA_PROD_SE", y_label = expression(bold(Delta ~ "Gross production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), y_lim = c(-15, 10), y_brk = seq(-15, 10, 5), add_zero = TRUE, add_lbl = TRUE)

main_plots_f2 <- ggarrange(
  grigg_cover, ncrmp_cover, delta_cover,
  grigg_prod,  ncrmp_prod,  delta_prod,
  ncol = 3, nrow = 2, align = "hv", labels = "AUTO",
  font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02
)

# ------------------------------------------------------------------------------
# Top Legend Construction (Aligned to 19-30 Latitude Range)
# ------------------------------------------------------------------------------
legend_nudges <- tribble(
  ~PLOT_LABEL,             ~leg_nx, ~leg_ny,
  "H\u014dlanik\u016b",  -0.3,    0.6,
  "Kuaihelani",          -0.1,    0.6,
  "Manawai",             -0.1,    0.6,
  "Kapou",               -0.3,    0.6,
  "Kamole",              -0.2,    0.6,
  "Kamokuokamohoali'i",  -0.2,    0.6,
  "P\u016bh\u0101honu",  -0.1,    0.6,
  "Lalo",                -0.3,    0.6,
  "Mokumanamana",        -0.2,    0.6,
  "Nihoa",               -0.2,    0.6,
  "Kaua'i",              -0.6,    0.6,
  "Ni'ihau",             -0.3,    0.6,
  "O'ahu",               -0.4,    0.6,
  "Moloka'i",            -0.2,    0.6,
  "L\u0101na'i",          0.0,    0.6,
  "Maui",                -0.3,    0.6,
  "Kaho'olawe",           0.0,    0.6,
  "Hawai'i",              0,      0.6
)

island_locs_f2 <- island_locs_f2 %>%
  left_join(legend_nudges, by = "PLOT_LABEL") %>%
  mutate(leg_nx = replace_na(leg_nx, 0), leg_ny = replace_na(leg_ny, 0.5))

p_legend_f2 <- ggplot() +
  geom_tile(
    data = data.frame(lat = seq(19, 30, length.out = 500), y = 0),
    aes(x = lat, y = y, fill = lat),
    height = 0.5
  ) +
  annotate("rect", xmin = 19, xmax = 30, ymin = -0.25, ymax = 0.25, color = "gray20", fill = NA, linewidth = 0.5) +
  
  # NWHI Bracket Bar 
  annotate("segment", x = 30, xend = 22.8, y = 4.8, yend = 4.8, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 30, xend = 30, y = 4.4, yend = 4.8, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 22.8, xend = 22.8, y = 4.4, yend = 4.8, color = "gray30", linewidth = 0.6) +
  annotate("text", x = 26.4, y = 5.2, label = "Northwestern Hawaiian Islands", fontface = "bold", size = 3.8, color = "gray20") +
  
  # MHI Bracket Bar 
  annotate("segment", x = 22.7, xend = 19.0, y = 4.8, yend = 4.8, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 22.7, xend = 22.7, y = 4.4, yend = 4.8, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 19.0, xend = 19.0, y = 4.4, yend = 4.8, color = "gray30", linewidth = 0.6) +
  annotate("text", x = 20.8, y = 5.2, label = "Main Hawaiian Islands", fontface = "bold", size = 3.8, color = "gray20") +
  
  annotate("text", x = 30.2, y = 0, label = "Latitude (°N)", fontface = "bold", size = 3.5, hjust = 1, color = "black") +
  
  scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
  scale_x_reverse(limits = c(30.2, 18.6), breaks = seq(19, 30, by = 1), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-0.6, 5.8), expand = c(0, 0)) + 
  coord_cartesian(clip = "off") +
  
  geom_point(data = island_locs_f2, aes(x = LATITUDE, y = 0, fill = LATITUDE), shape = 21, size = 2.0, color = "gray20", stroke = 0.6) +
  
  geom_text_repel(
    data = island_locs_f2, aes(x = LATITUDE, y = 0, label = PLOT_LABEL),
    nudge_x = island_locs_f2$leg_nx,
    nudge_y = island_locs_f2$leg_ny,
    angle = 45,             
    hjust = 0,               
    vjust = 0.5,               
    size = 2.8,
    fontface = "bold",
    color = "gray20",
    segment.color = "gray50",
    segment.size = 0.4,
    min.segment.length = 0,
    force = 0.5,
    box.padding = 0.1,
    max.overlaps = Inf
  ) +
  theme_void() +
  theme(
    axis.text.x       = element_text(size = 8, face = "bold", color = "gray20", margin = margin(t = 4)),
    axis.ticks.x      = element_line(color = "gray20"),
    axis.ticks.length = unit(0.15, "cm"),
    plot.margin       = margin(t = 10, r = 5, b = 5, l = 75, unit = "pt")
  )

# ------------------------------------------------------------------------------
# Figure Assembly & Save
# ------------------------------------------------------------------------------
# Increased legend height proportion from 0.35 to 0.42 to accommodate raised brackets
fig2_final <- ggarrange(p_legend_f2, main_plots_f2, nrow = 2, heights = c(0.42, 1))

ggsave(file.path(main_dir, "Figure2_Grigg_vs_NCRMP.png"), fig2_final, width = 28, height = 16.0, units = "cm", dpi = 300, bg = "white") 
ggsave(file.path(main_dir, "Figure2_Grigg_vs_NCRMP.tiff"), fig2_final, width = 28, height = 16.0, units = "cm", dpi = 300, compression = "lzw", bg = "white")