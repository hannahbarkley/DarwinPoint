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

pal_lat_grad <- c("#D73027", "#F46D43", "#FDAE61", "#FFFFBF", "#74ADD1", "#4575B4")

# ------------------------------------------------------------------------------
# Helpers: Naming & IO
# ------------------------------------------------------------------------------
# Full Hawaiian Archipelago dictionary mapping English/codes to Hawaiian names
island_lookup <- c(
  # Northwestern Hawaiian Islands (Papahānaumokuākea)
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
  
  # Main Hawaiian Islands
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

# Vectorized bootstrap for standard error of 90th percentile
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
      plot.margin       = margin(t = 10, r = 4, b = 2, l = 22, unit = "pt")
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

# Retain ALL unique islands across both datasets for the sidebar legend
island_locs_f2 <- bind_rows(
  select(df_grigg, PLOT_LABEL, LATITUDE),
  select(strs_sw_p90, PLOT_LABEL, LATITUDE)
) %>%
  filter(!is.na(LATITUDE) & !is.na(PLOT_LABEL)) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

# ------------------------------------------------------------------------------
# Label Position Offsets
# ------------------------------------------------------------------------------
nudge_lookup <- tribble(
  ~PLOT_LABEL,        ~cov_nx, ~cov_ny, ~prod_nx, ~prod_ny,
  "H\u014dlanik\u016b",   2.5,    -2.0,      1.2,     -2.5,
  "Manawai",             -1.0,   -12.0,     -1.0,      2.5,
  "Kapou",               -2.2,     8.0,     -2.2,      1.6,
  "Lalo",                -1.4,   -12.0,     -1.6,     -1.3,
  "Kaua'i",               1.0,     8.0,      1.0,      2.0,
  "Ni'ihau",              1.2,    -8.0,      1.2,     -2.5,
  "O'ahu",               -2.0,     6.0,     -1.2,     -2.8,
  "Moloka'i",            -1.0,     8.0,     -1.0,      2.5,
  "L\u0101na'i",         -1.2,    -8.0,     -1.2,     -2.5,
  "Maui",                 1.6,     8.0,      1.6,      2.4,
  "Hawai'i",              1.8,    14.0,      1.6,      2.5
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
    scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
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
  
  p <- p + geom_point(aes(fill = LATITUDE), shape = 21, color = "gray30", stroke = 0.5, size = 2.8, alpha = 0.85) +
    scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
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
  ncol = 3, nrow = 2, align = "hv", labels = "auto",
  font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02
)

# Right Sidebar Legend: plots all Hawaiian names with ample horizontal clearance
p_legend_f2 <- ggplot() +
  geom_tile(
    data = data.frame(lat = seq(18, 30, length.out = 500), x = 0),
    aes(x = x, y = lat, fill = lat),
    width = 0.5
  ) +
  annotate("rect", xmin = -0.25, xmax = 0.25, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  scale_x_continuous(limits = c(-0.25, 4.8), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(18, 30, by = 2), expand = c(0.02, 0.02)) +
  coord_cartesian(ylim = c(18, 30), clip = "off") +
  geom_point(data = island_locs_f2, aes(x = 0.25, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2, color = "gray30", stroke = 0.5) +
  geom_text_repel(
    data = island_locs_f2, aes(x = 0.25, y = LATITUDE, label = PLOT_LABEL),
    direction = "y", hjust = 0, nudge_x = 0.45, segment.color = "gray50",
    segment.size = 0.4, size = 2.4, fontface = "bold", color = "gray20", force = 3,
    max.overlaps = Inf
  ) +
  theme_void() +
  theme(
    plot.title        = element_text(size = 8, face = "bold", margin = margin(b = 4)),
    axis.text.y       = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 5)),
    axis.ticks.y      = element_line(color = "gray20"),
    axis.ticks.length = unit(0.15, "cm"),
    plot.margin       = margin(20, 0, 20, 10)
  ) +
  labs(title = "Latitude (°N)")

# ------------------------------------------------------------------------------
# Figure Assembly & Save
# ------------------------------------------------------------------------------
fig2_final <- ggarrange(main_plots_f2, p_legend_f2, ncol = 2, widths = c(1, 0.25))

ggsave(file.path(main_dir, "Figure2_Grigg_vs_NCRMP.png"), fig2_final, width = 28, height = 12.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(main_dir, "Figure2_Grigg_vs_NCRMP.tiff"), fig2_final, width = 28, height = 12.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")
cat("✅ Figure 2 saved successfully.\n")