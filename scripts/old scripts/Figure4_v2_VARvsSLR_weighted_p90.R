# ==============================================================================
# Script:  Figure4_VARvsSLR_weighted_p90.R
# Purpose: Vertical Reef Accretion Rate (VAR) 90th Percentiles vs. 2050 SLR.
#          Calculates the 90th percentile (P90) of VAR across sites using the 
#          weighted dataset, applies bootstrap SEs (as in script 4b), and 
#          plots island-level P90s against IPCC sea-level rise scenarios.
# Inputs:  - data/strs data/processed data/HA_2024_StRS_SiteLevel_Weighted.csv
#          - data/supporting data/IPCC_AR6_SLR_Rates_2040_2060.csv
# Outputs: - outputs/main figures/Figure4_VARvsSLR_weighted_p90.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
})

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# Helpers: Bootstrapping, Names & Themes
# ------------------------------------------------------------------------------
boot_q90_se <- function(x, n_boot = 1000) {
  x_clean <- na.omit(x)
  if (length(x_clean) < 2) return(0)
  boot_reps <- replicate(n_boot, {
    s <- sample(x_clean, size = length(x_clean), replace = TRUE)
    quantile(s, probs = 0.90, na.rm = TRUE, names = FALSE)
  })
  sd(boot_reps)
}

island_lookup <- c(
  "kure"                  = "Hōlanikū",
  "midway"                = "Kuaihelani",
  "pearl & hermes"        = "Manawai",
  "pearl and hermes"      = "Manawai",
  "lisianski"             = "Kapou",
  "laysan"                = "Kamole",
  "maro"                  = "Kamokuokamohoali'i",
  "gardner"               = "Pūhāhonu",
  "french frigate"        = "Lalo",
  "french frigate shoals" = "Lalo",
  "ffs"                   = "Lalo",
  "necker"                = "Mokumanamana",
  "nihoa"                 = "Nihoa",
  "kauai"                 = "Kaua'i",
  "niihau"                = "Ni'ihau",
  "oahu"                  = "O'ahu",
  "molokai"               = "Moloka'i",
  "lanai"                 = "Lāna'i",
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
      plot.margin       = margin(t = 12, r = 12, b = 4, l = 16, unit = "pt")
    )
}

saferead_csv <- function(candidates) {
  search_dirs <- c(here("data", "supporting data"), here("data", "strs data", "processed data"), here("data"), here())
  for (f in candidates) {
    for (d in search_dirs) {
      target <- file.path(d, f)
      if (file.exists(target)) return(read_csv(target, show_col_types = FALSE))
    }
  }
  stop(paste("None of the candidate files found:", paste(candidates, collapse = " | ")))
}

# ------------------------------------------------------------------------------
# 1. INGEST & PREPARE DATA
# ------------------------------------------------------------------------------
# SLR Projections
slr_raw <- saferead_csv(c("IPCC_AR6_SLR_Rates_2040_2060.csv", "ipcc_ar6_sea_level_projection_psmsl_id_523.csv"))

if ("Median_q50" %in% names(slr_raw)) {
  df_slr <- slr_raw %>% filter(!grepl("119|1-1.9", Scenario))
} else {
  df_slr <- slr_raw %>%
    filter(process == "totalrates", confidence == "medium", scenario != "ssp119") %>%
    select(scenario, quantile, rate = `2050`) %>%
    pivot_wider(names_from = quantile, values_from = rate, names_prefix = "q") %>%
    mutate(
      Scenario = case_match(scenario, "ssp126" ~ "SSP1-2.6", "ssp245" ~ "SSP2-4.5", "ssp370" ~ "SSP3-7.0", "ssp585" ~ "SSP5-8.5", .default = scenario),
      Median_q50 = q50, Lower_17th_q17 = q17, Upper_83rd_q83 = q83
    )
}

slr_box_min <- min(df_slr$Lower_17th_q17, na.rm = TRUE)
slr_box_max <- max(df_slr$Upper_83rd_q83, na.rm = TRUE)

rate_ssp126 <- df_slr %>% filter(grepl("126|1-2.6", Scenario)) %>% pull(Median_q50)
rate_ssp245 <- df_slr %>% filter(grepl("245|2-4.5", Scenario)) %>% pull(Median_q50)
rate_ssp370 <- df_slr %>% filter(grepl("370|3-7.0", Scenario)) %>% pull(Median_q50)
rate_ssp585 <- df_slr %>% filter(grepl("585|5-8.5", Scenario)) %>% pull(Median_q50)

# Site-Level Data & P90 Aggregation
df_site <- saferead_csv("HA_2024_StRS_SiteLevel_Weighted.csv") %>%
  mutate(
    LOCATION = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                        if ("ISLAND" %in% names(.)) ISLAND else NULL)
  )

set.seed(42)
df_strs_p90 <- df_site %>%
  group_by(LOCATION) %>%
  summarize(
    LATITUDE                    = mean(LATITUDE, na.rm = TRUE),
    VERTICAL_ACCRETION_MM_YR    = quantile(VERTICAL_ACCRETION_MM_YR, probs = 0.90, na.rm = TRUE, names = FALSE),
    se.VERTICAL_ACCRETION_MM_YR = boot_q90_se(VERTICAL_ACCRETION_MM_YR),
    .groups = "drop"
  ) %>%
  mutate(
    PLOT_LABEL = format_islands(LOCATION),
    is_kapou   = ifelse(PLOT_LABEL == "Kapou", "Kapou", "Other"),
    nudge_x = case_when(
      PLOT_LABEL == "Hawai'i" ~ 0.8, PLOT_LABEL == "Kaho'olawe" ~ 0.5, PLOT_LABEL == "Maui" ~ 1.0, 
      PLOT_LABEL == "Lāna'i" ~ -0.2, PLOT_LABEL == "Moloka'i" ~ -1.0, PLOT_LABEL == "O'ahu" ~ 0.3, 
      PLOT_LABEL == "Kaua'i" ~ -0.8, PLOT_LABEL == "Ni'ihau" ~ -0.6, PLOT_LABEL == "Lalo" ~ -0.8, 
      PLOT_LABEL == "Kapou" ~ -0.5, PLOT_LABEL == "Manawai" ~ 0.2, PLOT_LABEL == "Hōlanikū" ~ 0.8, 
      TRUE ~ 0
    ),
    nudge_y = case_when(
      PLOT_LABEL == "Hawai'i" ~ 0.5, PLOT_LABEL == "Kaho'olawe" ~ 2.5, PLOT_LABEL == "Maui" ~ -0.2, 
      PLOT_LABEL == "Lāna'i" ~ 1.5, PLOT_LABEL == "Moloka'i" ~ 0.8, PLOT_LABEL == "O'ahu" ~ -1.5, 
      PLOT_LABEL == "Kaua'i" ~ 0.8, PLOT_LABEL == "Ni'ihau" ~ -1.2, PLOT_LABEL == "Lalo" ~ 1.0, 
      PLOT_LABEL == "Kapou" ~ 1.5, PLOT_LABEL == "Manawai" ~ 1.5, PLOT_LABEL == "Hōlanikū" ~ -1.5, 
      TRUE ~ 0
    )
  )

# ------------------------------------------------------------------------------
# 2. FIT TRENDLINE (Full Dataset Only)
# ------------------------------------------------------------------------------
fit_full <- lm(VERTICAL_ACCRETION_MM_YR ~ LATITUDE, data = df_strs_p90)
lats     <- seq(min(df_strs_p90$LATITUDE, na.rm=T), max(df_strs_p90$LATITUDE, na.rm=T), length.out = 100)
pred_full<- predict(fit_full, newdata = data.frame(LATITUDE = lats), se.fit = TRUE)

df_trend_full <- data.frame(
  LATITUDE = lats, 
  fit      = pred_full$fit, 
  ymin     = pred_full$fit - 1.96 * pred_full$se.fit, 
  ymax     = pred_full$fit + 1.96 * pred_full$se.fit
)

# ------------------------------------------------------------------------------
# 3. PLOT CONSTRUCTION
# ------------------------------------------------------------------------------
pal_lat_grad <- c("#D73027", "#F46D43", "#FDAE61", "#FFFFBF", "#74ADD1", "#4575B4")

col_gray_ssp126 <- "#71717A"
col_gray_ssp245 <- "#52525B"
col_gray_ssp370 <- "#3F3F46"
col_gray_ssp585 <- "#18181B"

fig4_final <- ggplot(df_strs_p90, aes(x = LATITUDE, y = VERTICAL_ACCRETION_MM_YR)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = slr_box_min, ymax = slr_box_max, fill = "#F4F4F5", alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.55) +
  geom_hline(yintercept = rate_ssp585, linetype = "dotted", color = col_gray_ssp585, linewidth = 0.6) +
  geom_hline(yintercept = rate_ssp370, linetype = "dotted", color = col_gray_ssp370, linewidth = 0.6) +
  geom_hline(yintercept = rate_ssp245, linetype = "dotted", color = col_gray_ssp245, linewidth = 0.6) +
  geom_hline(yintercept = rate_ssp126, linetype = "dotted", color = col_gray_ssp126, linewidth = 0.6) +
  
  annotate("text", x = 29.8, y = 11.0, label = sprintf("2050 Sea Level Rise Projections for Hawai'i (%.1f\u2013%.1f mm yr\u207B\u00B9)", slr_box_min, slr_box_max), color = "#18181B", size = 2.8, fontface = "bold", hjust = 0) +
  
  annotate("label", x = 29.8, y = rate_ssp585, label = "SSP5-8.5", color = col_gray_ssp585, size = 2.4, fontface = "bold", hjust = 0, fill = "#F4F4F5", label.size = NA, label.padding = unit(0.12, "lines")) +
  annotate("label", x = 29.8, y = rate_ssp370, label = "SSP3-7.0", color = col_gray_ssp370, size = 2.4, fontface = "bold", hjust = 0, fill = "#F4F4F5", label.size = NA, label.padding = unit(0.12, "lines")) +
  annotate("label", x = 29.8, y = rate_ssp245, label = "SSP2-4.5", color = col_gray_ssp245, size = 2.4, fontface = "bold", hjust = 0, fill = "#F4F4F5", label.size = NA, label.padding = unit(0.12, "lines")) +
  annotate("label", x = 29.8, y = rate_ssp126, label = "SSP1-2.6", color = col_gray_ssp126, size = 2.4, fontface = "bold", hjust = 0, fill = "#F4F4F5", label.size = NA, label.padding = unit(0.12, "lines")) +
  
  # Trendline (Full dataset only)
  geom_ribbon(data = df_trend_full, aes(x = LATITUDE, ymin = ymin, ymax = ymax), fill = "#64748B", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = df_trend_full, aes(x = LATITUDE, y = fit), color = "#334155", linewidth = 0.75, linetype = "solid", inherit.aes = FALSE) +
  
  # Points & Error bars
  geom_errorbar(aes(ymin = VERTICAL_ACCRETION_MM_YR - se.VERTICAL_ACCRETION_MM_YR, ymax = VERTICAL_ACCRETION_MM_YR + se.VERTICAL_ACCRETION_MM_YR), color = "gray55", width = 0, linewidth = 0.45) +
  geom_point(aes(fill = LATITUDE, shape = is_kapou), color = "gray20", stroke = 0.5, size = 3.4, alpha = 0.95) +
  scale_shape_manual(values = c("Kapou" = 24, "Other" = 21), guide = "none") +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  
  geom_text_repel(aes(label = PLOT_LABEL), nudge_x = df_strs_p90$nudge_x, nudge_y = df_strs_p90$nudge_y, size = 2.8, fontface = "bold", color = "black", box.padding = 0, point.padding = 0.2, force = 0.5, max.overlaps = Inf, seed = 42, segment.color = "gray40", segment.alpha = 0.75, segment.linewidth = 0.35, min.segment.length = 0) +
  
  # Capped Y-Axis Limits
  scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 2)) +
  coord_cartesian(ylim = c(-4.5, 12.8), clip = "off") +
  theme_publication() +
  theme(plot.margin = margin(12, 16, 10, 14, "pt"), axis.title.y = element_text(margin = margin(r = 4)), axis.title.x = element_text(margin = margin(t = 6))) +
  labs(x = "Latitude (°N)", y = expression(bold("Vertical accretion rate 90th percentile (mm" ~ yr^{ -1 } * ")")))

# ------------------------------------------------------------------------------
# 4. EXPORT
# ------------------------------------------------------------------------------
ggsave(file.path(main_dir, "Figure4_VARvsSLR_weighted_p90.png"), fig4_final, width = 16.0, height = 11.5, units = "cm", dpi = 300, bg = "white")
cat("✅ Figure 4 P90 version saved successfully using weighted data (capped at 12)!\n")