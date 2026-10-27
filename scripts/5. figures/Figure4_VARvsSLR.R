# ==============================================================================
# Figure4_VARvsSLR.R
# Vertical Accretion Rate (VAR) vs. 2050 Sea-Level Rise (SLR)
# Formatted to PNAS 2-column specifications (17.8 cm width, 600 DPI)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(survey)
  library(scales)
})

options(survey.lonely.psu = "adjust")

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. COLOR PALETTES & STYLING
# ------------------------------------------------------------------------------

# Reversed Paul Tol Sunset: 19N = red, 30N = blue
pal_lat_grad <- c(
  "#A50026", "#DD3D2D", "#F67E4B", "#FDB366", "#FEDA8B",
  "#EAECCC", "#C2E4EF", "#98CAE1", "#6EA6CD", "#4A7BB7", "#364B9A"
)

pal_lat_vals <- c(
  seq(0, 0.3409, length.out = 6)[1:5],
  0.3409,
  seq(0.3409, 1, length.out = 6)[2:6]
)

# Scenario & element styling
col_ssp126       <- "#52525B"
col_ssp245       <- "#52525B"
col_ssp370       <- "#52525B"
col_ssp585       <- "#52525B"

col_site_outline <- "#52525B"
col_mean_outline <- "#09090B"
col_trend        <- "#0F172A"
col_trend_band   <- "#64748B"
slr_band_fill    <- "#D8D1E8"

# ------------------------------------------------------------------------------
# 2. ISLAND NOMENCLATURE
# ------------------------------------------------------------------------------

island_lookup <- c(
  "kure"                  = "Hōlanikū",
  "kure atoll"            = "Hōlanikū",
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
  "gardner"               = "Pūhāhonu",
  "gardner pinnacles"     = "Pūhāhonu",
  "gardner_pinnacles"     = "Pūhāhonu",
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
  "lanai"                 = "Lāna'i",
  "maui"                  = "Maui",
  "kahoolawe"             = "Kaho'olawe",
  "hawaii"                = "Hawai'i"
)

format_islands <- function(x) {
  clean_key <- tolower(trimws(as.character(x)))
  clean_key <- gsub("_", " ", clean_key)
  mapped <- island_lookup[clean_key]
  coalesce(unname(mapped), as.character(x))
}

# ------------------------------------------------------------------------------
# 3. PNAS THEME
# ------------------------------------------------------------------------------

theme_pnas <- function(base_size = 8, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      panel.background  = element_rect(fill = "#FFFFFF", color = NA),
      plot.background   = element_rect(fill = "#FFFFFF", color = NA),
      axis.line         = element_line(color = "#18181B", linewidth = 0.5),
      axis.ticks        = element_line(color = "#18181B", linewidth = 0.5),
      axis.ticks.length = unit(0.12, "cm"),
      axis.text         = element_text(color = "#27272A", size = 7),
      axis.title        = element_text(color = "#09090B", size = 8.5, face = "bold"),
      axis.title.y      = element_text(margin = margin(r = 5)),
      axis.title.x      = element_text(margin = margin(t = 5)),
      plot.margin       = margin(t = 8, r = 10, b = 6, l = 8, unit = "pt")
    )
}

# ------------------------------------------------------------------------------
# 4. DATA INGESTION
# ------------------------------------------------------------------------------

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "supporting data"),
    here("data", "strs data", "processed data"),
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

slr_raw <- saferead_csv(c(
  "IPCC_AR6_SLR_Rates_2040_2060.csv",
  "ipcc_ar6_sea_level_projection_psmsl_id_523.csv"
))

if ("Median_q50" %in% names(slr_raw)) {
  df_slr <- slr_raw %>% filter(!grepl("119|1-1.9", Scenario))
} else {
  df_slr <- slr_raw %>%
    filter(process == "totalrates", confidence == "medium", scenario != "ssp119") %>%
    select(scenario, quantile, rate = `2050`) %>%
    pivot_wider(names_from = quantile, values_from = rate, names_prefix = "q") %>%
    mutate(
      Scenario = case_match(
        scenario,
        "ssp126" ~ "SSP1-2.6",
        "ssp245" ~ "SSP2-4.5",
        "ssp370" ~ "SSP3-7.0",
        "ssp585" ~ "SSP5-8.5",
        .default = scenario
      ),
      Median_q50     = q50,
      Lower_17th_q17 = q17,
      Upper_83rd_q83 = q83
    )
}

slr_box_min <- min(df_slr$Lower_17th_q17, na.rm = TRUE)
slr_box_max <- max(df_slr$Upper_83rd_q83, na.rm = TRUE)

rate_ssp126 <- df_slr %>% filter(grepl("126|1-2.6", Scenario)) %>% pull(Median_q50) %>% first()
rate_ssp245 <- df_slr %>% filter(grepl("245|2-4.5", Scenario)) %>% pull(Median_q50) %>% first()
rate_ssp370 <- df_slr %>% filter(grepl("370|3-7.0", Scenario)) %>% pull(Median_q50) %>% first()
rate_ssp585 <- df_slr %>% filter(grepl("585|5-8.5", Scenario)) %>% pull(Median_q50) %>% first()

# Island-level data
df_strs_f4 <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>%
  mutate(
    LOCATION = coalesce(
      if ("LOCATION" %in% names(.)) LOCATION else NULL,
      if ("ISLAND" %in% names(.)) ISLAND else NULL
    ),
    PLOT_LABEL = format_islands(LOCATION),
    VERTICAL_ACCRETION_MM_YR = coalesce(
      if ("VERTICAL_ACCRETION_MM_YR" %in% names(.)) VERTICAL_ACCRETION_MM_YR else NULL, 0
    ),
    VISUAL_LAT = case_when(
      PLOT_LABEL == "Lāna'i"  ~ LATITUDE + 0.14,
      PLOT_LABEL == "Maui"    ~ LATITUDE - 0.05,
      PLOT_LABEL == "Kaua'i"  ~ LATITUDE + 0.08,
      PLOT_LABEL == "Ni'ihau" ~ LATITUDE - 0.08,
      TRUE ~ LATITUDE
    )
  )

# Site-level data
df_strs_site <- saferead_csv("HA_2024_StRS_SiteLevel_Weighted.csv") %>%
  mutate(
    LOCATION = coalesce(
      if ("LOCATION" %in% names(.)) LOCATION else NULL,
      if ("ISLAND" %in% names(.)) ISLAND else NULL
    ),
    PLOT_LABEL = format_islands(LOCATION),
    VERTICAL_ACCRETION_MM_YR = coalesce(
      if ("VERTICAL_ACCRETION_MM_YR" %in% names(.)) VERTICAL_ACCRETION_MM_YR else NULL, 0
    )
  ) %>%
  select(-any_of("LATITUDE")) %>%
  left_join(df_strs_f4 %>% select(PLOT_LABEL, LATITUDE), by = "PLOT_LABEL") %>%
  filter(!is.na(LATITUDE)) %>%
  mutate(
    VISUAL_LAT = case_when(
      PLOT_LABEL == "Lāna'i"  ~ LATITUDE + 0.14,
      PLOT_LABEL == "Maui"    ~ LATITUDE - 0.09,
      PLOT_LABEL == "Kaua'i"  ~ LATITUDE + 0.08,
      PLOT_LABEL == "Ni'ihau" ~ LATITUDE - 0.08,
      TRUE ~ LATITUDE
    )
  )

# ------------------------------------------------------------------------------
# 5. SURVEY DESIGN & GLM MODEL
# ------------------------------------------------------------------------------

des_file <- here("data", "strs data", "processed data", "HA_2024_StRS_SurveyDesign.rds")
des <- if (file.exists(des_file)) readRDS(des_file) else NULL

if (!is.null(des)) {
  mean_lat <- mean(des$variables$LOCATION_LAT, na.rm = TRUE)
  des <- update(des, LOCATION_LAT_c = LOCATION_LAT - mean_lat)
}

get_model_fit <- function(y_var, des_obj, df_island, mean_lat_ref = mean_lat) {
  if (is.null(des_obj) || !y_var %in% names(des_obj$variables)) return(NULL)
  
  mhi_min  <- min(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  mhi_max  <- max(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  nwhi_min <- min(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  nwhi_max <- max(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  
  mhi_lats  <- seq(mhi_min, mhi_max, length.out = 150)
  nwhi_lats <- seq(nwhi_min, nwhi_max, length.out = 180)
  reg_levels <- if (is.factor(des_obj$variables$REGIONCODE)) levels(des_obj$variables$REGIONCODE) else unique(des_obj$variables$REGIONCODE)
  
  fit_inter <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c * REGIONCODE")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  p_inter <- 1
  if (!is.null(fit_inter)) p_inter <- as.numeric(regTermTest(fit_inter, ~LOCATION_LAT_c:REGIONCODE)$p)
  
  if (!is.na(p_inter) && p_inter < 0.05) {
    pred_mhi  <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = mhi_lats - mean_lat_ref, REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = nwhi_lats - mean_lat_ref, REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    return(bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))
  }
  
  fit_main <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  p_main <- 1
  if (!is.null(fit_main)) p_main <- as.numeric(regTermTest(fit_main, ~LOCATION_LAT_c)$p)
  
  if (!is.na(p_main) && round(p_main, 3) <= 0.055) {
    arch_lats <- seq(min(df_island$LATITUDE, na.rm = TRUE), max(df_island$LATITUDE, na.rm = TRUE), length.out = 150)
    pred_arch <- predict(fit_main, newdata = data.frame(LOCATION_LAT_c = arch_lats - mean_lat_ref), se.fit = TRUE)
    return(data.frame(
      LATITUDE = arch_lats, fit = as.numeric(pred_arch), se = as.numeric(SE(pred_arch)), group = "Archipelago"
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))
  }
  return(NULL)
}

fit_df_full <- get_model_fit("VERTICAL_ACCRETION_MM_YR", des, df_strs_f4)

# ------------------------------------------------------------------------------
# 6. FIGURE COMPOSITION
# ------------------------------------------------------------------------------

fig4_final <- ggplot() +
  
  # A. Sea-Level Rise Shaded Reference Envelope
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = slr_box_min, ymax = slr_box_max,
    fill = slr_band_fill, alpha = 0.55
  ) +
  # Header Label: Top-left of lavender envelope
  annotate(
    "text", x = 30.1, y = slr_box_max - 0.5,
    label = "Range of sea level rise projections by 2050",
    color = "#52525B", size = 2.4, fontface = "italic", hjust = 0
  ) +
  
  # B. Accretion Baseline (0 mm/yr dashed line)
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#71717A",
    linewidth = 0.45
  ) +
  
  # C. IPCC AR6 Scenario Benchmarks (Dotted Lines in line with text)
  # Left segments (left axis to label start)
  annotate("segment", x = 30.2, xend = 29.85, y = rate_ssp585, yend = rate_ssp585, linetype = "dotted", color = col_ssp585, linewidth = 0.40) +
  annotate("segment", x = 30.2, xend = 29.85, y = rate_ssp370, yend = rate_ssp370, linetype = "dotted", color = col_ssp370, linewidth = 0.40) +
  annotate("segment", x = 30.2, xend = 29.85, y = rate_ssp245, yend = rate_ssp245, linetype = "dotted", color = col_ssp245, linewidth = 0.40) +
  annotate("segment", x = 30.2, xend = 29.85, y = rate_ssp126, yend = rate_ssp126, linetype = "dotted", color = col_ssp126, linewidth = 0.40) +
  
  # Right segments (cleanly matched to 28.95)
  annotate("segment", x = 28.95, xend = 18.6, y = rate_ssp585, yend = rate_ssp585, linetype = "dotted", color = col_ssp585, linewidth = 0.40) +
  annotate("segment", x = 28.95, xend = 18.6, y = rate_ssp370, yend = rate_ssp370, linetype = "dotted", color = col_ssp370, linewidth = 0.40) +
  annotate("segment", x = 28.95, xend = 18.6, y = rate_ssp245, yend = rate_ssp245, linetype = "dotted", color = col_ssp245, linewidth = 0.40) +
  annotate("segment", x = 28.95, xend = 18.6, y = rate_ssp126, yend = rate_ssp126, linetype = "dotted", color = col_ssp126, linewidth = 0.40) +
  
  # D. Scenario Text Labels (In line with the dotted lines, no box, no fill)
  annotate(
    "text", x = 29.8, y = rate_ssp585, label = "SSP5–8.5", color = col_ssp585,
    size = 2.3, fontface = "bold", hjust = 0, vjust = 0.5
  ) +
  annotate(
    "text", x = 29.8, y = rate_ssp370, label = "SSP3–7.0", color = col_ssp370,
    size = 2.3, fontface = "bold", hjust = 0, vjust = 0.5
  ) +
  annotate(
    "text", x = 29.8, y = rate_ssp245, label = "SSP2–4.5", color = col_ssp245,
    size = 2.3, fontface = "bold", hjust = 0, vjust = 0.5
  ) +
  annotate(
    "text", x = 29.8, y = rate_ssp126, label = "SSP1–2.6", color = col_ssp126,
    size = 2.3, fontface = "bold", hjust = 0, vjust = 0.5
  )

# E. Statistical Trendline & 95% Confidence Interval Ribbon
if (!is.null(fit_df_full)) {
  fig4_final <- fig4_final +
    geom_ribbon(
      data = fit_df_full,
      aes(x = LATITUDE, ymin = ymin, ymax = ymax, group = group),
      fill = col_trend_band, alpha = 0.16, inherit.aes = FALSE
    ) +
    geom_line(
      data = fit_df_full,
      aes(x = LATITUDE, y = fit, group = group),
      color = col_trend, linewidth = 0.70, inherit.aes = FALSE
    )
}

# F. Data Points: Site Jitter & Survey-Weighted Means
fig4_final <- fig4_final +
  # 1. Sites below SLR envelope (Context: semi-transparent, size = 1.25)
  geom_jitter(
    data = df_strs_site %>% filter(VERTICAL_ACCRETION_MM_YR < slr_box_min),
    aes(x = VISUAL_LAT, y = VERTICAL_ACCRETION_MM_YR, fill = LATITUDE),
    shape = 21, color = col_site_outline, alpha = 0.35, size = 1.5, stroke = 0.25,
    width = 0.04, height = 0, seed = 42
  ) +
  # 2. Sites inside lavender SLR envelope (High visibility: opaque, sharp border, same size = 1.25)
  geom_jitter(
    data = df_strs_site %>% filter(VERTICAL_ACCRETION_MM_YR >= slr_box_min),
    aes(x = VISUAL_LAT, y = VERTICAL_ACCRETION_MM_YR, fill = LATITUDE),
    shape = 21, color = "#18181B", alpha = 0.85, size = 1.5, stroke = 0.25,
    width = 0.04, height = 0, seed = 42
  ) +
  # 3. Survey-weighted island means (Foreground diamonds)
  geom_point(
    data = df_strs_f4,
    aes(x = VISUAL_LAT, y = VERTICAL_ACCRETION_MM_YR, fill = LATITUDE),
    shape = 23, size = 3.0, color = col_mean_outline, stroke = 0.45, alpha = 1.0
  ) +
  scale_fill_gradientn(
    colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30),
    oob = scales::squish, guide = "none"
  ) +
  
  # G. PNAS Coordinates & Typography
  scale_x_reverse(
    limits = c(30.2, 18.6),
    breaks = seq(19, 30, by = 1),
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(
    breaks = breaks_width(2)
  ) +
  coord_cartesian(clip = "off") +
  theme_pnas() +
  labs(
    x = "Latitude (°N)",
    y = expression(bold("Vertical accretion rate (mm" ~ yr^{-1} * ")"))
  )

# Display in RStudio
print(fig4_final)

# ------------------------------------------------------------------------------
# 7. EXPORT (PNAS 2-COLUMN STANDARD)
# ------------------------------------------------------------------------------

out_png  <- file.path(main_dir, "Figure4_VARvsSLR.png")
out_tiff <- file.path(main_dir, "Figure4_VARvsSLR.tiff")

ggsave(out_png,  fig4_final, width = 17.8, height = 9.8, units = "cm", dpi = 600, bg = "white")
ggsave(out_tiff, fig4_final, width = 17.8, height = 9.8, units = "cm", dpi = 600, compression = "lzw", bg = "white")

cat("\nFigure 4 exported successfully to PNG and TIFF.\n")