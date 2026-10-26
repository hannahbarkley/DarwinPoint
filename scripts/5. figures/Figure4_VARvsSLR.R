# ==============================================================================
# Script:  Figure4_VARvsSLR.R
# Purpose: Vertical Reef Accretion Rate (VAR) vs. 2050 Sea-Level Rise (SLR)
#          Formatted to Proceedings of the National Academy of Sciences (PNAS)
#          specifications (2-column width: 17.8 cm, Helvetica/Arial typography,
#          0.5 pt line hierarchy, and official IPCC AR6 scenario palettes).
#          Statistical trendlines follow the survey GLM model hierarchy:
#            1. Separate slopes for significant Lat x Region interactions.
#            2. Parallel slopes for significant additive Lat + Region effects.
#            3. Regional baselines if only Region is significant.
#            4. Continuous slope if only archipelago-wide Latitude is significant.
# Inputs:  - data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          - HA_2024_StRS_SiteLevel_Weighted.csv
#          - data/strs data/processed data/HA_2024_StRS_SurveyDesign.rds
#          - data/supporting data/IPCC_AR6_SLR_Rates_2040_2060.csv
# Outputs: - outputs/main figures/Figure4_VARvsSLR.tiff (600 DPI LZW raster)
#          - outputs/main figures/Figure4_VARvsSLR.png
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(survey)
  library(scales)
})

options(survey.lonely.psu = "adjust")

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

# --- 1. PNAS COLOR PALETTES & IPCC STANDARDS ---
pal_lat_grad <- c("#D73027", "#F46D43", "#FDAE61", "#FFFFBF", "#74ADD1", "#4575B4")

# IPCC AR6 Working Group I Official Scenario Colors
col_ssp126    <- "#003466"  # SSP1-2.6 (Deep Navy)
col_ssp245    <- "#F69320"  # SSP2-4.5 (Amber / Gold)
col_ssp370    <- "#DF0000"  # SSP3-7.0 (Red)
col_ssp585    <- "#980002"  # SSP5-8.5 (Dark Crimson)
slr_band_fill <- "#EBF4FA"  # Pale ocean-blue envelope

# ------------------------------------------------------------------------------
# 2. HELPERS: ARCHIPELAGO NOMENCLATURE & PNAS THEME
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

# ------------------------------------------------------------------------------
# 3. DATA INGESTION & GEOGRAPHIC DODGING
# ------------------------------------------------------------------------------
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

slr_box_min <- min(df_slr$Lower_17th_q17, na.rm = TRUE)  # 4.2 mm/yr
slr_box_max <- max(df_slr$Upper_83rd_q83, na.rm = TRUE)  # 11.8 mm/yr

rate_ssp126 <- df_slr %>% filter(grepl("126|1-2.6", Scenario)) %>% pull(Median_q50) %>% first()
rate_ssp245 <- df_slr %>% filter(grepl("245|2-4.5", Scenario)) %>% pull(Median_q50) %>% first()
rate_ssp370 <- df_slr %>% filter(grepl("370|3-7.0", Scenario)) %>% pull(Median_q50) %>% first()
rate_ssp585 <- df_slr %>% filter(grepl("585|5-8.5", Scenario)) %>% pull(Median_q50) %>% first()

# Island-level data
df_strs_darwin <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>%
  mutate(
    LOCATION                 = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                        if ("ISLAND" %in% names(.)) ISLAND else NULL),
    PLOT_LABEL               = format_islands(LOCATION),
    VERTICAL_ACCRETION_MM_YR = coalesce(if ("VERTICAL_ACCRETION_MM_YR" %in% names(.)) VERTICAL_ACCRETION_MM_YR else NULL, 0)
  )

# Exact customized repelling coordinate offsets
df_strs_f4 <- df_strs_darwin %>% mutate(
  nudge_x = case_when(
    PLOT_LABEL == "Hawai'i"            ~  0.6,  
    PLOT_LABEL == "Kaho'olawe"         ~  0.2, 
    PLOT_LABEL == "Maui"               ~ -0.1,  # Shift label slightly left with box
    PLOT_LABEL == "L\u0101na'i"        ~ -0.6, 
    PLOT_LABEL == "Moloka'i"           ~  0.1, 
    PLOT_LABEL == "O'ahu"              ~ -0.1,
    PLOT_LABEL == "Kaua'i"             ~ -0.6,  
    PLOT_LABEL == "Ni'ihau"            ~ -0.4, 
    PLOT_LABEL == "Lalo"               ~ -0.6,  
    PLOT_LABEL == "Kapou"              ~ -0.6,  
    PLOT_LABEL == "Manawai"            ~  0.5,  
    PLOT_LABEL == "H\u014dlanik\u016b" ~ -0.5, 
    TRUE ~ 0
  ),
  nudge_y = case_when(
    PLOT_LABEL == "Hawai'i"            ~ -0.5,  
    PLOT_LABEL == "Kaho'olawe"         ~ -2.9, 
    PLOT_LABEL == "Maui"               ~ -2.5,  
    PLOT_LABEL == "L\u0101na'i"        ~  2.0,  
    PLOT_LABEL == "Moloka'i"           ~ -3.4,  
    PLOT_LABEL == "O'ahu"              ~ -2.2,  
    PLOT_LABEL == "Kaua'i"             ~  1.2,  
    PLOT_LABEL == "Ni'ihau"            ~ -1.8,  
    PLOT_LABEL == "Lalo"               ~  0.6,  
    PLOT_LABEL == "Kapou"              ~  0.6,  
    PLOT_LABEL == "Manawai"            ~ -1.5,  
    PLOT_LABEL == "H\u014dlanik\u016b" ~ -1.5, 
    TRUE ~ 0
  )
)

# Site-level data with latitudinal alignment
df_strs_site <- saferead_csv("HA_2024_StRS_SiteLevel_Weighted.csv") %>%
  mutate(
    LOCATION                 = coalesce(if ("LOCATION" %in% names(.)) LOCATION else NULL,
                                        if ("ISLAND" %in% names(.)) ISLAND else NULL),
    PLOT_LABEL               = format_islands(LOCATION),
    VERTICAL_ACCRETION_MM_YR = coalesce(if ("VERTICAL_ACCRETION_MM_YR" %in% names(.)) VERTICAL_ACCRETION_MM_YR else NULL, 0)
  ) %>%
  select(-any_of("LATITUDE")) %>%
  left_join(df_strs_f4 %>% select(PLOT_LABEL, LATITUDE), by = "PLOT_LABEL") %>%
  filter(!is.na(LATITUDE))

# Geographic dodge to visually separate islands with overlapping latitudes
# Maui adjusted to LATITUDE - 0.05 to clear Kaho'olawe and stay spaced from Lana'i
df_strs_site <- df_strs_site %>%
  mutate(
    VISUAL_LAT = case_when(
      PLOT_LABEL == "L\u0101na'i" ~ LATITUDE + 0.14,
      PLOT_LABEL == "Maui"        ~ LATITUDE - 0.09,
      PLOT_LABEL == "Kaua'i"      ~ LATITUDE + 0.08,
      PLOT_LABEL == "Ni'ihau"     ~ LATITUDE - 0.08,
      TRUE                        ~ LATITUDE
    )
  )

df_strs_f4 <- df_strs_f4 %>%
  mutate(
    VISUAL_LAT = case_when(
      PLOT_LABEL == "L\u0101na'i" ~ LATITUDE + 0.14,
      PLOT_LABEL == "Maui"        ~ LATITUDE - 0.05,
      PLOT_LABEL == "Kaua'i"      ~ LATITUDE + 0.08,
      PLOT_LABEL == "Ni'ihau"     ~ LATITUDE - 0.08,
      TRUE                        ~ LATITUDE
    )
  )

des_file <- here("data", "strs data", "processed data", "HA_2024_StRS_SurveyDesign.rds")
des <- if (file.exists(des_file)) readRDS(des_file) else NULL

if (!is.null(des)) {
  mean_lat <- mean(des$variables$LOCATION_LAT, na.rm = TRUE)
  des <- update(des, LOCATION_LAT_c = LOCATION_LAT - mean_lat)
}

# ------------------------------------------------------------------------------
# 4. STATISTICAL MODEL FITS (SURVEY-WEIGHTED GLM HIERARCHY)
# ------------------------------------------------------------------------------
get_model_fit <- function(y_var, des_obj, df_island, mean_lat_ref = mean_lat) {
  if (is.null(des_obj) || !y_var %in% names(des_obj$variables)) return(NULL)
  
  mhi_min  <- min(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  mhi_max  <- max(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  nwhi_min <- min(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  nwhi_max <- max(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  
  mhi_lats   <- seq(mhi_min, mhi_max, length.out = 150)
  nwhi_lats  <- seq(nwhi_min, nwhi_max, length.out = 180)
  reg_levels <- if (is.factor(des_obj$variables$REGIONCODE)) levels(des_obj$variables$REGIONCODE) else unique(des_obj$variables$REGIONCODE)
  
  # Step 1: Interaction Model (y ~ LOCATION_LAT_c * REGIONCODE)
  fit_inter <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c * REGIONCODE")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  
  p_inter <- 1
  if (!is.null(fit_inter)) {
    wald_int <- regTermTest(fit_inter, ~LOCATION_LAT_c:REGIONCODE)
    p_inter  <- as.numeric(wald_int$p)
  }
  
  if (!is.na(p_inter) && p_inter < 0.05) {
    cat(sprintf("   ✓ %s: Significant REGIONAL INTERACTION (p = %.3f). Fitting separate slopes.\n", y_var, p_inter))
    pred_mhi  <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = mhi_lats - mean_lat_ref, REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = nwhi_lats - mean_lat_ref, REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    return(bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))
  }
  
  # Step 2: Additive Model (y ~ LOCATION_LAT_c + REGIONCODE)
  fit_add <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c + REGIONCODE")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  
  p_lat_add <- 1
  p_reg_add <- 1
  if (!is.null(fit_add)) {
    wald_lat  <- regTermTest(fit_add, ~LOCATION_LAT_c)
    wald_reg  <- regTermTest(fit_add, ~REGIONCODE)
    p_lat_add <- as.numeric(wald_lat$p)
    p_reg_add <- as.numeric(wald_reg$p)
  }
  
  # Condition 2A: Significant Parallel Regional Slopes
  if (!is.na(p_lat_add) && p_lat_add < 0.05 && !is.na(p_reg_add) && p_reg_add < 0.05) {
    cat(sprintf("   ✓ %s: Significant PARALLEL REGIONAL SLOPES (lat p = %.3e, reg p = %.3e).\n", y_var, p_lat_add, p_reg_add))
    pred_mhi  <- predict(fit_add, newdata = data.frame(LOCATION_LAT_c = mhi_lats - mean_lat_ref, REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_add, newdata = data.frame(LOCATION_LAT_c = nwhi_lats - mean_lat_ref, REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    return(bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))
  }
  
  # Condition 2B: Significant Region Step-Shift Only
  if (!is.na(p_reg_add) && p_reg_add < 0.05 && (is.na(p_lat_add) || p_lat_add >= 0.05)) {
    cat(sprintf("   ✓ %s: Significant REGIONAL STEP-SHIFT ONLY (reg p = %.3e). Fitting horizontal baselines.\n", y_var, p_reg_add))
    fit_reg   <- svyglm(as.formula(paste(y_var, "~ REGIONCODE")), design = des_obj, family = gaussian())
    pred_mhi  <- predict(fit_reg, newdata = data.frame(REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_reg, newdata = data.frame(REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    return(bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))
  }
  
  # Step 3: Single Continuous Model (y ~ LOCATION_LAT_c)
  fit_main <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  
  p_main <- 1
  if (!is.null(fit_main)) {
    wald_main <- regTermTest(fit_main, ~LOCATION_LAT_c)
    p_main    <- as.numeric(wald_main$p)
  }
  
  if (!is.na(p_main) && round(p_main, 3) <= 0.055) {
    cat(sprintf("   ✓ %s: Significant ARCHIPELAGO SLOPE (p = %.3f). Fitting continuous line.\n", y_var, p_main))
    arch_lats <- seq(min(df_island$LATITUDE, na.rm = TRUE), max(df_island$LATITUDE, na.rm = TRUE), length.out = 150)
    pred_arch <- predict(fit_main, newdata = data.frame(LOCATION_LAT_c = arch_lats - mean_lat_ref), se.fit = TRUE)
    
    return(data.frame(
      LATITUDE = arch_lats,
      fit      = as.numeric(pred_arch),
      se       = as.numeric(SE(pred_arch)),
      group    = "Archipelago"
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se))
  }
  
  cat(sprintf("   - %s: No significant latitudinal or regional trend. No fit added.\n", y_var))
  return(NULL)
}

cat("Evaluating survey GLM fit for Vertical Accretion...\n")
fit_df_full <- get_model_fit("VERTICAL_ACCRETION_MM_YR", des, df_strs_f4)

# ------------------------------------------------------------------------------
# 5. FIGURE COMPOSITION
# ------------------------------------------------------------------------------
fig4_final <- ggplot() +
  
  # A. Sea-Level Rise Shaded Reference Envelope
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = slr_box_min, ymax = slr_box_max, 
    fill = slr_band_fill, alpha = 0.95
  ) +
  
  # B. Net Neutral Accretion Baseline (0 mm/yr)
  geom_hline(yintercept = 0, linetype = "dashed", color = "#71717A", linewidth = 0.45) +
  
  # C. Edge-to-Edge IPCC AR6 Scenario Projection Benchmarks
  geom_hline(yintercept = rate_ssp585, linetype = "dotted", color = col_ssp585, linewidth = 0.55) +
  geom_hline(yintercept = rate_ssp370, linetype = "dotted", color = col_ssp370, linewidth = 0.55) +
  geom_hline(yintercept = rate_ssp245, linetype = "dotted", color = col_ssp245, linewidth = 0.55) +
  geom_hline(yintercept = rate_ssp126, linetype = "dotted", color = col_ssp126, linewidth = 0.55) +
  
  # D. Scenario Badges
  annotate("label", x = 29.8, y = rate_ssp585, 
           label = "SSP5\u20138.5", color = col_ssp585, 
           size = 2.3, fontface = "bold", hjust = 0, 
           fill = slr_band_fill, label.size = NA, label.padding = unit(0.12, "lines")) +
  
  annotate("label", x = 29.8, y = rate_ssp370, 
           label = "SSP3\u20137.0", color = col_ssp370, 
           size = 2.3, fontface = "bold", hjust = 0, 
           fill = slr_band_fill, label.size = NA, label.padding = unit(0.12, "lines")) +
  
  annotate("label", x = 29.8, y = rate_ssp245, 
           label = "SSP2\u20134.5", color = col_ssp245, 
           size = 2.3, fontface = "bold", hjust = 0, 
           fill = slr_band_fill, label.size = NA, label.padding = unit(0.12, "lines")) +
  
  annotate("label", x = 29.8, y = rate_ssp126, 
           label = "SSP1\u20132.6", color = col_ssp126, 
           size = 2.3, fontface = "bold", hjust = 0, 
           fill = slr_band_fill, label.size = NA, label.padding = unit(0.12, "lines"))

# E. Statistical Trend & 95% Confidence Interval Band
if (!is.null(fit_df_full)) {
  fig4_final <- fig4_final +
    geom_ribbon(
      data = fit_df_full,
      aes(x = LATITUDE, ymin = ymin, ymax = ymax, group = group),
      fill = "#64748B", alpha = 0.18, inherit.aes = FALSE
    ) +
    geom_line(
      data = fit_df_full,
      aes(x = LATITUDE, y = fit, group = group),
      color = "#1E293B", linewidth = 0.65, linetype = "solid", inherit.aes = FALSE
    )
}

# F. Site-Level Stratified Variance & Island Weighted Means
fig4_final <- fig4_final +
  geom_boxplot(
    data = df_strs_site,
    aes(x = VISUAL_LAT, y = VERTICAL_ACCRETION_MM_YR, group = PLOT_LABEL, fill = LATITUDE),
    color = "#3F3F46", alpha = 0.85, width = 0.18, linewidth = 0.38, position = "identity",
    fatten = 0, outlier.shape = 21, outlier.size = 1.0, outlier.color = "#3F3F46", outlier.alpha = 0.7
  ) +
  stat_summary(
    data = df_strs_site,
    aes(x = VISUAL_LAT, y = VERTICAL_ACCRETION_MM_YR, group = PLOT_LABEL, fill = LATITUDE),
    fun = mean, geom = "point", shape = 23, color = "#09090B", size = 2.6, stroke = 0.55,
    position = "identity"
  ) +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  geom_text_repel(
    data = df_strs_f4,
    aes(x = VISUAL_LAT, y = VERTICAL_ACCRETION_MM_YR, label = PLOT_LABEL), 
    nudge_x = df_strs_f4$nudge_x, 
    nudge_y = df_strs_f4$nudge_y, 
    size = 2.6, fontface = "bold", color = "#09090B", 
    box.padding = 0.25, point.padding = 0.35, force = 3, max.overlaps = Inf, seed = 42, 
    segment.color = "#71717A", segment.alpha = 0.75, segment.linewidth = 0.35, min.segment.length = 0
  ) +
  
  # G. PNAS Axis Scales & SI Typography
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
    y = expression(bold("Vertical accretion rate (mm" ~ yr^{ -1 } * ")"))
  )

# ------------------------------------------------------------------------------
# 6. EXPORT (PNAS 2-COLUMN STANDARD SPECIFICATIONS)
# ------------------------------------------------------------------------------
out_png  <- file.path(main_dir, "Figure4_VARvsSLR.png")
out_tiff <- file.path(main_dir, "Figure4_VARvsSLR.tiff")

ggsave(out_png,  fig4_final, width = 17.8, height = 9.8, units = "cm", dpi = 600, bg = "white")
ggsave(out_tiff, fig4_final, width = 17.8, height = 9.8, units = "cm", dpi = 600, compression = "lzw", bg = "white")

cat("\n✅ Figure 4 generated and exported to PNG and TIFF.\n")

# ------------------------------------------------------------------------------
# 7. SEA-LEVEL RISE ACCRETION DEFICIT CALCULATIONS
# ------------------------------------------------------------------------------
slr_min <- slr_box_min
pct_unweighted <- mean(df_strs_site$VERTICAL_ACCRETION_MM_YR > slr_min, na.rm = TRUE) * 100

cat(sprintf("\n🌊 2050 SEA-LEVEL RISE DEFICIT (Threshold = %.2f mm/yr):\n", slr_min))
cat(sprintf("  • Sites exceeding SLR (unweighted): %.2f%%\n", pct_unweighted))