# ==============================================================================
# Script:  Figure3_FullStRS.R
# Purpose: Archipelago-Wide StRS Metric Scaling across Latitude with 
#          Survey-Weighted Statistically Significant Regression Fits.
#          Hierarchical model selection:
#            1. Separate regional slopes if Lat x Region interaction is significant (p < 0.05).
#            2. Parallel regional slopes if additive terms (Lat + Region) are significant (p < 0.05).
#            3. Regional horizontal baselines if Region is significant but Lat is not (p_reg < 0.05, p_lat >= 0.05).
#            4. Single continuous slope if only Latitude is significant across the archipelago.
#            5. No line if non-significant.
#          Error bars, ribbons, and fits are bounded flush to plot limits.
# Inputs:  - data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          - data/strs data/processed data/HA_2024_StRS_SurveyDesign.rds
# Outputs: - outputs/main figures/Figure3_FullStRS.png
#          - outputs/main figures/Figure3_FullStRS.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(ggpubr)
  library(survey)
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
# Helpers: Hawaiian Archipelago Names & Theme
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

theme_publication <- function(base_size = 9, y_margin_r = 4, x_margin_t = 5,
                              t_margin = 12, r_margin = 16, b_margin = 16, l_margin = 25) {
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
      plot.margin       = margin(t = t_margin, r = r_margin, b = b_margin, l = l_margin, unit = "pt")
    )
}

# ------------------------------------------------------------------------------
# 1. DATA INGESTION & DESIGN CENTERING
# ------------------------------------------------------------------------------
island_file <- here("data", "strs data", "processed data", "HA_2024_StRS_IslandLevel.csv")
des_file    <- here("data", "strs data", "processed data", "HA_2024_StRS_SurveyDesign.rds")

stopifnot("HA_2024_StRS_IslandLevel.csv not found!" = file.exists(island_file))

df_strs_darwin <- read_csv(island_file, show_col_types = FALSE) %>%
  mutate(PLOT_LABEL = format_islands(LOCATION))

island_locs_f3 <- df_strs_darwin %>%
  select(PLOT_LABEL, LATITUDE) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

des <- if (file.exists(des_file)) readRDS(des_file) else NULL

if (!is.null(des)) {
  mean_lat <- mean(des$variables$LOCATION_LAT, na.rm = TRUE)
  des <- update(des, LOCATION_LAT_c = LOCATION_LAT - mean_lat)
}

# ------------------------------------------------------------------------------
# 2. FIT SURVEY GLM REGRESSIONS (Hierarchical Model Selection)
# ------------------------------------------------------------------------------
get_model_fit <- function(y_var, des_obj = des, df_island = df_strs_darwin, mean_lat_ref = mean_lat) {
  
  if (is.null(des_obj) || !y_var %in% names(des_obj$variables)) {
    return(list(fit_data = NULL))
  }
  
  mhi_min  <- min(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  mhi_max  <- max(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  nwhi_min <- min(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  nwhi_max <- max(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  arch_min <- min(df_island$LATITUDE, na.rm = TRUE)
  arch_max <- max(df_island$LATITUDE, na.rm = TRUE)
  
  mhi_lats   <- seq(mhi_min, mhi_max, length.out = 150)
  nwhi_lats  <- seq(nwhi_min, nwhi_max, length.out = 180)
  reg_levels <- if (is.factor(des_obj$variables$REGIONCODE)) levels(des_obj$variables$REGIONCODE) else unique(des_obj$variables$REGIONCODE)
  
  # ----------------------------------------------------------------------------
  # Step 1: Full Interaction Model (y ~ LOCATION_LAT_c * REGIONCODE)
  # ----------------------------------------------------------------------------
  fit_inter <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c * REGIONCODE")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  
  p_inter <- 1
  if (!is.null(fit_inter)) {
    wald_int <- regTermTest(fit_inter, ~LOCATION_LAT_c:REGIONCODE)
    p_inter  <- as.numeric(wald_int$p)
  }
  
  # Condition 1: Significant Interaction -> Separate Slopes per Region
  if (!is.na(p_inter) && p_inter < 0.05) {
    cat(sprintf("   ✓ %s: Significant REGIONAL INTERACTION (p = %.3f). Fitting separate slopes.\n", y_var, p_inter))
    
    pred_mhi  <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = mhi_lats - mean_lat_ref, REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = nwhi_lats - mean_lat_ref, REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    fit_df <- bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se)
    
    return(list(fit_data = fit_df))
  }
  
  # ----------------------------------------------------------------------------
  # Step 2: Additive Model (y ~ LOCATION_LAT_c + REGIONCODE)
  # ----------------------------------------------------------------------------
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
  
  # Condition 2A: Significant Latitude AND Region -> Parallel Regional Slopes
  if (!is.na(p_lat_add) && p_lat_add < 0.05 && !is.na(p_reg_add) && p_reg_add < 0.05) {
    cat(sprintf("   ✓ %s: Significant PARALLEL REGIONAL SLOPES (lat p = %.3e, reg p = %.3e).\n", y_var, p_lat_add, p_reg_add))
    
    pred_mhi  <- predict(fit_add, newdata = data.frame(LOCATION_LAT_c = mhi_lats - mean_lat_ref, REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_add, newdata = data.frame(LOCATION_LAT_c = nwhi_lats - mean_lat_ref, REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    fit_df <- bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se)
    
    return(list(fit_data = fit_df))
  }
  
  # Condition 2B: Significant Region ONLY (Latitude non-significant) -> Regional Baseline Step-Shift
  if (!is.na(p_reg_add) && p_reg_add < 0.05 && (is.na(p_lat_add) || p_lat_add >= 0.05)) {
    cat(sprintf("   ✓ %s: Significant REGIONAL STEP-SHIFT ONLY (reg p = %.3e, lat p = %.3f). Fitting horizontal regional baselines.\n", y_var, p_reg_add, p_lat_add))
    
    fit_reg <- svyglm(as.formula(paste(y_var, "~ REGIONCODE")), design = des_obj, family = gaussian())
    pred_mhi  <- predict(fit_reg, newdata = data.frame(REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi <- predict(fit_reg, newdata = data.frame(REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    fit_df <- bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se)
    
    return(list(fit_data = fit_df))
  }
  
  # ----------------------------------------------------------------------------
  # Step 3: Single Continuous Model (y ~ LOCATION_LAT_c, evaluated only if region is not significant)
  # ----------------------------------------------------------------------------
  fit_main <- tryCatch(
    svyglm(as.formula(paste(y_var, "~ LOCATION_LAT_c")), design = des_obj, family = gaussian()),
    error = function(e) NULL
  )
  
  p_main <- 1
  if (!is.null(fit_main)) {
    wald_main <- regTermTest(fit_main, ~LOCATION_LAT_c)
    p_main    <- as.numeric(wald_main$p)
  }
  
  # Condition 3: Significant Archipelago-Wide Latitude -> Continuous Line
  if (!is.na(p_main) && round(p_main, 3) <= 0.055) {
    cat(sprintf("   ✓ %s: Significant ARCHIPELAGO SLOPE (p = %.3f). Fitting continuous line.\n", y_var, p_main))
    
    arch_lats <- seq(arch_min, arch_max, length.out = 150)
    pred_arch <- predict(fit_main, newdata = data.frame(LOCATION_LAT_c = arch_lats - mean_lat_ref), se.fit = TRUE)
    
    fit_df <- data.frame(
      LATITUDE = arch_lats,
      fit      = as.numeric(pred_arch),
      se       = as.numeric(SE(pred_arch)),
      group    = "Archipelago"
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se)
    
    return(list(fit_data = fit_df))
  }
  
  # Condition 4: Non-Significant -> No Line
  cat(sprintf("   - %s: No significant latitudinal or regional trend. No fit added.\n", y_var))
  return(list(fit_data = NULL))
}

cat("Evaluating statistically justified survey GLM fits...\n")
fit_cov   <- get_model_fit("HARD_CORAL_COVER_PCT_MEAN")
fit_gprod <- get_model_fit("GROSS_CARB_PROD_KG_M2_YR_MEAN")
fit_geros <- get_model_fit("GROSS_EROSION_KG_M2_YR_MEAN")
fit_nprod <- get_model_fit("NET_CARB_PROD_KG_M2_YR_MEAN")

# ------------------------------------------------------------------------------
# 3. PLOT CONSTRUCTION
# ------------------------------------------------------------------------------
build_plot_f3 <- function(data, y_var, y_se, y_label, y_lim, y_brk, add_zero = FALSE, fit_res = NULL, show_x = TRUE) {
  
  p <- ggplot(data, aes(x = LATITUDE, y = .data[[y_var]])) +
    scale_x_reverse(limits = c(30.2, 18.6), breaks = seq(19, 30, by = 1), expand = c(0.01, 0.01)) +
    scale_y_continuous(breaks = y_brk) +
    coord_cartesian(ylim = y_lim, clip = "on") +
    theme_publication()
  
  if (!show_x) {
    p <- p + labs(x = NULL) + theme(axis.title.x = element_blank())
  } else {
    p <- p + labs(x = "Latitude (°N)")
  }
  
  if (!is.null(fit_res$fit_data)) {
    # 1. Shading dataframe: bounded cleanly flush to y_lim (e.g. 0)
    fit_ribbon <- fit_res$fit_data %>%
      mutate(
        ymin = pmax(y_lim[1], ymin),
        ymax = pmax(ymin, pmin(y_lim[2], ymax))
      )
    
    # 2. Central fit line: filter out points below y_lim[1] so line terminates 
    #    naturally without bending horizontally along the axis
    fit_line <- fit_res$fit_data %>%
      filter(fit >= y_lim[1] - 0.02, fit <= y_lim[2] + 0.02)
    
    p <- p +
      geom_ribbon(
        data = fit_ribbon,
        aes(x = LATITUDE, ymin = ymin, ymax = ymax, group = group),
        fill = "#94A3B8",
        alpha = 0.15,
        color = NA,
        inherit.aes = FALSE
      ) +
      geom_line(
        data = fit_line,
        aes(x = LATITUDE, y = fit, group = group),
        color = "#334155",
        linewidth = 0.85,
        inherit.aes = FALSE
      )
  }
  
  if (add_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6)
  }
  
  # Error bars clamped within axis boundaries to prevent margin overshoots
  p <- p + 
    geom_errorbar(
      aes(
        ymin = pmax(y_lim[1], .data[[y_var]] - .data[[y_se]]), 
        ymax = pmin(y_lim[2], .data[[y_var]] + .data[[y_se]])
      ), 
      color = "gray55", width = 0, linewidth = 0.5
    ) +
    geom_point(aes(fill = LATITUDE), shape = 21, size = 2.4, color = "gray20", stroke = 0.5, alpha = 0.85) +
    scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
    labs(y = y_label)
  
  return(p)
}

lbl_cov   <- "Coral cover (%)"
lbl_gprod <- expression(bold("Gross production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
lbl_geros <- expression(bold("Gross erosion rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
lbl_nprod <- expression(bold("Net production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))

p_strs_cov   <- build_plot_f3(df_strs_darwin, "HARD_CORAL_COVER_PCT_MEAN", "se.HARD_CORAL_COVER_PCT_MEAN", lbl_cov, c(0, 30), seq(0, 30, 5), fit_res = fit_cov, show_x = FALSE)
p_strs_gprod <- build_plot_f3(df_strs_darwin, "GROSS_CARB_PROD_KG_M2_YR_MEAN", "se.GROSS_CARB_PROD_KG_M2_YR_MEAN", lbl_gprod, c(0, 6), seq(0, 6, 1), fit_res = fit_gprod, show_x = FALSE)
p_strs_geros <- build_plot_f3(df_strs_darwin, "GROSS_EROSION_KG_M2_YR_MEAN", "se.GROSS_EROSION_KG_M2_YR_MEAN", lbl_geros, c(0, 4.5), seq(0, 4, 1), fit_res = fit_geros, show_x = TRUE)
p_strs_nprod <- build_plot_f3(df_strs_darwin, "NET_CARB_PROD_KG_M2_YR_MEAN", "se.NET_CARB_PROD_KG_M2_YR_MEAN", lbl_nprod, c(-4.5, 4.2), seq(-4, 4, by = 2), add_zero = TRUE, fit_res = fit_nprod, show_x = TRUE)

main_plots_f3 <- ggarrange(
  p_strs_cov,   p_strs_gprod,
  p_strs_geros, p_strs_nprod,
  ncol = 2, nrow = 2, align = "hv", labels = "AUTO",
  font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02
)

# ------------------------------------------------------------------------------
# 4. TOP HORIZONTAL LEGEND (Balanced Spacing)
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
  "Hawai'i",              0.0,    0.6
)

island_locs_f3 <- island_locs_f3 %>%
  left_join(legend_nudges, by = "PLOT_LABEL") %>%
  mutate(leg_nx = replace_na(leg_nx, 0), leg_ny = replace_na(leg_ny, 0.5))

p_legend_f3 <- ggplot() +
  geom_tile(
    data = data.frame(lat = seq(19, 30, length.out = 500), y = 0),
    aes(x = lat, y = y, fill = lat),
    height = 0.5
  ) +
  annotate("rect", xmin = 19, xmax = 30, ymin = -0.25, ymax = 0.25, color = "gray20", fill = NA, linewidth = 0.5) +
  
  # NWHI Bracket Bar 
  annotate("segment", x = 30.0, xend = 22.8, y = 3.6, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 30.0, xend = 30, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 22.8, xend = 22.8, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("text", x = 26.4, y = 4.0, label = "Northwestern Hawaiian Islands", fontface = "bold", size = 3.8, color = "gray20") +
  
  # MHI Bracket Bar
  annotate("segment", x = 22.7, xend = 19.0, y = 3.6, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 22.7, xend = 22.7, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 19.0, xend = 19.0, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("text", x = 20.8, y = 4.0, label = "Main Hawaiian Islands", fontface = "bold", size = 3.8, color = "gray20") +
  
  annotate("text", x = 30.2, y = 0, label = "Latitude (°N)", fontface = "bold", size = 3.5, hjust = 1, color = "black") +
  
  scale_fill_gradientn(colors = pal_lat_grad, values = pal_lat_vals, limits = c(19, 30), guide = "none") +
  scale_x_reverse(limits = c(30.2, 18.6), breaks = seq(19, 30, by = 1), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-0.6, 4.5), expand = c(0, 0)) + 
  coord_cartesian(clip = "off") +
  
  geom_point(data = island_locs_f3, aes(x = LATITUDE, y = 0, fill = LATITUDE), shape = 21, size = 2.0, color = "gray20", stroke = 0.6) +
  
  geom_text_repel(
    data = island_locs_f3, aes(x = LATITUDE, y = 0, label = PLOT_LABEL),
    nudge_x = island_locs_f3$leg_nx,
    nudge_y = island_locs_f3$leg_ny,
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
    plot.margin       = margin(t = 8, r = 5, b = 18, l = 75, unit = "pt")
  )

# ------------------------------------------------------------------------------
# 5. FIGURE ASSEMBLY & EXPORT
# ------------------------------------------------------------------------------
fig3_final <- ggarrange(p_legend_f3, main_plots_f3, nrow = 2, heights = c(0.31, 1))

ggsave(file.path(main_dir, "Figure3_FullStRS.png"), fig3_final, width = 28, height = 17.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(main_dir, "Figure3_FullStRS.tiff"), fig3_final, width = 28, height = 17.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")