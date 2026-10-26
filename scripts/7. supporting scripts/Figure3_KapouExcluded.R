# ==============================================================================
# Script:  FigureS7_KapouExcluded.R
# Purpose: Archipelago-Wide StRS Metric Scaling across Latitude with 
#          Survey-Weighted Statistically Significant Regression Fits 
#          EXCLUDING Kapou (Lisianski) for sensitivity analysis.
#          Aesthetics and layout harmonized to match Figure 3.
# Inputs:  - data/strs data/processed data/HA_2024_StRS_IslandLevel.csv
#          - data/strs data/processed data/HA_2024_StRS_SurveyDesign.rds
# Outputs: - outputs/supporting figures/Figure3_KapouExcluded.png
#          - outputs/supporting figures/Figure3_KapouExcluded.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggrepel)
  library(ggpubr)
  library(survey)
})

options(survey.lonely.psu = "adjust")

supp_dir <- here("outputs", "supporting figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

pal_lat_grad <- c("#D73027", "#F46D43", "#FDAE61", "#FFFFBF", "#74ADD1", "#4575B4")

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
# 1. DATA INGESTION, CENTERING, & SUBSETTING KAPOU
# ------------------------------------------------------------------------------
island_file <- here("data", "strs data", "processed data", "HA_2024_StRS_IslandLevel.csv")
des_file    <- here("data", "strs data", "processed data", "HA_2024_StRS_SurveyDesign.rds")

stopifnot("HA_2024_StRS_IslandLevel.csv not found!" = file.exists(island_file))

# Full dataset used for PLOTTING all island points
df_strs_darwin_full <- read_csv(island_file, show_col_types = FALSE) %>%
  mutate(
    PLOT_LABEL = format_islands(LOCATION),
    is_kapou   = ifelse(PLOT_LABEL == "Kapou", "Kapou", "Other")
  )

# Excluded dataset used strictly for MODELING and dynamic line bounds
df_strs_no_kapou <- df_strs_darwin_full %>%
  filter(!toupper(as.character(LOCATION)) %in% c("LIS", "LISIANSKI"))

island_locs_s6 <- df_strs_darwin_full %>%
  select(PLOT_LABEL, LATITUDE, is_kapou) %>%
  distinct(PLOT_LABEL, .keep_all = TRUE) %>%
  arrange(LATITUDE)

des_full <- if (file.exists(des_file)) readRDS(des_file) else NULL

if (!is.null(des_full)) {
  mean_lat <- mean(des_full$variables$LOCATION_LAT, na.rm = TRUE)
  des_full <- update(des_full, LOCATION_LAT_c = LOCATION_LAT - mean_lat)
  
  des_no_kapou <- subset(
    des_full,
    if ("ISLAND" %in% names(des_full$variables)) {
      !toupper(as.character(ISLAND)) %in% c("LIS", "LISIANSKI")
    } else {
      !toupper(as.character(LOCATION)) %in% c("LIS", "LISIANSKI")
    }
  )
}

# ------------------------------------------------------------------------------
# 2. FIT SURVEY GLM REGRESSIONS (Tight Data Bounding Without Kapou)
# ------------------------------------------------------------------------------
get_model_fit <- function(y_var, des_obj = des_no_kapou, df_island = df_strs_no_kapou, mean_lat_ref = mean_lat) {
  
  if (is.null(des_obj) || !y_var %in% names(des_obj$variables)) {
    return(list(fit_data = NULL))
  }
  
  mhi_min  <- min(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  mhi_max  <- max(df_island$LATITUDE[df_island$LATITUDE < 23], na.rm = TRUE)
  nwhi_min <- min(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  nwhi_max <- max(df_island$LATITUDE[df_island$LATITUDE > 23], na.rm = TRUE)
  arch_min <- min(df_island$LATITUDE, na.rm = TRUE)
  arch_max <- max(df_island$LATITUDE, na.rm = TRUE)
  
  # 1. Test Regional Interaction: LOCATION_LAT_c * REGIONCODE
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
    cat(sprintf("  ✓ %s: Significant REGIONAL INTERACTION (p = %.3f). Fitting separate MHI/NWHI lines.\n", y_var, p_inter))
    
    mhi_lats  <- seq(mhi_min, mhi_max, length.out = 60)
    nwhi_lats <- seq(nwhi_min, nwhi_max, length.out = 80)
    
    reg_levels <- levels(des_obj$variables$REGIONCODE)
    pred_mhi   <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = mhi_lats - mean_lat_ref, REGIONCODE = factor("MHI", levels = reg_levels)), se.fit = TRUE)
    pred_nwhi  <- predict(fit_inter, newdata = data.frame(LOCATION_LAT_c = nwhi_lats - mean_lat_ref, REGIONCODE = factor("NWHI", levels = reg_levels)), se.fit = TRUE)
    
    fit_df <- bind_rows(
      data.frame(LATITUDE = mhi_lats,  fit = as.numeric(pred_mhi),  se = as.numeric(SE(pred_mhi)),  group = "MHI"),
      data.frame(LATITUDE = nwhi_lats, fit = as.numeric(pred_nwhi), se = as.numeric(SE(pred_nwhi)), group = "NWHI")
    ) %>% mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se)
    
    return(list(fit_data = fit_df))
  }
  
  # 2. Test Archipelago-Wide Main Effect: LOCATION_LAT_c
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
    cat(sprintf("  ✓ %s: Significant ARCHIPELAGO SLOPE (p = %.3f). Fitting continuous line.\n", y_var, p_main))
    
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
  
  cat(sprintf("  - %s: No significant latitudinal trend (p = %.3f). No fit added.\n", y_var, p_main))
  return(list(fit_data = NULL))
}

cat("Evaluating Kapou-excluded survey GLM fits...\n")
fit_cov   <- get_model_fit("HARD_CORAL_COVER_PCT_MEAN")
fit_gprod <- get_model_fit("GROSS_CARB_PROD_KG_M2_YR_MEAN")
fit_geros <- get_model_fit("GROSS_EROSION_KG_M2_YR_MEAN")
fit_nprod <- get_model_fit("NET_CARB_PROD_KG_M2_YR_MEAN")

# ------------------------------------------------------------------------------
# 3. PLOT CONSTRUCTION
# ------------------------------------------------------------------------------
build_plot_s6 <- function(data, y_var, y_se, y_label, y_lim, y_brk, add_zero = FALSE, fit_res = NULL, show_x = TRUE) {
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
  
  # 1. 95% Confidence Interval Ribbons and Fitted Lines (Excluded Kapou Models)
  if (!is.null(fit_res$fit_data)) {
    p <- p +
      geom_ribbon(
        data = fit_res$fit_data,
        aes(x = LATITUDE, ymin = ymin, ymax = ymax, group = group),
        fill = "#94A3B8",
        alpha = 0.15,
        color = NA,
        inherit.aes = FALSE
      ) +
      geom_line(
        data = fit_res$fit_data,
        aes(x = LATITUDE, y = fit, group = group),
        color = "#334155",
        linewidth = 0.85,
        inherit.aes = FALSE
      )
  }
  
  # 2. Net Production Zero Baseline
  if (add_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6)
  }
  
  # 3. All Points Displayed (Kapou as Triangle shape 24, Others as Circle shape 21)
  p <- p + 
    geom_errorbar(aes(ymin = .data[[y_var]] - .data[[y_se]], 
                      ymax = .data[[y_var]] + .data[[y_se]]), 
                  color = "gray55", width = 0, linewidth = 0.5) +
    geom_point(aes(fill = LATITUDE, shape = is_kapou), size = 2.4, color = "gray20", stroke = 0.5, alpha = 0.85) +
    scale_shape_manual(values = c("Kapou" = 24, "Other" = 21), guide = "none") +
    scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
    labs(y = y_label)
  
  return(p)
}

lbl_cov   <- "Coral cover (%)"
lbl_gprod <- expression(bold("Gross production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
lbl_geros <- expression(bold("Gross erosion rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
lbl_nprod <- expression(bold("Net production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))

p_strs_cov   <- build_plot_s6(df_strs_darwin_full, "HARD_CORAL_COVER_PCT_MEAN", "se.HARD_CORAL_COVER_PCT_MEAN", lbl_cov, c(0, 30), seq(0, 30, 5), fit_res = fit_cov, show_x = FALSE)
p_strs_gprod <- build_plot_s6(df_strs_darwin_full, "GROSS_CARB_PROD_KG_M2_YR_MEAN", "se.GROSS_CARB_PROD_KG_M2_YR_MEAN", lbl_gprod, c(0, 6), seq(0, 6, 1), fit_res = fit_gprod, show_x = FALSE)
p_strs_geros <- build_plot_s6(df_strs_darwin_full, "GROSS_EROSION_KG_M2_YR_MEAN", "se.GROSS_EROSION_KG_M2_YR_MEAN", lbl_geros, c(0, 4), seq(0, 4, 1), fit_res = fit_geros, show_x = TRUE)
p_strs_nprod <- build_plot_s6(df_strs_darwin_full, "NET_CARB_PROD_KG_M2_YR_MEAN", "se.NET_CARB_PROD_KG_M2_YR_MEAN", lbl_nprod, c(-4.5, 4.2), seq(-4, 4, by = 2), add_zero = TRUE, fit_res = fit_nprod, show_x = TRUE)

main_plots_s6 <- ggarrange(
  p_strs_cov,   p_strs_gprod,
  p_strs_geros, p_strs_nprod,
  ncol = 2, nrow = 2, align = "hv", labels = "AUTO",
  font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02
)

# ------------------------------------------------------------------------------
# 4. TOP HORIZONTAL LEGEND (Figure 3 Matching Regional Brackets)
# ------------------------------------------------------------------------------
legend_nudges <- tribble(
  ~PLOT_LABEL,           ~leg_nx, ~leg_ny,
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

island_locs_s6 <- island_locs_s6 %>%
  left_join(legend_nudges, by = "PLOT_LABEL") %>%
  mutate(leg_nx = replace_na(leg_nx, 0), leg_ny = replace_na(leg_ny, 0.5))

p_legend_s6 <- ggplot() +
  geom_tile(
    data = data.frame(lat = seq(19, 30, length.out = 500), y = 0),
    aes(x = lat, y = y, fill = lat),
    height = 0.5
  ) +
  annotate("rect", xmin = 19, xmax = 30, ymin = -0.25, ymax = 0.25, color = "gray20", fill = NA, linewidth = 0.5) +
  
  # NWHI Bracket Bar
  annotate("segment", x = 30.0, xend = 22.8, y = 3.6, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 30.0, xend = 30.0, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 22.8, xend = 22.8, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("text", x = 26.4, y = 4.0, label = "Northwestern Hawaiian Islands", fontface = "bold", size = 3.8, color = "gray20") +
  
  # MHI Bracket Bar
  annotate("segment", x = 22.7, xend = 19.0, y = 3.6, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 22.7, xend = 22.7, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("segment", x = 19.0, xend = 19.0, y = 3.2, yend = 3.6, color = "gray30", linewidth = 0.6) +
  annotate("text", x = 20.8, y = 4.0, label = "Main Hawaiian Islands", fontface = "bold", size = 3.8, color = "gray20") +
  
  annotate("text", x = 30.2, y = 0, label = "Latitude (°N)", fontface = "bold", size = 3.5, hjust = 1, color = "black") +
  
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  scale_x_reverse(limits = c(30.2, 18.6), breaks = seq(19, 30, by = 1), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-0.6, 4.5), expand = c(0, 0)) + 
  coord_cartesian(clip = "off") +
  
  geom_point(data = island_locs_s6, aes(x = LATITUDE, y = 0, fill = LATITUDE, shape = is_kapou), size = 2.0, color = "gray20", stroke = 0.6) +
  scale_shape_manual(values = c("Kapou" = 24, "Other" = 21), guide = "none") +
  
  geom_text_repel(
    data = island_locs_s6, aes(x = LATITUDE, y = 0, label = PLOT_LABEL),
    nudge_x = island_locs_s6$leg_nx,
    nudge_y = island_locs_s6$leg_ny,
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
figS7_final <- ggarrange(p_legend_s6, main_plots_s6, nrow = 2, heights = c(0.31, 1))

ggsave(file.path(supp_dir, "Figure3_KapouExcluded.png"), figS7_final, width = 28, height = 17.5, units = "cm", dpi = 300, bg = "white")
