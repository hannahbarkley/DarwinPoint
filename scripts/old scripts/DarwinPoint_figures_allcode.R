# ==============================================================================
# 2024 Hawaiian Darwin Point & Carbonate Budget Figures (PNAS Edition)
# Overview: Combines Hawaiian StRS survey design weighting, survey GLMs, 
#           and multi-panel figure production for Figure 2, Figure 4, and Figure 5.
# Outputs: 
#   1) outputs/Fig2_Grigg_vs_NCRMP_Horizontal.png
#   2) outputs/Fig2_Grigg_vs_NCRMP_Vertical.png
#   3) outputs/Fig4_HI_2024_StRS_4Panel_Vertical.png
#   4) outputs/Fig5_NCRMP_Accretion_SLR.png
#   5) outputs/2024_HI_StRS_Islandlevel_DarwinPoint.csv
# ==============================================================================

library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(survey)


# --- ACCRETION & SLR CONSTANTS ---
rho_aragonite <- 2.89  # Density of solid aragonite (g/cm^3)
reef_porosity <- 0.50  # Estimated framework porosity (50%)
hi_slr_min    <- 6     # Predicted min SLR rate by 2050 (mm/yr)
hi_slr_max    <- 10    # Predicted max SLR rate by 2050 (mm/yr)


# 1. LOAD DATA -------------------------------------------------------------
strs    <- read_csv(here("data", "HA_2024_StRS_net.csv"), show_col_types = FALSE)
sm      <- read_csv(here("data", "SURVEY MASTER.csv"), show_col_types = FALSE)
sectors <- read_csv(here("data", "Sectors-Strata-Areas.csv"), show_col_types = FALSE)
strs_sw <- read_csv(here("outputs", "HA_2024_StRS_net_sw.csv"), show_col_types = FALSE)
grigg   <- read_csv(here("data", "Grigg_1982_data.csv"), show_col_types = FALSE)


# 2. DATA CLEANING, MERGING & SURVEY WEIGHTS -------------------------------

format_islands <- function(x) {
  case_match(
    x,
    "Oahu"                  ~ "O'ahu", 
    "Hawaii"                ~ "Hawai'i", 
    "Lanai"                 ~ "L\u0101na'i",
    "Molokai"               ~ "Moloka'i", 
    "Kauai"                 ~ "Kaua'i", 
    "Niihau"                ~ "Ni'ihau",
    "Kahoolawe"             ~ "Kaho'olawe", 
    "FFS"                   ~ "Lalo",
    "French Frigate Shoals" ~ "Lalo",
    "French Frigate"        ~ "Lalo",
    "Lisianski"             ~ "Kapou",
    "Pearl & Hermes"        ~ "Manawai",
    "Pearl and Hermes"      ~ "Manawai",
    "Kure"                  ~ "H\u014dlanik\u016b",
    # Added mappings for historical Grigg locations
    "Midway"                ~ "Kuaihelani",
    "Laysan"                ~ "Kamole",
    "Maro Reef"             ~ "Kamokuokamohoali'i",
    "Maro"                  ~ "Kamokuokamohoali'i",
    "Gardner Pinnacles"     ~ "P\u016bh\u0101honu",
    "Gardner"               ~ "P\u016bh\u0101honu",
    "Necker"                ~ "Mokumanamana",
    .default = x
  )
}

# Clean & Merge StRS Site Data
strs_clean <- strs %>%
  rename(OLD_SITE = SITE) %>%
  mutate(SITE = sub("-", "-0", OLD_SITE)) %>%
  select(SEC_NAME, SITE, OBS_YEAR, REGION, ISLAND, REEF_ZONE, DEPTH_BIN, 
         GROSS_PROD_KG_M2_YR, NET_PROD_KG_M2_YR, FISH_EROSION_KG_M2_YR_ALL, 
         URCHIN_EROSION_KG_M2_YR_MEAN, CORAL_COVER_PERC) %>%
  left_join(sm %>% select(SITE, LATITUDE_SV), by = "SITE")

# Unique site metadata
survey_site <- distinct(strs_clean, OBS_YEAR, REGION, ISLAND, SEC_NAME, SITE, REEF_ZONE, DEPTH_BIN, LATITUDE_SV)
if (any(table(survey_site$SITE) > 1)) stop("Duplicate sites found in metadata!")

meta <- left_join(survey_site, sectors, by = join_by(SEC_NAME, REEF_ZONE, DEPTH_BIN))

new_strs <- left_join(strs_clean, meta) %>%
  mutate(STRATANAME  = paste(SEC_NAME, REEF_ZONE, DEPTH_BIN, sep = "_"),
         REGION_YEAR = paste(REGION, OBS_YEAR, sep = "_"),
         DB_RZ       = paste(REEF_ZONE, DEPTH_BIN, sep = "_"))

# Calculate Sampling Weights & Island Latitudes
site_sw <- new_strs %>%
  mutate(n = n_distinct(SITE), sw = NH / n, .by = c(OBS_YEAR, SEC_NAME, REEF_ZONE, DEPTH_BIN, NH)) %>%
  mutate(ISLAND_LAT = mean(LATITUDE_SV, na.rm = TRUE), .by = ISLAND) %>%
  filter(n > 1) %>%
  mutate(across(c(OBS_YEAR, ISLAND), as.factor),
         Strat_conc = paste(OBS_YEAR, REGION, ISLAND, SEC_NAME, DB_RZ, sep = "_"),
         # Compute Gross Erosion (positive)
         URCHIN_EROS = if_else(is.na(URCHIN_EROSION_KG_M2_YR_MEAN), 0, URCHIN_EROSION_KG_M2_YR_MEAN),
         FISH_EROS   = if_else(is.na(FISH_EROSION_KG_M2_YR_ALL), 0, FISH_EROSION_KG_M2_YR_ALL),
         GROSS_EROSION_KG_M2_YR = URCHIN_EROS + FISH_EROS,
         # Compute Vertical Accretion
         VERTICAL_ACCRETION_MM_YR = NET_PROD_KG_M2_YR / (rho_aragonite * (1 - reef_porosity)))

# Construct Survey Design Object
des <- svydesign(id = ~1, strata = ~Strat_conc, weights = ~sw, data = site_sw)


# 3. SURVEY-WEIGHTED MEANS & MODELS ---------------------------------------

metrics <- c("GROSS_PROD_KG_M2_YR", "NET_PROD_KG_M2_YR", "VERTICAL_ACCRETION_MM_YR", "GROSS_EROSION_KG_M2_YR", "CORAL_COVER_PERC")

# Weighted Island Means
spatial_Imean <- svyby(reformulate(metrics), ~ISLAND, des, svymean, na.rm = TRUE) %>%
  left_join(distinct(site_sw, ISLAND, ISLAND_LAT), by = "ISLAND") %>% 
  rename(LATITUDE = ISLAND_LAT) %>%
  mutate(PLOT_LABEL = format_islands(as.character(ISLAND)))

if (!dir.exists(here("outputs"))) dir.create(here("outputs"))
write_csv(spatial_Imean, here("outputs", "2024_HI_StRS_Islandlevel_DarwinPoint.csv"))

# Fit Survey GLMs
mod_cov   <- svyglm(CORAL_COVER_PERC ~ ISLAND_LAT, design = des, family = gaussian())
mod_gprod <- svyglm(GROSS_PROD_KG_M2_YR ~ ISLAND_LAT, design = des, family = gaussian())
mod_geros <- svyglm(GROSS_EROSION_KG_M2_YR ~ ISLAND_LAT, design = des, family = gaussian())
mod_nprod <- svyglm(NET_PROD_KG_M2_YR ~ ISLAND_LAT, design = des, family = gaussian())
mod_vacc  <- svyglm(VERTICAL_ACCRETION_MM_YR ~ ISLAND_LAT, design = des, family = gaussian())


# 4. PREPARE FIGURE 2 DATA (GRIGG VS NCRMP) --------------------------------

# Added vertical accretion mapping for 1982 Grigg Data
grigg <- grigg %>% mutate(
  PLOT_LABEL = format_islands(LOCATION),
  VERTICAL_ACCRETION_MM_YR = GROSS_PROD / (rho_aragonite * (1 - reef_porosity))
)

strs_sw_means <- strs_sw %>%
  mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
  group_by(PLOT_LABEL) %>%
  summarize(
    LATITUDE                      = mean(LATITUDE, na.rm = TRUE),
    CORAL_COVER_SE                = sd(CORAL_COVER_PERC, na.rm = TRUE) / sqrt(sum(!is.na(CORAL_COVER_PERC))),
    CORAL_COVER_PERC              = mean(CORAL_COVER_PERC, na.rm = TRUE),
    GROSS_CARB_PROD_KG_M2_YR_SE   = sd(GROSS_PROD_KG_M2_YR, na.rm = TRUE) / sqrt(sum(!is.na(GROSS_PROD_KG_M2_YR))),
    GROSS_CARB_PROD_KG_M2_YR_MEAN = mean(GROSS_PROD_KG_M2_YR, na.rm = TRUE),
    .groups                       = "drop"
  )

shared_deltas <- inner_join(
  grigg %>% select(PLOT_LABEL, GRIGG_LAT = LATITUDE, GRIGG_COV = CORAL_COVER_PERC, GRIGG_PROD = GROSS_PROD),
  strs_sw_means %>% select(PLOT_LABEL, NCRMP_COV = CORAL_COVER_PERC, NCRMP_PROD = GROSS_CARB_PROD_KG_M2_YR_MEAN),
  by = "PLOT_LABEL"
) %>%
  mutate(
    LATITUDE   = GRIGG_LAT,
    DELTA_COV  = NCRMP_COV - GRIGG_COV,
    DELTA_PROD = NCRMP_PROD - GRIGG_PROD
  )


# 5. PLOTTING UTILITIES & THEME --------------------------------------------

lbl_cov        <- "Coral cover (%)"
lbl_prod       <- expression(bold("Gross production rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_delta_cov  <- "\u0394 Coral cover (%)"
lbl_delta_prod <- expression(bold(Delta~"Gross production rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_geros      <- expression(bold("Gross erosion rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_nprod      <- expression(bold("Net production rate (kg"~m^{-2}~yr^{-1}*")"))
lbl_vacc       <- expression(bold("Vertical accretion rate (mm"~yr^{-1}*")"))

theme_transparent <- function(txt_col = "black") {
  theme_bw() + theme(
    panel.grid        = element_blank(),
    panel.background  = element_rect(fill = "transparent", color = NA),
    plot.background   = element_rect(fill = "transparent", color = NA),
    plot.title        = element_text(color = txt_col, size = 10, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    axis.text         = element_text(color = txt_col, size = 8, face = "bold"),
    axis.title        = element_text(color = txt_col, size = 8.5, face = "bold"),
    panel.border      = element_rect(color = txt_col, linewidth = 0.6, fill = NA),
    axis.ticks        = element_line(color = txt_col, linewidth = 0.4),
    axis.title.x      = element_text(margin = margin(t = 6)),
    axis.title.y      = element_text(margin = margin(r = 8)),
    axis.ticks.length = unit(0.12, "cm"),
    plot.margin       = margin(t = 16, r = 8, b = 4, l = 12, unit = "pt")
  )
}

no_x_title <- function(p) {
  p + labs(x = NULL) + theme(axis.title.x = element_blank())
}

# Helper to generate prediction grids and p-value annotations from survey GLMs
get_glm_trend <- function(model, lat_range = seq(19, 30, length.out = 100)) {
  pred_grid <- data.frame(ISLAND_LAT = lat_range)
  preds     <- predict(model, newdata = pred_grid, se.fit = TRUE)
  
  pred_df <- data.frame(
    LATITUDE = lat_range,
    fit      = as.numeric(preds),
    se       = sqrt(as.numeric(attr(preds, "var")))
  ) %>%
    mutate(ymin = fit - 1.96 * se, ymax = fit + 1.96 * se)
  
  coefs   <- summary(model)$coefficients
  beta    <- coefs["ISLAND_LAT", "Estimate"]
  se_beta <- coefs["ISLAND_LAT", "Std. Error"]
  t_stat  <- coefs["ISLAND_LAT", "t value"]
  p_val   <- coefs["ISLAND_LAT", "Pr(>|t|)"]
  df_res  <- model$df.residual
  
  p_str <- ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val))
  
  stat_text <- sprintf("\u03B2 = %.2f \u00B1 %.2f\nt(%d) = %.2f, %s", 
                       beta, se_beta, df_res, t_stat, p_str)
  
  list(pred_df = pred_df, p_val = p_val, stat_text = stat_text)
}

# General Plot Builder
build_plot <- function(data, y_var, y_se = NULL, txt_col = "black", y_label, 
                       title_text = NULL, y_lim, y_brk, add_zero_line = FALSE, model = NULL) {
  seg_col   <- ifelse(txt_col == "#EEEBE6", "#EEEBE6", "gray50")
  err_alpha <- ifelse(txt_col == "#EEEBE6", 0.6, 1)
  
  p <- ggplot(data, aes(x = LATITUDE, y = .data[[y_var]])) +
    scale_x_reverse(limits = c(30, 19), breaks = seq(19, 30, by = 1)) +
    theme_transparent(txt_col)
  
  # Render Survey GLM Trendline & Confidence Interval if provided
  if (!is.null(model)) {
    m_info    <- get_glm_trend(model)
    line_type <- ifelse(m_info$p_val < 0.05, "solid", "dotted")
    
    p <- p + 
      geom_ribbon(data = m_info$pred_df, aes(x = LATITUDE, ymin = ymin, ymax = ymax),
                  fill = "gray85", alpha = 0.5, inherit.aes = FALSE) +
      geom_line(data = m_info$pred_df, aes(x = LATITUDE, y = fit),
                color = "gray30", linewidth = 0.6, linetype = line_type, inherit.aes = FALSE)
    
    y_stat_pos <- y_lim[2] - (y_lim[2] - y_lim[1]) * 0.05
    
    p <- p + annotate("text", x = 29.5, y = y_stat_pos, label = m_info$stat_text,
                      hjust = 0, vjust = 1, size = 2.5, fontface = "bold", 
                      color = "black", lineheight = 0.95)
  }
  
  if (!is.null(y_se)) {
    p <- p + geom_errorbar(
      aes(ymin = .data[[y_var]] - .data[[y_se]], ymax = .data[[y_var]] + .data[[y_se]]), 
      color = txt_col, alpha = err_alpha, width = 0, linewidth = 0.4
    )
  }
  
  p <- p + scale_y_continuous(limits = y_lim, breaks = y_brk, expand = c(0, 0))
  
  if (add_zero_line || (y_lim[1] < 0 && !is.null(title_text) && title_text == "2024 \u2212 1982")) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5)
  }
  
  p <- p + geom_point(color = txt_col, size = 1.8)
  
  p <- p + geom_text_repel(
    aes(label = PLOT_LABEL), 
    size              = 2.4, 
    fontface          = "bold", 
    color             = "black", 
    box.padding       = 0.35,   
    point.padding     = 0.25,   
    force             = 3,      
    max.overlaps      = Inf, 
    seed              = 123,
    segment.color     = seg_col, 
    segment.alpha     = 0.6, 
    segment.linewidth = 0.3, 
    min.segment.length= 0       
  ) +
    labs(title = title_text, x = "Latitude (°N)", y = y_label)
  
  return(p)
}


# 6. BUILD & EXPORT FIGURE 2 (GRIGG VS NCRMP) -----------------------------

grigg_cover <- build_plot(grigg, "CORAL_COVER_PERC", y_label = lbl_cov, title_text = "1982", y_lim = c(0, 100), y_brk = seq(0, 100, 20))
grigg_prod  <- build_plot(grigg, "GROSS_PROD", y_label = lbl_prod, title_text = "1982", y_lim = c(-2, 20), y_brk = seq(-2, 20, 4))

ncrmp_cover_sw <- build_plot(strs_sw_means, "CORAL_COVER_PERC", "CORAL_COVER_SE", y_label = lbl_cov, title_text = "2024", y_lim = c(0, 40), y_brk = seq(0, 40, 10))
ncrmp_prod_sw  <- build_plot(strs_sw_means, "GROSS_CARB_PROD_KG_M2_YR_MEAN", "GROSS_CARB_PROD_KG_M2_YR_SE", y_label = lbl_prod, title_text = "2024", y_lim = c(-2, 8), y_brk = seq(-2, 8, 2))

delta_cover <- build_plot(shared_deltas, "DELTA_COV", y_label = lbl_delta_cov, title_text = "2024 \u2212 1982", y_lim = c(-80, 20), y_brk = seq(-80, 20, 20))
delta_prod  <- build_plot(shared_deltas, "DELTA_PROD", y_label = lbl_delta_prod, title_text = "2024 \u2212 1982", y_lim = c(-15, 5), y_brk = seq(-15, 5, 5))

# Composites
grigg_ncrmp_plot_h <- ggarrange(
  no_x_title(grigg_cover), no_x_title(ncrmp_cover_sw), no_x_title(delta_cover),
  grigg_prod,              ncrmp_prod_sw,              delta_prod,
  ncol = 3, nrow = 2, align = "hv", labels = "auto",
  font.label = list(size = 10, face = "bold"), label.x = 0.01, label.y = 0.99
)

grigg_ncrmp_plot_v <- ggarrange(
  no_x_title(grigg_cover),    no_x_title(grigg_prod), 
  no_x_title(ncrmp_cover_sw), no_x_title(ncrmp_prod_sw),
  delta_cover,                delta_prod,
  ncol = 2, nrow = 3, align = "hv", labels = "auto",
  font.label = list(size = 10, face = "bold"), label.x = 0.01, label.y = 0.99
)

ggsave(here("outputs", "Fig2_Grigg_vs_NCRMP_Horizontal.png"), grigg_ncrmp_plot_h, width = 25, height = 12.5, units = "cm", dpi = 300, bg = "white")
ggsave(here("outputs", "Fig2_Grigg_vs_NCRMP_Vertical.png"), grigg_ncrmp_plot_v, width = 20, height = 21.5, units = "cm", dpi = 300, bg = "white")


# 7. BUILD & EXPORT FIGURE 4 (WITH SURVEY GLM RESULTS) -------------------

# 1) Coral Cover (mod_cov)
p_strs_cov <- build_plot(
  spatial_Imean, "CORAL_COVER_PERC", "se.CORAL_COVER_PERC", 
  y_label = lbl_cov, y_lim = c(0, 40), y_brk = seq(0, 40, 10), 
  model = mod_cov
)

# 2) Gross Production (mod_gprod)
p_strs_gprod <- build_plot(
  spatial_Imean, "GROSS_PROD_KG_M2_YR", "se.GROSS_PROD_KG_M2_YR", 
  y_label = lbl_prod, y_lim = c(-1, 8), y_brk = seq(-1, 8, 2), 
  add_zero_line = TRUE, model = mod_gprod
)

# 3) Gross Erosion (mod_geros)
p_strs_geros <- build_plot(
  spatial_Imean, "GROSS_EROSION_KG_M2_YR", "se.GROSS_EROSION_KG_M2_YR", 
  y_label = lbl_geros, y_lim = c(0, 6), y_brk = seq(0, 6, 1), 
  model = mod_geros
)

# 4) Net Production (mod_nprod) - Bottom plot
p_strs_nprod <- build_plot(
  spatial_Imean, "NET_PROD_KG_M2_YR", "se.NET_PROD_KG_M2_YR", 
  y_label = lbl_nprod, y_lim = c(-6, 6), y_brk = seq(-6, 6, 2), 
  add_zero_line = TRUE, model = mod_nprod
)

# Arrange 4x1 Vertical Panel
fig4_strs_4panel <- ggarrange(
  no_x_title(p_strs_cov),
  no_x_title(p_strs_gprod),
  no_x_title(p_strs_geros),
  p_strs_nprod, # Bottom plot keeps x-axis title
  ncol = 1, nrow = 4, align = "v",
  labels = "auto",
  font.label = list(size = 10, face = "bold"),
  label.x = 0.01, label.y = 0.99
)

# Export Figure 4
ggsave(here("outputs", "Fig4_HI_2024_StRS_4Panel_Vertical.png"), fig4_strs_4panel, 
       width = 15.0, height = 22.0, units = "cm", dpi = 300, bg = "white")


# 8. BUILD & EXPORT FIGURE 5 (NCRMP VERTICAL ACCRETION PLOT) -------------

# Predicted SLR range benchmark layer (6 to 10 mm/yr by 2050)
slr_layer <- list(
  # Clean shaded band without boundary lines
  annotate("rect", xmin = 19, xmax = 30, ymin = hi_slr_min, ymax = hi_slr_max,
           fill = "#D55E00", alpha = 0.12),
  
  # Text label centered vertically within the shaded band (y = 8.0)
  annotate("text", x = 29.5, y = (hi_slr_min + hi_slr_max) / 2, 
           label = paste0("Predicted SLR by 2050 (", hi_slr_min, "\u2013", hi_slr_max, " mm/yr)"),
           color = "#B33F00", size = 2.7, fontface = "bold", hjust = 0)
)

# Build single-panel NCRMP plot using build_plot & mod_vacc GLM trend line
p_ncrmp_vacc <- build_plot(
  data          = spatial_Imean, 
  y_var         = "VERTICAL_ACCRETION_MM_YR", 
  y_se          = "se.VERTICAL_ACCRETION_MM_YR", 
  y_label       = lbl_vacc,  
  y_lim         = c(-4, 13), 
  y_brk         = seq(-4, 12, 2), 
  add_zero_line = TRUE,
  model         = mod_vacc
) + slr_layer

# Export Figure 5
ggsave(here("outputs", "Fig5_NCRMP_Accretion_SLR.png"), p_ncrmp_vacc,
       width = 15.0, height = 11.5, units = "cm", dpi = 300, bg = "white")