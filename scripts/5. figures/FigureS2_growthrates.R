# ==============================================================================
# Script:  FigureS2_growthrates.R
# Purpose: Coral Growth Rates and Estimated Production vs. Latitude
# Inputs:  data/dbase/prod_dbase_ha_strs_coeffs.csv
#          data/dbase/prod_dbase_ha_strs.csv (or prod_dbase_ha.csv)
# Outputs: outputs/supplementary figures/FigureS2_growthrates.png
#          outputs/supplementary figures/FigureS2_growthrates.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(ggpubr)
})

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

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
      plot.margin       = margin(t = 10, r = 8, b = 4, l = 14, unit = "pt")
    )
}

no_x_title <- function(p) p + labs(x = NULL) + theme(axis.title.x = element_blank())

saferead_csv <- function(candidates) {
  search_dirs <- c(here("data", "dbase"), here("data"), here())
  for (f in candidates) {
    for (d in search_dirs) {
      target <- file.path(d, f)
      if (file.exists(target)) return(read_csv(target, show_col_types = FALSE))
    }
  }
  stop(paste("None of the candidate files found:", paste(candidates, collapse = " | ")))
}

# ------------------------------------------------------------------------------
# 1. DATA INGESTION & PROCESSING
# ------------------------------------------------------------------------------
coeffs_dbase_s2 <- saferead_csv("prod_dbase_ha_strs_coeffs.csv")
raw_dbase_s2    <- saferead_csv(c("prod_dbase_ha_strs.csv", "prod_dbase_ha.csv"))

target_codes <- c("MOSP-EM", "POCS-BR", "POSP-BR", "POSP-EM", "POSP-MD")

group_labels_s2 <- c(
  expression("Encrusting "*italic("Montipora")), 
  expression("Branching "*italic("Pocillopora")), 
  expression("Branching "*italic("Porites")), 
  expression("Encrusting "*italic("Porites")), 
  expression("Massive "*italic("Porites"))
)

median_lengths <- c(
  "MOSP-EM" = 8.56,
  "POCS-BR" = 11.90,
  "POSP-BR" = 16.64,
  "POSP-EM" = 5.66,
  "POSP-MD" = 11.23
)

growth_df_s2 <- coeffs_dbase_s2 %>%
  mutate(
    SUBSTRATE_CODE = if_else(SUBSTRATE_CODE %in% c("MOEN", "MONTI", "MCAP"), "MOSP", SUBSTRATE_CODE),
    MORPHOLOGYCODE = if_else(MORPHOLOGYCODE %in% c("ENC", "EN"), "EM", MORPHOLOGYCODE),
    CORAL_GROUP    = paste0(SUBSTRATE_CODE, "-", MORPHOLOGYCODE)
  ) %>%
  filter(CORAL_GROUP %in% target_codes) %>%
  left_join(
    raw_dbase_s2 %>% 
      mutate(
        SUBSTRATE_CODE = if_else(SUBSTRATE_CODE %in% c("MOEN", "MONTI", "MCAP"), "MOSP", SUBSTRATE_CODE),
        MORPHOLOGYCODE = if_else(MORPHOLOGYCODE %in% c("ENC", "EN"), "EM", MORPHOLOGYCODE)
      ) %>% 
      select(LOCATIONCODE, SUBSTRATE_CODE, MORPHOLOGYCODE, LATITUDE) %>% 
      distinct(),
    by = c("LOCATIONCODE", "SUBSTRATE_CODE", "MORPHOLOGYCODE")
  ) %>%
  group_by(CORAL_GROUP, MORPHOLOGYCODE, LOCATIONCODE, LATITUDE) %>%
  summarize(
    Extension_cm_yr = mean(EXTENSION_CM_YR, na.rm = TRUE),
    Density_g_cm3   = mean(DENSITY_G_CM3, na.rm = TRUE),
    Coef_Mean       = mean(COEFFICIENT_MEAN, na.rm = TRUE),
    Intercept_Mean  = mean(INTERCEPT_MEAN, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  mutate(
    x           = median_lengths[as.character(CORAL_GROUP)],
    Colony_Prod = Coef_Mean * x + Intercept_Mean
  )

growth_df_s2$CORAL_GROUP <- factor(growth_df_s2$CORAL_GROUP, levels = target_codes)

# Preserved original color palette
group_pal_s2 <- c("#003087", "#E41A1C", "#1EBEC7", "#0085CA", "#FF8400")

growth_df_s2 <- growth_df_s2 %>%
  mutate(
    Ext_plot  = Extension_cm_yr + c("MOSP-EM" = -0.02, "POSP-EM" = 0.0, "POCS-BR" = 0.00, "POSP-BR" = 0.0, "POSP-MD" = 0.02)[as.character(CORAL_GROUP)],
    Den_plot  = Density_g_cm3   + c("MOSP-EM" = -0.02, "POSP-EM" = 0.0, "POCS-BR" = 0.02, "POSP-BR" = 0.0, "POSP-MD" = 0.02)[as.character(CORAL_GROUP)],
    Prod_plot = Colony_Prod
  )

# ------------------------------------------------------------------------------
# 2. SUBPLOTS
# ------------------------------------------------------------------------------

# Panel A: Linear Extension Rate
p_ext_s2 <- ggplot(growth_df_s2, aes(x = LATITUDE, y = Ext_plot)) +
  geom_smooth(aes(color = CORAL_GROUP), method = "loess", se = FALSE, span = 0.8, linewidth = 0.75) +
  geom_point(aes(fill = CORAL_GROUP), shape = 21, color = "gray20", stroke = 0.6, size = 2.8, alpha = 0.85) +
  scale_color_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_fill_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_x_reverse(
    limits = c(30.2, 18.6), 
    breaks = seq(19, 30, by = 1), 
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2.4, 0.4)) +
  coord_cartesian(ylim = c(0, 2.4), clip = "off") +
  labs(x = "Latitude (°N)", y = expression(bold("Linear extension rate (cm" ~ yr^{ -1 } * ")"))) +
  theme_publication() +
  theme(legend.position = "none")

# Panel B: Skeletal Density
p_den_s2 <- ggplot(growth_df_s2, aes(x = LATITUDE, y = Den_plot)) +
  geom_smooth(aes(color = CORAL_GROUP), method = "loess", se = FALSE, span = 0.8, linewidth = 0.75) +
  geom_point(aes(fill = CORAL_GROUP), shape = 21, color = "gray20", stroke = 0.6, size = 2.8, alpha = 0.85) +
  scale_color_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_fill_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_x_reverse(
    limits = c(30.2, 18.6), 
    breaks = seq(19, 30, by = 1), 
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(limits = c(1.2, 1.8), breaks = seq(1.2, 1.8, 0.1)) +
  coord_cartesian(ylim = c(1.2, 1.8), clip = "off") +
  labs(x = "Latitude (°N)", y = expression(bold("Skeletal density (g" ~ cm^{ -3 } * ")"))) +
  theme_publication() +
  theme(legend.position = "none")

# Panel C: Colony Production Rate
p_prod_s2 <- ggplot(growth_df_s2, aes(x = LATITUDE, y = Prod_plot)) +
  geom_smooth(aes(color = CORAL_GROUP), method = "loess", se = FALSE, span = 0.8, linewidth = 0.75) +
  geom_point(aes(fill = CORAL_GROUP), shape = 21, color = "gray20", stroke = 0.6, size = 2.8, alpha = 0.85) +
  scale_color_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_fill_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_x_reverse(
    limits = c(30.2, 18.6), 
    breaks = seq(19, 30, by = 1), 
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0, 25, 5)) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  labs(x = "Latitude (°N)", y = expression(bold("Colony production rate (g" ~ yr^{ -1 } * ")"))) +
  theme_publication() +
  theme(legend.position = "none")

# Shared Right-Side Legend
p_legend_s2_growth <- get_legend(
  ggplot(growth_df_s2, aes(x = LATITUDE, y = Den_plot, fill = CORAL_GROUP)) +
    geom_point(shape = 21, color = "gray20", stroke = 0.6, size = 2.8) +
    scale_fill_manual(values = group_pal_s2, labels = group_labels_s2, name = NULL) +
    guides(fill = guide_legend(
      ncol = 1, 
      override.aes = list(shape = 21, size = 2.8, stroke = 0.6, color = "gray20")
    )) +
    theme_publication() +
    theme(
      legend.position = "right", 
      legend.title    = element_blank(), 
      legend.text     = element_text(size = 8.5, color = "gray20")
    )
)

# ------------------------------------------------------------------------------
# 3. ASSEMBLY & EXPORT
# ------------------------------------------------------------------------------
main_plots_s2 <- ggarrange(
  no_x_title(p_ext_s2),
  no_x_title(p_den_s2),
  p_prod_s2,
  ncol       = 1, 
  nrow       = 3, 
  align      = "v",
  labels     = c("A", "B", "C"),
  font.label = list(size = 11, face = "bold"),
  label.x    = 0.01, 
  label.y    = 1.02
)

figS2_final <- ggarrange(
  main_plots_s2, 
  p_legend_s2_growth, 
  ncol   = 2, 
  nrow   = 1, 
  widths = c(1, 0.3)
)

ggsave(file.path(supp_dir, "FigureS2_growthrates.png"), figS2_final, width = 17.5, height = 21.0, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS2_growthrates.tiff"), figS2_final, width = 17.5, height = 21.0, units = "cm", dpi = 300, compression = "lzw", bg = "white")
