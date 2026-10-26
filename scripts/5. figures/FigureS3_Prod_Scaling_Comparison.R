# ==============================================================================
# Script:  FigureS3_Prod_Scaling_Comparison.R
# Purpose: Figure S3 - Sensitivity Analysis Comparing Scaled (Darwin Point) vs. 
#          Unscaled Gross Carbonate Production Estimates across Hawaiian Islands
# Inputs:  - data/strs data/field data/HA_2024_StRS_Prod_SiteLevel.csv
#          - data/strs data/field data/HA_2024_StRS_Prod_SiteLevel_notscaled.csv
#          - data/supporting data/Sectors-Strata-Areas.csv
# Outputs: - outputs/supplementary figures/FigureS3_Prod_Scaling_Comparison.png
#          - outputs/supplementary figures/FigureS3_Prod_Scaling_Comparison.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(survey)
})

# ------------------------------------------------------------------------------
# 1. DIRECTORIES & PNAS AESTHETIC DEFINITIONS 
# ------------------------------------------------------------------------------
supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

color_unscaled <- "#D73027"
color_scaled   <- "#4575B4"
point_border   <- "#09090B"

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

theme_pnas <- function(base_size = 8, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "#CBD5E1", linewidth = 0.5),
      panel.background   = element_rect(fill = "#FFFFFF", color = NA),
      plot.background    = element_rect(fill = "#FFFFFF", color = NA),
      axis.line          = element_line(color = "#18181B", linewidth = 0.5),
      axis.ticks         = element_line(color = "#18181B", linewidth = 0.5),
      axis.ticks.length  = unit(0.12, "cm"),
      axis.text          = element_text(color = "#27272A", size = 7),
      axis.text.x        = element_text(color = "#27272A", size = 7.5, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title         = element_text(color = "#09090B", size = 8.5, face = "bold"),
      axis.title.y       = element_text(margin = margin(r = 5)),
      axis.title.x       = element_blank(),
      legend.title       = element_blank(),
      legend.text        = element_text(size = 7.5, color = "#09090B", face = "bold"),
      legend.position    = c(0.86, 0.88),
      legend.background  = element_blank(),
      legend.key.size    = unit(0.35, "cm"),
      plot.margin        = margin(t = 8, r = 10, b = 6, l = 8, unit = "pt")
    )
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "strs data", "field data"),
    here("data", "strs data", "processed data"),
    here("data", "supporting data"),
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
# 2. DATA INGESTION & ISLAND-LEVEL SURVEY-WEIGHTED SUMMARY
# ------------------------------------------------------------------------------

# Load raw datasets and sector strata areas
df_scaled   <- saferead_csv("HA_2024_StRS_Prod_SiteLevel.csv")
df_unscaled <- saferead_csv("HA_2024_StRS_Prod_SiteLevel_notscaled.csv")
sectors     <- saferead_csv("Sectors-Strata-Areas.csv") %>%
  select(SEC_NAME, REEF_ZONE, DEPTH_BIN, NH) %>%
  distinct()

# Helper function to compute design-weighted island means and standard errors (SE)
calc_survey_stats <- function(df, model_name) {
  weighted_df <- df %>%
    left_join(sectors, by = c("SEC_NAME", "REEF_ZONE", "DEPTH_BIN")) %>%
    mutate(
      n = n_distinct(SITE),
      sw = NH / n,
      .by = c(OBS_YEAR, SEC_NAME, REEF_ZONE, DEPTH_BIN, NH)
    ) %>%
    filter(!is.na(sw) & n > 1) %>%
    mutate(
      Strat_conc = paste(OBS_YEAR, ISLAND, SEC_NAME, REEF_ZONE, DEPTH_BIN, sep = "_"),
      ISLAND = as.character(ISLAND)
    )
  
  des <- svydesign(id = ~1, strata = ~Strat_conc, weights = ~sw, data = weighted_df)
  
  # svyby automatically returns PROD (mean) and se
  svyby(~PROD, ~ISLAND, des, svymean, na.rm = TRUE) %>%
    as_tibble() %>%
    rename(Mean = PROD, SE = se) %>%
    mutate(Model = model_name)
}

# Calculate weighted stats for both
stats_scaled   <- calc_survey_stats(df_scaled, "Latitude-Scaled")
stats_unscaled <- calc_survey_stats(df_unscaled, "Unscaled")

# Combine and format
island_stats <- bind_rows(stats_unscaled, stats_scaled) %>%
  mutate(PLOT_LABEL = format_islands(ISLAND))

# Extract latitudes for ordering (North to South)
island_lats <- df_scaled %>%
  group_by(ISLAND) %>%
  summarize(LATITUDE = mean(LATITUDE, na.rm = TRUE), .groups = "drop") %>%
  mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
  arrange(desc(LATITUDE))

island_order <- island_lats$PLOT_LABEL

island_stats <- island_stats %>%
  mutate(
    Model      = factor(Model, levels = c("Unscaled", "Latitude-Scaled")),
    PLOT_LABEL = factor(PLOT_LABEL, levels = island_order)
  )

# ------------------------------------------------------------------------------
# 3. FIGURE S3 PLOT GENERATION
# ------------------------------------------------------------------------------
fig_s3 <- ggplot(island_stats, aes(x = PLOT_LABEL, y = Mean, fill = Model)) +
  geom_col(
    position  = position_dodge(0.75),
    width     = 0.65,
    color     = point_border,
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = pmax(0, Mean - SE), ymax = Mean + SE),
    position  = position_dodge(0.75),
    width     = 0.22,
    linewidth = 0.5,
    color     = point_border
  ) +
  scale_y_continuous(
    limits = c(0, 6.0),
    breaks = seq(0, 6, 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = c(
      "Unscaled"        = color_unscaled,
      "Latitude-Scaled" = color_scaled
    )
  ) +
  labs(
    x = NULL,
    y = expression(bold("Gross production rate (kg" ~ m^{-2} ~ yr^{-1} * ")"))
  ) +
  coord_cartesian(clip = "off") +
  theme_pnas()

# ------------------------------------------------------------------------------
# 4. EXPORT (PNAS 2-COLUMN STANDARD SPECIFICATIONS: 17.8 CM WIDTH)
# ------------------------------------------------------------------------------
out_png  <- file.path(supp_dir, "FigureS3_Prod_Scaling_Comparison.png")
out_tiff <- file.path(supp_dir, "FigureS3_Prod_Scaling_Comparison.tiff")

ggsave(out_png,  fig_s3, width = 17.8, height = 9.8, units = "cm", dpi = 600, bg = "white")
ggsave(out_tiff, fig_s3, width = 17.8, height = 9.8, units = "cm", dpi = 600, compression = "lzw", bg = "white")

cat(sprintf("✅ Figure S3 successfully saved to:\n • %s\n • %s\n", out_png, out_tiff))