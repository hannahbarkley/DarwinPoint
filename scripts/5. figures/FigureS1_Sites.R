# ==============================================================================
# Script:  FigureS1_Sites.R
# Purpose: Sampling Site Spatial Distribution across the Hawaiian Archipelago
#          - Geometry indicates survey type:
#              • Circles  (shape 21): StRS survey stations
#              • Triangles(shape 24): Fixed OCC survey stations
#          - Color indicates environmental subset:
#              • Orange (#FF8400): Grigg SW Fore-Reef subset
#              • Blue   (#1D70B8): Regional / Other stations
# Inputs:  data/gis layers/islands.shp
#          data/gis layers/srtm30v5h_30contr.shp
#          data/gis layers/srtm30v5h_100contr.shp
#          data/gis layers/gebco_2026_smoothed_360.tif (if present)
#          data/grigg data/HA_2024_Grigg_Equivalent_Sites.csv
#          data/strs data/processed data/HA_2024_StRS_Net_SiteLevel_cnet.csv
#          data/fixed data/processed data/HA_2024_Fixed_Net_SiteLevel_cnet.csv
# Outputs: outputs/supplementary figures/FigureS1_Sites.png
#          outputs/supplementary figures/FigureS1_Sites.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sf)
  library(terra)
  library(patchwork)
  library(ggspatial)
})

sf_use_s2(FALSE)
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

supp_dir <- here("outputs", "supplementary figures")
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# --- 1. UNIFIED COLOR & SYMBOL DEFINITIONS ---
color_grigg    <- "#FF8400"  # Vibrant Orange (Grigg SW fore-reef sites)
color_acc_blue <- "#1D70B8"  # High-contrast Cobalt Blue (Regional / Other sites)
point_border   <- "#0F172A"  # Deep slate outline for high contrast on all surfaces

# --- 2. FILE HELPERS ---
load_360_sf <- function(file_name) {
  file_path <- here("data", "gis layers", file_name)
  if (file.exists(file_path)) {
    sf_obj <- st_read(file_path, quiet = TRUE) %>% 
      st_zm(drop = TRUE, what = "ZM") %>% 
      st_transform(4326) %>% 
      st_make_valid()
    bb <- st_bbox(sf_obj)
    if (bb["xmin"] < 0 || bb["xmax"] < 0)
      sf_obj <- st_as_sf(terra::shift(terra::vect(sf_obj), dx = 360))
    if (any(grepl("LINE", st_geometry_type(sf_obj)))) {
      poly_sf <- suppressWarnings(st_polygonize(st_union(sf_obj)))
      poly_sf <- st_collection_extract(poly_sf, "POLYGON")
      if (length(poly_sf) > 0) sf_obj <- poly_sf
    }
    return(sf_obj)
  }
  return(NULL)
}

saferead_csv <- function(candidates) {
  search_dirs <- c(
    here("data", "grigg data"),
    here("data", "strs data", "processed data"),
    here("data", "fixed data", "processed data"),
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

# --- 3. DATA INGESTION & SUBSET PARTITIONING ---
islands_sf   <- load_360_sf("islands.shp")
contour_30m  <- load_360_sf("srtm30v5h_30contr.shp")
contour_100m <- load_360_sf("srtm30v5h_100contr.shp")

# 1. Grigg-Equivalent master dataset (both StRS and Fixed stations)
df_grigg_raw <- saferead_csv("HA_2024_Grigg_Equivalent_Sites.csv") %>%
  mutate(
    SITE      = coalesce(if ("SITE" %in% names(.)) SITE else NULL,
                         if ("OCC_SITEID" %in% names(.)) OCC_SITEID else NULL),
    LONGITUDE = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE)
  )

# Separate Grigg stations by survey type
df_strs_grigg  <- df_grigg_raw %>% 
  filter(SURVEY_TYPE == "StRS" | (!"SURVEY_TYPE" %in% names(.) & !grepl("^OCC", SITE)))

df_fixed_grigg <- df_grigg_raw %>% 
  filter(SURVEY_TYPE == "Fixed" | (!"SURVEY_TYPE" %in% names(.) & grepl("^OCC", SITE)))

# 2. Regional StRS stations (all StRS, isolated to non-Grigg)
df_strs_raw <- saferead_csv("HA_2024_StRS_Net_SiteLevel_cnet.csv") %>%
  mutate(
    SITE      = coalesce(if ("SITE" %in% names(.)) SITE else NULL,
                         if ("OCC_SITEID" %in% names(.)) OCC_SITEID else NULL),
    LONGITUDE = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE)
  )

df_strs_other <- df_strs_raw %>% filter(!SITE %in% df_strs_grigg$SITE)

# 3. Regional Fixed stations (all Fixed, isolated to non-Grigg)
df_fixed_raw <- saferead_csv("HA_2024_Fixed_Net_SiteLevel_cnet.csv") %>%
  mutate(
    SITE      = coalesce(if ("SITE" %in% names(.)) SITE else NULL,
                         if ("OCC_SITEID" %in% names(.)) OCC_SITEID else NULL),
    LONGITUDE = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE)
  )

df_fixed_other <- df_fixed_raw %>% filter(!SITE %in% df_fixed_grigg$SITE)

# --- 4. BATHYMETRY & BASEMAP PREPARATION ---
mhi_bbox    <- st_bbox(c(xmin = 198.0, ymin = 18.5, xmax = 206.0, ymax = 23.0), crs = st_crs(islands_sf))
mhi_islands <- st_crop(islands_sf, mhi_bbox)

mhi_bathy_sf    <- NULL
smooth_tif_path <- here("data", "gis layers", "gebco_2026_smoothed_360.tif")
if (file.exists(smooth_tif_path)) {
  bathy_smooth        <- rast(smooth_tif_path)
  mhi_ext_360         <- ext(198.0, 206.0, 18.0, 23.0)
  bathy_mhi           <- crop(bathy_smooth, mhi_ext_360)
  mhi_islands_v       <- vect(mhi_islands)
  mhi_immediate_mask  <- buffer(mhi_islands_v, width = 0.50)
  bathy_mhi_smooth    <- mask(bathy_mhi, mhi_immediate_mask)
  bathy_mhi_smooth[bathy_mhi_smooth > 0 & !is.na(bathy_mhi_smooth)] <- 0
  
  rcl <- matrix(c(-12000, -3000, 1, -3000, -1500, 2, -1500, -400, 3, -400, 0, 4), ncol = 3, byrow = TRUE)
  bathy_poly_v <- as.polygons(classify(bathy_mhi_smooth, rcl, include.lowest = TRUE))
  if (length(bathy_poly_v) > 0) {
    mhi_bathy_sf <- st_as_sf(bathy_poly_v)
    attr_col     <- setdiff(names(mhi_bathy_sf), attr(mhi_bathy_sf, "sf_column"))[1]
    if (!is.na(attr_col)) {
      names(mhi_bathy_sf)[names(mhi_bathy_sf) == attr_col] <- "depth_bin"
      mhi_bathy_sf$depth_bin <- as.character(mhi_bathy_sf$depth_bin)
    }
  }
}

# Base theme configured with automatic tag positioning inside the top-left of each panel
light_theme_base <- theme_minimal(base_size = 8, base_family = "sans") + 
  theme(
    plot.background   = element_rect(fill = "#FFFFFF", color = NA),
    panel.background  = element_rect(fill = "#FFFFFF", color = NA),
    panel.border      = element_rect(color = "#0F172A", fill = NA, linewidth = 0.6),
    panel.grid        = element_blank(),
    axis.text         = element_blank(),
    axis.title        = element_blank(),
    axis.ticks        = element_blank(),
    legend.position   = "none",
    plot.margin       = margin(0, 0, 0, 0, "pt"),
    plot.tag.position = c(0.03, 0.95),
    plot.tag          = element_text(face = "bold", size = 12, color = "#0F172A", hjust = 0, vjust = 1)
  )

# --- 5. INSET BUILDER ---
# target_ar = 1.0 (square insets) fills the 1.0-width grid cell perfectly
create_inset <- function(lon_c, lat_c, span = 0.20, l_30 = NULL, l_100 = NULL, target_ar = 1.0) {
  span_lon <- (span * target_ar) / cos(lat_c * pi / 180)
  
  p <- ggplot()
  
  if (!is.null(l_100)) p <- p + geom_sf(data = l_100, fill = "#F1F5F9", color = "#CBD5E1", linewidth = 0.15)
  if (!is.null(l_30))  p <- p + geom_sf(data = l_30,  fill = "#E2E8F0", color = "#94A3B8", linewidth = 0.15)
  
  p + 
    geom_sf(data = islands_sf, fill = "#1E293B", color = "#0F172A", linewidth = 0.2) +
    
    # 1. Base Layer: StRS Circles
    geom_point(data = df_strs_other, aes(x = LONGITUDE, y = LATITUDE), 
               shape = 21, fill = color_acc_blue, color = point_border, stroke = 0.35, size = 1.9, alpha = 0.90) +
    geom_point(data = df_strs_grigg, aes(x = LONGITUDE, y = LATITUDE), 
               shape = 21, fill = color_grigg,    color = point_border, stroke = 0.38, size = 1.9, alpha = 0.95) +
    
    # 2. Top Layer: Fixed Triangles
    geom_point(data = df_fixed_other, aes(x = LONGITUDE, y = LATITUDE), 
               shape = 24, fill = color_acc_blue, color = point_border, stroke = 0.40, size = 2.4, alpha = 0.95) +
    geom_point(data = df_fixed_grigg, aes(x = LONGITUDE, y = LATITUDE), 
               shape = 24, fill = color_grigg,    color = point_border, stroke = 0.42, size = 2.4, alpha = 0.95) +
    
    annotation_scale(location = "bl", width_hint = 0.35, style = "bar", bar_cols = c("#0F172A", "white"), 
                     text_col = "#0F172A", text_cex = 0.58, height = unit(0.04, "cm")) +
    coord_sf(xlim = c(lon_c - span_lon, lon_c + span_lon), 
             ylim = c(lat_c - span, lat_c + span), 
             expand = FALSE, datum = NA) + 
    light_theme_base
}

inset_kure    <- create_inset(181.66, 28.42, 0.11, contour_30m, contour_100m)
inset_manawai <- create_inset(184.13, 27.84, 0.20, contour_30m, contour_100m)
inset_lisian  <- create_inset(186.05, 26.05, 0.23, contour_30m, contour_100m)
inset_ffs     <- create_inset(193.79, 23.75, 0.20, contour_30m, contour_100m)

# --- 6. MAIN MHI MAP ---
main_map_s1 <- ggplot()

if (!is.null(mhi_bathy_sf)) {
  main_map_s1 <- main_map_s1 + 
    geom_sf(data = mhi_bathy_sf, aes(fill = depth_bin), color = NA) +
    scale_fill_manual(values = c("1" = "#FFFFFF", "2" = "#F8FAFC", "3" = "#F1F5F9", "4" = "#E2E8F0"), na.value = "transparent")
}

main_map_s1 <- main_map_s1 + 
  geom_sf(data = mhi_islands, fill = "#CBD5E1", color = "#B6BDC3", linewidth = 0.3) +
  
  # 1. Base Layer: StRS Circles
  geom_point(data = subset(df_strs_other, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), 
             shape = 21, fill = color_acc_blue, color = point_border, stroke = 0.40, size = 2.2, alpha = 0.90) +
  geom_point(data = subset(df_strs_grigg, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), 
             shape = 21, fill = color_grigg,    color = point_border, stroke = 0.42, size = 2.2, alpha = 0.95) +
  
  # 2. Top Layer: Fixed Triangles
  geom_point(data = subset(df_fixed_other, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), 
             shape = 24, fill = color_acc_blue, color = point_border, stroke = 0.45, size = 2.8, alpha = 0.95) +
  geom_point(data = subset(df_fixed_grigg, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), 
             shape = 24, fill = color_grigg,    color = point_border, stroke = 0.48, size = 2.8, alpha = 0.95) +
  
  geom_text(data = data.frame(
    lbl     = c("Kaua'i", "Ni'ihau", "O'ahu", "Moloka'i", "L\u0101na'i", "Maui", "Kaho'olawe", "Hawai'i"), 
    lbl_lon = c(200.15, 199.40, 202.05, 203.10, 202.70, 204.05, 202.90, 203.60), 
    lbl_lat = c(22.35,  21.75,  21.88,  21.38,  20.75,  21.05,  20.40,  19.25)
  ), aes(x = lbl_lon, y = lbl_lat, label = lbl), size = 3.6, color = "#0F172A", fontface = "bold") +
  annotation_scale(location = "bl", width_hint = 0.18, style = "bar", bar_cols = c("#0F172A", "white"), text_cex = 0.70, height = unit(0.05, "cm")) +
  coord_sf(xlim = c(198.9, 205.5), ylim = c(18.8, 22.6), expand = FALSE, datum = NA) + 
  light_theme_base

# --- 7. FINAL COMPOSITE LAYOUT & EXPORT ---
gap_w <- 0.06
gap_h <- 0.06
h_e   <- (4 + 3 * gap_w) / 1.6247

design_s1 <- c(
  patchwork::area(1, 1, 1, 1), # Inset A
  patchwork::area(1, 3, 1, 3), # Inset B
  patchwork::area(1, 5, 1, 5), # Inset C
  patchwork::area(1, 7, 1, 7), # Inset D
  patchwork::area(3, 1, 3, 7)  # Main map E
)

figS1_final <- inset_kure + inset_manawai + inset_lisian + inset_ffs + main_map_s1 + 
  plot_layout(
    design  = design_s1,
    widths  = c(1, gap_w, 1, gap_w, 1, gap_w, 1),
    heights = c(1, gap_h, h_e)
  ) + 
  plot_annotation(tag_levels = 'A')

out_png  <- file.path(supp_dir, "FigureS1_Sites.png")
out_tiff <- file.path(supp_dir, "FigureS1_Sites.tiff")

ggsave(out_png,  figS1_final, width = 17.8, height = 15.5, units = "cm", dpi = 600, bg = "white")
ggsave(out_tiff, figS1_final, width = 17.8, height = 15.5, units = "cm", dpi = 600, compression = "lzw", bg = "white")
