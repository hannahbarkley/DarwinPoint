# ==============================================================================
# Script: FigureS4_Porosity_Map.R
# Hawaiian Archipelago - Site-Level Framework Porosity Spatial Map
# Aesthetic: Refined Marine Blue Bathymetry, High-Contrast Light Land Mass, Plasma Ramp
# Outputs: data/2024_Site_Estimated_Porosity_3Morph.csv
#          outputs/FigureS4_Porosity_Map.png & outputs/FigureS4_Porosity_Map.tiff
# ==============================================================================

library(terra)
library(ggplot2)
library(tidyterra)
library(scales)
library(sf)
library(patchwork)        
library(ggtext)
library(here)
library(ggspatial)
library(tidyverse)
library(ggnewscale)

# Restore missing .shx files automatically if only .shp is present
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
sf_use_s2(FALSE)

# Directory Routing: CSVs to data/, Figures to outputs/
data_dir <- here("data")
out_dir  <- here("outputs")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(out_dir))  dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. SPATIAL & BATHYMETRY DATA PREPARATION
# ==============================================================================

load_360_sf <- function(file_name) {
  file_path <- file.path(data_dir, "gis layers", file_name)
  if (file.exists(file_path)) {
    sf_obj <- st_read(file_path, quiet = TRUE) %>% 
      st_zm(drop = TRUE, what = "ZM") %>% 
      st_transform(4326) %>% 
      st_make_valid()
    
    bb <- st_bbox(sf_obj)
    if (bb["xmin"] < 0 || bb["xmax"] < 0) {
      v <- terra::vect(sf_obj)
      v_shifted <- terra::shift(v, dx = 360)
      sf_obj <- st_as_sf(v_shifted)
    }
    
    if (any(grepl("LINE", st_geometry_type(sf_obj)))) {
      poly_sf <- suppressWarnings(st_polygonize(st_union(sf_obj)))
      poly_sf <- st_collection_extract(poly_sf, "POLYGON")
      if (length(poly_sf) > 0) {
        sf_obj <- poly_sf
      }
    }
    return(sf_obj)
  } else {
    warning(paste(file_name, "not found in data/gis layers/"))
    return(NULL)
  }
}

islands_sf   <- load_360_sf("islands.shp")
contour_30m  <- load_360_sf("srtm30v5h_30contr.shp")
contour_100m <- load_360_sf("srtm30v5h_100contr.shp")

# ------------------------------------------------------------------------------
# MHI Bathymetry Pipeline: Classified Vector Contours
# ------------------------------------------------------------------------------
asc_path <- file.path(data_dir, "gis layers", "gebco_2026_n30.0_s18.0_w-180.0_e-153.0_ascii.asc")
tif_path <- file.path(data_dir, "gis layers", "gebco_2026_cropped_native.tif")

if (file.exists(tif_path)) {
  bathy_raw <- rast(tif_path)
} else if (file.exists(asc_path)) {
  bathy_raster <- rast(asc_path)
  crs(bathy_raster) <- "EPSG:4326"
  bathy_raw <- crop(bathy_raster, ext(-183.0, -152.0, 16.5, 32.5))
} else {
  bathy_raw <- NULL
}

mhi_bathy_sf <- NULL
mhi_bbox <- st_bbox(c(xmin = 198.0, ymin = 18.5, xmax = 206.0, ymax = 23.0), crs = st_crs(islands_sf))
mhi_islands <- st_crop(islands_sf, mhi_bbox)

if (!is.null(bathy_raw)) {
  bathy_360 <- shift(bathy_raw, dx = 360)
  mhi_ext_360 <- ext(198.0, 206.0, 18.0, 23.0)
  bathy_mhi <- crop(bathy_360, mhi_ext_360)
  
  mhi_islands_v <- vect(mhi_islands)
  mhi_immediate_mask <- buffer(mhi_islands_v, width = 0.50)
  
  bathy_masked <- mask(bathy_mhi, mhi_immediate_mask)
  bathy_masked[bathy_masked > 0] <- 0
  
  bathy_filled <- focal(bathy_masked, w = 11, fun = mean, na.policy = "only", na.rm = TRUE)
  bathy_smooth <- focal(bathy_filled, w = 9,  fun = mean, na.rm = TRUE)
  bathy_smooth <- mask(bathy_smooth, mhi_immediate_mask)
  bathy_smooth[bathy_smooth > 0 & !is.na(bathy_smooth)] <- 0
  
  rcl <- matrix(c(
    -12000, -3000, 1,
    -3000, -1500, 2,
    -1500,  -400, 3,
    -400,      0, 4
  ), ncol = 3, byrow = TRUE)
  
  bathy_binned <- classify(bathy_smooth, rcl, include.lowest = TRUE)
  names(bathy_binned) <- "depth_bin"
  
  bathy_poly_v <- as.polygons(bathy_binned)
  
  if (length(bathy_poly_v) > 0) {
    mhi_bathy_sf <- st_as_sf(bathy_poly_v)
    
    attr_col <- setdiff(names(mhi_bathy_sf), attr(mhi_bathy_sf, "sf_column"))[1]
    if (!is.na(attr_col)) {
      names(mhi_bathy_sf)[names(mhi_bathy_sf) == attr_col] <- "depth_bin"
      mhi_bathy_sf$depth_bin <- as.character(mhi_bathy_sf$depth_bin)
    }
  }
}

# ==============================================================================
# 2. LOAD & CALCULATE SITE ESTIMATED POROSITY
# ==============================================================================

cnet_csv   <- file.path(data_dir, "2024_Hawaii_Cnet_cover.csv")
lookup_csv <- file.path(data_dir, "Hawaii_Photoquad_codes.csv")

if (!file.exists(cnet_csv))   cnet_csv   <- here("2024_Hawaii_Cnet_cover.csv")
if (!file.exists(lookup_csv)) lookup_csv <- here("Hawaii_Photoquad_codes.csv")

if (!file.exists(cnet_csv))   stop("Dataset 2024_Hawaii_Cnet_cover.csv not found.")
if (!file.exists(lookup_csv)) stop("Lookup table Hawaii_Photoquad_codes.csv not found.")

cover_df  <- read_csv(cnet_csv, show_col_types = FALSE)
lookup_df <- read_csv(lookup_csv, show_col_types = FALSE) %>%
  select(CODE, TIER_2) %>%
  distinct(CODE, .keep_all = TRUE)

extra_lookup <- tibble(
  CODE   = c("HAL", "HYCO"),
  TIER_2 = c("UPMA", "HYCO")
)
lookup_df <- bind_rows(lookup_df, extra_lookup)

meta_cols <- c("SITEVISITID", "SITE", "ISLAND", "REGION", "REEF_ZONE", "DEPTH_BIN", "LATITUDE", "LONGITUDE")
meta_present <- intersect(meta_cols, names(cover_df))
code_cols    <- setdiff(names(cover_df), c(meta_present, "OBS_YEAR", "HABITAT_CODE", "DATE_", "MISSIONID", "OCC_SITEID", "ANALYSIS_YEAR", "ANALYSIS_SCHEME", "SEC_NAME", "EXCLUDE_FLAG", "TRANSECT_PHOTOS", "Oceanography", "new_MIN_DEPTH_M", "new_MAX_DEPTH_M", "METHOD", "PERM_SITE", "CLIMATE_STATION_YN", "N", "new.N"))

# Model End-Member Parameters
phi_massive    <- 0.30
phi_encrusting <- 0.25
phi_branching  <- 0.60
phi_default    <- 0.40

porosity_data <- cover_df %>%
  pivot_longer(cols = all_of(code_cols), names_to = "CODE", values_to = "COVER_PERC") %>%
  left_join(lookup_df, by = "CODE") %>%
  filter(TIER_2 %in% c("BR", "MASS", "ENC")) %>%
  group_by(across(all_of(meta_present)), TIER_2) %>%
  summarize(Total_Cover = sum(COVER_PERC, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = TIER_2, values_from = Total_Cover, values_fill = 0) %>%
  rename(Cover_Branching = BR, Cover_Massive = MASS, Cover_Encrusting = ENC) %>%
  mutate(
    Total_Coral = Cover_Branching + Cover_Massive + Cover_Encrusting,
    Rel_BR   = if_else(Total_Coral > 0, Cover_Branching / Total_Coral, 0),
    Rel_MASS = if_else(Total_Coral > 0, Cover_Massive / Total_Coral, 0),
    Rel_ENC  = if_else(Total_Coral > 0, Cover_Encrusting / Total_Coral, 0),
    Est_Porosity = if_else(
      Total_Coral > 0,
      (Rel_MASS * phi_massive) + (Rel_ENC * phi_encrusting) + (Rel_BR * phi_branching),
      phi_default
    ),
    LONGITUDE = if_else(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE)
  )

# Export CSV dataset directly to data/ directory
write_csv(porosity_data, file.path(data_dir, "2024_Site_Estimated_Porosity_3Morph.csv"))

# ==============================================================================
# 3. BASE THEMES & AESTHETICS
# ==============================================================================

mhi_label_data <- data.frame(
  lbl = c("Kaua'i", "Ni'ihau", "O'ahu", "Moloka'i", "L\u0101na'i", "Maui", "Kaho'olawe", "Hawai'i"),
  lbl_lon = c(200.28, 199.55, 202.05, 203.10, 202.66, 203.89, 202.95, 203.75), 
  lbl_lat = c(22.40,  21.65,  21.83,  21.39,  20.62,  21.12,  20.39,  19.30)
)

light_theme_base <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    plot.background      = element_rect(fill = "#FFFFFF", color = NA),
    panel.background     = element_rect(fill = "#FFFFFF", color = NA),
    panel.border         = element_rect(color = "#0F172A", fill = NA, linewidth = 0.6), 
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.text            = element_blank(),
    axis.title           = element_blank(),
    axis.ticks           = element_blank(),
    legend.position      = "none",
    plot.margin          = margin(2, 2, 2, 2, "pt")
  )

# ==============================================================================
# 4. BUILD INSET MAPS
# ==============================================================================

create_inset <- function(lon_center_360, lat_center, span = 0.20, layer_30m = NULL, layer_100m = NULL, tag_label = "A") {
  p <- ggplot()
  
  if (!is.null(layer_100m)) p <- p + geom_sf(data = layer_100m, fill = "#BAE6FD", color = "#BAE6FD", linewidth = 0.10)
  if (!is.null(layer_30m))  p <- p + geom_sf(data = layer_30m,  fill = "#7DD3FC", color = "#7DD3FC", linewidth = 0.10)
  
  p <- p +
    geom_sf(data = islands_sf, fill = "#E5E7EB", color = "#1F2937", linewidth = 0.30) +
    
    geom_point(
      data = porosity_data, 
      aes(x = LONGITUDE, y = LATITUDE, fill = Est_Porosity), 
      shape = 21, color = "#0F172A", stroke = 0.35, size = 2.0, alpha = 0.95
    ) +
    
    scale_fill_viridis_c(
      option = "plasma", 
      limits = c(0.25, 0.60), 
      labels = scales::percent_format(accuracy = 1),
      name = expression(bold("Framework Porosity ("*phi*")"))
    ) +
    
    annotate(
      "text", 
      x = lon_center_360 - span + (2 * span * 0.05), 
      y = lat_center + span - (2 * span * 0.05), 
      label = tag_label, fontface = "bold", size = 4.2, color = "#0F172A", hjust = 0, vjust = 1
    ) +
    
    annotation_scale(
      location = "bl", width_hint = 0.35, 
      style = "bar", bar_cols = c("#0F172A", "white"),
      text_col = "#0F172A", text_cex = 0.58, text_family = "sans",
      height = unit(0.04, "cm"), pad_x = unit(0.05, "npc"), pad_y = unit(0.05, "npc")
    ) +
    
    coord_sf(
      xlim   = c(lon_center_360 - span, lon_center_360 + span), 
      ylim   = c(lat_center - span, lat_center + span), 
      expand = FALSE, datum = NA 
    ) +
    theme_void(base_size = 8, base_family = "sans") +
    theme(
      plot.background  = element_rect(fill = "#FFFFFF", color = NA), 
      panel.background = element_rect(fill = "#F0F9FF", color = NA),
      panel.border     = element_rect(color = "#0F172A", fill = NA, linewidth = 0.6),
      legend.position  = "none",
      plot.margin      = margin(2, 2, 2, 2, "pt")
    )
  
  return(p)
}

inset_kure    <- create_inset(181.66, 28.42, span = 0.11, layer_30m = contour_30m, layer_100m = contour_100m, tag_label = "A")
inset_manawai <- create_inset(184.13, 27.84, span = 0.20, layer_30m = contour_30m, layer_100m = contour_100m, tag_label = "B")
inset_lisian  <- create_inset(186.05, 26.05, span = 0.23, layer_30m = contour_30m, layer_100m = contour_100m, tag_label = "C") 
inset_ffs     <- create_inset(193.79, 23.75, span = 0.20, layer_30m = contour_30m, layer_100m = contour_100m, tag_label = "D")

# ==============================================================================
# 5. BUILD MAIN MAP & ASSEMBLE
# ==============================================================================

mhi_porosity <- subset(porosity_data, LONGITUDE > 198.0)

main_map <- ggplot()

if (!is.null(mhi_bathy_sf)) {
  main_map <- main_map +
    geom_sf(data = mhi_bathy_sf, aes(fill = depth_bin), color = NA) +
    scale_fill_manual(
      values = c(
        "1" = "#FFFFFF",
        "2" = "#F0F9FF",
        "3" = "#E0F2FE",
        "4" = "#BAE6FD"
      ),
      guide = "none"
    ) +
    new_scale_fill()
}

main_map <- main_map +
  geom_sf(data = mhi_islands, fill = "#E5E7EB", color = "#1F2937", linewidth = 0.35) +
  
  geom_point(
    data = mhi_porosity, 
    aes(x = LONGITUDE, y = LATITUDE, fill = Est_Porosity), 
    shape = 21, color = "#0F172A", stroke = 0.35, size = 2.1, alpha = 0.95
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    limits = c(0.25, 0.60),
    breaks = seq(0.25, 0.60, 0.05),
    labels = c("25%", "30%", "35%", "40%", "45%", "50%", "55%", "60%"),
    name   = expression(bold("Estimated Framework Porosity ("*phi*")")),
    guide  = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5.6, "cm"),  # Calibrated width for exact label enclosure
      barheight      = unit(0.35, "cm"),
      frame.colour   = "#1F2937",
      frame.linewidth= 0.5,
      ticks.colour   = "#1F2937"
    )
  ) +
  
  geom_text(
    data = mhi_label_data, aes(x = lbl_lon, y = lbl_lat, label = lbl),
    size = 3.4, color = "#0F172A", fontface = "bold", family = "sans"
  ) +
  
  annotate(
    "text", 
    x = 198.9 + (6.6 * 0.02), 
    y = 22.6 - (3.8 * 0.03), 
    label = "E", fontface = "bold", size = 4.2, color = "#0F172A", hjust = 0, vjust = 1
  ) +
  
  annotation_scale(
    location = "bl", width_hint = 0.18, style = "bar", bar_cols = c("#0F172A", "white"),
    text_family = "sans", text_cex = 0.70, text_col = "#0F172A",
    height = unit(0.05, "cm"), pad_x = unit(0.04, "npc"), pad_y = unit(0.04, "npc")
  ) +
  
  coord_sf(
    xlim   = c(198.9, 205.5), 
    ylim   = c(18.8, 22.6),   
    expand = FALSE, datum = NA
  ) +
  light_theme_base +
  theme(
    legend.position      = c(0.78, 0.88),  # Shifted rightward into open ocean
    legend.direction     = "horizontal",
    legend.background    = element_rect(fill = alpha("white", 0.95), color = "#CBD5E1", linewidth = 0.5),
    legend.margin        = margin(t = 6, r = 18, b = 6, l = 18, unit = "pt"), # Expanded margins enclose 25% and 60% labels
    legend.title         = element_text(size = 8, face = "bold"),
    legend.text          = element_text(size = 7.5, face = "bold"),
    plot.margin          = margin(2, 2, 2, 2, "pt")
  )

# Layout design for patchwork
layout_design <- c(
  area(t = 1, l = 1, b = 1, r = 1),
  area(t = 1, l = 2, b = 1, r = 2),
  area(t = 1, l = 3, b = 1, r = 3),
  area(t = 1, l = 4, b = 1, r = 4),
  area(t = 2, l = 1, b = 2, r = 4)
)

final_figure <- inset_kure + inset_manawai + inset_lisian + inset_ffs + main_map +
  plot_layout(design = layout_design, heights = c(1, 2.5))

# ==============================================================================
# 6. EXPORT FIGURES
# ==============================================================================
pnas_width  <- 17.8
pnas_height <- 15.0

# PNG Export (600 DPI)
ggsave(
  filename = file.path(out_dir, "HA_Porosity_Map.png"),
  plot     = final_figure, width = pnas_width, height = pnas_height,
  units    = "cm", dpi = 600, bg = "white"
)

# TIFF Export (600 DPI, LZW compression)
ggsave(
  filename    = file.path(out_dir, "HA_Porosity_Map.tiff"),
  plot        = final_figure, width = pnas_width, height = pnas_height,
  units       = "cm", dpi = 600, compression = "lzw", bg = "white"
)

final_figure