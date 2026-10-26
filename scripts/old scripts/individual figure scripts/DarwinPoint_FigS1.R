# ==============================================================================
# Script: FigureS1_Sites.R
# Hawaiian Archipelago - Top Row Insets & Main Hawaiian Islands
# Aesthetic: High-Contrast Palette (Cobalt Blue, Warm Coral, Radiant Gold)
# Outputs: outputs/FigureS1_Sites.png & outputs/FigureS1_Sites.tiff
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

# Restore missing .shx files automatically if only .shp is present
Sys.setenv(SHAPE_RESTORE_SHX = "YES")
sf_use_s2(FALSE)

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. SPATIAL & BATHYMETRY DATA PREPARATION
# ==============================================================================

load_360_sf <- function(file_name) {
  file_path <- here("data", "gis layers", file_name)
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
# Robust MHI Bathymetry Pipeline: Classified Vector Contours
# ------------------------------------------------------------------------------
asc_path <- here("data", "gis layers", "gebco_2026_n30.0_s18.0_w-180.0_e-153.0_ascii.asc")
tif_path <- here("data", "gis layers", "gebco_2026_cropped_native.tif")

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
  
  # Bin depths into discrete classifications
  rcl <- matrix(c(
    -12000, -3000, 1,
    -3000, -1500, 2,
    -1500,  -400, 3,
    -400,     0, 4
  ), ncol = 3, byrow = TRUE)
  
  bathy_binned <- classify(bathy_smooth, rcl, include.lowest = TRUE)
  names(bathy_binned) <- "depth_bin"
  
  # Convert to vector polygons for robust rendering
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
# 2. LOAD & PREPARE DATASETS
# ==============================================================================
prod_csv <- here("data", "HA_2024_StRS_net.csv")
if (!file.exists(prod_csv)) prod_csv <- here("HA_2024_StRS_net.csv")
if (file.exists(prod_csv)) {
  prod_data <- read.csv(prod_csv)
  prod_data$LONGITUDE <- ifelse(prod_data$LONGITUDE < 0, prod_data$LONGITUDE + 360, prod_data$LONGITUDE)
} else {
  prod_data <- data.frame(LONGITUDE = numeric(), LATITUDE = numeric())
}

golden_csv <- here("data", "HA_2024_StRS_net_sw.csv")
if (!file.exists(golden_csv)) golden_csv <- here("HA_2024_StRS_net_sw.csv")
if (file.exists(golden_csv)) {
  golden_data <- read.csv(golden_csv)
  golden_data$LONGITUDE <- ifelse(golden_data$LONGITUDE < 0, golden_data$LONGITUDE + 360, golden_data$LONGITUDE)
} else {
  golden_data <- data.frame(LONGITUDE = numeric(), LATITUDE = numeric())
}

fixed_csv <- here("data", "HA_2024_fixed_net.csv")
if (!file.exists(fixed_csv)) fixed_csv <- here("HA_2024_fixed_net.csv")
if (file.exists(fixed_csv)) {
  fixed_data <- read.csv(fixed_csv)
  fixed_data$LONGITUDE <- ifelse(fixed_data$LONGITUDE < 0, fixed_data$LONGITUDE + 360, fixed_data$LONGITUDE)
} else {
  fixed_data <- data.frame(LONGITUDE = numeric(), LATITUDE = numeric())
}

# ------------------------------------------------------------------------------
# High-Contrast Palette:
# - Cobalt Royal Blue (#2563EB): All StRS Sampling Sites
# - Warm Coral Orange (#ED6925): Shallow-Water StRS Subset (Matches Fig 2/S3 P90 accent)
# - Radiant Gold/Yellow (#FACC15): Fixed Monitoring Sites (Matches Fig 2/S2 scale)
# ------------------------------------------------------------------------------
color_std   <- "#2563EB"  # Cobalt Royal Blue
color_sw    <- "#ED6925"  # Warm Coral Orange
color_fixed <- "#FACC15"  # Radiant Gold/Yellow

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
  
  # Bathymetry contours
  if (!is.null(layer_100m)) p <- p + geom_sf(data = layer_100m, fill = "#CBD5E1", color = "#CBD5E1", linewidth = 0.10)
  if (!is.null(layer_30m))  p <- p + geom_sf(data = layer_30m,  fill = "#94A3B8", color = "#94A3B8", linewidth = 0.10)
  
  p <- p +
    # Dark slate land fill
    geom_sf(data = islands_sf, fill = "#1E293B", color = NA) +
    
    # White background halos for triangles ONLY
    geom_point(data = fixed_data,  aes(x = LONGITUDE, y = LATITUDE), shape = 17, color = "#FFFFFF", size = 2.0, alpha = 0.50) +
    
    # Primary sampling points
    geom_point(data = prod_data,   aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_std,   color = "#020617", stroke = 0.45, size = 1.25, alpha = 0.90) +
    geom_point(data = golden_data, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_sw,    color = "#020617", stroke = 0.45, size = 1.25, alpha = 1.0) +
    geom_point(data = fixed_data,  aes(x = LONGITUDE, y = LATITUDE), shape = 24, fill = color_fixed, color = "#020617", stroke = 0.45, size = 1.25, alpha = 1.0) +
    
    # Subfigure tag inside top-left corner
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
      panel.background = element_rect(fill = "#FFFFFF", color = NA),
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
mhi_prod   <- subset(prod_data, LONGITUDE > 198.0)
mhi_golden <- subset(golden_data, LONGITUDE > 198.0)
mhi_fixed  <- subset(fixed_data, LONGITUDE > 198.0)

main_map <- ggplot()

# Render strictly classified bathymetry
if (!is.null(mhi_bathy_sf)) {
  main_map <- main_map +
    geom_sf(data = mhi_bathy_sf, aes(fill = depth_bin), color = NA) +
    scale_fill_manual(
      values = c(
        "1" = "#FFFFFF", # > 3,000 m
        "2" = "#F1F5F9", # 1,500 – 3,000 m
        "3" = "#CBD5E1", # 400 – 1,500 m
        "4" = "#94A3B8"  # < 400 m
      ),
      na.value = "transparent"
    )
}

main_map <- main_map +
  # Dark slate land fill
  geom_sf(data = mhi_islands, fill = "#1E293B", color = NA) +
  
  # White background halos for triangles ONLY
  geom_point(data = mhi_fixed,  aes(x = LONGITUDE, y = LATITUDE), shape = 17, color = "#FFFFFF", size = 2.5, alpha = 0.50) +
  
  # Primary sampling points
  geom_point(data = mhi_prod,   aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_std,   color = "#020617", stroke = 0.45, size = 1.6, alpha = 0.90) +
  geom_point(data = mhi_golden, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_sw,    color = "#020617", stroke = 0.45, size = 1.6, alpha = 1.0) +
  geom_point(data = mhi_fixed,  aes(x = LONGITUDE, y = LATITUDE), shape = 24, fill = color_fixed, color = "#020617", stroke = 0.45, size = 1.6, alpha = 1.0) +
  
  geom_text(
    data = mhi_label_data, aes(x = lbl_lon, y = lbl_lat, label = lbl),
    size = 3.4, color = "#0F172A", fontface = "bold", family = "sans"
  ) +
  
  # Subfigure tag E inside top-left corner
  annotate(
    "text", 
    x = 198.9 + (6.6 * 0.02), 
    y = 22.6 - (3.8 * 0.03), 
    label = "E", fontface = "bold", size = 4.2, color = "#0F172A", hjust = 0, vjust = 1
  ) +
  
  # Scale Bar at Bottom-Left
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
    legend.position   = "none",
    plot.margin       = margin(2, 2, 2, 2, "pt")
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
# 6. EXPORT PNAS FORMATS
# ==============================================================================
pnas_width  <- 17.8
pnas_height <- 15.0

# PNG Export (600 DPI)
ggsave(
  filename = file.path(out_dir, "FigureS1_Sites.png"),
  plot     = final_figure, width = pnas_width, height = pnas_height,
  units    = "cm", dpi = 600, bg = "white"
)

# TIFF Export (600 DPI, LZW compression)
ggsave(
  filename = file.path(out_dir, "FigureS1_Sites.tiff"),
  plot     = final_figure, width = pnas_width, height = pnas_height,
  units    = "cm", dpi = 600, compression = "lzw", bg = "white"
)