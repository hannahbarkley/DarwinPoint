# ==============================================================================
# Script: Figure1_LightNetProduction_PNAS.R
# Hawaiian Archipelago - MHI Right-Aligned, Large Widescreen Insets
# Aesthetic: Balanced Zoom (Clear Inset Boundaries & Distinct MHI Labels)
# Outputs:  outputs/Figure1_LightMap_PNAS.tiff & .pdf
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
# 1. SPATIAL DATA PREPARATION
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

# ==============================================================================
# 2. LOAD & PREPARE DATASETS
# ==============================================================================
prod_csv <- here("data", "HA_2024_StRS_net.csv")
if (file.exists(prod_csv)) {
  prod_data <- read.csv(prod_csv)
  prod_data$LONGITUDE <- ifelse(prod_data$LONGITUDE < 0, prod_data$LONGITUDE + 360, prod_data$LONGITUDE)
} else {
  prod_data <- data.frame(LONGITUDE = numeric(), LATITUDE = numeric())
}

golden_csv <- here("data", "HA_2024_StRS_net_sw.csv")
if (file.exists(golden_csv)) {
  golden_data <- read.csv(golden_csv)
  golden_data$LONGITUDE <- ifelse(golden_data$LONGITUDE < 0, golden_data$LONGITUDE + 360, golden_data$LONGITUDE)
} else {
  golden_data <- data.frame(LONGITUDE = numeric(), LATITUDE = numeric())
}

# ==============================================================================
# 3. BASE THEMES & AESTHETICS
# ==============================================================================

# Adjusted Label Coordinates for MHI
mhi_label_data <- data.frame(
  lbl = c("Kaua'i", "Ni'ihau", "O'ahu", "Moloka'i", "L\u0101na'i", "Maui", "Kaho'olawe", "Hawai'i"),
  lbl_lon = c(200.0, 199.35, 201.8, 203.0, 202.4, 203.9, 203.1, 203.6), 
  lbl_lat = c(22.45, 21.75,  21.85, 21.45, 20.70, 21.15, 20.05, 19.30)
)

light_theme_base <- theme_minimal(base_size = 8, base_family = "sans") +
  theme(
    plot.background      = element_rect(fill = "#FFFFFF", color = NA),
    panel.background     = element_rect(fill = "#FFFFFF", color = NA),
    panel.border         = element_rect(color = "#334155", fill = NA, linewidth = 0.5), 
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    axis.text            = element_blank(),
    axis.title           = element_blank(),
    axis.ticks           = element_blank(),
    legend.position      = "none",
    plot.margin          = margin(4, 4, 4, 4, "pt")
  )

# ==============================================================================
# 4. BUILD INSET MAPS
# ==============================================================================
create_inset <- function(lon_center_360, lat_center, span = 0.25, inset_title = "", layer_30m = NULL, layer_100m = NULL) {
  p <- ggplot()
  
  if (!is.null(layer_100m)) p <- p + geom_sf(data = layer_100m, fill = "#CBD5E1", color = "#CBD5E1", linewidth = 0.15)
  if (!is.null(layer_30m))  p <- p + geom_sf(data = layer_30m,  fill = "#94A3B8", color = "#94A3B8", linewidth = 0.15)
  
  p <- p +
    geom_sf(data = islands_sf, fill = "#F1F5F9", color = "#334155", linewidth = 0.2) +
    
    geom_point(data = prod_data, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = "#3B82F6", color = "#0F172A", stroke = 0.35, size = 0.8, alpha = 0.90) +
    geom_point(data = golden_data, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = "#FFC300", color = "#0F172A", stroke = 0.35, size = 0.8, alpha = 1.0) +
    
    annotation_scale(
      location = "bl", width_hint = 0.25, 
      text_col = "transparent", line_width = 0.5, height = unit(0.03, "cm"),
      pad_x = unit(0.04, "npc"), pad_y = unit(0.04, "npc")
    ) +
    
    coord_sf(
      xlim  = c(lon_center_360 - span, lon_center_360 + span), 
      ylim  = c(lat_center - span, lat_center + span), 
      expand = FALSE, datum = NA 
    ) +
    labs(title = inset_title) +
    theme_void(base_size = 8, base_family = "sans") +
    theme(
      plot.background  = element_rect(fill = "#FFFFFF", color = NA), 
      panel.background = element_rect(fill = "#FFFFFF", color = NA),
      plot.title       = element_markdown(lineheight = 0.10, hjust = 0.5, margin = margin(t = 2, b = 2, l = 0, r = 0)),
      panel.border     = element_rect(color = "#0F172A", fill = NA, linewidth = 0.5),
      plot.margin      = margin(0, 0, 0, 0, "pt")
    )
  
  return(p)
}

# Compact Inset Titles (Hawaiian 6.5pt, English 5.0pt)
title_kure    <- "<span style='font-size: 6.5pt; color:#0F172A;'><b>H\u014dlanik\u016b</b></span><br><span style='font-size: 5.0pt; color:#475569;'><i>(Kure Atoll)</i></span>"
title_manawai <- "<span style='font-size: 6.5pt; color:#0F172A;'><b>Manawai</b></span><br><span style='font-size: 5.0pt; color:#475569;'><i>(Pearl & Hermes Atoll)</i></span>"
title_kapou   <- "<span style='font-size: 6.5pt; color:#0F172A;'><b>Kapou</b></span><br><span style='font-size: 5.0pt; color:#475569;'><i>(Lisianski Island)</i></span>"
title_lalo    <- "<span style='font-size: 6.5pt; color:#0F172A;'><b>Lalo</b></span><br><span style='font-size: 5.0pt; color:#475569;'><i>(French Frigate Shoals)</i></span>"

inset_kure    <- create_inset(181.64, 28.42, span = 0.15, inset_title = title_kure,    layer_30m = contour_30m, layer_100m = contour_100m)
inset_manawai <- create_inset(184.12, 27.84, span = 0.24, inset_title = title_manawai, layer_30m = contour_30m, layer_100m = contour_100m)
inset_lisian  <- create_inset(186.05, 26.05, span = 0.28, inset_title = title_kapou,   layer_30m = contour_30m, layer_100m = contour_100m) 
inset_ffs     <- create_inset(193.80, 23.75, span = 0.24, inset_title = title_lalo,    layer_30m = contour_30m, layer_100m = contour_100m)

# ==============================================================================
# 5. BUILD MAIN MAP & ASSEMBLE
# ==============================================================================
mhi_bbox <- st_bbox(c(xmin = 198.5, ymin = 15.0, xmax = 210.0, ymax = 28.0), crs = st_crs(islands_sf))
mhi_islands <- st_crop(islands_sf, mhi_bbox)
mhi_prod <- subset(prod_data, LONGITUDE > 198.5)
mhi_golden <- subset(golden_data, LONGITUDE > 198.5)

main_map <- ggplot() +
  geom_sf(data = mhi_islands, fill = "#F1F5F9", color = "#334155", linewidth = 0.25) +
  
  geom_point(data = mhi_prod, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = "#3B82F6", color = "#0F172A", stroke = 0.35, size = 0.85, alpha = 0.90) +
  geom_point(data = mhi_golden, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = "#FFC300", color = "#0F172A", stroke = 0.35, size = 0.85, alpha = 1.0) +
  
  geom_text(
    data = mhi_label_data, aes(x = lbl_lon, y = lbl_lat, label = lbl),
    size = 2.4, color = "#0F172A", fontface = "bold", family = "sans"
  ) +
  
  annotation_scale(
    location = "br", width_hint = 0.08, text_family = "sans", text_cex = 0.60,
    line_width = 0.5, height = unit(0.06, "cm"),
    pad_x = unit(0.04, "npc"), pad_y = unit(0.04, "npc")
  ) +
  
  # Balanced coordinate bounds giving room for Kaua'i and Ni'ihau clear of the insets
  coord_sf(
    xlim   = c(191, 206), 
    ylim   = c(18.3, 24.0),   
    expand = FALSE, datum = NA, clip = "on"
  ) +
  light_theme_base

# Assemble with slightly tightened inset bounds (right = 0.51) to keep clear of Ni'ihau
final_figure <- main_map + 
  inset_element(inset_kure,    left = 0.01, right = 0.15, bottom = 0.64, top = 0.96, align_to = "full") +
  inset_element(inset_manawai, left = 0.13, right = 0.27, bottom = 0.60, top = 0.92, align_to = "full") +
  inset_element(inset_lisian,  left = 0.25, right = 0.39, bottom = 0.56, top = 0.88, align_to = "full") +
  inset_element(inset_ffs,     left = 0.37, right = 0.51, bottom = 0.52, top = 0.84, align_to = "full")

# ==============================================================================
# 6. EXPORT PNAS FORMATS
# ==============================================================================
pnas_width  <- 17.8
pnas_height <- 8.9

ggsave(
  filename = file.path(out_dir, "FigureS1.tiff"),
  plot     = final_figure, width = pnas_width, height = pnas_height,
  units    = "cm", dpi = 600, compression = "lzw", bg = "white"
)

ggsave(
  filename = file.path(out_dir, "FigureS1.pdf"),
  plot     = final_figure, width = pnas_width, height = pnas_height,
  units    = "cm", device = cairo_pdf, bg = "white"
)

final_figure