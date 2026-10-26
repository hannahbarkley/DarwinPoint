# ==============================================================================
# Script: Hawaiian_Darwin_Point_Figures.R
# Hawaiian Darwin Point Bathymetry, Map Figures & Coral Reef Comparisons
# Outputs: 
#   - outputs/Figure1_HA_map.png & .tiff
# ==============================================================================

# --- LOAD LIBRARIES ---
library(terra)
library(tidyverse)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(patchwork)
library(ggtext)
library(ggrepel)
library(ggpubr)
library(here)

# --- GLOBAL SETTINGS & DIRECTORY SETUP ---
sf_use_s2(FALSE) # Disable s2 spherical geometry to eliminate topology loop errors

out_dir <- here("outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ==============================================================================
# SECTION 1: FIGURE 1 - FULL MONUMENT MAP & BATHYMETRY
# ==============================================================================

# Base File Paths & Variables
asc_path        <- here("data", "gis layers", "gebco_2026_n30.0_s18.0_w-180.0_e-153.0_ascii.asc")
tif_path        <- here("data", "gis layers", "gebco_2026_cropped_native.tif")
smooth_tif_path <- here("data", "gis layers", "gebco_2026_smoothed_360.tif") # Cached raster

depth_breaks <- c(-12000, -3000, -1500, -400, Inf)
pmnm_label_2lines <- "Papahānaumokuākea National \nMarine Sanctuary Boundary"

# Science-Grade Orthographic Globe Base Setup
ortho_crs <- "+proj=ortho +lat_0=12 +lon_0=-160"

world_pacific <- ne_countries(scale = "small", returnclass = "sf") %>% 
  st_make_valid()

vis_box_east <- st_as_sfc(st_bbox(c(xmin = 110.5, ymin = -89.0, xmax = 180.0, ymax = 89.0), crs = 4326))
vis_box_west <- st_as_sfc(st_bbox(c(xmin = -180.0, ymin = -89.0, xmax = -70.5, ymax = 89.0), crs = 4326))

world_east <- suppressWarnings(st_crop(world_pacific, vis_box_east))
world_west <- suppressWarnings(st_crop(world_pacific, vis_box_west))

world_visible_wgs84 <- rbind(world_east, world_west) %>% st_make_valid()

ocean_circle <- st_sfc(st_point(c(0, 0)), crs = ortho_crs) %>% 
  st_buffer(6370000)

world_visible_ortho <- suppressWarnings(
  st_transform(world_visible_wgs84, crs = ortho_crs) %>% st_make_valid()
)
world_visible_ortho <- suppressWarnings(st_intersection(world_visible_ortho, ocean_circle))

# Load Islands
islands_sf <- st_read(here("data", "gis layers", "islands.shp"), quiet = TRUE) %>% 
  st_transform(4326) %>% 
  st_make_valid()

# Process Bathymetry (Check cache for computationally heavy focal smoothing)
if (file.exists(smooth_tif_path)) {
  bathy_smooth <- rast(smooth_tif_path)
} else {
  crop_box <- ext(-183.0, -152.0, 16.5, 32.5) 
  
  if (!file.exists(tif_path)) {
    bathy_raster <- rast(asc_path)
    crs(bathy_raster) <- "EPSG:4326"
    bathy_cropped <- crop(bathy_raster, crop_box)
    writeRaster(bathy_cropped, tif_path, overwrite = TRUE)
  } else {
    bathy_cropped <- rast(tif_path)
  }
  
  # Island & Feature targets for masking
  islands_v      <- vect(islands_sf)
  islands_simple <- simplifyGeom(islands_v, tolerance = 0.01)
  
  # Temp data structure for masking spatial feature points
  pts_df <- data.frame(
    lon = c(-178.3, -177.3, -175.8, -173.9, -171.7, -170.5, -167.9, -166.1, -164.7, -161.9, -159.5, -160.1, -158.0, -157.0, -156.9, -156.3, -156.6, -155.5),
    lat = c(28.4, 28.2, 27.8, 26.0, 25.7, 25.4, 25.0, 23.8, 23.5, 23.0, 22.1, 21.9, 21.4, 21.1, 20.8, 20.8, 20.5, 19.5)
  )
  feature_pts_v <- vect(pts_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  targets_v <- rbind(buffer(islands_simple, width = 0.25), buffer(feature_pts_v, width = 0.25))
  
  bathy_coarse  <- aggregate(bathy_cropped, fact = 5, fun = "max")
  ridge_binary  <- bathy_coarse >= -3200
  ridge_patches <- patches(ridge_binary, zeroAsNA = TRUE)
  
  intersect_ids <- unique(terra::extract(ridge_patches, targets_v, ID = FALSE)[[1]])
  intersect_ids <- intersect_ids[!is.na(intersect_ids)]
  
  ridge_mask_coarse <- match(ridge_patches, intersect_ids) > 0
  ridge_mask_vect   <- buffer(aggregate(as.polygons(ridge_mask_coarse, values = TRUE)), width = 0.06)
  
  bathy_cropped <- mask(bathy_cropped, ridge_mask_vect)
  bathy_cropped[bathy_cropped > 0] <- 0
  
  bathy_filled <- focal(bathy_cropped, w = 11, fun = mean, na.policy = "only", na.rm = TRUE)
  bathy_smooth <- focal(bathy_filled, w = 9, fun = mean, na.rm = TRUE)
  bathy_smooth[bathy_smooth > 0 & !is.na(bathy_smooth)] <- 0
  
  # Shift raster coordinates into 0..360 space & cache to disk
  bathy_smooth <- shift(bathy_smooth, dx = 360)
  writeRaster(bathy_smooth, smooth_tif_path, overwrite = TRUE)
}

# Shift islands into 0..360 space directly in sf
st_geometry(islands_sf) <- st_geometry(islands_sf) + c(360, 0)
st_crs(islands_sf)      <- 4326

# Load PMNM Extended Boundary
monument_raw <- st_read(here("data", "gis layers", "PMNM_Boundaries.shp"), quiet = TRUE) %>% 
  st_transform(4326) %>% 
  st_make_valid()

monument_sf <- monument_raw[toupper(monument_raw$Site_ID) == "MNM1", ]

# Direct vertex cleaning in 0..360 space
coords_raw <- st_coordinates(st_boundary(st_make_valid(monument_sf)))
coords_raw[, "X"] <- ifelse(coords_raw[, "X"] < 0, coords_raw[, "X"] + 360, coords_raw[, "X"])

segment_list <- list()
if (!is.null(coords_raw) && nrow(coords_raw) >= 2) {
  for (j in 1:(nrow(coords_raw) - 1)) {
    if ("L1" %in% colnames(coords_raw) && coords_raw[j, "L1"] != coords_raw[j + 1, "L1"]) next
    if ("L2" %in% colnames(coords_raw) && coords_raw[j, "L2"] != coords_raw[j + 1, "L2"]) next
    
    x1 <- coords_raw[j, "X"];     y1 <- coords_raw[j, "Y"]
    x2 <- coords_raw[j + 1, "X"]; y2 <- coords_raw[j + 1, "Y"]
    
    is_horizontal_25 <- (abs(y1 - 25.33) < 0.25) && (abs(y2 - 25.33) < 0.25) && (abs(x1 - x2) > 0.05)
    is_horizontal_31 <- (abs(y1 - 31.47) < 0.15) && (abs(y2 - 31.47) < 0.15) && (abs(x1 - x2) > 0.5)
    
    if (!is_horizontal_25 && !is_horizontal_31) {
      segment_list[[length(segment_list) + 1]] <- st_linestring(matrix(c(x1, x2, y1, y2), ncol = 2))
    }
  }
}

monument_sf_clean <- st_sf(geometry = st_sfc(st_line_merge(st_union(st_sfc(segment_list, crs = 4326))), crs = 4326))

# Define Label Data (Standardized Hawaiian Orthography with \u02BB ʻOkina)
label_data <- data.frame(
  lon = c(-178.3, -177.3, -175.8, -173.9, -171.7, -170.5, -167.9, -166.1, -164.7, -161.9, 
          -159.5, -160.1, -158.0, -157.0, -156.9, -156.3, -156.6, -155.5),
  lat = c(28.4,   28.2,   27.8,   26.0,   25.7,   25.4,   25.0,   23.8,   23.5,   23.0,   
          22.1,   21.9,   21.4,   21.1,   20.8,   20.8,   20.5,   19.5),
  
  lbl_lon = c(-178.3, -177.3, -175.8, -174.5, -171.7, -170.3, -167.9, -166.1, -164.9, -161.9, 
              -159.5, -160.7, -158.0, -156.3, -157.7, -155.5, -157.6, -156.8),
  lbl_lat = c(29.65,  27.1,   29.15,  25.0,   26.8,   23.8,   26.55,  22.6,   24.60,  24.20,   
              22.7,   21.5,   22.0,   21.5,   20.75,  21.1,   20.3,   19.4),
  
  seg_xend = c(-178.3, -177.3, -175.8, -174.4, -171.7, -170.3, -167.9, -166.1, -164.9, -161.9, 
               NA, NA, NA, NA, NA, NA, NA, NA),
  seg_yend = c(29.15,  27.60,  28.65,  25.50,  26.30,  24.30,  26.05,  23.10,  24.10,  23.70,
               NA, NA, NA, NA, NA, NA, NA, NA),
  
  haw = c("H\u014dlanik\u016b", "Kuaihelani", "Manawai", "Kapou", "Kamole", "Kamokuokamohoali\u02BBi", 
          "\u02BB\u014cn\u016bnui, \u02BB\u014cn\u016biki", "Lalo", "Mokumanamana", "Nihoa", 
          "Kaua\u02BBi", "Ni\u02BBihau", "O\u02BBahu", "Moloka\u02BBi", "L\u0101na\u02BBi", "Maui", "Kaho\u02BBolawe", "Hawai\u02BBi"),
  
  eng = c("Kure Atoll", "Midway Atoll", "Pearl & Hermes Atoll", "Lisianski Island", 
          "Laysan Island", "Maro Reef", "Gardner Pinnacles", "French Frigate Shoals", 
          "Necker Island", "", "", "", "", "", "", "", "", ""),
  
  draw_line = c(rep(TRUE, 10), rep(FALSE, 8))
) %>%
  mutate(
    full_label = ifelse(
      eng != "", 
      paste0("<span style='font-size: 5.1pt; color:#0F172A;'><b>", haw, "</b></span><br><span style='font-size: 4.2pt; color:#475569;'><i>(", eng, ")</i></span>"), 
      paste0("<span style='font-size: 5.1pt; color:#0F172A;'><b>", haw, "</b></span>")
    ),
    lon      = ifelse(lon < 0, lon + 360, lon),
    lbl_lon  = ifelse(lbl_lon < 0, lbl_lon + 360, lbl_lon),
    seg_xend = ifelse(!is.na(seg_xend) & seg_xend < 0, seg_xend + 360, seg_xend)
  )

# Globe inset
highlight_bbox    <- st_as_sfc(st_bbox(c(xmin = -183.0, ymin = 17.0, xmax = -152.0, ymax = 32.5), crs = 4326))
highlight_ortho   <- suppressWarnings(st_transform(highlight_bbox, crs = ortho_crs) %>% st_make_valid())
highlight_visible <- suppressWarnings(st_intersection(highlight_ortho, ocean_circle))

globe_inset <- ggplot() +
  geom_sf(data = ocean_circle, fill = "#FAFAFA", color = "#888888", linewidth = 0.05) +
  geom_sf(data = world_visible_ortho, fill = "#1E293B", color = NA) +
  geom_sf(data = highlight_visible, fill = "#0F172A", color = "#020617", alpha = 0.35, linewidth = 0.20) +
  coord_sf(crs = ortho_crs) +
  theme_void() + 
  theme(
    panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.15),
    plot.background  = element_blank(),
    panel.background = element_blank(),
    plot.margin      = margin(0, 0, 0, 0)
  )

main_map <- ggplot() +
  geom_spatraster_contour_filled(
    data   = bathy_smooth, 
    breaks = depth_breaks
  ) +
  scale_fill_manual(
    name   = "Depth (m)",
    values = c("#FFFFFF", "#F1F5F9", "#CBD5E1", "#94A3B8"),
    labels = c("> 3,000", "1,500 – 3,000", "400 – 1,500", "< 400"),
    na.value = "transparent",
    guide  = guide_legend(
      order          = 1,
      title.position = "top",
      title.hjust    = 0,
      ncol           = 1,
      byrow          = TRUE,
      reverse        = TRUE,
      override.aes   = list(color = "#0F172A", linewidth = 0.25)
    )
  ) +
  geom_sf(
    data = monument_sf_clean, 
    aes(color = pmnm_label_2lines), 
    linewidth = 0.40,
    linetype  = "dotted"
  ) +
  scale_color_manual(
    name   = NULL,
    values = setNames("#1E293B", pmnm_label_2lines),
    guide  = guide_legend(
      order        = 2,
      override.aes = list(linetype = "dotted", linewidth = 0.45, color = "#1E293B")
    )
  ) +
  geom_sf(data = islands_sf, fill = "#0F172A", color = "#0F172A", linewidth = 0.05) +
  geom_segment(
    data = subset(label_data, draw_line),
    aes(x = lon, y = lat, xend = seg_xend, yend = seg_yend), 
    color = "#334155", linewidth = 0.22, alpha = 0.90, lineend = "round"
  ) +
  geom_richtext(
    data = label_data, 
    aes(x = lbl_lon, y = lbl_lat, label = full_label),
    size = 1.9, lineheight = 0.88, color = "black", 
    family = "sans", hjust = 0.5, 
    fill = NA, label.color = NA,
    label.padding = unit(1, "pt"), label.r = unit(1, "pt")
  ) +
  coord_sf(
    xlim   = c(177.0, 208.0), 
    ylim   = c(17.0, 32.5), 
    expand = FALSE, 
    crs    = "+proj=longlat +datum=WGS84 +over",
    clip   = "on"
  ) +
  scale_x_continuous(
    breaks   = seq(180, 208, by = 4),
    labels   = function(x) ifelse(x == 180, "180°", ifelse(x > 180, paste0(360 - x, "°W"), paste0(x, "°E"))),
    sec.axis = dup_axis(name = NULL)
  ) +
  scale_y_continuous(
    breaks   = seq(18, 32, by = 2),
    labels   = function(x) paste0(x, "°N"),
    sec.axis = dup_axis(name = NULL)
  ) +
  theme_minimal(base_size = 6.5, base_family = "sans") +
  theme(
    plot.background      = element_rect(fill = "#FFFFFF", color = NA),
    panel.background     = element_rect(fill = "#FFFFFF", color = NA),
    panel.border         = element_rect(color = "#0F172A", fill = NA, linewidth = 0.45),
    axis.line            = element_blank(),
    axis.ticks           = element_line(color = "#0F172A", linewidth = 0.30),
    axis.ticks.length    = unit(0.08, "cm"),
    panel.grid.major     = element_line(color = "#EEF2F6", linewidth = 0.20, linetype = "solid"),
    panel.grid.minor     = element_blank(),
    legend.position      = c(0.040, 0.035),
    legend.justification = c(0, 0),
    legend.box           = "vertical",
    legend.box.background= element_rect(fill = alpha("#FFFFFF", 0.90), color = "#0F172A", linewidth = 0.35),
    legend.background    = element_blank(),
    legend.box.margin    = margin(4, 6, 4, 6, "pt"),
    legend.margin        = margin(0, 0, 0, 0, "pt"),
    legend.key.size      = unit(0.30, "cm"),
    legend.key           = element_rect(fill = NA, color = NA),
    legend.title         = element_text(size = 5.5, face = "bold", color = "#0F172A"),
    legend.text          = element_text(size = 5.0, color = "#0F172A", lineheight = 0.85),
    legend.spacing.y     = unit(1.5, "pt"),
    axis.text            = element_text(color = "#0F172A", size = 6.0, face = "bold"),
    axis.title           = element_blank(),
    plot.margin          = margin(8, 24, 8, 28, "pt")
  )

fig1_full_monument <- main_map + inset_element(globe_inset, left = 0.76, bottom = 0.68, right = 1.00, top = 0.98)

# Export Figure 1 (Full Monument Version)
ggsave(file.path(out_dir, "Figure1_HA_map.png"), fig1_full_monument, width = 17.8, height = 9.6, units = "cm", dpi = 600)
ggsave(file.path(out_dir, "Figure1_HA_map.tiff"), fig1_full_monument, width = 17.8, height = 9.6, units = "cm", dpi = 600, compression = "lzw")
