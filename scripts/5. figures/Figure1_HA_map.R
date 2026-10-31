# ==============================================================================
# Script:  Figure1_HA_map.R
# Purpose: Full Hawaiian Archipelago Bathymetric Map and Global Inset
# Inputs:  data/gis layers/islands.shp
#          data/gis layers/PMNM_Boundaries.shp
#          data/gis layers/gebco_2026_n30.0_s18.0_w-180.0_e-153.0_ascii.asc
#          data/gis layers/gebco_2026_smoothed_360.tif (preprocessed cache, if present)
#          data/gis layers/gebco_2026_cropped_native.tif (intermediate raster, if present)
# Outputs: outputs/main figures/Figure1_HA_map.png
#          outputs/main figures/Figure1_HA_map.tiff
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(sf)
  library(terra)
  library(tidyterra)
  library(rnaturalearth)
  library(ggtext)
  library(patchwork)
  library(geomtextpath)
})

sf_use_s2(FALSE)
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

main_dir <- here("outputs", "main figures")
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

ortho_crs <- "+proj=ortho +lat_0=12 +lon_0=-160"
world_pacific <- ne_countries(scale = "small", returnclass = "sf") %>% st_make_valid()
vis_box_east <- st_as_sfc(st_bbox(c(xmin = 110.5, ymin = -89.0, xmax = 180.0, ymax = 89.0), crs = 4326))
vis_box_west <- st_as_sfc(st_bbox(c(xmin = -180.0, ymin = -89.0, xmax = -70.5, ymax = 89.0), crs = 4326))
world_visible_wgs84 <- rbind(
  suppressWarnings(st_crop(world_pacific, vis_box_east)),
  suppressWarnings(st_crop(world_pacific, vis_box_west))
) %>% st_make_valid()
ocean_circle <- st_sfc(st_point(c(0, 0)), crs = ortho_crs) %>% st_buffer(6370000)
world_visible_ortho <- suppressWarnings(st_intersection(st_transform(world_visible_wgs84, crs = ortho_crs) %>% st_make_valid(), ocean_circle))

islands_f1_sf <- st_read(here("data", "gis layers", "islands.shp"), quiet = TRUE) %>% st_transform(4326) %>% st_make_valid()

asc_path <- here("data", "gis layers", "gebco_2026_n30.0_s18.0_w-180.0_e-153.0_ascii.asc")
tif_path <- here("data", "gis layers", "gebco_2026_cropped_native.tif")
smooth_tif_path <- here("data", "gis layers", "gebco_2026_smoothed_360.tif")

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
  
  pts_df <- data.frame(
    lon = c(-178.3, -177.3, -175.8, -173.9, -171.7, -170.5, -167.9, -166.1, -164.7, -161.9, -159.5, -160.1, -158.0, -157.0, -156.9, -156.3, -156.6, -155.5),
    lat = c(28.4, 28.2, 27.8, 26.0, 25.7, 25.4, 25.0, 23.8, 23.5, 23.0, 22.1, 21.9, 21.4, 21.1, 20.8, 20.8, 20.5, 19.5)
  )
  feature_pts_v <- vect(pts_df, geom = c("lon", "lat"), crs = "EPSG:4326")
  targets_v <- rbind(buffer(simplifyGeom(vect(islands_f1_sf), tolerance = 0.01), width = 0.25),
                     buffer(feature_pts_v, width = 0.25))
  
  bathy_coarse <- aggregate(bathy_cropped, fact = 5, fun = "max")
  ridge_patches <- patches(bathy_coarse >= -3200, zeroAsNA = TRUE)
  intersect_ids <- unique(terra::extract(ridge_patches, targets_v, ID = FALSE)[[1]])
  intersect_ids <- intersect_ids[!is.na(intersect_ids)]
  
  bathy_cropped <- mask(bathy_cropped, buffer(aggregate(as.polygons(match(ridge_patches, intersect_ids) > 0, values = TRUE)), width = 0.06))
  bathy_cropped[bathy_cropped > 0] <- 0
  bathy_smooth <- shift(focal(focal(bathy_cropped, w = 11, fun = mean, na.policy = "only", na.rm = TRUE), w = 9, fun = mean, na.rm = TRUE), dx = 360)
  writeRaster(bathy_smooth, smooth_tif_path, overwrite = TRUE)
}

st_geometry(islands_f1_sf) <- st_geometry(islands_f1_sf) + c(360, 0)
st_crs(islands_f1_sf) <- 4326

monument_raw <- st_read(here("data", "gis layers", "PMNM_Boundaries.shp"), quiet = TRUE) %>% st_transform(4326) %>% st_make_valid()
coords_raw <- st_coordinates(st_boundary(st_make_valid(monument_raw[toupper(monument_raw$Site_ID) == "MNM1", ])))
coords_raw[, "X"] <- ifelse(coords_raw[, "X"] < 0, coords_raw[, "X"] + 360, coords_raw[, "X"])

segment_list <- list()
for (j in 1:(nrow(coords_raw) - 1)) {
  if (coords_raw[j, "L1"] == coords_raw[j + 1, "L1"]) {
    x1 <- coords_raw[j, "X"]; y1 <- coords_raw[j, "Y"]
    x2 <- coords_raw[j + 1, "X"]; y2 <- coords_raw[j + 1, "Y"]
    if (!((abs(y1 - 25.33) < 0.25 && abs(y2 - 25.33) < 0.25 && abs(x1 - x2) > 0.05) ||
          (abs(y1 - 31.47) < 0.15 && abs(y2 - 31.47) < 0.15 && abs(x1 - x2) > 0.5))) {
      segment_list[[length(segment_list) + 1]] <- st_linestring(matrix(c(x1, x2, y1, y2), ncol = 2))
    }
  }
}
monument_sf_clean <- st_sf(geometry = st_sfc(st_line_merge(st_union(st_sfc(segment_list, crs = 4326))), crs = 4326))

label_data <- data.frame(
  lon = c(-178.3, -177.3, -175.8, -173.9, -171.7, -170.5, -167.9, -166.1, -164.7, -161.9, -159.5, -160.1, -158.0, -157.0, -156.9, -156.3, -156.6, -155.5),
  lat = c(28.4, 28.2, 27.8, 26.0, 25.7, 25.4, 25.0, 23.8, 23.5, 23.0, 22.1, 21.9, 21.4, 21.1, 20.8, 20.8, 20.5, 19.5),
  lbl_lon = c(-179.5, -177.3, -176.8, -174.5, -171.7, -171.3, -167.9, -165.0, -164.9, -161.9, -159.5, -160.7, -158.0, -156.3, -157.7, -155.5, -157.6, -156.8),
  lbl_lat = c(27.5, 29.5, 26.50, 25.0, 26.8, 23.8, 23.20, 25.0, 22.0, 24.20, 22.7, 21.5, 22.0, 21.5, 20.75, 21.1, 20.3, 19.4),
  seg_xend = c(-179.5, -177.3, -176.8, -174.4, -171.7, -170.8, -167.9, -165.0, -164.9, -161.9, rep(NA, 8)),
  seg_yend = c(28, 29.0, 27.0, 25.5, 26.30, 24.30, 23.70, 24.5, 22.5, 23.70, rep(NA, 8)),
  haw = c("H\u014dlanik\u016b", "Kuaihelani", "Manawai", "Kapou", "Kamole", "Kamokuokamohoali\u02BBi", "\u02BB\u014Cn\u016Bnui, \u02BB\u014Cn\u016Biki", "Lalo", "Mokumanamana", "Nihoa", "Kaua\u02BBi", "Ni\u02BBihau", "O\u02BBahu", "Moloka\u02BBi", "L\u0101na\u02BBi", "Maui", "Kaho\u02BBolawe", "Hawai\u02BBi"),
  eng = c("Kure Atoll", "Midway Atoll", "Pearl & Hermes Atoll", "Lisianski Island", "Laysan Island", "Maro Reef", "Gardner Pinnacles", "French Frigate Shoals", "Necker Island", rep("", 9)),
  draw_line = c(rep(TRUE, 10), rep(FALSE, 8))
) %>% mutate(
  full_label = ifelse(eng != "", paste0("<span style='font-size: 5.1pt; color:#0F172A;'><b>", haw, "</b></span><br><span style='font-size: 4.2pt; color:#475569;'><i>(", eng, ")</i></span>"), paste0("<span style='font-size: 5.1pt; color:#0F172A;'><b>", haw, "</b></span>")),
  lon = ifelse(lon < 0, lon + 360, lon),
  lbl_lon = ifelse(lbl_lon < 0, lbl_lon + 360, lbl_lon),
  seg_xend = ifelse(!is.na(seg_xend) & seg_xend < 0, seg_xend + 360, seg_xend)
)

# Straight label for NWHI
nwhi_label <- data.frame(
  region = "Northwestern Hawaiian Islands",
  lon    = 189.5,
  lat    = 27.9,
  angle  = -18
)

# Curved arc path for the Main Hawaiian Islands
mhi_curve <- tibble(
  lon = seq(200.5, 206.2, length.out = 30)
) %>%
  mutate(
    lat = 23.2 - 0.35 * (lon - 201.8) - 0.035 * (lon - 201.8)^2,
    label = "Main Hawaiian Islands"
  )

globe_inset <- ggplot() +
  geom_sf(data = ocean_circle, fill = "#FAFAFA", color = "#888888", linewidth = 0.05) +
  geom_sf(data = world_visible_ortho, fill = "#1E293B", color = NA) +
  geom_sf(data = suppressWarnings(st_intersection(suppressWarnings(st_transform(st_as_sfc(st_bbox(c(xmin = -183.0, ymin = 17.0, xmax = -152.0, ymax = 32.5), crs = 4326)), crs = ortho_crs) %>% st_make_valid()), ocean_circle)), fill = "#0F172A", color = "#020617", alpha = 0.35, linewidth = 0.20) +
  coord_sf(crs = ortho_crs) + theme_void() + theme(panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.15))

fig1_map <- ggplot() +
  geom_spatraster_contour_filled(data = bathy_smooth, breaks = c(-12000, -3000, -1500, -400, Inf)) +
  scale_fill_manual(values = c("#FFFFFF", "#F1F5F9", "#CBD5E1", "#94A3B8"), labels = c("> 3,000", "1,500 – 3,000", "400 – 1,500", "< 400"), na.value = "transparent", guide = guide_legend(order = 1, title = "Depth (m)", title.position = "top", ncol = 1, reverse = TRUE, override.aes = list(color = "black", linewidth = 0.40))) +
  geom_sf(data = monument_sf_clean, aes(color = "Papahānaumokuākea National Marine Sanctuary/ \nPapahānaumokuākea Marine National Monument"), linewidth = 0.40, linetype = "dotted") +
  scale_color_manual(name = NULL, values = "black", guide = guide_legend(order = 2, override.aes = list(linetype = "dotted", linewidth = 0.45))) +
  geom_sf(data = islands_f1_sf, fill = "#0F172A", color = "#0F172A", linewidth = 0.05) +
  geom_segment(data = subset(label_data, draw_line), aes(x = lon, y = lat, xend = seg_xend, yend = seg_yend), color = "#334155", linewidth = 0.22, alpha = 0.90, lineend = "round") +
  geom_richtext(data = label_data, aes(x = lbl_lon, y = lbl_lat, label = full_label), size = 1.9, lineheight = 0.88, color = "black", family = "sans", hjust = 0.5, fill = NA, label.color = NA) +
  geom_text(
    data     = nwhi_label,
    aes(x = lon, y = lat, label = region, angle = angle),
    color    = "#334155",
    fontface = "bold",
    size     = 2.7,
    family   = "sans"
  ) +
  geom_textpath(
    data        = mhi_curve,
    aes(x = lon, y = lat, label = label),
    text_only   = TRUE,
    color       = "#334155",
    fontface    = "bold",
    size        = 2.7,
    family      = "sans",
    vjust       = -0.2,
    upright     = TRUE
  ) +
  
  coord_sf(xlim = c(177.0, 208.0), ylim = c(17.0, 32.5), expand = FALSE, crs = "+proj=longlat +datum=WGS84 +over", clip = "on") +
  scale_x_continuous(breaks = seq(180, 208, by = 4), labels = function(x) ifelse(x == 180, "180°", ifelse(x > 180, paste0(360 - x, "°W"), paste0(x, "°E"))), sec.axis = dup_axis(name = NULL)) +
  scale_y_continuous(breaks = seq(18, 32, by = 2), labels = function(x) paste0(x, "°N"), sec.axis = dup_axis(name = NULL)) +
  theme_minimal(base_size = 6.5, base_family = "sans") +
  theme(
    plot.background      = element_rect(fill = "#FFFFFF", color = NA),
    panel.background     = element_rect(fill = "#FFFFFF", color = NA),
    panel.border         = element_rect(color = "black", fill = NA, linewidth = 0.45),
    axis.line            = element_blank(),
    axis.title           = element_blank(),
    axis.ticks           = element_line(color = "black", linewidth = 0.30),
    axis.ticks.length    = unit(0.08, "cm"),
    panel.grid.major     = element_line(color = "#EEF2F6", linewidth = 0.20),
    panel.grid.minor     = element_blank(),
    legend.position      = c(0.03, 0.03),
    legend.justification = c(0, 0),
    legend.key.height    = unit(2.7, "mm"),
    legend.key.width     = unit(4.5, "mm"),
    legend.spacing.y     = unit(0.5, "pt"),
    legend.margin        = margin(t = 2.5, r = 5, b = 2.5, l = 5, unit = "pt"),
    legend.box.margin    = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.background = element_rect(fill = alpha("#FFFFFF", 0.95), color = "black", linewidth = 0.45),
    legend.title         = element_text(size = 5.5, face = "bold", color = "black", margin = margin(b = 1)),
    legend.text          = element_text(size = 4.8, color = "black"),
    plot.margin          = margin(8, 24, 8, 28, "pt")
  )

fig1_final <- fig1_map + inset_element(globe_inset, left = 0.76, bottom = 0.68, right = 1.00, top = 0.98)

ggsave(file.path(main_dir, "Figure1_HA_map.png"), fig1_final, width = 17.8, height = 9.6, units = "cm", dpi = 600)
ggsave(file.path(main_dir, "Figure1_HA_map.tiff"), fig1_final, width = 17.8, height = 9.6, units = "cm", dpi = 600, compression = "lzw")