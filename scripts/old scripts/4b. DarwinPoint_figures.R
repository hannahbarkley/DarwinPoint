# ==============================================================================
# Script:  Master_Hawaiian_Darwin_Point_Figures_and_Tables.R
# Purpose: Generates all publication-ready main figures (1-4), supplementary 
#          figures (S1-S5), and Table 1 for the Hawaiian Darwin Point study. 
#          Compares historical (1982) carbonate budgets with modern (2024) 
#          NCRMP surveys across the archipelago's latitudinal gradient.
# Inputs:  - StRS Processed CSVs (in data/strs data/processed data/): 
#                  HA_2024_StRS_IslandLevel.csv, HA_2024_StRS_Net_SiteLevel.csv, 
#                  HA_2024_StRS_Net_SiteLevel_GriggSubset.csv
#          - StRS Field CSVs (in data/strs data/field data/):
#                  HA_2024_StRS_Prod_SiteLevel.csv
#          - Fixed CSVs (in fixed data/data/processed data/): 
#                  HA_2024_Fixed_Net_SiteLevel.csv
#          - Grigg CSV (in data/grigg data/): 
#                  HA_1982_Grigg.csv
#          - Other CSVs (in data/dbase): 
#                  prod_dbase_ha_strs.csv
#          - GIS (in data/gis layers/): 
#                  islands.shp, srtm30v5h_30contr.shp, srtm30v5h_100contr.shp,
#                  PMNM_Boundaries.shp, GEBCO bathymetry rasters (ASC/TIF)
# Outputs: - outputs/main figures/ (Figures 1-4 as PNG/TIFF, Table 1 as DOCX)
#          - outputs/supplementary figures/ (Figures S1-S5 as PNG/TIFF)
# ==============================================================================

# --- 0. ENVIRONMENT & LIBRARIES ---
library(tidyverse)
library(here)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(patchwork)
library(ggtext)
library(ggspatial)
library(scales)
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(knitr)
library(gt)

# Disable S2 spherical geometry for Pacific cropping
sf_use_s2(FALSE)
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

# Output directory setup
out_dir <- here("outputs")
main_dir <- file.path(out_dir, "main figures")
supp_dir <- file.path(out_dir, "supplementary figures")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)
if (!dir.exists(supp_dir)) dir.create(supp_dir, recursive = TRUE)

# --- 1. GLOBAL UNIFIED COLOR PALETTE ---
color_acc_red    <- "#D73027"
color_acc_orange <- "#F46D43"
color_acc_yellow <- "#FDAE61"
color_acc_lblue  <- "#74ADD1"
color_acc_blue   <- "#4575B4"
color_acc_gray   <- "#94A3B8"

pal_lat_grad <- c(
  color_acc_red,
  color_acc_orange,
  color_acc_yellow,
  "#FFFFBF",
  color_acc_lblue,
  color_acc_blue
)

# --- 2. SHARED HELPER FUNCTIONS ---

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

format_islands_with_english <- function(x, html = FALSE) {
  sep <- if (html) "<br>" else " "
  case_match(
    x,
    "H\u014dlanik\u016b"        ~ paste0("H\u014dlanik\u016b", sep, "(Kure Atoll)"),
    "Kuaihelani"               ~ paste0("Kuaihelani", sep, "(Midway Atoll)"),
    "Manawai"                  ~ paste0("Manawai", sep, "(Pearl & Hermes)"),
    "Kapou"                    ~ paste0("Kapou", sep, "(Lisianski Island)"),
    "Kamole"                   ~ paste0("Kamole", sep, "(Laysan Island)"),
    "Kamokuokamohoali'i"       ~ paste0("Kamokuokamohoali'i", sep, "(Maro Reef)"),
    "P\u016bh\u0101honu"       ~ paste0("P\u016bh\u0101honu", sep, "(Gardner Pinnacles)"),
    "Lalo"                     ~ paste0("Lalo", sep, "(French Frigate Shoals)"),
    "Mokumanamana"             ~ paste0("Mokumanamana", sep, "(Necker Island)"),
    .default = x
  )
}

theme_publication <- function(base_size = 9, y_margin_r = 2, x_margin_t = 4) {
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
      plot.margin       = margin(t = 10, r = 4, b = 2, l = 22, unit = "pt")
    )
}

no_x_title <- function(p) {
  p + labs(x = NULL) + theme(axis.title.x = element_blank())
}

boot_q90_se <- function(x, B = 1000) {
  x <- na.omit(x)
  if (length(x) <= 1) return(0)
  boot_vals <- replicate(B, quantile(sample(x, replace = TRUE), probs = 0.90, na.rm = TRUE))
  sd(boot_vals)
}

fmt_val_se <- function(val, se = NULL, digits = 2) {
  map2_chr(val, if (is.null(se)) list(NULL) else se, function(v, s) {
    if (is.na(v)) return("—")
    v_str <- sprintf(paste0("%.", digits, "f"), v)
    if (is.null(s) || is.na(s)) {
      return(v_str)
    } else {
      s_str <- sprintf(paste0("%.", digits, "f"), s)
      return(paste0(v_str, " (", s_str, ")"))
    }
  })
}

load_360_sf <- function(file_name) {
  file_path <- here("data", "gis layers", file_name)
  if (file.exists(file_path)) {
    sf_obj <- st_read(file_path, quiet = TRUE) %>% st_zm(drop = TRUE, what = "ZM") %>% st_transform(4326) %>% st_make_valid()
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

saferead_csv <- function(file_name) {
  # Scans all expected directories for datasets, including field data
  paths_to_check <- c(
    here("data", "strs data", "processed data", file_name),
    here("data", "strs data", "field data", file_name),      # Added field data path
    here("fixed data", "data", "processed data", file_name),
    here("data", "fixed data", "processed data", file_name), 
    here("data", "grigg data", file_name),   
    here("data", "dbase", file_name),
    here("data", file_name),
    here("outputs", file_name),
    here(file_name)
  )
  
  for (p in paths_to_check) {
    if (file.exists(p)) return(read_csv(p, show_col_types = FALSE))
  }
  stop(paste("Dataset", file_name, "not found in any expected directory."))
}

# --- 3. GLOBAL DATA INGESTION ---
df_strs_darwin   <- saferead_csv("HA_2024_StRS_IslandLevel.csv") %>% mutate(PLOT_LABEL = format_islands(as.character(ISLAND)))
df_strs_net      <- saferead_csv("HA_2024_StRS_Net_SiteLevel.csv") %>% mutate(LONGITUDE = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE))

# Load SW/Grigg Subset and standardize column names if using CoralNet naming (TOTAL/PROD)
df_strs_sw       <- saferead_csv("HA_2024_StRS_Net_SiteLevel_GriggSubset.csv") %>% mutate(LONGITUDE = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE))

if (!"CORAL_COVER_PERC" %in% names(df_strs_sw) && "TOTAL" %in% names(df_strs_sw)) {
  df_strs_sw$CORAL_COVER_PERC <- df_strs_sw$TOTAL
}
if (!"GROSS_PROD_KG_M2_YR" %in% names(df_strs_sw) && "PROD" %in% names(df_strs_sw)) {
  df_strs_sw$GROSS_PROD_KG_M2_YR <- df_strs_sw$PROD
}

df_fixed_net     <- saferead_csv("HA_2024_Fixed_Net_SiteLevel.csv") %>% mutate(LONGITUDE = ifelse(LONGITUDE < 0, LONGITUDE + 360, LONGITUDE))
df_grigg         <- saferead_csv("HA_1982_Grigg.csv") %>% mutate(PLOT_LABEL = format_islands(LOCATION))
df_coralnet_prod <- saferead_csv("HA_2024_StRS_Prod_SiteLevel.csv")

# NOTE: For 'prod_dbase_ha.csv' (used in S2), the saferead_csv will find it in data/ or root.
prod_dbase_s2    <- saferead_csv("prod_dbase_ha_strs.csv")

islands_sf   <- load_360_sf("islands.shp")
contour_30m  <- load_360_sf("srtm30v5h_30contr.shp")
contour_100m <- load_360_sf("srtm30v5h_100contr.shp")

# ==============================================================================
# FIGURE 1: Full Monument Map & Bathymetry
# ==============================================================================
cat("Rendering Figure 1...\n")

ortho_crs <- "+proj=ortho +lat_0=12 +lon_0=-160"
world_pacific <- ne_countries(scale = "small", returnclass = "sf") %>% st_make_valid()
vis_box_east <- st_as_sfc(st_bbox(c(xmin = 110.5, ymin = -89.0, xmax = 180.0, ymax = 89.0), crs = 4326))
vis_box_west <- st_as_sfc(st_bbox(c(xmin = -180.0, ymin = -89.0, xmax = -70.5, ymax = 89.0), crs = 4326))
world_visible_wgs84 <- rbind(suppressWarnings(st_crop(world_pacific, vis_box_east)), suppressWarnings(st_crop(world_pacific, vis_box_west))) %>% st_make_valid()
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
  } else
    bathy_cropped <- rast(tif_path)
  
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
  lbl_lon = c(-178.3, -177.3, -175.8, -174.5, -171.7, -170.3, -167.9, -166.1, -164.9, -161.9, -159.5, -160.7, -158.0, -156.3, -157.7, -155.5, -157.6, -156.8),
  lbl_lat = c(29.65, 27.1, 29.15, 25.0, 26.8, 23.8, 26.55, 22.6, 24.60, 24.20, 22.7, 21.5, 22.0, 21.5, 20.75, 21.1, 20.3, 19.4),
  seg_xend = c(-178.3, -177.3, -175.8, -174.4, -171.7, -170.3, -167.9, -166.1, -164.9, -161.9, rep(NA, 8)),
  seg_yend = c(29.15, 27.60, 28.65, 25.50, 26.30, 24.30, 26.05, 23.10, 24.10, 23.70, rep(NA, 8)),
  haw = c("H\u014dlanik\u016b", "Kuaihelani", "Manawai", "Kapou", "Kamole", "Kamokuokamohoali\u02BBi", "P\u016bh\u0101honu", "Lalo", "Mokumanamana", "Nihoa", "Kaua\u02BBi", "Ni\u02BBihau", "O\u02BBahu", "Moloka\u02BBi", "L\u0101na\u02BBi", "Maui", "Kaho\u02BBolawe", "Hawai\u02BBi"),
  eng = c("Kure Atoll", "Midway Atoll", "Pearl & Hermes Atoll", "Lisianski Island", "Laysan Island", "Maro Reef", "Gardner Pinnacles", "French Frigate Shoals", "Necker Island", rep("", 9)),
  draw_line = c(rep(TRUE, 10), rep(FALSE, 8))
) %>% mutate(
  full_label = ifelse(eng != "", paste0("<span style='font-size: 5.1pt; color:#0F172A;'><b>", haw, "</b></span><br><span style='font-size: 4.2pt; color:#475569;'><i>(", eng, ")</i></span>"), paste0("<span style='font-size: 5.1pt; color:#0F172A;'><b>", haw, "</b></span>")),
  lon = ifelse(lon < 0, lon + 360, lon),
  lbl_lon = ifelse(lbl_lon < 0, lbl_lon + 360, lbl_lon),
  seg_xend = ifelse(!is.na(seg_xend) & seg_xend < 0, seg_xend + 360, seg_xend)
)

globe_inset <- ggplot() +
  geom_sf(data = ocean_circle, fill = "#FAFAFA", color = "#888888", linewidth = 0.05) +
  geom_sf(data = world_visible_ortho, fill = "#1E293B", color = NA) +
  geom_sf(data = suppressWarnings(st_intersection(suppressWarnings(st_transform(st_as_sfc(st_bbox(c(xmin = -183.0, ymin = 17.0, xmax = -152.0, ymax = 32.5), crs = 4326)), crs = ortho_crs) %>% st_make_valid()), ocean_circle)), fill = "#0F172A", color = "#020617", alpha = 0.35, linewidth = 0.20) +
  coord_sf(crs = ortho_crs) + theme_void() + theme(panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.15))

fig1_map <- ggplot() +
  geom_spatraster_contour_filled(data = bathy_smooth, breaks = c(-12000, -3000, -1500, -400, Inf)) +
  scale_fill_manual(values = c("#FFFFFF", "#F1F5F9", "#CBD5E1", "#94A3B8"), labels = c("> 3,000", "1,500 – 3,000", "400 – 1,500", "< 400"), na.value = "transparent", guide = guide_legend(order = 1, title = "Depth (m)", title.position = "top", ncol = 1, reverse = TRUE, override.aes = list(color = "black", linewidth = 0.40))) +
  geom_sf(data = monument_sf_clean, aes(color = "Papahānaumokuākea National\nMarine Sanctuary Boundary"), linewidth = 0.40, linetype = "dotted") +
  scale_color_manual(name = NULL, values = "black", guide = guide_legend(order = 2, override.aes = list(linetype = "dotted", linewidth = 0.45))) +
  geom_sf(data = islands_f1_sf, fill = "#0F172A", color = "#0F172A", linewidth = 0.05) +
  geom_segment(data = subset(label_data, draw_line), aes(x = lon, y = lat, xend = seg_xend, yend = seg_yend), color = "#334155", linewidth = 0.22, alpha = 0.90, lineend = "round") +
  geom_richtext(data = label_data, aes(x = lbl_lon, y = lbl_lat, label = full_label), size = 1.9, lineheight = 0.88, color = "black", family = "sans", hjust = 0.5, fill = NA, label.color = NA) +
  coord_sf(xlim = c(177.0, 208.0), ylim = c(17.0, 32.5), expand = FALSE, crs = "+proj=longlat +datum=WGS84 +over", clip = "on") +
  scale_x_continuous(breaks = seq(180, 208, by = 4), labels = function(x) ifelse(x == 180, "180°", ifelse(x > 180, paste0(360 - x, "°W"), paste0(x, "°E"))), sec.axis = dup_axis(name = NULL)) +
  scale_y_continuous(breaks = seq(18, 32, by = 2), labels = function(x) paste0(x, "°N"), sec.axis = dup_axis(name = NULL)) +
  theme_minimal(base_size = 6.5, base_family = "sans") +
  theme(
    plot.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.background = element_rect(fill = "#FFFFFF", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.30),
    axis.ticks.length = unit(0.08, "cm"),
    panel.grid.major = element_line(color = "#EEF2F6", linewidth = 0.20),
    panel.grid.minor = element_blank(),
    legend.position = c(0.03, 0.03),
    legend.justification = c(0, 0),
    legend.key.height = unit(2.7, "mm"),
    legend.key.width = unit(4.5, "mm"),
    legend.spacing.y = unit(0.5, "pt"),
    legend.margin = margin(t = 2.5, r = 5, b = 2.5, l = 5, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.background = element_rect(fill = alpha("#FFFFFF", 0.95), color = "black", linewidth = 0.45),
    legend.title = element_text(size = 5.5, face = "bold", color = "black", margin = margin(b = 1)),
    legend.text = element_text(size = 4.8, color = "black"),
    plot.margin = margin(8, 24, 8, 28, "pt")
  )

fig1_final <- fig1_map + inset_element(globe_inset, left = 0.76, bottom = 0.68, right = 1.00, top = 0.98)
ggsave(file.path(main_dir, "Figure1_HA_map.png"), fig1_final, width = 17.8, height = 9.6, units = "cm", dpi = 600)
ggsave(file.path(main_dir, "Figure1_HA_map.tiff"), fig1_final, width = 17.8, height = 9.6, units = "cm", dpi = 600, compression = "lzw")


# ==============================================================================
# FIGURE 2: Historical (Grigg 1982) vs Modern (NCRMP) P90 Comparison
# ==============================================================================
cat("Rendering Figure 2...\n")

df_grigg <- df_grigg %>% mutate(VERTICAL_ACCRETION_MM_YR = GROSS_PROD / (2.89 * (1 - 0.50)))

set.seed(42)
strs_sw_p90 <- df_strs_sw %>% mutate(PLOT_LABEL = format_islands(ISLAND)) %>% group_by(PLOT_LABEL) %>%
  summarize(
    LATITUDE = mean(LATITUDE, na.rm = TRUE),
    CORAL_COVER_P90 = quantile(CORAL_COVER_PERC, probs = 0.90, na.rm = TRUE),
    CORAL_COVER_P90_SE = boot_q90_se(CORAL_COVER_PERC),
    GROSS_CARB_PROD_KG_M2_YR_P90 = quantile(GROSS_PROD_KG_M2_YR, probs = 0.90, na.rm = TRUE),
    GROSS_CARB_PROD_KG_M2_YR_P90_SE = boot_q90_se(GROSS_PROD_KG_M2_YR),
    .groups = "drop"
  )

shared_p90_deltas <- strs_sw_p90 %>% inner_join(df_grigg %>% select(PLOT_LABEL, GRIGG_LAT = LATITUDE, GRIGG_COV = CORAL_COVER_PERC, GRIGG_PROD = GROSS_PROD), by = "PLOT_LABEL") %>%
  mutate(
    LATITUDE = GRIGG_LAT,
    DELTA_COV_P90 = CORAL_COVER_P90 - GRIGG_COV,
    DELTA_COV_SE = CORAL_COVER_P90_SE,
    DELTA_PROD_P90 = GROSS_CARB_PROD_KG_M2_YR_P90 - GRIGG_PROD,
    DELTA_PROD_SE = GROSS_CARB_PROD_KG_M2_YR_P90_SE
  )

island_locs_f2 <- bind_rows(select(df_grigg, PLOT_LABEL, LATITUDE), select(strs_sw_p90, PLOT_LABEL, LATITUDE)) %>% distinct(PLOT_LABEL, .keep_all = TRUE) %>% arrange(LATITUDE)

delta_cov_data <- shared_p90_deltas %>% mutate(
  nudge_x = case_when(PLOT_LABEL == "H\u014dlanik\u016b" ~ 2.5, PLOT_LABEL == "Manawai" ~ -1.0, PLOT_LABEL == "Kapou" ~ -2.2, PLOT_LABEL == "Lalo" ~ -1.4, PLOT_LABEL == "Kaua'i" ~ 1.0, PLOT_LABEL == "Ni'ihau" ~ 1.2, PLOT_LABEL == "O'ahu" ~ -2.0, PLOT_LABEL == "Moloka'i" ~ -1.0, PLOT_LABEL == "L\u0101na'i" ~ -1.2, PLOT_LABEL == "Maui" ~ 1.6, PLOT_LABEL == "Hawai'i" ~ 1.8, TRUE ~ 1.0),
  nudge_y = case_when(PLOT_LABEL == "H\u014dlanik\u016b" ~ -2.0, PLOT_LABEL == "Manawai" ~ -12.0, PLOT_LABEL == "Kapou" ~ 8.0, PLOT_LABEL == "Lalo" ~ -12.0, PLOT_LABEL == "Kaua'i" ~ 8.0, PLOT_LABEL == "Ni'ihau" ~ -8.0, PLOT_LABEL == "O'ahu" ~ 6.0, PLOT_LABEL == "Moloka'i" ~ 8.0, PLOT_LABEL == "L\u0101na'i" ~ -8.0, PLOT_LABEL == "Maui" ~ 8.0, PLOT_LABEL == "Hawai'i" ~ 14.0, TRUE ~ 0.0)
)
delta_prod_data <- shared_p90_deltas %>% mutate(
  nudge_x = case_when(PLOT_LABEL == "H\u014dlanik\u016b" ~ 1.2, PLOT_LABEL == "Manawai" ~ -1.0, PLOT_LABEL == "Kapou" ~ -2.2, PLOT_LABEL == "Lalo" ~ -1.6, PLOT_LABEL == "Kaua'i" ~ 1.0, PLOT_LABEL == "Ni'ihau" ~ 1.2, PLOT_LABEL == "O'ahu" ~ -1.2, PLOT_LABEL == "Moloka'i" ~ -1.0, PLOT_LABEL == "L\u0101na'i" ~ -1.2, PLOT_LABEL == "Maui" ~ 1.6, PLOT_LABEL == "Hawai'i" ~ 1.6, TRUE ~ 1.0),
  nudge_y = case_when(PLOT_LABEL == "H\u014dlanik\u016b" ~ -2.5, PLOT_LABEL == "Manawai" ~ 2.5, PLOT_LABEL == "Kapou" ~ 1.6, PLOT_LABEL == "Lalo" ~ -1.3, PLOT_LABEL == "Kaua'i" ~ 2.0, PLOT_LABEL == "Ni'ihau" ~ -2.5, PLOT_LABEL == "O'ahu" ~ -2.8, PLOT_LABEL == "Moloka'i" ~ 2.5, PLOT_LABEL == "L\u0101na'i" ~ -2.5, PLOT_LABEL == "Maui" ~ 2.4, PLOT_LABEL == "Hawai'i" ~ 2.5, TRUE ~ 0.0)
)

build_plot_f2 <- function(data, y_var, y_se = NULL, y_label, y_lim, y_brk, add_zero = FALSE, add_lbl = FALSE) {
  p <- ggplot(data, aes(x = LATITUDE, y = .data[[y_var]])) + scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) + theme_publication()
  if (!is.null(y_se))
    p <- p + geom_errorbar(aes(ymin = .data[[y_var]] - .data[[y_se]], ymax = .data[[y_var]] + .data[[y_se]]), color = "gray55", width = 0, linewidth = 0.5)
  p <- p + scale_y_continuous(breaks = y_brk) + coord_cartesian(ylim = y_lim, clip = "off")
  if (add_zero)
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60", linewidth = 0.6)
  
  p <- p + geom_point(aes(fill = LATITUDE), shape = 21, color = "gray30", stroke = 0.5, size = 2.8, alpha = 0.85) +
    scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
    labs(x = "Latitude (°N)", y = y_label)
  
  if (add_lbl) {
    nx <- if ("nudge_x" %in% names(data)) data$nudge_x else 0
    ny <- if ("nudge_y" %in% names(data)) data$nudge_y else 0
    p <- p + geom_text_repel(aes(label = PLOT_LABEL), size = 2.4, fontface = "bold", color = "gray20", box.padding = 0.25, point.padding = 0.35, nudge_x = nx, nudge_y = ny, force = 1, seed = 42, segment.color = "gray60", segment.alpha = 0.8, segment.linewidth = 0.35, min.segment.length = 0)
  }
  return(p)
}

lbl_cov <- "Coral cover (%)"
lbl_prod <- expression(bold("Gross production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
grigg_cover <- build_plot_f2(df_grigg, "CORAL_COVER_PERC", y_label = lbl_cov, y_lim = c(0, 100), y_brk = seq(0, 100, 20))
ncrmp_cover <- build_plot_f2(strs_sw_p90, "CORAL_COVER_P90", "CORAL_COVER_P90_SE", y_label = lbl_cov, y_lim = c(0, 100), y_brk = seq(0, 100, 20))
grigg_prod  <- build_plot_f2(df_grigg, "GROSS_PROD", y_label = lbl_prod, y_lim = c(0, 20), y_brk = seq(0, 20, 5))
ncrmp_prod  <- build_plot_f2(strs_sw_p90, "GROSS_CARB_PROD_KG_M2_YR_P90", "GROSS_CARB_PROD_KG_M2_YR_P90_SE", y_label = lbl_prod, y_lim = c(0, 20), y_brk = seq(0, 20, 5))
delta_cover <- build_plot_f2(delta_cov_data, "DELTA_COV_P90", "DELTA_COV_SE", y_label = "\u0394 Coral cover (%)", y_lim = c(-80, 40), y_brk = seq(-80, 40, 20), TRUE, TRUE)
delta_prod  <- build_plot_f2(delta_prod_data, "DELTA_PROD_P90", "DELTA_PROD_SE", y_label = expression(bold(Delta ~ "Gross production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), y_lim = c(-15, 10), y_brk = seq(-15, 10, 5), TRUE, TRUE)

main_plots_f2 <- ggarrange(no_x_title(grigg_cover), no_x_title(ncrmp_cover), no_x_title(delta_cover), grigg_prod, ncrmp_prod, delta_prod, ncol = 3, nrow = 2, align = "hv", labels = "auto", font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02)
p_legend_f2 <- ggplot() + geom_tile(data = data.frame(lat = seq(18, 30, len = 500), x = 0), aes(x = x, y = lat, fill = lat), width = 0.5) + annotate("rect", xmin = -0.25, xmax = 0.25, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) + scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") + scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.02, 0.02)) + scale_x_continuous(limits = c(-0.25, 3.8), expand = c(0, 0)) + geom_point(data = island_locs_f2, aes(x = 0.25, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2, color = "gray30", stroke = 0.5) + geom_text_repel(data = island_locs_f2, aes(x = 0.25, y = LATITUDE, label = PLOT_LABEL), direction = "y", hjust = 0, nudge_x = 0.4, segment.color = "gray50", segment.linewidth = 0.4, size = 2.6, fontface = "bold", color = "gray20", force = 2) + theme_void() + theme(plot.title = element_text(size = 8, face = "bold", margin = margin(b = 4)), axis.text.y = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 5)), axis.ticks.y = element_line(color = "gray20"), axis.ticks.length = unit(0.15, "cm"), plot.margin = margin(20, 0, 20, 10)) + labs(title = "Latitude (°N)")

fig2_final <- ggarrange(main_plots_f2, p_legend_f2, ncol = 2, widths = c(1, 0.22))
ggsave(file.path(main_dir, "Figure2_Grigg_vs_NCRMP.png"), fig2_final, width = 28, height = 12.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(main_dir, "Figure2_Grigg_vs_NCRMP.tiff"), fig2_final, width = 28, height = 12.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")


# ==============================================================================
# FIGURE 3: Full StRS Metrics across Latitude
# ==============================================================================
cat("Rendering Figure 3...\n")

island_locs_f3 <- df_strs_darwin %>% select(PLOT_LABEL, LATITUDE) %>% distinct(PLOT_LABEL, .keep_all = TRUE) %>% arrange(LATITUDE)

p_strs_cov   <- build_plot_f2(df_strs_darwin, "CORAL_COVER_PERC", "se.CORAL_COVER_PERC", y_label = lbl_cov, y_lim = c(0, 30), y_brk = seq(0, 30, 5))
p_strs_gprod <- build_plot_f2(df_strs_darwin, "GROSS_PROD_KG_M2_YR", "se.GROSS_PROD_KG_M2_YR", y_label = lbl_prod, y_lim = c(0, 6), y_brk = seq(0, 6, 1), FALSE)
p_strs_geros <- build_plot_f2(df_strs_darwin, "GROSS_EROSION_KG_M2_YR", "se.GROSS_EROSION_KG_M2_YR", y_label = expression(bold("Gross erosion rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), y_lim = c(0, 4), y_brk = seq(0, 4, 1))
p_strs_nprod <- build_plot_f2(df_strs_darwin, "NET_PROD_KG_M2_YR", "se.NET_PROD_KG_M2_YR", y_label = expression(bold("Net production rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), y_lim = c(-4.5, 4.2), y_brk = seq(-4, 4, by = 2), TRUE)

main_plots_f3 <- ggarrange(no_x_title(p_strs_cov), no_x_title(p_strs_gprod), p_strs_geros, p_strs_nprod, ncol = 2, nrow = 2, align = "hv", labels = "auto", font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02)
p_legend_f3 <- ggplot() + geom_tile(data = data.frame(lat = seq(18, 30, len = 500), x = 0), aes(x = x, y = lat, fill = lat), width = 0.5) + annotate("rect", xmin = -0.25, xmax = 0.25, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) + scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") + scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.02, 0.02)) + scale_x_continuous(limits = c(-0.25, 3.8), expand = c(0, 0)) + geom_point(data = island_locs_f3, aes(x = 0.25, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2, color = "gray30", stroke = 0.5) + geom_text_repel(data = island_locs_f3, aes(x = 0.25, y = LATITUDE, label = PLOT_LABEL), direction = "y", hjust = 0, nudge_x = 0.4, segment.color = "gray50", segment.linewidth = 0.4, size = 2.6, fontface = "bold", color = "gray20", force = 3, box.padding = unit(0.3, "lines"), min.segment.length = 0) + coord_cartesian(clip = "off") + theme_void() + theme(plot.title = element_text(size = 8, face = "bold", margin = margin(b = 4)), axis.text.y = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 5)), axis.ticks.y = element_line(color = "gray20"), axis.ticks.length = unit(0.15, "cm"), plot.margin = margin(15, 0, 15, 10)) + labs(title = "Latitude (°N)")

fig3_final <- ggarrange(main_plots_f3, p_legend_f3, ncol = 2, widths = c(1, 0.22))
ggsave(file.path(main_dir, "Figure3_FullStRS.png"), fig3_final, width = 24.0, height = 14.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(main_dir, "Figure3_FullStRS.tiff"), fig3_final, width = 24.0, height = 14.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")


# ==============================================================================
# FIGURE 4: VAR vs SLR
# ==============================================================================
cat("Rendering Figure 4...\n")

se_col_f4 <- if ("se.VERTICAL_ACCRETION_MM_YR" %in% names(df_strs_darwin)) {
  "se.VERTICAL_ACCRETION_MM_YR"
} else {
  "VERTICAL_ACCRETION_MM_YR_SE"
}
df_strs_f4 <- df_strs_darwin %>% mutate(
  nudge_x = case_when(PLOT_LABEL == "Hawai'i" ~ 0.8, PLOT_LABEL == "Kaho'olawe" ~ 0.5, PLOT_LABEL == "Maui" ~ 1.0, PLOT_LABEL == "L\u0101na'i" ~ -0.2, PLOT_LABEL == "Moloka'i" ~ -1.0, PLOT_LABEL == "O'ahu" ~ 0.3, PLOT_LABEL == "Kaua'i" ~ -0.8, PLOT_LABEL == "Ni'ihau" ~ -0.6, PLOT_LABEL == "Lalo" ~ -0.8, PLOT_LABEL == "Kapou" ~ -0.5, PLOT_LABEL == "Manawai" ~ 0.2, PLOT_LABEL == "H\u014dlanik\u016b" ~ 0.8, TRUE ~ 0),
  nudge_y = case_when(PLOT_LABEL == "Hawai'i" ~ 0.5, PLOT_LABEL == "Kaho'olawe" ~ 2.5, PLOT_LABEL == "Maui" ~ -0.2, PLOT_LABEL == "L\u0101na'i" ~ 1.5, PLOT_LABEL == "Moloka'i" ~ 0.8, PLOT_LABEL == "O'ahu" ~ -1.5, PLOT_LABEL == "Kaua'i" ~ 0.8, PLOT_LABEL == "Ni'ihau" ~ -1.2, PLOT_LABEL == "Lalo" ~ 1.0, PLOT_LABEL == "Kapou" ~ 1.5, PLOT_LABEL == "Manawai" ~ 1.5, PLOT_LABEL == "H\u014dlanik\u016b" ~ -1.5, TRUE ~ 0)
)

fig4_final <- ggplot(df_strs_f4, aes(x = LATITUDE, y = VERTICAL_ACCRETION_MM_YR)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 6, ymax = 10, fill = "gray80", alpha = 0.15) +
  geom_hline(yintercept = c(6, 10), linetype = "dotted", color = "gray80", linewidth = 0.45) +
  annotate("text", x = 29.8, y = 8, label = "Predicted sea level rise by 2050", color = "gray20", size = 2.8, fontface = "bold", hjust = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.55) +
  geom_errorbar(aes(ymin = VERTICAL_ACCRETION_MM_YR - .data[[se_col_f4]], ymax = VERTICAL_ACCRETION_MM_YR + .data[[se_col_f4]]), color = "gray55", width = 0, linewidth = 0.45) +
  geom_point(aes(fill = LATITUDE), shape = 21, color = "gray20", stroke = 0.5, size = 3.4, alpha = 0.95) +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  geom_text_repel(aes(label = PLOT_LABEL), nudge_x = df_strs_f4$nudge_x, nudge_y = df_strs_f4$nudge_y, size = 2.8, fontface = "bold", color = "black", box.padding = 0, point.padding = 0.2, force = 0.5, max.overlaps = Inf, seed = 42, segment.color = "gray40", segment.alpha = 0.75, segment.linewidth = 0.35, min.segment.length = 0) +
  scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
  scale_y_continuous(breaks = seq(-4, 12, by = 2)) +
  coord_cartesian(ylim = c(-4.5, 12.5), clip = "off") +
  theme_publication() +
  theme(plot.margin = margin(15, 20, 10, 15, "pt"), axis.title.y = element_text(margin = margin(r = 4)), axis.title.x = element_text(margin = margin(t = 6))) +
  labs(x = "Latitude (°N)", y = expression(bold("Vertical accretion rate (mm" ~ yr^{ -1 } * ")")))

ggsave(file.path(main_dir, "Figure4_VARvsSLR.png"), fig4_final, width = 16.0, height = 11.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(main_dir, "Figure4_VARvsSLR.tiff"), fig4_final, width = 16.0, height = 11.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")


# ==============================================================================
# FIGURE S1: Insets & Main Hawaiian Islands Map
# ==============================================================================
cat("Rendering Figure S1...\n")

mhi_bbox <- st_bbox(c(xmin = 198.0, ymin = 18.5, xmax = 206.0, ymax = 23.0), crs = st_crs(islands_sf))
mhi_islands <- st_crop(islands_sf, mhi_bbox)

mhi_bathy_sf <- NULL
if (exists("bathy_smooth")) {
  mhi_ext_360 <- ext(198.0, 206.0, 18.0, 23.0)
  bathy_mhi <- crop(bathy_smooth, mhi_ext_360)
  mhi_islands_v <- vect(mhi_islands)
  mhi_immediate_mask <- buffer(mhi_islands_v, width = 0.50)
  bathy_mhi_smooth <- mask(bathy_mhi, mhi_immediate_mask)
  bathy_mhi_smooth[bathy_mhi_smooth > 0 & !is.na(bathy_mhi_smooth)] <- 0
  
  rcl <- matrix(c(-12000, -3000, 1, -3000, -1500, 2, -1500, -400, 3, -400, 0, 4), ncol = 3, byrow = TRUE)
  bathy_poly_v <- as.polygons(classify(bathy_mhi_smooth, rcl, include.lowest = TRUE))
  if (length(bathy_poly_v) > 0) {
    mhi_bathy_sf <- st_as_sf(bathy_poly_v)
    attr_col <- setdiff(names(mhi_bathy_sf), attr(mhi_bathy_sf, "sf_column"))[1]
    if (!is.na(attr_col)) {
      names(mhi_bathy_sf)[names(mhi_bathy_sf) == attr_col] <- "depth_bin"
      mhi_bathy_sf$depth_bin <- as.character(mhi_bathy_sf$depth_bin)
    }
  }
}

light_theme_base <- theme_minimal(base_size = 8, base_family = "sans") + theme(
  plot.background = element_rect(fill = "#FFFFFF", color = NA),
  panel.background = element_rect(fill = "#FFFFFF", color = NA),
  panel.border = element_rect(color = "#0F172A", fill = NA, linewidth = 0.6),
  panel.grid = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "none",
  plot.margin = margin(2, 2, 2, 2, "pt")
)

create_inset <- function(lon_c, lat_c, span = 0.20, l_30 = NULL, l_100 = NULL, tag = "A") {
  p <- ggplot()
  if (!is.null(l_100))
    p <- p + geom_sf(data = l_100, fill = "#CBD5E1", color = "#CBD5E1", linewidth = 0.1)
  if (!is.null(l_30))
    p <- p + geom_sf(data = l_30, fill = "#94A3B8", color = "#94A3B8", linewidth = 0.1)
  
  # AESTHETICS: Bumped point sizes to 1.4/1.3 for better visibility 
  p + geom_sf(data = islands_sf, fill = "#1E293B", color = "#0F172A", linewidth = 0.2) +
    geom_point(data = df_strs_net, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_acc_blue, color = "#FFFFFF", stroke = 0.25, size = 1.4, alpha = 0.85) +
    geom_point(data = df_strs_sw, aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_acc_orange, color = "#FFFFFF", stroke = 0.25, size = 1.4, alpha = 0.90) +
    geom_point(data = df_fixed_net, aes(x = LONGITUDE, y = LATITUDE), shape = 24, fill = color_acc_yellow, color = "#0F172A", stroke = 0.30, size = 1.3, alpha = 0.95) +
    annotate("text", x = lon_c - span + (2 * span * 0.05), y = lat_c + span - (2 * span * 0.05), label = tag, fontface = "bold", size = 4.2, color = "#0F172A", hjust = 0, vjust = 1) +
    annotation_scale(location = "bl", width_hint = 0.35, style = "bar", bar_cols = c("#0F172A", "white"), text_col = "#0F172A", text_cex = 0.58, height = unit(0.04, "cm")) +
    coord_sf(xlim = c(lon_c - span, lon_c + span), ylim = c(lat_c - span, lat_c + span), expand = FALSE, datum = NA) + light_theme_base
}

inset_kure    <- create_inset(181.66, 28.42, 0.11, contour_30m, contour_100m, "A")
inset_manawai <- create_inset(184.13, 27.84, 0.20, contour_30m, contour_100m, "B")
inset_lisian  <- create_inset(186.05, 26.05, 0.23, contour_30m, contour_100m, "C")
inset_ffs     <- create_inset(193.79, 23.75, 0.20, contour_30m, contour_100m, "D")

main_map_s1 <- ggplot()
if (!is.null(mhi_bathy_sf))
  main_map_s1 <- main_map_s1 + geom_sf(data = mhi_bathy_sf, aes(fill = depth_bin), color = NA) +
  scale_fill_manual(values = c("1" = "#FFFFFF", "2" = "#F1F5F9", "3" = "#CBD5E1", "4" = "#94A3B8"), na.value = "transparent")

# Main Map AESTHETICS: Bumped point sizes to 1.4/1.3 for better visibility 
main_map_s1 <- main_map_s1 + 
  geom_sf(data = mhi_islands, fill = "#1E293B", color = "#0F172A", linewidth = 0.2) +
  geom_point(data = subset(df_strs_net, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_acc_blue, color = "#FFFFFF", stroke = 0.25, size = 1.4, alpha = 0.85) +
  geom_point(data = subset(df_strs_sw, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), shape = 21, fill = color_acc_orange, color = "#FFFFFF", stroke = 0.25, size = 1.4, alpha = 0.90) +
  geom_point(data = subset(df_fixed_net, LONGITUDE > 198.0), aes(x = LONGITUDE, y = LATITUDE), shape = 24, fill = color_acc_yellow, color = "#0F172A", stroke = 0.30, size = 1.3, alpha = 0.95) +
  
  geom_text(data = data.frame(
    lbl = c("Kaua'i", "Ni'ihau", "O'ahu", "Moloka'i", "L\u0101na'i", "Maui", "Kaho'olawe", "Hawai'i"), 
    lbl_lon = c(200.15, 199.40, 202.05, 203.10, 202.70, 204.05, 202.90, 203.60), 
    lbl_lat = c(22.35,  21.75,  21.88,  21.38,  20.75,  21.05,  20.40,  19.25)
  ), aes(x = lbl_lon, y = lbl_lat, label = lbl), size = 3.6, color = "#0F172A", fontface = "bold") +
  
  annotate("text", x = 198.9 + (6.6 * 0.02), y = 22.6 - (3.8 * 0.03), label = "E", fontface = "bold", size = 4.2, color = "#0F172A", hjust = 0, vjust = 1) +
  annotation_scale(location = "bl", width_hint = 0.18, style = "bar", bar_cols = c("#0F172A", "white"), text_cex = 0.70, height = unit(0.05, "cm")) +
  coord_sf(xlim = c(198.9, 205.5), ylim = c(18.8, 22.6), expand = FALSE, datum = NA) + light_theme_base

figS1_final <- inset_kure + inset_manawai + inset_lisian + inset_ffs + main_map_s1 +
  plot_layout(design = c(patchwork::area(1, 1, 1, 1), patchwork::area(1, 2, 1, 2), patchwork::area(1, 3, 1, 3), patchwork::area(1, 4, 1, 4), patchwork::area(2, 1, 2, 4)), heights = c(1, 2.5))

ggsave(file.path(supp_dir, "FigureS1_Sites.png"), figS1_final, width = 17.8, height = 15.0, units = "cm", dpi = 600, bg = "white")
ggsave(file.path(supp_dir, "FigureS1_Sites.tiff"), figS1_final, width = 17.8, height = 15.0, units = "cm", dpi = 600, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE S2: Coral Growth Rates & Estimated Production
# ==============================================================================
cat("Rendering Figure S2 (Growth Rates & Production)...\n")

# Load pre-computed database with coefficients and raw dbase with latitudes
coeffs_dbase_s2 <- saferead_csv("prod_dbase_ha_strs_coeffs.csv")
raw_dbase_s2    <- prod_dbase_s2 # Loaded during global data ingestion

# Target substrate-morphology codes sorted ALPHABETICALLY BY GENUS NAME:
target_codes <- c("MOSP-EM", "POCS-BR", "POSP-BR", "POSP-EM", "POSP-MD")

# Native Plotmath Labels
group_labels_s2 <- c(
  expression("Encrusting "*italic("Montipora")), 
  expression("Branching "*italic("Pocillopora")), 
  expression("Branching "*italic("Porites")), 
  expression("Encrusting "*italic("Porites")), 
  expression("Massive "*italic("Porites"))
)

# Empirical median colony lengths (cm) from HA_2025_StRS_demography.csv:
#   MOSP-EM: 8.56 cm
#   POCS-BR: 11.90 cm
#   POSP-BR: 16.64 cm
#   POSP-EM:  5.66 cm
#   POSP-MD: 11.23 cm
median_lengths <- c(
  "MOSP-EM" = 8.56,
  "POCS-BR" = 11.90,
  "POSP-BR" = 16.64,
  "POSP-EM" = 5.66,
  "POSP-MD" = 11.23
)

# Join LATITUDE from raw_dbase_s2 back onto coeffs_dbase_s2
growth_df_s2 <- coeffs_dbase_s2 %>%
  mutate(
    SUBSTRATE_CODE = case_when(
      SUBSTRATE_CODE %in% c("MOEN", "MONTI", "MCAP") ~ "MOSP",
      TRUE ~ SUBSTRATE_CODE
    ),
    MORPHOLOGYCODE = case_when(
      MORPHOLOGYCODE %in% c("ENC", "EN") ~ "EM",
      TRUE ~ MORPHOLOGYCODE
    ),
    CORAL_GROUP = paste0(SUBSTRATE_CODE, "-", MORPHOLOGYCODE)
  ) %>%
  filter(CORAL_GROUP %in% target_codes) %>%
  left_join(
    raw_dbase_s2 %>% 
      mutate(
        SUBSTRATE_CODE = case_when(
          SUBSTRATE_CODE %in% c("MOEN", "MONTI", "MCAP") ~ "MOSP",
          TRUE ~ SUBSTRATE_CODE
        ),
        MORPHOLOGYCODE = case_when(
          MORPHOLOGYCODE %in% c("ENC", "EN") ~ "EM",
          TRUE ~ MORPHOLOGYCODE
        )
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
    # Taxon-specific median colony length (cm) from demography data
    x = median_lengths[as.character(CORAL_GROUP)],
    # Annual colony production in g/yr for median-sized colony
    Colony_Prod = Coef_Mean * x + Intercept_Mean
  )

# Convert to ordered factor matching alphabetical genus sequence
growth_df_s2$CORAL_GROUP <- factor(growth_df_s2$CORAL_GROUP, levels = target_codes)

# Publication palette
group_pal_s2 <- c(
  "#003087", # Encrusting Montipora
  "#E41A1C", # Branching Pocillopora
  "#1EBEC7", # Branching Porites
  "#0085CA", # Encrusting Porites
  "#FF8400"  # Massive Porites
)

# Apply vertical visual offsets
growth_df_s2 <- growth_df_s2 %>%
  mutate(
    Ext_plot  = Extension_cm_yr + c("MOSP-EM" = -0.02, "POSP-EM" = 0.0, "POCS-BR" = 0.00, "POSP-BR" = 0.0, "POSP-MD" = 0.02)[as.character(CORAL_GROUP)],
    Den_plot  = Density_g_cm3   + c("MOSP-EM" = -0.02, "POSP-EM" = 0.0, "POCS-BR" = 0.02, "POSP-BR" = 0.0, "POSP-MD" = 0.02)[as.character(CORAL_GROUP)],
    Prod_plot = Colony_Prod
  )

# Panel A: Linear Extension Rate vs. Latitude
p_ext_s2 <- ggplot(growth_df_s2, aes(x = LATITUDE, y = Ext_plot)) +
  geom_smooth(aes(color = CORAL_GROUP), method = "loess", se = FALSE, span = 0.8, linewidth = 0.8) +
  geom_point(aes(fill = CORAL_GROUP), shape = 21, color = "gray30", stroke = 0.5, size = 2.8, alpha = 0.85) +
  scale_color_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_fill_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
  scale_y_continuous(limits = c(0, 2.4), breaks = seq(0, 2.4, 0.4)) +
  coord_cartesian(ylim = c(0, 2.4), clip = "off") +
  labs(x = "Latitude (°N)", y = expression(bold("Linear extension rate (cm" ~ yr^{ -1 } * ")"))) +
  theme_publication() +
  theme(legend.position = "none")

# Panel B: Skeletal Bulk Density vs. Latitude
p_den_s2 <- ggplot(growth_df_s2, aes(x = LATITUDE, y = Den_plot)) +
  geom_smooth(aes(color = CORAL_GROUP), method = "loess", se = FALSE, span = 0.8, linewidth = 0.8) +
  geom_point(aes(fill = CORAL_GROUP), shape = 21, color = "gray30", stroke = 0.5, size = 2.8, alpha = 0.85) +
  scale_color_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_fill_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
  scale_y_continuous(limits = c(1.2, 1.8), breaks = seq(1.2, 1.8, 0.1)) +
  coord_cartesian(ylim = c(1.2, 1.8), clip = "off") +
  labs(x = "Latitude (°N)", y = expression(bold("Skeletal density (g" ~ cm^{ -3 } * ")"))) +
  theme_publication() +
  theme(legend.position = "none")

# Panel C: Estimated Median Colony Production Rate vs. Latitude
p_prod_s2 <- ggplot(growth_df_s2, aes(x = LATITUDE, y = Prod_plot)) +
  geom_smooth(aes(color = CORAL_GROUP), method = "loess", se = FALSE, span = 0.8, linewidth = 0.8) +
  geom_point(aes(fill = CORAL_GROUP), shape = 21, color = "gray30", stroke = 0.5, size = 2.8, alpha = 0.85) +
  scale_color_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_fill_manual(values = group_pal_s2, labels = group_labels_s2) +
  scale_x_reverse(limits = c(30.5, 17.5), breaks = seq(18, 30, by = 2)) +
  scale_y_continuous(limits = c(0, 25), breaks = seq(0,25, 5)) +
  coord_cartesian(ylim = c(0, 25), clip = "off") +
  labs(x = "Latitude (°N)", y = expression(bold("Colony production rate (g" ~ yr^{ -1 } * ")"))) +
  theme_publication() +
  theme(legend.position = "none")

# Single Shared Taxon Legend
p_legend_s2_growth <- get_legend(
  ggplot(growth_df_s2, aes(x = LATITUDE, y = Den_plot, fill = CORAL_GROUP)) +
    geom_point(shape = 21, color = "gray30", stroke = 0.5, size = 3) +
    scale_fill_manual(values = group_pal_s2, labels = group_labels_s2, name = NULL) +
    guides(fill = guide_legend(ncol = 1, override.aes = list(shape = 21, size = 3, stroke = 0.5, color = "gray30"))) +
    theme_publication() +
    theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 8.5, color = "gray20"))
)

# Arrange 3 stacked panels
main_plots_s2 <- ggarrange(
  no_x_title(p_ext_s2),
  no_x_title(p_den_s2),
  p_prod_s2,
  ncol = 1,
  nrow = 3,
  align = "v",
  labels = c("a", "b", "c"),
  font.label = list(size = 12, face = "bold"),
  label.x = 0.01,
  label.y = 1.02
)

# Final Assembly & Export
figS2_final <- ggarrange(
  main_plots_s2,
  p_legend_s2_growth,
  ncol = 2,
  nrow = 1,
  widths = c(1, 0.3)
)

ggsave(file.path(supp_dir, "FigureS2_growthrates.png"), figS2_final, width = 17.5, height = 21.0, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS2_growthrates.tiff"), figS2_final, width = 17.5, height = 21.0, units = "cm", dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# FIGURE S3: SfM vs CoralNet
# ==============================================================================
cat("Rendering Figure S3...\n")

matched_s3 <- df_fixed_net %>% select(OCC_SITEID, LOCATION, LATITUDE, LONGITUDE, SFM_PROD = GROSS_CARB_PROD_KG_M2_YR_MEAN, SFM_SE = GROSS_CARB_PROD_KG_M2_YR_SE) %>%
  inner_join(
    df_coralnet_prod %>% select(SITE, CORALNET_PROD = PROD, CORALNET_LWR = PROD_lwr, CORALNET_UPR = PROD_upr),
    by = c("OCC_SITEID" = "SITE")
  ) %>%
  mutate(
    ISLAND_LABEL = format_islands(LOCATION),
    CORALNET_SE = (CORALNET_UPR - CORALNET_LWR) / (2 * 1.96)
  )

island_locs_s3 <- matched_s3 %>% group_by(ISLAND_LABEL) %>% summarize(LATITUDE = mean(LATITUDE), .groups = "drop") %>% rename(PLOT_LABEL = ISLAND_LABEL) %>% arrange(LATITUDE)
max_sfm_val <- max(matched_s3$SFM_PROD)

p_main_s3 <- ggplot(matched_s3, aes(x = SFM_PROD, y = CORALNET_PROD)) +
  geom_segment(x = 0, y = 0, xend = max_sfm_val, yend = max_sfm_val, linetype = "dashed", color = "gray60", linewidth = 0.55) +
  geom_smooth(method = "lm", formula = y ~ x, color = "gray15", linewidth = 0.75, se = TRUE, fill = "gray80", alpha = 0.35) +
  geom_errorbar(aes(ymin = CORALNET_PROD - CORALNET_SE, ymax = CORALNET_PROD + CORALNET_SE), color = "gray70", width = 0, linewidth = 0.35) +
  geom_errorbarh(aes(xmin = SFM_PROD - SFM_SE, xmax = SFM_PROD + SFM_SE), color = "gray70", height = 0, linewidth = 0.35) +
  geom_point(aes(fill = LATITUDE), shape = 21, color = "gray20", stroke = 0.55, size = 3.3, alpha = 0.90) +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 2), expand = c(0.01, 0.01)) +
  coord_fixed(ratio = 1, xlim = c(0, 14), ylim = c(0, 14), clip = "off") +
  theme_publication(base_size = 9.5, y_margin_r = 6, x_margin_t = 6) +
  theme(plot.margin = margin(10, 4, 10, 10, "pt")) +
  labs(x = expression(bold("SfM-based Gross Production (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")),
       y = expression(bold("CoralNet-based Gross Production (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")))

p_legend_s3 <- ggplot() +
  geom_tile(data = data.frame(lat = seq(18, 30, len = 500), x = 0), aes(x = x, y = lat, fill = lat), width = 0.45) +
  annotate("rect", xmin = -0.225, xmax = 0.225, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.01, 0.01)) +
  scale_x_continuous(limits = c(-0.25, 3.2), expand = c(0, 0)) +
  geom_point(data = island_locs_s3, aes(x = 0.225, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2.1, color = "gray30", stroke = 0.5) +
  geom_text_repel(data = island_locs_s3, aes(x = 0.225, y = LATITUDE, label = PLOT_LABEL), direction = "y", hjust = 0, nudge_x = 0.35, segment.color = "gray50", segment.linewidth = 0.4, size = 2.7, fontface = "bold", color = "gray20", force = 2, min.segment.length = 0) +
  theme_void() +
  theme(
    plot.title = element_text(size = 8.5, face = "bold", margin = margin(b = 4)),
    axis.text.y = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 4)),
    axis.ticks.y = element_line(color = "gray20"),
    axis.ticks.length = unit(0.12, "cm"),
    plot.margin = margin(12, 0, 12, 0)
  ) +
  labs(title = "Latitude (°N)")

figS3_final <- ggarrange(p_main_s3, p_legend_s3, ncol = 2, widths = c(1, 0.28), align = "h")
ggsave(file.path(supp_dir, "FigureS3_SfM_vs_CoralNet.png"), figS3_final, width = 17.5, height = 13.5, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS3_SfM_vs_CoralNet.tiff"), figS3_final, width = 17.5, height = 13.5, units = "cm", dpi = 300, compression = "lzw", bg = "white")


# ==============================================================================
# FIGURE S4: Latitudinal Sensitivity & P90 Comparisons (Combined 2x2)
# ==============================================================================
cat("Rendering Figure S4 (Sensitivity)...\n")

prep_sens <- function(sw_data, full_data, var_name) {
  full_domain <- full_data %>%
    mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
    select(PLOT_LABEL, `All Sites` = .data[[var_name]])
  
  sw_data %>%
    mutate(PLOT_LABEL = format_islands(ISLAND)) %>%
    group_by(PLOT_LABEL) %>%
    summarize(
      LATITUDE  = mean(LATITUDE, na.rm = TRUE),
      `SW Mean` = mean(.data[[var_name]], na.rm = TRUE),
      `SW 75th` = quantile(.data[[var_name]], probs = 0.75, na.rm = TRUE),
      `SW 90th` = quantile(.data[[var_name]], probs = 0.90, na.rm = TRUE),
      `SW Max`  = max(.data[[var_name]], na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    left_join(full_domain, by = "PLOT_LABEL") %>%
    arrange(desc(LATITUDE)) %>%
    mutate(PLOT_LABEL = factor(PLOT_LABEL, levels = PLOT_LABEL))
}

coral_stats <- prep_sens(df_strs_sw, df_strs_darwin, "CORAL_COVER_PERC")
prod_stats  <- prep_sens(df_strs_sw, df_strs_darwin, "GROSS_PROD_KG_M2_YR")

build_sens_plot <- function(df_wide, y_lab, y_lim, y_brk) {
  df_long <- df_wide %>%
    pivot_longer(
      cols = c("All Sites", "SW Mean", "SW 75th", "SW 90th", "SW Max"),
      names_to = "Statistic",
      values_to = "Value"
    ) %>%
    mutate(Statistic = factor(Statistic, levels = c("All Sites", "SW Mean", "SW 75th", "SW 90th", "SW Max")))
  
  ggplot() +
    geom_segment(data = df_wide, aes(x = PLOT_LABEL, xend = PLOT_LABEL, y = pmin(`All Sites`, `SW Mean`, na.rm = TRUE), yend = `SW Max`), color = "gray60", linetype = "dashed", alpha = 0.7, linewidth = 0.6) +
    geom_point(data = df_long, aes(x = PLOT_LABEL, y = Value, fill = Statistic), shape = 21, color = "black", size = 2.8, stroke = 0.5, alpha = 0.85) +
    scale_y_continuous(limits = y_lim, breaks = y_brk) +
    scale_fill_manual(values = c("All Sites" = "#1E293B", "SW Mean" = color_acc_blue, "SW 75th" = color_acc_lblue, "SW 90th" = color_acc_orange, "SW Max" = color_acc_red)) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
      axis.text.y        = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x        = element_text(color = "black", size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title         = element_text(color = "black", size = 11, face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.6),
      axis.ticks         = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length  = unit(0.15, "cm"),
      legend.title = element_blank(),
      plot.margin        = margin(10, 14, 10, 14, "pt")
    )
}

build_dumbbell_plot <- function(df_wide, y_lab, y_lim, y_brk) {
  ggplot(df_wide) +
    geom_segment(aes(x = PLOT_LABEL, xend = PLOT_LABEL, y = `All Sites`, yend = `SW 90th`), color = "gray50", linewidth = 0.8) +
    geom_point(aes(x = PLOT_LABEL, y = `All Sites`, fill = "All Sites"), shape = 21, color = "black", size = 3.2, stroke = 0.5, alpha = 0.90) +
    geom_point(aes(x = PLOT_LABEL, y = `SW 90th`, fill = "SW 90th"), shape = 21, color = "black", size = 3.2, stroke = 0.5, alpha = 0.90) +
    scale_y_continuous(limits = y_lim, breaks = y_brk) +
    scale_fill_manual(values = c("All Sites" = "#1E293B", "SW 90th" = color_acc_orange)) +
    labs(x = NULL, y = y_lab) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "gray70"),
      axis.text.y        = element_text(color = "black", size = 10, face = "bold"),
      axis.text.x        = element_text(color = "black", size = 10, face = "bold", angle = 45, hjust = 1, vjust = 1),
      axis.title         = element_text(color = "black", size = 11, face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.6),
      axis.ticks         = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length  = unit(0.15, "cm"),
      legend.title = element_blank(),
      plot.margin        = margin(10, 14, 10, 14, "pt")
    )
}

p_sens_cov  <- build_sens_plot(coral_stats, "Coral Cover (%)", c(0, 100), seq(0, 100, 20))
p_sens_prod <- build_sens_plot(prod_stats, expression(bold("Gross Production Rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), c(0, 20), seq(0, 20, 5))

p_dumb_cov  <- build_dumbbell_plot(coral_stats, "Coral Cover (%)", c(0, 100), seq(0, 100, 20))
p_dumb_prod <- build_dumbbell_plot(prod_stats, expression(bold("Gross Production Rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")")), c(0, 20), seq(0, 20, 5))

figS4_combined <- ggarrange(
  no_x_title(p_sens_cov),
  no_x_title(p_dumb_cov),
  p_sens_prod,
  p_dumb_prod,
  ncol = 2,
  nrow = 2,
  align = "hv",
  labels = "auto",
  font.label = list(size = 12, face = "bold"),
  common.legend = TRUE,
  legend = "top"
)

ggsave(file.path(supp_dir, "FigureS4_Sensitivity.png"), figS4_combined, width = 28.0, height = 20.0, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS4_Sensitivity.tiff"), figS4_combined, width = 28.0, height = 20.0, units = "cm", dpi = 300, compression = "lzw", bg = "white")


# ==============================================================================
# FIGURE S5: Bioerosion Rates (Urchin & Fish)
# ==============================================================================
cat("Rendering Figure S5...\n")

df_urchin_island <- df_fixed_net %>%
  mutate(PLOT_LABEL = format_islands(LOCATION)) %>%
  group_by(PLOT_LABEL) %>%
  summarize(
    LATITUDE = mean(LATITUDE, na.rm = TRUE),
    URCHIN_EROSION_AVG = mean(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE),
    se.URCHIN_EROSION  = sd(URCHIN_EROSION_KG_M2_YR_MEAN, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

lbl_urchin <- expression(bold("Urchin erosion rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))
lbl_fish   <- expression(bold("Parrotfish erosion rate (kg" ~ m^{ -2 } ~ yr^{ -1 } * ")"))

p_urchin <- build_plot_f2(data = df_urchin_island, y_var = "URCHIN_EROSION_AVG", y_se = "se.URCHIN_EROSION", y_label = lbl_urchin, y_lim = c(0, 4), y_brk = seq(0, 4, 1), add_zero = FALSE)
p_fish   <- build_plot_f2(data = df_strs_darwin, y_var = "FISH_EROSION_KG_M2_YR_ALL", y_se = "se.FISH_EROSION_KG_M2_YR_ALL", y_label = lbl_fish, y_lim = c(0, 2), y_brk = seq(0, 2, 0.5), add_zero = FALSE)

main_plots_s5 <- ggarrange(no_x_title(p_urchin), p_fish, ncol = 1, nrow = 2, align = "v", labels = "auto", font.label = list(size = 12, face = "bold"), label.x = 0.01, label.y = 1.02)

island_locs_s5 <- df_strs_darwin %>% select(PLOT_LABEL, LATITUDE) %>% distinct(PLOT_LABEL, .keep_all = TRUE) %>% arrange(LATITUDE)

p_legend_s5 <- ggplot() +
  geom_tile(data = data.frame(lat = seq(18, 30, len = 500), x = 0), aes(x = x, y = lat, fill = lat), width = 0.5) +
  annotate("rect", xmin = -0.25, xmax = 0.25, ymin = 18, ymax = 30, color = "gray20", fill = NA, linewidth = 0.5) +
  scale_fill_gradientn(colors = pal_lat_grad, limits = c(18, 30), guide = "none") +
  scale_y_continuous(limits = c(18, 30), breaks = seq(18, 30, by = 2), expand = c(0.02, 0.02)) +
  scale_x_continuous(limits = c(-0.25, 3.8), expand = c(0, 0)) +
  geom_point(data = island_locs_s5, aes(x = 0.25, y = LATITUDE, fill = LATITUDE), shape = 21, size = 2, color = "gray30", stroke = 0.5) +
  geom_text_repel(data = island_locs_s5, aes(x = 0.25, y = LATITUDE, label = PLOT_LABEL), direction = "y", hjust = 0, nudge_x = 0.4, segment.color = "gray50", segment.linewidth = 0.4, size = 2.6, fontface = "bold", color = "gray20", force = 3, box.padding = unit(0.3, "lines")) +
  theme_void() +
  theme(
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 4)),
    axis.text.y = element_text(size = 7.5, face = "bold", color = "gray20", margin = margin(r = 5)),
    axis.ticks.y = element_line(color = "gray20"),
    axis.ticks.length = unit(0.15, "cm"),
    plot.margin = margin(15, 0, 15, 10)
  ) +
  labs(title = "Latitude (°N)")

figS5_final <- ggarrange(main_plots_s5, p_legend_s5, ncol = 2, widths = c(1, 0.35))

ggsave(file.path(supp_dir, "FigureS5_Bioerosion.png"), figS5_final, width = 18.0, height = 15.0, units = "cm", dpi = 300, bg = "white")
ggsave(file.path(supp_dir, "FigureS5_Bioerosion.tiff"), figS5_final, width = 18.0, height = 15.0, units = "cm", dpi = 300, compression = "lzw", bg = "white")

# ==============================================================================
# EXPORT TABLE 1
# ==============================================================================
target_docx <- file.path(main_dir, "Table1_IslandSummary.docx")

# 1. Delete pre-existing file so Pandoc does not throw overwrite error 22
if (file.exists(target_docx)) {
  unlink(target_docx, force = TRUE)
}

# 2. Save using explicit filename and path arguments to avoid CLI path truncation
table_saved <- FALSE
tryCatch(
  {
    gtsave(gt_journal_table, filename = "Table1_IslandSummary.docx", path = main_dir)
    table_saved <- file.exists(target_docx) && file.info(target_docx)$size > 0
  },
  error = function(e) {
    message("Notice: Pandoc DOCX export encountered an issue; switching to native Word RTF.")
  }
)

# 3. Native Word RTF Fallback (requires 0 Pandoc dependencies, opens directly in MS Word)
if (!table_saved) {
  target_rtf <- file.path(main_dir, "Table1_IslandSummary.rtf")
  if (file.exists(target_rtf)) unlink(target_rtf, force = TRUE)
  gtsave(gt_journal_table, filename = "Table1_IslandSummary.rtf", path = main_dir)
}
