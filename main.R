# ============================
# Adaptive project launcher
# ============================

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg)),
                         winslash = "/", mustWork = TRUE))
  }
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx <- tryCatch(rstudioapi::getActiveDocumentContext(), error = function(e) NULL)
    if (!is.null(ctx) && nzchar(ctx$path)) {
      return(normalizePath(dirname(ctx$path), winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ROOT_DIR <- get_script_dir()
message("Project root: ", ROOT_DIR)
setwd(ROOT_DIR)

source_utf8 <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  source(path, encoding = "UTF-8", local = FALSE)
}

# ----------------------------
# Module loader
# ----------------------------
source_modules <- function() {
  candidates <- c(
    file.path(ROOT_DIR, "R"),
    file.path(ROOT_DIR, "modules"),
    ROOT_DIR
  )
  mod_dir <- NULL
  for (d in candidates) {
    if (dir.exists(d)) {
      test <- list.files(d, pattern = "^[0-9]{2}_.+\\.R$", full.names = TRUE)
      if (length(test) > 0) { mod_dir <- d; break }
    }
  }
  if (is.null(mod_dir)) {
    stop(
      "No module folder found.\n",
      "Create one of these and put your module files inside:\n",
      "  - ", file.path(ROOT_DIR, "modules"), "\n",
      "  - ", file.path(ROOT_DIR, "R"), "\n",
      "Or place 00_config.R, 01_utils_*.R, etc. directly in:\n",
      "  - ", ROOT_DIR, "\n"
    )
  }
  files <- list.files(mod_dir, pattern = "^[0-9]{2}_.+\\.R$", full.names = TRUE)
  files <- files[order(basename(files))]
  message("Loading modules from: ", mod_dir)
  for (f in files) {
    message("  - ", basename(f))
    source_utf8(f)
  }
  invisible(files)
}

# ----------------------------
# Packages used across modules
# ----------------------------
library(sf)
library(dplyr)
library(ggplot2)
library(httr)
library(readr)
library(grid)

# ----------------------------
# Load your modular pipeline
# ----------------------------
source_modules()
message("plot_final_map loaded from: ", environmentName(environment(plot_final_map)))
print(args(plot_final_map))

# ----------------------------
# Load your transect designer module (separate file)
# ----------------------------
tdm <- file.path(ROOT_DIR, "TransectDesignerModule.R")
if (!file.exists(tdm)) stop("Missing TransectDesignerModule.R at: ", tdm)
source_utf8(tdm)

# ----------------------------
# Load terre
# ----------------------------
terre_obj <- load_terre(TERRE_PATH)
terre     <- terre_obj$terre
terre_crs <- terre_obj$crs

# ----------------------------
# Ensure exports root (fallback if modules didn't set it)
# ----------------------------
if (!exists("TRANSECT_EXPORT_DIR") || !nzchar(TRANSECT_EXPORT_DIR)) {
  TRANSECT_EXPORT_DIR <- file.path(ROOT_DIR, "exports")
  dir.create(TRANSECT_EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)
}

# ----------------------------
# STEP 0A: AOI shapefile (optional; takes precedence over Landmark)
# ----------------------------
cat("STEP 0A: Optional AOI shapefile\n")
aoi_path <- get_input_safe("AOI polygon path (.shp/.gpkg) [Enter to skip]: ", "")
aoi_path <- trimws(aoi_path)
has_aoi  <- nzchar(aoi_path)

province_code <- NULL
aoi_ll <- NULL
aoi_terre <- NULL

if (has_aoi) {
  message("Loading AOI: ", aoi_path)
  # Prefer helper if present; else read directly
  aoi <- if (exists("read_aoi_polygon")) {
    read_aoi_polygon(aoi_path)
  } else {
    suppressWarnings(sf::st_read(aoi_path, quiet = TRUE)) |>
      sf::st_make_valid() |>
      suppressWarnings(sf::st_collection_extract("POLYGON")) |>
      (\(x) { if (nrow(x) == 0) stop("AOI contains no POLYGON"); sf::st_union(sf::st_geometry(x)) |> sf::st_as_sf() })()
  }
  
  aoi_ll    <- sf::st_transform(aoi, 4326)
  aoi_terre <- sf::st_transform(aoi, terre_crs)
  
  bb_final <- sf::st_bbox(aoi_ll)
  span     <- bbox_span_km(bb_final)
  bb_txt   <- sprintf(
    "AOI: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f  (~%.1f km x %.1f km)",
    bb_final["xmin"], bb_final["ymin"], bb_final["xmax"], bb_final["ymax"],
    span["width_km"], span["height_km"]
  )
  
  # Crop terre & limits with padding
  bb_terre   <- sf::st_bbox(aoi_terre)
  terre_crop <- sf::st_crop(terre, bb_terre)
  pad_x <- as.numeric(bb_terre["xmax"] - bb_terre["xmin"]) * 0.05
  pad_y <- as.numeric(bb_terre["ymax"] - bb_terre["ymin"]) * 0.05
  xlim_terre <- c(as.numeric(bb_terre["xmin"]) - pad_x, as.numeric(bb_terre["xmax"]) + pad_x)
  ylim_terre <- c(as.numeric(bb_terre["ymin"]) - pad_y, as.numeric(bb_terre["ymax"]) + pad_y)
  
  grat_final <- sf::st_graticule(bbox = bb_final, crs = sf::st_crs(4326), ndx = 8, ndy = 8)
  grat_final <- sf::st_transform(sf::st_sf(geometry = grat_final), terre_crs)
  
  # Landmark substitute = AOI centroid via UTM to avoid lon/lat centroid warning
  bb      <- sf::st_bbox(aoi_ll)
  mid_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  mid_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  zone    <- floor((mid_lon + 180) / 6) + 1
  crs_m   <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
  pt_m    <- sf::st_centroid(sf::st_union(sf::st_transform(aoi_ll, crs_m)))
  pt      <- sf::st_transform(sf::st_as_sf(pt_m), 4326)
}

# ----------------------------
# STEP 0: Landmark -> BBox (SKIPPED if AOI provided)
# ----------------------------
if (!has_aoi) {
  cat("STEP 0: Landmark search and bbox selection\n")
  landmark <- get_input_safe("Landmark search (e.g., Perce, Gaspe): ", "Perce")
  province_code <- get_input_safe("Province code (Enter to skip; Quebec=24): ", "24")
  province_code <- trimws(province_code); if (province_code == "") province_code <- NULL
  
  cand <- geonames_search(landmark, province = province_code, num = 20)
  if (nrow(cand) == 0 && !is.null(province_code)) cand <- geonames_search(landmark, province = NULL, num = 20)
  pt <- choose_geoname(cand)
  
  bb_final <- refine_bbox(pt, terre, initial_buffer_km = 25, grid_step_deg = 0.01)
  
  span   <- bbox_span_km(bb_final)
  bb_txt <- sprintf(
    "BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f  (~%.1f km x %.1f km)",
    bb_final["xmin"], bb_final["ymin"], bb_final["xmax"], bb_final["ymax"],
    span["width_km"], span["height_km"]
  )
  
  bb_terre   <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb_final), terre_crs))
  terre_crop <- sf::st_crop(terre, bb_terre)
  xlim_terre <- c(as.numeric(bb_terre["xmin"]), as.numeric(bb_terre["xmax"]))
  ylim_terre <- c(as.numeric(bb_terre["ymin"]), as.numeric(bb_terre["ymax"]))
  
  grat_final <- sf::st_graticule(bbox = bb_final, crs = sf::st_crs(4326), ndx = 8, ndy = 8)
  grat_final <- sf::st_transform(sf::st_sf(geometry = grat_final), terre_crs)
}

# ----------------------------
# STEP 1: Airports
# ----------------------------
air_sf_all <- load_airports(AIRPORTS_CSV)
air_sf     <- bbox_filter_points(air_sf_all, bb_final) |> make_airport_labels() |> sf::st_transform(terre_crs)

if (has_aoi) {
  crs_m  <- "EPSG:3857"
  aoi_buf <- sf::st_transform(aoi_terre, crs_m) |> sf::st_buffer(25000) |> sf::st_transform(terre_crs)
  air_sf  <- air_sf[sf::st_intersects(air_sf, aoi_buf, sparse = FALSE), , drop = FALSE]
}

# ----------------------------
# STEP 2: Municipalities (NRCan GeoNames)
# ----------------------------
muni_raw <- geonames_bbox_paged(bb_final, province = province_code, concise = "CITY")
muni_ll  <- to_points_safe(muni_raw, bb_final) |> bbox_filter_points(bb_final)
muni_sf  <- sf::st_transform(muni_ll, terre_crs)

if (has_aoi) {
  crs_m  <- "EPSG:3857"
  aoi_buf <- sf::st_transform(aoi_terre, crs_m) |> sf::st_buffer(25000) |> sf::st_transform(terre_crs)
  muni_sf <- muni_sf[sf::st_intersects(muni_sf, aoi_buf, sparse = FALSE), , drop = FALSE]
}

# ----------------------------
# STEP 3: Transects (interactive)
# ----------------------------
transects_ll <- refine_transects(
  bb_ll = bb_final,
  terre = terre,
  pt_landmark_ll = pt,
  air_sf = air_sf,
  muni_sf = muni_sf,
  n = 10,
  len_km = 10,
  sp_km = 1,
  brg_deg = 90,
  clip_to_bbox = TRUE,
  show_full_terre = FALSE
)
transects_sf <- sf::st_transform(transects_ll, terre_crs)

# If AOI is provided, clip transects to AOI polygon (post design)
if (has_aoi) {
  tr_before    <- nrow(transects_sf)
  transects_sf <- clip_transects_to_aoi(transects_sf, aoi_terre)  # helper from modules/02_aoi_support.R
  if (nrow(transects_sf) == 0) {
    warning("All transects were outside AOI after clipping.")
  } else if (nrow(transects_sf) < tr_before) {
    message("Transects clipped to AOI: ", tr_before, " -> ", nrow(transects_sf))
  }
}

# ---- Normalize transect geometries after any clipping ----
# Make sure we only keep LINESTRING/MULTILINESTRING and drop empties
if (inherits(transects_sf, "sf")) {
  # Drop empty
  transects_sf <- transects_sf[!sf::st_is_empty(transects_sf), , drop = FALSE]
  
  # Extract LINESTRING parts from collections (if any)
  # This is robust when st_intersection() produced collections
  transects_sf <- suppressWarnings(sf::st_collection_extract(transects_sf, "LINESTRING"))
  
  # If still mixed, cast MULTILINESTRING -> LINESTRING parts
  # (This will split multi-lines into individual line pieces; IDs are preserved)
  if (!all(sf::st_geometry_type(transects_sf) %in% c("LINESTRING", "MULTILINESTRING"))) {
    transects_sf <- suppressWarnings(sf::st_cast(transects_sf, "LINESTRING"))
  }
  
  # Final guard: keep only non-zero length lines
  len_m <- as.numeric(sf::st_length(sf::st_transform(transects_sf, sf::st_crs(transects_sf))))
  transects_sf <- transects_sf[len_m > 0, , drop = FALSE]
}

# Optional: print to confirm
cat("transects_sf normalized to: "); print(unique(sf::st_geometry_type(transects_sf)))

# Labels at east end (recompute after potential AOI clip)
tran_labels_sf <- transect_labels_at_east_end(transects_sf, bb_final)
tran_labels_sf$label <- paste0("T", tran_labels_sf$id)

# ----------------------------
# STEP 4: Crop + graticule + limits (only recompute when NOT in AOI mode)
# ----------------------------
if (!has_aoi) {
  bb_terre   <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb_final), terre_crs))
  terre_crop <- sf::st_crop(terre, bb_terre)
  xlim_terre <- c(as.numeric(bb_terre["xmin"]), as.numeric(bb_terre["xmax"]))
  ylim_terre <- c(as.numeric(bb_terre["ymin"]), as.numeric(bb_terre["ymax"]))
  grat_final <- sf::st_graticule(bbox = bb_final, crs = sf::st_crs(4326), ndx = 8, ndy = 8)
  grat_final <- sf::st_transform(sf::st_sf(geometry = grat_final), terre_crs)
}

# ----------------------------
# OUTPUT NAMING & FOLDERS
# ----------------------------
out_name <- trimws(get_input_safe("Transect shapefile name (e.g., Perce_T_2026): ", "transects"))

OUT_DIR_TRANSECTS <- file.path(TRANSECT_EXPORT_DIR, paste0(out_name, "_transects"))
OUT_DIR_FLIGHTS   <- file.path(TRANSECT_EXPORT_DIR, paste0(out_name, "_flights"))
dir.create(OUT_DIR_TRANSECTS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_FLIGHTS,   showWarnings = FALSE, recursive = TRUE)

OUT_DIR_TRANSECTS <- normalizePath(OUT_DIR_TRANSECTS, winslash = "/", mustWork = FALSE)
OUT_DIR_FLIGHTS   <- normalizePath(OUT_DIR_FLIGHTS,   winslash = "/", mustWork = FALSE)
out_png_transects <- file.path(OUT_DIR_TRANSECTS, paste0(out_name, "_transects.png"))

# ----------------------------
# STEP 5: Plot (base transects map) – single call
# ----------------------------
gg_base <- plot_final_map(
  terre_crop = terre_crop,
  grat_final = grat_final,
  transects_sf = transects_sf,
  tran_labels_sf = tran_labels_sf,
  muni_sf = muni_sf,
  air_sf = air_sf,
  xlim_terre = xlim_terre,
  ylim_terre = ylim_terre,
  bb_txt = bb_txt,
  out_png = out_png_transects,   # save into the transects-only folder
  terre_crs = terre_crs,
  tran_label_offset_frac = 0.012
)
if (has_aoi) {
  gg_base <- gg_base + ggplot2::geom_sf(data = aoi_terre, fill = NA, color = "black", linewidth = 0.7)
}
if (inherits(gg_base, "ggplot")) print(gg_base)

# ----------------------------
# STEP 6: Export transects to ESRI Shapefile
# ----------------------------
export_transects_shp(transects_sf, OUT_DIR_TRANSECTS, paste0(out_name, "_transects"))

# ----------------------------
# STEP 7 (Optional): Optimize flight routes
# ----------------------------
opt_ans <- tolower(get_input_safe("Run flight route optimizer now? (y/N): ", "n"))
if (startsWith(opt_ans, "y")) {
  ap_start <- trimws(get_input_safe("Start airport code/name (Enter=nearest): ", ""))
  if (ap_start == "") ap_start <- NULL
  
  ap_end <- trimws(get_input_safe("End at a different airport? (code/name or Enter to return to start): ", ""))
  if (ap_end == "") ap_end <- NULL
  
  multi_in <- trimws(get_input_safe("Multi-airport mode? (Enter=No, 'auto 2' for K=2, or list 'CYGP,CYBC'): ", ""))
  
  speed_kn <- as.numeric(get_input_safe("Aircraft speed in knots [54]: ", "54"))
  max_hrs  <- as.numeric(get_input_safe("Max endurance per trip [2]: ", "2"))
  
  air_all_for_opt <- sf::st_transform(air_sf_all, terre_crs)
  dir.create(OUT_DIR_FLIGHTS, showWarnings = FALSE, recursive = TRUE)
  
  multi_airports <- NULL; multi_k <- 2
  if (nzchar(multi_in)) {
    if (grepl("^auto\\s+\\d+$", multi_in, ignore.case = TRUE)) {
      multi_k <- as.integer(sub(".*\\s+", "", multi_in))
      multi_airports <- "auto"
    } else {
      multi_airports <- strsplit(multi_in, "\\s*,\\s*")[[1]]
    }
  }
  
  message("Running optimizer...")
  opt_res <- run_transect_optimizer(
    transects_sf = transects_sf,
    air_sf       = air_all_for_opt,
    terre_crs    = terre_crs,
    airport_code = ap_start,
    end_airport_code = ap_end,
    multi_airports   = multi_airports,
    multi_k          = multi_k,
    speed_knots  = speed_kn,
    max_hours    = max_hrs,
    export_dir   = OUT_DIR_FLIGHTS,
    out_basename = out_name,
    make_plot    = TRUE,
    terre_crop   = terre_crop,
    xlim_terre   = xlim_terre,
    ylim_terre   = ylim_terre,
    verbose      = TRUE
  )
  
  gg_routes <- plot_final_map(
    terre_crop = terre_crop, grat_final = grat_final,
    transects_sf = transects_sf, tran_labels_sf = tran_labels_sf,
    muni_sf = muni_sf, air_sf = air_sf,
    xlim_terre = xlim_terre, ylim_terre = ylim_terre,
    bb_txt = paste0(bb_txt, " | Optimized routes"),
    out_png = file.path(TRANSECT_EXPORT_DIR, paste0(out_name, "_map_with_routes.png")),
    terre_crs = terre_crs,
    tran_label_offset_frac = 0.012,
    routes_sf = opt_res$routes_sf
  )
  if (has_aoi) {
    gg_routes <- gg_routes + ggplot2::geom_sf(data = aoi_terre, fill = NA, color = "black", linewidth = 0.7)
  }
  if (inherits(gg_routes, "ggplot")) print(gg_routes)
}