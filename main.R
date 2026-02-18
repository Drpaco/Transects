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
# STEP 0A: AOI shapefile (optional; NAD83-first; Landmark skipped if AOI provided)
# ----------------------------
cat("STEP 0A: Optional AOI shapefile\n")
aoi_path <- get_input_safe("AOI polygon path (.shp/.gpkg) [Enter to skip]: ", "")
aoi_path <- trimws(aoi_path)
has_aoi  <- nzchar(aoi_path)

# We keep the DISPLAY/EXTENTS in NAD83 degrees (EPSG:4269)
# Measurements (km, m) will be done TEMPORARILY in a local UTM and converted back.

province_code <- NULL
aoi_ll    <- NULL    # AOI in lon/lat NAD83 (EPSG:4269)
aoi_terre <- NULL    # Alias: AOI in terre_crs (should also be 4269 in your setup)

# Helper: derive a local UTM for metrics (no permanent state)
.local_crs_m_from_bbox <- function(x_ll_4269) {
  bb <- sf::st_bbox(x_ll_4269)
  mid_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  mid_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  zone <- floor((mid_lon + 180) / 6) + 1
  if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
}

if (has_aoi) {
  message("Loading AOI: ", aoi_path)
  
  # 1) Read raw AOI
  aoi_raw <- suppressWarnings(sf::st_read(aoi_path, quiet = TRUE))
  aoi_raw <- sf::st_make_valid(aoi_raw)
  
  # 2) If CRS is missing, ask to ASSIGN the correct EPSG (do NOT transform yet)
  if (is.na(sf::st_crs(aoi_raw))) {
    epsg_in <- as.integer(get_input_safe(
      "AOI has no CRS. Enter EPSG to ASSIGN (e.g., 32198 = NAD83 / Quebec Lambert): ", "32198"))
    if (is.finite(epsg_in)) sf::st_crs(aoi_raw) <- epsg_in
  }
  
  # 3) Keep ONLY polygonal AOI; dissolve multipart
  aoi_poly <- suppressWarnings(sf::st_collection_extract(aoi_raw, "POLYGON"))
  if (nrow(aoi_poly) == 0) stop("AOI file contains no POLYGON geometries: ", aoi_path)
  aoi_poly <- sf::st_union(sf::st_geometry(aoi_poly)) |> sf::st_as_sf()
  
  # 4) Transform AOI to NAD83 (EPSG:4269) for ALL map logic/limits
  #    (Your terre_crs should be 4269 per your debug; if not, we still target 4269 for extents)
  aoi_ll <- sf::st_transform(aoi_poly, 4269)   # lon/lat NAD83 (degrees)
  
  # 5) Sanity: does the lon/lat NAD83 bbox look like Québec?
  bb_ll <- sf::st_bbox(aoi_ll)
  looks_like_quebec <- (as.numeric(bb_ll["xmin"]) < -40 && as.numeric(bb_ll["xmax"]) > -90 &&
                          as.numeric(bb_ll["ymin"]) < 70  && as.numeric(bb_ll["ymax"]) > 40)
  if (!looks_like_quebec) {
    message(sprintf(
      "AOI in NAD83(4269) looks off: [lon %.3f..%.3f, lat %.3f..%.3f]",
      bb_ll["xmin"], bb_ll["xmax"], bb_ll["ymin"], bb_ll["ymax"]
    ))
    epsg_in <- as.integer(get_input_safe(
      "Likely wrong ASSIGNED CRS. Enter correct EPSG to ASSIGN (e.g., 32198): ", "32198"))
    if (is.finite(epsg_in)) {
      sf::st_crs(aoi_poly) <- epsg_in   # assign, not transform
      aoi_ll <- sf::st_transform(aoi_poly, 4269)
      bb_ll  <- sf::st_bbox(aoi_ll)
    }
  }
  
  # 6) For SPANS (km) and CENTROID computations ONLY: hop to a local UTM (metric), then back to NAD83
  crs_m  <- .local_crs_m_from_bbox(aoi_ll)
  aoi_m  <- sf::st_transform(aoi_ll, crs_m)
  # AOI centroid back to NAD83 (will be used as landmark substitute in NAD83 degrees)
  pt_m   <- sf::st_centroid(sf::st_union(sf::st_geometry(aoi_m)))
  pt     <- sf::st_transform(sf::st_as_sf(pt_m), 4269)   # NAD83: keep everything in 4269
  
  # Compute span in km (AOI bbox measured in local UTM)
  bb_m   <- sf::st_bbox(aoi_m)
  span_km_x <- (as.numeric(bb_m["xmax"]) - as.numeric(bb_m["xmin"])) / 1000
  span_km_y <- (as.numeric(bb_m["ymax"]) - as.numeric(bb_m["ymin"])) / 1000
  
  # 7) Build NAD83 bbox/labels and all map artifacts in NAD83
  bb_final <- sf::st_bbox(aoi_ll)   # NAD83 degrees
  bb_txt <- sprintf(
    "AOI (NAD83): lon=%.4f..%.4f  lat=%.4f..%.4f  (~%.1f km x %.1f km)",
    bb_final["xmin"], bb_final["xmax"], bb_final["ymin"], bb_final["ymax"],
    span_km_x, span_km_y
  )
  
  # 8) Crop terre in NAD83 and build NAD83 limits with a small pad
  #    NOTE: terre is already in NAD83 per your debug; if not, we transform once here.
  if (!identical(sf::st_crs(terre), sf::st_crs(4269))) {
    terre <- sf::st_transform(terre, 4269)
  }
  aoi_terre <- aoi_ll
  bb_terre  <- sf::st_bbox(aoi_terre)    # NAD83 degrees
  terre_crop <- sf::st_crop(terre, bb_terre)
  
  pad_x <- as.numeric(bb_terre["xmax"] - bb_terre["xmin"]) * 0.05
  pad_y <- as.numeric(bb_terre["ymax"] - bb_terre["ymin"]) * 0.05
  xlim_terre <- c(as.numeric(bb_terre["xmin"]) - pad_x, as.numeric(bb_terre["xmax"]) + pad_x)
  ylim_terre <- c(as.numeric(bb_terre["ymin"]) - pad_y, as.numeric(bb_terre["ymax"]) + pad_y)
  
  # Graticule built in 4269, kept in 4269
  grat_final <- sf::st_graticule(bbox = bb_final, crs = sf::st_crs(4269), ndx = 8, ndy = 8)
  grat_final <- sf::st_sf(geometry = grat_final)  # already in 4269
  
  message("AOI mode: ON (Landmark skipped). Map/display locked to NAD83 (EPSG:4269).")
} else {
  message("AOI mode: OFF (will use Landmark flow). Map/display locked to NAD83 (EPSG:4269).")
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
# Use filtered airports for the main map (as you already do)
air_sf <- bbox_filter_points(air_sf_all, bb_final)
air_sf <- make_airport_labels(air_sf)
air_sf <- st_transform(air_sf, terre_crs)

if (has_aoi) {
  crs_m  <- "EPSG:3857"
  aoi_buf <- sf::st_transform(aoi_terre, crs_m) |> sf::st_buffer(25000) |> sf::st_transform(terre_crs)
  air_sf  <- air_sf[sf::st_intersects(air_sf, aoi_buf, sparse = FALSE), , drop = FALSE]
}

cat("\n--- DEBUG AIRPORTS ---\n")
print( AIRPORTS_CSV )
print( names(air_sf_all) )
# Show the key fields if they exist (will print NAs for missing columns, which is fine)
print( head(air_sf_all[, c("ident","icao","iata","code","airport_label","name")], 10) )

cat("Exact match ident == 'CYGP':\n")
print( subset(air_sf_all, ident == "CYGP") )

cat("Partial search for 'CYGP' in ident:\n")
print( air_sf_all[grepl("CYGP", air_sf_all$ident, ignore.case = TRUE), c("ident","name")] )

cat("--- END DEBUG AIRPORTS ---\n")

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
  show_full_terre = FALSE,            # crop terre to bbox during preview
  aoi_ll = if (has_aoi) aoi_ll else NULL
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
cat("\n--- DEBUG before STEP 5 ---\n")
print(sf::st_crs(terre_crs))
cat("terre_crop: n=", if (inherits(terre_crop,"sf")) nrow(terre_crop) else NA, "\n"); print(sf::st_crs(terre_crop))
cat("grat_final: n=", if (inherits(grat_final,"sf")) nrow(grat_final) else NA, "\n"); print(sf::st_crs(grat_final))
cat("transects_sf: n=", if (inherits(transects_sf,"sf")) nrow(transects_sf) else NA, "\n"); print(sf::st_crs(transects_sf))
cat("muni_sf: n=", if (inherits(muni_sf,"sf")) nrow(muni_sf) else NA, "\n"); print(sf::st_crs(muni_sf))
cat("air_sf: n=", if (inherits(air_sf,"sf")) nrow(air_sf) else NA, "\n"); print(sf::st_crs(air_sf))
if (exists("has_aoi") && has_aoi) { cat("aoi_terre: n=", nrow(aoi_terre), "\n"); print(sf::st_crs(aoi_terre)) }
cat("xlim_terre: ", paste(xlim_terre, collapse=" , "), "\n")
cat("ylim_terre: ", paste(ylim_terre, collapse=" , "), "\n")
cat("--- END DEBUG ---\n")

# Sanity: keep only finite limits
if (any(!is.finite(as.numeric(xlim_terre))) || any(!is.finite(as.numeric(ylim_terre)))) {
  warning("Invalid x/y limits passed to plot_final_map; deriving from terre_crop bbox.")
  bb <- sf::st_bbox(terre_crop)
  xlim_terre <- c(as.numeric(bb["xmin"]), as.numeric(bb["xmax"]))
  ylim_terre <- c(as.numeric(bb["ymin"]), as.numeric(bb["ymax"]))
}

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
  tran_label_offset_frac = 0.012,
  routes_sf = NULL,
  aoi_sf = if (has_aoi) aoi_terre else NULL
  
)
# if (has_aoi) {
#   gg_base <- gg_base + ggplot2::geom_sf(data = aoi_terre, fill = NA, color = "black", linewidth = 0.7)
# }
if (inherits(gg_base, "ggplot")) print(gg_base)

# ----------------------------
# STEP 6: Export transects to ESRI Shapefile
# ----------------------------
export_transects_shp(transects_sf, OUT_DIR_TRANSECTS, paste0(out_name, "_transects"))

# ----------------------------
# STEP 7 (Optional): Optimize flight routes
# ----------------------------

cat("lines:", nrow(transects_sf), "  types:", paste(unique(sf::st_geometry_type(transects_sf)), collapse=","), "\n")

# ---- Prepare a clean transect set for the optimizer ----
# Use NAD83 for display but compute lengths in a temporary UTM to drop zero-length pieces.
derive_crs_m <- function(x_ll_4269) {
  x_ll <- sf::st_transform(x_ll_4269, 4326)
  bb <- sf::st_bbox(x_ll)
  mid_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  mid_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  zone <- floor((mid_lon + 180) / 6) + 1
  if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
}

transects_opt <- transects_sf
transects_opt <- suppressWarnings(sf::st_collection_extract(transects_opt, "LINESTRING"))
transects_opt <- transects_opt[!sf::st_is_empty(transects_opt), , drop = FALSE]

if (nrow(transects_opt) == 0) {
  stop("No transects remain after AOI clipping. Increase transect length or expand AOI.")
}

crs_m_chk <- derive_crs_m(transects_opt)
len_m <- as.numeric(sf::st_length(sf::st_transform(transects_opt, crs_m_chk)))
# Keep strictly positive segments (>= 1 m)
transects_opt <- transects_opt[len_m >= 1, , drop = FALSE]

if (nrow(transects_opt) == 0) {
  stop("All transects are sub-meter after clipping. Increase 'len_km' or reduce clipping.")
}

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
  
  cat("OPT start/end (raw):", ap_start, "|", ap_end, "\n")
  
  message("Running optimizer...")
  opt_res <- run_transect_optimizer(
    transects_sf = transects_opt,   # <--- use the cleaned set
    air_sf       = air_all_for_opt, # full airports (not bbox-filtered)
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
    routes_sf = opt_res$routes_sf,
    aoi_sf = if (has_aoi) aoi_terre else NULL
  )
  # if (has_aoi) {
  #   gg_routes <- gg_routes + ggplot2::geom_sf(data = aoi_terre, fill = NA, color = "black", linewidth = 0.7)
  # }
  if (inherits(gg_routes, "ggplot")) print(gg_routes)
}