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

# Load modules from:
# 1) ROOT_DIR/R (if exists), else
# 2) ROOT_DIR/modules (if exists), else
# 3) ROOT_DIR itself (00_*.R etc)
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
      "Or place files named like 00_config.R, 01_utils.R directly in:\n",
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
terre <- terre_obj$terre
terre_crs <- terre_obj$crs

# ----------------------------
# STEP 0: Landmark -> BBox
# ----------------------------
cat("STEP 0: Landmark search and bbox selection\n")
landmark <- get_input_safe("Landmark search (e.g., Perce, Gaspe): ", "Perce")
province_code <- get_input_safe("Province code (Enter to skip; Quebec=24): ", "24")
province_code <- trimws(province_code)
if (province_code == "") province_code <- NULL

cand <- geonames_search(landmark, province = province_code, num = 20)
if (nrow(cand) == 0 && !is.null(province_code)) {
  cand <- geonames_search(landmark, province = NULL, num = 20)
}
pt <- choose_geoname(cand)

bb_final <- refine_bbox(pt, terre, initial_buffer_km = 25, grid_step_deg = 0.01)

span <- bbox_span_km(bb_final)
bb_txt <- sprintf(
  "BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f  (~%.1f km x %.1f km)",
  bb_final["xmin"], bb_final["ymin"], bb_final["xmax"], bb_final["ymax"],
  span["width_km"], span["height_km"]
)

# ----------------------------
# STEP 1: Airports
# ----------------------------
air_sf_all <- load_airports(AIRPORTS_CSV)   # keep the full set (unfiltered)
# Use filtered airports for the main map (as you already do)
air_sf <- bbox_filter_points(air_sf_all, bb_final)
air_sf <- make_airport_labels(air_sf)
air_sf <- st_transform(air_sf, terre_crs)

# ----------------------------
# STEP 2: Municipalities (NRCan GeoNames)
# ----------------------------
muni_raw <- geonames_bbox_paged(bb_final, province = province_code, concise = "CITY")
muni_ll  <- to_points_safe(muni_raw, bb_final)
muni_ll  <- bbox_filter_points(muni_ll, bb_final)
muni_sf  <- st_transform(muni_ll, terre_crs)

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

transects_sf <- st_transform(transects_ll, terre_crs)

# Create labels at each transect’s own east end (T1, T2, ...)
tran_labels_sf <- transect_labels_at_east_end(transects_sf, bb_final)
tran_labels_sf$label <- paste0("T", tran_labels_sf$id)

# ----------------------------
# STEP 4: Crop + graticule + limits
# ----------------------------
bb_terre <- st_bbox(st_transform(st_as_sfc(bb_final), terre_crs))
terre_crop <- st_crop(terre, bb_terre)

xlim_terre <- c(as.numeric(bb_terre["xmin"]), as.numeric(bb_terre["xmax"]))
ylim_terre <- c(as.numeric(bb_terre["ymin"]), as.numeric(bb_terre["ymax"]))

grat_final <- st_graticule(bbox = bb_final, crs = st_crs(4326), ndx = 8, ndy = 8)
grat_final <- st_transform(st_sf(geometry = grat_final), terre_crs)


# ============================
# OUTPUT NAMING & FOLDERS
# ============================
# Ask once and reuse for all outputs
out_name <- get_input_safe("Transect shapefile name (e.g., Perce_T_2026): ", "transects")
out_name <- trimws(out_name)

# Build two dedicated folders:
#  - Transects-only outputs (base shapefile + base map PNG)
#  - Flights (optimizer) outputs (routes shp/gpkg + trips CSV + routes PNG)
OUT_DIR_TRANSECTS <- file.path(TRANSECT_EXPORT_DIR, paste0(out_name, "_transects"))
OUT_DIR_FLIGHTS   <- file.path(TRANSECT_EXPORT_DIR, paste0(out_name, "_flights"))

# Ensure folders exist (helps GDAL and ggsave)
dir.create(OUT_DIR_TRANSECTS, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR_FLIGHTS,   showWarnings = FALSE, recursive = TRUE)

# Normalize paths (Windows/network paths and accents)
OUT_DIR_TRANSECTS <- normalizePath(OUT_DIR_TRANSECTS, winslash = "/", mustWork = FALSE)
OUT_DIR_FLIGHTS   <- normalizePath(OUT_DIR_FLIGHTS,   winslash = "/", mustWork = FALSE)

# Base transects PNG path (used in STEP 5)
out_png_transects <- file.path(OUT_DIR_TRANSECTS, paste0(out_name, "_transects.png"))

# ----------------------------
# STEP 5: Plot (base transects map)
# ----------------------------
gg_base <- plot_final_map(
  terre_crop = terre_crop,
  grat_final = grat_final,
  transects_sf = transects_sf,
  tran_labels_sf = tran_labels_sf,  # if you meant tran_labels_sf, keep as is; else typo fix below
  muni_sf = muni_sf,
  air_sf = air_sf,
  xlim_terre = xlim_terre,
  ylim_terre = ylim_terre,
  bb_txt = bb_txt,
  out_png = out_png_transects,  # save into the transects-only folder
  terre_crs = terre_crs,
  tran_label_offset_frac = 0.012
)
if (inherits(gg_base, "ggplot")) print(gg_base)


# ----------------------------
# STEP 6: Export transects to ESRI Shapefile
# ----------------------------
# out_name and OUT_DIR_TRANSECTS already defined above
export_transects_shp(transects_sf, OUT_DIR_TRANSECTS, paste0(out_name, "_transects"))


# ----------------------------
# STEP 7 (Optional): Optimize flight routes
# ----------------------------
opt_ans <- tolower(get_input_safe("Run flight route optimizer now? (y/N): ", "n"))
if (startsWith(opt_ans, "y")) {
  
  ap_start <- get_input_safe("Start airport code/name (Enter=nearest): ", "")
  ap_start <- trimws(ap_start); if (ap_start == "") ap_start <- NULL
  
  ap_end <- get_input_safe("End at a different airport? (code/name or Enter to return to start): ", "")
  ap_end <- trimws(ap_end); if (ap_end == "") ap_end <- NULL
  
  multi_in <- get_input_safe("Multi-airport mode? (Enter=No, 'auto 2' for K=2, or list e.g. 'CYGP,CYBC'): ", "")
  multi_in <- trimws(multi_in)
  
  speed_kn <- as.numeric(get_input_safe("Aircraft speed in knots [54]: ", "54"))
  max_hrs  <- as.numeric(get_input_safe("Max endurance per trip [2]: ", "2"))
  
  air_all_for_opt <- st_transform(air_sf_all, terre_crs)
  
  # flights folder
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
    airport_code = ap_start,             # start
    end_airport_code = ap_end,           # optional end (single-base mode)
    multi_airports   = multi_airports,   # NULL | "auto" | character vector
    multi_k          = multi_k,
    speed_knots  = speed_kn,
    max_hours    = max_hrs,
    export_dir   = OUT_DIR_FLIGHTS,
    out_basename = out_name,
    make_plot    = TRUE,
    terre_crop   = terre_crop,
    xlim_terre   = xlim_terre,
    ylim_terre   = ylim_terre,
    verbose      = TRUE                  # shows timer + progress bar
  )
  
  # Viewer + overlay map (combined)
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
  if (inherits(gg_routes, "ggplot")) print(gg_routes)
  # (plot_final_map auto-prints in Viewer with the earlier patch)
}
