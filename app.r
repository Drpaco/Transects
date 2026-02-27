# app.R

needed <- c(
  "shiny","shinyWidgets","shinycssloaders","leaflet","leaflet.extras","DT",
  "sf","dplyr","ggplot2","readr","grid","zip"
)
to_install <- setdiff(needed, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(leaflet)
library(leaflet.extras)
library(DT)
library(sf)
library(dplyr)
library(ggplot2)
library(readr)
library(grid)
library(zip)

# Accept up to 200 MB uploads (adjust as needed)
options(shiny.maxRequestSize = 200 * 1024^2)

# ---- Cache location (local vs shinyapps.io) --
is_shinyapps <- nzchar(Sys.getenv("SHINY_PORT"))
MUNI_CACHE_FILE <- if (is_shinyapps) {
  file.path(tempdir(), "municipalities_cache.csv")
} else {
  file.path("data", "cache", "municipalities_cache.csv")
}
dir.create(dirname(MUNI_CACHE_FILE), recursive = TRUE, showWarnings = FALSE)

# --- Runtime options (no debugger halt unless explicitly enabled) -------------
if (interactive() && isTRUE(getOption("shiny.debug", FALSE))) {
  options(shiny.fullstacktrace = TRUE, shiny.error = browser)
} else {
  options(shiny.error = NULL)
}

# ---- Offline municipalities cache -------------------------------------------
MUNI_CACHE_FILE <- file.path("data", "cache", "municipalities_cache.csv")
dir.create(dirname(MUNI_CACHE_FILE), recursive = TRUE, showWarnings = FALSE)
# Use GEOS (planar) to avoid strict s2 errors during bbox/crop/intersects
sf::sf_use_s2(FALSE)

# --- Basemap selection (topo-friendly choices) --------------------------------
TOPO_CHOICES <- intersect(
  c("Esri.WorldTopoMap", "OpenTopoMap", "CartoDB.Voyager", "CartoDB.Positron", "Esri.WorldImagery"),
  names(leaflet::providers)
)
BASEMAP_DEFAULT <- if (length(TOPO_CHOICES)) TOPO_CHOICES[[1]] else "OpenStreetMap"

# --- Load config & your modules ------------------------------------------------
source("00_config.R", local = TRUE)

cat("WORKDIR:", normalizePath(getwd(), winslash="/"), "\n")
cat("TERRE_PATH:", TERRE_PATH, "\n")
cat("AIRPORTS_CSV:", AIRPORTS_CSV, "\n")
cat("TRANSECT_EXPORT_DIR:", TRANSECT_EXPORT_DIR, "\n")

safe_source <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  cat("SOURCING:", path, "\n")
  if (!file.exists(path)) stop("File not found: ", path)
  tryCatch(
    sys.source(path, envir = globalenv()),
    error = function(e) {
      stop("Error while sourcing '", path, "': ", conditionMessage(e))
    }
  )
}

safe_source("modules/01_utils_bbox.R")
safe_source("modules/02_terre.R")
safe_source("modules/03_geonames.R")
safe_source("modules/04_bbox_selector.R")
safe_source("modules/05_airports.R")
safe_source("modules/07_plot.R")
safe_source("modules/08_transect_optimizer.R")
safe_source("TransectDesignerModule.R")

# Use TERRE_PATH from 00_config.R (do not override it here)
terre_obj <- load_terre(TERRE_PATH)
terre     <- suppressWarnings(sf::st_make_valid(terre_obj$terre))
terre_crs <- terre_obj$crs   # keep NAD83 display in degrees

# --- Helpers -------------------------------------------------------------------
# ---- Municipalities offline cache helpers ------------------------------------

# Write a POINT sf (EPSG:4326) to CSV with stable columns
save_muni_cache <- function(muni_ll) {
  stopifnot(inherits(muni_ll, "sf"))
  muni_ll <- sf::st_transform(muni_ll, 4326)
  xy <- sf::st_coordinates(muni_ll)
  df <- cbind(
    data.frame(
      name     = as.character(muni_ll$name %||% muni_ll$NAME %||% muni_ll$title %||% ""),
      province = as.character(muni_ll$province %||% muni_ll$PR %||% muni_ll$prov %||% ""),
      lon      = xy[, 1],
      lat      = xy[, 2],
      stringsAsFactors = FALSE
    )
  )
  dir.create(dirname(MUNI_CACHE_FILE), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, MUNI_CACHE_FILE, row.names = FALSE, na = "")
  cat("DEBUG[save_muni_cache] wrote", nrow(df), "rows to", normalizePath(MUNI_CACHE_FILE, mustWork = FALSE), "\n")
}

# Read CSV -> POINT sf (EPSG:4326)
load_muni_cache_sf <- function() {
  if (!file.exists(MUNI_CACHE_FILE)) {
    cat("DEBUG[load_muni_cache_sf] no cache at", normalizePath(MUNI_CACHE_FILE, mustWork = FALSE), "\n")
    return(NULL)
  }
  df <- try(utils::read.csv(MUNI_CACHE_FILE, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(df, "try-error") || !nrow(df)) {
    cat("DEBUG[load_muni_cache_sf] cache read error or zero rows\n")
    return(NULL)
  }
  has_cols <- all(c("lon","lat") %in% names(df))
  if (!has_cols) {
    cat("DEBUG[load_muni_cache_sf] cache missing lon/lat columns\n")
    return(NULL)
  }
  sfobj <- sf::st_as_sf(df, coords = c("lon","lat"), crs = 4326, remove = FALSE)
  cat("DEBUG[load_muni_cache_sf] loaded", nrow(sfobj), "rows from cache\n")
  sfobj
}

# Offline-first with robust online fallback + ALWAYS save once on success
get_muni_points_offline_first <- function(query_bbox_ll, region_terre, terre_crs,
                                          province = NULL, concise = NULL, theme = "POPULATED_PLACE") {
  # 1) Try cache first
  cache_ll <- load_muni_cache_sf()
  if (!is.null(cache_ll)) {
    cache_terre <- try(sf::st_transform(cache_ll, terre_crs), silent = TRUE)
    if (!inherits(cache_terre, "try-error") && !is.null(region_terre)) {
      keep <- try(sf::st_intersects(cache_terre, region_terre, sparse = FALSE)[, 1], silent = TRUE)
      if (!inherits(keep, "try-error")) {
        n_total <- nrow(cache_terre); n_keep <- sum(keep)
        cat("DEBUG[get_muni] cache used: total =", n_total, " in region =", n_keep, "\n")
        return(cache_terre[keep, , drop = FALSE])
      }
    }
  }
  
  # 2) Primary online fetch (theme, broader than CITY only)
  cat("DEBUG[get_muni] cache missing/unusable; trying online (theme)\n")
  muni_raw <- try(
    geonames_bbox_paged(query_bbox_ll, province = province, concise = concise, theme = theme, num = 1000),
    silent = TRUE
  )
  rows_primary <- if (!inherits(muni_raw, "try-error") && !is.null(muni_raw)) nrow(muni_raw) else 0
  
  # 3) Fallback: union of concise classes if primary is empty
  if (rows_primary == 0) {
    cat("DEBUG[get_muni] primary returned 0; trying concise union with theme=NULL\n")
    concise_set <- c("CITY","TOWN","VILLAGE","HAMLET","LOCALITY")
    lst <- lapply(concise_set, function(cz)
      geonames_bbox_paged(
        bb_ll    = query_bbox_ll,
        province = NULL,
        concise  = cz,
        theme    = NULL,
        num      = 1000
      )
    )
    muni_raw <- suppressWarnings(do.call(rbind, Filter(Negate(is.null), lst)))
    if (is.null(muni_raw)) muni_raw <- sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
    cat("DEBUG[get_muni] union rows =", nrow(muni_raw), "\n")
  }
  
  # 4) Convert to POINTS in 4326
  if (is.null(muni_raw) || nrow(muni_raw) == 0) {
    cat("DEBUG[get_muni] no municipalities after fallback; returning NULL\n")
    return(NULL)
  }
  muni_ll <- to_points_safe(muni_raw, query_bbox_ll)
  if (is.null(muni_ll) || nrow(muni_ll) == 0) {
    cat("DEBUG[get_muni] to_points_safe produced no rows; returning NULL\n")
    return(NULL)
  }
  
  # 5) SAVE the cache once (POINTS in 4326)
  tryCatch(
    { save_muni_cache(muni_ll) },   # writes data/cache/municipalities_cache.csv (locally)
    error = function(e) cat("DEBUG[save_muni_cache] failed:", conditionMessage(e), "\n")
  )
  
  # 6) Transform to terre_crs and intersect with filter region
  muni_terre <- try(sf::st_transform(muni_ll, terre_crs), silent = TRUE)
  if (inherits(muni_terre, "try-error")) return(NULL)
  keep <- try(sf::st_intersects(muni_terre, region_terre, sparse = FALSE)[, 1], silent = TRUE)
  if (inherits(keep, "try-error")) return(muni_terre[0, , drop = FALSE])
  muni_terre[keep, , drop = FALSE]
}


# cat("DEBUG[muni] n_airports =", if (!is.null(rv$air_sf)) nrow(rv$air_sf) else 0,
#     " n_municipalities =", if (!is.null(rv$muni_sf)) nrow(rv$muni_sf) else 0, "\n")

# Build the effective filter region (in terre_crs) from UI state
# - If limit_to_aoi_buf == TRUE and AOI exists: AOI buffered by buf_m
# - Else: bbox polygon buffered by buf_m
# buf_km = 0 -> use auto-scaled buffer from bbox (scaled_buffer_m)
make_filter_region <- function(bb_final, aoi_terre, terre_crs, limit_to_aoi_buf, buf_km) {
  buf_m <- if (buf_km > 0) buf_km * 1000 else scaled_buffer_m(bb_final)
  
  if (isTRUE(limit_to_aoi_buf) && !is.null(aoi_terre)) {
    region_3857 <- sf::st_transform(aoi_terre, 3857) |>
      sf::st_buffer(buf_m)
    region <- sf::st_transform(region_3857, terre_crs)
  } else {
    bb_poly <- sf::st_as_sfc(bb_final)                 # 4326
    bb_3857 <- sf::st_transform(bb_poly, 3857) |>
      sf::st_buffer(buf_m)
    region <- sf::st_transform(bb_3857, terre_crs)
  }
  region
}

# Return a lon/lat (EPSG:4326) bbox for API queries, guaranteed finite
# Falls back to bb_final (w/ tiny padding) if region bbox is NA/Inf/empty
safe_query_bbox_ll <- function(region, bb_final) {
  ok <- function(x) is.numeric(x) && all(is.finite(x)) && !any(is.na(x))
  bb_try <- try({
    bb_ll <- region |>
      sf::st_transform(4326) |>
      sf::st_bbox()
    c(xmin = as.numeric(bb_ll["xmin"]), ymin = as.numeric(bb_ll["ymin"]),
      xmax = as.numeric(bb_ll["xmax"]), ymax = as.numeric(bb_ll["ymax"]))
  }, silent = TRUE)
  
  if (inherits(bb_try, "try-error") || !ok(bb_try)) {
    # Fallback: pad bb_final slightly so it's non-degenerate
    padx <- max(1e-4, 0.001 * abs(as.numeric(bb_final["xmax"] - bb_final["xmin"])))
    pady <- max(1e-4, 0.001 * abs(as.numeric(bb_final["ymax"] - bb_final["ymin"])))
    fb <- c(
      xmin = as.numeric(bb_final["xmin"]) - padx,
      ymin = as.numeric(bb_final["ymin"]) - pady,
      xmax = as.numeric(bb_final["xmax"]) + padx,
      ymax = as.numeric(bb_final["ymax"]) + pady
    )
    # DEBUG
    cat("DEBUG[query_bbox] Using fallback bb_final padded:",
        sprintf("xmin=%.6f ymin=%.6f xmax=%.6f ymax=%.6f\n", fb["xmin"], fb["ymin"], fb["xmax"], fb["ymax"]))
    return(sf::st_bbox(fb, crs = sf::st_crs(4326)))
  }
  
  # DEBUG
  cat("DEBUG[query_bbox] From region:",
      sprintf("xmin=%.6f ymin=%.6f xmax=%.6f ymax=%.6f\n", bb_try["xmin"], bb_try["ymin"], bb_try["xmax"], bb_try["ymax"]))
  sf::st_bbox(bb_try, crs = sf::st_crs(4326))
}

compute_limits_from_bbox <- function(bb_final, terre, terre_crs) {
  old_s2 <- sf::sf_use_s2()
  on.exit(sf::sf_use_s2(old_s2), add = TRUE)
  sf::sf_use_s2(FALSE)
  
  bb_terre   <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb_final), terre_crs))
  terre_valid <- tryCatch(sf::st_make_valid(terre), error = function(e) terre)
  
  terre_crop <- tryCatch(
    sf::st_crop(terre_valid, bb_terre),
    error = function(e) {
      warning("st_crop(terre, bb_terre) failed: ", conditionMessage(e))
      NULL
    }
  )
  
  pad_x <- as.numeric(bb_terre["xmax"] - bb_terre["xmin"]) * 0.05
  pad_y <- as.numeric(bb_terre["ymax"] - bb_terre["ymin"]) * 0.05
  xlim_terre <- c(as.numeric(bb_terre["xmin"]) - pad_x, as.numeric(bb_terre["xmax"]) + pad_x)
  ylim_terre <- c(as.numeric(bb_terre["ymin"]) - pad_y, as.numeric(bb_terre["ymax"]) + pad_y)
  
  list(bb_terre = bb_terre, terre_crop = terre_crop,
       xlim_terre = xlim_terre, ylim_terre = ylim_terre)
}

scaled_buffer_m <- function(bb_final, min_km = 15, frac = 0.08) {
  mid_lon <- as.numeric((bb_final["xmin"] + bb_final["xmax"]) / 2)
  mid_lat <- as.numeric((bb_final["ymin"] + bb_final["ymax"]) / 2)
  zone <- floor((mid_lon + 180) / 6) + 1
  crs_m <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
  p1 <- sf::st_sfc(sf::st_point(c(bb_final["xmin"], bb_final["ymin"])), crs = 4326) |> sf::st_transform(crs_m)
  p2 <- sf::st_sfc(sf::st_point(c(bb_final["xmax"], bb_final["ymax"])), crs = 4326) |> sf::st_transform(crs_m)
  diag_m <- as.numeric(sf::st_distance(p1, p2))
  max(min_km*1000, frac * diag_m)
}

filter_points_after_bbox <- function(pts_sf_all, bb_final, terre_crs, aoi_terre = NULL, use_aoi_buf = FALSE, buf_m = 25000) {
  pts <- bbox_filter_points(pts_sf_all, bb_final)
  pts <- sf::st_transform(pts, terre_crs)
  if (!is.null(aoi_terre) && use_aoi_buf) {
    aoi_buf <- sf::st_transform(aoi_terre, "EPSG:3857") |> sf::st_buffer(buf_m) |> sf::st_transform(terre_crs)
    keep <- sf::st_intersects(pts, aoi_buf, sparse = FALSE)
    pts  <- pts[keep, , drop = FALSE]
  }
  pts
}


# Convert a region (terre_crs) into a lon/lat query bbox for APIs (e.g., GeoNames)
region_to_query_bbox_ll <- function(region, terre_crs) {
  bb_ll <- region |> sf::st_transform(4326) |> sf::st_bbox()
  # Ensure bbox elements are named the usual way
  sf::st_bbox(c(xmin = bb_ll["xmin"], ymin = bb_ll["ymin"], xmax = bb_ll["xmax"], ymax = bb_ll["ymax"]), crs = sf::st_crs(4326))
}



# ---- SAFE reader: returns list(obj=sf, error=char) and never stops the app
read_vector_upload <- function(input_file) {
  stopifnot(!is.null(input_file))
  
  nm  <- input_file$name
  f   <- input_file$datapath
  ext <- tolower(tools::file_ext(nm))
  
  ok   <- function(sfobj, prj_missing = FALSE) list(obj = sfobj, error = NULL, prj_missing = prj_missing)
  fail <- function(msg)                       list(obj = NULL,  error = msg,  prj_missing = FALSE)
  
  cat("DEBUG[read_vector_upload] name:", nm, "ext:", ext, "size:", input_file$size, "\n")
  
  if (ext == "gpkg") {
    lyr_info <- try(sf::st_layers(f), silent = TRUE)
    if (inherits(lyr_info, "try-error")) {
      return(fail(paste("Unable to read GPKG layers:", as.character(lyr_info))))
    }
    geomtype_up <- toupper(lyr_info$geomtype)
    poly_idx    <- which(grepl("POLYGON", geomtype_up, fixed = TRUE))
    layer_name  <- if (length(poly_idx) >= 1) lyr_info$name[poly_idx[1]] else lyr_info$name[1]
    cat("DEBUG[read_vector_upload] gpkg layer chosen:", layer_name, "\n")
    
    sfobj <- try(sf::st_read(f, layer = layer_name, quiet = TRUE, options = "ENCODING=UTF-8"), silent = TRUE)
    if (inherits(sfobj, "try-error")) {
      return(fail(paste("Reading GPKG layer failed:", as.character(sfobj))))
    }
    return(ok(sfobj))
    
  } else if (ext == "zip") {
    td <- tempfile("unz"); dir.create(td, showWarnings = FALSE, recursive = TRUE)
    unzip(f, exdir = td)
    
    # Check for PRJ presence (just for user awareness)
    prj_files <- list.files(td, pattern = "\\.prj$", full.names = TRUE, recursive = TRUE)
    prj_missing <- length(prj_files) == 0
    
    shp <- list.files(td, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
    if (length(shp) == 0) {
      return(fail("No .shp found in ZIP. Zip the .shp, .shx, .dbf, .prj (and .cpg) together."))
    }
    cat("DEBUG[read_vector_upload] shp found:", shp[1], "\n")
    
    sfobj <- try(sf::st_read(shp[1], quiet = TRUE, options = "ENCODING=UTF-8"), silent = TRUE)
    if (inherits(sfobj, "try-error")) {
      return(fail(paste("Reading shapefile failed:", as.character(sfobj))))
    }
    return(ok(sfobj, prj_missing = prj_missing))
    
  } else if (ext %in% c("shp","dbf","shx","prj","cpg")) {
    return(fail("Upload shapefiles as a single ZIP (.zip) containing .shp, .shx, .dbf, .prj (and .cpg)."))
    
  } else {
    return(fail("Unsupported file. Upload a .gpkg or a zipped shapefile (.zip)."))
  }
}

# Landmark helpers: ensure a POINT in 4326, and outline polygon if available
landmark_as_point_ll <- function(x) {
  stopifnot(inherits(x, "sf"))
  x <- sf::st_make_valid(x)
  g <- sf::st_geometry(x)
  g <- sf::st_zm(g, drop = TRUE, what = "ZM")
  if (length(g) == 0 || all(sf::st_is_empty(g))) stop("Selected landmark has no geometry.")
  if (inherits(g, "sfc_POINT")) return(sf::st_transform(x, 4326))
  g_pt <- try(suppressWarnings(sf::st_point_on_surface(g)), silent = TRUE)
  if (inherits(g_pt, "try-error")) g_pt <- sf::st_centroid(g)
  sf::st_set_crs(g_pt, sf::st_crs(x)) |>
    sf::st_transform(4326) |>
    (\(geom) sf::st_set_geometry(x, geom))()
}
landmark_polygon_ll_or_null <- function(x) {
  stopifnot(inherits(x, "sf"))
  x <- sf::st_make_valid(x)
  gtype <- unique(as.character(sf::st_geometry_type(x, by_geometry = TRUE)))
  if (!any(grepl("POLYGON", gtype, fixed = TRUE))) return(NULL)
  suppressWarnings(sf::st_transform(x, 4326))
}

# --- UI ------------------------------------------------------------------------
ui <- navbarPage(
  
  title = tagList(
    span("Transect Planner"),
    tags$small(style="margin-left:6px; color:#6c757d;", em("by François Bolduc, Canadian Wildlife Service, Quebec Region"))
  ),
  id = "main_tabs",
  tabPanel("Area & BBox",
           fluidRow(
             column(4,
                    h4("1) Choose Area"),
                    # Basemap selector
                    selectInput("basemap", "Basemap style", choices = TOPO_CHOICES, selected = BASEMAP_DEFAULT),
                    
                    radioButtons("area_mode", "Mode", choices = c("AOI shapefile" = "aoi", "Landmark search" = "lm"), inline = TRUE),
                    conditionalPanel("input.area_mode == 'aoi'",
                                     fileInput(
                                       "aoi_file",
                                       "Upload AOI (.gpkg or zipped shapefile)",
                                       accept = c(".gpkg", ".zip")   # ⟵ only allow valid packages
                                     ),
                                     helpText("If shapefile, zip .shp/.shx/.dbf/.prj/.cpg and upload the .zip.")
                    ),
                    conditionalPanel("input.area_mode == 'lm'",
                                     textInput("lm_query", "Landmark search (e.g., Perce, Gaspe)", "Perce"),
                                     textInput("lm_province", "Province code (e.g., 24=QC, blank=Any)", ""),
                                     sliderInput("lm_buf_km", "Initial bbox buffer (km)", min = 0, max = 250, value = 10, step = 10),
                                     actionButton("lm_search", "Search Landmark")
                    ),
                    hr(),
                    h4("2) Refine BBox"),
                    prettySwitch("draw_rect_mode", "Draw rectangle", value = FALSE, status = "primary"),
                    actionButton("use_view_bbox", "Use current map view as BBox", icon = icon("crop")),
                    helpText("Pan/zoom or draw a rectangle → click 'Use current map view as BBox'."),
                    hr(),
                    actionButton("reset_bbox", "Start over (reset all)", icon = icon("undo"), class = "btn-warning"),
                    hr(),
                    verbatimTextOutput("bbox_text")
             ),
             column(8,
                    leafletOutput("map_bbox", height = 520) %>% withSpinner(color="#2C3E50")
             )
           )
  ),
  tabPanel("Points (Airports & Municipalities)",
           fluidRow(
             column(4,
                    h4("Points"),
                    prettySwitch("limit_to_aoi_buf", "Limit to AOI + buffer", value = TRUE, status = "primary"),
                    sliderInput("buf_km", "Buffer distance (km) — 0 = auto-scale", min = 0, max = 100, value = 0, step = 5),
                    actionButton("refresh_points", "Apply Filter", icon = icon("filter")),
                    hr(),
                    actionButton("reset_points", "Reset points", icon = icon("undo"), class = "btn-warning"),
                    hr(),
                    verbatimTextOutput("pts_counts")
             ),
             column(8,
                    leafletOutput("map_points", height = 520) %>% withSpinner(color="#2C3E50")
             )
           )
  ),
  tabPanel("Transects",
           fluidRow(
             column(4,
                    
                    # Scoped CSS for this tab: bold step headers, normal-weight input labels
                    tags$head(tags$style(HTML(
                      "#transects_tab .control-label, 
         #transects_tab label, 
         #transects_tab .help-block { font-weight: 400 !important; }
         #transects_tab summary { font-weight: 700 !important; font-size: 16px; }
         #transects_tab summary::-webkit-details-marker { color: #555; }"
                    ))),
                    
                    # All controls are inside this wrapper so CSS scoping works
                    div(id = "transects_tab",
                        
                        # --- Step 1 ------------------------------------------------------------
                        tags$details(
                                     tags$summary("1) Use existing transects"),
                                     br(),
                                     fileInput("tr_file", "Upload transects (.gpkg or zipped shapefile)",
                                               accept = c(".zip",".gpkg",".shp",".dbf",".shx",".prj",".cpg"))
                        ),
                        hr(),
                        
                        # --- Step 2 ------------------------------------------------------------
                        tags$details(
                                     tags$summary("2) Or build transects"),
                                     br(),
                                     numericInput("n_lines", "Number of lines", 10, min = 1, step = 1),
                                     numericInput("len_km", "Length (km)", 10, min = 0.1, step = 0.5),
                                     numericInput("sp_km", "Spacing (km)", 1, min = 0.1, step = 0.1),
                                     numericInput("brg_deg", "Bearing (deg) 0=N, 90=E", 90, min = 0, max = 359, step = 1),
                                     helpText("Anchor defaults to bbox center. Click on map to set anchor."),
                                     actionButton("build_tr", "Build/Refresh Transects", icon = icon("sliders"))
                        ),
                        hr(),
                        
                        # --- Step 3 ------------------------------------------------------------
                        tags$details(
                          tags$summary("3) Move transects"),
                          br(),
                          prettySwitch("enable_move", "Enable drag-to-move (anchor)", value = FALSE, status = "info"),
                          checkboxInput("rebuild_on_drag", "Rebuild lines when dragging (instead of shift)", value = FALSE),
                          actionButton("reset_anchor", "Reset anchor to bbox center", icon = icon("crosshairs"))
                        ),
                        hr(),
                        
                        # --- Step 4 ------------------------------------------------------------
                        tags$details(
                          tags$summary("4) Save transects"),
                          br(),
                          downloadButton("dl_transects_gpkg", "Download Transects (GPKG)"),
                          downloadButton("dl_transects_shp",  "Download Transects (Shapefile ZIP)")
                        ),
                        hr(),
                        
                        # --- Step 5 ------------------------------------------------------------
                        tags$details(
                          tags$summary("5) Subset transects"),
                          br(),
                          DTOutput("transects_table")
                        ),
                        hr(),
                        
                        # Reset button for this tab
                        actionButton("reset_transects", "Reset transects", icon = icon("undo"), class = "btn-warning")
                    )
             ),
             
             column(8,
                    leafletOutput("map_transects", height = 520) %>% withSpinner(color="#2C3E50")
             )
           )
  ),
  tabPanel("Optimize",
           fluidRow(
             column(4,
                    h4("Optimizer Settings"),
                    textInput("start_ap", "Start airport (code/name, Enter=nearest)", ""),
                    textInput("end_ap",   "End airport (code/name, blank=return to start)", ""),
                    textInput("multi_mode", "Multi-airports: blank=No, 'auto 3' or list 'CYGP,CYBC'", ""),
                    numericInput("speed_kn", "Speed (knots)", 54, min = 10, step = 1),
                    numericInput("max_hrs",  "Max endurance per trip (hours)", 2, min = 0.1, step = 0.1),
                    textInput("out_basename", "Output base name", "routes"),
                    actionButton("run_opt", "Run Optimizer", icon = icon("route")),
                    hr(),
                    # Downloads for final routes
                    downloadButton("dl_routes_gpkg", "Download Routes (GPKG)"),
                    downloadButton("dl_routes_shp",  "Download Routes (Shapefile ZIP)"),
                    downloadButton("dl_trips_csv",   "Download Trips CSV"),
                    downloadButton("dl_manifest_csv","Download Manifest CSV"),
                    hr(),
                    actionButton("reset_opt", "Reset optimizer outputs", icon = icon("undo"), class = "btn-warning")
             ),
             column(8,
                    # Leaflet map for optimizer (replaces static plot)
                    leafletOutput("map_opt", height = 520) %>% withSpinner(color="#2C3E50"),
                    hr(),
                    h5("Trips Summary"),
                    DTOutput("trips_table"),
                    hr(),
                    h5("Flight Manifest"),
                    DTOutput("manifest_table")
             )
           )
  ),
  tabPanel("Guide",
           fluidRow(
             column(
               12,
               # --- Scoped styles for readability ---
               tags$head(tags$style(HTML("
        #guide h3 { margin-top: 0.6rem; }
        #guide h4 { margin-top: 1.2rem; }
        #guide .callout { padding: 10px 12px; border-left: 4px solid #2C3E50; background: #f7f9fb; margin: 12px 0; }
        #guide code, #guide samp { background: #f1f3f5; padding: 0 4px; border-radius: 3px; }
        #guide details { margin: 10px 0 6px 0; }
        #guide summary { font-weight: 700; cursor: pointer; }
        #guide li { margin: 4px 0; }
        #guide .step { font-weight: 700; }
        #guide .ok { color: #2c7a7b; }
        #guide .warn { color: #b7791f; }
        #guide .err { color: #c53030; }
        #guide .kbd { padding:0 5px; border:1px solid #ccc; border-bottom-width:2px; border-radius:3px; background:#fff; }
      "))),
               div(
                 id = "guide",
                 
                 h3("How to use the Transect Planner"),
                 # --- Note (goes right after h3(...)) ---
                 div(class = "callout",
                     p(
                       "The Transect Planner builds ", strong("parallel transects only"), 
                       ", with a focus on eastern Canada. Use it to select your study area; load municipalities and airports; build or upload transects; subset transects for partial surveys; and generate an optimized route to survey the selected transects, based on your chosen start and end airports."
                     )                 ),
                 div(class = "callout",
                     p(strong("Quick start")),
                     tags$ol(
                       tags$li("Open ", strong("Area & BBox"), " → upload an AOI (Area of Interest)", em("(or)"),
                               " search a landmark → set your bounding box."),
                       tags$li("Switch to ", strong("Points"), " → set the buffer (0–250 km) and click ",
                               strong("Apply Filter"), " to load airports and municipalities."),
                       tags$li("Open ", strong("Transects"), " → either upload existing transects ",
                               em("or"), " set parameters and ", strong("Build/Refresh Transects"), "."),
                       tags$li("Drag the ", em("anchor"), " to move transects if needed (Step 3), then ",
                               strong("Save Transects"), " (Step 4)."),
                       tags$li("Go to ", strong("Optimize"), " → set flight settings → ", strong("Run Optimizer"),
                               " → view colored routes, legend, and download outputs.")
                     )
                 ),
                 
                 h4("Tabs overview"),
                 
                 tags$details(open = "open",
                              tags$summary("1) Area & BBox"),
                              tags$ul(
                                tags$li(strong("Basemap"), ": choose a background map from the dropdown."),
                                tags$li(strong("AOI shapefile/GeoPackage"), ": upload a polygon AOI. If the AOI has no CRS, a small prompt asks you to enter an EPSG (e.g., 32198 for Québec Lambert)."),
                                tags$li(strong("Landmark search"), ": search by name and (optionally) province code; pick a result to seed the BBox."),
                                tags$li(strong("Refine BBox"), ": pan/zoom or draw a rectangle, then click ", code("Use current map view as BBox"), "."),
                                tags$li(strong("What you see on the map"),
                                        tags$ul(
                                          tags$li("AOI outline (black), if provided"),
                                          tags$li("Bounding box (magenta)")
                                        )
                                ),
                                tags$li(strong("Reset"), ": ", code("Start over (reset all)"), " clears state and maps.")
                              )
                 ),
                 
                 tags$details(open = "open",
                              tags$summary("2) Points (Airports & Municipalities)"),
                              tags$ul(
                                tags$li(strong("Filter region"), ": controlled by the toggle ",
                                        code("Limit to AOI + buffer"),
                                        " and the buffer distance slider (", strong("0–250 km"), "). ",
                                        "When ON, the app uses ", em("AOI + buffer"), "; when OFF, it uses ", em("BBox + buffer"), "."),
                                tags$li(strong("Apply Filter"), ": loads airports and municipalities, then draws them with always-on labels."),
                                tags$li(strong("Filter outline"), ": a thin grey outline of the current filter region is visible on this tab only."),
                                tags$li(strong("Counts"), ": the panel shows how many airports and municipalities were retained."),
                                tags$li(strong("Reset"), ": ", code("Reset points"), " clears airports/municipalities from the map.")
                              ),
                              div(class="callout",
                                  p(strong("Offline tip (municipalities)")),
                                  p("If the GeoGratis GeoNames service is temporarily unavailable, the app can use an ",
                                    "offline CSV cache of municipalities (same points and labels). If enabled by your admin, ",
                                    "the first successful fetch writes:"),
                                  tags$ul(
                                    tags$li("Cache path: ", code(textOutput("guide_cache_path", inline = TRUE)))
                                  ),
                                  p("Future runs read this CSV and filter to your region—no web call needed. ",
                                    em("If you don’t see the CSV, it simply means the cache wasn’t enabled or saved yet."))
                              )
                 ),
                 
                 tags$details(open = "open",
                   tags$summary("3) Transects (5 steps)"),
                   p("This tab is organized as five numbered, collapsible steps. Steps 3–5 are hidden by default."),
                   tags$ol(
                     tags$li(span(class="step", "1) Use existing transects (optional)"), ": upload a ",
                             code(".gpkg"), " or zipped shapefile. The map draws them immediately."),
                     tags$li(span(class="step", "2) Or build transects"), ": set ",
                             code("Number of lines"), ", ", code("Length (km)"), ", ", code("Spacing (km)"),
                             ", ", code("Bearing (deg)"), "; click ", strong("Build/Refresh Transects"),
                             ". The anchor defaults to BBox center; click the map to set it elsewhere."),
                     tags$li(span(class="step", "3) Move transects"), ": toggle ",
                             code("Enable drag-to-move (anchor)"),
                             " to show a draggable anchor:",
                             tags$ul(
                               tags$li(code("Rebuild lines when dragging"), ": when checked, transects are recomputed using the current parameters; when unchecked, the current set is translated."),
                               tags$li(code("Reset anchor"), ": send the anchor back to BBox center.")
                             )),
                     tags$li(span(class="step", "4) Save transects"), ": download ",
                             code("GPKG"), " or a ", code("Shapefile ZIP"), "."),
                     tags$li(span(class="step", "5) Subset remaining transects"), ": select rows in the table to remove individual lines from the current set.")
                   ),
                   tags$ul(
                     tags$li(strong("Clipping to AOI"), ": if an AOI is present, newly built (and rebuild‑on‑drag) transects are clipped at the AOI boundary."),
                     tags$li(strong("Labels"), ": each transect shows an always‑visible ID label near its midpoint."),
                     tags$li(strong("Reset"), ": ", code("Reset transects"), " clears transects and the anchor.")
                   )
                 ),
                 
                 tags$details(open = "open",
                   tags$summary("4) Optimize"),
                   tags$ul(
                     tags$li(strong("Settings"), ": choose start/end airports (or leave end blank to return to start), multi‑airports, speed (knots), endurance (hours), and output base name."),
                     tags$li(strong("Run Optimizer"), ": shows a progress dialog; when finished, routes are drawn on the map in ", em("distinct colors per trip"), " with a legend."),
                     tags$li(strong("Downloads"), ": ",
                             code("Routes (GPKG, Shapefile ZIP)"), ", ",
                             code("Trips CSV"), ", ", code("Flight Manifest CSV"), "."),
                     tags$li(strong("Reset"), ": ", code("Reset optimizer outputs"), " clears routes and tables.")
                   )
                 ),
                 
                 tags$details(open = "open",
                   tags$summary("Troubleshooting & FAQs"),
                   tags$ul(
                     tags$li(strong("AOI upload closes the app? "), span(class="err","(rare)"),
                             " Your AOI may be corrupt or missing critical files. The app now guards against hard read errors and shows a notification instead—try zipping the entire shapefile set ",
                             em("(.shp, .shx, .dbf, .prj, .cpg)"), " and upload the ", code(".zip"), "."),
                     tags$li(strong("AOI has no CRS (.prj missing)"), ": a small CRS prompt appears; enter the correct EPSG (e.g., ", code("32198"), " for Québec Lambert) and the app will proceed and zoom to the AOI."),
                     tags$li(strong("Municipalities show 0"), ": the filter region may be too tight (try turning OFF ",
                             code("Limit to AOI + buffer"), " and increase buffer). If the GeoGratis service is down, ",
                             "municipalities won’t load unless you have the offline CSV cache."),
                     tags$li(strong("I don’t see labels"), ": zoom in a bit; also check the layer control—labels are separate overlay groups that can be toggled."),
                     tags$li(strong("No map on a tab"), ": each map renders immediately; if you still see a blank widget, switch tabs once or resize the window—rendering then triggers and layers appear.")
                   )
                 ),
                 
                 tags$details(open = "open",
                   tags$summary("Map tips"),
                   tags$ul(
                     tags$li("Use the layer control (top‑right) to toggle AOI, BBox, Airports, Municipalities, Labels, Transects, Routes, etc."),
                     tags$li("The scale bar (bottom‑left) shows distances in km."),
                     tags$li("Mouse wheel = zoom; drag = pan; click in Transects tab sets the anchor when move mode is enabled.")
                   )
                 )               )
             )
           )
  )
  
)

# --- SERVER --------------------------------------------------------------------
server <- function(input, output, session) {
  rv <- reactiveValues(
    aoi_ll = NULL, aoi_terre = NULL,
    bb_final = NULL, bb_terre = NULL,
    xlim_terre = NULL, ylim_terre = NULL,
    terre_crop = NULL, grat_final = NULL,
    air_all = NULL, air_sf = NULL,
    muni_sf = NULL, muni_raw = NULL,
    lm_results = NULL,
    anchor_lonlat = NULL,
    transects_sf = NULL,
    routes_sf = NULL, trip_summary = NULL, manifest_df = NULL
  )
  
  # --- Diagnostics: show the cache path and whether it exists
  observe({
    cat("MUNI_CACHE_FILE =", normalizePath(MUNI_CACHE_FILE, winslash = "/", mustWork = FALSE), "\n")
    cat("  exists? ", file.exists(MUNI_CACHE_FILE), "\n")
  })
  
  # ---- Reactive aliases (centralize reads of rv$ and input$) ----
  bb_final_rx <- reactive({ rv$bb_final })                       # bbox as a reactive
  aoi_terre_rx <- reactive({ rv$aoi_terre })                     # AOI in terre CRS as a reactive
  limit_to_aoi_rx <- reactive({ isTRUE(input$limit_to_aoi_buf) })
  buf_km_rx <- reactive({ input$buf_km })
  
  # Effective filter region (in terre_crs) as a reactive
  region_rx <- reactive({
    req(bb_final_rx())  # ensures not NULL
    make_filter_region(
      bb_final         = bb_final_rx(),
      aoi_terre        = aoi_terre_rx(),
      terre_crs        = terre_crs,
      limit_to_aoi_buf = limit_to_aoi_rx(),
      buf_km           = buf_km_rx()
    )
  })
  
  # Robust lon/lat query bbox for APIs (never NA/Inf)
  query_bbox_ll_rx <- reactive({
    req(bb_final_rx(), region_rx())
    safe_query_bbox_ll(region_rx(), bb_final_rx())
  })
  
  # Holds an AOI needing a CRS assignment
  rv$aoi_pending <- NULL
  
  # UI modal to ask for EPSG when AOI has no CRS
  show_aoi_crs_modal <- function() {
    showModal(modalDialog(
      title = "AOI: Coordinate Reference System (CRS) required",
      size = "m",
      tagList(
        p("The uploaded AOI has no CRS (projection). Please enter the EPSG code."),
        p(em("Tips for Québec data:"), " common CRS include ", code("EPSG:4326 (WGS84 lon/lat)"),
          ", ", code("EPSG:26920 (NAD83 / UTM zone 20N)"),
          ", ", code("EPSG:32188 (NAD83 / MTM zone 8)"), ", ", code("EPSG:32189 (MTM 9)"), "."),
        numericInput("aoi_epsg_input", "EPSG code",
                     value = 32198,  # default to Québec Lambert
                     min = 2000, max = 999999, step = 1)      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("aoi_set_epsg", "Apply CRS", class = "btn-primary")
      )
    ))
  }
  
  # Apply EPSG chosen in the modal
  observeEvent(input$aoi_set_epsg, {
    req(rv$aoi_pending)
    epsg <- as.integer(input$aoi_epsg_input)
    removeModal()
    if (is.na(epsg)) {
      showNotification("Invalid EPSG code.", type = "error", duration = 6)
      return(invisible(NULL))
    }
    try({
      sf::st_crs(rv$aoi_pending) <- epsg
    }, silent = TRUE)
    if (is.na(sf::st_crs(rv$aoi_pending))) {
      showNotification(paste("Could not set CRS to EPSG:", epsg), type = "error", duration = 6)
      return(invisible(NULL))
    }
    # Proceed with AOI now that CRS is set
    process_aoi(rv$aoi_pending)
    rv$aoi_pending <- NULL
  })
  
  # Utility: check map exists before proxy ops
  map_exists <- function(id) {
    !is.null(input[[paste0(id, "_zoom")]]) || !is.null(input[[paste0(id, "_bounds")]])
  }
  
  # --- Basemap switcher (safe; only updates maps that already exist) ----------
  observeEvent(input$basemap, ignoreInit = TRUE, {
    req(input$basemap)
    if (!(input$basemap %in% names(leaflet::providers))) return()
    for (map_id in c("map_bbox", "map_points", "map_transects", "map_opt")) {
      if (map_exists(map_id)) {
        try({
          leafletProxy(map_id, session = session) %>%
            clearTiles() %>%
            addProviderTiles(input$basemap)
        }, silent = TRUE)
      }
    }
  })
  
  # Base bbox map
  output$map_bbox <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(if (is.null(input$basemap)) BASEMAP_DEFAULT else input$basemap) |>
      addScaleBar(position = "bottomleft") |>
      addDrawToolbar(
        targetGroup = "draw",
        polygonOptions = FALSE, circleOptions = FALSE,
        circleMarkerOptions = FALSE, markerOptions = FALSE,
        polylineOptions = FALSE, rectangleOptions = TRUE,
        editOptions = editToolbarOptions()
      ) |>
      setView(lng = -64.5, lat = 48.8, zoom = 7)
  })
  # Ensure the BBox map renders even if the tab loads later
  outputOptions(output, "map_bbox", suspendWhenHidden = FALSE)
  
  # Initialize bb_final from the BBox map view as soon as bounds are available
  observeEvent(input$map_bbox_bounds, {
    if (is.null(rv$bb_final)) {
      b <- input$map_bbox_bounds
      rv$bb_final <- sf::st_bbox(
        c(xmin = b$west, ymin = b$south, xmax = b$east, ymax = b$north),
        crs = sf::st_crs(4326)
      )
      lims <- compute_limits_from_bbox(rv$bb_final, terre, terre_crs)
      rv$bb_terre   <- lims$bb_terre
      rv$terre_crop <- lims$terre_crop
      rv$xlim_terre <- lims$xlim_terre
      rv$ylim_terre <- lims$ylim_terre
      rv$grat_final <- sf::st_graticule(bbox = rv$bb_final, crs = sf::st_crs(4326), ndx = 8, ndy = 8) |>
        sf::st_transform(terre_crs) |>
        sf::st_sf()
    }
  }, ignoreInit = FALSE)
  
  # AOI upload
  # Processes a valid AOI sf into app state and redraws/zooms maps
  process_aoi <- function(aoi_sf) {
    aoi_raw <- sf::st_make_valid(aoi_sf)
    
    # Polygon-only
    aoi_poly <- suppressWarnings(sf::st_collection_extract(aoi_raw, "POLYGON"))
    if (nrow(aoi_poly) == 0) {
      showNotification("Uploaded layer has no POLYGON geometries.", type = "error", duration = 6)
      cat("DEBUG[AOI] No polygon features\n")
      return(invisible(NULL))
    }
    
    # Store AOI (lon/lat & projected)
    rv$aoi_ll    <- sf::st_transform(aoi_poly, 4269)     # NAD83 lon/lat
    rv$aoi_terre <- sf::st_transform(rv$aoi_ll, terre_crs)
    
    # Derive bbox + dependents
    rv$bb_final <- sf::st_bbox(rv$aoi_ll)
    lims <- compute_limits_from_bbox(rv$bb_final, terre, terre_crs)
    rv$bb_terre   <- lims$bb_terre
    rv$terre_crop <- lims$terre_crop
    rv$xlim_terre <- lims$xlim_terre
    rv$ylim_terre <- lims$ylim_terre
    rv$grat_final <- sf::st_graticule(bbox = rv$bb_final, crs = sf::st_crs(4269), ndx = 8, ndy = 8) |>
      sf::st_sf()
    
    # Draw AOI + bbox on the BBox map and FIT
    if (!is.null(input$map_bbox_bounds) || !is.null(input$map_bbox_zoom)) {
      leafletProxy("map_bbox", session = session) %>%
        clearShapes() %>%
        addPolygons(data = sf::st_transform(rv$aoi_ll, 4326),
                    color = "black", weight = 2, fill = FALSE) %>%
        addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                      lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                      color = "magenta", weight = 2) %>%
        fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                  rv$bb_final["xmax"], rv$bb_final["ymax"])
    }
    
    # Also trigger the other maps to draw/fit using existing observers/helpers
    # (They already listen to rv$bb_final and call redraw_* helpers)
    try(redraw_points_map(),    silent = TRUE)
    try(redraw_transects_map(), silent = TRUE)
    try(redraw_opt_map(),       silent = TRUE)
    
    showNotification("AOI loaded and zoomed.", type = "message", duration = 3)
  }
  
  # AOI upload observer
  observeEvent(input$aoi_file, ignoreInit = TRUE, {
    cat("DEBUG[AOI] upload name:", input$aoi_file$name, "size:", input$aoi_file$size, "\n")
    res <- read_vector_upload(input$aoi_file)
    
    if (!is.null(res$error)) {
      showNotification(paste("AOI upload failed:", res$error), type = "error", duration = 10)
      cat("DEBUG[AOI] ERROR:", res$error, "\n")
      return(invisible(NULL))
    }
    
    aoi_raw <- res$obj
    
    # If shapefile had no .prj, warn early (helps users pick a CRS correctly)
    if (isTRUE(res$prj_missing)) {
      showNotification("Warning: the shapefile ZIP contains no .prj. CRS is likely missing.",
                       type = "warning", duration = 8)
    }
    
    # If CRS is missing, ask user which EPSG to apply
    if (is.na(sf::st_crs(aoi_raw))) {
      rv$aoi_pending <- aoi_raw
      show_aoi_crs_modal()
      return(invisible(NULL))
    }
    
    # CRS present → process and zoom
    process_aoi(aoi_raw)
  })
  
  
  # Landmark search + modal picker (robust)
  observeEvent(input$lm_search, {
    req(nzchar(trimws(input$lm_query)))
    
    withProgress(message = "Searching landmarks…", value = 0, {
      incProgress(0.2)
      prov <- trimws(input$lm_province); if (prov == "") prov <- NULL
      
      cand <- tryCatch(
        geonames_search(input$lm_query, province = prov, num = 20),
        error = function(e) {
          showNotification(paste("Landmark search failed:", conditionMessage(e)),
                           type = "error", duration = 6)
          return(NULL)
        }
      )
      if (is.null(cand)) return(invisible(NULL))
      if (!nrow(cand)) {
        showNotification("No landmark results. Try a different name or clear the province code.",
                         type = "warning", duration = 5)
        return(invisible(NULL))
      }
      
      rv$lm_results <- cand
      showModal(modalDialog(
        title = "Select a landmark",
        size = "l",
        DTOutput("lm_table"),
        footer = tagList(actionButton("lm_pick_ok", "Use selected", class = "btn-primary"),
                         modalButton("Cancel"))
      ))
    })
  })
  
  output$lm_table <- renderDT({
    req(rv$lm_results)
    cols <- intersect(names(rv$lm_results), c("name","concise",'theme',"province","latitude","longitude","id","key"))
    dat <- as.data.frame(sf::st_drop_geometry(rv$lm_results[, cols, drop = FALSE]))
    datatable(dat, selection = "single", options = list(pageLength = 10))
  })
  
  observeEvent(input$lm_pick_ok, {
    req(rv$lm_results)
    sel <- input$lm_table_rows_selected
    if (length(sel) != 1) {
      showNotification("Select a single landmark.", type = "message")
      return(invisible(NULL))
    }
    removeModal()
    
    try({
      lm_feat    <- rv$lm_results[sel, , drop = FALSE]
      lm_pt_ll   <- landmark_as_point_ll(lm_feat)          # POINT for bbox + marker
      lm_poly_ll <- landmark_polygon_ll_or_null(lm_feat)   # outline if polygon
      
      # Build bbox from point
      rv$bb_final <- bbox_from_point_km(lm_pt_ll, buffer_km = input$lm_buf_km)
      
      # Recompute derived items
      lims <- compute_limits_from_bbox(rv$bb_final, terre, terre_crs)
      rv$bb_terre   <- lims$bb_terre
      rv$terre_crop <- lims$terre_crop
      rv$xlim_terre <- lims$xlim_terre
      rv$ylim_terre <- lims$ylim_terre
      rv$grat_final <- sf::st_graticule(bbox = rv$bb_final, crs = sf::st_crs(4326), ndx = 8, ndy = 8) |>
        sf::st_transform(terre_crs) |>
        sf::st_sf()
      
      p <- leafletProxy("map_bbox", session = session) %>% clearShapes()
      if (!is.null(lm_poly_ll)) {
        p <- p %>% addPolygons(data = lm_poly_ll, color = "#333333", weight = 2, fill = FALSE, opacity = 0.8)
      }
      p %>%
        addMarkers(data = lm_pt_ll) %>%
        addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                      lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                      color = "magenta", weight = 2) %>%
        fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                  rv$bb_final["xmax"], rv$bb_final["ymax"])
    }, silent = TRUE)
  })
  
  # Use current map view as bbox
  observeEvent(input$use_view_bbox, {
    b <- input$map_bbox_bounds
    validate(need(!is.null(b), "No map bounds available"))
    rv$bb_final <- sf::st_bbox(c(xmin = b$west, ymin = b$south, xmax = b$east, ymax = b$north), crs = sf::st_crs(4326))
    
    lims <- compute_limits_from_bbox(rv$bb_final, terre, terre_crs)
    rv$bb_terre   <- lims$bb_terre
    rv$terre_crop <- lims$terre_crop
    rv$xlim_terre <- lims$xlim_terre
    rv$ylim_terre <- lims$ylim_terre
    rv$grat_final <- sf::st_graticule(bbox = rv$bb_final, crs = sf::st_crs(4269), ndx = 8, ndy = 8) |> sf::st_sf()
    
    leafletProxy("map_bbox", session = session) %>%
      clearGroup("bbox") %>%
      addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                    lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                    color = "magenta", weight = 2, group = "bbox")
  })
  
  output$bbox_text <- renderText({
    req(rv$bb_final)
    sprintf("BBox: lon=%.4f..%.4f, lat=%.4f..%.4f",
            rv$bb_final["xmin"], rv$bb_final["xmax"], rv$bb_final["ymin"], rv$bb_final["ymax"])
  })
  
  # --- RESET: AREA/BBOX (Start over) ------------------------------------------
  observeEvent(input$reset_bbox, {
    # Clear all reactive state
    rv$aoi_ll <- rv$aoi_terre <- NULL
    rv$bb_final <- rv$bb_terre <- NULL
    rv$xlim_terre <- rv$ylim_terre <- NULL
    rv$terre_crop <- rv$grat_final <- NULL
    rv$air_sf <- rv$muni_sf <- NULL
    rv$transects_sf <- NULL
    rv$anchor_lonlat <- NULL
    rv$routes_sf <- rv$trip_summary <- rv$manifest_df <- NULL
    
    # Clear maps safely
    for (map_id in c("map_bbox","map_points","map_transects","map_opt")) {
      if (map_exists(map_id)) {
        try({
          leafletProxy(map_id, session = session) %>%
            clearShapes() %>% clearMarkers() %>%
            clearGroup("aoi") %>% clearGroup("bbox") %>%
            clearGroup("airports") %>% clearGroup("municipalities") %>%
            clearGroup("airports_labels") %>% clearGroup("municipalities_labels") %>%
            clearGroup("transects") %>% clearGroup("transect_labels") %>%
            clearGroup("anchor") %>% clearGroup("routes") %>%
            clearGroup("filter_region") %>%
            setView(lng = -64.5, lat = 48.8, zoom = 6)
        }, silent = TRUE)
      }
    }
    showNotification("Application reset.", type = "message", duration = 4)
  })
  
  # --- POINTS TAB --------------------------------------------------------------
  output$map_points <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(if (is.null(input$basemap)) BASEMAP_DEFAULT else input$basemap) |>
      addScaleBar(position = "bottomleft") |>
      addLayersControl(
        baseGroups = NULL,
        overlayGroups = c(
          "aoi", "bbox",
          "airports", "municipalities",
          "airports_labels", "municipalities_labels",
          "filter_region"          
        ),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      setView(lng = -64.5, lat = 48.8, zoom = 6)
  })
  
  outputOptions(output, "map_points", suspendWhenHidden = FALSE)
  
  # Safely add always-on text labels (works even on older leaflet versions)
  safe_add_label_only <- function(p, data, label, group, textsize = "11px", offset = c(10, -10)) {
    if ("addLabelOnlyMarkers" %in% getNamespaceExports("leaflet")) {
      p %>% leaflet::addLabelOnlyMarkers(
        data = data,
        label = label,
        labelOptions = labelOptions(
          noHide = TRUE, textOnly = TRUE, direction = "auto",
          textsize = textsize, offset = offset, opacity = 0.9
        ),
        group = group
      )
    } else {
      p %>% addCircleMarkers(
        data = data,
        radius = 1, stroke = FALSE, fillOpacity = 0,
        label = label,
        labelOptions = labelOptions(
          noHide = TRUE, direction = "auto",
          textsize = textsize, offset = offset, opacity = 0.9
        ),
        group = group
      )
    }
  }
  
  redraw_points_map <- function() {
    if (!map_exists("map_points")) return(invisible(NULL))
    
    p <- leafletProxy("map_points", session = session) %>%
      clearGroup("aoi") %>%
      clearGroup("bbox") %>%
      clearGroup("airports") %>%
      clearGroup("municipalities") %>%
      clearGroup("airports_labels") %>%
      clearGroup("municipalities_labels") %>%
      clearGroup("filter_region")   # outline is only on the POINTS map
    
    # AOI outline
    if (!is.null(rv$aoi_ll)) {
      p <- p %>% addPolygons(
        data  = sf::st_transform(rv$aoi_ll, 4326),
        color = "black", weight = 2, fill = FALSE, group = "aoi"
      )
    }
    
    # BBox rectangle
    if (!is.null(bb_final_rx())) {
      bb <- bb_final_rx()
      p <- p %>% addRectangles(
        lng1 = bb["xmin"], lat1 = bb["ymin"],
        lng2 = bb["xmax"], lat2 = bb["ymax"],
        color = "magenta", weight = 2, group = "bbox"
      )
    }
    
    # Airports + labels
    if (!is.null(rv$air_sf) && nrow(rv$air_sf) > 0) {
      air_ll <- sf::st_transform(rv$air_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = air_ll,
        radius = 4, color = "blue",
        label  = ~airport_label,
        group  = "airports"
      )
      p <- safe_add_label_only(
        p, data = air_ll, label = ~airport_label,
        group = "airports_labels", textsize = "11px", offset = c(10, -10)
      )
    }
    
    # Municipalities + labels
    if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0) {
      muni_ll <- sf::st_transform(rv$muni_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = muni_ll,
        radius = 3, color = "red",
        label  = ~as.character(name),
        group  = "municipalities"
      )
      p <- safe_add_label_only(
        p, data = muni_ll, label = ~as.character(name),
        group = "municipalities_labels", textsize = "10px", offset = c(8, -8)
      )
    }
    
    # ---------- FILTER OUTLINE (only on the POINTS map) ----------
    if (!is.null(bb_final_rx())) {
      region <- region_rx()  # uses your reactive alias
      p <- p %>% addPolygons(
        data  = sf::st_transform(region, 4326),
        color = "#666", weight = 1, fill = FALSE,
        group = "filter_region"
      )
    }
    # -------------------------------------------------------------
    
    # Fit to bbox if available
    if (!is.null(bb_final_rx())) {
      bb <- bb_final_rx()
      p %>% fitBounds(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
    }
  }
  
  # Fit & draw when bbox changes (Points map)
  observeEvent(rv$bb_final, {
    if (map_exists("map_points")) {
      leafletProxy("map_points", session = session) %>%
        fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                  rv$bb_final["xmax"], rv$bb_final["ymax"])
      try(redraw_points_map(), silent = TRUE)
    }
  }, ignoreInit = TRUE)
  
  # Apply Filter (robust)
  observeEvent(input$refresh_points, ignoreInit = TRUE, {
    req(rv$bb_final)  # We at least need an initial bbox to scale a default buffer
    
    withProgress(message = "Applying point filters…", value = 0, {
      on.exit({ incProgress(1) }, add = TRUE)
      
      incProgress(0.05, detail = "Loading airports")
      if (is.null(rv$air_all)) {
        rv$air_all <- tryCatch(
          load_airports(AIRPORTS_CSV),
          error = function(e) {
            showNotification(
              paste("Failed to load airports from", AIRPORTS_CSV, ":", conditionMessage(e)),
              type = "error", duration = 8
            )
            return(NULL)
          }
        )
      }
      if (is.null(rv$air_all)) return(invisible(NULL))
      
      # ---- Build the effective filter region + a robust lon/lat query bbox ----
      region <- region_rx()
      query_bbox_ll <- safe_query_bbox_ll(region, bb_final_rx())
      
      cat(sprintf("DEBUG[bbox_ll] west=%.6f south=%.6f east=%.6f north=%.6f\n",
                  query_bbox_ll["xmin"], query_bbox_ll["ymin"], query_bbox_ll["xmax"], query_bbox_ll["ymax"]))
      
      # ---- Airports: transform & intersect with the same region ----
      incProgress(0.35, detail = "Filtering airports")
      air_terre <- tryCatch(sf::st_transform(rv$air_all, terre_crs), error = function(e) NULL)
      if (is.null(air_terre)) return(invisible(NULL))
      
      keep_air <- tryCatch(
        sf::st_intersects(air_terre, region, sparse = FALSE)[, 1],
        error = function(e) rep(FALSE, nrow(air_terre))
      )
      rv$air_sf <- air_terre[keep_air, , drop = FALSE]
      
      if (nrow(rv$air_sf) > 0) {
        rv$air_sf <- tryCatch(make_airport_labels(rv$air_sf), error = function(e) rv$air_sf)
      }
      
      # # ---- Municipalities: OFFLINE-FIRST (CSV cache) -------------------------------
      # incProgress(0.65, detail = "Preparing municipalities (offline-first)")
      # rv$muni_sf <- tryCatch(
      #   get_muni_points_offline_first(
      #     query_bbox_ll = query_bbox_ll,     # from your safe_query_bbox_ll(region, bb)
      #     region_terre  = region,            # same region used for airports filtering
      #     terre_crs     = terre_crs,
      #     province      = NULL,
      #     concise       = NULL,
      #     theme    = "POPULATED_PLACE"
      #   ),
      #   error = function(e) {
      #     showNotification(paste("Municipality fetch (offline-first) failed:", conditionMessage(e)),
      #                      type = "error", duration = 8)
      #     NULL
      #   }
      # )
      # 
      
      # ---- Municipalities: DIAGNOSTIC (union of concise classes; online-only) ----
      incProgress(0.65, detail = "Municipalities (diagnostic union)")
      
      # 1) Pull several concise classes and union (theme=NULL for this test)
      concise_set <- c("CITY","TOWN","VILLAGE","HAMLET","LOCALITY")
      lst <- lapply(concise_set, function(cz)
        geonames_bbox_paged(
          bb_ll    = query_bbox_ll,
          province = NULL,
          concise  = cz,
          theme    = NULL,   # IMPORTANT: remove theme for the test
          num      = 1000
        )
      )
      
      # Keep only non-null pages
      lst <- Filter(Negate(is.null), lst)
      
      # Safely row-bind sf pages (preserves geometry). If there are none, make an empty sf with a geometry column.
      muni_raw <- if (length(lst)) {
        suppressWarnings(do.call(rbind, lst))
      } else {
        sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
      }
      
      cat("DEBUG[test] union rows =", nrow(muni_raw), "\n")
      
      # 2) Convert to POINTS, transform to terre_crs, intersect with region
      #    (Start with a valid empty sf so we never fail on assignment.)
      rv$muni_sf <- sf::st_sf(geometry = sf::st_sfc(), crs = 4326)
      
      if (!is.null(muni_raw) && nrow(muni_raw) > 0) {
        muni_ll <- to_points_safe(muni_raw, query_bbox_ll)  # returns POINT sf in 4326
        if (!is.null(muni_ll) && nrow(muni_ll) > 0) {
          muni_terre <- sf::st_transform(muni_ll, terre_crs)
          keep <- try(sf::st_intersects(muni_terre, region, sparse = FALSE)[,1], silent = TRUE)
          if (!inherits(keep, "try-error") && any(keep)) {
            rv$muni_sf <- muni_terre[keep, , drop = FALSE]
          } else {
            rv$muni_sf <- muni_terre[0, , drop = FALSE]
          }
        }
      }
      
      # Ensure a 'name' column exists for labels
      if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0 && !("name" %in% names(rv$muni_sf))) {
        cand <- c("GEONAME","NAME","title")
        hit  <- intersect(cand, names(rv$muni_sf))
        rv$muni_sf$name <- if (length(hit)) as.character(rv$muni_sf[[hit[1]]]) else ""
      }
      
      cat("DEBUG[test] rv$muni_sf rows =", nrow(rv$muni_sf), "\n")
      
      
      
      
      
      cat("DEBUG[muni] rv$muni_sf rows = ", if (is.null(rv$muni_sf)) 0 else nrow(rv$muni_sf), "\n")
      
      if (file.exists(MUNI_CACHE_FILE)) {
        info <- file.info(MUNI_CACHE_FILE)
        cat("DEBUG[muni-cache] exists, size(bytes) =", info$size, "\n")
      } else {
        cat("DEBUG[muni-cache] does NOT exist\n")
      }
      
      # --- ALWAYS save a cache copy when we have municipalities ---
      if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0) {
        # Ensure a 'name' column (for labels and for cache)
        if (!("name" %in% names(rv$muni_sf))) {
          cand <- c("GEONAME","NAME","title")
          hit  <- intersect(cand, names(rv$muni_sf))
          rv$muni_sf$name <- if (length(hit)) as.character(rv$muni_sf[[hit[1]]]) else ""
        }
        # Ensure a 'province' column (cache schema)
        if (!("province" %in% names(rv$muni_sf))) {
          rv$muni_sf$province <- ""
        }
        
        # Save cache as POINTS in EPSG:4326 (the cache schema)
        muni_ll_for_cache <- try(sf::st_transform(rv$muni_sf, 4326), silent = TRUE)
        if (!inherits(muni_ll_for_cache, "try-error")) {
          tryCatch(
            {
              save_muni_cache(muni_ll_for_cache)   # writes MUNI_CACHE_FILE
              cat("DEBUG[cache-write] wrote file ->", normalizePath(MUNI_CACHE_FILE, mustWork = FALSE), "\n")
              cat("DEBUG[cache-write] exists? ", file.exists(MUNI_CACHE_FILE), "\n")
            },
            error = function(e) {
              cat("DEBUG[cache-write] failed:", conditionMessage(e), "\n")
            }
          )
        } else {
          cat("DEBUG[cache] transform to 4326 failed; no cache written\n")
        }
      }
      
      # ---- Redraw the Points map (proxy) + feedback ----
      redraw_points_map()
      n_air <- if (!is.null(rv$air_sf)) nrow(rv$air_sf) else 0
      n_mun <- if (!is.null(rv$muni_sf)) nrow(rv$muni_sf) else 0
      showNotification(sprintf("Filtered: %d airports, %d municipalities.", n_air, n_mun),
                       type = "message", duration = 4)
    })
  })
  
  
  # bb_dbg <- sf::st_bbox(sf::st_transform(region, 4326))
  # cat("DEBUG[region bbox] xmin=", bb_dbg["xmin"], " ymin=", bb_dbg["ymin"],
  #     " xmax=", bb_dbg["xmax"], " ymax=", bb_dbg["ymax"], "\n")
  
  output$pts_counts <- renderText({
    a <- if (!is.null(rv$air_sf)) nrow(rv$air_sf) else 0
    m <- if (!is.null(rv$muni_sf)) nrow(rv$muni_sf) else 0
    if (a + m == 0) "Click 'Apply Filter'" else sprintf("Airports: %d | Municipalities: %d", a, m)
  })
  
  
  
  # Reset points
  observeEvent(input$reset_points, {
    rv$air_sf <- rv$muni_sf <- NULL
    redraw_points_map()
    showNotification("Points reset.", type = "message", duration = 3)
  })
  
  # --- TRANSECTS TAB -----------------------------------------------------------
  output$map_transects <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(if (is.null(input$basemap)) BASEMAP_DEFAULT else input$basemap) |>
      addScaleBar(position = "bottomleft") |>
      addLayersControl(
        baseGroups = NULL,
        overlayGroups = c("aoi", "bbox", "airports", "municipalities",
                          "airports_labels", "municipalities_labels",
                          "transects", "transect_labels", "anchor"),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      setView(lng = -64.5, lat = 48.8, zoom = 6)
  })
  outputOptions(output, "map_transects", suspendWhenHidden = FALSE)
  
  # Redraw helper for transects map (includes AOI, bbox, points, transects, labels, anchor)
  redraw_transects_map <- function(fit = FALSE) {
    if (!map_exists("map_transects")) return(invisible(NULL))
    
    p <- leafletProxy("map_transects", session = session) %>%
      clearGroup("aoi") %>%
      clearGroup("bbox") %>%
      clearGroup("airports") %>%
      clearGroup("municipalities") %>%
      clearGroup("airports_labels") %>%
      clearGroup("municipalities_labels") %>%
      clearGroup("transects") %>%
      clearGroup("transect_labels") %>%
      clearGroup("anchor") 
    
    
    # AOI
    if (!is.null(rv$aoi_ll)) {
      p <- p %>% addPolygons(
        data  = sf::st_transform(rv$aoi_ll, 4326),
        color = "black", weight = 2, fill = FALSE, group = "aoi"
      )
    }
    
    # BBox
    if (!is.null(rv$bb_final)) {
      p <- p %>% addRectangles(
        lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
        lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
        color = "magenta", weight = 2, group = "bbox"
      )
    }
    
    # Airports + labels (if available from Points tab)
    if (!is.null(rv$air_sf) && nrow(rv$air_sf) > 0) {
      air_ll <- sf::st_transform(rv$air_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = air_ll, radius = 4, color = "blue",
        label  = ~airport_label, group  = "airports"
      )
      p <- safe_add_label_only(
        p, data = air_ll, label = ~airport_label,
        group = "airports_labels", textsize = "11px", offset = c(10, -10)
      )
    }
    
    # Municipalities + labels
    if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0) {
      muni_ll <- sf::st_transform(rv$muni_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = muni_ll, radius = 3, color = "red",
        label  = ~as.character(name), group  = "municipalities"
      )
      p <- safe_add_label_only(
        p, data = muni_ll, label = ~as.character(name),
        group = "municipalities_labels", textsize = "10px", offset = c(8, -8)
      )
    }
    
    # Municipalities: points + always-visible labels
    if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0) {
      muni_ll <- sf::st_transform(rv$muni_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = muni_ll,
        radius = 3, color = "red",
        label  = ~as.character(name),
        group  = "municipalities"
      )
      p <- safe_add_label_only(
        p, data = muni_ll, label = ~as.character(name),
        group = "municipalities_labels", textsize = "10px", offset = c(8, -8)
      )
    }
    
    # Keep current bbox in view if available
    if (!is.null(bb_final_rx())) {
      bb <- bb_final_rx()
      p %>% fitBounds(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"])
    }
    
    # Transects + always-visible labels (robust midpoint with fallback)
    if (!is.null(rv$transects_sf) && nrow(rv$transects_sf) > 0) {
      tr_ll <- sf::st_transform(rv$transects_sf, 4326)
      p <- p %>% addPolylines(
        data  = tr_ll, color = "red", weight = 2,
        label = ~paste0("T", id), group = "transects"
      )
      
      midpts <- NULL
      # Try robust midpoint on merged lines
      midpts <- tryCatch({
        gm <- sf::st_line_merge(sf::st_geometry(tr_ll))
        sf::st_line_sample(gm, sample = 0.5) |> sf::st_cast("POINT")
      }, error = function(e) NULL)
      # Fallback to centroids if needed
      if (is.null(midpts) || length(midpts) != nrow(tr_ll)) {
        midpts <- tryCatch(sf::st_centroid(sf::st_geometry(tr_ll)), error = function(e) NULL)
      }
      
      if (!is.null(midpts)) {
        n <- min(length(midpts), nrow(tr_ll))
        lab_sf <- sf::st_as_sf(
          data.frame(id = tr_ll$id[seq_len(n)]),
          geometry = sf::st_sfc(midpts)[seq_len(n)], crs = 4326
        )
        p <- safe_add_label_only(
          p, data = lab_sf, label = ~paste0("T", id),
          group = "transect_labels", textsize = "11px", offset = c(8, -8)
        )
      }
    }
    
    # Draggable anchor marker (if enabled)
    if (isTRUE(input$enable_move)) {
      if (is.null(rv$anchor_lonlat) && !is.null(rv$bb_final)) {
        rv$anchor_lonlat <- c(
          as.numeric((rv$bb_final["xmin"] + rv$bb_final["xmax"]) / 2),
          as.numeric((rv$bb_final["ymin"] + rv$bb_final["ymax"]) / 2)
        )
      }
      if (!is.null(rv$anchor_lonlat)) {
        p <- p %>% addMarkers(
          lng = rv$anchor_lonlat[1], lat = rv$anchor_lonlat[2],
          layerId = "anchor", group = "anchor",
          popup = "Drag me to move transects",
          options = markerOptions(draggable = TRUE)
        )
      }
    }
    
    if (fit && !is.null(rv$bb_final)) {
      p %>% fitBounds(
        rv$bb_final["xmin"], rv$bb_final["ymin"],
        rv$bb_final["xmax"], rv$bb_final["ymax"]
      )
    }
  }
  
  # Fit & draw when bbox changes (Transects)
  observeEvent(rv$bb_final, {
    if (map_exists("map_transects")) {
      leafletProxy("map_transects", session = session) %>%
        fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                  rv$bb_final["xmax"], rv$bb_final["ymax"])
      try(redraw_transects_map(), silent = TRUE)
    }
  }, ignoreInit = TRUE)
  
  observeEvent(rv$transects_sf, { try(redraw_transects_map(), silent = TRUE) }, ignoreInit = TRUE)
  observeEvent(input$enable_move, { try(redraw_transects_map(), silent = TRUE) }, ignoreInit = TRUE)
  
  # Click to set anchor
  observeEvent(input$map_transects_click, {
    rv$anchor_lonlat <- c(input$map_transects_click$lng, input$map_transects_click$lat)
    try(redraw_transects_map(), silent = TRUE)
  })
  
  # Drag anchor: rebuild or translate
  observeEvent(input$map_transects_marker_dragend, {
    info <- input$map_transects_marker_dragend
    req(!is.null(info))
    if (!is.null(info$id) && info$id != "anchor") return()
    
    new_lonlat <- c(info$lng, info$lat)
    old_lonlat <- rv$anchor_lonlat
    rv$anchor_lonlat <- new_lonlat
    
    if (isTRUE(input$rebuild_on_drag)) {
      try({
        tr_ll <- make_transects_ll(
          anchor_lonlat = rv$anchor_lonlat,
          n = input$n_lines, len_km = input$len_km, sp_km = input$sp_km,
          brg_deg = input$brg_deg, clip_to_bbox = TRUE, bbox_ll = rv$bb_final
        )
        tr_sf <- sf::st_transform(tr_ll, terre_crs)

        if (!is.null(rv$aoi_terre)) {
          aoi_u <- sf::st_make_valid(sf::st_union(rv$aoi_terre))
          tr_sf <- tryCatch(sf::st_intersection(tr_sf, aoi_u), error = function(e) tr_sf)
          tr_sf <- tr_sf[!sf::st_is_empty(tr_sf), , drop = FALSE]
        }
        
        tr_sf <- clean_lines_for_optimizer(tr_sf)
        if (!("id" %in% names(tr_sf))) tr_sf$id <- seq_len(nrow(tr_sf))
        rv$transects_sf <- tr_sf
      }, silent = TRUE)
    } else {
      if (!is.null(old_lonlat) && !is.null(rv$transects_sf) && nrow(rv$transects_sf) > 0) {
        old_pt <- sf::st_sfc(sf::st_point(old_lonlat), crs = 4326) |> sf::st_transform(terre_crs)
        new_pt <- sf::st_sfc(sf::st_point(new_lonlat), crs = 4326) |> sf::st_transform(terre_crs)
        delta  <- sf::st_coordinates(new_pt) - sf::st_coordinates(old_pt)
        sf::st_geometry(rv$transects_sf) <- sf::st_geometry(rv$transects_sf) + c(delta[1], delta[2])
      }
    }
    try(redraw_transects_map(), silent = TRUE)
  }, ignoreInit = TRUE)
  
  # Reset anchor
  observeEvent(input$reset_anchor, {
    req(rv$bb_final)
    new_anchor <- c(
      as.numeric((rv$bb_final["xmin"] + rv$bb_final["xmax"]) / 2),
      as.numeric((rv$bb_final["ymin"] + rv$bb_final["ymax"]) / 2)
    )
    old <- rv$anchor_lonlat
    rv$anchor_lonlat <- new_anchor
    
    if (isTRUE(input$rebuild_on_drag)) {
      try({
        tr_ll <- make_transects_ll(
          anchor_lonlat = rv$anchor_lonlat,
          n = input$n_lines, len_km = input$len_km, sp_km = input$sp_km,
          brg_deg = input$brg_deg, clip_to_bbox = TRUE, bbox_ll = rv$bb_final
        )
        tr_sf <- sf::st_transform(tr_ll, terre_crs)
        tr_sf <- clean_lines_for_optimizer(tr_sf)
        if (!("id" %in% names(tr_sf))) tr_sf$id <- seq_len(nrow(tr_sf))
        rv$transects_sf <- tr_sf
      }, silent = TRUE)
    } else if (!is.null(old) && !is.null(rv$transects_sf) && nrow(rv$transects_sf) > 0) {
      old_pt <- sf::st_sfc(sf::st_point(old), crs = 4326) |> sf::st_transform(terre_crs)
      new_pt <- sf::st_sfc(sf::st_point(new_anchor), crs = 4326) |> sf::st_transform(terre_crs)
      delta  <- sf::st_coordinates(new_pt) - sf::st_coordinates(old_pt)
      sf::st_geometry(rv$transects_sf) <- sf::st_geometry(rv$transects_sf) + c(delta[1], delta[2])
    }
    
    try(redraw_transects_map(fit = TRUE), silent = TRUE)
  }, ignoreInit = TRUE)
  
  # Load existing transects
  observeEvent(input$tr_file, {
    tr_raw <- read_vector_upload(input$tr_file)
    tr_raw <- sf::st_make_valid(tr_raw) |> sf::st_transform(terre_crs)
    tr_sf  <- clean_lines_for_optimizer(tr_raw)
    if (!("id" %in% names(tr_sf))) tr_sf$id <- seq_len(nrow(tr_sf))
    rv$transects_sf <- tr_sf
    try(redraw_transects_map(), silent = TRUE)
  })
  
  # Build transects from controls
  observeEvent(input$build_tr, {
    req(rv$bb_final)
    withProgress(message = "Building transects…", value = 0, {
      incProgress(0.25, detail = "Computing lines")
      out <- tryCatch({
        if (is.null(rv$anchor_lonlat)) {
          rv$anchor_lonlat <- c(
            as.numeric((rv$bb_final["xmin"] + rv$bb_final["xmax"]) / 2),
            as.numeric((rv$bb_final["ymin"] + rv$bb_final["ymax"]) / 2)
          )
        }
        tr_ll <- make_transects_ll(
          anchor_lonlat = rv$anchor_lonlat,
          n = input$n_lines, len_km = input$len_km, sp_km = input$sp_km,
          brg_deg = input$brg_deg, clip_to_bbox = TRUE, bbox_ll = rv$bb_final
        )
        incProgress(0.6, detail = "Transform & clean")
        tr_sf <- sf::st_transform(tr_ll, terre_crs)

        # --- Clip to AOI if available ---
        if (!is.null(rv$aoi_terre)) {
          # union AOI to ensure a single polygon; make valid to avoid topology errors
          aoi_u <- sf::st_make_valid(sf::st_union(rv$aoi_terre))
          tr_sf <- tryCatch(
            sf::st_intersection(tr_sf, aoi_u),
            error = function(e) {
              showNotification(paste("Transect clipping to AOI failed:", conditionMessage(e)),
                               type = "warning", duration = 8)
              tr_sf
            }
          )
          # drop empty geometries if any
          tr_sf <- tr_sf[!sf::st_is_empty(tr_sf), , drop = FALSE]
        }
        
        tr_sf <- clean_lines_for_optimizer(tr_sf)
        if (!("id" %in% names(tr_sf))) tr_sf$id <- seq_len(nrow(tr_sf))
        tr_sf
      }, error = function(e) {
        showNotification(paste("Failed to build transects:", conditionMessage(e)),
                         type = "error", duration = 8)
        return(NULL)
      })
      if (is.null(out)) return(invisible(NULL))
      rv$transects_sf <- out
      incProgress(1, detail = "Drawing")
      try(redraw_transects_map(), silent = TRUE)
      showNotification(sprintf("Built %d transects.", nrow(rv$transects_sf)),
                       type = "message", duration = 4)
    })
  })
  
  # Subset remaining transects interactively
  output$transects_table <- renderDT({
    req(rv$transects_sf)
    datatable(data.frame(id = rv$transects_sf$id), selection = "multiple",
              options = list(pageLength = 10))
  })
  observeEvent(input$transects_table_rows_selected, {
    req(rv$transects_sf)
    sel <- input$transects_table_rows_selected
    if (length(sel) > 0) {
      keep <- setdiff(seq_len(nrow(rv$transects_sf)), sel)
      rv$transects_sf <- rv$transects_sf[keep, , drop = FALSE]
      try(redraw_transects_map(), silent = TRUE)
    }
  })
  
  # Save transects
  output$dl_transects_gpkg <- downloadHandler(
    filename = function() sprintf("%s_transects.gpkg", input$out_basename),
    content = function(file) {
      req(rv$transects_sf)
      sf::st_write(rv$transects_sf, dsn = file, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
    }
  )
  output$dl_transects_shp <- downloadHandler(
    filename = function() sprintf("%s_transects_shp.zip", input$out_basename),
    content = function(file) {
      req(rv$transects_sf)
      tmpdir <- tempfile("shp"); dir.create(tmpdir)
      layer  <- "transects"
      suppressWarnings(sf::st_write(rv$transects_sf, dsn = tmpdir, layer = layer,
                                    driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE))
      zip(zipfile = file, files = list.files(tmpdir, full.names = TRUE), mode = "cherry-pick")
    }
  )
  
  # Reset transects
  observeEvent(input$reset_transects, {
    rv$transects_sf <- NULL
    rv$anchor_lonlat <- NULL
    redraw_transects_map()
    showNotification("Transects reset.", type = "message", duration = 3)
  })
  
  # --- OPTIMIZER TAB (Leaflet-based) ------------------------------------------
  output$map_opt <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(if (is.null(input$basemap)) BASEMAP_DEFAULT else input$basemap) |>
      addScaleBar(position = "bottomleft") |>
      addLayersControl(
        baseGroups = NULL,
        overlayGroups = c("aoi", "bbox", "airports", "municipalities",
                          "airports_labels", "municipalities_labels",
                          "transects", "transect_labels", "routes"),
        options = layersControlOptions(collapsed = FALSE)
      ) |>
      setView(lng = -64.5, lat = 48.8, zoom = 6)
  })
  outputOptions(output, "map_opt", suspendWhenHidden = FALSE)
  
  # Redraw helper for optimizer map (routes + everything else)
  redraw_opt_map <- function(fit = FALSE) {
    if (!map_exists("map_opt")) return(invisible(NULL))
    
    p <- leafletProxy("map_opt", session = session) %>%
      clearGroup("aoi") %>%
      clearGroup("bbox") %>%
      clearGroup("airports") %>%
      clearGroup("municipalities") %>%
      clearGroup("airports_labels") %>%
      clearGroup("municipalities_labels") %>%
      clearGroup("transects") %>%
      clearGroup("transect_labels") %>%
      clearGroup("routes")
    
    # AOI
    if (!is.null(rv$aoi_ll)) {
      p <- p %>% addPolygons(
        data  = sf::st_transform(rv$aoi_ll, 4326),
        color = "black", weight = 2, fill = FALSE, group = "aoi"
      )
    }
    
    # BBox
    if (!is.null(rv$bb_final)) {
      p <- p %>% addRectangles(
        lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
        lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
        color = "magenta", weight = 2, group = "bbox"
      )
    }
    
    # Airports + labels
    if (!is.null(rv$air_sf) && nrow(rv$air_sf) > 0) {
      air_ll <- sf::st_transform(rv$air_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = air_ll,
        radius = 4, color = "blue", label = ~airport_label,
        group  = "airports"
      )
      p <- safe_add_label_only(
        p, data = air_ll, label = ~airport_label,
        group = "airports_labels", textsize = "11px", offset = c(10, -10)
      )
    }
    
    # Municipalities + labels
    if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0) {
      muni_ll <- sf::st_transform(rv$muni_sf, 4326)
      p <- p %>% addCircleMarkers(
        data   = muni_ll, radius = 3, color = "red", label = ~as.character(name),
        group  = "municipalities"
      )
      p <- safe_add_label_only(
        p, data = muni_ll, label = ~as.character(name),
        group = "municipalities_labels", textsize = "10px", offset = c(8, -8)
      )
    }
    
    # Transects + labels
    if (!is.null(rv$transects_sf) && nrow(rv$transects_sf) > 0) {
      tr_ll <- sf::st_transform(rv$transects_sf, 4326)
      p <- p %>% addPolylines(
        data  = tr_ll, color = "red", weight = 2,
        label = ~paste0("T", id), group = "transects"
      )
      midpts <- NULL
      midpts <- tryCatch({
        gm <- sf::st_line_merge(sf::st_geometry(tr_ll))
        sf::st_line_sample(gm, sample = 0.5) |> sf::st_cast("POINT")
      }, error = function(e) NULL)
      if (is.null(midpts) || length(midpts) != nrow(tr_ll)) {
        midpts <- tryCatch(sf::st_centroid(sf::st_geometry(tr_ll)), error = function(e) NULL)
      }
      if (!is.null(midpts)) {
        n <- min(length(midpts), nrow(tr_ll))
        lab_sf <- sf::st_as_sf(data.frame(id = tr_ll$id[seq_len(n)]),
                               geometry = sf::st_sfc(midpts)[seq_len(n)], crs = 4326)
        p <- safe_add_label_only(
          p, data = lab_sf, label = ~paste0("T", id),
          group = "transect_labels", textsize = "11px", offset = c(8, -8)
        )
      }
    }
    
    # Routes (once optimizer finishes) — DISTINCT COLORS PER TRIP
    if (!is.null(rv$routes_sf) && nrow(rv$routes_sf) > 0) {
      routes_sf <- rv$routes_sf
      if (is.na(sf::st_crs(routes_sf))) {
        sf::st_crs(routes_sf) <- terre_crs
      }
      routes_ll <- sf::st_transform(routes_sf, 4326)
      
      # Trip grouping
      if ("trip_id" %in% names(routes_ll)) {
        routes_ll$trip_group <- as.character(routes_ll$trip_id)
      } else if ("route_id" %in% names(routes_ll)) {
        routes_ll$trip_group <- as.character(routes_ll$route_id)
      } else if ("id" %in% names(routes_ll)) {
        routes_ll$trip_group <- as.character(routes_ll$id)
      } else {
        routes_ll$trip_group <- as.character(seq_len(nrow(routes_ll)))
      }
      
      routes_ll$label_txt  <- routes_ll$trip_group
      routes_ll$trip_group <- factor(routes_ll$trip_group)
      
      pal <- leaflet::colorFactor(
        palette = "Set2",
        domain  = levels(routes_ll$trip_group),
        na.color = "#888888"
      )
      
      p <- p %>% addPolylines(
        data   = routes_ll,
        color  = ~pal(trip_group),
        weight = 3, opacity = 0.9,
        label  = ~label_txt,
        group  = "routes"
      ) %>% showGroup("routes")
      
      # Legend (remove older one if supported, then add)
      try(p <- p %>% removeControl("routes_legend"), silent = TRUE)
      p <- p %>% addLegend(
        position = "topright",
        pal      = pal,
        values   = routes_ll$trip_group,
        title    = "Trip",
        opacity  = 0.9,
        layerId  = "routes_legend"
      )
    }
    
    if (fit && !is.null(rv$bb_final)) {
      p %>% fitBounds(
        rv$bb_final["xmin"], rv$bb_final["ymin"],
        rv$bb_final["xmax"], rv$bb_final["ymax"]
      )
    }
  }
  
  # Fit & draw when bbox changes (Optimizer map)
  observeEvent(rv$bb_final, {
    if (map_exists("map_opt")) {
      leafletProxy("map_opt", session = session) %>%
        fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                  rv$bb_final["xmax"], rv$bb_final["ymax"])
      try(redraw_opt_map(), silent = TRUE)
    }
  }, ignoreInit = TRUE)
  
  # Redraw optimizer map when inputs change
  observeEvent(list(rv$aoi_ll, rv$air_sf, rv$muni_sf, rv$transects_sf, rv$routes_sf), {
    try(redraw_opt_map(), silent = TRUE)
  }, ignoreInit = TRUE)
  
  # Also (re)draw when opening the Optimize tab
  observeEvent(input$main_tabs, {
    if (identical(input$main_tabs, "Optimize") && map_exists("map_opt")) {
      try(redraw_opt_map(), silent = TRUE)
    }
  }, ignoreInit = TRUE)
  
  # Run optimizer
  observeEvent(input$run_opt, {
    req(rv$transects_sf)
    
    withProgress(message = "Optimizing routes…", value = 0, {
      showNotification("Optimizer started…", type = "message", duration = 3)
      incProgress(0.15, detail = "Preparing data")
      
      air_all_for_opt <- sf::st_transform(load_airports(AIRPORTS_CSV), terre_crs)
      
      multi_airports <- NULL; multi_k <- 2
      mi <- trimws(input$multi_mode)
      if (nzchar(mi)) {
        if (grepl("^auto\\s+\\d+$", mi, ignore.case = TRUE)) {
          multi_k <- as.integer(sub(".*\\s+", "", mi)); multi_airports <- "auto"
        } else {
          multi_airports <- strsplit(mi, "\\s*,\\s*")[[1]]
        }
      }
      
      incProgress(0.45, detail = "Running solver")
      res <- tryCatch({
        run_transect_optimizer(
          transects_sf = rv$transects_sf,
          air_sf       = air_all_for_opt,
          terre_crs    = terre_crs,
          airport_code = if (nzchar(input$start_ap)) input$start_ap else NULL,
          end_airport_code = if (nzchar(input$end_ap)) input$end_ap else NULL,
          multi_airports   = multi_airports,
          multi_k          = multi_k,
          speed_knots  = input$speed_kn,
          max_hours    = input$max_hrs,
          export_dir   = TRANSECT_EXPORT_DIR,
          out_basename = input$out_basename,
          make_plot    = FALSE,        # we render with Leaflet now
          terre_crop   = rv$terre_crop,
          xlim_terre   = rv$xlim_terre,
          ylim_terre   = rv$ylim_terre,
          air_sf_plot  = rv$air_sf,
          muni_sf_plot = rv$muni_sf,
          aoi_sf_plot  = rv$aoi_terre
        )
      }, error = function(e) {
        showNotification(paste("Optimizer failed:", conditionMessage(e)),
                         type = "error", duration = 8)
        return(NULL)
      })
      
      if (is.null(res)) return(invisible(NULL))
      incProgress(0.85, detail = "Finalizing outputs")
      
      rv$routes_sf    <- res$routes_sf
      rv$trip_summary <- res$trip_summary
      rv$manifest_df  <- res$flight_manifest
      
      if (is.null(rv$routes_sf) || nrow(rv$routes_sf) == 0) {
        incProgress(1, detail = "No routes produced")
        showNotification("Optimizer finished, but produced no routes for the given settings.",
                         type = "warning", duration = 8)
      } else {
        try(redraw_opt_map(fit = TRUE), silent = TRUE)
        incProgress(1, detail = "Routes drawn")
        showNotification(
          sprintf("Optimizer finished. %d route segment(s) drawn.", nrow(rv$routes_sf)),
          type = "message", duration = 6
        )
      }
    })
  })
  
  output$trips_table <- renderDT({
    req(rv$trip_summary); datatable(rv$trip_summary, options = list(pageLength = 5))
  })
  output$manifest_table <- renderDT({
    req(rv$manifest_df); datatable(rv$manifest_df, options = list(pageLength = 5))
  })
  
  # Downloads (routes)
  output$dl_routes_gpkg <- downloadHandler(
    filename = function() sprintf("%s_flight_routes.gpkg", input$out_basename),
    content = function(file) {
      req(rv$routes_sf)
      sf::st_write(rv$routes_sf, dsn = file, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
    }
  )
  output$dl_routes_shp <- downloadHandler(
    filename = function() sprintf("%s_flight_routes_shp.zip", input$out_basename),
    content = function(file) {
      req(rv$routes_sf)
      tmpdir <- tempfile("shp"); dir.create(tmpdir)
      layer  <- "routes"
      suppressWarnings(sf::st_write(rv$routes_sf, dsn = tmpdir, layer = layer,
                                    driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE))
      zip(zipfile = file, files = list.files(tmpdir, full.names = TRUE), mode = "cherry-pick")
    }
  )
  output$dl_trips_csv <- downloadHandler(
    filename = function() sprintf("%s_flight_trips.csv", input$out_basename),
    content = function(file) {
      req(rv$trip_summary); write.csv(rv$trip_summary, file, row.names = FALSE, na = "")
    }
  )
  output$dl_manifest_csv <- downloadHandler(
    filename = function() sprintf("%s_flight_manifest.csv", input$out_basename),
    content = function(file) {
      req(rv$manifest_df); write.csv(rv$manifest_df, file, row.names = FALSE, na = "")
    }
  )
  
  # Reset optimizer outputs
  observeEvent(input$reset_opt, {
    rv$routes_sf <- rv$trip_summary <- rv$manifest_df <- NULL
    redraw_opt_map()
    showNotification("Optimizer outputs reset.", type = "message", duration = 3)
  })
}

shinyApp(ui, server)