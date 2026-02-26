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

# --- Runtime options (no debugger halt unless explicitly enabled) -------------
if (interactive() && isTRUE(getOption("shiny.debug", FALSE))) {
  options(shiny.fullstacktrace = TRUE, shiny.error = browser)
} else {
  options(shiny.error = NULL)
}

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

read_vector_upload <- function(input_file) {
  req(input_file)
  f <- input_file$datapath
  ext <- tolower(tools::file_ext(input_file$name))
  if (ext == "gpkg") {
    return(suppressWarnings(sf::st_read(f, quiet = TRUE)))
  } else if (ext == "zip") {
    td <- tempfile("unz"); dir.create(td)
    unzip(f, exdir = td)
    shp <- list.files(td, pattern = "\\.shp$", full.names = TRUE)
    validate(need(length(shp) >= 1, "No .shp found in zip"))
    return(suppressWarnings(sf::st_read(shp[1], quiet = TRUE)))
  } else if (ext %in% c("shp","dbf","shx","prj","cpg")) {
    return(suppressWarnings(sf::st_read(f, quiet = TRUE)))
  } else {
    validate(need(FALSE, "Unsupported file: upload .gpkg or a zipped shapefile (.zip)"))
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
  title = "Pelagic Transects Planner (Shiny)",
  id = "main_tabs",
  tabPanel("Area & BBox",
           fluidRow(
             column(4,
                    h4("1) Choose Area"),
                    # Basemap selector
                    selectInput("basemap", "Basemap style", choices = TOPO_CHOICES, selected = BASEMAP_DEFAULT),
                    
                    radioButtons("area_mode", "Mode", choices = c("AOI shapefile" = "aoi", "Landmark search" = "lm"), inline = TRUE),
                    conditionalPanel("input.area_mode == 'aoi'",
                                     fileInput("aoi_file", "Upload AOI (.gpkg or zipped shapefile)", accept = c(".zip",".gpkg",".shp",".dbf",".shx",".prj",".cpg")),
                                     helpText("If shapefile, zip .shp/.shx/.dbf/.prj/.cpg and upload the .zip.")
                    ),
                    conditionalPanel("input.area_mode == 'lm'",
                                     textInput("lm_query", "Landmark search (e.g., Perce, Gaspe)", "Perce"),
                                     textInput("lm_province", "Province code (e.g., 24=QC, blank=Any)", "24"),
                                     sliderInput("lm_buf_km", "Initial bbox buffer (km)", min = 5, max = 100, value = 25, step = 5),
                                     actionButton("lm_search", "Search Landmark")
                    ),
                    hr(),
                    h4("2) Refine BBox"),
                    prettySwitch("draw_rect_mode", "Draw rectangle", value = FALSE, status = "primary"),
                    actionButton("use_view_bbox", "Use current map view as BBox", icon = icon("crop")),
                    helpText("Pan/zoom or draw a rectangle → click 'Use current map view as BBox'."),
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
                    h4("Filtering"),
                    prettySwitch("limit_to_aoi_buf", "Limit to AOI + buffer", value = TRUE, status = "primary"),
                    sliderInput("buf_km", "Buffer distance (km) — 0 = auto-scale", min = 0, max = 100, value = 0, step = 5),
                    actionButton("refresh_points", "Apply Filter", icon = icon("filter")),
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
                    h4("Use existing transects (optional)"),
                    fileInput("tr_file", "Upload transects (.gpkg or zipped shapefile)", accept = c(".zip",".gpkg",".shp",".dbf",".shx",".prj",".cpg")),
                    hr(),
                    
                    h4("Or build transects"),
                    numericInput("n_lines", "Number of lines", 10, min = 1, step = 1),
                    numericInput("len_km", "Length (km)", 10, min = 0.1, step = 0.5),
                    numericInput("sp_km", "Spacing (km)", 1, min = 0.1, step = 0.1),
                    numericInput("brg_deg", "Bearing (deg) 0=N, 90=E", 90, min = 0, max = 359, step = 1),
                    helpText("Anchor defaults to bbox center. Click on map to set anchor."),
                    actionButton("build_tr", "Build/Refresh Transects", icon = icon("sliders")),
                    hr(),
                    
                    # Move transects (above subset)
                    h4("Move transects"),
                    prettySwitch("enable_move", "Enable drag-to-move (anchor)", value = FALSE, status = "info"),
                    checkboxInput("rebuild_on_drag", "Rebuild lines when dragging (instead of shift)", value = FALSE),
                    actionButton("reset_anchor", "Reset anchor to bbox center", icon = icon("crosshairs")),
                    hr(),
                    
                    # Save transects
                    h4("Save transects"),
                    downloadButton("dl_transects_gpkg", "Download Transects (GPKG)"),
                    downloadButton("dl_transects_shp",  "Download Transects (Shapefile ZIP)"),
                    hr(),
                    
                    h4("Subset remaining transects"),
                    DTOutput("transects_table")
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
                    downloadButton("dl_manifest_csv","Download Manifest CSV")
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
  
  # --- Basemap switcher (safe; only updates maps that already exist) ----------
  observeEvent(input$basemap, ignoreInit = TRUE, {
    req(input$basemap)
    if (!(input$basemap %in% names(leaflet::providers))) return()
    
    map_exists <- function(id) {
      !is.null(input[[paste0(id, "_zoom")]]) || !is.null(input[[paste0(id, "_bounds")]])
    }
    
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
  observeEvent(input$aoi_file, ignoreInit = TRUE, {
    aoi_raw <- read_vector_upload(input$aoi_file)
    aoi_raw <- sf::st_make_valid(aoi_raw)
    aoi_poly <- suppressWarnings(sf::st_collection_extract(aoi_raw, "POLYGON"))
    validate(need(nrow(aoi_poly) > 0, "AOI has no polygon geometries"))
    
    rv$aoi_ll    <- sf::st_transform(aoi_poly, 4269)
    rv$aoi_terre <- sf::st_transform(rv$aoi_ll, terre_crs)
    
    rv$bb_final <- sf::st_bbox(rv$aoi_ll)
    lims <- compute_limits_from_bbox(rv$bb_final, terre, terre_crs)
    rv$bb_terre   <- lims$bb_terre
    rv$terre_crop <- lims$terre_crop
    rv$xlim_terre <- lims$xlim_terre
    rv$ylim_terre <- lims$ylim_terre
    rv$grat_final <- sf::st_graticule(bbox = rv$bb_final, crs = sf::st_crs(4269), ndx = 8, ndy = 8) |> sf::st_sf()
    
    leafletProxy("map_bbox", session = session) %>%
      clearShapes() %>%
      addPolygons(data = sf::st_transform(rv$aoi_ll, 4326), color = "black", weight = 2, fill = FALSE) %>%
      addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                    lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                    color = "magenta", weight = 2) %>%
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"], rv$bb_final["xmax"], rv$bb_final["ymax"])
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
    cols <- intersect(names(rv$lm_results), c("name","concise","province","latitude","longitude","id","key"))
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
  
  # --- POINTS TAB --------------------------------------------------------------
  # Render the map immediately; add layers control; we will draw layers via proxy
  output$map_points <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(if (is.null(input$basemap)) BASEMAP_DEFAULT else input$basemap) |>
      addScaleBar(position = "bottomleft") |>
      addLayersControl(
        baseGroups = NULL,
        overlayGroups = c("aoi", "bbox", "airports", "municipalities",
                          "airports_labels", "municipalities_labels"),
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
  
  # Helper to (re)draw points via proxy — defined INSIDE server so it sees session
  redraw_points_map <- function() {
    # Only update if the widget exists client-side
    if (is.null(input$map_points_bounds) && is.null(input$map_points_zoom)) return(invisible(NULL))
    
    p <- leafletProxy("map_points", session = session) %>%
      clearGroup("aoi") %>%
      clearGroup("bbox") %>%
      clearGroup("airports") %>%
      clearGroup("municipalities") %>%
      clearGroup("airports_labels") %>%
      clearGroup("municipalities_labels")
    
    # AOI outline
    if (!is.null(rv$aoi_ll)) {
      p <- p %>% addPolygons(
        data  = sf::st_transform(rv$aoi_ll, 4326),
        color = "black", weight = 2, fill = FALSE, group = "aoi"
      )
    }
    
    # BBox rectangle
    if (!is.null(rv$bb_final)) {
      p <- p %>% addRectangles(
        lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
        lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
        color = "magenta", weight = 2, group = "bbox"
      )
    }
    
    # Airports: points + always-visible labels
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
    
    # Fit if bbox available
    if (!is.null(rv$bb_final)) {
      p %>% fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                      rv$bb_final["xmax"], rv$bb_final["ymax"])
    }
  }
  
  # When BBox becomes available/changes, fit & draw on the Points map
  observeEvent(rv$bb_final, {
    if (is.null(input$map_points_bounds) && is.null(input$map_points_zoom)) return(invisible(NULL))
    leafletProxy("map_points", session = session) %>%
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                rv$bb_final["xmax"], rv$bb_final["ymax"])
    try(redraw_points_map(), silent = TRUE)
  }, ignoreInit = TRUE)
  
  # Apply Filter (robust)
  observeEvent(input$refresh_points, ignoreInit = TRUE, {
    req(rv$bb_final)
    
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
      
      buf_m <- if (input$buf_km > 0) input$buf_km * 1000 else scaled_buffer_m(rv$bb_final)
      
      incProgress(0.35, detail = "Filtering airports")
      rv$air_sf <- tryCatch(
        filter_points_after_bbox(
          rv$air_all, rv$bb_final, terre_crs,
          aoi_terre   = rv$aoi_terre,
          use_aoi_buf = isTRUE(input$limit_to_aoi_buf),
          buf_m       = buf_m
        ),
        error = function(e) {
          showNotification(paste("Airport filtering failed:", conditionMessage(e)),
                           type = "error", duration = 8)
          return(sf::st_sf())
        }
      )
      if (nrow(rv$air_sf) > 0) {
        rv$air_sf <- tryCatch(
          make_airport_labels(rv$air_sf),
          error = function(e) {
            showNotification(paste("Could not create airport labels:", conditionMessage(e)),
                             type = "warning", duration = 6)
            rv$air_sf
          }
        )
      }
      
      incProgress(0.65, detail = "Fetching municipalities")
      muni_raw <- tryCatch(
        geonames_bbox_paged(rv$bb_final, province = NULL, concise = "CITY"),
        error = function(e) {
          showNotification(paste("Municipality lookup failed:", conditionMessage(e)),
                           type = "error", duration = 8)
          return(NULL)
        }
      )
      
      incProgress(0.8, detail = "Preparing municipality points")
      if (!is.null(muni_raw)) {
        muni_ll <- tryCatch(
          to_points_safe(muni_raw, rv$bb_final) |> bbox_filter_points(rv$bb_final),
          error = function(e) {
            showNotification(paste("Municipality point conversion failed:", conditionMessage(e)),
                             type = "warning", duration = 6)
            return(NULL)
          }
        )
        if (!is.null(muni_ll) && is.na(sf::st_crs(muni_ll))) {
          muni_ll <- sf::st_set_crs(muni_ll, 4326)
        }
        rv$muni_sf <- if (!is.null(muni_ll)) sf::st_transform(muni_ll, terre_crs) else NULL
      } else {
        rv$muni_sf <- NULL
      }
      
      redraw_points_map()
      
      n_air <- if (!is.null(rv$air_sf)) nrow(rv$air_sf) else 0
      n_mun <- if (!is.null(rv$muni_sf)) nrow(rv$muni_sf) else 0
      showNotification(sprintf("Filtered: %d airports, %d municipalities.", n_air, n_mun),
                       type = "message", duration = 4)
    })
  })
  
  output$pts_counts <- renderText({
    a <- if (!is.null(rv$air_sf)) nrow(rv$air_sf) else 0
    m <- if (!is.null(rv$muni_sf)) nrow(rv$muni_sf) else 0
    if (a + m == 0) "Click 'Apply Filter'" else sprintf("Airports: %d | Municipalities: %d", a, m)
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
    if (is.null(input$map_transects_bounds) && is.null(input$map_transects_zoom)) return(invisible(NULL))
    
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
    
    # Airports + labels (reuse from points tab state if available)
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
        data   = muni_ll,
        radius = 3, color = "red", label = ~as.character(name),
        group  = "municipalities"
      )
      p <- safe_add_label_only(
        p, data = muni_ll, label = ~as.character(name),
        group = "municipalities_labels", textsize = "10px", offset = c(8, -8)
      )
    }
    
    # Transects + labels at midpoints
    if (!is.null(rv$transects_sf) && nrow(rv$transects_sf) > 0) {
      tr_ll <- sf::st_transform(rv$transects_sf, 4326)
      p <- p %>% addPolylines(
        data  = tr_ll, color = "red", weight = 2,
        label = ~paste0("T", id), group = "transects"
      )
      midpts <- tryCatch({
        sf::st_line_sample(sf::st_geometry(tr_ll), sample = 0.5) |> sf::st_cast("POINT")
      }, error = function(e) NULL)
      
      if (!is.null(midpts)) {
        ids <- rv$transects_sf$id
        lab_sf <- sf::st_as_sf(data.frame(id = ids), geometry = midpts, crs = 4326)
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
  
  # Fit & draw when bbox changes
  observeEvent(rv$bb_final, {
    if (is.null(input$map_transects_bounds) && is.null(input$map_transects_zoom)) return(invisible(NULL))
    leafletProxy("map_transects", session = session) %>%
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                rv$bb_final["xmax"], rv$bb_final["ymax"])
    try(redraw_transects_map(), silent = TRUE)
  }, ignoreInit = TRUE)
  
  # Redraw when transects or move toggle change
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
    if (is.null(input$map_opt_bounds) && is.null(input$map_opt_zoom)) return(invisible(NULL))
    
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
        data   = muni_ll,
        radius = 3, color = "red", label = ~as.character(name),
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
      midpts <- tryCatch({
        sf::st_line_sample(sf::st_geometry(tr_ll), sample = 0.5) |> sf::st_cast("POINT")
      }, error = function(e) NULL)
      if (!is.null(midpts)) {
        ids <- rv$transects_sf$id
        lab_sf <- sf::st_as_sf(data.frame(id = ids), geometry = midpts, crs = 4326)
        p <- safe_add_label_only(
          p, data = lab_sf, label = ~paste0("T", id),
          group = "transect_labels", textsize = "11px", offset = c(8, -8)
        )
      }
    }
    
    # Routes (once optimizer finishes) — DISTINCT COLORS PER TRIP
    if (!is.null(rv$routes_sf) && nrow(rv$routes_sf) > 0) {
      routes_sf <- rv$routes_sf
      
      # Ensure we have a CRS before transforming for Leaflet
      if (is.na(sf::st_crs(routes_sf))) {
        sf::st_crs(routes_sf) <- terre_crs
      }
      routes_ll <- sf::st_transform(routes_sf, 4326)
      
      # Build a robust trip grouping and label column
      if ("trip_id" %in% names(routes_ll)) {
        routes_ll$trip_group <- as.character(routes_ll$trip_id)
      } else if ("route_id" %in% names(routes_ll)) {
        routes_ll$trip_group <- as.character(routes_ll$route_id)
      } else if ("id" %in% names(routes_ll)) {
        routes_ll$trip_group <- as.character(routes_ll$id)
      } else {
        routes_ll$trip_group <- as.character(seq_len(nrow(routes_ll)))
      }
      
      # Label text (kept simple & robust)
      routes_ll$label_txt <- routes_ll$trip_group
      
      # Make it a factor so the palette has stable levels
      routes_ll$trip_group <- factor(routes_ll$trip_group)
      
      # A pleasant qualitative palette for categories
      pal <- leaflet::colorFactor(
        palette = "Set2",                    # try "Dark2" or "Set1" if you prefer
        domain  = levels(routes_ll$trip_group),
        na.color = "#888888"
      )
      
      # Draw polylines colored by trip
      p <- p %>% addPolylines(
        data   = routes_ll,
        color  = ~pal(trip_group),
        weight = 3, opacity = 0.9,
        label  = ~label_txt,
        group  = "routes"
      ) %>% showGroup("routes")
      
      # Update the legend (remove an older one if present, then add a fresh one)
      # If your leaflet version supports layerId on legends:
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
    if (is.null(input$map_opt_bounds) && is.null(input$map_opt_zoom)) return(invisible(NULL))
    leafletProxy("map_opt", session = session) %>%
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"],
                rv$bb_final["xmax"], rv$bb_final["ymax"])
    try(redraw_opt_map(), silent = TRUE)
  }, ignoreInit = TRUE)
  
  # Redraw optimizer map when inputs change
  observeEvent(list(rv$aoi_ll, rv$air_sf, rv$muni_sf, rv$transects_sf, rv$routes_sf), {
    try(redraw_opt_map(), silent = TRUE)
  }, ignoreInit = TRUE)
  
  observeEvent(input$main_tabs, {
    if (identical(input$main_tabs, "Optimize")) {
      # Only act if the widget exists
      if (!is.null(input$map_opt_bounds) || !is.null(input$map_opt_zoom)) {
        try(redraw_opt_map(), silent = TRUE)
      }
    }
  }, ignoreInit = TRUE)
  
  # Run optimizer
  `%||%` <- function(a, b) if (!is.null(a)) a else b
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
          make_plot    = FALSE,        # we now render with Leaflet, not static PNG
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
}

shinyApp(ui, server)