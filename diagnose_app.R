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

options(shiny.fullstacktrace = TRUE, shiny.error = browser)

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
    source(path, encoding = "UTF-8", local = TRUE),
    error = function(e) {
      stop("Error while sourcing '", path, "': ", conditionMessage(e))
    }
  )
}

# Use YOUR module names and files
safe_source("modules/01_utils_bbox.R")
safe_source("modules/02_terre.R")
safe_source("modules/03_geonames.R")
safe_source("modules/04_bbox_selector.R")
safe_source("modules/05_airports.R")
safe_source("modules/07_plot.R")
safe_source("modules/08_transect_optimizer.R")
safe_source("TransectDesignerModule.R")

# Transect builder (uses make_transects_ll(); not refine_transects())
source_utf8("TransectDesignerModule.R")

# If your repo defines load_terre() in another module, source it.
# Otherwise, uncomment a minimal fallback here:
# if (!exists("load_terre")) {
#   stop("Missing load_terre(). Please source the module that defines load_terre(TERRE_PATH).")
#   # # Minimal fallback (uncomment if needed and you know your terre layer name):
#   # load_terre <- function(path) {
#   #   x <- suppressWarnings(sf::st_read(path, quiet = TRUE))
#   #   list(terre = x, crs = sf::st_crs(x))
#   # }
}
TERRE_PATH <- "data_downloads/background/TerreNAD83.shp"
terre_obj <- load_terre(TERRE_PATH)
terre     <- terre_obj$terre
terre_crs <- terre_obj$crs   # keep NAD83 display in degrees

# --- Helpers -------------------------------------------------------------------
compute_limits_from_bbox <- function(bb_final, terre, terre_crs) {
  bb_terre   <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb_final), terre_crs))
  terre_crop <- sf::st_crop(terre, bb_terre)
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

# --- UI ------------------------------------------------------------------------
ui <- navbarPage(
  title = "Pelagic Transects Planner (Shiny)",
  id = "main_tabs",
  tabPanel("Area & BBox",
           fluidRow(
             column(4,
                    h4("1) Choose Area"),
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
                    downloadButton("dl_routes_gpkg", "Download GPKG"),
                    downloadButton("dl_routes_shp",  "Download Shapefile (ZIP)"),
                    downloadButton("dl_trips_csv",   "Download Trips CSV"),
                    downloadButton("dl_manifest_csv","Download Manifest CSV"),
                    downloadButton("dl_routes_png",  "Download Routes PNG")
             ),
             column(8,
                    plotOutput("map_routes_png", height = 520) %>% withSpinner(color="#2C3E50"),
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
  
  # Base bbox map
  output$map_bbox <- renderLeaflet({
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(providers$Esri.WorldLightGrayBase) |>
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
    
    leafletProxy("map_bbox") %>%
      clearShapes() %>%
      addPolygons(data = sf::st_transform(rv$aoi_ll, 4326), color = "black", weight = 2, fill = FALSE) %>%
      addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                    lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                    color = "magenta", weight = 2) %>%
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"], rv$bb_final["xmax"], rv$bb_final["ymax"])
  })
  
  # Landmark search + table picker (no console)
  observeEvent(input$lm_search, {
    prov <- trimws(input$lm_province); if (prov == "") prov <- NULL
    cand <- geonames_search(input$lm_query, province = prov, num = 20)
    validate(need(nrow(cand) > 0, "No landmark results"))
    
    rv$lm_results <- cand
    showModal(modalDialog(
      title = "Select a landmark",
      size = "l",
      DTOutput("lm_table"),
      footer = tagList(actionButton("lm_pick_ok", "Use selected"), modalButton("Cancel"))
    ))
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
    validate(need(length(sel) == 1, "Select one landmark"))
    removeModal()
    
    pt <- rv$lm_results[sel, , drop = FALSE]
    rv$bb_final <- bbox_from_point_km(pt, buffer_km = input$lm_buf_km)
    
    lims <- compute_limits_from_bbox(rv$bb_final, terre, terre_crs)
    rv$bb_terre   <- lims$bb_terre
    rv$terre_crop <- lims$terre_crop
    rv$xlim_terre <- lims$xlim_terre
    rv$ylim_terre <- lims$ylim_terre
    rv$grat_final <- sf::st_graticule(bbox = rv$bb_final, crs = sf::st_crs(4326), ndx = 8, ndy = 8) |>
      sf::st_transform(terre_crs) |> sf::st_sf()
    
    leafletProxy("map_bbox") %>%
      clearShapes() %>%
      addMarkers(data = sf::st_transform(pt, 4326)) %>%
      addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                    lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                    color = "magenta", weight = 2) %>%
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"], rv$bb_final["xmax"], rv$bb_final["ymax"])
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
    
    leafletProxy("map_bbox") %>%
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
  
  # Points (Airports & Municipalities)
  observeEvent(input$refresh_points, {
    req(rv$bb_final)
    rv$air_all <- load_airports(AIRPORTS_CSV)
    buf_m <- if (input$buf_km > 0) input$buf_km*1000 else scaled_buffer_m(rv$bb_final)
    
    rv$air_sf <- filter_points_after_bbox(rv$air_all, rv$bb_final, terre_crs,
                                          aoi_terre = rv$aoi_terre,
                                          use_aoi_buf = isTRUE(input$limit_to_aoi_buf),
                                          buf_m = buf_m)
    rv$air_sf <- make_airport_labels(rv$air_sf)
    
    muni_raw <- geonames_bbox_paged(rv$bb_final, province = NULL, concise = "CITY")
    muni_ll  <- to_points_safe(muni_raw, rv$bb_final) |> bbox_filter_points(rv$bb_final)
    rv$muni_sf <- sf::st_transform(muni_ll, terre_crs)
  }, ignoreInit = TRUE)
  
  output$pts_counts <- renderText({
    if (is.null(rv$air_sf) || is.null(rv$muni_sf)) return("Click 'Apply Filter'")
    sprintf("Airports: %d | Municipalities: %d", nrow(rv$air_sf), nrow(rv$muni_sf))
  })
  
  output$map_points <- renderLeaflet({
    req(rv$bb_final, rv$terre_crop)
    p <- leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(providers$Esri.WorldLightGrayBase) |>
      addScaleBar(position = "bottomleft")
    if (!is.null(rv$aoi_ll)) p <- p |> addPolygons(data = sf::st_transform(rv$aoi_ll, 4326),
                                                   color = "black", weight = 2, fill = FALSE)
    if (!is.null(rv$air_sf) && nrow(rv$air_sf) > 0)
      p <- p |> addCircleMarkers(data = sf::st_transform(rv$air_sf, 4326),
                                 radius = 4, color = "blue", label = ~airport_label)
    if (!is.null(rv$muni_sf) && nrow(rv$muni_sf) > 0)
      p <- p |> addCircleMarkers(data = sf::st_transform(rv$muni_sf, 4326),
                                 radius = 3, color = "red", label = ~as.character(name))
    p |> fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"], rv$bb_final["xmax"], rv$bb_final["ymax"])
  })
  
  # Transects
  output$map_transects <- renderLeaflet({
    req(rv$bb_final)
    leaflet(options = leafletOptions(worldCopyJump = FALSE)) |>
      addProviderTiles(providers$Esri.WorldLightGrayBase) |>
      addScaleBar(position = "bottomleft") |>
      fitBounds(rv$bb_final["xmin"], rv$bb_final["ymin"], rv$bb_final["xmax"], rv$bb_final["ymax"])
  })
  
  observeEvent(input$map_transects_click, {
    rv$anchor_lonlat <- c(input$map_transects_click$lng, input$map_transects_click$lat)
  })
  
  # Load existing transects
  observeEvent(input$tr_file, {
    tr_raw <- read_vector_upload(input$tr_file)
    tr_raw <- sf::st_make_valid(tr_raw) |> sf::st_transform(terre_crs)
    tr_sf  <- clean_lines_for_optimizer(tr_raw)
    if (!("id" %in% names(tr_sf))) tr_sf$id <- seq_len(nrow(tr_sf))
    rv$transects_sf <- tr_sf
  })
  
  # Build transects from controls
  observeEvent(input$build_tr, {
    req(rv$bb_final)
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
    tr_sf <- sf::st_transform(tr_ll, terre_crs)
    tr_sf <- clean_lines_for_optimizer(tr_sf)
    if (!("id" %in% names(tr_sf))) tr_sf$id <- seq_len(nrow(tr_sf))
    rv$transects_sf <- tr_sf
  })
  
  observe({
    req(rv$bb_final, rv$transects_sf)
    p <- leafletProxy("map_transects") %>% clearShapes() %>% clearMarkers()
    if (!is.null(rv$aoi_ll)) {
      p <- p %>% addPolygons(data = sf::st_transform(rv$aoi_ll, 4326), color = "black", weight = 2, fill = FALSE)
    }
    p %>%
      addRectangles(lng1 = rv$bb_final["xmin"], lat1 = rv$bb_final["ymin"],
                    lng2 = rv$bb_final["xmax"], lat2 = rv$bb_final["ymax"],
                    color = "magenta", weight = 2) %>%
      addPolylines(data = sf::st_transform(rv$transects_sf, 4326), color = "red", weight = 2,
                   label = ~paste0("T", id))
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
    }
  })
  
  # Optimize
  observeEvent(input$run_opt, {
    req(rv$transects_sf)
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
    
    res <- run_transect_optimizer(
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
      make_plot    = TRUE,
      terre_crop   = rv$terre_crop,
      xlim_terre   = rv$xlim_terre,
      ylim_terre   = rv$ylim_terre,
      # overlays for optimizer quick PNG:
      air_sf_plot  = rv$air_sf,
      muni_sf_plot = rv$muni_sf,
      aoi_sf_plot  = rv$aoi_terre
    )
    rv$routes_sf    <- res$routes_sf
    rv$trip_summary <- res$trip_summary
    rv$manifest_df  <- res$flight_manifest
  })
  
  output$map_routes_png <- renderPlot({
    req(rv$routes_sf, rv$transects_sf)
    plot_final_map(
      terre_crop = rv$terre_crop, grat_final = rv$grat_final,
      transects_sf = rv$transects_sf, tran_labels_sf = NULL,
      muni_sf = rv$muni_sf, air_sf = rv$air_sf,
      xlim_terre = rv$xlim_terre, ylim_terre = rv$ylim_terre,
      bb_txt = "Optimized routes",
      out_png = NULL,
      terre_crs = terre_crs,
      tran_label_offset_frac = 0.012,
      routes_sf = rv$routes_sf,
      aoi_sf = rv$aoi_terre
    )
  })
  
  output$trips_table <- renderDT({
    req(rv$trip_summary); datatable(rv$trip_summary, options = list(pageLength = 5))
  })
  output$manifest_table <- renderDT({
    req(rv$manifest_df); datatable(rv$manifest_df, options = list(pageLength = 5))
  })
  
  # Downloads
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
  output$dl_routes_png <- downloadHandler(
    filename = function() sprintf("%s_map_with_routes.png", input$out_basename),
    content = function(file) {
      req(rv$routes_sf, rv$transects_sf)
      tmp_png <- tempfile(fileext = ".png")
      plot_final_map(
        terre_crop = rv$terre_crop, grat_final = rv$grat_final,
        transects_sf = rv$transects_sf, tran_labels_sf = NULL,
        muni_sf = rv$muni_sf, air_sf = rv$air_sf,
        xlim_terre = rv$xlim_terre, ylim_terre = rv$ylim_terre,
        bb_txt = "Optimized routes",
        out_png = tmp_png,
        terre_crs = terre_crs,
        tran_label_offset_frac = 0.012,
        routes_sf = rv$routes_sf,
        aoi_sf = rv$aoi_terre
      )
      file.copy(tmp_png, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)