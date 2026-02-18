# ============================================================
# TRANSECT DESIGNER (interactive)
# Creates N parallel lines of length (km), spacing (km), at bearing (deg)
# Anchor point is moved with n/s/e/w (km), like panning.
# ============================================================

# Compute UTM EPSG from lon/lat (north hemisphere)
utm_epsg_from_lonlat <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  # assume northern hemisphere for your area
  32600 + zone
}

# Build transects in UTM (meters), then transform back to lon/lat (4326)
make_transects_ll <- function(anchor_lonlat, n = 10, len_km = 20, sp_km = 2,
                              brg_deg = 90, clip_to_bbox = TRUE, bbox_ll = NULL) {
  stopifnot(length(anchor_lonlat) == 2)
  lon <- anchor_lonlat[1]; lat <- anchor_lonlat[2]
  utm <- utm_epsg_from_lonlat(lon, lat)
  
  # Anchor point as sf
  anchor_pt_ll  <- sf::st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon","lat"), crs = 4326)
  anchor_pt_utm <- sf::st_transform(anchor_pt_ll, utm)
  
  # Direction unit vectors in UTM for given bearing (0 = North, 90 = East)
  theta <- brg_deg * pi / 180
  dx <- sin(theta); dy <- cos(theta)
  
  # Perpendicular vector (to the left)
  pdx <- -dy; pdy <- dx
  
  L <- len_km * 1000
  S <- sp_km  * 1000
  
  # Center the stack around the anchor (so anchor is middle line)
  # offsets = ..., -2S, -S, 0, S, 2S, ...
  offsets <- ((seq_len(n) - (n + 1)/2) * S)
  
  # Create lines in UTM
  axy   <- sf::st_coordinates(anchor_pt_utm)[1, ]
  lines <- lapply(seq_len(n), function(i) {
    off <- offsets[i]
    x0 <- axy["X"] + off * pdx
    y0 <- axy["Y"] + off * pdy
    
    # line endpoints centered at anchor, length L along direction (dx,dy)
    x1 <- x0 - (L/2) * dx
    y1 <- y0 - (L/2) * dy
    x2 <- x0 + (L/2) * dx
    y2 <- y0 + (L/2) * dy
    
    sf::st_linestring(matrix(c(x1,y1,x2,y2), ncol = 2, byrow = TRUE))
  })
  
  tran_utm <- sf::st_sf(
    id = seq_len(n),
    geometry = sf::st_sfc(lines, crs = utm)
  )
  
  # Optionally clip to bbox (bbox_ll in EPSG:4326)
  if (clip_to_bbox && !is.null(bbox_ll)) {
    bbox_poly_ll  <- sf::st_as_sfc(validate_bbox(bbox_ll)) |> sf::st_set_crs(4326)
    bbox_poly_utm <- sf::st_transform(bbox_poly_ll, utm)
    tran_utm      <- sf::st_intersection(tran_utm, bbox_poly_utm)
  }
  
  sf::st_transform(tran_utm, 4326)
}

# Plot preview: full terre context, bbox outline, and transects
# Plot preview: terre base, bbox outline, transects, airports, municipalities, and optional AOI
plot_transect_preview <- function(bb_ll, terre, pt_landmark_ll, tran_ll,
                                  air_sf = NULL, muni_sf = NULL,
                                  aoi_ll = NULL,                # <--- NEW: AOI in EPSG:4326 (optional)
                                  show_full_terre = TRUE) {
  
  bb_ll <- validate_bbox(bb_ll)
  
  # Numeric bbox limits (prevents ggplot from auto-expanding)
  xlim <- c(as.numeric(bb_ll[["xmin"]]), as.numeric(bb_ll[["xmax"]]))
  ylim <- c(as.numeric(bb_ll[["ymin"]]), as.numeric(bb_ll[["ymax"]]))
  
  # bbox polygon overlay (4326)
  bb_poly <- sf::st_as_sfc(bb_ll) |> sf::st_set_crs(4326)
  
  # transform layers to lon/lat for consistent preview
  terre_ll <- sf::st_transform(terre, 4326)
  tran_ll  <- sf::st_transform(tran_ll, 4326)
  
  # draw only cropped terre for speed OR full terre if requested
  terre_draw <- if (show_full_terre) terre_ll else sf::st_crop(terre_ll, sf::st_bbox(bb_poly))
  
  p <- ggplot2::ggplot() +
    # base land
    ggplot2::geom_sf(data = terre_draw, fill = "grey97", color = "grey85", linewidth = 0.15) +
    # AOI outline (if provided)
    { if (!is.null(aoi_ll)) ggplot2::geom_sf(data = aoi_ll, fill = NA, color = "black", linewidth = 0.7) } +
    # bbox outline
    ggplot2::geom_sf(data = sf::st_sf(geometry = bb_poly), fill = NA, color = "magenta", linewidth = 1) +
    # transects + labels
    ggplot2::geom_sf(data = tran_ll, color = "red", linewidth = 1) +
    ggplot2::geom_sf_text(data = tran_ll, ggplot2::aes(label = id), color = "red4", size = 3,
                          check_overlap = TRUE, show.legend = FALSE) +
    # landmark point (optional)
    { if (!is.null(pt_landmark_ll)) ggplot2::geom_sf(data = sf::st_transform(pt_landmark_ll, 4326),
                                                     color = "blue", size = 3) } +
    # airports (optional)
    { if (!is.null(air_sf) && nrow(air_sf) > 0)
      ggplot2::geom_sf(data = sf::st_transform(air_sf, 4326), color = "blue", size = 1.1) } +
    # municipalities (optional)
    { if (!is.null(muni_sf) && nrow(muni_sf) > 0)
      ggplot2::geom_sf(data = sf::st_transform(muni_sf, 4326), color = "red", size = 1.1) } +
    # lock the view to bbox
    ggplot2::coord_sf(crs = sf::st_crs(4326), xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("Transect preview — BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f",
                                  xlim[1], ylim[1], xlim[2], ylim[2]))
  
  print(p)
}

# ------------------------------------------------------------
# Interactive BBox refiner starting from an existing bbox (e.g., AOI extent)
# Commands:
#   ok                   -> accept bbox
#   zin / zout           -> zoom in/out by zoom_factor (default 1.5x)
#   pan n|s|e|w [km]     -> pan by km (default pan_step_km)
#   set xmin|xmax|ymin|ymax <val>  -> set side explicitly (degrees)
#   show full|crop       -> toggle full/cropped 'terre' in preview
#   help                 -> show commands
#
# Inputs/Outputs: bbox in EPSG:4326 (named: xmin ymin xmax ymax)
# ------------------------------------------------------------
refine_bbox_from_bbox <- function(bb_ll,
                                  terre,
                                  aoi_ll = NULL,             # optional AOI outline in 4326
                                  initial_zoom_factor = 1.0,
                                  pan_step_km = 10,
                                  zoom_factor = 1.5,
                                  show_full_terre = FALSE) {
  
  # --- normalize bbox ---
  stopifnot(all(c("xmin","ymin","xmax","ymax") %in% names(bb_ll)))
  bb_ll <- sf::st_bbox(c(
    xmin = as.numeric(bb_ll["xmin"]),
    ymin = as.numeric(bb_ll["ymin"]),
    xmax = as.numeric(bb_ll["xmax"]),
    ymax = as.numeric(bb_ll["ymax"])
  ), crs = sf::st_crs(4326))
  
  # initial zoom if requested
  if (is.finite(initial_zoom_factor) && initial_zoom_factor != 1) {
    cx <- (bb_ll["xmin"] + bb_ll["xmax"]) / 2
    cy <- (bb_ll["ymin"] + bb_ll["ymax"]) / 2
    w  <- (bb_ll["xmax"] - bb_ll["xmin"]) / initial_zoom_factor
    h  <- (bb_ll["ymax"] - bb_ll["ymin"]) / initial_zoom_factor
    bb_ll["xmin"] <- cx - w/2; bb_ll["xmax"] <- cx + w/2
    bb_ll["ymin"] <- cy - h/2; bb_ll["ymax"] <- cy + h/2
  }
  
  # internal preview
  .plot_bbox_preview <- function(bb_ll, terre, show_full_terre, aoi_ll = NULL) {
    bb_poly  <- sf::st_as_sfc(bb_ll) |> sf::st_set_crs(4326)
    terre_ll <- sf::st_transform(terre, 4326)
    terre_draw <- if (isTRUE(show_full_terre)) terre_ll else sf::st_crop(terre_ll, sf::st_bbox(bb_poly))
    xlim <- c(as.numeric(bb_ll["xmin"]), as.numeric(bb_ll["xmax"]))
    ylim <- c(as.numeric(bb_ll["ymin"]), as.numeric(bb_ll["ymax"]))
    p <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = terre_draw, fill = "grey97", color = "grey85", linewidth = 0.2) +
      { if (!is.null(aoi_ll)) ggplot2::geom_sf(data = aoi_ll, fill = NA, color = "black", linewidth = 0.6) } +
      ggplot2::geom_sf(data = sf::st_sf(geometry = bb_poly), fill = NA, color = "magenta", linewidth = 1) +
      ggplot2::coord_sf(crs = sf::st_crs(4326), xlim = xlim, ylim = ylim, expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = sprintf("BBox refine — lon: %.4f..%.4f  lat: %.4f..%.4f",
                        xlim[1], xlim[2], ylim[1], ylim[2])
      )
    print(p)
  }
  
  .help <- function() {
    cat("\nCommands:\n",
        "  ok                    -> accept current bbox\n",
        "  zin                   -> zoom in (", zoom_factor, "x)\n",
        "  zout                  -> zoom out (", zoom_factor, "x)\n",
        "  pan <dir> [km]        -> pan n/s/e/w by km (default ", pan_step_km, ")\n",
        "  set xmin <lon>        -> set xmin (degrees)\n",
        "  set xmax <lon>        -> set xmax (degrees)\n",
        "  set ymin <lat>        -> set ymin (degrees)\n",
        "  set ymax <lat>        -> set ymax (degrees)\n",
        "  show full|crop        -> draw full terre or cropped (faster)\n",
        "  help                  -> show commands\n\n", sep = "")
  }
  
  .help()
  repeat {
    .plot_bbox_preview(bb_ll, terre, show_full_terre, aoi_ll = aoi_ll)
    ans <- tolower(trimws(get_input_safe("BBox refine command: ", "ok")))
    if (ans %in% c("ok","accept")) break
    if (ans %in% c("help","?")) { .help(); next }
    
    tok <- strsplit(ans, "\\s+")[[1]]
    if (length(tok) == 0) next
    cmd <- tok[1]
    
    if (cmd == "zin") {
      cx <- (bb_ll["xmin"] + bb_ll["xmax"]) / 2
      cy <- (bb_ll["ymin"] + bb_ll["ymax"]) / 2
      w  <- (bb_ll["xmax"] - bb_ll["xmin"]) / zoom_factor
      h  <- (bb_ll["ymax"] - bb_ll["ymin"]) / zoom_factor
      bb_ll["xmin"] <- cx - w/2; bb_ll["xmax"] <- cx + w/2
      bb_ll["ymin"] <- cy - h/2; bb_ll["ymax"] <- cy + h/2
      next
    }
    
    if (cmd == "zout") {
      cx <- (bb_ll["xmin"] + bb_ll["xmax"]) / 2
      cy <- (bb_ll["ymin"] + bb_ll["ymax"]) / 2
      w  <- (bb_ll["xmax"] - bb_ll["xmin"]) * zoom_factor
      h  <- (bb_ll["ymax"] - bb_ll["ymin"]) * zoom_factor
      bb_ll["xmin"] <- cx - w/2; bb_ll["xmax"] <- cx + w/2
      bb_ll["ymin"] <- cy - h/2; bb_ll["ymax"] <- cy + h/2
      next
    }
    
    if (cmd == "show" && length(tok) >= 2) {
      if (tok[2] == "full") show_full_terre <- TRUE
      if (tok[2] == "crop") show_full_terre <- FALSE
      next
    }
    
    if (cmd == "pan" && length(tok) >= 2) {
      dir <- tok[2]
      step_km <- if (length(tok) >= 3) suppressWarnings(as.numeric(tok[3])) else pan_step_km
      if (!is.finite(step_km) || step_km <= 0) step_km <- pan_step_km
      
      # derive a local UTM from current bbox center
      b <- bb_ll
      mid_lon <- as.numeric((b["xmin"] + b["xmax"]) / 2)
      mid_lat <- as.numeric((b["ymin"] + b["ymax"]) / 2)
      zone <- floor((mid_lon + 180) / 6) + 1
      crs_m <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
      
      bb_poly <- sf::st_as_sfc(bb_ll) |> sf::st_set_crs(4326)
      bb_m    <- sf::st_transform(bb_poly, crs_m)
      bb_env  <- sf::st_bbox(bb_m)
      
      dx <- 0; dy <- 0
      if (dir %in% c("n","north")) dy <-  step_km*1000
      if (dir %in% c("s","south")) dy <- -step_km*1000
      if (dir %in% c("e","east"))  dx <-  step_km*1000
      if (dir %in% c("w","west"))  dx <- -step_km*1000
      
      # move bbox corners in meters then back to 4326
      p1 <- sf::st_sfc(sf::st_point(c(bb_env["xmin"] + dx, bb_env["ymin"] + dy)), crs = crs_m)
      p2 <- sf::st_sfc(sf::st_point(c(bb_env["xmax"] + dx, bb_env["ymax"] + dy)), crs = crs_m)
      bb2 <- sf::st_bbox(sf::st_union(p1, p2) |> sf::st_as_sfc())
      bb2_ll <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb2), 4326))
      bb_ll <- sf::st_bbox(c(xmin = bb2_ll["xmin"], ymin = bb2_ll["ymin"],
                             xmax = bb2_ll["xmax"], ymax = bb2_ll["ymax"]),
                           crs = sf::st_crs(4326))
      next
    }
    
    if (cmd == "set" && length(tok) >= 3) {
      fld <- tok[2]; v <- suppressWarnings(as.numeric(tok[3]))
      if (fld %in% c("xmin","xmax","ymin","ymax") && is.finite(v)) {
        bb_ll[fld] <- v
      } else {
        message("Usage: set xmin <lon> | set xmax <lon> | set ymin <lat> | set ymax <lat>")
      }
      next
    }
    
    message("Unknown command. Type 'help' for options.")
  }
  
  bb_ll
}

# Interactive editor for transects

refine_transects <- function(bb_ll, terre, pt_landmark_ll,
                             air_sf = NULL, muni_sf = NULL,
                             n = 10, len_km = 20, sp_km = 2, brg_deg = 90,
                             anchor_lonlat = NULL,
                             clip_to_bbox = TRUE,
                             show_full_terre = TRUE,
                             aoi_ll = NULL)
                             {
  
  bb_ll <- validate_bbox(bb_ll)
  
  # Default anchor: bbox center
  if (is.null(anchor_lonlat)) {
    anchor_lonlat <- c(
      as.numeric((bb_ll["xmin"] + bb_ll["xmax"]) / 2),
      as.numeric((bb_ll["ymin"] + bb_ll["ymax"]) / 2)
    )
  }
  
  repeat {
    tran_ll <- make_transects_ll(anchor_lonlat, n, len_km, sp_km, brg_deg,
                                 clip_to_bbox = clip_to_bbox, bbox_ll = bb_ll)
    
    cat("\nTransects settings:\n")
    cat(sprintf("  n=%d  len_km=%.2f  sp_km=%.2f  brg_deg=%.1f\n", n, len_km, sp_km, brg_deg))
    cat(sprintf("  anchor lon/lat = %.5f, %.5f\n", anchor_lonlat[1], anchor_lonlat[2]))
    cat("Commands:\n")
    cat("  ok                         -> accept transects\n")
    cat("  n <int>                    -> number of lines (e.g., n 17)\n")
    cat("  len <km>                   -> length in km (e.g., len 25)\n")
    cat("  sp <km>                    -> spacing in km (e.g., sp 1.5)\n")
    cat("  brg <deg>                  -> bearing degrees (0=N,90=E) (e.g., brg 45)\n")
    cat("  anchor <lon> <lat>         -> set anchor (e.g., anchor -64.2 48.8)\n")
    cat("  move <dir> [km]            -> move anchor (n/s/e/w or north/south/east/west)\n")
    cat("                                 e.g., move east 20\n")
    cat("  clip on|off                -> clip transects to bbox\n\n")
    
    
    plot_transect_preview(bb_ll, terre, pt_landmark_ll, tran_ll,
                          air_sf = air_sf, muni_sf = muni_sf,
                          aoi_ll = aoi_ll,                  # <--- pass AOI (shows outline if present)
                          show_full_terre = show_full_terre)
    
    
    ans  <- trimws(get_input_safe("Transect choice: ", "ok"))
    ans2 <- tolower(ans)
    if (ans2 == "ok") return(tran_ll)
    
    # split command tokens
    tok <- strsplit(ans, "\\s+")[[1]]
    
    # ---- COMMAND NORMALIZATION ----
    # Allow: move <dir> [km], where dir can be n/s/e/w or north/south/east/west
    if (length(tok) >= 2 && tok[1] == "move") {
      dir <- tolower(tok[2])
      
      # normalize long directions
      if (dir == "north") dir <- "n"
      if (dir == "south") dir <- "s"
      if (dir == "east")  dir <- "e"
      if (dir == "west")  dir <- "w"
      
      step_km <- if (length(tok) >= 3) suppressWarnings(as.numeric(tok[3])) else sp_km
      if (!is.finite(step_km) || step_km <= 0) step_km <- sp_km
      
      if (dir %in% c("n","s","e","w")) {
        utm  <- utm_epsg_from_lonlat(anchor_lonlat[1], anchor_lonlat[2])
        pt0  <- sf::st_as_sf(data.frame(lon=anchor_lonlat[1], lat=anchor_lonlat[2]),
                             coords=c("lon","lat"), crs=4326)
        pt0m <- sf::st_transform(pt0, utm)
        xy   <- sf::st_coordinates(pt0m)[1,]
        
        dx <- 0; dy <- 0
        if (dir == "n") dy <-  step_km * 1000
        if (dir == "s") dy <- -step_km * 1000
        if (dir == "e") dx <-  step_km * 1000
        if (dir == "w") dx <- -step_km * 1000
        
        pt1m <- sf::st_as_sf(data.frame(x = xy["X"] + dx, y = xy["Y"] + dy),
                             coords=c("x","y"), crs=utm)
        pt1  <- sf::st_transform(pt1m, 4326)
        ll   <- sf::st_coordinates(pt1)[1,]
        anchor_lonlat <- c(ll["X"], ll["Y"])
        
        cat(sprintf("Moved anchor %s by %.2f km -> lon/lat = %.5f, %.5f\n",
                    dir, step_km, anchor_lonlat[1], anchor_lonlat[2]))
        next
      } else {
        message("Move direction must be one of: n/s/e/w/north/south/east/west")
        next
      }
    }
    
    # parameter setters
    if (length(tok) >= 2 && tok[1] == "n") {
      n_new <- suppressWarnings(as.integer(tok[2]))
      if (is.finite(n_new) && n_new > 0) n <- n_new else message("Bad n")
      next
    }
    if (length(tok) >= 2 && tok[1] == "len") {
      v <- suppressWarnings(as.numeric(tok[2]))
      if (is.finite(v) && v > 0) len_km <- v else message("Bad len")
      next
    }
    if (length(tok) >= 2 && tok[1] == "sp") {
      v <- suppressWarnings(as.numeric(tok[2]))
      if (is.finite(v) && v > 0) sp_km <- v else message("Bad sp")
      next
    }
    if (length(tok) >= 2 && tok[1] == "brg") {
      v <- suppressWarnings(as.numeric(tok[2]))
      if (is.finite(v)) brg_deg <- v %% 360 else message("Bad brg")
      next
    }
    if (length(tok) >= 3 && tok[1] == "anchor") {
      lo <- suppressWarnings(as.numeric(tok[2]))
      la <- suppressWarnings(as.numeric(tok[3]))
      if (is.finite(lo) && is.finite(la)) anchor_lonlat <- c(lo, la) else message("Bad anchor")
      next
    }
    
    # clip toggle
    if (length(tok) >= 2 && tok[1] == "clip") {
      if (tok[2] %in% c("on","true","1"))  clip_to_bbox <- TRUE
      if (tok[2] %in% c("off","false","0")) clip_to_bbox <- FALSE
      next
    }
    
    # panning anchor in km; default step = spacing
    dir <- tok[1]
    if (dir %in% c("n","s","e","w")) {
      step_km <- if (length(tok) >= 2) suppressWarnings(as.numeric(tok[2])) else sp_km
      if (!is.finite(step_km) || step_km <= 0) step_km <- sp_km
      
      utm  <- utm_epsg_from_lonlat(anchor_lonlat[1], anchor_lonlat[2])
      pt0  <- sf::st_as_sf(data.frame(lon=anchor_lonlat[1], lat=anchor_lonlat[2]),
                           coords=c("lon","lat"), crs=4326)
      pt0m <- sf::st_transform(pt0, utm)
      xy   <- sf::st_coordinates(pt0m)[1,]
      dx <- 0; dy <- 0
      if (dir == "n") dy <-  step_km*1000
      if (dir == "s") dy <- -step_km*1000
      if (dir == "e") dx <-  step_km*1000
      if (dir == "w") dx <- -step_km*1000
      
      pt1m <- sf::st_as_sf(data.frame(x = xy["X"] + dx, y = xy["Y"] + dy), coords=c("x","y"), crs=utm)
      pt1  <- sf::st_transform(pt1m, 4326)
      ll   <- sf::st_coordinates(pt1)[1,]
      anchor_lonlat <- c(ll["X"], ll["Y"])
      next
    }
    
    message("Unknown command.")
  }
}
