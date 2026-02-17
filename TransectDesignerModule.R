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
  anchor_pt_ll <- st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon","lat"), crs = 4326)
  anchor_pt_utm <- st_transform(anchor_pt_ll, utm)
  
  # Direction unit vectors in UTM for given bearing
  # Bearing: 0 = North, 90 = East
  theta <- brg_deg * pi / 180
  dx <- sin(theta)
  dy <- cos(theta)
  
  # Perpendicular vector (to the left)
  pdx <- -dy
  pdy <- dx
  
  L <- len_km * 1000
  S <- sp_km * 1000
  
  # Center the stack around the anchor (so anchor is middle line)
  # offsets = ..., -2S, -S, 0, S, 2S, ...
  offsets <- ((seq_len(n) - (n + 1)/2) * S)
  
  # Create lines in UTM
  axy <- st_coordinates(anchor_pt_utm)[1, ]
  lines <- lapply(seq_len(n), function(i) {
    off <- offsets[i]
    x0 <- axy["X"] + off * pdx
    y0 <- axy["Y"] + off * pdy
    
    # line endpoints centered at anchor, length L along direction (dx,dy)
    x1 <- x0 - (L/2) * dx
    y1 <- y0 - (L/2) * dy
    x2 <- x0 + (L/2) * dx
    y2 <- y0 + (L/2) * dy
    
    st_linestring(matrix(c(x1,y1,x2,y2), ncol = 2, byrow = TRUE))
  })
  
  tran_utm <- st_sf(
    id = seq_len(n),
    geometry = st_sfc(lines, crs = utm)
  )
  
  # Optionally clip to bbox (bbox_ll in EPSG:4326)
  if (clip_to_bbox && !is.null(bbox_ll)) {
    bbox_poly_ll <- st_as_sfc(validate_bbox(bbox_ll)) |> st_set_crs(4326)
    bbox_poly_utm <- st_transform(bbox_poly_ll, utm)
    tran_utm <- st_intersection(tran_utm, bbox_poly_utm)
  }
  
  st_transform(tran_utm, 4326)
}

# Plot preview: full terre context, bbox outline, and transects
plot_transect_preview <- function(bb_ll, terre, pt_landmark_ll, tran_ll,
                                  air_sf = NULL, muni_sf = NULL,
                                  show_full_terre = TRUE) {
  
  bb_ll <- validate_bbox(bb_ll)
  
  # Force numeric bbox limits (prevents ggplot from ignoring them)
  xlim <- c(as.numeric(bb_ll[["xmin"]]), as.numeric(bb_ll[["xmax"]]))
  ylim <- c(as.numeric(bb_ll[["ymin"]]), as.numeric(bb_ll[["ymax"]]))
  
  # bbox polygon overlay
  bb_poly <- st_as_sfc(bb_ll) |> st_set_crs(4326)
  
  # transform layers to lon/lat for consistent preview
  terre_ll <- st_transform(terre, 4326)
  tran_ll  <- st_transform(tran_ll, 4326)
  
  # draw only cropped terre for speed OR full terre if requested
  terre_draw <- if (show_full_terre) terre_ll else st_crop(terre_ll, st_bbox(bb_poly))
  
  p <- ggplot() +
    geom_sf(data = terre_draw, fill = "grey97", color = "grey85", linewidth = 0.15) +
    geom_sf(data = st_sf(geometry = bb_poly), fill = NA, color = "magenta", linewidth = 1) +
    geom_sf(data = tran_ll, color = "red", linewidth = 1) +
    geom_sf_text(data = tran_ll, aes(label = id), color = "red4", size = 3,
                 check_overlap = TRUE, show.legend = FALSE) +
    { if (!is.null(pt_landmark_ll)) geom_sf(data = st_transform(pt_landmark_ll, 4326),
                                            color = "blue", size = 3) } +
    # Optional overlays (points)
    { if (!is.null(air_sf) && nrow(air_sf) > 0)
      geom_sf(data = st_transform(air_sf, 4326), color = "blue", size = 1.1) } +
    { if (!is.null(muni_sf) && nrow(muni_sf) > 0)
      geom_sf(data = st_transform(muni_sf, 4326), color = "red", size = 1.1) } +
    # IMPORTANT: force bbox view and CRS
    coord_sf(crs = st_crs(4326), xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_minimal() +
    labs(title = sprintf("Transect preview — BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f",
                         xlim[1], ylim[1], xlim[2], ylim[2]))
  
  print(p)
}

# Interactive editor for transects

refine_transects <- function(bb_ll, terre, pt_landmark_ll,
                             air_sf = NULL, muni_sf = NULL,
                             n = 10, len_km = 20, sp_km = 2, brg_deg = 90,
                             anchor_lonlat = NULL,
                             clip_to_bbox = TRUE,
                             show_full_terre = TRUE) {
  
  
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
    cat("  move <dir> [km]            -> move anchor (dir = n/s/e/w or north/south/east/west)\n")
    cat("                               e.g., move east 20\n")
    cat("  clip on|off                -> clip transects to bbox\n\n")
    
    plot_transect_preview(bb_ll, terre, pt_landmark_ll, tran_ll,
                          air_sf = air_sf, muni_sf = muni_sf,
                          show_full_terre = show_full_terre)
    
    ans <- trimws(get_input_safe("Transect choice: ", "ok"))
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
        utm <- utm_epsg_from_lonlat(anchor_lonlat[1], anchor_lonlat[2])
        
        pt0  <- st_as_sf(data.frame(lon=anchor_lonlat[1], lat=anchor_lonlat[2]),
                         coords=c("lon","lat"), crs=4326)
        pt0m <- st_transform(pt0, utm)
        xy   <- st_coordinates(pt0m)[1,]
        
        dx <- 0; dy <- 0
        if (dir == "n") dy <-  step_km * 1000
        if (dir == "s") dy <- -step_km * 1000
        if (dir == "e") dx <-  step_km * 1000
        if (dir == "w") dx <- -step_km * 1000
        
        pt1m <- st_as_sf(data.frame(x = xy["X"] + dx, y = xy["Y"] + dy),
                         coords=c("x","y"), crs=utm)
        pt1  <- st_transform(pt1m, 4326)
        ll   <- st_coordinates(pt1)[1,]
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
      if (tok[2] %in% c("on","true","1")) clip_to_bbox <- TRUE
      if (tok[2] %in% c("off","false","0")) clip_to_bbox <- FALSE
      next
    }
    
    # panning anchor in km; default step = spacing
    dir <- tok[1]
    if (dir %in% c("n","s","e","w")) {
      step_km <- if (length(tok) >= 2) suppressWarnings(as.numeric(tok[2])) else sp_km
      if (!is.finite(step_km) || step_km <= 0) step_km <- sp_km
      
      # move anchor in UTM meters then convert back to lon/lat
      utm <- utm_epsg_from_lonlat(anchor_lonlat[1], anchor_lonlat[2])
      pt0 <- st_as_sf(data.frame(lon=anchor_lonlat[1], lat=anchor_lonlat[2]),
                      coords=c("lon","lat"), crs=4326)
      pt0m <- st_transform(pt0, utm)
      xy <- st_coordinates(pt0m)[1,]
      dx <- 0; dy <- 0
      if (dir == "n") dy <- step_km*1000
      if (dir == "s") dy <- -step_km*1000
      if (dir == "e") dx <- step_km*1000
      if (dir == "w") dx <- -step_km*1000
      
      pt1m <- st_as_sf(data.frame(x = xy["X"] + dx, y = xy["Y"] + dy), coords=c("x","y"), crs=utm)
      pt1 <- st_transform(pt1m, 4326)
      ll <- st_coordinates(pt1)[1,]
      anchor_lonlat <- c(ll["X"], ll["Y"])
      next
    }
    
    message("Unknown command.")
  }
}