bbox_from_point_km <- function(pt_sf, buffer_km = 25) {
  pt_ll <- sf::st_transform(pt_sf, 4326)
  lon <- sf::st_coordinates(pt_ll)[1,1]
  zone <- floor((lon + 180) / 6) + 1
  utm <- 32600 + zone
  
  buf <- sf::st_buffer(sf::st_transform(pt_ll, utm), buffer_km * 1000)
  sf::st_bbox(sf::st_transform(sf::st_as_sfc(sf::st_bbox(buf)), 4326))
}

plot_bbox_preview <- function(bb_ll, terre, pt_sf = NULL, grid_step_deg = 0.01) {
  bb_ll <- validate_bbox(bb_ll)
  terre_ll <- sf::st_transform(terre, 4326)
  bb_poly <- sf::st_as_sfc(bb_ll) |> sf::st_set_crs(4326)
  
  grat <- sf::st_graticule(
    bbox = bb_ll, crs = sf::st_crs(4326),
    ndx = max(2, round((bb_ll["xmax"] - bb_ll["xmin"]) / grid_step_deg)),
    ndy = max(2, round((bb_ll["ymax"] - bb_ll["ymin"]) / grid_step_deg))
  )
  
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = terre_ll, fill = "grey97", color = "grey75", linewidth = 0.2) +
    ggplot2::geom_sf(data = sf::st_sf(geometry = grat), color = "grey88", linewidth = 0.2) +
    ggplot2::geom_sf(data = sf::st_sf(geometry = bb_poly), fill = NA, color = "magenta", linewidth = 1) +
    { if (!is.null(pt_sf)) ggplot2::geom_sf(data = sf::st_transform(pt_sf, 4326), color = "blue", size = 3) } +
    ggplot2::coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
                      ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
                      expand = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = sprintf("BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f",
                                  bb_ll["xmin"], bb_ll["ymin"], bb_ll["xmax"], bb_ll["ymax"]))
}

refine_bbox <- function(pt_sf, terre, initial_buffer_km = 25, grid_step_deg = 0.01) {
  bb <- validate_bbox(bbox_from_point_km(pt_sf, initial_buffer_km))
  
  zoom_out_step <- 1.25
  zoom_in_step  <- 0.80
  pan_step      <- 0.20
  
  repeat {
    span <- bbox_span_km(bb)
    cat(sprintf("\nCurrent window ~ %.1f km (W–E) x %.1f km (S–N)\n", span["width_km"], span["height_km"]))
    print(plot_bbox_preview(bb, terre, pt_sf, grid_step_deg))
    
    cat("\nCommands:\n",
        "  ok                        -> accept\n",
        "  xmin,ymin,xmax,ymax        -> set bbox manually\n",
        "  zoom (any number)          -> <1 in, >1 out (e.g., 0.9, 1.2)\n",
        "  + / out                    -> step zoom out (", zoom_out_step, ")\n",
        "  - / in                     -> step zoom in  (", zoom_in_step, ")\n",
        "  pan: n/s/e/w [fraction]    -> move without zoom (default ", pan_step, ")\n", sep = "")
    
    ans <- trimws(get_input_safe("Your choice: ", "ok"))
    ans2 <- tolower(trimws(ans))
    if (ans2 == "ok") break
    
    if (grepl(",", ans2)) {
      vals <- as.numeric(strsplit(ans2, ",")[[1]])
      if (length(vals) == 4 && all(is.finite(vals))) {
        bb <- validate_bbox(sf::st_bbox(c(xmin = vals[1], ymin = vals[2], xmax = vals[3], ymax = vals[4]),
                                        crs = sf::st_crs(4326)))
      }
      next
    }
    
    pan_match <- regexec("^([nsew])\\s*([0-9]*\\.?[0-9]+)?$", ans2, perl = TRUE)
    pm <- regmatches(ans2, pan_match)[[1]]
    if (length(pm) > 0) {
      dir <- pm[2]
      frac <- if (length(pm) >= 3 && nzchar(pm[3])) as.numeric(pm[3]) else pan_step
      if (!is.finite(frac) || frac <= 0) frac <- pan_step
      
      xmin <- as.numeric(bb["xmin"]); xmax <- as.numeric(bb["xmax"])
      ymin <- as.numeric(bb["ymin"]); ymax <- as.numeric(bb["ymax"])
      dx <- (xmax - xmin) * frac
      dy <- (ymax - ymin) * frac
      
      if (dir == "n") { ymin <- ymin + dy; ymax <- ymax + dy }
      if (dir == "s") { ymin <- ymin - dy; ymax <- ymax - dy }
      if (dir == "e") { xmin <- xmin + dx; xmax <- xmax + dx }
      if (dir == "w") { xmin <- xmin - dx; xmax <- xmax - dx }
      
      xmin <- max(-180, xmin); xmax <- min(180, xmax)
      ymin <- max(-90, ymin);  ymax <- min(90, ymax)
      
      bb <- validate_bbox(sf::st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = sf::st_crs(4326)))
      next
    }
    
    if (ans2 %in% c("+", "out")) z <- zoom_out_step else
      if (ans2 %in% c("-", "in")) z <- zoom_in_step else
        z <- suppressWarnings(as.numeric(ans2))
    
    if (is.finite(z) && z > 0) {
      xmin <- as.numeric(bb["xmin"]); xmax <- as.numeric(bb["xmax"])
      ymin <- as.numeric(bb["ymin"]); ymax <- as.numeric(bb["ymax"])
      cx <- (xmin + xmax) / 2
      cy <- (ymin + ymax) / 2
      w <- (xmax - xmin) * z / 2
      h <- (ymax - ymin) * z / 2
      bb <- validate_bbox(sf::st_bbox(c(xmin = cx - w, ymin = cy - h, xmax = cx + w, ymax = cy + h),
                                      crs = sf::st_crs(4326)))
    } else {
      message("Unrecognized command.")
    }
  }
  
  bb
}