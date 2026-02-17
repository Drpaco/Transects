get_input_safe <- function(prompt, default = "") {
  if (!interactive()) return(default)
  ans <- tryCatch(readline(prompt), error = function(e) default)
  ans <- trimws(ans)
  if (nchar(ans) == 0) default else ans
}

validate_bbox <- function(bb) {
  v <- as.numeric(bb[c("xmin","ymin","xmax","ymax")])
  if (any(!is.finite(v))) stop("BBox has NA/Inf values. Enter bbox manually or restart.")
  if (v[1] >= v[3] || v[2] >= v[4]) stop("BBox invalid (xmin>=xmax or ymin>=ymax).")
  bb
}

bbox_span_km <- function(bb_ll) {
  b <- validate_bbox(bb_ll)
  midlat <- as.numeric((b["ymin"] + b["ymax"]) / 2)
  dx_km <- as.numeric(b["xmax"] - b["xmin"]) * 111 * cos(midlat * pi / 180)
  dy_km <- as.numeric(b["ymax"] - b["ymin"]) * 111
  c(width_km = dx_km, height_km = dy_km)
}

bbox_filter_points <- function(sf_points, bb_ll) {
  bb_ll <- validate_bbox(bb_ll)
  pts_ll <- sf::st_transform(sf_points, 4326)
  xy <- sf::st_coordinates(pts_ll)
  keep <- xy[,1] >= bb_ll["xmin"] & xy[,1] <= bb_ll["xmax"] &
    xy[,2] >= bb_ll["ymin"] & xy[,2] <= bb_ll["ymax"]
  sf_points[keep, , drop = FALSE]
}

to_points_safe <- function(x, bb_ll) {
  bb_ll <- validate_bbox(bb_ll)
  
  if (!inherits(x, "sf")) {
    stop("to_points_safe() expected an sf object but got: ", paste(class(x), collapse = ", "))
  }
  
  if (nrow(x) == 0) {
    return(sf::st_sf(geometry = sf::st_sfc(crs = 4326)))
  }
  
  lon0 <- as.numeric((bb_ll["xmin"] + bb_ll["xmax"]) / 2)
  zone <- floor((lon0 + 180) / 6) + 1
  utm_epsg <- 32600 + zone
  
  x_utm <- sf::st_transform(x, utm_epsg)
  gt <- sf::st_geometry_type(x_utm)
  
  nonpt <- which(!(gt %in% c("POINT","MULTIPOINT")))
  if (length(nonpt) > 0) {
    sf::st_geometry(x_utm)[nonpt] <- sf::st_point_on_surface(sf::st_geometry(x_utm)[nonpt])
  }
  
  x_pt <- sf::st_cast(x_utm, "POINT", warn = FALSE)
  sf::st_transform(x_pt, 4326)
}