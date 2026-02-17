# 06_transects.R
# Transect utilities: label points + shapefile export

utm_epsg_from_bbox <- function(bb_ll) {
  bb_ll <- validate_bbox(bb_ll)
  lon0 <- as.numeric((bb_ll["xmin"] + bb_ll["xmax"]) / 2)
  zone <- floor((lon0 + 180) / 6) + 1
  32600 + zone
}

# Label points placed at EACH transect's EAST end (rightmost end).
# IMPORTANT: no offset here — offset is applied in plotting as a fraction of map width.
transect_labels_at_east_end <- function(transects_sf, bb_ll, utm_epsg = NULL) {
  if (!inherits(transects_sf, "sf") || nrow(transects_sf) == 0) {
    return(sf::st_sf(id = integer(0), geometry = sf::st_sfc(crs = sf::st_crs(transects_sf))))
  }
  
  bb_ll <- validate_bbox(bb_ll)
  if (is.null(utm_epsg)) utm_epsg <- utm_epsg_from_bbox(bb_ll)
  
  # Ensure id exists
  if (!("id" %in% names(transects_sf))) transects_sf$id <- seq_len(nrow(transects_sf))
  
  # Work in meters for stable "eastmost end" detection
  tr_utm <- sf::st_transform(transects_sf, utm_epsg)
  
  pts <- vector("list", nrow(tr_utm))
  for (i in seq_len(nrow(tr_utm))) {
    g <- sf::st_geometry(tr_utm[i, ])
    
    # Works even if MULTILINESTRING after clipping
    parts <- sf::st_cast(g, "LINESTRING", warn = FALSE)
    coords <- do.call(rbind, lapply(parts, sf::st_coordinates))
    
    idx <- which.max(coords[, "X"])  # eastmost vertex
    x_end <- coords[idx, "X"]
    y_end <- coords[idx, "Y"]
    
    pts[[i]] <- sf::st_point(c(x_end, y_end))
  }
  
  lbl_utm <- sf::st_sf(id = tr_utm$id, geometry = sf::st_sfc(pts, crs = utm_epsg))
  sf::st_transform(lbl_utm, sf::st_crs(transects_sf))
}

# Export transects as ESRI Shapefile with user-chosen name
export_transects_shp <- function(transects_sf, export_dir, filename_base) {
  if (!inherits(transects_sf, "sf") || nrow(transects_sf) == 0) {
    message("No transects to export.")
    return(invisible(NULL))
  }
  
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
  
  filename_base <- gsub("[^A-Za-z0-9_\\-]", "_", filename_base)
  shp_path <- file.path(export_dir, paste0(filename_base, ".shp"))
  
  # Simple fields for shapefile
  if (!("id" %in% names(transects_sf))) transects_sf$id <- seq_len(nrow(transects_sf))
  transects_sf$label <- paste0("T", transects_sf$id)
  
  # Remove existing components if overwriting
  base <- tools::file_path_sans_ext(shp_path)
  unlink(paste0(base, c(".shp",".shx",".dbf",".prj",".cpg",".sbn",".sbx")), force = TRUE)
  
  sf::st_write(transects_sf[, c("id","label","geometry")],
               shp_path, driver = "ESRI Shapefile",
               delete_layer = TRUE, quiet = TRUE)
  
  message("Transects exported to: ", shp_path)
  invisible(shp_path)
}