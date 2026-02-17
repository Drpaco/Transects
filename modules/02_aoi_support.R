# 02_aoi_support.R
# Helpers to load an AOI polygon, normalize CRS, and clip transects to AOI

read_aoi_polygon <- function(path) {
  if (!file.exists(path)) stop("AOI file not found: ", path)
  aoi <- sf::st_read(path, quiet = TRUE)
  if (!inherits(aoi, "sf")) stop("AOI is not an sf object: ", path)
  
  # keep only polygonal geometries; dissolve multipart
  aoi <- sf::st_make_valid(aoi)
  aoi <- sf::st_collection_extract(aoi, "POLYGON")
  if (nrow(aoi) == 0) stop("AOI file contains no POLYGON geometries.")
  aoi <- sf::st_union(sf::st_geometry(aoi)) |> sf::st_as_sf()
  aoi
}

# Clip any LINESTRING sf to AOI polygon; drop 0-length fragments
clip_transects_to_aoi <- function(transects_sf, aoi_sf) {
  stopifnot(inherits(transects_sf, "sf"), inherits(aoi_sf, "sf"))
  # Use a robust projected CRS for clipping
  # Derive local UTM from AOI bbox center
  aoi_ll <- sf::st_transform(aoi_sf, 4326)
  bb <- sf::st_bbox(aoi_ll)
  mid_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  mid_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  zone <- floor((mid_lon + 180) / 6) + 1
  crs_m <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
  
  tr_m  <- sf::st_transform(transects_sf, crs_m)
  aoi_m <- sf::st_transform(aoi_sf,      crs_m)
  
  clipped <- suppressWarnings(sf::st_intersection(tr_m, aoi_m))
  if (nrow(clipped) == 0) return(clipped)
  
  # Drop 0-length fragments
  len <- as.numeric(sf::st_length(clipped))
  clipped <- clipped[len > 0, , drop = FALSE]
  
  # Re‑project back to original CRS of transects
  sf::st_transform(clipped, sf::st_crs(transects_sf))
}

# cat("transects_sf geometry type(s):\n"); print(unique(sf::st_geometry_type(transects_sf)))
# cat("transects_sf sfc class: "); print(class(sf::st_geometry(transects_sf)))
# 
# cat("air_sf type(s): ");        print(unique(sf::st_geometry_type(air_sf)))
# cat("muni_sf type(s): ");       print(unique(sf::st_geometry_type(muni_sf)))
# cat("aoi_terre type(s): ");     print(unique(sf::st_geometry_type(aoi_terre)))