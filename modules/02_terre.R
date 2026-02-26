# modules/02_terre.R
# Minimal loader for Terre land layer (NAD83), returning both layer and CRS

load_terre <- function(path, target_epsg = 4269) {
  # Accept either a file (.gpkg / .shp) OR a folder that contains one
  if (dir.exists(path)) {
    cand <- list.files(path, pattern = "\\.(gpkg|shp)$", full.names = TRUE, ignore.case = TRUE)
    if (!length(cand)) {
      stop("No .gpkg or .shp found in directory: ", path)
    }
    f <- cand[1]
  } else {
    f <- path
  }
  
  # Read the vector layer
  x <- suppressWarnings(sf::st_read(f, quiet = TRUE))
  if (!inherits(x, "sf")) {
    stop("load_terre(): could not read sf from: ", f)
  }
  
  # Ensure CRS is present; Terre should be in NAD83. If missing, stop (safer).
  crs_in <- sf::st_crs(x)
  if (is.na(crs_in)) {
    stop("Terre layer has no CRS. Please ensure TerreNAD83.shp has a defined CRS (EPSG:4269).")
  }
  
  # Transform to NAD83 (EPSG:4269) for consistent display if needed
  if (!identical(crs_in, sf::st_crs(target_epsg))) {
    x <- sf::st_transform(x, target_epsg)
    crs_in <- sf::st_crs(target_epsg)
  }
  
  list(terre = x, crs = crs_in)
}