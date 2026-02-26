# ==============================================================================
# 00_config.R
# PROJECT CONFIGURATION
# ------------------------------------------------------------------------------

# EDIT PATHS BELOW TO MATCH YOUR LOCAL FILE LOCATIONS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1) Terre land layer (NAD83)
#    Your Terre shapefile lives in: data_downloads/background/TerreNAD83.shp
#    Keep it there — you do NOT need to move the file.
# ------------------------------------------------------------------------------
TERRE_PATH <- "data_downloads/background/TerreNAD83.shp"

# ------------------------------------------------------------------------------
# 2) OurAirports CSV cache
#    The Shiny app will download airports.csv into this path if missing.
# ------------------------------------------------------------------------------
AIRPORTS_CSV <- "data_downloads/airports.csv"

# ------------------------------------------------------------------------------
# 3) Export directory
#    All transects, optimizer outputs (shapefile/GPKG/CSV/PNG) go here.
# ------------------------------------------------------------------------------
TRANSECT_EXPORT_DIR <- "exports"

# ==============================================================================
# CREATE FOLDERS IF THEY DON'T EXIST
# ==============================================================================
dir.create(dirname(AIRPORTS_CSV), showWarnings = FALSE, recursive = TRUE)
dir.create(TRANSECT_EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# END OF CONFIG
# ==============================================================================