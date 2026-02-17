# Global config
options(timeout = 600)

# Avoid s2 issues with invalid polygons in terre
sf::sf_use_s2(FALSE)

# Paths (edit once here)
ROOT_DIR <- "U:/SOMEC/Aérien/InventaireParPhotosAeriennes/transectsPhoto"

TERRE_PATH <- file.path(ROOT_DIR, "data_downloads/background/terreNAD83.shp")

CACHE_DIR <- file.path(ROOT_DIR, "data_downloads")
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

AIRPORTS_CSV <- file.path(CACHE_DIR, "airports.csv")

OUTPUT_DIR <- file.path(ROOT_DIR, "data_downloads/outputs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_PNG <- file.path(OUTPUT_DIR, "map_bbox_preview.png")

TRANSECT_EXPORT_DIR <- file.path(ROOT_DIR, "data_downloads/transects_export")
dir.create(TRANSECT_EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)