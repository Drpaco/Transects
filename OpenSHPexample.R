# Install once if needed
#install.packages(c("sf", "ggplot2"))   # optional: "tmap", "dplyr"

# Load libraries
library(sf)
library(ggplot2)
# library(tmap)
# library(dplyr)

# Example: full path to the .shp file
shp_path <- "F:/U/Oiseaux souillés/Inventaires pélagiques/Golfe/Zones2/BaieChaleurs.shp"
# On macOS/Linux it might look like:
# shp_path <- "/Users/you/data/boundaries/municipal_boundaries.shp"

# Read the shapefile
# stringsAsFactors = FALSE prevents automatic factor conversion of attributes
# quiet = TRUE silences verbose output
g <- st_read(shp_path, stringsAsFactors = FALSE, quiet = TRUE)

# Inspect
print(g)           # basic info
st_geometry_type(g) # geometry types (POINT, LINESTRING, POLYGON, MULTIPOLYGON, etc.)
st_crs(g)          # coordinate reference system
head(g)            # attribute preview


plot(st_geometry(g), col = "lightblue", border = "gray40", main = "Shapefile preview")
