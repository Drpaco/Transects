setwd('U:/SOMEC/Aérien/InventaireParPhotosAeriennes/transectsPhoto')

library(sf)
library(dplyr)
library(ggplot2)
library(httr)
library(readr)
library(grid)   # unit()
options(timeout = 600)

# Avoid s2 polygon validity failures from terre
sf::sf_use_s2(FALSE)

# ----------------------------
# Paths
# ----------------------------
terre_path <- "U:/SOMEC/Aérien/InventaireParPhotosAeriennes/transectsPhoto/data_downloads/background/terreNAD83.shp"
dir.create("data_downloads", showWarnings = FALSE, recursive = TRUE)
air_csv <- "data_downloads/airports.csv"

# ----------------------------
# Helpers
# ----------------------------
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

bbox_filter_points <- function(sf_points, bb_ll) {
  bb_ll <- validate_bbox(bb_ll)
  pts_ll <- st_transform(sf_points, 4326)
  xy <- st_coordinates(pts_ll)
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

bbox_span_km <- function(bb_ll) {
  b <- validate_bbox(bb_ll)
  midlat <- as.numeric((b["ymin"] + b["ymax"]) / 2)
  dx_km <- as.numeric(b["xmax"] - b["xmin"]) * 111 * cos(midlat * pi / 180)
  dy_km <- as.numeric(b["ymax"] - b["ymin"]) * 111
  c(width_km = dx_km, height_km = dy_km)
}

# ----------------------------
# NEW: Transect label points at each transect's EAST end
# (rightmost end of each line), nudged east by offset_m so label sits beside the line.
# Works even if transects are MULTILINESTRING after clipping.
# ----------------------------
transect_labels_at_east_end <- function(transects_sf, bb_ll, offset_m = 300, utm_epsg = NULL) {
  if (!inherits(transects_sf, "sf") || nrow(transects_sf) == 0) {
    return(st_sf(id = integer(0), geometry = st_sfc(crs = st_crs(transects_sf))))
  }
  
  bb_ll <- validate_bbox(bb_ll)
  
  # Choose UTM zone from bbox center
  if (is.null(utm_epsg)) {
    lon0 <- as.numeric((bb_ll["xmin"] + bb_ll["xmax"]) / 2)
    zone <- floor((lon0 + 180) / 6) + 1
    utm_epsg <- 32600 + zone
  }
  
  # Ensure id exists
  if (!("id" %in% names(transects_sf))) transects_sf$id <- seq_len(nrow(transects_sf))
  
  tr_utm <- st_transform(transects_sf, utm_epsg)
  
  label_pts <- vector("list", nrow(tr_utm))
  
  for (i in seq_len(nrow(tr_utm))) {
    g <- st_geometry(tr_utm[i, ])
    
    # flatten to LINESTRING parts
    parts <- st_cast(g, "LINESTRING", warn = FALSE)
    coords <- do.call(rbind, lapply(parts, st_coordinates))
    
    # eastmost vertex
    idx <- which.max(coords[, "X"])
    x_end <- coords[idx, "X"]
    y_end <- coords[idx, "Y"]
    
    # nudge east for label placement
    label_pts[[i]] <- st_point(c(x_end + offset_m, y_end))
  }
  
  lbl_utm <- st_sf(
    id = tr_utm$id,
    geometry = st_sfc(label_pts, crs = utm_epsg)
  )
  
  st_transform(lbl_utm, st_crs(transects_sf))
}

# ----------------------------
# Read terre background
# ----------------------------
terre <- st_read(terre_path, quiet = TRUE)
terre_crs <- st_crs(terre)

# ----------------------------
# STEP 0: Landmark search (NRCan GeoNames API) + bbox refine
# GeoNames supports q/province/bbox/num/skip and GeoJSON output. 
# ----------------------------
geonames_search <- function(q, province = NULL, num = 20, lang = "en") {
  base <- sprintf("https://geogratis.gc.ca/services/geoname/%s/geonames.geojson", lang)  # 
  query <- list(q = q, num = num, exclude = "links")
  if (!is.null(province)) query$province <- province
  
  r <- GET(base, query = query, timeout(60))
  stop_for_status(r)
  txt <- content(r, as = "text", encoding = "UTF-8")
  st_read(txt, quiet = TRUE)
}

choose_geoname <- function(pts_sf, n_show = 12) {
  if (nrow(pts_sf) == 0) stop("No results returned. Try another term.")
  cols <- intersect(names(pts_sf), c("name","concise","province","latitude","longitude","id","key"))
  if (length(cols) == 0) cols <- names(pts_sf)[1:min(8, ncol(pts_sf))]
  print(pts_sf[1:min(n_show, nrow(pts_sf)), cols])
  
  idx <- as.integer(get_input_safe(sprintf("Select row (1..%d): ", min(n_show, nrow(pts_sf))), "1"))
  if (is.na(idx) || idx < 1 || idx > min(n_show, nrow(pts_sf))) stop("Invalid selection.")
  pts_sf[idx, ]
}

bbox_from_point_km <- function(pt_sf, buffer_km = 25) {
  pt_ll <- st_transform(pt_sf, 4326)
  lon <- st_coordinates(pt_ll)[1,1]
  zone <- floor((lon + 180) / 6) + 1
  utm <- 32600 + zone
  
  buf <- st_buffer(st_transform(pt_ll, utm), buffer_km * 1000)
  st_bbox(st_transform(st_as_sfc(st_bbox(buf)), 4326))
}

plot_bbox_preview <- function(bb_ll, pt_sf = NULL, grid_step_deg = 0.01) {
  bb_ll <- validate_bbox(bb_ll)
  
  terre_ll <- st_transform(terre, 4326)
  bb_poly <- st_as_sfc(bb_ll) |> st_set_crs(4326)
  
  grat <- st_graticule(
    bbox = bb_ll, crs = st_crs(4326),
    ndx = max(2, round((bb_ll["xmax"] - bb_ll["xmin"]) / grid_step_deg)),
    ndy = max(2, round((bb_ll["ymax"] - bb_ll["ymin"]) / grid_step_deg))
  )
  
  ggplot() +
    geom_sf(data = terre_ll, fill = "grey97", color = "grey75", linewidth = 0.2) +
    geom_sf(data = st_sf(geometry = grat), color = "grey88", linewidth = 0.2) +
    geom_sf(data = st_sf(geometry = bb_poly), fill = NA, color = "magenta", linewidth = 1) +
    { if (!is.null(pt_sf)) geom_sf(data = st_transform(pt_sf, 4326), color = "blue", size = 3) } +
    coord_sf(xlim = c(bb_ll["xmin"], bb_ll["xmax"]),
             ylim = c(bb_ll["ymin"], bb_ll["ymax"]),
             expand = FALSE) +
    theme_minimal() +
    labs(title = sprintf("BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f",
                         bb_ll["xmin"], bb_ll["ymin"], bb_ll["xmax"], bb_ll["ymax"]))
}

refine_bbox <- function(pt_sf, initial_buffer_km = 25, grid_step_deg = 0.01) {
  bb <- validate_bbox(bbox_from_point_km(pt_sf, initial_buffer_km))
  
  zoom_out_step <- 1.25
  zoom_in_step  <- 0.80
  pan_step      <- 0.20
  
  repeat {
    span <- bbox_span_km(bb)
    cat(sprintf("\nCurrent window ~ %.1f km (W–E) x %.1f km (S–N)\n", span["width_km"], span["height_km"]))
    print(plot_bbox_preview(bb, pt_sf, grid_step_deg))
    
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
        bb <- validate_bbox(st_bbox(c(xmin = vals[1], ymin = vals[2], xmax = vals[3], ymax = vals[4]),
                                    crs = st_crs(4326)))
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
      
      bb <- validate_bbox(st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), crs = st_crs(4326)))
      next
    }
    
    if (ans2 %in% c("+", "out")) {
      z <- zoom_out_step
    } else if (ans2 %in% c("-", "in")) {
      z <- zoom_in_step
    } else {
      z <- suppressWarnings(as.numeric(ans2))
    }
    
    if (is.finite(z) && z > 0) {
      xmin <- as.numeric(bb["xmin"]); xmax <- as.numeric(bb["xmax"])
      ymin <- as.numeric(bb["ymin"]); ymax <- as.numeric(bb["ymax"])
      cx <- (xmin + xmax) / 2
      cy <- (ymin + ymax) / 2
      w <- (xmax - xmin) * z / 2
      h <- (ymax - ymin) * z / 2
      bb <- validate_bbox(st_bbox(c(xmin = cx - w, ymin = cy - h, xmax = cx + w, ymax = cy + h),
                                  crs = st_crs(4326)))
    } else {
      message("Unrecognized command.")
    }
  }
  bb
}

# --- load transect module (must define refine_transects / make_transects_ll / etc) ---
source('TransectDesignerModule.R', encoding = "UTF-8")

cat("STEP 0: Landmark search and bbox selection\n")
landmark <- get_input_safe("Landmark search (e.g., Perce, Gaspe, Anticosti): ", "Perce")
province_code <- get_input_safe("Province code (Enter to skip; Quebec=24): ", "24")
province_code <- trimws(province_code); if (province_code == "") province_code <- NULL

cand <- geonames_search(landmark, province = province_code, num = 20)
if (nrow(cand) == 0 && !is.null(province_code)) cand <- geonames_search(landmark, province = NULL, num = 20)
pt <- choose_geoname(cand)

bb_final <- refine_bbox(pt, initial_buffer_km = 25, grid_step_deg = 0.01)

span <- bbox_span_km(bb_final)
bb_txt <- sprintf("BBox: xmin=%.4f ymin=%.4f xmax=%.4f ymax=%.4f  (~%.1f km x %.1f km)",
                  bb_final["xmin"], bb_final["ymin"], bb_final["xmax"], bb_final["ymax"],
                  span["width_km"], span["height_km"])

# ----------------------------
# STEP 1: Airports from OurAirports (cached) 
# ----------------------------
if (!file.exists(air_csv)) {
  download.file("https://ourairports.com/data/airports.csv", air_csv, mode = "wb")  # 
}
air <- read_csv(air_csv, show_col_types = FALSE) |>
  filter(iso_country == "CA", !is.na(longitude_deg), !is.na(latitude_deg))

air_sf <- st_as_sf(air, coords = c("longitude_deg", "latitude_deg"), crs = 4326)
air_sf <- bbox_filter_points(air_sf, bb_final)

if (nrow(air_sf) > 0) {
  air_sf$airport_label <- with(air_sf, {
    code <- ifelse(!is.na(iata_code) & iata_code != "", iata_code, ident)
    ifelse(!is.na(name) & name != "", paste0(code, " – ", name), code)
  })
} else {
  air_sf$airport_label <- character(0)
}
air_sf <- st_transform(air_sf, terre_crs)

# ----------------------------
# STEP 2: Main municipalities (NRCan GeoNames bbox query) 
# ----------------------------
geonames_bbox_paged <- function(bb_ll, province = "24", concise = "CITY", theme = NULL,
                                lang = "en", num = 1000) {
  bb_ll <- validate_bbox(bb_ll)
  base <- sprintf("https://geogratis.gc.ca/services/geoname/%s/geonames.geojson", lang)  # 
  bbox_param <- paste(bb_ll["xmin"], bb_ll["ymin"], bb_ll["xmax"], bb_ll["ymax"], sep = ",")
  
  chunks <- list(); skip <- 0
  repeat {
    query <- list(bbox = bbox_param, num = num, skip = skip, exclude = "links")
    if (!is.null(province)) query$province <- province
    if (!is.null(concise))  query$concise  <- concise
    if (!is.null(theme))    query$theme    <- theme
    
    r <- httr::GET(base, query = query, httr::timeout(60))
    httr::stop_for_status(r)
    
    txt <- httr::content(r, as = "text", encoding = "UTF-8")
    chunk <- sf::st_read(txt, quiet = TRUE)
    
    if (nrow(chunk) == 0) break
    chunks[[length(chunks) + 1]] <- chunk
    if (nrow(chunk) < num) break
    skip <- skip + num
  }
  
  if (length(chunks) == 0) return(sf::st_sf(geometry = sf::st_sfc(crs = 4326)))
  do.call(rbind, chunks)
}

muni_raw <- geonames_bbox_paged(bb_final, province = province_code, concise = "CITY")
muni_ll  <- to_points_safe(muni_raw, bb_final)
muni_ll  <- bbox_filter_points(muni_ll, bb_final)
muni_sf  <- st_transform(muni_ll, terre_crs)

# ----------------------------
# Transect designer call (module)
# ----------------------------
transects_ll <- refine_transects(
  bb_ll = bb_final,
  terre = terre,
  pt_landmark_ll = pt,
  air_sf = air_sf,
  muni_sf = muni_sf,
  n = 17,
  len_km = 20,
  sp_km = 2,
  brg_deg = 90,
  clip_to_bbox = TRUE,
  show_full_terre = FALSE
)
transects_sf <- st_transform(transects_ll, terre_crs)

# ----------------------------
# STEP 3: Crop terre and define plot limits
# ----------------------------
bb_terre <- st_bbox(st_transform(st_as_sfc(bb_final), terre_crs))
terre_crop <- st_crop(terre, bb_terre)

xlim_terre <- c(as.numeric(bb_terre["xmin"]), as.numeric(bb_terre["xmax"]))
ylim_terre <- c(as.numeric(bb_terre["ymin"]), as.numeric(bb_terre["ymax"]))

grat_final <- st_graticule(bbox = bb_final, crs = st_crs(4326), ndx = 8, ndy = 8)
grat_final <- st_transform(st_sf(geometry = grat_final), terre_crs)

# --- Transect labels BESIDE each transect at its own EAST end ---
tran_labels_sf <- transect_labels_at_east_end(transects_sf, bb_final, offset_m = 2000)
tran_labels_sf$label <- paste0("T", tran_labels_sf$id)
tran_labels_df <- cbind(st_coordinates(tran_labels_sf), st_drop_geometry(tran_labels_sf))

# ----------------------------
# STEP 4: Final plot
# ----------------------------
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(ggrepel)

p <- ggplot() +
  geom_sf(data = terre_crop, fill = "darkkhaki", color = "grey40") +
  geom_sf(data = grat_final, color = "grey70", linewidth = 0.25)

# Transects + labels beside east end
if (inherits(transects_sf, "sf") && nrow(transects_sf) > 0) {
  p <- p +
    geom_sf(data = transects_sf, aes(color = "Transects"), linewidth = 1.1) +
    geom_text(
      data = tran_labels_df,
      aes(x = X, y = Y, label = label),
      color = "red4",
      size = 3,
      fontface = "bold"
    )
}

# Municipalities
if (inherits(muni_sf, "sf") && nrow(muni_sf) > 0) {
  p <- p + geom_sf(data = muni_sf, aes(color = "Municipalities"), size = 1.6)
  
  name_field <- if ("name" %in% names(muni_sf)) "name" else names(muni_sf)[1]
  p <- p + ggrepel::geom_text_repel(
    data = cbind(st_coordinates(muni_sf), st_drop_geometry(muni_sf)),
    aes(x = X, y = Y, label = .data[[name_field]]),
    color = "red4", size = 3, min.segment.length = 0,
    box.padding = 0.25, point.padding = 0.15, max.overlaps = Inf
  )
}

# Airports
if (inherits(air_sf, "sf") && nrow(air_sf) > 0) {
  p <- p + geom_sf(data = air_sf, aes(color = "Airports"), size = 1.3)
  
  p <- p + ggrepel::geom_text_repel(
    data = cbind(st_coordinates(air_sf), st_drop_geometry(air_sf)),
    aes(x = X, y = Y, label = airport_label),
    color = "blue4", size = 2.8, min.segment.length = 0,
    box.padding = 0.25, point.padding = 0.15, max.overlaps = Inf
  )
}

p <- p +
  scale_color_manual(
    name = NULL,
    values = c("Municipalities" = "red",
               "Airports" = "blue",
               "Transects" = "red3")
  ) +
  coord_sf(xlim = xlim_terre, ylim = ylim_terre, expand = FALSE, clip = "off") +
  labs(title = bb_txt) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.height = unit(0.18, "lines"),
    legend.key.width  = unit(0.65, "lines"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 10),
    plot.margin = margin(5.5, 25, 2, 5.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), nrow = 1))

print(p)
ggsave("data_downloads/map_bbox_preview.png", p, width = 10, height = 8, dpi = 200)
cat("\nDone. Map saved to data_downloads/map_bbox_preview.png\n")

# ----------------------------
# FINAL STEP: Save transects as ESRI Shapefile (name of your choosing)
# ----------------------------

# Ensure transects exist
if (exists("transects_sf") && inherits(transects_sf, "sf") && nrow(transects_sf) > 0) {
  
  # Ensure simple attributes for shapefile (field names <= 10 chars)
  if (!("id" %in% names(transects_sf))) transects_sf$id <- seq_len(nrow(transects_sf))
  transects_sf$label <- paste0("T", transects_sf$id)  # T1, T2, ...
  
  # Ask user for output name (without extension is fine)
  out_name <- get_input_safe("Enter transect shapefile name (e.g., transects_Perce_2026): ", "transects")
  out_name <- gsub("[^A-Za-z0-9_\\-]", "_", out_name)  # sanitize a bit
  
  # Choose output folder (you can change this)
  out_dir <- file.path("data_downloads", "transects_export")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  shp_path <- file.path(out_dir, paste0(out_name, ".shp"))
  
  # If you want to overwrite without prompts:
  if (file.exists(shp_path)) {
    # Delete all companion files if already present
    base <- tools::file_path_sans_ext(shp_path)
    unlink(paste0(base, c(".shp",".shx",".dbf",".prj",".cpg",".sbn",".sbx")), force = TRUE)
  }
  
  # Write shapefile
  st_write(transects_sf[, c("id","label","geometry")],
           shp_path, driver = "ESRI Shapefile",
           delete_layer = TRUE, quiet = TRUE)
  
  cat("\nTransects saved to:\n", shp_path, "\n", sep = "")
  
} else {
  cat("\nNo transects to save (transects_sf is missing or empty).\n")
}