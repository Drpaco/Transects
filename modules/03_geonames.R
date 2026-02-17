# NRCan GeoNames API supports q/bbox/province/concise/paging; GeoJSON output.
geonames_search <- function(q, province = NULL, num = 20, lang = "en") {
  base <- sprintf("https://geogratis.gc.ca/services/geoname/%s/geonames.geojson", lang)  # 
  query <- list(q = q, num = num, exclude = "links")
  if (!is.null(province)) query$province <- province
  
  r <- httr::GET(base, query = query, httr::timeout(60))
  httr::stop_for_status(r)
  txt <- httr::content(r, as = "text", encoding = "UTF-8")
  sf::st_read(txt, quiet = TRUE)
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

geonames_bbox_paged <- function(bb_ll, province = "24", concise = "CITY", theme = NULL,
                                lang = "en", num = 1000) {
  bb_ll <- validate_bbox(bb_ll)
  base <- sprintf("https://geogratis.gc.ca/services/geoname/%s/geonames.geojson", lang)  # 
  bbox_param <- paste(bb_ll["xmin"], bb_ll["ymin"], bb_ll["xmax"], bb_ll["ymax"], sep = ",")
  
  chunks <- list()
  skip <- 0
  
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