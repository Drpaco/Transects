# 05_airports.R
# Airports from OurAirports (cached CSV)

# 05_airports.R
load_airports <- function(csv_path) {
  if (!file.exists(csv_path)) {
    download.file("https://ourairports.com/data/airports.csv", csv_path, mode = "wb")
  }
  
  air <- readr::read_csv(csv_path, show_col_types = FALSE) |>
    dplyr::filter(
      iso_country == "CA",
      !is.na(longitude_deg), !is.na(latitude_deg)
    ) |>
    # Normalize common fields for matching
    dplyr::mutate(
      icao = ident,             # OurAirports 'ident' is ICAO for CA (e.g., CYGP)
      iata = iata_code,         # alias
      local = local_code,       # alias
      code  = dplyr::if_else(!is.na(iata_code) & iata_code != "", iata_code, ident),
      airport_label = dplyr::if_else(
        !is.na(name) & name != "", paste0(code, " \u2013 ", name), code
      )
    )
  
  sf::st_as_sf(air, coords = c("longitude_deg", "latitude_deg"), crs = 4326)
}

make_airport_labels <- function(air_sf) {
  if (!inherits(air_sf, "sf")) stop("make_airport_labels() expects sf")
  if (nrow(air_sf) == 0) {
    air_sf$airport_label <- character(0)
    return(air_sf)
  }
  # Keep what load_airports produced; fill if missing
  if (!"airport_label" %in% names(air_sf)) {
    air_sf$airport_label <- with(air_sf, {
      code <- ifelse(!is.na(iata) & iata != "", iata,
                     ifelse(!is.na(ident) & ident != "", ident, ""))
      ifelse(!is.na(name) & name != "", paste0(code, " \u2013 ", name), code)
    })
  }
  air_sf
}