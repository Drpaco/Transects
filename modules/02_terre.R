load_terre <- function(path) {
  terre <- sf::st_read(path, quiet = TRUE)
  list(terre = terre, crs = sf::st_crs(terre))
}