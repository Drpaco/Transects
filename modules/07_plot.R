# 07_plot.R
# Final plotting (CRS-safe) + scale-adaptive transect label offset
# Now supports:
#   - out_png NULL-safe (returns ggplot when out_png is NULL)
#   - optional `routes_sf` layer (optimized flight routes)
#   - optional `aoi_sf` outline UNDERLAY

plot_final_map <- function(
    terre_crop, grat_final,
    transects_sf, tran_labels_sf,
    muni_sf, air_sf,
    xlim_terre, ylim_terre,
    bb_txt, out_png,
    terre_crs,
    tran_label_offset_frac = 0.01,
    routes_sf = NULL,
    aoi_sf = NULL
) {
  
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
  
  # ---- Force everything into terre_crs ----
  if (inherits(terre_crop, "sf")) terre_crop <- sf::st_transform(terre_crop, terre_crs)
  if (inherits(grat_final, "sf")) grat_final <- sf::st_transform(grat_final, terre_crs)
  
  if (inherits(transects_sf, "sf") && nrow(transects_sf) > 0) transects_sf <- sf::st_transform(transects_sf, terre_crs)
  if (inherits(tran_labels_sf, "sf") && nrow(tran_labels_sf) > 0) tran_labels_sf <- sf::st_transform(tran_labels_sf, terre_crs)
  
  if (inherits(muni_sf, "sf") && nrow(muni_sf) > 0) muni_sf <- sf::st_transform(muni_sf, terre_crs)
  if (inherits(air_sf, "sf")  && nrow(air_sf)  > 0) air_sf  <- sf::st_transform(air_sf,  terre_crs)
  
  if (inherits(routes_sf, "sf") && nrow(routes_sf) > 0) routes_sf <- sf::st_transform(routes_sf, terre_crs)
  
  # NEW: transform AOI underlay
  if (inherits(aoi_sf, "sf") && nrow(aoi_sf) > 0) aoi_sf <- sf::st_transform(aoi_sf, terre_crs)
  
  # ---- Build label data AFTER CRS align ----
  tran_labels_df <- data.frame(X=numeric(0), Y=numeric(0), label=character(0))
  if (inherits(tran_labels_sf, "sf") && nrow(tran_labels_sf) > 0) {
    if (!("label" %in% names(tran_labels_sf))) tran_labels_sf$label <- paste0("T", tran_labels_sf$id)
    xy <- sf::st_coordinates(tran_labels_sf)
    tran_labels_df <- data.frame(X = xy[,1], Y = xy[,2], label = as.character(tran_labels_sf$label), stringsAsFactors = FALSE)
  }
  
  muni_labels_df <- data.frame(X=numeric(0), Y=numeric(0), label=character(0))
  if (inherits(muni_sf, "sf") && nrow(muni_sf) > 0) {
    nm <- if ("name" %in% names(muni_sf)) as.character(muni_sf$name) else as.character(seq_len(nrow(muni_sf)))
    xy <- sf::st_coordinates(muni_sf)
    muni_labels_df <- data.frame(X = xy[,1], Y = xy[,2], label = nm, stringsAsFactors = FALSE)
  }
  
  air_labels_df <- data.frame(X=numeric(0), Y=numeric(0), label=character(0))
  if (inherits(air_sf, "sf") && nrow(air_sf) > 0) {
    if (!("airport_label" %in% names(air_sf))) {
      air_sf$airport_label <- if ("ident" %in% names(air_sf)) as.character(air_sf$ident) else ""
    }
    xy <- sf::st_coordinates(air_sf)
    air_labels_df <- data.frame(X = xy[,1], Y = xy[,2], label = as.character(air_sf$airport_label), stringsAsFactors = FALSE)
  }
  
  # ---- Scale-adaptive transect label offset (east) ----
  offset_x <- as.numeric(tran_label_offset_frac) * (as.numeric(xlim_terre[2]) - as.numeric(xlim_terre[1]))
  if (!is.finite(offset_x)) offset_x <- 0
  if (nrow(tran_labels_df) > 0) tran_labels_df$X <- tran_labels_df$X + offset_x
  
  # NEW: Validate/repair x/y limits if they look like degrees while we’re in projected CRS
  repair_limits <- function(xlim, ylim) {
    # invalid?
    inv <- any(!is.finite(as.numeric(xlim))) || any(!is.finite(as.numeric(ylim)))
    # suspiciously tiny ranges (e.g., -10..10) when using metres -> likely degrees passed in
    tiny <- (abs(diff(as.numeric(xlim))) <= 10 && abs(diff(as.numeric(ylim))) <= 10)
    # derive combined bbox from terra/lines/points present
    if (inv || tiny) {
      layers <- list(terre_crop, transects_sf, muni_sf, air_sf, routes_sf, aoi_sf)
      layers <- Filter(function(x) inherits(x,"sf") && nrow(x) > 0, layers)
      if (length(layers) > 0) {
        bb <- Reduce(sf::st_union, lapply(layers, sf::st_geometry)) |> sf::st_bbox()
        padx <- as.numeric(bb["xmax"] - bb["xmin"]) * 0.05
        pady <- as.numeric(bb["ymax"] - bb["ymin"]) * 0.05
        xlim <- c(as.numeric(bb["xmin"])-padx, as.numeric(bb["xmax"])+padx)
        ylim <- c(as.numeric(bb["ymin"])-pady, as.numeric(bb["ymax"])+pady)
      }
    }
    list(x=xlim, y=ylim)
  }
  lims <- repair_limits(xlim_terre, ylim_terre)
  xlim_terre <- lims$x
  ylim_terre <- lims$y
  
  # ---- Start plot ----
  p <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = terre_crop, fill = "darkkhaki", color = "grey40") +
    # AOI underlay (outline)
    { if (inherits(aoi_sf, "sf") && nrow(aoi_sf) > 0)
      ggplot2::geom_sf(data = aoi_sf, fill = NA, color = "black", linewidth = 0.7) } +
    ggplot2::geom_sf(data = grat_final, color = "grey70", linewidth = 0.25)
  
  # Optional routes
  if (inherits(routes_sf, "sf") && nrow(routes_sf) > 0) {
    if ("trip" %in% names(routes_sf)) {
      p <- p + ggplot2::geom_sf(data = routes_sf, ggplot2::aes(color = factor(trip)), linewidth = 1.1)
    } else {
      p <- p + ggplot2::geom_sf(data = routes_sf, ggplot2::aes(color = "Routes"), linewidth = 1.1)
    }
  }
  
  # Transects + labels
  if (inherits(transects_sf, "sf") && nrow(transects_sf) > 0) {
    p <- p + ggplot2::geom_sf(data = transects_sf, ggplot2::aes(color = "Transects"), linewidth = 1.1)
  }
  if (nrow(tran_labels_df) > 0) {
    p <- p + ggplot2::geom_text(
      data = tran_labels_df,
      ggplot2::aes(x = X, y = Y, label = label),
      color = "red4", size = 3, fontface = "bold",
      hjust = 0, inherit.aes = FALSE
    )
  }
  
  # Municipalities + labels
  if (inherits(muni_sf, "sf") && nrow(muni_sf) > 0) {
    p <- p + ggplot2::geom_sf(data = muni_sf, ggplot2::aes(color = "Municipalities"), size = 1.6)
  }
  if (nrow(muni_labels_df) > 0) {
    if (has_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data = muni_labels_df,
        ggplot2::aes(x = X, y = Y, label = label),
        color = "red4", size = 3,
        min.segment.length = 0,
        box.padding = 0.25, point.padding = 0.15,
        max.overlaps = Inf, inherit.aes = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = muni_labels_df,
        ggplot2::aes(x = X, y = Y, label = label),
        color = "red4", size = 3, inherit.aes = FALSE
      )
    }
  }
  
  # Airports + labels
  if (inherits(air_sf, "sf") && nrow(air_sf) > 0) {
    p <- p + ggplot2::geom_sf(data = air_sf, ggplot2::aes(color = "Airports"), size = 1.3)
  }
  if (nrow(air_labels_df) > 0) {
    if (has_ggrepel) {
      p <- p + ggrepel::geom_text_repel(
        data = air_labels_df,
        ggplot2::aes(x = X, y = Y, label = label),
        color = "blue4", size = 2.8,
        min.segment.length = 0,
        box.padding = 0.25, point.padding = 0.15,
        max.overlaps = Inf, inherit.aes = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = air_labels_df,
        ggplot2::aes(x = X, y = Y, label = label),
        color = "blue4", size = 2.8, inherit.aes = FALSE
      )
    }
  }
  
  # Legend + extent + CRS lock
  legend_values <- c(); legend_colors <- c()
  if (inherits(muni_sf, "sf") && nrow(muni_sf) > 0) { legend_values <- c(legend_values,"Municipalities"); legend_colors <- c(legend_colors,"red") }
  if (inherits(air_sf, "sf")  && nrow(air_sf)  > 0) { legend_values <- c(legend_values,"Airports");       legend_colors <- c(legend_colors,"blue") }
  if (inherits(transects_sf, "sf") && nrow(transects_sf) > 0) { legend_values <- c(legend_values,"Transects"); legend_colors <- c(legend_colors,"red3") }
  if (inherits(routes_sf, "sf") && nrow(routes_sf) > 0 && !("trip" %in% names(routes_sf))) {
    legend_values <- c(legend_values, "Routes"); legend_colors <- c(legend_colors, "dodgerblue3")
  }
  
  p <- p +
    ggplot2::coord_sf(crs = terre_crs, xlim = xlim_terre, ylim = ylim_terre, expand = FALSE, clip = "off") +
    ggplot2::labs(title = bb_txt) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.height = grid::unit(0.18, "lines"),
      legend.key.width  = grid::unit(0.65, "lines"),
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(5.5, 35, 2, 5.5)
    )
  
  if (length(legend_values) > 0) {
    names(legend_colors) <- legend_values
    p <- p + ggplot2::scale_color_manual(
      name = NULL, values = legend_colors, limits = legend_values, drop = FALSE
    )
  } else if (inherits(routes_sf, "sf") && "trip" %in% names(routes_sf)) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Trip"))
  }
  
  # Save or return
  if (!is.null(out_png) && nzchar(out_png)) {
    dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)
    if (interactive()) print(p)
    ggplot2::ggsave(filename = out_png, plot = p, width = 10, height = 8, dpi = 200)
    return(invisible(p))
  } else {
    if (interactive()) print(p)
    return(p)
  }
}