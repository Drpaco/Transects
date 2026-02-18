# 08_transect_optimizer.R
# Optional post-processing: optimize flight routes to cover transects within
# endurance limits, starting/ending at a chosen airport, with optional multi-base mode.

# Assumes these packages are already loaded by your launcher: sf, dplyr, ggplot2, grid

# ------------------------------
# Helpers
# ------------------------------
ensure_id <- function(lines_sf) {
  if (!("id" %in% names(lines_sf))) {
    lines_sf$id <- seq_len(nrow(lines_sf))
  }
  lines_sf
}

utm_epsg_from_lonlat <- function(lon, lat) {
  zone <- floor((lon + 180) / 6) + 1
  if (lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
}

# Distance in meters between two XY points
pt_dist <- function(x1,y1,x2,y2) sqrt((x2-x1)^2 + (y2-y1)^2)

# Pick nearest airport to transects centroid using a local projected CRS
nearest_airport <- function(transects_sf, air_sf) {
  if (nrow(air_sf) == 0) {
    stop("No airports available to compute nearest. Pass air_sf_all to the optimizer.")
  }
  tr_ll   <- sf::st_transform(transects_sf, 4326)
  bbox_ll <- sf::st_bbox(tr_ll)
  mid_lon <- as.numeric((bbox_ll["xmin"] + bbox_ll["xmax"]) / 2)
  mid_lat <- as.numeric((bbox_ll["ymin"] + bbox_ll["ymax"]) / 2)
  zone    <- floor((mid_lon + 180) / 6) + 1
  crs_m   <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
  
  tr_m  <- sf::st_transform(transects_sf, crs_m)
  air_m <- sf::st_transform(air_sf,      crs_m)
  
  center_m <- sf::st_centroid(sf::st_union(sf::st_geometry(tr_m)))
  d <- sf::st_distance(sf::st_geometry(air_m), center_m)
  idx <- which.min(as.numeric(d))
  air_sf[idx, , drop = FALSE]
}

# Robust airport chooser: fast exact ident -> exact other fields -> partial -> nearest
choose_airport <- function(air_sf, transects_sf, query = NULL) {
  if (!inherits(air_sf, "sf") || nrow(air_sf) == 0) {
    stop("No airports available. Pass an unfiltered airport layer (e.g., air_sf_all) to run_transect_optimizer.")
  }
  # If user didn’t enter anything, go straight to nearest
  if (is.null(query) || !nzchar(query)) return(nearest_airport(transects_sf, air_sf))
  
  # Normalize the query once
  key <- toupper(trimws(as.character(query)))
  key <- iconv(key, from = "", to = "ASCII//TRANSLIT")
  key <- gsub("\\s+", "", key)
  
  norm_col <- function(x) {
    v <- toupper(trimws(as.character(x)))
    v <- iconv(v, from = "", to = "ASCII//TRANSLIT")
    gsub("\\s+", "", v)
  }
  
  # ---- 0) FAST PATH: exact match on ident first (this is your ICAO in CA, e.g., CYGP) ----
  if ("ident" %in% names(air_sf)) {
    v <- norm_col(air_sf$ident)
    hit <- which(v == key)
    if (length(hit) > 0) return(air_sf[hit[1], , drop = FALSE])
  }
  
  # ---- 1) exact match across other code fields ----
  cols_priority_exact <- intersect(
    c("icao_code","gps_code","icao","iata","iata_code","local_code","local","code"),
    names(air_sf)
  )
  for (cn in cols_priority_exact) {
    v <- norm_col(air_sf[[cn]])
    hit <- which(v == key)
    if (length(hit) > 0) return(air_sf[hit[1], , drop = FALSE])
  }
  
  # ---- 2) partial match across common label/name fields ----
  cols_partial <- intersect(
    c("ident","icao_code","gps_code","icao","iata","iata_code","local_code","local","code","name","airport_label"),
    names(air_sf)
  )
  if (length(cols_partial) > 0) {
    mask <- Reduce(`|`, lapply(cols_partial, function(cn) {
      vv <- toupper(trimws(as.character(air_sf[[cn]])))
      grepl(key, vv, fixed = TRUE)
    }))
    cand <- air_sf[which(mask), , drop = FALSE]
    if (nrow(cand) >= 1) return(cand[1, , drop = FALSE])
  }
  
  message("No airport matched '", query, "'. Falling back to nearest airport.")
  nearest_airport(transects_sf, air_sf)
}

# --- Helper: choose a stable projected CRS from bbox center (UTM) ---
crs_m_from_bbox <- function(x_sf) {
  x_ll  <- sf::st_transform(x_sf, 4326)
  bb    <- sf::st_bbox(x_ll)
  mid_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  mid_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  zone  <- floor((mid_lon + 180) / 6) + 1
  if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
}

# --- Helper: readable label for an airport row ---
airport_label <- function(ap) {
  for (cn in c("ident","icao","iata","code","name","airport_label","label")) {
    if (cn %in% names(ap)) {
      val <- as.character(ap[[cn]][1])
      if (!is.na(val) && nzchar(val)) return(val)
    }
  }
  "AP"
}

# --- Build routes using per-trip (start_ap, end_ap) selection ---
build_routes_sf_per_trip_airports <- function(lines_sf, plan, crs_m, target_crs, trip_airports) {
  lines_m <- sf::st_transform(ensure_id(lines_sf), crs_m)
  
  get_line_xy <- function(id) {
    geom <- lines_m$geometry[which(lines_m$id == id)][[1]]
    sf::st_coordinates(geom)[, c("X","Y"), drop = FALSE]
  }
  maybe_reverse <- function(xy, reverse = FALSE) if (reverse) xy[nrow(xy):1, , drop = FALSE] else xy
  seg_xy <- function(x1,y1,x2,y2) matrix(c(x1,y1,x2,y2), ncol=2, byrow=TRUE)
  
  routes_sfc <- sf::st_sfc(
    lapply(seq_along(plan$trips), function(t) {
      ap_start_ll <- trip_airports$start[[t]]
      ap_end_ll   <- trip_airports$end[[t]]
      ap_start_m  <- sf::st_transform(ap_start_ll, crs_m)
      ap_end_m    <- sf::st_transform(ap_end_ll,   crs_m)
      as_xy <- sf::st_coordinates(ap_start_m)[1, ]
      ae_xy <- sf::st_coordinates(ap_end_m)[1, ]
      ax <- as_xy["X"]; ay <- as_xy["Y"]
      ex <- ae_xy["X"]; ey <- ae_xy["Y"]
      
      trip_ids <- plan$trips[[t]]
      trip_dir <- plan$trip_dirs_end_endpoint[[t]]
      
      route_xy <- matrix(c(ax, ay), ncol=2)
      for (k in seq_along(trip_ids)) {
        id <- trip_ids[k]
        rev_flag <- (trip_dir[k] == 1)   # 1: end=A -> traverse B->A
        xy_line <- get_line_xy(id)
        xy_line <- maybe_reverse(xy_line, rev_flag)
        
        last <- route_xy[nrow(route_xy), ]
        sx <- xy_line[1,1]; sy <- xy_line[1,2]
        if (last[1] != sx || last[2] != sy) {
          route_xy <- rbind(route_xy, seg_xy(last[1], last[2], sx, sy)[2, , drop=FALSE])
        }
        route_xy <- rbind(route_xy, xy_line[-1, , drop=FALSE])
      }
      
      last <- route_xy[nrow(route_xy), ]
      if (last[1] != ex || last[2] != ey) {
        route_xy <- rbind(route_xy, seg_xy(last[1], last[2], ex, ey)[2, , drop=FALSE])
      }
      sf::st_linestring(route_xy)
    }),
    crs = crs_m
  )
  
  routes_sf <- sf::st_sf(trip = seq_along(routes_sfc), geometry = routes_sfc)
  routes_sf <- sf::st_transform(routes_sf, target_crs)
  routes_m  <- sf::st_transform(routes_sf, crs_m)
  dist_m    <- as.numeric(sf::st_length(routes_m))
  routes_sf$dist_m  <- dist_m
  routes_sf$dist_nm <- dist_m / 1852
  routes_sf
}

# --- Manifest with inter-trip ferry legs + per-trip totals ---
build_manifest_with_airports <- function(lines_sf, plan, crs_m, speed_knots, trip_airports, dep_labels, arr_labels) {
  m_per_nm <- 1852
  lines_m  <- sf::st_transform(ensure_id(lines_sf), crs_m)
  
  # endpoints A/B in meters for each line id
  get_line_xy <- function(id) {
    geom <- lines_m$geometry[which(lines_m$id == id)][[1]]
    sf::st_coordinates(geom)[, c("X","Y"), drop = FALSE]
  }
  end_xy   <- function(xy, end_is_A) if (end_is_A) xy[1, ] else xy[nrow(xy), ]
  start_xy_for_dir <- function(xy, end_is_A) if (end_is_A) xy[nrow(xy), ] else xy[1, ]
  
  rows <- list()
  
  # add a "ferry" leg from previous trip's ARR AP to current DEP AP (t>1)
  add_ferry_leg <- function(t, from_ap_ll, to_ap_ll) {
    ap_from_m <- sf::st_transform(from_ap_ll, crs_m)
    ap_to_m   <- sf::st_transform(to_ap_ll,   crs_m)
    fxy <- sf::st_coordinates(ap_from_m)[1, ]
    txy <- sf::st_coordinates(ap_to_m)[1, ]
    ferry_m <- sqrt((txy["X"] - fxy["X"])^2 + (txy["Y"] - fxy["Y"])^2)
    ferry_nm <- ferry_m / m_per_nm
    ferry_min <- (ferry_nm / speed_knots) * 60
    data.frame(
      trip = t, leg_seq = 0L, leg_type = "ferry",
      from = airport_label(from_ap_ll), to = airport_label(to_ap_ll),
      transect_id = NA_integer_, direction = NA_character_,
      dist_m = ferry_m, dist_nm = ferry_nm, est_minutes = ferry_min,
      stringsAsFactors = FALSE
    )
  }
  
  for (t in seq_along(plan$trips)) {
    ap_start_ll <- trip_airports$start[[t]]
    ap_end_ll   <- trip_airports$end[[t]]
    
    # inter-trip ferry (t>1): previous ARR -> this DEP
    if (t > 1) {
      rows[[length(rows)+1]] <- add_ferry_leg(t, trip_airports$end[[t-1]], ap_start_ll)
    }
    
    trip_ids <- plan$trips[[t]]
    trip_dir <- plan$trip_dirs_end_endpoint[[t]]
    
    # set current point at start airport
    ap_start_m <- sf::st_transform(ap_start_ll, crs_m)
    cur <- sf::st_coordinates(ap_start_m)[1, ]
    leg_seq <- 1L
    
    for (k in seq_along(trip_ids)) {
      id <- trip_ids[k]
      xy <- get_line_xy(id)
      end_is_A <- (trip_dir[k] == 1)      # 1 => end at A
      sxy <- start_xy_for_dir(xy, end_is_A)
      exy <- end_xy(xy, end_is_A)
      
      # connector: cur -> start of line
      conn_m <- sqrt((sxy["X"] - cur["X"])^2 + (sxy["Y"] - cur["Y"])^2)
      if (conn_m > 0) {
        conn_nm  <- conn_m / m_per_nm
        conn_min <- (conn_nm / speed_knots) * 60
        rows[[length(rows)+1]] <- data.frame(
          trip = t, leg_seq = leg_seq, leg_type = "connector",
          from = if (leg_seq == 1) dep_labels[t] else paste0("T", trip_ids[k-1]),
          to   = paste0("T", id),
          transect_id = NA_integer_, direction = NA_character_,
          dist_m = conn_m, dist_nm = conn_nm, est_minutes = conn_min,
          stringsAsFactors = FALSE
        )
        leg_seq <- leg_seq + 1L
      }
      
      # line (transect) leg
      line_m <- as.numeric(sf::st_length(lines_m[which(lines_m$id == id), ]))
      line_nm <- line_m / m_per_nm
      line_min <- (line_nm / speed_knots) * 60
      rows[[length(rows)+1]] <- data.frame(
        trip = t, leg_seq = leg_seq, leg_type = "transect",
        from = paste0("T", id), to = paste0("T", id),
        transect_id = id, direction = if (end_is_A) "B->A" else "A->B",
        dist_m = line_m, dist_nm = line_nm, est_minutes = line_min,
        stringsAsFactors = FALSE
      )
      leg_seq <- leg_seq + 1L
      
      cur <- exy
    }
    
    # final connector to ARR AP
    ap_end_m <- sf::st_transform(ap_end_ll, crs_m)
    aexy <- sf::st_coordinates(ap_end_m)[1, ]
    back_m <- sqrt((aexy["X"] - cur["X"])^2 + (aexy["Y"] - cur["Y"])^2)
    if (back_m > 0) {
      back_nm  <- back_m / m_per_nm
      back_min <- (back_nm / speed_knots) * 60
      rows[[length(rows)+1]] <- data.frame(
        trip = t, leg_seq = leg_seq, leg_type = "connector",
        from = paste0("T", tail(trip_ids, 1)), to = arr_labels[t],
        transect_id = NA_integer_, direction = NA_character_,
        dist_m = back_m, dist_nm = back_nm, est_minutes = back_min,
        stringsAsFactors = FALSE
      )
    }
  }
  
  manifest <- do.call(rbind, rows)
  
  # add per-trip totals (same value repeated on all rows of the trip)
  totals <- manifest |>
    dplyr::group_by(trip) |>
    dplyr::summarise(trip_total_minutes = sum(est_minutes), .groups = "drop")
  manifest <- dplyr::left_join(manifest, totals, by = "trip")
  manifest
}

# Prepare line endpoints and airports in a metric CRS; supports different arrival airport
prep_lines <- function(lines_sf, airport_start_sf, crs_m, airport_end_sf = NULL) {
  # 1) Project everything to a metric CRS for planning
  stopifnot(inherits(lines_sf, "sf"))
  lines_m       <- sf::st_transform(lines_sf,       crs_m)
  airport_start <- sf::st_transform(airport_start_sf, crs_m)
  airport_end   <- if (is.null(airport_end_sf)) airport_start else sf::st_transform(airport_end_sf, crs_m)
  
  # 2) Normalize geometries:
  #    - make valid
  #    - explode MULTILINESTRING/collections into LINESTRING parts
  #    - drop empties
  lines_m <- sf::st_make_valid(lines_m)
  lines_m <- suppressWarnings(sf::st_cast(lines_m, "LINESTRING"))
  lines_m <- suppressWarnings(sf::st_collection_extract(lines_m, "LINESTRING"))
  lines_m <- lines_m[!sf::st_is_empty(lines_m), , drop = FALSE]
  if (nrow(lines_m) == 0) {
    stop("No usable transect line coordinates found (degenerate geometries).")
  }
  
  # 3) Drop ~0-length pieces (use >= 1 meter)
  len_m <- as.numeric(sf::st_length(lines_m))
  keep  <- is.finite(len_m) & (len_m >= 1)
  lines_m <- lines_m[keep, , drop = FALSE]
  len_m   <- len_m[keep]
  if (nrow(lines_m) == 0) {
    stop("No usable transect line coordinates found (degenerate geometries).")
  }
  
  # 4) Give each retained piece its own id (we're planning per usable segment)
  lines_m$id <- seq_len(nrow(lines_m))
  
  # 5) Extract endpoints A/B safely for each piece
  Ax <- Ay <- Bx <- By <- numeric(nrow(lines_m))
  for (i in seq_len(nrow(lines_m))) {
    coords <- sf::st_coordinates(lines_m[i, ])
    # For safety, if multiple parts slipped through, just take first/last row
    first <- coords[1, c("X", "Y")]
    last  <- coords[nrow(coords), c("X", "Y")]
    Ax[i] <- first[["X"]]; Ay[i] <- first[["Y"]]
    Bx[i] <- last[["X"]];  By[i] <- last[["Y"]]
  }
  
  # 6) Airports XY in meters
  a_start <- sf::st_coordinates(airport_start)[1, ]
  a_end   <- sf::st_coordinates(airport_end)[1,  ]
  
  list(
    ids      = lines_m$id,
    line_len = len_m,
    Ax = Ax, Ay = Ay, Bx = Bx, By = By,
    ax = a_start["X"], ay = a_start["Y"],   # departure
    bx = a_end["X"],   by = a_end["Y"]      # arrival
  )
}

# ------------------------------
# Core optimization logic
# ------------------------------
trip_cost_for_block <- function(dat, order_ids, i, j) {
  block <- order_ids[i:j]
  idx <- match(block, dat$ids)
  m <- length(idx)
  
  end_xy   <- function(ii, e) if (e == 1) c(dat$Ax[ii], dat$Ay[ii]) else c(dat$Bx[ii], dat$By[ii])
  start_xy <- function(ii, e_end) if (e_end == 1) c(dat$Bx[ii], dat$By[ii]) else c(dat$Ax[ii], dat$Ay[ii])
  
  dp  <- matrix(Inf, m, 2)
  prv <- matrix(0L,  m, 2)
  
  # first line: airport_start -> start + line length
  ii <- idx[1]
  for (e in 1:2) {
    s <- start_xy(ii, e)
    dp[1, e] <- pt_dist(dat$ax, dat$ay, s[1], s[2]) + dat$line_len[ii]
  }
  
  # transitions
  if (m >= 2) {
    for (t in 2:m) {
      ii <- idx[t]
      for (e in 1:2) {
        s <- start_xy(ii, e)
        best <- Inf; best_pe <- 0L
        for (pe in 1:2) {
          p <- end_xy(idx[t-1], pe)
          cand <- dp[t-1, pe] + pt_dist(p[1],p[2], s[1],s[2]) + dat$line_len[ii]
          if (cand < best) { best <- cand; best_pe <- pe }
        }
        dp[t, e] <- best
        prv[t, e] <- best_pe
      }
    }
  }
  
  # return to airport_end (may be same as start)
  last <- idx[m]
  best_total <- Inf; best_last_e <- 0L
  for (e in 1:2) {
    p <- end_xy(last, e)
    cand <- dp[m, e] + pt_dist(p[1], p[2], dat$bx, dat$by)
    if (cand < best_total) { best_total <- cand; best_last_e <- e }
  }
  
  # recover directions
  dir_end <- integer(m)
  dir_end[m] <- best_last_e
  if (m >= 2) for (t in (m-1):1) dir_end[t] <- prv[t+1, dir_end[t+1]]
  
  list(cost_m = best_total, end_endpoint = dir_end)
}

plan_trips_dp <- function(dat, order_ids, speed_knots, max_hours) {
  m_per_nm <- 1852
  n <- length(order_ids)
  max_m <- speed_knots * max_hours * m_per_nm
  
  cost <- matrix(Inf, n, n)
  dirs <- vector("list", n*n)
  
  for (i in 1:n) for (j in i:n) {
    tc <- trip_cost_for_block(dat, order_ids, i, j)
    if (tc$cost_m <= max_m) {
      cost[i,j] <- tc$cost_m
      dirs[[ (i-1)*n + j ]] <- tc$end_endpoint
    }
  }
  
  best <- rep(Inf, n+1)
  prev_cut <- rep(NA_integer_, n+1)
  best[1] <- 0
  
  for (k in 2:(n+1)) {
    j <- k-1
    for (i in 1:j) {
      if (is.finite(cost[i,j]) && is.finite(best[i])) {
        cand <- best[i] + cost[i,j]
        if (cand < best[k]) { best[k] <- cand; prev_cut[k] <- i }
      }
    }
  }
  if (!is.finite(best[n+1])) stop("No feasible partition under max_hours. Increase endurance or speed.")
  
  trips <- list(); trip_dirs <- list()
  k <- n+1
  while (k > 1) {
    i <- prev_cut[k]; j <- k-1
    trips[[length(trips)+1]] <- order_ids[i:j]
    trip_dirs[[length(trip_dirs)+1]] <- dirs[[ (i-1)*n + j ]]
    k <- i
  }
  trips <- rev(trips); trip_dirs <- rev(trip_dirs)
  
  total_m <- best[n+1]
  total_nm <- total_m / m_per_nm
  total_h  <- total_nm / speed_knots
  
  list(
    order_ids = order_ids,
    trips = trips,
    trip_dirs_end_endpoint = trip_dirs,
    total_m = total_m,
    total_nm = total_nm,
    total_hours = total_h,
    total_minutes = total_h * 60
  )
}

midpoint_xy <- function(dat) {
  mx <- (dat$Ax + dat$Bx) / 2
  my <- (dat$Ay + dat$By) / 2
  cbind(mx, my)
}

make_mid_dist <- function(dat) {
  M <- midpoint_xy(dat)
  n <- nrow(M)
  D <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    D[i,j] <- pt_dist(M[i,1], M[i,2], M[j,1], M[j,2])
  }
  D
}

nn_order <- function(D, start = 1) {
  n <- nrow(D)
  visited <- rep(FALSE, n)
  ord <- integer(n)
  ord[1] <- start
  visited[start] <- TRUE
  for (k in 2:n) {
    last <- ord[k-1]
    cand <- which(!visited)
    nxt <- cand[which.min(D[last, cand])]
    ord[k] <- nxt
    visited[nxt] <- TRUE
  }
  ord
}

two_opt <- function(ord, score_fun) {
  improved <- TRUE
  best_ord <- ord
  best_score <- score_fun(best_ord)
  while (improved) {
    improved <- FALSE
    n <- length(best_ord)
    for (i in 2:(n-2)) {
      for (k in (i+1):(n-1)) {
        new_ord <- best_ord
        new_ord[i:k] <- rev(new_ord[i:k])
        sc <- score_fun(new_ord)
        if (sc < best_score) {
          best_ord <- new_ord
          best_score <- sc
          improved <- TRUE
        }
      }
    }
  }
  best_ord
}

optimize_order_with_refuel <- function(dat, speed_knots, max_hours, n_starts = 20, verbose = TRUE) {
  Dmid <- make_mid_dist(dat)
  n <- length(dat$ids)
  
  best_plan <- NULL
  best_m <- Inf
  
  starts <- unique(pmin(n, pmax(1, round(seq(1, n, length.out = n_starts)))))
  
  score_fun <- function(ord_idx) {
    order_ids <- dat$ids[ord_idx]
    pl <- plan_trips_dp(dat, order_ids, speed_knots, max_hours)
    pl$total_m
  }
  
  t0 <- Sys.time()
  if (isTRUE(verbose)) {
    message(sprintf("Optimization started: %s  (n=%d, starts=%d)", format(t0), n, length(starts)))
    pb <- utils::txtProgressBar(min = 0, max = length(starts), style = 3)
  }
  
  for (i in seq_along(starts)) {
    s <- starts[i]
    ord0 <- nn_order(Dmid, start = s)
    ord1 <- two_opt(ord0, score_fun)
    order_ids <- dat$ids[ord1]
    pl <- plan_trips_dp(dat, order_ids, speed_knots, max_hours)
    
    if (pl$total_m < best_m) {
      best_m <- pl$total_m
      best_plan <- pl
    }
    
    if (isTRUE(verbose)) {
      utils::setTxtProgressBar(pb, i)
      if (i %% 3 == 0) {
        dt <- difftime(Sys.time(), t0, units = "secs")
        cat(sprintf("  elapsed: %.1f s\r", as.numeric(dt)))
      }
    }
  }
  
  if (isTRUE(verbose)) {
    close(pb)
    dt <- difftime(Sys.time(), t0, units = "secs")
    message(sprintf("\nOptimization finished in %.1f seconds.", as.numeric(dt)))
  }
  
  best_plan
}

# ------------------------------
# Build route linework from the plan (supports different arrival airport)
# ------------------------------
build_routes_sf <- function(lines_sf, airport_start_sf, plan, crs_m, target_crs, airport_end_sf = NULL) {
  lines_m    <- sf::st_transform(ensure_id(lines_sf), crs_m)
  ap_start_m <- sf::st_transform(airport_start_sf, crs_m)
  ap_end_m   <- if (is.null(airport_end_sf)) ap_start_m else sf::st_transform(airport_end_sf, crs_m)
  
  get_line_xy <- function(id) {
    geom <- lines_m$geometry[which(lines_m$id == id)][[1]]
    sf::st_coordinates(geom)[, c("X","Y"), drop = FALSE]
  }
  maybe_reverse <- function(xy, reverse = FALSE) if (reverse) xy[nrow(xy):1, , drop = FALSE] else xy
  seg_xy <- function(x1,y1,x2,y2) matrix(c(x1,y1,x2,y2), ncol=2, byrow=TRUE)
  
  a0 <- sf::st_coordinates(ap_start_m)[1, ]
  ax <- a0["X"]; ay <- a0["Y"]
  a1 <- sf::st_coordinates(ap_end_m)[1,  ]
  ex <- a1["X"]; ey <- a1["Y"]
  
  build_trip_route_xy <- function(trip_ids, trip_end_endpoint, is_last_trip) {
    route_xy <- matrix(c(ax, ay), ncol=2)
    for (t in seq_along(trip_ids)) {
      id  <- trip_ids[t]
      end <- trip_end_endpoint[t]
      reverse <- (end == 1)       # end at A => traverse B->A
      xy_line <- get_line_xy(id)
      xy_line <- maybe_reverse(xy_line, reverse)
      
      last <- route_xy[nrow(route_xy), ]
      sx <- xy_line[1,1]; sy <- xy_line[1,2]
      if (last[1] != sx || last[2] != sy) {
        route_xy <- rbind(route_xy, seg_xy(last[1], last[2], sx, sy)[2, , drop=FALSE])
      }
      route_xy <- rbind(route_xy, xy_line[-1, , drop=FALSE])
    }
    # return to airport_end (for last trip), otherwise to start (base-return)
    last <- route_xy[nrow(route_xy), ]
    if (isTRUE(is_last_trip)) {
      if (last[1] != ex || last[2] != ey) {
        route_xy <- rbind(route_xy, seg_xy(last[1], last[2], ex, ey)[2, , drop=FALSE])
      }
    } else {
      if (last[1] != ax || last[2] != ay) {
        route_xy <- rbind(route_xy, seg_xy(last[1], last[2], ax, ay)[2, , drop=FALSE])
      }
    }
    route_xy
  }
  
  routes_sfc <- sf::st_sfc(
    lapply(seq_along(plan$trips), function(t) {
      xy <- build_trip_route_xy(plan$trips[[t]], plan$trip_dirs_end_endpoint[[t]], is_last_trip = (t == length(plan$trips)))
      sf::st_linestring(xy)
    }),
    crs = crs_m
  )
  
  routes_sf <- sf::st_sf(trip = seq_along(routes_sfc), geometry = routes_sfc)
  routes_sf <- sf::st_transform(routes_sf, target_crs)
  
  routes_m <- sf::st_transform(routes_sf, crs_m)
  dist_m   <- as.numeric(sf::st_length(routes_m))
  routes_sf$dist_m  <- dist_m
  routes_sf$dist_nm <- dist_m / 1852
  
  routes_sf
}

make_route_arrows <- function(routes_sf, target_crs, crs_m, start_frac = 0.48, end_frac = 0.52) {
  routes_m <- sf::st_transform(routes_sf, crs_m)
  out <- lapply(seq_len(nrow(routes_m)), function(i) {
    pts_m <- sf::st_line_sample(routes_m[i, ], sample = c(start_frac, end_frac))
    pts_m <- sf::st_cast(pts_m, "POINT")
    if (length(pts_m) < 2) return(NULL)
    pts_p <- sf::st_transform(sf::st_as_sf(pts_m), target_crs)
    xy <- sf::st_coordinates(pts_p)
    data.frame(
      trip = routes_sf$trip[i],
      x = xy[1, "X"], y = xy[1, "Y"],
      xend = xy[2, "X"], yend = xy[2, "Y"]
    )
  })
  do.call(rbind, out)
}

# -------- Output helper: detailed flight manifest (legs) --------
build_flight_manifest <- function(dat, plan, speed_knots, dep_label = "DEP", arr_label = "DEP") {
  m_per_nm <- 1852
  manifest <- list()
  
  end_xy   <- function(ii, e) if (e == 1) c(dat$Ax[ii], dat$Ay[ii]) else c(dat$Bx[ii], dat$By[ii])
  start_xy <- function(ii, e_end) if (e_end == 1) c(dat$Bx[ii], dat$By[ii]) else c(dat$Ax[ii], dat$Ay[ii])
  
  for (t in seq_along(plan$trips)) {
    trip_ids <- plan$trips[[t]]
    trip_dir <- plan$trip_dirs_end_endpoint[[t]]
    
    curx <- dat$ax; cury <- dat$ay
    leg_seq <- 1L
    
    for (k in seq_along(trip_ids)) {
      id <- trip_ids[k]
      ii <- match(id, dat$ids)
      end_e   <- trip_dir[k]
      s_xy    <- start_xy(ii, end_e)
      e_xy    <- end_xy(ii, end_e)
      dir_txt <- if (end_e == 1) "B->A" else "A->B"
      
      conn_m <- sqrt((s_xy[1]-curx)^2 + (s_xy[2]-cury)^2)
      if (conn_m > 0) {
        conn_nm  <- conn_m / m_per_nm
        conn_min <- (conn_nm / speed_knots) * 60
        manifest[[length(manifest)+1]] <- data.frame(
          trip = t, leg_seq = leg_seq, leg_type = "connector",
          from = if (leg_seq == 1) dep_label else paste0("T", trip_ids[k-1]),
          to   = paste0("T", id),
          transect_id = NA_integer_, direction = NA_character_,
          dist_m = conn_m, dist_nm = conn_nm, est_minutes = conn_min,
          stringsAsFactors = FALSE
        )
        leg_seq <- leg_seq + 1L
      }
      
      line_m <- as.numeric(dat$line_len[ii])
      line_nm <- line_m / m_per_nm
      line_min <- (line_nm / speed_knots) * 60
      manifest[[length(manifest)+1]] <- data.frame(
        trip = t, leg_seq = leg_seq, leg_type = "transect",
        from = paste0("T", id), to = paste0("T", id),
        transect_id = id, direction = dir_txt,
        dist_m = line_m, dist_nm = line_nm, est_minutes = line_min,
        stringsAsFactors = FALSE
      )
      leg_seq <- leg_seq + 1L
      
      curx <- e_xy[1]; cury <- e_xy[2]
    }
    
    # back connector
    back_m <- sqrt((dat$bx - curx)^2 + (dat$by - cury)^2)
    if (back_m > 0) {
      back_nm  <- back_m / m_per_nm
      back_min <- (back_nm / speed_knots) * 60
      manifest[[length(manifest)+1]] <- data.frame(
        trip = t, leg_seq = leg_seq, leg_type = "connector",
        from = paste0("T", tail(trip_ids,1)), to = if (t == length(plan$trips)) arr_label else dep_label,
        transect_id = NA_integer_, direction = NA_character_,
        dist_m = back_m, dist_nm = back_nm, est_minutes = back_min,
        stringsAsFactors = FALSE
      )
    }
  }
  
  do.call(rbind, manifest)
}

# Clean lines for optimizer: keep valid LINESTRINGs with positive length (>= 1 m)
clean_lines_for_optimizer <- function(lines_sf) {
  stopifnot(inherits(lines_sf, "sf"))
  # In case Z/M crept in, drop them; keep safe types; make valid
  lines_sf <- sf::st_zm(lines_sf, drop = TRUE, what = "ZM")
  lines_sf <- sf::st_make_valid(lines_sf)
  lines_sf <- suppressWarnings(sf::st_collection_extract(lines_sf, "LINESTRING"))
  lines_sf <- lines_sf[!sf::st_is_empty(lines_sf), , drop = FALSE]
  if (nrow(lines_sf) == 0) return(lines_sf)
  # Compute length in a local metric CRS and drop ~0-length pieces
  lines_ll <- sf::st_transform(lines_sf, 4326)
  bb <- sf::st_bbox(lines_ll)
  mid_lon <- as.numeric((bb["xmin"] + bb["xmax"]) / 2)
  mid_lat <- as.numeric((bb["ymin"] + bb["ymax"]) / 2)
  zone    <- floor((mid_lon + 180) / 6) + 1
  crs_m   <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
  len_m   <- as.numeric(sf::st_length(sf::st_transform(lines_sf, crs_m)))
  lines_sf[len_m >= 1, , drop = FALSE]  # keep >= 1 m
}

# ------------------------------ SINGLE-BASE (supports different arrival) ------------------------------
run_transect_optimizer_single <- function(transects_sf,
                                          air_sf,
                                          terre_crs,
                                          start_airport_code = NULL,
                                          end_airport_code   = NULL,
                                          speed_knots = 54,
                                          max_hours = 2,
                                          export_dir = NULL,
                                          out_basename = "routes",
                                          make_plot = TRUE,
                                          terre_crop = NULL,
                                          xlim_terre = NULL,
                                          ylim_terre = NULL,
                                          verbose = TRUE) 
{
  
  
  # ---- NEW: clean the input transects right here ----
  transects_sf <- clean_lines_for_optimizer(transects_sf)
  if (nrow(transects_sf) == 0) stop(
    "No usable transect line coordinates found after cleaning. ",
    "Increase transect length or reduce clipping."
  )
  
  stopifnot(inherits(transects_sf, "sf"), inherits(air_sf, "sf"))
  transects_sf <- ensure_id(transects_sf)
  
  # pick start and end airports
  ap_start <- choose_airport(air_sf, transects_sf, start_airport_code)
  ap_end   <- if (is.null(end_airport_code)) ap_start else choose_airport(air_sf, transects_sf, end_airport_code)
  
  ap_start_ll <- sf::st_transform(ap_start, 4326)
  ap_end_ll   <- sf::st_transform(ap_end,   4326)
  ap_start_xy <- sf::st_coordinates(ap_start_ll)[1, ]
  crs_m <- utm_epsg_from_lonlat(ap_start_xy["X"], ap_start_xy["Y"])
  
  dat <- prep_lines(transects_sf, ap_start_ll, crs_m, airport_end_sf = ap_end_ll)
  n_starts <- max(5, min(20, length(dat$ids)))
  best <- optimize_order_with_refuel(dat, speed_knots, max_hours, n_starts = n_starts, verbose = verbose)
  
  routes_sf <- build_routes_sf(transects_sf, ap_start_ll, best, crs_m, terre_crs, airport_end_sf = ap_end_ll)
  routes_sf$est_hours   <- routes_sf$dist_nm / speed_knots
  routes_sf$est_minutes <- routes_sf$est_hours * 60
  
  # labels
  lab_from <- if ("ident" %in% names(ap_start)) as.character(ap_start$ident[1]) else "DEP"
  lab_to   <- if ("ident" %in% names(ap_end))   as.character(ap_end$ident[1])   else "ARR"
  
  manifest_df <- build_flight_manifest(dat, best, speed_knots, dep_label = lab_from, arr_label = lab_to)
  
  trip_summary <- data.frame(
    trip        = routes_sf$trip,
    dist_nm     = routes_sf$dist_nm,
    est_hours   = routes_sf$est_hours,
    est_minutes = routes_sf$est_minutes,
    start_airport = lab_from,
    end_airport   = lab_to
  )
  
  # ---------- Exports ----------
  if (!is.null(export_dir)) {
    dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
    export_dir <- normalizePath(export_dir, winslash = "/", mustWork = FALSE)
    
    layer_name <- paste0(out_basename, "_flight_routes")
    shp_path   <- file.path(export_dir, paste0(layer_name, ".shp"))
    csv_path   <- file.path(export_dir, paste0(out_basename, "_flight_trips.csv"))
    man_path   <- file.path(export_dir, paste0(out_basename, "_flight_manifest.csv"))
    
    # Shapefile write using dsn/layer (Windows/network friendly)
    exts <- c(".shp",".dbf",".shx",".prj",".cpg")
    for (ext in exts) {
      fp <- file.path(export_dir, paste0(layer_name, ext))
      if (file.exists(fp)) unlink(fp)
    }
    
    write_ok <- TRUE
    tryCatch({
      suppressWarnings(suppressMessages(
        sf::st_write(routes_sf, dsn = export_dir, layer = layer_name,
                     driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE)
      ))
    }, error = function(e) {
      message("Shapefile write failed (", conditionMessage(e), "). Falling back to GeoPackage.")
      write_ok <<- FALSE
    })
    if (!write_ok || !file.exists(shp_path)) {
      gpkg_path <- file.path(export_dir, paste0(out_basename, "_flight_routes.gpkg"))
      suppressWarnings(suppressMessages(
        sf::st_write(routes_sf, dsn = gpkg_path, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
      ))
    }
    
    utils::write.csv(trip_summary, csv_path, row.names = FALSE)
    utils::write.csv(manifest_df, man_path, row.names = FALSE)
  }
  
  # ---------- Optional quick PNG ----------
  if (isTRUE(make_plot)) {
    plot_crs <- terre_crs
    lines_p   <- sf::st_transform(transects_sf, plot_crs)
    routes_p  <- sf::st_transform(routes_sf,  plot_crs)
    ap_start_p<- sf::st_transform(ap_start_ll, plot_crs)
    ap_end_p  <- sf::st_transform(ap_end_ll,   plot_crs)
    
    arrows_mid <- rbind(
      make_route_arrows(routes_sf, plot_crs, crs_m, 0.18, 0.22),
      make_route_arrows(routes_sf, plot_crs, crs_m, 0.48, 0.52),
      make_route_arrows(routes_sf, plot_crs, crs_m, 0.78, 0.82)
    )
    
    geoms <- c(sf::st_geometry(lines_p), sf::st_geometry(routes_p), sf::st_geometry(ap_start_p))
    bb <- sf::st_bbox(geoms)
    pad  <- 0.08
    xpad <- as.numeric(bb["xmax"] - bb["xmin"]) * pad
    ypad <- as.numeric(bb["ymax"] - bb["ymin"]) * pad
    xlim <- if (!is.null(xlim_terre)) xlim_terre else c(bb["xmin"] - xpad, bb["xmax"] + xpad)
    ylim <- if (!is.null(ylim_terre)) ylim_terre else c(bb["ymin"] - ypad, bb["ymax"] + ypad)
    
    gg <- ggplot2::ggplot() +
      { if (!is.null(terre_crop)) ggplot2::geom_sf(data = sf::st_transform(terre_crop, plot_crs),
                                                   fill = "cornsilk", color = "grey60") } +
      ggplot2::geom_sf(data = routes_p, ggplot2::aes(color = factor(trip)), linewidth = 1.1) +
      ggplot2::geom_segment(
        data = arrows_mid,
        ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = factor(trip)),
        arrow = grid::arrow(length = grid::unit(0.15, "inches"), type = "closed"),
        linewidth = 0.9
      ) +
      ggplot2::geom_sf(data = lines_p, color = "red", linewidth = 0.7) +
      ggplot2::geom_sf(data = ap_start_p, color = "black", size = 2.5) +
      ggplot2::geom_sf_text(data = ap_start_p, ggplot2::aes(label = "DEP"), nudge_x = 0.01, size = 3) +
      { if (!all(sf::st_equals(ap_start_p, ap_end_p)[[1]])) ggplot2::geom_sf(data = ap_end_p, color = "black", size = 2.5) } +
      { if (!all(sf::st_equals(ap_start_p, ap_end_p)[[1]])) ggplot2::geom_sf_text(data = ap_end_p, ggplot2::aes(label = "ARR"), nudge_x = 0.01, size = 3) } +
      ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Trip"))
    
    if (!is.null(export_dir)) {
      png_path <- file.path(export_dir, paste0(out_basename, "_flight_routes.png"))
      dir.create(dirname(png_path), showWarnings = FALSE, recursive = TRUE)
      ggplot2::ggsave(filename = png_path, plot = gg, width = 12, height = 10, dpi = 300)
    } else {
      print(gg)  # show in Viewer when not writing files
    }
  }
  
  list(
    plan = best,
    routes_sf = routes_sf,
    trip_summary = trip_summary,
    flight_manifest = manifest_df,
    start_airport = ap_start,
    end_airport = ap_end
  )
}

# ------------------------------ MULTI-BASE (optional) ------------------------------
select_candidate_airports <- function(air_sf, transects_sf, mode = c("auto","list"), k = 2, codes = NULL) {
  mode <- match.arg(mode)
  if (mode == "list") {
    if (is.null(codes) || length(codes) == 0) stop("codes must be provided for mode='list'")
    cols <- intersect(names(air_sf), c("ident","icao","iata","code","name","label","airport_label"))
    mask <- Reduce(`|`, lapply(cols, function(cn) air_sf[[cn]] %in% codes))
    cand <- air_sf[which(mask), , drop = FALSE]
    if (nrow(cand) == 0) stop("None of the provided airport codes/names matched.")
    return(cand)
  } else {
    tr_ll <- sf::st_transform(transects_sf, 4326)
    bbox_ll <- sf::st_bbox(tr_ll)
    mid_lon <- as.numeric((bbox_ll["xmin"] + bbox_ll["xmax"]) / 2)
    mid_lat <- as.numeric((bbox_ll["ymin"] + bbox_ll["ymax"]) / 2)
    crs_m <- if (mid_lat >= 0) paste0("EPSG:", 32600 + floor((mid_lon + 180)/6) + 1) else paste0("EPSG:", 32700 + floor((mid_lon + 180)/6) + 1)
    tr_m  <- sf::st_transform(transects_sf, crs_m)
    air_m <- sf::st_transform(air_sf, crs_m)
    center_m <- sf::st_centroid(sf::st_union(sf::st_geometry(tr_m)))
    d <- sf::st_distance(sf::st_geometry(air_m), center_m)
    ord <- order(as.numeric(d))
    air_m[ord[seq_len(min(k, length(ord)))], , drop = FALSE] |> sf::st_transform(sf::st_crs(air_sf))
  }
}

assign_lines_to_airports <- function(lines_sf, cand_air_sf) {
  # Determine a local projected CRS (UTM) from lines bbox center
  lines_ll <- sf::st_transform(lines_sf, 4326)
  bbox_ll  <- sf::st_bbox(lines_ll)
  mid_lon  <- as.numeric((bbox_ll["xmin"] + bbox_ll["xmax"]) / 2)
  mid_lat  <- as.numeric((bbox_ll["ymin"] + bbox_ll["ymax"]) / 2)
  zone     <- floor((mid_lon + 180) / 6) + 1
  crs_m    <- if (mid_lat >= 0) paste0("EPSG:", 32600 + zone) else paste0("EPSG:", 32700 + zone)
  
  # Transform to projected CRS **before** sampling
  lines_m <- sf::st_transform(lines_sf, crs_m)
  air_m   <- sf::st_transform(cand_air_sf, crs_m)
  
  # Sample line midpoints in meters (projected)
  lines_mid_m <- sf::st_line_sample(lines_m, sample = 0.5)
  lines_mid_m <- sf::st_cast(lines_mid_m, "POINT")
  # Convert to sf POINT with same crs_m
  lm <- sf::st_as_sf(lines_mid_m)
  
  # Distance matrix: airports x line midpoints
  D <- sf::st_distance(sf::st_geometry(air_m), sf::st_geometry(lm))
  nearest_idx <- apply(as.matrix(D), 2, which.min)  # nearest airport for each line
  
  # Split original lines (in original CRS) by assigned airport index
  splits <- split(seq_len(nrow(lines_sf)), nearest_idx)
  
  lapply(names(splits), function(i) {
    idx <- splits[[i]]
    list(
      airport = cand_air_sf[as.integer(i), , drop = FALSE],
      lines   = lines_sf[idx, , drop = FALSE]
    )
  })
}

run_transect_optimizer_multi <- function(transects_sf,
                                         air_sf,
                                         terre_crs,
                                         candidate_mode = c("auto","list"),
                                         candidate_k = 2,
                                         candidate_codes = NULL,
                                         speed_knots = 54,
                                         max_hours = 2,
                                         export_dir = NULL,
                                         out_basename = "routes",
                                         make_plot = TRUE,
                                         terre_crop = NULL,
                                         xlim_terre = NULL,
                                         ylim_terre = NULL,
                                         verbose = TRUE,
                                         # NEW: anchors passed from wrapper
                                         start_airport_code = NULL,
                                         end_airport_code   = NULL) {
  
  # ---- NEW: clean the input transects right here ----
  transects_sf <- clean_lines_for_optimizer(transects_sf)
  if (nrow(transects_sf) == 0) stop(
    "No usable transect line coordinates found after cleaning. ",
    "Increase transect length or reduce clipping."
  )
  
  candidate_mode <- match.arg(candidate_mode)
  # Candidate set
  cand_air <- select_candidate_airports(air_sf, transects_sf,
                                        mode = candidate_mode, k = candidate_k, codes = candidate_codes)
  
  # --- NEW: make sure typed start/end (if any) are INCLUDED in candidate set (even if far) ---
  force_codes <- unique(na.omit(c(start_airport_code, end_airport_code)))
  if (length(force_codes) > 0) {
    force_idx <- integer(0)
    for (q in force_codes) {
      if (!nzchar(q)) next
      key <- toupper(gsub("\\s+","", trimws(as.character(q))))
      idx <- integer(0)
      # search in full air_sf
      if ("ident" %in% names(air_sf)) {
        v <- toupper(gsub("\\s+","", trimws(as.character(air_sf$ident))))
        idx <- c(idx, which(v == key))
      }
      for (cn in intersect(c("icao_code","gps_code","icao","iata","iata_code","local_code","local","code"), names(air_sf))) {
        v <- toupper(gsub("\\s+","", trimws(as.character(air_sf[[cn]]))))
        idx <- c(idx, which(v == key))
      }
      idx <- unique(idx)
      if (length(idx) > 0) force_idx <- c(force_idx, idx[1])
    }
    if (length(force_idx) > 0) {
      cand_air <- suppressWarnings(dplyr::bind_rows(
        cand_air,
        air_sf[unique(force_idx), , drop = FALSE]
      )) |>
        dplyr::distinct(ident, .keep_all = TRUE)  # de-dup by ident (or pick a stable key)
    }
  }
  
  if (nrow(cand_air) == 0) stop("No candidate airports available for multi-base.")
  
  
  
  # Metric CRS for planning
  crs_m <- crs_m_from_bbox(transects_sf)
  
  # A global plan (order/partition) does not depend on airport choices
  # Use nearest (or forced start) just to build the DP 'dat0'
  ap_start <- if (!is.null(start_airport_code) && nzchar(start_airport_code))
    choose_airport(cand_air, transects_sf, start_airport_code) else
      nearest_airport(transects_sf, cand_air)
  ap_end   <- if (!is.null(end_airport_code) && nzchar(end_airport_code))
    choose_airport(cand_air, transects_sf, end_airport_code) else ap_start
  
  ap_start_ll <- sf::st_transform(ap_start, 4326)
  ap_end_ll   <- sf::st_transform(ap_end,   4326)
  
  dat0 <- prep_lines(transects_sf, ap_start_ll, crs_m, airport_end_sf = ap_end_ll)
  n_starts <- max(5, min(20, length(dat0$ids)))
  if (isTRUE(verbose)) message("Computing global plan...")
  best <- optimize_order_with_refuel(dat0, speed_knots, max_hours, n_starts = n_starts, verbose = verbose)
  
  # Precompute candidate airport coords in meters
  cand_ll <- sf::st_transform(cand_air, 4326)
  cand_m  <- sf::st_transform(cand_ll,  crs_m)
  ap_coords <- sf::st_coordinates(cand_m)
  
  # Utility: choose the best (start,end) for each trip, accounting for ferry from prev end
  pick_trip_airports <- function(prev_end_idx, first_xy, last_xy, is_first, is_last) {
    forced_start_idx <- if (!is.null(start_airport_code) && nzchar(start_airport_code) && is_first) {
      cols <- intersect(names(cand_ll), c("ident","icao","iata","code","name","airport_label","label"))
      which(Reduce(`|`, lapply(cols, function(cn) tolower(as.character(cand_ll[[cn]])) == tolower(start_airport_code))))[1]
    } else NA_integer_
    
    forced_end_idx <- if (!is.null(end_airport_code) && nzchar(end_airport_code) && is_last) {
      cols <- intersect(names(cand_ll), c("ident","icao","iata","code","name","airport_label","label"))
      which(Reduce(`|`, lapply(cols, function(cn) tolower(as.character(cand_ll[[cn]])) == tolower(end_airport_code))))[1]
    } else NA_integer_
    
    best_cost <- Inf; best_pair <- c(NA_integer_, NA_integer_)
    for (si in seq_len(nrow(cand_ll))) {
      if (!is.na(forced_start_idx) && si != forced_start_idx) next
      
      # ferry from prev end -> this start (0 if first trip or prev_end_idx is NA)
      ferry_m <- 0
      if (!is.na(prev_end_idx)) {
        prev_xy <- ap_coords[prev_end_idx, ]
        this_xy <- ap_coords[si, ]
        ferry_m <- sqrt((this_xy["X"] - prev_xy["X"])^2 + (this_xy["Y"] - prev_xy["Y"])^2)
      }
      # start AP -> first line start
      sxy <- ap_coords[si, ]
      start_conn_m <- sqrt((first_xy["X"] - sxy["X"])^2 + (first_xy["Y"] - sxy["Y"])^2)
      
      for (ei in seq_len(nrow(cand_ll))) {
        if (!is.na(forced_end_idx) && ei != forced_end_idx) next
        exy <- ap_coords[ei, ]
        end_conn_m <- sqrt((last_xy["X"] - exy["X"])^2 + (last_xy["Y"] - exy["Y"])^2)
        total_m <- ferry_m + start_conn_m + end_conn_m
        if (total_m < best_cost) { best_cost <- total_m; best_pair <- c(si, ei) }
      }
    }
    best_pair
  }
  
  # Get first/last line XY for each trip (given global 'best' directions)
  lines_m <- sf::st_transform(ensure_id(transects_sf), crs_m)
  get_line_xy <- function(id) {
    geom <- lines_m$geometry[which(lines_m$id == id)][[1]]
    sf::st_coordinates(geom)[, c("X","Y"), drop = FALSE]
  }
  trip_endpoints_xy <- function(t) {
    ids <- best$trips[[t]]; dirs <- best$trip_dirs_end_endpoint[[t]]
    # first
    xy1 <- get_line_xy(ids[1]); if (dirs[1] == 1) xy1 <- xy1[nrow(xy1):1, , drop=FALSE]
    sxy <- xy1[1, ]
    # last
    xyn <- get_line_xy(tail(ids, 1)); if (tail(dirs,1) == 1) xyn <- xyn[nrow(xyn):1, , drop=FALSE]
    exy <- xyn[nrow(xyn), ]
    list(start_xy = sxy, end_xy = exy)
  }
  
  # Choose per-trip airports in sequence (this allows trips to hop between airports)
  Tn <- length(best$trips)
  trip_airports <- list(start = vector("list", Tn), end = vector("list", Tn))
  dep_labels <- character(Tn); arr_labels <- character(Tn)
  prev_end_idx <- NA_integer_
  
  for (t in seq_len(Tn)) {
    ep <- trip_endpoints_xy(t)
    pair <- pick_trip_airports(prev_end_idx, ep$start_xy, ep$end_xy, is_first = (t==1), is_last = (t==Tn))
    si <- pair[1]; ei <- pair[2]
    trip_airports$start[[t]] <- cand_ll[si, , drop = FALSE]
    trip_airports$end[[t]]   <- cand_ll[ei, , drop = FALSE]
    dep_labels[t] <- airport_label(cand_ll[si, , drop = FALSE])
    arr_labels[t] <- airport_label(cand_ll[ei, , drop = FALSE])
    prev_end_idx <- ei
  }
  
  # Build outputs
  routes_sf <- build_routes_sf_per_trip_airports(transects_sf, best, crs_m, terre_crs, trip_airports)
  # overwrite time cols with the actual speed
  routes_sf$est_hours   <- routes_sf$dist_nm / speed_knots
  routes_sf$est_minutes <- routes_sf$est_hours * 60
  
  manifest_df <- build_manifest_with_airports(transects_sf, best, crs_m, speed_knots, trip_airports,
                                              dep_labels = dep_labels, arr_labels = arr_labels)
  
  # Export combined set
  if (!is.null(export_dir)) {
    dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
    export_dir <- normalizePath(export_dir, winslash = "/", mustWork = FALSE)
    
    layer_name <- paste0(out_basename, "_flight_routes")
    shp_path   <- file.path(export_dir, paste0(layer_name, ".shp"))
    trips_path <- file.path(export_dir, paste0(out_basename, "_flight_trips.csv"))
    man_path   <- file.path(export_dir, paste0(out_basename, "_flight_manifest.csv"))
    
    # Shapefile write using dsn/layer (Windows/network friendly)
    exts <- c(".shp",".dbf",".shx",".prj",".cpg")
    for (ext in exts) {
      fp <- file.path(export_dir, paste0(layer_name, ext))
      if (file.exists(fp)) unlink(fp)
    }
    write_ok <- TRUE
    tryCatch({
      suppressWarnings(suppressMessages(
        sf::st_write(routes_sf, dsn = export_dir, layer = layer_name,
                     driver = "ESRI Shapefile", delete_layer = TRUE, quiet = TRUE)
      ))
    }, error = function(e) {
      message("Combined Shapefile write failed (", conditionMessage(e), "). Falling back to GeoPackage.")
      write_ok <<- FALSE
    })
    if (!write_ok || !file.exists(shp_path)) {
      gpkg_path <- file.path(export_dir, paste0(out_basename, "_flight_routes.gpkg"))
      suppressWarnings(suppressMessages(
        sf::st_write(routes_sf, dsn = gpkg_path, driver = "GPKG", delete_dsn = TRUE, quiet = TRUE)
      ))
    }
    
    # Per-trip totals table (trip 1..N in execution order)
    trips_all <- data.frame(
      trip = routes_sf$trip,
      dist_nm = routes_sf$dist_nm,
      est_hours = routes_sf$est_hours,
      est_minutes = routes_sf$est_minutes,
      start_airport = dep_labels[routes_sf$trip],
      end_airport   = arr_labels[routes_sf$trip],
      stringsAsFactors = FALSE
    )
    utils::write.csv(trips_all, trips_path, row.names = FALSE)
    utils::write.csv(manifest_df, man_path, row.names = FALSE)
    
    if (isTRUE(make_plot)) {
      plot_crs <- terre_crs
      lines_p  <- sf::st_transform(transects_sf, plot_crs)
      routes_p <- sf::st_transform(routes_sf,  plot_crs)
      geoms <- c(sf::st_geometry(lines_p), sf::st_geometry(routes_p))
      bb <- sf::st_bbox(geoms); pad <- 0.08
      xpad <- as.numeric(bb["xmax"] - bb["xmin"]) * pad
      ypad <- as.numeric(bb["ymax"] - bb["ymin"]) * pad
      xlim <- if (!is.null(xlim_terre)) xlim_terre else c(bb["xmin"] - xpad, bb["xmax"] + xpad)
      ylim <- if (!is.null(ylim_terre)) ylim_terre else c(bb["ymin"] - ypad, bb["ymax"] + ypad)
      gg <- ggplot2::ggplot() +
        { if (!is.null(terre_crop)) ggplot2::geom_sf(data = sf::st_transform(terre_crop, plot_crs),
                                                     fill = "cornsilk", color = "grey60") } +
        ggplot2::geom_sf(data = routes_p, ggplot2::aes(color = factor(trip)), linewidth = 1.1) +
        ggplot2::geom_sf(data = lines_p, color = "red", linewidth = 0.7) +
        ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
        ggplot2::theme_minimal() +
        ggplot2::guides(color = ggplot2::guide_legend(title = "Trip"))
      png_path <- file.path(export_dir, paste0(out_basename, "_flight_routes.png"))
      ggplot2::ggsave(filename = png_path, plot = gg, width = 12, height = 10, dpi = 300)
    }
  }
  
  list(
    routes_sf = routes_sf,
    trip_summary = data.frame(
      trip = routes_sf$trip,
      dist_nm = routes_sf$dist_nm,
      est_hours = routes_sf$est_hours,
      est_minutes = routes_sf$est_minutes,
      start_airport = dep_labels[routes_sf$trip],
      end_airport   = arr_labels[routes_sf$trip],
      stringsAsFactors = FALSE
    ),
    flight_manifest = manifest_df
  )
}

# ------------------------------ Backward-compatible wrapper ------------------------------
run_transect_optimizer <- function(transects_sf,
                                   air_sf,
                                   terre_crs,
                                   airport_code = NULL,         # start
                                   speed_knots = 54,
                                   max_hours = 2,
                                   export_dir = NULL,
                                   out_basename = "routes",
                                   make_plot = TRUE,
                                   terre_crop = NULL,
                                   xlim_terre = NULL,
                                   ylim_terre = NULL,
                                   end_airport_code = NULL,
                                   multi_airports = NULL,       # NULL | "auto" | character vector
                                   multi_k = 2,
                                   verbose = TRUE) {
  
  if (is.null(multi_airports)) {
    return(
      run_transect_optimizer_single(
        transects_sf = transects_sf,
        air_sf       = air_sf,
        terre_crs    = terre_crs,
        start_airport_code = airport_code,
        end_airport_code   = end_airport_code,
        speed_knots = speed_knots,
        max_hours   = max_hours,
        export_dir  = export_dir,
        out_basename = out_basename,
        make_plot   = make_plot,
        terre_crop  = terre_crop,
        xlim_terre  = xlim_terre,
        ylim_terre  = ylim_terre,
        verbose     = verbose
      )
    )
  }
  
  mode  <- if (identical(multi_airports, "auto")) "auto" else "list"
  codes <- if (mode == "list") multi_airports else NULL
  run_transect_optimizer_multi(
    transects_sf = transects_sf,
    air_sf       = air_sf,
    terre_crs    = terre_crs,
    candidate_mode = mode,
    candidate_k    = multi_k,
    candidate_codes = codes,
    speed_knots = speed_knots,
    max_hours   = max_hours,
    export_dir  = export_dir,
    out_basename = out_basename,
    make_plot   = make_plot,
    terre_crop  = terre_crop,
    xlim_terre  = xlim_terre,
    ylim_terre  = ylim_terre,
    verbose     = verbose,
    start_airport_code = airport_code,
    end_airport_code   = end_airport_code
  )
}

