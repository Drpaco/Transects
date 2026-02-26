# ============================================================
# GSD + Ground Coverage with Excel export & persistent user catalog
# - Camera model lookup from built-ins + camera_catalog_user.csv (if present)
# - Start/End/By for focal length (mm) & altitude (ft)
# - Output: GSD (cm/pixel), Coverage ("Width x Length" in ft or m) + Area (ft^2 or m^2)
# - Excel export (.xlsx): one workbook, multiple sheets
# - Heatmaps: interactive loop; try any metric until "OK/Done"
# - NEW: After choosing an output type that includes coverage, ask for units (ft/m)
# To disable autorun: options(gsd.autorun = FALSE) before source().
# ============================================================

# ---------- Robust interactive input helper ----------
get_input_safe <- function(prompt, default = "") {
  if (!interactive()) return(default)
  ans <- tryCatch(readline(prompt), error = function(e) default)
  ans <- trimws(ans)
  if (nchar(ans) == 0) default else ans
}

# ---------- Small utilities ----------
read_num_safe <- function(prompt, default, min_val = -Inf, max_val = Inf) {
  repeat {
    s <- get_input_safe(prompt, as.character(default))
    v <- suppressWarnings(as.numeric(s))
    if (!is.na(v) && v >= min_val && v <= max_val) return(v)
    if (!interactive()) {
      warning(sprintf("Invalid numeric input. Using default %s.", default))
      return(default)
    }
    cat(sprintf("  ! Invalid input. Enter a number between %g and %g.\n", min_val, max_val)); flush.console()
  }
}

read_seq_params_safe <- function(name, unit, defaults) {
  cat(sprintf("\n-- %s (%s) --\n", name, unit)); flush.console()
  start_val <- read_num_safe(sprintf("  Start %s (%s) [%s]: ", name, unit, defaults$start),
                             default = defaults$start)
  end_val   <- read_num_safe(sprintf("  End %s (%s) [%s]: ",   name, unit, defaults$end),
                             default = defaults$end)
  by_val    <- read_num_safe(sprintf("  Step %s (%s) [%s]: ",  name, unit, defaults$by),
                             default = defaults$by, min_val = .Machine$double.eps)
  
  if (end_val < start_val) {
    cat("  ! End < Start. Swapping them.\n"); flush.console()
    tmp <- start_val; start_val <- end_val; end_val <- tmp
  }
  if (by_val <= 0) {
    cat("  ! Step must be > 0. Using default.\n"); flush.console()
    by_val <- defaults$by
  }
  list(start = start_val, end = end_val, by = by_val)
}

make_seq_inclusive <- function(start, end, by) {
  if (by <= 0) stop("`by` must be > 0")
  seq_vals <- seq(from = start, to = end, by = by)
  if (length(seq_vals) == 0) {
    seq_vals <- c(start, end)
  } else {
    last <- tail(seq_vals, 1)
    if (end > last && (end - last) > (abs(by) * 1e-7)) {
      seq_vals <- c(seq_vals, end)
    }
  }
  sort(unique(seq_vals))
}

# ---------- Camera catalog (built-in) ----------
# Built-ins: DJI Air 2S, Mavic 3 Classic, Sony a7R IV, Canon 5D Mark IV, Nikon D850, Sony RX100 VII.
camera_catalog_builtin <- data.frame(
  model = c(
    "DJI Air 2S",
    "DJI Mavic 3 Classic",
    "Sony a7R IV (ILCE-7RM4A)",
    "Canon EOS 5D Mark IV",
    "Nikon D850",
    "Sony RX100 VII"
  ),
  mp = c(20.0, 20.0, 61.0, 30.4, 45.7, 20.1),
  sensor_w_mm = c(13.2, 17.3, 35.7, 36.0, 35.9, 13.2),
  sensor_h_mm = c( 8.8, 13.0, 23.8, 24.0, 23.9,  8.8),
  aka = c(
    "Air2S|Air 2 S|Air-2S",
    "Mavic 3|Mavic3 Classic|M3 Classic",
    "A7R IV|α7R IV|ILCE-7RM4|ILCE-7RM4A",
    "5D4|5D Mk IV|5DM4",
    "D-850|Nikon D 850",
    "RX100M7|RX100 VII|DSC-RX100M7"
  ),
  stringsAsFactors = FALSE
)

# ---------- Persistent user catalog: CSV merge utilities ----------
.user_catalog_path <- function() "camera_catalog_user.csv"

load_camera_catalog <- function() {
  cat_path <- .user_catalog_path()
  if (file.exists(cat_path)) {
    user <- tryCatch({
      read.csv(cat_path, stringsAsFactors = FALSE, check.names = FALSE)
    }, error = function(e) {
      warning(sprintf("Could not read '%s': %s", cat_path, e$message))
      NULL
    })
  } else {
    user <- NULL
  }
  want <- c("model","mp","sensor_w_mm","sensor_h_mm","aka")
  if (!is.null(user)) {
    missing <- setdiff(want, names(user))
    for (m in missing) user[[m]] <- if (m == "aka") "" else NA_real_
    user <- user[ , want]
  }
  if (!is.null(user) && nrow(user)) {
    user <- user[!duplicated(tolower(user$model), fromLast=TRUE), ]
    bi <- camera_catalog_builtin[!tolower(camera_catalog_builtin$model) %in% tolower(user$model), ]
    out <- rbind(bi, user)
  } else {
    out <- camera_catalog_builtin
  }
  rownames(out) <- NULL
  out
}

save_camera_catalog_user <- function(df_user) {
  write.csv(df_user, .user_catalog_path(), row.names = FALSE)
  invisible(TRUE)
}

list_known_cameras <- function(camera_catalog) {
  cat("\nKnown camera presets:\n")
  for (i in seq_len(nrow(camera_catalog))) {
    cat(sprintf("  - %s  [%.1f MP, %.2f x %.2f mm]\n",
                camera_catalog$model[i], camera_catalog$mp[i],
                camera_catalog$sensor_w_mm[i], camera_catalog$sensor_h_mm[i]))
  }
}

lookup_camera <- function(query, camera_catalog) {
  if (is.null(query) || trimws(query) == "") return(NULL)
  q <- trimws(query)
  if (tolower(q) %in% c("list","help","?")) { list_known_cameras(camera_catalog); return(NULL) }
  hits <- which(grepl(q, camera_catalog$model, ignore.case = TRUE) |
                  grepl(q, camera_catalog$aka,   ignore.case = TRUE))
  if (length(hits) == 0) return(NULL)
  if (length(hits) == 1) return(camera_catalog[hits, , drop = FALSE])
  if (!interactive()) return(camera_catalog[hits[1], , drop = FALSE])
  cat("\nMultiple matches:\n")
  for (i in seq_along(hits)) {
    r <- camera_catalog[hits[i], ]
    cat(sprintf("  [%d] %s  (%.1f MP, %.2f x %.2f mm)\n",
                i, r$model, r$mp, r$sensor_w_mm, r$sensor_h_mm))
  }
  repeat {
    s <- get_input_safe(sprintf("Choose 1-%d or 0 to cancel [1]: ", length(hits)), default = "1")
    k <- suppressWarnings(as.integer(s))
    if (!is.na(k) && k == 0) return(NULL)
    if (!is.na(k) && k >= 1 && k <= length(hits)) return(camera_catalog[hits[k], , drop = FALSE])
    cat("  ! Invalid selection.\n")
  }
}

# Optional: interactive add/update; persists to CSV layer
add_or_update_camera_interactive <- function() {
  path <- .user_catalog_path()
  user <- if (file.exists(path)) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  else data.frame(model=character(), mp=numeric(), sensor_w_mm=numeric(),
                  sensor_h_mm=numeric(), aka=character(), stringsAsFactors = FALSE)
  
  cat("\n=== Add/Update camera preset ===\n")
  model <- get_input_safe("Model name (unique key): ", default = "")
  if (model == "") { cat("  ! Cancelled.\n"); return(invisible(FALSE)) }
  
  already <- which(tolower(user$model) == tolower(model))
  mp        <- suppressWarnings(as.numeric(get_input_safe("Effective megapixels (e.g., 20.0): ", "20")))
  sw        <- suppressWarnings(as.numeric(get_input_safe("Sensor width (mm): ", "13.2")))
  sh        <- suppressWarnings(as.numeric(get_input_safe("Sensor height (mm): ", "8.8")))
  aka       <- get_input_safe("Aliases (optional, '|' separated): ", "")
  
  if (any(is.na(c(mp,sw,sh))) || mp <= 0 || sw <= 0 || sh <= 0) {
    cat("  ! Invalid numeric values. Aborting.\n"); return(invisible(FALSE))
  }
  
  row <- data.frame(model=model, mp=mp, sensor_w_mm=sw, sensor_h_mm=sh, aka=aka, stringsAsFactors = FALSE)
  
  if (length(already)) {
    user[already[1], ] <- row
    cat("  -> Updated existing preset.\n")
  } else {
    user <- rbind(user, row)
    cat("  -> Added new preset.\n")
  }
  save_camera_catalog_user(user)
  cat(sprintf("  -> Saved to %s\n", path))
  invisible(TRUE)
}

# ---------- Excel writer (openxlsx preferred; fallback writexl) ----------
.save_excel <- function(path, sheets) {
  have_openxlsx <- requireNamespace("openxlsx", quietly = TRUE)
  have_writexl  <- requireNamespace("writexl",  quietly = TRUE)
  if (!have_openxlsx && !have_writexl) {
    stop("Please install an Excel writer: install.packages('openxlsx') or install.packages('writexl')")
  }
  if (have_openxlsx) {
    wb <- openxlsx::createWorkbook()
    for (nm in names(sheets)) {
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, nm, sheets[[nm]])
    }
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  } else {
    writexl::write_xlsx(sheets, path = path)
  }
  invisible(TRUE)
}

# ---------- Heatmap helper (robust zlim handling) ----------
.plot_heatmap <- function(Z_alt_by_focal, main_title,
                          xlab = "Focal length (mm)", ylab = "Altitude (ft)",
                          focal_labels = NULL, alt_labels = NULL, zlim = NULL) {
  Z <- t(as.matrix(Z_alt_by_focal))  # [focal x alt]
  X <- seq_len(nrow(Z)); Y <- seq_len(ncol(Z))
  if (is.null(focal_labels)) focal_labels <- colnames(Z_alt_by_focal)
  if (is.null(alt_labels))   alt_labels   <- rownames(Z_alt_by_focal)
  focal_labels <- sub(" mm$", "", focal_labels)
  alt_labels   <- sub(" ft.*$", "", alt_labels)
  if (any(!is.finite(Z))) Z[!is.finite(Z)] <- NA_real_
  
  if (is.null(zlim) || length(zlim) != 2 || any(!is.finite(zlim)) || diff(zlim) <= 0) {
    rng <- range(Z, na.rm = TRUE, finite = TRUE)
    if (all(is.finite(rng)) && diff(rng) > 0) {
      zlim <- rng
    } else {
      val <- suppressWarnings(as.numeric(stats::median(Z, na.rm = TRUE)))
      zlim <- if (is.finite(val)) c(val - 1e-9, val + 1e-9) else c(0, 1)
    }
  }
  
  image(x = X, y = Y, z = Z,
        xlab = xlab, ylab = ylab,
        col = hcl.colors(32, "YlOrRd", rev = TRUE),
        axes = FALSE, useRaster = TRUE,
        zlim = zlim)
  axis(1, at = X, labels = focal_labels, las = 2, cex.axis = 0.85)
  axis(2, at = Y, labels = alt_labels,   las = 1, cex.axis = 0.85)
  box()
  contour(x = X, y = Y, z = Z, add = TRUE, drawlabels = TRUE, col = "gray30")
  title(main = main_title)
}

# ---------- Looping heatmap menu ----------
.show_heatmap_loop <- function(show_gsd, show_cov, gsd_cm, width_out, length_out, area_out, cov_units_label) {
  if (nrow(gsd_cm) < 2 || ncol(gsd_cm) < 2) {
    warning(sprintf("Heatmap needs at least 2 altitude classes and 2 focal classes. You have %d x %d.",
                    nrow(gsd_cm), ncol(gsd_cm)))
    return(invisible())
  }
  
  repeat {
    if (show_gsd && !show_cov) {
      cat("\nHeatmap metric:\n  1 = GSD (cm/pixel)\n  0 = OK / Done\n")
      s <- tolower(get_input_safe("Choose [1]: ", default = "1"))
      if (s %in% c("0","ok","done","q","quit","exit")) break
      old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
      par(mar = c(6, 6, 3, 1))
      .plot_heatmap(gsd_cm, "Ground Sample Distance (cm/pixel)")
      cat("  Heatmap drawn.\n")
      
    } else if (!show_gsd && show_cov) {
      cat("\nHeatmap metric (coverage):\n")
      cat(sprintf("  1 = Width (%s)\n", cov_units_label))
      cat(sprintf("  2 = Length (%s)\n", cov_units_label))
      cat(sprintf("  3 = Area (%s^2)\n", cov_units_label))
      cat("  4 = ALL Coverage (Width, Length, Area)\n")
      cat("  5 = Width & Length together (side-by-side)\n")
      cat("  0 = OK / Done\n")
      s <- tolower(get_input_safe("Choose [4]: ", default = "4"))
      if (s %in% c("0","ok","done","q","quit","exit")) break
      choice <- suppressWarnings(as.integer(s))
      
      old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
      if (choice == 1) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(width_out,  sprintf("Ground Coverage Width (%s)", cov_units_label))
      } else if (choice == 2) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(length_out, sprintf("Ground Coverage Length (%s)", cov_units_label))
      } else if (choice == 3) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(area_out,   sprintf("Ground Coverage Area (%s^2)", cov_units_label))
      } else if (choice == 5) {
        zlim <- range(c(width_out, length_out), na.rm = TRUE, finite = TRUE)
        if (!all(is.finite(zlim)) || diff(zlim) <= 0) {
          zlim <- range(width_out, na.rm = TRUE, finite = TRUE)
          if (!all(is.finite(zlim)) || diff(zlim) <= 0) zlim <- c(0, 1)
        }
        par(mfrow = c(1, 2), mar = c(6, 6, 3, 1))
        .plot_heatmap(width_out,  sprintf("Width (%s)",  cov_units_label), zlim = zlim)
        .plot_heatmap(length_out, sprintf("Length (%s)", cov_units_label), zlim = zlim)
        par(mfrow = c(1, 1))
      } else {
        par(mfrow = c(1, 3), mar = c(6, 6, 3, 1))
        .plot_heatmap(width_out,  sprintf("Width (%s)",  cov_units_label))
        .plot_heatmap(length_out, sprintf("Length (%s)", cov_units_label))
        .plot_heatmap(area_out,   sprintf("Area (%s^2)", cov_units_label))
        par(mfrow = c(1, 1))
      }
      cat("  Heatmap(s) drawn.\n")
      
    } else { # show_gsd && show_cov
      cat("\nHeatmap metric:\n")
      cat("  1 = GSD (cm/pixel)\n")
      cat(sprintf("  2 = Coverage: Width (%s)\n",  cov_units_label))
      cat(sprintf("  3 = Coverage: Length (%s)\n", cov_units_label))
      cat(sprintf("  4 = Coverage: Area (%s^2)\n", cov_units_label))
      cat("  5 = ALL (GSD, Width, Length, Area)\n")
      cat("  6 = Coverage: Width & Length together (side-by-side)\n")
      cat("  0 = OK / Done\n")
      s <- tolower(get_input_safe("Choose [5]: ", default = "5"))
      if (s %in% c("0","ok","done","q","quit","exit")) break
      choice <- suppressWarnings(as.integer(s))
      
      old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
      if (choice == 1) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(gsd_cm,   "Ground Sample Distance (cm/pixel)")
      } else if (choice == 2) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(width_out,  sprintf("Ground Coverage Width (%s)", cov_units_label))
      } else if (choice == 3) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(length_out, sprintf("Ground Coverage Length (%s)", cov_units_label))
      } else if (choice == 4) {
        par(mar = c(6, 6, 3, 1)); .plot_heatmap(area_out,   sprintf("Ground Coverage Area (%s^2)", cov_units_label))
      } else if (choice == 6) {
        zlim <- range(c(width_out, length_out), na.rm = TRUE, finite = TRUE)
        if (!all(is.finite(zlim)) || diff(zlim) <= 0) {
          zlim <- range(width_out, na.rm = TRUE, finite = TRUE)
          if (!all(is.finite(zlim)) || diff(zlim) <= 0) zlim <- c(0, 1)
        }
        par(mfrow = c(1, 2), mar = c(6, 6, 3, 1))
        .plot_heatmap(width_out,  sprintf("Width (%s)",  cov_units_label), zlim = zlim)
        .plot_heatmap(length_out, sprintf("Length (%s)", cov_units_label), zlim = zlim)
        par(mfrow = c(1, 1))
      } else {
        par(mfrow = c(2, 2), mar = c(6, 6, 3, 1))
        .plot_heatmap(gsd_cm,     "GSD (cm/pixel)")
        .plot_heatmap(width_out,  sprintf("Width (%s)",  cov_units_label))
        .plot_heatmap(length_out, sprintf("Length (%s)", cov_units_label))
        .plot_heatmap(area_out,   sprintf("Area (%s^2)", cov_units_label))
        par(mfrow = c(1, 1))
      }
      cat("  Heatmap(s) drawn.\n")
    }
  }
  
  cat("  OK — exiting heatmap menu.\n"); flush.console()
  invisible(TRUE)
}

# ---------- Main function ----------
run_gsd_ranges <- function(
    defaults = list(
      mp = 20,
      sensor_width_mm  = 13.2,
      sensor_height_mm = 8.8,
      focal_start_mm = 24, focal_end_mm = 85, focal_by_mm = 5,
      alt_start_ft   = 100, alt_end_ft   = 400, alt_by_ft   = 50
    ),
    ask_save = TRUE,
    ask_plot = TRUE
) {
  cat("=== Camera & Sensor Parameters ===\n"); flush.console()
  
  # Load catalog (built-ins + user CSV override)
  camera_catalog <- load_camera_catalog()
  
  # Optional quick edit before starting
  edit_presets <- tolower(get_input_safe("\nEdit camera presets first? (y/n) [n]: ", "n"))
  if (edit_presets %in% c("y","yes")) {
    add_or_update_camera_interactive()
    camera_catalog <- load_camera_catalog()  # reload after edit
  }
  
  # Camera lookup OR manual
  cat("\nType a camera model to auto-fill (e.g., 'DJI Air 2S', 'Sony a7R IV'),\n")
  cat("or type 'list' to see presets, or press Enter to input manually.\n")
  cam_query <- get_input_safe("Camera model [manual]: ", default = "")
  cam_hit <- lookup_camera(cam_query, camera_catalog)
  
  if (!is.null(cam_hit)) {
    cat(sprintf("\nFound: %s\n  -> %.1f MP, sensor %.2f x %.2f mm\n",
                cam_hit$model, cam_hit$mp, cam_hit$sensor_w_mm, cam_hit$sensor_h_mm)); flush.console()
    use_it <- tolower(get_input_safe("Use these values? (y/n) [y]: ", default = "y"))
    if (use_it %in% c("y","yes")) {
      mp <- cam_hit$mp
      sensor_width_mm  <- cam_hit$sensor_w_mm
      sensor_height_mm <- cam_hit$sensor_h_mm
    } else {
      mp <- read_num_safe("Total pixels (megapixels): ", default = defaults$mp, min_val = 0.1)
      sensor_width_mm  <- read_num_safe("Sensor width (mm): ", default = defaults$sensor_width_mm,  min_val = 1e-6)
      sensor_height_mm <- read_num_safe("Sensor height (mm): ", default = defaults$sensor_height_mm, min_val = 1e-6)
    }
  } else {
    mp <- read_num_safe("Total pixels (megapixels, e.g., 20): ", default = defaults$mp, min_val = 0.1)
    sensor_width_mm  <- read_num_safe("Sensor width (mm), e.g., 13.2: ", default = defaults$sensor_width_mm,  min_val = 1e-6)
    sensor_height_mm <- read_num_safe("Sensor height (mm), e.g., 8.8: ", default = defaults$sensor_height_mm, min_val = 1e-6)
  }
  
  # Pixel pitch (assumes square pixels)
  total_pixels    <- mp * 1e6
  sensor_area_mm2 <- sensor_width_mm * sensor_height_mm
  pixel_area_mm2  <- sensor_area_mm2 / total_pixels
  pixel_pitch_mm  <- sqrt(pixel_area_mm2)
  pixel_pitch_um  <- pixel_pitch_mm * 1000
  
  cat(sprintf("\nComputed pixel pitch ≈ %.3f µm (%.6f mm/pixel)\n",
              pixel_pitch_um, pixel_pitch_mm)); flush.console()
  
  # Ranged inputs
  focal_params <- read_seq_params_safe(
    name = "Focal length", unit = "mm",
    defaults = list(start = defaults$focal_start_mm,
                    end   = defaults$focal_end_mm,
                    by    = defaults$focal_by_mm)
  )
  alt_params <- read_seq_params_safe(
    name = "Altitude", unit = "ft",
    defaults = list(start = defaults$alt_start_ft,
                    end   = defaults$alt_end_ft,
                    by    = defaults$alt_by_ft)
  )
  
  focal_list_mm <- make_seq_inclusive(focal_params$start, focal_params$end, focal_params$by)
  alt_list_ft   <- make_seq_inclusive(alt_params$start,   alt_params$end,   alt_params$by)
  
  if (length(focal_list_mm) < 2) stop("Focal list must have at least 2 classes.")
  if (length(alt_list_ft)   < 2) stop("Altitude list must have at least 2 classes.")
  
  # Convert units
  ft_to_m    <- 0.3048
  m_to_ft    <- 3.28083989501312
  alt_list_m <- alt_list_ft * ft_to_m
  
  # ============ COMPUTATIONS ============
  # 1) GSD (cm/pixel)
  gsd_m  <- outer(alt_list_m, focal_list_mm, function(H, f) H * (pixel_pitch_mm / f))
  gsd_cm <- gsd_m * 100
  
  # 2) Ground coverage (whole frame footprint)
  width_m   <- outer(alt_list_m, focal_list_mm, function(H, f) H * (sensor_width_mm  / f))
  length_m  <- outer(alt_list_m, focal_list_mm, function(H, f) H * (sensor_height_mm / f))
  area_m2   <- width_m * length_m
  
  width_ft  <- width_m  * m_to_ft
  length_ft <- length_m * m_to_ft
  area_ft2  <- width_ft * length_ft
  
  # Labels (rows = altitude, cols = focal)
  rn <- paste0(sprintf("%.2f", alt_list_ft), " ft (", sprintf("%.2f", alt_list_m), " m)")
  cn <- paste0(sprintf("%.2f", focal_list_mm), " mm")
  dimnames(gsd_cm)    <- list(rn, cn)
  dimnames(width_m)   <- list(rn, cn)
  dimnames(length_m)  <- list(rn, cn)
  dimnames(area_m2)   <- list(rn, cn)
  dimnames(width_ft)  <- list(rn, cn)
  dimnames(length_ft) <- list(rn, cn)
  dimnames(area_ft2)  <- list(rn, cn)
  
  # ============ OUTPUT CHOICE ============
  cat("\nOutput type:\n  1 = GSD (cm/pixel)\n  2 = Coverage (Width x Length) + Area\n  3 = Both\n")
  choice <- suppressWarnings(as.integer(read_num_safe("Choose [3]: ", default = 3, min_val = 1, max_val = 3)))
  show_gsd <- choice %in% c(1, 3)
  show_cov <- choice %in% c(2, 3)
  
  # NEW: choose coverage units (only if coverage is requested)
  cov_units <- "ft"   # default
  if (show_cov) {
    u <- tolower(get_input_safe("\nCoverage unit? (ft/m) [ft]: ", default = "ft"))
    if (u %in% c("m","meter","meters","metre","metres")) cov_units <- "m"
  }
  cov_units_label <- if (cov_units == "m") "m" else "ft"
  
  # Pick matrices by unit
  if (cov_units == "m") {
    width_out  <- width_m
    length_out <- length_m
    area_out   <- area_m2
    cov_title_units <- "meters"
    area_title_units <- "m^2"
    cov_suffix <- "_m"
    area_suffix <- "_m2"
    round_wl <- 2; round_area <- 2
  } else {
    width_out  <- width_ft
    length_out <- length_ft
    area_out   <- area_ft2
    cov_title_units <- "feet"
    area_title_units <- "ft^2"
    cov_suffix <- "_ft"
    area_suffix <- "_ft2"
    round_wl <- 2; round_area <- 1
  }
  
  # One "Width x Length" character matrix (in chosen units)
  cov_wxl <- matrix(
    paste0(round(width_out, round_wl), " x ", round(length_out, round_wl)),
    nrow = nrow(width_out), ncol = ncol(width_out),
    dimnames = list(rn, cn)
  )
  
  # Display
  if (show_gsd) {
    cat("\n=== Ground Sample Distance (cm/pixel) ===\n"); flush.console()
    print(round(gsd_cm, 3))
    cat(sprintf("\nGSD range: %.3f to %.3f cm/pixel\n", min(gsd_cm), max(gsd_cm))); flush.console()
  }
  if (show_cov) {
    cat(sprintf("\n=== Ground Coverage (Width x Length, %s) ===\n", cov_title_units)); flush.console()
    print(cov_wxl, quote = FALSE)
    cat(sprintf("\n-- Total Area (%s) --\n", area_title_units)); print(round(area_out, round_area))
    cat(sprintf("\nArea range: %.1f to %.1f %s\n",
                min(area_out), max(area_out), area_title_units)); flush.console()
  }
  
  # ============ EXCEL EXPORT ============
  if (ask_save) {
    save_xlsx <- tolower(get_input_safe("\nSave to Excel (.xlsx)? (y/n) [n]: ", default = "n"))
    if (save_xlsx %in% c("y", "yes")) {
      base <- get_input_safe("  Base file name (no extension, e.g., gsd_ranges): ",
                             default = "gsd_ranges")
      path <- paste0(base, ".xlsx")
      
      sheets <- list()
      sheets[["Parameters"]] <- data.frame(
        Camera_Model     = if (!is.null(cam_hit)) cam_hit$model else "Manual / Custom",
        Megapixels       = mp,
        Sensor_Width_mm  = sensor_width_mm,
        Sensor_Height_mm = sensor_height_mm,
        Pixel_Pitch_mm   = round(pixel_pitch_mm, 6),
        Altitude_Input   = "feet",
        Coverage_Units   = cov_units_label,
        stringsAsFactors = FALSE
      )
      if (show_gsd) {
        df_gsd <- data.frame(
          Altitude_ft = alt_list_ft, Altitude_m = alt_list_m,
          round(gsd_cm, 4), check.names = FALSE
        )
        names(df_gsd)[-(1:2)] <- paste0("GSD_cm_per_px @ ", sprintf("%.2f", focal_list_mm), " mm")
        sheets[["GSD_cm_per_px"]] <- df_gsd
      }
      if (show_cov) {
        df_cov <- data.frame(
          Altitude_ft = alt_list_ft, Altitude_m  = alt_list_m,
          as.data.frame(cov_wxl, stringsAsFactors = FALSE), check.names = FALSE
        )
        names(df_cov)[-(1:2)] <- paste0("Width x Length (", cov_units_label, ") @ ", sprintf("%.2f", focal_list_mm), " mm")
        sheets[[paste0("Coverage_WidthxLength", cov_suffix)]] <- df_cov
        
        df_w <- data.frame(Altitude_ft = alt_list_ft, Altitude_m = alt_list_m,
                           round(width_out, round_wl), check.names = FALSE)
        names(df_w)[-(1:2)] <- paste0("Width_", cov_units_label, " @ ", sprintf("%.2f", focal_list_mm), " mm")
        sheets[[paste0("Width", cov_suffix)]] <- df_w
        
        df_l <- data.frame(Altitude_ft = alt_list_ft, Altitude_m = alt_list_m,
                           round(length_out, round_wl), check.names = FALSE)
        names(df_l)[-(1:2)] <- paste0("Length_", cov_units_label, " @ ", sprintf("%.2f", focal_list_mm), " mm")
        sheets[[paste0("Length", cov_suffix)]] <- df_l
        
        df_a <- data.frame(Altitude_ft = alt_list_ft, Altitude_m = alt_list_m,
                           round(area_out, round_area), check.names = FALSE)
        names(df_a)[-(1:2)] <- paste0("Area", area_suffix, " @ ", sprintf("%.2f", focal_list_mm), " mm")
        sheets[[paste0("Area", area_suffix)]] <- df_a
      }
      
      .save_excel(path, sheets)
      cat(sprintf("  -> Saved Excel workbook: %s\n", path)); flush.console()
    }
  }
  
  # ============ HEATMAPS (interactive loop) ============
  if (ask_plot) {
    open_menu <- tolower(get_input_safe("\nOpen heatmap menu now? (y/n) [n]: ", default = "n"))
    if (open_menu %in% c("y","yes")) {
      .show_heatmap_loop(show_gsd, show_cov, gsd_cm, width_out, length_out, area_out, cov_units_label)
    }
  }
  
  cat("\nDone.\n"); flush.console()
  
  invisible(list(
    camera_hit       = if (!is.null(cam_hit)) cam_hit$model else NULL,
    pixel_pitch_mm   = pixel_pitch_mm,
    pixel_pitch_um   = pixel_pitch_um,
    focal_list_mm    = focal_list_mm,
    alt_list_ft      = alt_list_ft,
    alt_list_m       = alt_list_m,
    gsd_cm           = gsd_cm,
    # coverage outputs in the unit chosen by the user:
    coverage_units   = cov_units_label,
    width            = width_out,
    length           = length_out,
    area             = area_out,
    coverage_wxl     = cov_wxl
  ))
}

# ---------- Auto-run control ----------
if (isTRUE(getOption("gsd.autorun", TRUE))) {
  invisible(run_gsd_ranges())
}
