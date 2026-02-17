# ============================
# Global config (dynamic paths; no hard-coding)
# ============================

options(timeout = 600)

# Avoid s2 issues with invalid polygons in 'terre'
sf::sf_use_s2(FALSE)

# ---- helpers ---------------------------------------------------------------
.norm_path <- function(p) {
  if (!nzchar(p)) return("")
  normalizePath(p, winslash = "/", mustWork = FALSE)
}
.get_env <- function(var, default = "") {
  val <- Sys.getenv(var, unset = "")
  if (nzchar(val)) val else default
}

# ---- ROOT_DIR --------------------------------------------------------------
# Precedence:
# 1) Keep any ROOT_DIR pre-set (e.g., by main.r)
# 2) ENV ROOT_DIR
# 3) Fallback to current working directory

if (!exists("ROOT_DIR") || !nzchar(ROOT_DIR)) {
  env_root <- .get_env("ROOT_DIR")
  if (nzchar(env_root)) {
    ROOT_DIR <- .norm_path(env_root)
  } else {
    ROOT_DIR <- .norm_path(getwd())
  }
}

# ---- Base data/cache/outputs ------------------------------------------------
# Allow ENV overrides for each; else default relative to ROOT_DIR

# TERRE_PATH (can be file or directory; you use a .shp file)
TERRE_PATH <- if (!exists("TERRE_PATH") || !nzchar(TERRE_PATH)) {
  env_terre <- .get_env("TERRE_PATH")
  if (nzchar(env_terre)) .norm_path(env_terre) else
    file.path(ROOT_DIR, "data_downloads/background/terreNAD83.shp")
} else TERRE_PATH
TERRE_PATH <- .norm_path(TERRE_PATH)

# CACHE_DIR
CACHE_DIR <- if (!exists("CACHE_DIR") || !nzchar(CACHE_DIR)) {
  env_cache <- .get_env("CACHE_DIR")
  if (nzchar(env_cache)) .norm_path(env_cache) else
    file.path(ROOT_DIR, "data_downloads")
} else CACHE_DIR
CACHE_DIR <- .norm_path(CACHE_DIR)
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# AIRPORTS_CSV
AIRPORTS_CSV <- if (!exists("AIRPORTS_CSV") || !nzchar(AIRPORTS_CSV)) {
  env_air_csv <- .get_env("AIRPORTS_CSV")
  if (nzchar(env_air_csv)) .norm_path(env_air_csv) else
    file.path(CACHE_DIR, "airports.csv")
} else AIRPORTS_CSV
AIRPORTS_CSV <- .norm_path(AIRPORTS_CSV)

# OUTPUT_DIR (for generic quick outputs)
OUTPUT_DIR <- if (!exists("OUTPUT_DIR") || !nzchar(OUTPUT_DIR)) {
  env_out <- .get_env("OUTPUT_DIR")
  if (nzchar(env_out)) .norm_path(env_out) else
    file.path(ROOT_DIR, "data_downloads/outputs")
} else OUTPUT_DIR
OUTPUT_DIR <- .norm_path(OUTPUT_DIR)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# OUT_PNG (preview map path)
OUT_PNG <- if (!exists("OUT_PNG") || !nzchar(OUT_PNG)) {
  env_out_png <- .get_env("OUT_PNG")
  if (nzchar(env_out_png)) .norm_path(env_out_png) else
    file.path(OUTPUT_DIR, "map_bbox_preview.png")
} else OUT_PNG
OUT_PNG <- .norm_path(OUT_PNG)

# TRANSECT_EXPORT_DIR (root under which _transects/ and _flights/ folders are created)
# Precedence:
# 1) ENV TRANSECT_EXPORT_DIR
# 2) Keep existing in session
# 3) Default to ROOT_DIR/exports (simple, predictable)
TRANSECT_EXPORT_DIR <- (function(current_value) {
  env_dir <- .get_env("TRANSECT_EXPORT_DIR")
  if (nzchar(env_dir)) {
    .norm_path(env_dir)
  } else if (!missing(current_value) && !is.null(current_value) && nzchar(current_value)) {
    .norm_path(current_value)
  } else {
    file.path(ROOT_DIR, "exports")
  }
})(if (exists("TRANSECT_EXPORT_DIR")) TRANSECT_EXPORT_DIR else "")

TRANSECT_EXPORT_DIR <- .norm_path(TRANSECT_EXPORT_DIR)
dir.create(TRANSECT_EXPORT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Optional: one-time debug (comment out if noisy) -----------------------
 message("CONFIG >> ROOT_DIR: ", ROOT_DIR)
 message("CONFIG >> TERRE_PATH: ", TERRE_PATH)
 message("CONFIG >> CACHE_DIR: ", CACHE_DIR)
 message("CONFIG >> OUTPUT_DIR: ", OUTPUT_DIR)
 message("CONFIG >> TRANSECT_EXPORT_DIR: ", TRANSECT_EXPORT_DIR)