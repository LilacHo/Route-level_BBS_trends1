# ==========================================================================
# tests/check_sdm.R
# ==========================================================================
# Diagnostic for the rcp45 x10000 scaling bug found while building the
# visualization dashboard (see conversation / PR notes): 5 species show
# rcp45 values of 0, 10000, 20000, ... 60000 in
# output/species_routes_covariates/per_species_sdm/*_route_trends_sdm.csv
# instead of the expected 0-7 SDM classified-change codes. rcp85 for the
# same species looks fine (plain 0-7), and the range-shift map PNGs for
# these species (output/species_maps/) also look fine, which suggests the
# problem is specific to the rcp45 .tif source rasters for these species,
# not to 3c_add_SDM_covariates.R's extraction code (which does no
# rescaling itself: extract() -> vals[, 2], straight assignment).
#
# This script inspects the raw rcp45 (and rcp85, as a control) rasters for
# the affected species -- plus one known-good species from the same group
# as a baseline -- and reports:
#   - datatype, scale/offset (scoff), whether the layer is categorical
#     (is.factor / cats / coltab), min/max
#   - a frequency table of raw cell values
#   - a re-run of the exact extract() logic 3c uses, on that species' own
#     BBS route coordinates, so we can see whether the x10000 values come
#     out of extract() (i.e. really are the raw pixel values) or only
#     appear later in the CSV (which would point elsewhere).
#
# No files are written; read-only. Run from the project root:
#   Rscript tests/check_sdm.R
# ==========================================================================

library(here)
library(terra)
library(dplyr)

here::i_am("tests/check_sdm.R")

# ---------------------------------------------------------------------------
# Species to inspect: the 5 flagged as having corrupted rcp45 values, plus
# one control species (Yellow-rumped Warbler, MYWA -- swap for any western_forests
# species you know is clean) from the same group for comparison.
# ---------------------------------------------------------------------------
suspects <- tribble(
  ~code,  ~group,           ~label,
  "ACWO", "western_forests", "Acorn Woodpecker (suspect)",
  "AMDI", "western_forests", "American Dipper (suspect)",
  "BHGR", "western_forests", "Black-headed Grosbeak (suspect)",
  "WILL", "waterbirds",      "Willet (suspect)",
  "WODU", "waterbirds",      "Wood Duck (suspect)"
)

control <- tribble(
  ~code,  ~group,           ~label,
  "BRCR", "western_forests", "control (western_forests)",
  "MALL", "waterbirds",      "control (waterbirds)"
)

sdm_dir <- function(rcp, group) here::here("data", paste0(rcp, "_", group))
sdm_file <- function(rcp, group, code) {
  rcp_num <- sub("^rcp", "", rcp)
  file.path(sdm_dir(rcp, group), code,
            paste0(group, "_", code, "_breeding_2025_", rcp_num, "_ENSEMBLE_classifiedchange.tif"))
}

# ---------------------------------------------------------------------------
# Per-raster metadata + raw value frequency table
# ---------------------------------------------------------------------------
inspect_raster <- function(path, rcp, code, label) {
  cat("\n----", label, "|", rcp, "----\n")
  cat("File:", path, "\n")
  if (!file.exists(path)) {
    cat("  ** FILE NOT FOUND **\n")
    return(invisible(NULL))
  }
  r <- rast(path)
  cat("  dim:", paste(dim(r), collapse = " x "),
      "| res:", paste(round(res(r), 2), collapse = ", "),
      "| crs:", crs(r, describe = TRUE)$name, "\n")
  cat("  datatype:", datatype(r), "\n")
  so <- scoff(r)
  cat("  scale/offset:", if (is.null(so)) "none" else paste(so, collapse = " / "), "\n")
  cat("  is.factor:", is.factor(r), "\n")
  if (is.factor(r)) {
    cat("  levels:\n")
    print(levels(r))
  }
  ct <- coltab(r)
  cat("  has color table:", !is.null(ct[[1]]), "\n")
  mm <- minmax(r)
  cat("  min/max:", paste(mm, collapse = " / "), "\n")

  cat("  frequency table of raw cell values:\n")
  ft <- tryCatch(freq(r), error = function(e) { cat("    freq() failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ft)) print(as.data.frame(ft))

  invisible(r)
}

# ---------------------------------------------------------------------------
# Reproduce 3c_add_SDM_covariates.R's extract() step for one species/rcp,
# using that species' own route coordinates from the pre-SDM CSV
# (output/species_routes_covariates/per_species/<name>_<model>_route_trends.csv).
# Falls back to the post-SDM CSV's lat/lon if the pre-SDM one isn't found.
# ---------------------------------------------------------------------------
reproduce_extract <- function(r, code, rcp) {
  pre_dir <- here::here("output", "species_routes_covariates", "per_species")
  candidates <- list.files(pre_dir, pattern = paste0("_base_route_trends\\.csv$"), full.names = TRUE)
  # Match by species_code inside the file rather than guessing the filename
  # (names have apostrophes/hyphens that don't map 1:1 to the code).
  hit <- NULL
  for (f in candidates) {
    first_line <- tryCatch(read.csv(f, nrows = 1), error = function(e) NULL)
    if (!is.null(first_line) && !is.null(first_line$species_code) &&
        first_line$species_code == code) {
      hit <- f
      break
    }
  }
  if (is.null(hit)) {
    cat("  (could not find pre-SDM route CSV for", code, "-- skipping extract() reproduction)\n")
    return(invisible(NULL))
  }
  routes <- read.csv(hit)
  routes_sv <- vect(routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  routes_sv_proj <- project(routes_sv, crs(r))
  vals <- extract(r, routes_sv_proj)
  cat("  extract() on", nrow(routes), "of", code, "'s own routes -- value table:\n")
  print(table(vals[, 2], useNA = "ifany"))
  invisible(vals)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
run_one <- function(code, group, label) {
  cat("\n==========================================================\n")
  cat(label, "(", code, "/", group, ")\n")
  cat("==========================================================\n")
  for (rcp in c("rcp45", "rcp85")) {
    path <- sdm_file(rcp, group, code)
    r <- inspect_raster(path, rcp, code, label)
    if (!is.null(r)) reproduce_extract(r, code, rcp)
  }
}

for (i in seq_len(nrow(suspects))) {
  run_one(suspects$code[i], suspects$group[i], suspects$label[i])
}
for (i in seq_len(nrow(control))) {
  run_one(control$code[i], control$group[i], control$label[i])
}

cat("\n\n=== Done. Compare suspects' rcp45 frequency tables/extract() output against",
    "their own rcp85 and against the control species' rcp45. ===\n")
