# 3c_add_SDM_covariates.R
#
# Adds climate-suitability columns to each per-species-per-model route CSV
# from 2c_generate_route_trend_csvs_covariates.R, by extracting species
# distribution model (SDM) raster values at each route's location. Each
# species' raster folder is resolved from its own bird group, so this runs
# across all 12 groups automatically. A species with a missing raster
# folder/file is skipped (with a warning) rather than stopping the run.
#
# Input: output/species_routes_covariates/per_species/*_route_trends.csv
# (from 2c), data/spp_names_codes_group_aou.csv (to resolve each species'
# group), and the SDM rasters in data/rcp45_<group>/ and data/rcp85_<group>/
# for every group (file naming:
# <group>_<species_code>_breeding_2025_<45|85>_ENSEMBLE_classifiedchange.tif).
#
# Output: output/species_routes_covariates/per_species_sdm/<species>_<tag>_route_trends_sdm.csv
# -- same columns as the input, plus rcp45 and rcp85 (raster values at each
# route, reprojected to the raster CRS).
#
# Raster value legend (Bateman et al. 2020):
#   0 = never suitable   1 = extirpation      2 = worsening (-50 to -100%)
#   3 = slightly worsening (-25 to -50%)      4 = neutral (-25 to +25%)
#   5 = slightly improving (25 to 50%)        6 = improving (50%+)
#   7 = colonization

library(here)
library(tidyverse)
library(terra)
library(sf)

here::i_am("3c_add_SDM_covariates.R")

# Read species name/code lookup table --------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"))

# List all species-per-model route CSVs (written by 2c) ---------------------
species_routes_dir <- here::here("output", "species_routes_covariates", "per_species")
species_files <- list.files(species_routes_dir, pattern = "_route_trends\\.csv$",
                            full.names = TRUE)

cat("Found", length(species_files), "species-per-model route files\n\n")

# Create output directory --------------------------------------------------
out_dir <- here::here("output", "species_routes_covariates", "per_species_sdm")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_skipped_missing_group   <- 0
n_skipped_missing_rasters <- 0
n_processed               <- 0

# Process each species-per-model file ---------------------------------------
for (f in species_files) {

  # Skip if output already exists
  out_file <- file.path(out_dir, sub("_route_trends\\.csv$", "_route_trends_sdm.csv", basename(f)))
  if (file.exists(out_file)) {
    cat("Skipping (already exists):", basename(out_file), "\n")
    next
  }

  sp_routes <- read.csv(f)

  # Get the species code from the CSV (should be unique per file)
  abbr <- unique(sp_routes$species_code)
  if (length(abbr) != 1) {
    message("WARNING: Multiple species codes in ", basename(f), " — skipping")
    next
  }

  if (!abbr %in% spp_df$Code) {
    message("WARNING: Code '", abbr, "' not found in spp_names_codes_group_aou.csv — skipping")
    n_skipped_missing_group <- n_skipped_missing_group + 1
    next
  }

  # Species' bird group drives the raster directory + file prefix
  bird_group <- unique(spp_df$Group[spp_df$Code == abbr])
  if (length(bird_group) != 1 || is.na(bird_group) || bird_group == "") {
    message("WARNING: No unique Group for code '", abbr, "' — skipping")
    n_skipped_missing_group <- n_skipped_missing_group + 1
    next
  }

  cat("Processing:", abbr, "[", bird_group, "] (", basename(f), ")\n")

  sp_routes$rcp45 <- NA
  sp_routes$rcp85 <- NA

  rcp45_dir  <- here::here("data", paste0("rcp45_", bird_group))
  rcp85_dir  <- here::here("data", paste0("rcp85_", bird_group))
  rcp45_file <- file.path(rcp45_dir, abbr,
                          paste0(bird_group, "_", abbr, "_breeding_2025_45_ENSEMBLE_classifiedchange.tif"))
  rcp85_file <- file.path(rcp85_dir, abbr,
                          paste0(bird_group, "_", abbr, "_breeding_2025_85_ENSEMBLE_classifiedchange.tif"))

  # Skip (don't stop the whole run) when either raster is missing.
  missing <- c(
    if (!dir.exists(rcp45_dir))   paste("rcp45 folder not found:", rcp45_dir),
    if (dir.exists(rcp45_dir) && !file.exists(rcp45_file)) paste("rcp45 raster not found:", rcp45_file),
    if (!dir.exists(rcp85_dir))   paste("rcp85 folder not found:", rcp85_dir),
    if (dir.exists(rcp85_dir) && !file.exists(rcp85_file)) paste("rcp85 raster not found:", rcp85_file)
  )
  if (length(missing) > 0) {
    message("WARNING: Skipping ", abbr, " [", bird_group, "] — ", paste(missing, collapse = "; "))
    n_skipped_missing_rasters <- n_skipped_missing_rasters + 1
    next
  }

  # Extract rcp45 values
  rcp45_rast <- rast(rcp45_file)
  routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  routes_sv_proj <- project(routes_sv, crs(rcp45_rast))
  vals_45 <- extract(rcp45_rast, routes_sv_proj)
  sp_routes$rcp45 <- vals_45[, 2]

  # Extract rcp85 values
  rcp85_rast <- rast(rcp85_file)
  routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
  routes_sv_proj <- project(routes_sv, crs(rcp85_rast))
  vals_85 <- extract(rcp85_rast, routes_sv_proj)
  sp_routes$rcp85 <- vals_85[, 2]

  write.csv(sp_routes, out_file, row.names = FALSE)
  n_processed <- n_processed + 1

  cat("  Wrote:", basename(out_file),
      "| rcp45 non-NA:", sum(!is.na(sp_routes$rcp45)),
      "| rcp85 non-NA:", sum(!is.na(sp_routes$rcp85)), "\n")
}

# Summary ------------------------------------------------------------------
cat("\n=== SDM extraction done (all groups, base + anthro) ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files written this run:", n_processed, "\n")
cat("Skipped (no Group match):", n_skipped_missing_group, "\n")
cat("Skipped (missing raster folder/file):", n_skipped_missing_rasters, "\n")
cat("Total SDM CSVs on disk now:", length(list.files(out_dir, pattern = "\\.csv$")), "\n")
