# 3c_add_SDM_covariates.R
#
# Adds climate-suitability columns to each per-species-per-model route CSV
# from 2c, by extracting species distribution model (SDM) raster values at
# each route's location.
#
# Input: output/species_routes_covariates/per_species/*_route_trends.csv
# (from 2c), data/spp_names_codes_group_aou.csv (to resolve each species'
# land-cover group), and the SDM rasters in data/rcp45_grasslands/ and
# data/rcp85_grasslands/ (file naming:
# grasslands_<species_code>_breeding_2025_<45|85>_ENSEMBLE_classifiedchange.tif).
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
    next
  }

  # Species' land cover group drives the raster directory + file prefix
  land_cover <- unique(spp_df$Group[spp_df$Code == abbr])
  if (length(land_cover) != 1 || is.na(land_cover) || land_cover == "") {
    message("WARNING: No unique Group for code '", abbr, "' — skipping")
    next
  }

  cat("Processing:", abbr, "[", land_cover, "] (", basename(f), ")\n")

  sp_routes$rcp45 <- NA
  sp_routes$rcp85 <- NA

  rcp45_dir  <- here::here("data", paste0("rcp45_", land_cover))
  rcp85_dir  <- here::here("data", paste0("rcp85_", land_cover))
  rcp45_file <- file.path(rcp45_dir, abbr,
                          paste0(land_cover, "_", abbr, "_breeding_2025_45_ENSEMBLE_classifiedchange.tif"))
  rcp85_file <- file.path(rcp85_dir, abbr,
                          paste0(land_cover, "_", abbr, "_breeding_2025_85_ENSEMBLE_classifiedchange.tif"))

  # Require both raster files; stop with a specific message if either is missing.
  if (!dir.exists(rcp45_dir)) {
    stop("rcp45 raster folder not found: ", rcp45_dir,
         " (species ", abbr, ", group '", land_cover, "')")
  }
  if (!file.exists(rcp45_file)) {
    stop("rcp45 raster not found for ", abbr, ": ", rcp45_file)
  }
  if (!dir.exists(rcp85_dir)) {
    stop("rcp85 raster folder not found: ", rcp85_dir,
         " (species ", abbr, ", group '", land_cover, "')")
  }
  if (!file.exists(rcp85_file)) {
    stop("rcp85 raster not found for ", abbr, ": ", rcp85_file)
  }

  # Extract rcp45 values
  if (file.exists(rcp45_file)) {
    rcp45_rast <- rast(rcp45_file)
    routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    routes_sv_proj <- project(routes_sv, crs(rcp45_rast))
    vals_45 <- extract(rcp45_rast, routes_sv_proj)
    sp_routes$rcp45 <- vals_45[, 2]
  }

  # Extract rcp85 values
  if (file.exists(rcp85_file)) {
    rcp85_rast <- rast(rcp85_file)
    routes_sv <- vect(sp_routes, geom = c("longitude", "latitude"), crs = "EPSG:4326")
    routes_sv_proj <- project(routes_sv, crs(rcp85_rast))
    vals_85 <- extract(rcp85_rast, routes_sv_proj)
    sp_routes$rcp85 <- vals_85[, 2]
  }

  write.csv(sp_routes, out_file, row.names = FALSE)

  cat("  Wrote:", basename(out_file),
      "| rcp45 non-NA:", sum(!is.na(sp_routes$rcp45)),
      "| rcp85 non-NA:", sum(!is.na(sp_routes$rcp85)), "\n")
}

# Summary ------------------------------------------------------------------
cat("\n=== SDM extraction done (covariate pipeline) ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files written:", length(list.files(out_dir, pattern = "\\.csv$")), "\n")
