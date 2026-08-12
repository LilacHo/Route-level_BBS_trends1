# 3c_add_SDM_covariates.R
#
# Covariate-pipeline counterpart of 3_add_SDM.R — identical logic, different
# input/output directories, so it can run against the per-species-per-model
# CSVs from 2c_generate_route_trend_csvs_covariates.R without touching
# 3_add_SDM.R or colliding with its output/species_routes /
# output/species_routes_sdm directories (which belong to the original,
# non-covariate pipeline: 1_species_iCAR_2010_2025.R ->
# 2_generate_route_trend_csvs.R -> 3_add_SDM.R).
#
# Add rcp45 and rcp85 climate-scenario columns to each per-species-per-model
# route CSV written by 2c_generate_route_trend_csvs_covariates.R
# (output/species_routes_covariates/per_species/) by extracting SDM
# classified-change raster values at each route location, and write the
# revised CSVs to output/species_routes_covariates/per_species_sdm/.
#
# Each input file already has exactly one species_code (2c writes one file
# per species-per-model-tag), which is what this script requires — the
# combined all_route_trends_*.csv from 2c (many species AND many models in
# one file) would NOT work here directly, hence the per-species files.
# Files also carry "group" and "model" columns through unchanged (this
# script doesn't look at "model", just passes it along like every other
# pre-existing column; "group" is used below to resolve raster paths, and
# also carried through so it doesn't need to be re-joined later).
#
# Since 1c_species_iCAR_covariates.R / 2c now fit and combine EVERY species
# across all 12 groups in one run (not one bird_group at a time), this
# script will encounter files from every group in a single run too, and
# needs rcp45/rcp85 rasters for whichever groups' species it's asked to
# process. A missing raster for one species/group is common when only some
# groups' SDM rasters have been generated so far — this script SKIPS that
# species (with a warning) rather than aborting the whole run, so it can
# make progress on whatever groups do have rasters available while you
# backfill the rest.
#
# Workflow (per *_route_trends.csv file):
#   1. Skip the file if its *_route_trends_sdm.csv output already exists.
#   2. Read the species code (species_code column, one code per file) and
#      group (group column, written directly by 2c — no lookup needed).
#   3. Build the rcp45/rcp85 raster paths from the group, which drives both the
#      directory (rcp45_<group>/, rcp85_<group>/) and the file name prefix
#      (<group>_<code>_breeding_2025_<45|85>_ENSEMBLE_classifiedchange.tif).
#   4. If either the group's raster folder or the species' raster file is
#      missing, warn and skip this file (leaves rcp45/rcp85 as NA rather than
#      writing anything) — does NOT stop the whole run.
#   5. Otherwise extract raster values at each route's longitude/latitude
#      (reprojected to the raster CRS) into the rcp45 and rcp85 columns and
#      write the result to output/species_routes_covariates/per_species_sdm/.
#
# Raster value legend (from Bateman et al. 2020):
#   0 = never suitable
#   1 = extirpation
#   2 = worsening (-50 to -100% change)
#   3 = slightly worsening (-25 to -50% change)
#   4 = neutral (-25 to +25% change)
#   5 = slightly improving (25 to 50% change)
#   6 = improving (50 to Inf % change)
#   7 = colonization

library(here)
library(tidyverse)
library(terra)
library(sf)

here::i_am("3c_add_SDM_covariates.R")

# List all species-per-model route CSVs -------------------------------------
# (written by 2c_generate_route_trend_csvs_covariates.R)
species_routes_dir <- here::here("output", "species_routes_covariates", "per_species")
species_files <- list.files(species_routes_dir, pattern = "_route_trends\\.csv$",
                            full.names = TRUE)

cat("Found", length(species_files), "species-per-model route files\n\n")

# Create output directory --------------------------------------------------
out_dir <- here::here("output", "species_routes_covariates", "per_species_sdm")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

n_skipped_missing_raster <- 0

# Process each species-per-model file ---------------------------------------
for (f in species_files) {

  # Skip if output already exists
  out_file <- file.path(out_dir, sub("_route_trends\\.csv$", "_route_trends_sdm.csv", basename(f)))
  if (file.exists(out_file)) {
    cat("Skipping (already exists):", basename(out_file), "\n")
    next
  }

  sp_routes <- read.csv(f)

  # Get the species abbreviation from the CSV (should be unique per file)
  abbr <- unique(sp_routes$species_code)
  if (length(abbr) != 1) {
    message("WARNING: Multiple species codes in ", basename(f), " — skipping")
    next
  }

  # Land-cover group is written directly by 2c_generate_route_trend_csvs_covariates.R
  # (the "group" column) -- no separate lookup needed.
  if (!"group" %in% names(sp_routes)) {
    message("WARNING: No 'group' column in ", basename(f),
            " — re-run 2c_generate_route_trend_csvs_covariates.R to regenerate it. Skipping.")
    next
  }
  bird_group <- unique(sp_routes$group)
  if (length(bird_group) != 1 || is.na(bird_group) || bird_group == "") {
    message("WARNING: No unique group value in ", basename(f), " — skipping")
    next
  }

  cat("Processing:", abbr, "[", bird_group, "] (", basename(f), ")\n")

  # Initialize new columns

  sp_routes$rcp45 <- NA
  sp_routes$rcp85 <- NA

  # Build raster directory + file paths
  rcp45_dir  <- here::here("data", paste0("rcp45_", bird_group))
  rcp85_dir  <- here::here("data", paste0("rcp85_", bird_group))
  rcp45_file <- file.path(rcp45_dir, abbr,
                          paste0(bird_group, "_", abbr, "_breeding_2025_45_ENSEMBLE_classifiedchange.tif"))
  rcp85_file <- file.path(rcp85_dir, abbr,
                          paste0(bird_group, "_", abbr, "_breeding_2025_85_ENSEMBLE_classifiedchange.tif"))

  # Require both raster files; WARN and skip this file (not stop the whole
  # run) if either is missing -- with ~600 species across 12 groups, it's
  # expected that only some groups' SDM rasters exist at any given time.
  # Distinguish "whole group folder missing" from "just this species' file
  # missing" so the message points at the right thing to fix.
  missing_msg <- NULL
  if (!dir.exists(rcp45_dir)) {
    missing_msg <- paste0("rcp45 raster folder not found: ", rcp45_dir,
                          " (species ", abbr, ", group '", bird_group, "')")
  } else if (!file.exists(rcp45_file)) {
    missing_msg <- paste0("rcp45 raster not found for ", abbr, ": ", rcp45_file)
  } else if (!dir.exists(rcp85_dir)) {
    missing_msg <- paste0("rcp85 raster folder not found: ", rcp85_dir,
                          " (species ", abbr, ", group '", bird_group, "')")
  } else if (!file.exists(rcp85_file)) {
    missing_msg <- paste0("rcp85 raster not found for ", abbr, ": ", rcp85_file)
  }
  if (!is.null(missing_msg)) {
    message("WARNING: ", missing_msg, " — skipping this species/model.")
    n_skipped_missing_raster <- n_skipped_missing_raster + 1
    next
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

  # Write revised CSV to output/species_routes_covariates/per_species_sdm/
  write.csv(sp_routes, out_file, row.names = FALSE)

  cat("  Wrote:", basename(out_file),
      "| rcp45 non-NA:", sum(!is.na(sp_routes$rcp45)),
      "| rcp85 non-NA:", sum(!is.na(sp_routes$rcp85)), "\n")
}

# Summary ------------------------------------------------------------------
cat("\n=== SDM extraction done (covariate pipeline) ===\n")
cat("Output directory:", out_dir, "\n")
cat("Files written:", length(list.files(out_dir, pattern = "\\.csv$")), "\n")
cat("Files skipped for missing rcp45/rcp85 rasters:", n_skipped_missing_raster,
    "(re-run this script once those groups' rasters are available)\n")
