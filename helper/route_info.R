## =============================================================================
## helper/route_info.R
##
## Extract the lightweight route lookup (route, routeF, latitude, longitude)
## that 2c_generate_route_trend_csvs_covariates.R (and 2_generate_route_trend_csvs.R)
## need, directly from an EXISTING data/stan_data/*_stan_data.RData file —
## without re-running 1c_species_iCAR_covariates.R or
## 1_species_iCAR_2010_2025.R, i.e. without re-fitting any models.
##
## Use this to backfill data/route_info/*.rds for any stan_data.RData that
## doesn't yet have a matching route_info.rds — e.g. runs from before 1c
## started saving route_info.rds directly, or any run where 1c's
## stan_data.RData save succeeded but its route_info.rds save didn't (see the
## unlink()-vs-save() ordering issue discussed for
## 1c_species_iCAR_covariates.R's fit_one_covariate_model()).
##
## Two known stan_data.RData shapes, both handled:
##   - 1c_species_iCAR_covariates.R / tests/test_16covariates_pilot.R:
##     saves `new_data`, a plain data.frame with route/routeF/latitude/
##     longitude columns directly.
##   - 1_species_iCAR_2010_2025.R: saves `route_map`, an sf object (route,
##     routeF, strat, geometry) reprojected to the strata map's CRS — needs
##     st_transform(4326) + st_coordinates() to recover lat/lon, the same way
##     2_generate_route_trend_csvs.R already does.
##
## Usage:
##   - source("helper/route_info.R") to get the extract_route_info() function
##     for use on a single file, e.g.:
##       extract_route_info(here::here("data", "stan_data", "X_stan_data.RData"),
##                          here::here("data", "route_info", "X_route_info.rds"))
##   - Rscript helper/route_info.R (or source() from the console) also runs
##     the batch section below, which backfills every
##     data/stan_data/*_stan_data.RData file that's missing a
##     data/route_info/*_route_info.rds counterpart. Safe to re-run — it
##     skips files that already have a route_info.rds unless overwrite = TRUE.
## =============================================================================

library(dplyr)
library(sf)
library(here)

here::i_am("helper/route_info.R")

#' Extract route_info (route, routeF, latitude, longitude) from a saved
#' stan_data.RData file, without touching the caller's global environment and
#' without re-running any model fitting.
#'
#' @param stan_data_file path to a *_stan_data.RData file
#' @param route_info_file path to write the resulting route_info .rds to
#' @param overwrite if FALSE (default), skip if route_info_file already exists
#' @return TRUE if a route_info.rds was written (or already existed and
#'   overwrite = FALSE), FALSE if extraction failed.
extract_route_info <- function(stan_data_file, route_info_file, overwrite = FALSE) {
  if (file.exists(route_info_file) && !overwrite) {
    cat("  Skipping (already exists):", basename(route_info_file), "\n")
    return(TRUE)
  }
  if (!file.exists(stan_data_file)) {
    message("  [WARN] stan_data file not found: ", stan_data_file)
    return(FALSE)
  }

  # Load into an isolated environment — never touches the caller's global
  # environment, and stan_data.RData can contain very large objects.
  e <- new.env()
  loaded <- tryCatch(load(stan_data_file, envir = e), error = function(err) {
    message("  [ERROR] Failed to load ", basename(stan_data_file), ": ", conditionMessage(err))
    NULL
  })
  if (is.null(loaded)) return(FALSE)

  route_info <- NULL

  if ("new_data" %in% loaded && is.data.frame(e$new_data) &&
      all(c("route", "routeF", "latitude", "longitude") %in% names(e$new_data))) {
    # 1c_species_iCAR_covariates.R / test_16covariates_pilot.R shape.
    route_info <- e$new_data %>% distinct(route, routeF, latitude, longitude)

  } else if ("route_map" %in% loaded && inherits(e$route_map, "sf")) {
    # 1_species_iCAR_2010_2025.R shape — route_map is an sf object already
    # reprojected to the strata map's CRS (not WGS84), so convert back the
    # same way 2_generate_route_trend_csvs.R does.
    route_info <- tryCatch({
      e$route_map %>%
        st_transform(4326) %>%
        mutate(longitude = st_coordinates(.)[, 1],
               latitude  = st_coordinates(.)[, 2]) %>%
        st_drop_geometry() %>%
        distinct(route, routeF, latitude, longitude)
    }, error = function(err) {
      message("  [ERROR] Failed to derive route_info from route_map in ",
              basename(stan_data_file), ": ", conditionMessage(err))
      NULL
    })
  }

  if (is.null(route_info)) {
    message("  [WARN] No usable new_data/route_map (with route/routeF/latitude/longitude) found in ",
            basename(stan_data_file), " — objects loaded: ", paste(loaded, collapse = ", "))
    return(FALSE)
  }

  dir.create(dirname(route_info_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(route_info, route_info_file)
  cat("  Wrote:", basename(route_info_file), "(", nrow(route_info), "routes )\n")
  TRUE
}

## ==========================================================================
## Batch mode: backfill every stan_data.RData in data/stan_data/ that's
## missing a route_info.rds in data/route_info/. Runs automatically whenever
## this file is executed top-to-bottom (Rscript helper/route_info.R, or
## source("helper/route_info.R")).
## ==========================================================================
stan_data_dir  <- here::here("data", "stan_data")
route_info_dir <- here::here("data", "route_info")

stan_data_files <- list.files(stan_data_dir, pattern = "_stan_data\\.RData$", full.names = TRUE)
cat("=== route_info backfill ===\n")
cat("Found", length(stan_data_files), "stan_data.RData file(s) in", stan_data_dir, "\n\n")

n_ok   <- 0
n_fail <- 0
for (f in stan_data_files) {
  route_info_file <- file.path(route_info_dir,
                               sub("_stan_data\\.RData$", "_route_info.rds", basename(f)))
  cat("[", basename(f), "]\n")
  ok <- extract_route_info(f, route_info_file)
  if (ok) n_ok <- n_ok + 1 else n_fail <- n_fail + 1
}

cat("\n=== Done ===\n")
cat("OK:", n_ok, " | Failed:", n_fail, "\n")
cat("route_info .rds files in:", route_info_dir, "\n")
