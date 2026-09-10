## =============================================================================
## functions/covariate_model_fitting.R
##
## Shared data-prep + fitting functions for the iCAR route-level slope model
## (base vs. anthro covariate), used by every script that fits or refits
## these models:
##   - 1c_species_iCAR_covariates.R        (production: all species x base+anthro)
##   - 1d_refit_nonconverged_species.R     (targeted refit of non-converged fits)
##   - tests/test_green_heron_debug.R      (single-species debug re-run)
##
## Consolidated here (previously hand-duplicated in all three files, which
## risked drifting out of sync) so there is exactly ONE implementation of
## the data-prep pipeline and the fitting call. Differences between
## production/refit/debug use are expressed entirely through function
## ARGUMENTS (chains, iterations, adapt_delta, show_exceptions, output
## directories) -- not through separate copies of the function bodies.
##
## Provides:
##   species_to_f(sp)          -- file-safe species name
##   load_covariate(...)       -- read one route-year covariate CSV
##   prepare_species_data(...) -- BBS pull + covariate join + NA-drop +
##                                 spatial (Voronoi) neighbours for one species
##   fit_one_covariate_model(...) -- fit "base" or a named covariate model tag
##                                 for one species, given prepare_species_data()'s
##                                 output; saves stanfit/summ/route_info/
##                                 stan_data to the directories passed in, and
##                                 returns a one-row diagnostics data.frame.
##
## fit_one_covariate_model() also includes a chain-completion gate: right
## after stan_model$sample() returns, it checks stanfit$return_codes()
## BEFORE touching $summary()/$time(). If any chain failed (nonzero return
## code), it prints a clear "N of 4 chains completed, chain X failed"
## message. This exists because a partially-failed chain otherwise surfaces
## only as a cryptic downstream R-level error from $summary()/$time() (e.g.
## "arguments imply differing number of rows: 4, 3") with no indication of
## which chain crashed or why -- this gate makes that diagnosable from
## console output alone, for every caller (production runs included, not
## just debug runs).
##
## This file sources functions/neighbours_define_voronoi.R and
## functions/posterior_summary_functions.R itself, so callers no longer need
## to source those separately before sourcing this file.
## =============================================================================

library(dplyr)
library(here)

source(here::here("functions", "neighbours_define_voronoi.R"))
source(here::here("functions", "posterior_summary_functions.R"))

# ---------------------------------------------------------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# ---------------------------------------------------------------------------
#' Read one route-year covariate CSV (one row per BBS route per year) and
#' build a "<StateNum>-<Route>" route_key to match new_data$route in
#' prepare_species_data() below.
#'
#' @param file CSV filename under data/ (e.g. "Anthro.csv")
#' @param value_col column name in the CSV holding the covariate value
#' @param out_name name to give the value column in the returned lookup
#'   (also doubles as the model_tag used by fit_one_covariate_model())
#' @return data.frame(route_key, year, <out_name>), one row per route-year,
#'   NA values dropped, deduplicated on route_key/year
load_covariate <- function(file, value_col, out_name) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              !!out_name := .data[[value_col]]) %>%
    filter(!is.na(.data[[out_name]])) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

# ---------------------------------------------------------------------------
#' Data-prep pipeline for one species: BBS data pull, covariate join,
#' NA-drop, spatial (Voronoi) neighbours. Generalized to an arbitrary NAMED
#' LIST of covariate lookups so it isn't tied to a fixed covariate set.
#' Counts with a missing value for ANY covariate in covariate_lookups are
#' dropped up front, so every model tag fit from the same prepped object
#' sees identical rows ("same data for a fair comparison").
#'
#' @param species Common.Name used for messages/plot labeling
#' @param species_bbs bbs_english name used by bbsBayes2::stratify()
#' @param strat stratification, e.g. "bcr"
#' @param firstYear,lastYear year range
#' @param covariate_lookups named list of load_covariate() outputs; names
#'   double as both the joined column name and the model_tag
#'   fit_one_covariate_model() will look for
#' @return list(new_data, route_map, realized_strata_map, car_stan_dat,
#'   sd_alpha_prior, base_stan_data, n_dropped)
prepare_species_data <- function(species, species_bbs, strat,
                                 firstYear, lastYear, covariate_lookups) {

  cat("    Preparing data for:", species, "\n")

  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data <- data_pkg[["raw_data"]]
  strata_map <- load_map(strat)

  # Continental US only -- the land-cover covariate CSVs only cover the
  # continental US, and state_num == 3 is excluded there too.
  raw_data <- raw_data %>%
    filter(country_num == 840) %>%
    filter(state_num != 3)

  new_data <- data.frame(
    strat_name = raw_data$strata_name,
    strat      = raw_data$strata,
    route      = raw_data$route,
    latitude   = raw_data$latitude,
    longitude  = raw_data$longitude,
    count      = raw_data$count,
    year       = raw_data$year_num,
    firstyr    = raw_data$first_year,
    ObsN       = raw_data$observer,
    r_year     = raw_data$year
  )

  if (nrow(new_data) < 50) {
    stop("Too few observations (", nrow(new_data), ")")
  }

  # Join every covariate by route + calendar year ---------------------------
  n_before <- nrow(new_data)
  cov_names <- names(covariate_lookups)
  for (cn in cov_names) {
    new_data <- new_data %>%
      left_join(covariate_lookups[[cn]], by = c("route" = "route_key", "r_year" = "year"))
  }

  match_rate <- mean(complete.cases(new_data[, cov_names, drop = FALSE]))
  cat("    Covariate match rate (all of", paste(cov_names, collapse = ", "), "present):",
      round(100 * match_rate, 1), "%\n")

  cat("    Sample new_data$route values:      ",
      paste(head(unique(new_data$route), 6), collapse = ", "), "\n")
  cat("    Sample covariate route_key values: ",
      paste(head(unique(covariate_lookups[[cov_names[1]]]$route_key), 6), collapse = ", "), "\n")

  if (match_rate < 0.5) {
    stop("Less than half of observations matched all of ", paste(cov_names, collapse = ", "),
         " for ", species, " — new_data$route doesn't look like '<StateNum>-<Route>'. ",
         "Compare the sample values printed above. This previously happened because strat ",
         "was 'bbs_usgs' instead of 'bcr'; if strat is already 'bcr', check the ",
         "continental-US filter and the key construction in load_covariate() instead.")
  }

  # Drop rows missing any covariate -- same reduced-dataset principle used
  # throughout this project, so every model tag (including "base") fits on
  # identical rows.
  new_data <- new_data %>% filter(complete.cases(new_data[, cov_names, drop = FALSE]))
  n_dropped <- n_before - nrow(new_data)
  cat("    Dropped", n_dropped, "of", n_before,
      "observations with missing covariate data\n")

  if (nrow(new_data) < 50) {
    stop("Too few observations after dropping missing-covariate rows (",
         nrow(new_data), ")")
  }

  realized_strata_map <- strata_map %>%
    filter(strata_name %in% unique(new_data$strat_name))

  # Spatial neighbours -----------------------------------------------------
  new_data$routeF <- as.integer(factor(new_data$route))

  route_map <- unique(data.frame(
    route     = new_data$route,
    routeF    = new_data$routeF,
    strat     = new_data$strat_name,
    latitude  = new_data$latitude,
    longitude = new_data$longitude
  ))

  dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
  while (length(dups) > 0) {
    route_map[dups, "latitude"]  <- route_map[dups, "latitude"]  + 0.01
    route_map[dups, "longitude"] <- route_map[dups, "longitude"] + 0.01
    dups <- which(duplicated(route_map[, c("latitude", "longitude")]))
  }

  route_map <- st_as_sf(route_map, coords = c("longitude", "latitude"))
  st_crs(route_map) <- 4326
  route_map <- st_transform(route_map, crs = st_crs(strata_map))

  car_stan_dat <- neighbours_define_voronoi(
    real_point_map  = route_map,
    species         = species,
    strat_indicator = "routeF",
    strata_map      = realized_strata_map,
    concavity       = 1,
    save_plot_data  = TRUE,
    plot_dir        = here::here("data", "maps")
  )

  cat("    Routes:", max(new_data$routeF),
      " | Obs:", nrow(new_data),
      " | Edges:", car_stan_dat$N_edges, "\n")

  tmp <- data.frame(route = new_data$routeF, count = new_data$count) %>%
    group_by(route) %>%
    summarise(mean_count = mean(log(count + 1)))
  sd_alpha_prior <- sd(tmp$mean_count)

  base_stan_data <- list(
    count      = new_data$count,
    ncounts    = nrow(new_data),
    strat      = new_data$strat,
    route      = new_data$routeF,
    year       = new_data$year,
    firstyr    = new_data$firstyr,
    fixedyear  = floor(mean(1:max(new_data$year))),
    nyears     = max(new_data$year),
    observer   = as.integer(factor(new_data$ObsN)),
    nobservers = length(unique(new_data$ObsN)),
    N_edges    = car_stan_dat$N_edges,
    node1      = car_stan_dat$node1,
    node2      = car_stan_dat$node2,
    nroutes    = max(new_data$routeF),
    sd_alpha_prior = sd_alpha_prior
  )

  if (car_stan_dat$N != base_stan_data[["nroutes"]]) {
    stop("Some routes are missing from adjacency matrix")
  }

  list(new_data = new_data, route_map = route_map,
       realized_strata_map = realized_strata_map, car_stan_dat = car_stan_dat,
       sd_alpha_prior = sd_alpha_prior, base_stan_data = base_stan_data,
       n_dropped = n_dropped)
}

# ---------------------------------------------------------------------------
#' Fit one model_tag ("base" or a covariate name present in
#' prepped$new_data, e.g. "anthro") for one species, given
#' prepare_species_data()'s output. Saves the stanfit, its posterior
#' summary, a lightweight route lookup, and the stan_data/new_data to the
#' directories passed in -- always at the SAME filenames
#' (<species_f>_iCAR_<model_tag>_<firstYear>_<lastYear>_*), regardless of
#' caller, so a refit or debug run overwrites/replaces a production run's
#' output in place rather than creating a parallel set of files.
#'
#' Includes a chain-completion gate: immediately after stan_model$sample()
#' returns, checks stanfit$return_codes() before calling $summary()/$time().
#' If fewer chains completed than requested, prints which chain(s) failed --
#' this is the most common root cause when $summary()/$time() throw a
#' cryptic "arguments imply differing number of rows" error (fewer completed
#' chains than the number configured).
#'
#' @param species,species_f Common.Name and file-safe name
#' @param model_tag "base" or a covariate name in prepped$new_data
#' @param prepped output of prepare_species_data()
#' @param firstYear,lastYear year range (used only for output filenames)
#' @param model_base,model_single compiled CmdStanModel objects for the
#'   no-covariate and single-covariate Stan programs (pass NULL for whichever
#'   isn't needed by the tags actually being fit)
#' @param rds_dir,route_info_dir,cmdstanr_output_dir,stan_data_dir output
#'   directories -- pass the project's production dirs to overwrite/extend a
#'   real run, or separate *_debug/*_test dirs to keep a run's output fully
#'   isolated
#' @param chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth,
#'   show_exceptions passed straight through to
#'   CmdStanModel$sample() -- production defaults match
#'   1c_species_iCAR_covariates.R's original settings; pass different values
#'   for a refit (more iterations, higher adapt_delta) or a debug run
#'   (show_exceptions = TRUE)
#' @param overwrite_route_info if TRUE (default), always (re)write
#'   route_info.rds; if FALSE, skip writing when it already exists (minor
#'   efficiency option for a large first-time production run where the file
#'   can't yet exist anyway)
#' @return a one-row diagnostics data.frame (species, model, max_rhat,
#'   min_ess, converged, gamma1_mean, n_obs, n_dropped_missing_cov, fit_min)
fit_one_covariate_model <- function(species, species_f, model_tag, prepped,
                                    firstYear, lastYear,
                                    model_base, model_single,
                                    rds_dir, route_info_dir, cmdstanr_output_dir,
                                    stan_data_dir = here::here("data", "stan_data"),
                                    chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                                    adapt_delta = 0.8, max_treedepth = 10,
                                    show_exceptions = FALSE,
                                    overwrite_route_info = TRUE) {

  new_data       <- prepped$new_data
  base_stan_data <- prepped$base_stan_data

  if (model_tag == "base") {
    stan_model <- model_base
    stan_data  <- base_stan_data
  } else {
    if (!model_tag %in% names(new_data)) {
      stop("model_tag '", model_tag, "' has no matching covariate column in new_data — ",
           "check covariate_specs.")
    }
    stan_model <- model_single
    stan_data  <- c(base_stan_data, list(covariate = new_data[[model_tag]]))
  }

  cat("    Fitting model:", model_tag,
      "(chains =", chains, ", iter_warmup =", iter_warmup, ", iter_sampling =", iter_sampling,
      ", adapt_delta =", adapt_delta, ")\n")

  # Per-species/model_tag output subfolder so raw cmdstan CSVs stay isolated
  # and can be deleted right after this fit is saved to .rds.
  run_output_dir <- file.path(cmdstanr_output_dir, paste0(species_f, "_", model_tag))
  if (!dir.exists(run_output_dir)) dir.create(run_output_dir, recursive = TRUE)

  stanfit <- stan_model$sample(
    data = stan_data,
    refresh = 400,
    chains = chains, iter_sampling = iter_sampling, iter_warmup = iter_warmup,
    parallel_chains = chains,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    show_exceptions = show_exceptions,
    output_dir = run_output_dir
  )

  # --- Chain-completion gate -- see function doc above. ---------------------
  codes <- tryCatch(stanfit$return_codes(), error = function(e) NA)
  if (!all(is.na(codes)) && any(codes != 0, na.rm = TRUE)) {
    n_ok <- sum(codes == 0, na.rm = TRUE)
    cat("    [WARNING] Only", n_ok, "of", length(codes), "chains completed successfully for",
        species, "(", model_tag, ") -- chain(s)",
        paste(which(codes != 0), collapse = ", "), "failed (nonzero return code).\n")
    cat("    This is very likely the root cause if summary()/time() error out next.",
        "Re-run with show_exceptions = TRUE to see the chain's own crash/exception text.\n")
  }

  summ <- tryCatch(stanfit$summary(), error = function(e) {
    stop("stanfit$summary() failed for ", species, " (", model_tag, "): ", conditionMessage(e),
         " -- see chain return codes printed above.")
  })
  fit_time <- tryCatch(round(stanfit$time()[["total"]] / 60, 1), error = function(e) {
    stop("stanfit$time() failed for ", species, " (", model_tag, "): ", conditionMessage(e),
         " -- see chain return codes printed above.")
  })
  cat("    Fit time:", fit_time, "minutes\n")

  out_base <- paste0(species_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear)
  stanfit$save_object(file.path(rds_dir, paste0(out_base, "_stanfit.rds")))
  saveRDS(summ, file.path(rds_dir, paste0(out_base, "_summ_fit.rds")))

  # Lightweight route lookup (route, routeF, latitude, longitude) — saved on
  # its own, before cleanup below, so 2c can get what it needs without
  # depending on the full stan_data.RData save further down.
  route_info_file <- file.path(route_info_dir,
                               paste0(species_f, "_", model_tag, "_",
                                      firstYear, "_", lastYear, "_route_info.rds"))
  if (overwrite_route_info || !file.exists(route_info_file)) {
    route_info <- new_data %>% distinct(route, routeF, latitude, longitude)
    saveRDS(route_info, route_info_file)
  }

  # Raw per-chain CSVs are no longer needed once the fit is saved above —
  # delete them now instead of letting them accumulate in %TEMP%.
  unlink(run_output_dir, recursive = TRUE)

  if (!dir.exists(stan_data_dir)) dir.create(stan_data_dir, recursive = TRUE)
  sp_data_file <- file.path(stan_data_dir,
                            paste0(species_f, "_", model_tag, "_",
                                   firstYear, "_", lastYear, "_stan_data.RData"))
  save(list = c("stan_data", "new_data", "model_tag"), file = sp_data_file)

  # Convergence diagnostics + covariate effect -------------------------------
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess  <- min(summ$ess_bulk, na.rm = TRUE)

  gamma1_mean <- if ("gamma1" %in% summ$variable) {
    summ$mean[summ$variable == "gamma1"]
  } else {
    NA_real_
  }

  # Whole-model convergence: "species are accepted in models converged based
  # on Rhat < 1.01 and bulk effective sample sizes > 400" -- same criterion
  # helper/model_convergence.R applies when reading a saved fit back later;
  # computed here too so it's visible immediately in this run's own console
  # output and diagnostics row, not just on a later separate pass.
  converged <- is.finite(max_rhat) && is.finite(min_ess) &&
    max_rhat < 1.01 && min_ess > 400

  cat("    Max Rhat:", round(max_rhat, 4),
      " | Min ESS:", round(min_ess, 0),
      " | Converged (Rhat<1.01 & ESS>400):", converged,
      if (!is.na(gamma1_mean)) paste(" | gamma1:", round(gamma1_mean, 4)) else "",
      "\n")

  data.frame(
    species     = species,
    model       = model_tag,
    max_rhat    = round(max_rhat, 4),
    min_ess     = round(min_ess, 0),
    converged   = converged,
    gamma1_mean = round(gamma1_mean, 4),
    n_obs       = nrow(new_data),
    n_dropped_missing_cov = prepped$n_dropped,
    fit_min     = fit_time
  )
}
