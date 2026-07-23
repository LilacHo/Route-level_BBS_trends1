## =============================================================================
## Grassland birds: iCAR route-level slope model + land-cover covariates
## Time period: 2010-2025
##
## Fits FOUR models per species on the SAME reduced dataset (only route-years
## where grassland_habitat, grassland_anthro, AND grassland_habitat_to_anthro
## are all non-NA — i.e. routes without full covariate coverage are dropped
## entirely, not just the missing observations), so all four are directly
## comparable to each other:
##
##   base                      models/slope_iCAR_route_NB_New.stan
##     log(lambda[r,t]) = alpha[r] + beta[r]*t
##   grassland_habitat          models/slope_iCAR_route_NB_New_covariate.stan
##     log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*GrasslandHabitat[r,t]
##   grassland_anthro           models/slope_iCAR_route_NB_New_covariate.stan
##     log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*GrasslandAnthro[r,t]
##   grassland_habitat_to_anthro  models/slope_iCAR_route_NB_New_covariate.stan
##     log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*GrasslandHabitatToAnthro[r,t]
##
## NOTE: "base" here is a SEPARATE fit from 1_species_iCAR_2010_2025.R's own
## base fit. That script fits on the FULL dataset (all routes); this script's
## "base" tag refits the same no-covariate model but restricted to the
## covariate-reduced dataset, specifically so it's comparable to the other
## three models on identical data. Both fits are kept (different output file
## names), since they answer different questions.
##
## Covariate data (each one row per BBS route per year):
##   data/grassland habitat.csv          — "grassland habitat" column.
##     Confirmed equal to LndCov_71_GrasslandHerbaceous + LndCov_81_PastureHay
##     (natural grassland + hayed/pastured grassland), range [0, 0.99].
##   data/grassland anthro.csv           — "grassland anthro" column.
##     Confirmed equal to the sum of the four Developed* NLCD classes (open
##     space, low/medium/high intensity) plus CultivatedCrops, range [0, 1.0].
##     This REPLACES the earlier data/Anthro.csv as the "anthro" covariate —
##     Anthro.csv's definition couldn't be reconstructed from the NLCD
##     land-cover components and is no longer used here.
##   data/grassland habitat to anthro.csv — "grassland habitat to anthro"
##     column, a year-over-year grassland-habitat-to-anthro transition
##     metric, range [0, ~0.15], mean ~0.00095 (population sd/max ratio
##     ~0.015 — even more zero-inflated/narrow than the earlier grasslands-
##     to-Anthro metric it replaces). Rescaled by its own observed max below
##     for the same reason grass_to_anthro was rescaled previously (see
##     fit_one_covariate_model()/ghta-style comment there): raw values only
##     span a sliver of [0,1], far narrower than gamma1 ~ normal(0,1)
##     implicitly assumes, and standardizing (rather than max-rescaling) was
##     shown empirically to collapse this class of covariate's effect toward
##     zero, so max-rescaling is used again here.
##
## grassland_habitat and grassland_anthro are NOT rescaled — both already
## span a wide, non-degenerate range within [0,1] (sd/max ~0.20 and ~0.27
## respectively), the same situation as the earlier grassland/anthro
## covariates, which were also used raw.
##
## Route-key format: new_data$route (from bbsBayes2::prepare_data()) is
## formatted "<StateNum>-<Route>" — but ONLY when stratified by "bcr"; the
## "bbs_usgs" stratification renames routes differently and gave a ~38%
## match rate against the covariate CSVs. Confirmed against
## 1_baseline_NB_fit.R, a working reference script that joins these same
## covariate CSVs successfully using strat = "bcr" plus a continental-US
## filter (country_num == 840, state_num != 3). The join below still prints
## sample values from both sides and checks its own match rate as a
## safety net.
##
## Per user decision (2026-07-13 conversation):
##   - Counts with a missing (NA) covariate value are dropped from that
##     species' stan_data (reduces ncounts) — applied once, up front, so all
##     four models (including "base") see identical data.
##   - The 2010 gap in the covariate CSVs (they start at 2011) is NOT
##     specially handled here — 2010 counts simply end up with NA covariates
##     and are dropped by the rule above. No back-filling or year-shifting.
##
## This script fits the models and saves, per species and per model tag
## ("base", "grassland_habitat", "grassland_anthro",
## "grassland_habitat_to_anthro"):
##   output/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds
##   output/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData
##
## No leave-future-out CV here: models/slope_iCAR_route_NB_cv.stan (used for
## CV in 1_species_iCAR_2010_2025.R) was not extended with covariate terms,
## so CV is out of scope for this script.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("1c_species_iCAR_covariates.R")
source("functions/neighbours_define_voronoi.R")
source("functions/posterior_summary_functions.R")

# Ensure cmdstanr writes output to CSV (default) and set output_dir for temp
# files. This avoids the "No CmdStan config files found" error by using
# output_dir so cmdstanr can find its files.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ----------------------------------------------------------------
land_cover <- "grasslands"   # target group: must match a value in the
                             # 'Group' column of spp_names_codes_group_aou.csv
firstYear  <- 2010
lastYear   <- 2025
dt         <- lastYear - firstYear

strat <- "bcr"   # NOT "bbs_usgs" — bbsBayes2 renames/reformats route IDs
                 # differently per stratification scheme, and "bcr" is what
                 # produces the "<StateNum>-<Route>" route IDs the covariate
                 # CSVs are keyed by (confirmed by 1_baseline_NB_fit.R, a
                 # working reference that joins these same covariate files
                 # successfully using this stratification).

# Re-fit control: if TRUE, ignore existing fits and re-run the model (needed
# after changing the model, priors, or data thresholds).
force_refit <- TRUE

# Output directories
output_dir <- here::here("output")
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(here::here("data"))) dir.create(here::here("data"))
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)

# Land-cover covariate lookups ---------------------------------------------
# One row per route per year; build a "<StateNum>-<Route>" key to match
# new_data$route below, and keep just (route_key, year, value).
# check.names = FALSE preserves value_col exactly as it appears in the CSV
# header (needed since all three new covariate names contain spaces).
load_covariate <- function(file, value_col) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              value     = .data[[value_col]]) %>%
    filter(!is.na(value)) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

grass_habitat_lookup <- load_covariate("grassland habitat.csv", "grassland habitat") %>%
  rename(grassland_habitat = value)
grass_anthro_lookup  <- load_covariate("grassland anthro.csv", "grassland anthro") %>%
  rename(grassland_anthro = value)
ghta_lookup           <- load_covariate("grassland habitat to anthro.csv", "grassland habitat to anthro") %>%
  rename(grassland_habitat_to_anthro = value)

# grassland habitat to anthro.csv is a year-over-year TRANSITION rate
# (fraction of a route's area converting from grassland habitat to
# anthropogenic cover that year), not a standing-stock proportion like
# grassland_habitat/grassland_anthro — most routes see zero or near-zero
# conversion in a given year, so raw values only span [0, ~0.15] with a mean
# near 0.00095. gamma1 ~ normal(0,1) in
# slope_iCAR_route_NB_New_covariate.stan implicitly assumes a covariate that
# swings across roughly [0,1], so left un-rescaled this covariate could
# barely move the fit regardless of the true relationship. Rescaling by its
# own observed max (a pure linear rescale — order/ratios between routes are
# unchanged, zeros stay zero) stretches it back out to [0,1], matching the
# scale grassland_habitat/grassland_anthro already use and keeping the same
# bounded single-covariate Stan model/prior valid for it.
ghta_scale <- max(ghta_lookup$grassland_habitat_to_anthro, na.rm = TRUE)
ghta_lookup <- ghta_lookup %>% mutate(grassland_habitat_to_anthro = grassland_habitat_to_anthro / ghta_scale)
cat("Grassland-habitat-to-Anthro covariate rescaled by its observed max (",
    round(ghta_scale, 6), ") to span [0,1]\n", sep = "")

cat("Grassland habitat covariate:           ", nrow(grass_habitat_lookup), " route-year rows (",
    length(unique(grass_habitat_lookup$route_key)), " routes)\n", sep = "")
cat("Grassland anthro covariate:            ", nrow(grass_anthro_lookup), " route-year rows (",
    length(unique(grass_anthro_lookup$route_key)), " routes)\n", sep = "")
cat("Grassland habitat to anthro covariate:  ", nrow(ghta_lookup), " route-year rows (",
    length(unique(ghta_lookup$route_key)), " routes)\n", sep = "")

# Target group species list -----------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

# Optional external override — set `species_filter <- c("Some Species", ...)`
# in the calling environment before sourcing this script to restrict the run
# to just those species (used by tests/test_single_species.R).
if (exists("species_filter")) {
  target_spp <- target_spp %>% filter(Common.Name %in% species_filter)
  cat("species_filter applied — running for:", paste(species_filter, collapse = ", "), "\n")
}

cat("=== iCAR route-level slope model + land-cover covariates ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat(land_cover, "species to fit (n =", nrow(target_spp), "):\n")
print(target_spp %>% select(Common.Name, Code, bbs_english))

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Compile the Stan models once (reused across species) ---------------------
model_base   <- cmdstan_model("models/slope_iCAR_route_NB_New.stan",
                              stanc_options = list("O1"))
model_single <- cmdstan_model("models/slope_iCAR_route_NB_New_covariate.stan",
                              stanc_options = list("O1"))

# Four model tags fit per species, on the same prepared (reduced) dataset --
model_tags <- c("base", "grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro")

# ==========================================================================
# Function: data-prep pipeline for one species (BBS data, covariate join,
# NA-drop, spatial neighbours). Shared across all four model tags so they
# fit on identical data. Returns a list with new_data, route_map,
# realized_strata_map, car_stan_dat, sd_alpha_prior, base_stan_data (the
# model-agnostic Stan data fields), and n_dropped.
# ==========================================================================
prepare_species_data <- function(species, species_bbs, strat,
                                 firstYear, lastYear,
                                 grass_habitat_lookup, grass_anthro_lookup, ghta_lookup) {

  cat("    Preparing data for:", species, "\n")

  # BBS data ---------------------------------------------------------------
  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data <- data_pkg[["raw_data"]]
  strata_map <- load_map(strat)

  # Continental US only (matches 1_baseline_NB_fit.R) — the land-cover
  # covariate rasters/CSVs only cover the continental US, and state_num == 3
  # is excluded there too.
  raw_data <- raw_data %>%
    filter(country_num == 840) %>%
    filter(state_num != 3)

  # Use all available routes (no further spatial route filtering) ----------
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

  # Join land-cover covariates by route + calendar year ---------------------
  n_before <- nrow(new_data)
  new_data <- new_data %>%
    left_join(grass_habitat_lookup, by = c("route" = "route_key", "r_year" = "year")) %>%
    left_join(grass_anthro_lookup,  by = c("route" = "route_key", "r_year" = "year")) %>%
    left_join(ghta_lookup,          by = c("route" = "route_key", "r_year" = "year"))

  match_rate <- mean(!is.na(new_data$grassland_habitat) & !is.na(new_data$grassland_anthro) &
                      !is.na(new_data$grassland_habitat_to_anthro))
  cat("    Covariate match rate (grassland_habitat, grassland_anthro & grassland_habitat_to_anthro all present):",
      round(100 * match_rate, 1), "%\n")

  # Diagnostic: always print a few sample route IDs from each side of the
  # join, so a format mismatch is visible even when match_rate is nonzero
  # (e.g. only some states/strata happen to line up under the current
  # "<StateNum>-<Route>" assumption).
  cat("    Sample new_data$route values:      ",
      paste(head(unique(new_data$route), 6), collapse = ", "), "\n")
  cat("    Sample covariate route_key values: ",
      paste(head(unique(grass_habitat_lookup$route_key), 6), collapse = ", "), "\n")

  if (match_rate < 0.5) {
    stop("Less than half of observations matched a grassland_habitat/grassland_anthro/",
         "grassland_habitat_to_anthro value for ", species, " — new_data$route doesn't ",
         "look like '<StateNum>-<Route>'. Compare the sample values printed ",
         "above. This previously happened because strat was 'bbs_usgs' ",
         "instead of 'bcr' (bbsBayes2 renames routes differently per ",
         "stratification); if strat is already 'bcr', check the ",
         "continental-US filter and the key construction in ",
         "load_covariate() instead.")
  }

  # Per user decision: drop counts with a missing covariate (rather than
  # impute or special-case the pre-2011 gap in the covariate CSVs). Dropped
  # once here, up front, so grassland_habitat-only, grassland_anthro-only,
  # grassland_habitat_to_anthro-only, and base all fit on the same rows.
  new_data <- new_data %>% filter(!is.na(grassland_habitat), !is.na(grassland_anthro),
                                   !is.na(grassland_habitat_to_anthro))
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

  # Model-agnostic Stan data fields (covariate fields added per model tag) --
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

# ==========================================================================
# Function: fit one of the four models (base, grassland_habitat,
# grassland_anthro, grassland_habitat_to_anthro) for one species, given the
# shared prepared data from prepare_species_data(). Returns a one-row
# diagnostics data.frame (or NULL on failure, handled by the caller).
# ==========================================================================
fit_one_covariate_model <- function(species, species_f, model_tag, prepped,
                                    firstYear, lastYear) {

  new_data       <- prepped$new_data
  base_stan_data <- prepped$base_stan_data

  # Covariate field + Stan model for this tag --------------------------------
  if (model_tag == "base") {
    stan_model <- model_base
    stan_data  <- base_stan_data
  } else if (model_tag == "grassland_habitat") {
    stan_model <- model_single
    stan_data  <- c(base_stan_data, list(covariate = new_data$grassland_habitat))
  } else if (model_tag == "grassland_anthro") {
    stan_model <- model_single
    stan_data  <- c(base_stan_data, list(covariate = new_data$grassland_anthro))
  } else if (model_tag == "grassland_habitat_to_anthro") {
    stan_model <- model_single
    stan_data  <- c(base_stan_data, list(covariate = new_data$grassland_habitat_to_anthro))
  } else {
    stop("Unknown model_tag: ", model_tag)
  }

  cat("    Fitting model:", model_tag, "\n")

  stanfit <- stan_model$sample(
    data = stan_data,
    refresh = 400,
    chains = 4, iter_sampling = 2000, iter_warmup = 2000,
    parallel_chains = 4,
    adapt_delta = 0.8,
    max_treedepth = 10,
    show_exceptions = FALSE,
    output_dir = cmdstanr_output_dir
  )

  summ <- stanfit$summary()
  fit_time <- round(stanfit$time()[["total"]] / 60, 1)
  cat("    Fit time:", fit_time, "minutes\n")

  out_base <- paste0(species_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear)
  stanfit$save_object(file.path(output_dir, paste0(out_base, "_stanfit.rds")))
  saveRDS(summ, file.path(output_dir, paste0(out_base, "_summ_fit.rds")))

  sp_data_file <- here::here("data", "stan_data",
                             paste0(species_f, "_", model_tag, "_",
                                    firstYear, "_", lastYear, "_stan_data.RData"))
  save(list = c("stan_data", "new_data", "model_tag"), file = sp_data_file)

  # Convergence diagnostics + covariate effect -------------------------------
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess  <- min(summ$ess_bulk, na.rm = TRUE)

  # "base" has no gamma1 parameter at all; the three covariate models each
  # have exactly one (gamma1). Extract defensively so a missing variable
  # (rather than an unexpected one) never errors here.
  gamma1_mean <- if ("gamma1" %in% summ$variable) {
    summ$mean[summ$variable == "gamma1"]
  } else {
    NA_real_
  }

  cat("    Max Rhat:", round(max_rhat, 4),
      " | Min ESS:", round(min_ess, 0),
      if (!is.na(gamma1_mean)) paste(" | gamma1:", round(gamma1_mean, 4)) else "",
      "\n")

  data.frame(
    species     = species,
    model       = model_tag,
    max_rhat    = round(max_rhat, 4),
    min_ess     = round(min_ess, 0),
    gamma1_mean = round(gamma1_mean, 4),
    n_obs       = nrow(new_data),
    n_dropped_missing_cov = prepped$n_dropped,
    fit_min     = fit_time
  )
}

# ==========================================================================
# Main loop: for each species, prepare data once, then fit all four models
# ==========================================================================
results_list <- list()
diagnostics_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]
  sp_bbs  <- target_spp$bbs_english[i]   # name used to look up species in BBS

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(target_spp), "]", sp, "\n")
  cat("================================================================\n")

  # Skip data prep entirely if all four model outputs already exist
  all_exist <- all(vapply(model_tags, function(tag) {
    out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
    file.exists(summ_file)
  }, logical(1)))

  if (all_exist && !force_refit) {
    cat("  Skipping (all", length(model_tags), "models already fitted for", sp, ")\n")
    next
  }

  # Prepare data once (shared across all four model fits) --------------------
  prepped <- tryCatch(
    prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                         firstYear = firstYear, lastYear = lastYear,
                         grass_habitat_lookup = grass_habitat_lookup,
                         grass_anthro_lookup = grass_anthro_lookup,
                         ghta_lookup = ghta_lookup),
    error = function(e) {
      message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(prepped)) next

  # Fit each of the four models -----------------------------------------------
  for (tag in model_tags) {
    out_base       <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file      <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
    stan_data_file <- here::here("data", "stan_data",
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_stan_data.RData"))
    if (file.exists(summ_file) && file.exists(stan_data_file) && !force_refit) {
      cat("  Skipping (already fitted):", basename(summ_file), "\n")
      next
    }

    diagnostics <- tryCatch(
      fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                              prepped = prepped, firstYear = firstYear,
                              lastYear = lastYear),
      error = function(e) {
        message("  [ERROR] Model '", tag, "' failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(diagnostics)) next

    results_list[[paste(sp, tag, sep = " | ")]] <- summ_file
    diagnostics_list[[paste(sp, tag, sep = " | ")]] <- diagnostics

    cat("  Done —", tag, "fit + diagnostics saved for", sp, "\n")
  }
}

# Combined diagnostics CSV -------------------------------------------------
if (length(diagnostics_list) > 0) {
  diagnostics_all <- bind_rows(diagnostics_list)
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_covariates_", land_cover, "_",
                               firstYear, "_", lastYear, ".csv"))
  write.csv(diagnostics_all, diag_csv, row.names = FALSE)
  cat("\nDiagnostics written to:", diag_csv, "\n")
  print(diagnostics_all)
}

cat("\n=== Summary ===\n")
cat("Species x model fits this run:", length(results_list), "\n")
cat("Fits saved in:", output_dir, "\n")
