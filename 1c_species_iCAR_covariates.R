## =============================================================================
## Forest birds: iCAR route-level slope model + land-cover covariates
## Time period: 2010-2025
##
## Fits ONE no-covariate model ("base") plus several single-covariate models
## per species, for a RANDOM n_species_per_group-species panel drawn from
## each of four forest bird groups (boreal_forests, eastern_forests,
## subtropical_forests, western_forests). The group being fit determines
## which covariates are used:
##
##   boreal_forests / western_forests (SIX models each):
##     base, forests, evergreen_and_mixed_forests, anthro,
##     forests_to_anthro, evergreen_and_mixed_forests_to_anthro
##   eastern_forests (SIX models):
##     base, forests, eastern_forests, anthro,
##     forests_to_anthro, eastern_forests_to_anthro
##   subtropical_forests (FOUR models):
##     base, forests, anthro, forests_to_anthro
##     (the user's request named only forests/anthro/forests_to_anthro but
##     said "four models" -- "base" is assumed to be the implicit 4th, to
##     match the pattern used by every other group here. Flag if that's
##     wrong.)
##
## All non-base models share the same single-covariate Stan model
## (models/slope_iCAR_route_NB_New_covariate.stan):
##   log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*Covariate[r,t]
## "base" uses models/slope_iCAR_route_NB_New.stan (no covariate term).
##
## For each group, all of that group's models are fit on the SAME reduced
## dataset (only route-years where every one of that group's covariates is
## non-NA), so they're directly comparable to each other -- same principle
## used by the grasslands/marshlands versions of this script.
##
## Covariate data -- these CSVs now exist in data/ (dropped in after this
## script was first drafted), in the same "one row per BBS route per year"
## shape as every other covariate CSV in this project
## (CountryNum/StateNum/Route/.../year/<value column>), loaded with
## load_covariate() exactly like grassland habitat.csv / wetlands.csv /
## Anthro.csv were:
##   data/forests.csv                                     -- "forests" column
##   data/evergreen_and_mixed_forests.csv                  -- "evergreen_and_mixed_forests" column
##   data/eastern_forests.csv                              -- "eastern_forests" column
##   data/Anthro.csv                                       -- "Anthro" column (already used by the
##                                                             marshlands version of this script)
##   data/forests to Anthro.csv                            -- "forests to Anthro" column
##   data/evergreen_and_mixed_forests to Anthro.csv         -- "evergreen_and_mixed_forests to Anthro" column
##   data/eastern_forests to Anthro.csv                     -- "eastern_forests to Anthro" column
## Naming convention actually used by the files as provided: the multi-word
## habitat covariates ("evergreen_and_mixed_forests", "eastern_forests")
## use UNDERSCORES in both the filename and column header, matching their
## model_tag names exactly -- but the " to Anthro" suffix on the transition
## covariates keeps a literal SPACE (e.g. "eastern_forests to Anthro.csv"),
## same as grassland habitat to anthro.csv / wetlands to Anthro.csv before
## it. NOTE the tag "eastern_forests" (a covariate, for the eastern_forests
## group only) is spelled identically to the bird group land_cover value
## "eastern_forests" -- they mean different things (one is a Group value,
## the other is a model_tag/covariate name) but share a name by coincidence;
## watch for this if grepping the codebase later.
##
## Rescaling: the three "*_to_anthro" transition covariates are expected to
## be zero-inflated/narrow like grassland_habitat_to_anthro and
## wetlands_to_anthro were (most routes near-zero conversion in a given
## year), so they're rescaled by their own observed max the same way, right
## after loading. "forests", "evergreen_and_mixed_forests", "eastern_forests"
## and "anthro" are standing-stock proportions and are NOT rescaled, matching
## grassland_habitat/grassland_anthro/wetlands/anthro previously.
##
## Route-key format: new_data$route (from bbsBayes2::prepare_data()) is
## formatted "<StateNum>-<Route>" -- but ONLY when stratified by "bcr"; see
## the grasslands/marshlands versions of this script for the full
## explanation. strat = "bcr" is used here for the same reason.
##
## Per user decision (carried over from the grasslands/marshlands versions):
##   - Counts with a missing (NA) covariate value (for ANY of that group's
##     non-base covariates) are dropped from that species' stan_data, applied
##     once up front, so every model tag for a group (including "base") sees
##     identical data.
##
## Random species selection: n_species_per_group species are drawn at random
## (not hand-picked, not alphabetical) from each group's in_bbs == TRUE pool
## via a fixed random seed (RANDOM_SEED below), so the panel is reproducible
## across re-runs of this script without needing to hard-code species names.
##
## This script fits the models and saves, per species and per model tag:
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds
##   data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData
## Plus, per group, a combined diagnostics CSV and (via
## helper/gamma_lookup.R) a gamma1 lookup table -- both written once per
## group so the four groups' outputs don't overwrite each other.
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
# files. Each species/model_tag combo gets its own subfolder under here
# (created/deleted inside fit_one_covariate_model()) so raw per-chain CSVs
# don't pile up in %TEMP% across a whole run of many species x model tags.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ----------------------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025
dt         <- lastYear - firstYear

strat <- "bcr"   # NOT "bbs_usgs" -- see grasslands/marshlands versions of
                 # this script for the full explanation.

# Re-fit control: if TRUE, ignore existing fits and re-run the model.
force_refit <- TRUE

forest_groups       <- c("boreal_forests", "eastern_forests",
                         "subtropical_forests", "western_forests")
n_species_per_group <- 10
RANDOM_SEED          <- 2026   # fixed so the random species panels below are
                               # reproducible across re-runs

# Output directories --------------------------------------------------------
output_dir <- here::here("output")
if (!dir.exists(output_dir)) dir.create(output_dir)
rds_dir <- here::here("output", "rds")
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
if (!dir.exists(here::here("data"))) dir.create(here::here("data"))
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)
route_info_dir <- here::here("data", "route_info")
if (!dir.exists(route_info_dir)) dir.create(route_info_dir, recursive = TRUE)

# Covariate loader ----------------------------------------------------------
# One row per route per year; build a "<StateNum>-<Route>" key to match
# new_data$route below. check.names = FALSE preserves value_col exactly as
# it appears in the CSV header (several of these covariate names contain
# spaces). out_name becomes both the joined column name in new_data AND the
# model_tag used to select it in fit_one_covariate_model() (see below).
load_covariate <- function(file, value_col, out_name) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              !!out_name := .data[[value_col]]) %>%
    filter(!is.na(.data[[out_name]])) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

# Load one covariate per the spec above, optionally rescaling by its own
# observed max (for the zero-inflated "*_to_anthro" transition covariates --
# see header comment).
load_group_covariate <- function(file, value_col, out_name, rescale) {
  lk <- load_covariate(file, value_col, out_name)
  if (rescale) {
    cscale <- max(lk[[out_name]], na.rm = TRUE)
    lk[[out_name]] <- lk[[out_name]] / cscale
    cat("  ", out_name, "covariate rescaled by its observed max (", round(cscale, 6), ") to span [0,1]\n")
  }
  cat("  ", out_name, "covariate:", nrow(lk), "route-year rows (",
      length(unique(lk$route_key)), "routes)\n")
  lk
}

# Per-group model tags + covariate specs -------------------------------------
# "forests"/"anthro" style tags below double as both the model_tag AND the
# joined column name in new_data (see prepare_species_data()/
# fit_one_covariate_model() below, which key off this directly instead of a
# per-tag if/else chain).
conifer_group_covariates <- list(
  forests                              = list(file = "forests.csv", value_col = "forests", rescale = FALSE),
  evergreen_and_mixed_forests          = list(file = "evergreen_and_mixed_forests.csv", value_col = "evergreen_and_mixed_forests", rescale = FALSE),
  anthro                               = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE),
  forests_to_anthro                    = list(file = "forests to Anthro.csv", value_col = "forests to Anthro", rescale = TRUE),
  evergreen_and_mixed_forests_to_anthro = list(file = "evergreen_and_mixed_forests to Anthro.csv", value_col = "evergreen_and_mixed_forests to Anthro", rescale = TRUE)
)

group_config <- list(
  boreal_forests = list(
    model_tags = c("base", "forests", "evergreen_and_mixed_forests", "anthro",
                   "forests_to_anthro", "evergreen_and_mixed_forests_to_anthro"),
    covariates = conifer_group_covariates
  ),
  western_forests = list(
    model_tags = c("base", "forests", "evergreen_and_mixed_forests", "anthro",
                   "forests_to_anthro", "evergreen_and_mixed_forests_to_anthro"),
    covariates = conifer_group_covariates
  ),
  eastern_forests = list(
    model_tags = c("base", "forests", "eastern_forests", "anthro",
                   "forests_to_anthro", "eastern_forests_to_anthro"),
    covariates = list(
      forests                   = list(file = "forests.csv", value_col = "forests", rescale = FALSE),
      eastern_forests           = list(file = "eastern_forests.csv", value_col = "eastern_forests", rescale = FALSE),
      anthro                    = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE),
      forests_to_anthro         = list(file = "forests to Anthro.csv", value_col = "forests to Anthro", rescale = TRUE),
      eastern_forests_to_anthro = list(file = "eastern_forests to Anthro.csv", value_col = "eastern_forests to Anthro", rescale = TRUE)
    )
  ),
  subtropical_forests = list(
    # Request said "four models: forests, Anthro, forests to Anthro" (three
    # named) -- "base" added here as the assumed implicit 4th, matching
    # every other group's pattern (see header comment).
    model_tags = c("base", "forests", "anthro", "forests_to_anthro"),
    covariates = list(
      forests           = list(file = "forests.csv", value_col = "forests", rescale = FALSE),
      anthro            = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE),
      forests_to_anthro = list(file = "forests to Anthro.csv", value_col = "forests to Anthro", rescale = TRUE)
    )
  )
)

# Species list + random per-group panels -------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

set.seed(RANDOM_SEED)
species_by_group <- list()
for (g in forest_groups) {
  pool <- spp_df %>%
    filter(Group == g, in_bbs == TRUE) %>%
    distinct(Common.Name, Code, .keep_all = TRUE)
  species_by_group[[g]] <- pool %>%
    slice_sample(n = min(n_species_per_group, nrow(pool))) %>%
    arrange(Common.Name)
}

cat("=== Random", n_species_per_group, "-species panels per forest group (seed =", RANDOM_SEED, ") ===\n")
for (g in forest_groups) {
  cat(" ", g, "(", nrow(species_by_group[[g]]), "species ):",
      paste(species_by_group[[g]]$Common.Name, collapse = ", "), "\n")
}

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Compile the Stan models once (reused across every group/species) ---------
model_base   <- cmdstan_model("models/slope_iCAR_route_NB_New.stan",
                              stanc_options = list("O1"))
model_single <- cmdstan_model("models/slope_iCAR_route_NB_New_covariate.stan",
                              stanc_options = list("O1"))

# ==========================================================================
# Function: data-prep pipeline for one species (BBS data, covariate join,
# NA-drop, spatial neighbours). Generalized to an arbitrary NAMED LIST of
# covariate lookups (one per non-base model_tag for the current group)
# instead of a fixed 3-covariate signature, so the same function serves
# every group's covariate set. Shared across all of a group's model tags for
# a species so they fit on identical data.
# ==========================================================================
prepare_species_data <- function(species, species_bbs, strat,
                                 firstYear, lastYear, covariate_lookups) {

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

  # Join every covariate for this group by route + calendar year ------------
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

  # Drop rows missing ANY of this group's covariates -- same reduced-dataset
  # principle used throughout this project, so every model tag (including
  # "base") fits on identical rows.
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

# ==========================================================================
# Function: fit one model_tag for one species, given the shared prepared
# data from prepare_species_data(). Generalized: any non-"base" tag is
# assumed to name a column already joined onto new_data (see
# prepare_species_data() above, which joins each covariate under exactly its
# config name) -- no per-tag if/else chain needed. Returns a one-row
# diagnostics data.frame (or NULL on failure, handled by the caller).
# ==========================================================================
fit_one_covariate_model <- function(species, species_f, model_tag, prepped,
                                    firstYear, lastYear) {

  new_data       <- prepped$new_data
  base_stan_data <- prepped$base_stan_data

  if (model_tag == "base") {
    stan_model <- model_base
    stan_data  <- base_stan_data
  } else {
    if (!model_tag %in% names(new_data)) {
      stop("model_tag '", model_tag, "' has no matching covariate column in new_data — ",
           "check group_config's covariates list for this group.")
    }
    stan_model <- model_single
    stan_data  <- c(base_stan_data, list(covariate = new_data[[model_tag]]))
  }

  cat("    Fitting model:", model_tag, "\n")

  # Per-species/model_tag output subfolder so raw cmdstan CSVs stay isolated
  # and can be deleted right after this fit is saved to .rds.
  run_output_dir <- file.path(cmdstanr_output_dir, paste0(species_f, "_", model_tag))
  if (!dir.exists(run_output_dir)) dir.create(run_output_dir, recursive = TRUE)

  stanfit <- stan_model$sample(
    data = stan_data,
    refresh = 400,
    chains = 4, iter_sampling = 2000, iter_warmup = 2000,
    parallel_chains = 4,
    adapt_delta = 0.8,
    max_treedepth = 10,
    show_exceptions = FALSE,
    output_dir = run_output_dir
  )

  summ <- stanfit$summary()
  fit_time <- round(stanfit$time()[["total"]] / 60, 1)
  cat("    Fit time:", fit_time, "minutes\n")

  out_base <- paste0(species_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear)
  stanfit$save_object(file.path(rds_dir, paste0(out_base, "_stanfit.rds")))
  saveRDS(summ, file.path(rds_dir, paste0(out_base, "_summ_fit.rds")))

  # Lightweight route lookup (route, routeF, latitude, longitude) — saved on
  # its own, before cleanup below, so 2c can get what it needs without
  # depending on the full stan_data.RData save further down. Skip if a file
  # of the same name already exists.
  route_info_file <- file.path(route_info_dir,
                               paste0(species_f, "_", model_tag, "_",
                                      firstYear, "_", lastYear, "_route_info.rds"))
  if (!file.exists(route_info_file)) {
    route_info <- new_data %>% distinct(route, routeF, latitude, longitude)
    saveRDS(route_info, route_info_file)
  }

  # Raw per-chain CSVs are no longer needed once the fit is saved above —
  # delete them now instead of letting them accumulate in %TEMP%.
  unlink(run_output_dir, recursive = TRUE)

  sp_data_file <- here::here("data", "stan_data",
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
# Outer loop: one forest group at a time -- load that group's covariates,
# fit its species panel x model_tags, write its diagnostics + gamma1 lookup.
# ==========================================================================
# gamma_lookup_table() is soft-coded (see helper/gamma_lookup.R) -- source
# it once here (skipping its auto-run, since target_spp/model_tags/
# land_cover aren't fixed for the whole script — they change every group
# below) and call gamma_lookup_table() explicitly per group instead.
gamma_lookup_skip_autorun <- TRUE
source(here::here("helper", "gamma_lookup.R"))

results_list_all <- list()

for (land_cover in forest_groups) {

  cfg         <- group_config[[land_cover]]
  model_tags  <- cfg$model_tags
  target_spp  <- species_by_group[[land_cover]]

  cat("\n\n################################################################\n")
  cat("# FOREST GROUP:", land_cover, " | models:", paste(model_tags, collapse = ", "), "\n")
  cat("################################################################\n")

  # Load this group's covariate lookups (as-if the CSVs already exist --
  # see header comment) ------------------------------------------------------
  cat("Loading covariates for", land_cover, "...\n")
  covariate_lookups <- list()
  for (cn in setdiff(model_tags, "base")) {
    spec <- cfg$covariates[[cn]]
    covariate_lookups[[cn]] <- load_group_covariate(spec$file, spec$value_col, cn, spec$rescale)
  }

  results_list <- list()
  diagnostics_list <- list()

  for (i in seq_len(nrow(target_spp))) {
    sp      <- target_spp$Common.Name[i]
    sp_f    <- species_to_f(sp)
    sp_code <- target_spp$Code[i]
    sp_bbs  <- target_spp$bbs_english[i]

    cat("\n================================================================\n")
    cat("  [", land_cover, i, "/", nrow(target_spp), "]", sp, "\n")
    cat("================================================================\n")

    all_exist <- all(vapply(model_tags, function(tag) {
      out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
      summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
      file.exists(summ_file)
    }, logical(1)))

    if (all_exist && !force_refit) {
      cat("  Skipping (all", length(model_tags), "models already fitted for", sp, ")\n")
      next
    }

    prepped <- tryCatch(
      prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                           firstYear = firstYear, lastYear = lastYear,
                           covariate_lookups = covariate_lookups),
      error = function(e) {
        message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(prepped)) next

    for (tag in model_tags) {
      out_base       <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
      summ_file      <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
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

  # Combined diagnostics CSV for this group -----------------------------------
  if (length(diagnostics_list) > 0) {
    diagnostics_all <- bind_rows(diagnostics_list)
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_covariates_", land_cover, "_",
                                 firstYear, "_", lastYear, ".csv"))
    write.csv(diagnostics_all, diag_csv, row.names = FALSE)
    cat("\nDiagnostics for", land_cover, "written to:", diag_csv, "\n")
    print(diagnostics_all)
  }

  # gamma1 lookup table for this group -----------------------------------------
  gamma_lookup_table(target_spp = target_spp, model_tags = model_tags,
                     firstYear = firstYear, lastYear = lastYear,
                     rds_dir = rds_dir, land_cover = land_cover)

  cat("\n=== Summary for", land_cover, "===\n")
  cat("Species x model fits this run:", length(results_list), "\n")

  results_list_all[[land_cover]] <- results_list
}

cat("\n\n=== ALL FOREST GROUPS DONE ===\n")
cat("Total species x model fits this run:", sum(lengths(results_list_all)), "\n")
cat("Fits saved in:", output_dir, "\n")
