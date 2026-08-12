## =============================================================================
## Full production run: iCAR route-level slope model, base vs. anthro
## Time period: 2010-2025
##
## Fits TWO models for EVERY BBS species in this project's species list
## (data/spp_names_codes_group_aou.csv, in_bbs == TRUE -- all groups:
## aridlands, boreal_forests, eastern_forests, western_forests,
## subtropical_forests, grasslands, marshlands, coastal, waterbirds,
## generalists, arctic, urban_suburban -- currently 607 species):
##
##   base   -- no covariate    (models/slope_iCAR_route_NB_New.stan)
##   anthro -- data/Anthro.csv  (models/slope_iCAR_route_NB_New_covariate.stan)
##     log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*Anthro[r,t]
##
## This replaces the earlier per-group, many-covariates design (random
## N-species panels x several habitat/anthro/transition covariates per
## forest group). That pilot work (see output/files/aridlands_2010_2025.txt
## and the gamma_lookup/beta_compare write-ups) found: (1) the "*_to_anthro"
## transition covariates never produced a credible gamma1 in any tested
## species -- too rare/sparse to be identifiable at this route/year
## resolution; (2) "anthro" is the one covariate that consistently produced
## a real, credible, well-powered effect; and (3) group-specific habitat
## covariates aren't available/consistent across every one of the 12 groups
## the way Anthro.csv is (it's already used project-wide). So this version
## drops the per-group covariate machinery entirely and just runs anthro,
## for everyone, everywhere.
##
## Both models are fit on the SAME reduced dataset per species (only
## route-years where Anthro is non-NA), same "same data for a fair
## comparison" principle used throughout this project.
##
## Covariate data: data/Anthro.csv, one row per BBS route per year
## (CountryNum/StateNum/Route/.../year/Anthro), loaded with load_covariate()
## the same way every other covariate in this project has been. Not
## rescaled (standing-stock proportion already spans most of [0,1]).
##
## Route-key format: new_data$route (from bbsBayes2::prepare_data()) is
## formatted "<StateNum>-<Route>" -- but ONLY when stratified by "bcr" (see
## every earlier version of this script for the full explanation). strat =
## "bcr" is used here for the same reason.
##
## Per user decision (carried over from every earlier version of this
## script): counts with a missing (NA) Anthro value are dropped from that
## species' stan_data, applied once up front, so both model tags see
## identical data.
##
## Species: EVERY species in the project's list, not a random panel --
## sorted alphabetically for a stable, resumable run order.
##
## force_refit defaults to FALSE here (unlike the earlier pilot versions of
## this script, which defaulted to TRUE) -- this is a full ~600-species
## production run, not a small re-run-every-time pilot, so by default this
## script SKIPS any species/tag that already has a saved summ_fit.rds and
## simply resumes where a previous, possibly-interrupted run left off. Set
## force_refit <- TRUE to force a clean re-fit of everything instead.
##
## This is a LONG run: ~600 species x 2 models x (BBS data pull + Voronoi
## neighbours + 4-chain/4000-iteration Stan fit) each. Expect this to take
## many hours to days depending on hardware -- that's expected, not a bug.
##
## This script fits the models and saves, per species and per model tag:
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_stanfit.rds
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds
##   data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData
## Plus one combined diagnostics CSV and (via helper/gamma_lookup.R) one
## gamma1 lookup table across all species.
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
# don't pile up in %TEMP% across a run of ~600 species x 2 model tags.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ----------------------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025
dt         <- lastYear - firstYear

strat <- "bcr"   # NOT "bbs_usgs" -- see earlier versions of this script for
                 # the full explanation.

# Re-fit control: FALSE by default for this full-species run -- see header
# comment. Set TRUE to force a clean re-fit of everything.
force_refit <- FALSE

model_tags <- c("base", "anthro")

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
# it appears in the CSV header. out_name becomes both the joined column name
# in new_data AND the model_tag used to select it in
# fit_one_covariate_model() (see below).
load_covariate <- function(file, value_col, out_name) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              !!out_name := .data[[value_col]]) %>%
    filter(!is.na(.data[[out_name]])) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

# Load one covariate per the spec above, optionally rescaling by its own
# observed max. Kept as a thin wrapper (rather than calling load_covariate()
# directly) so this script stays easy to extend with another covariate
# later without restructuring the loop below.
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

# Covariate spec -- "anthro" doubles as both the model_tag AND the joined
# column name in new_data (see prepare_species_data()/
# fit_one_covariate_model() below, which key off this directly instead of a
# per-tag if/else chain). Anthro.csv is used project-wide, independent of
# bird group, so there is no per-group covariate config here anymore.
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# Species list -- EVERY species, not a random panel -----------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Full run: ", nrow(target_spp), " species x ", length(model_tags),
    " models (", paste(model_tags, collapse = ", "), ") ===\n", sep = "")
cat("Groups represented:", paste(sort(unique(target_spp$Group)), collapse = ", "), "\n")

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Compile the Stan models once (reused across every species) ---------------
model_base   <- cmdstan_model("models/slope_iCAR_route_NB_New.stan",
                              stanc_options = list("O1"))
model_single <- cmdstan_model("models/slope_iCAR_route_NB_New_covariate.stan",
                              stanc_options = list("O1"))

# Load the (single) covariate once, up front, for every species ------------
cat("Loading covariates...\n")
covariate_lookups <- list()
for (cn in setdiff(model_tags, "base")) {
  spec <- covariate_specs[[cn]]
  covariate_lookups[[cn]] <- load_group_covariate(spec$file, spec$value_col, cn, spec$rescale)
}

# ==========================================================================
# Function: data-prep pipeline for one species (BBS data, covariate join,
# NA-drop, spatial neighbours). Generalized to an arbitrary NAMED LIST of
# covariate lookups (currently just "anthro") instead of a fixed signature,
# so this still supports adding another covariate back later without
# restructuring. Shared across both model tags for a species so they fit on
# identical data.
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

  # Drop rows missing anthro -- same reduced-dataset principle used
  # throughout this project, so both model tags (including "base") fit on
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

# ==========================================================================
# Function: fit one model_tag ("base" or "anthro") for one species, given
# the shared prepared data from prepare_species_data(). Returns a one-row
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
           "check covariate_specs.")
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
# Main loop: every species, both model tags. Flat -- no per-group looping
# anymore, since anthro applies the same way to every species regardless of
# its Group.
# ==========================================================================
gamma_lookup_skip_autorun <- TRUE
source(here::here("helper", "gamma_lookup.R"))

results_list <- list()
diagnostics_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp       <- target_spp$Common.Name[i]
  sp_f     <- species_to_f(sp)
  sp_code  <- target_spp$Code[i]
  sp_bbs   <- target_spp$bbs_english[i]
  sp_group <- target_spp$Group[i]

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(target_spp), "]", sp, " (", sp_group, ")\n")
  cat("================================================================\n")

  all_exist <- all(vapply(model_tags, function(tag) {
    out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    file.exists(summ_file)
  }, logical(1)))

  if (all_exist && !force_refit) {
    cat("  Skipping (both", paste(model_tags, collapse = " and "), "already fitted for", sp, ")\n")
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

    diagnostics$group <- sp_group
    results_list[[paste(sp, tag, sep = " | ")]] <- summ_file
    diagnostics_list[[paste(sp, tag, sep = " | ")]] <- diagnostics

    cat("  Done —", tag, "fit + diagnostics saved for", sp, "\n")
  }

  # Write/refresh the combined diagnostics CSV every 25 species, not just at
  # the very end -- with a ~600-species, multi-hour/day run, losing all
  # progress info to an interruption right before the final write would be
  # painful. Cheap to do (append-in-memory, rewrite whole file).
  if (length(diagnostics_list) > 0 && i %% 25 == 0) {
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_covariates_all_species_anthro_",
                                 firstYear, "_", lastYear, ".csv"))
    write.csv(bind_rows(diagnostics_list), diag_csv, row.names = FALSE)
    cat("  [checkpoint] diagnostics written to:", diag_csv, "(", i, "/", nrow(target_spp), "species processed )\n")
  }
}

# Final combined diagnostics CSV ---------------------------------------------
if (length(diagnostics_list) > 0) {
  diagnostics_all <- bind_rows(diagnostics_list)
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_covariates_all_species_anthro_",
                               firstYear, "_", lastYear, ".csv"))
  write.csv(diagnostics_all, diag_csv, row.names = FALSE)
  cat("\nDiagnostics written to:", diag_csv, "\n")
}

# gamma1 lookup table across every species -----------------------------------
gamma_lookup_table(target_spp = target_spp, model_tags = model_tags,
                   firstYear = firstYear, lastYear = lastYear,
                   rds_dir = rds_dir, bird_group = "all_species_anthro")

cat("\n\n=== DONE ===\n")
cat("Species x model fits this run:", length(results_list), "\n")
cat("Fits saved in:", output_dir, "\n")
