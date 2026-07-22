## =============================================================================
## Pilot test: iCAR route-level slope model + 16 simultaneous land-cover
## covariates (models/slope_iCAR_route_NB_New_16covariates.stan), across a
## 10-species panel of grassland birds.
##
## Unlike 1c_species_iCAR_covariates.R / tests/test_multi_species.R (which fit
## FIVE separate models built from 3 covariates: base, grassland, anthro,
## grassland_anthro, grass_to_anthro), this script fits ONE model per species
## that includes all 16 NLCD land-cover covariates simultaneously:
##
##   log(lambda[r,t]) = alpha[r] + beta[r]*t + sum_k gamma[k] * X[r,t,k]
##
## Covariate data: data/LndCov_<NLCD code>_<category>.csv, one file per land-
## cover class (16 files total: OpenWater, PerennialIceSnow, DevelopedOpenSpace,
## DevelopedLowIntensity, DevelopedMediumIntensity, DevelopedHighIntensity,
## BarrenLand, DeciduousForest, EvergreenForest, MixedForest, ShrubScrub,
## GrasslandHerbaceous, PastureHay, CultivatedCrops, WoodyWetlands,
## EmergentHerbaceousWetlands). Each file has one row per BBS route per year,
## columns CountryNum/StateNum/Route/.../year/<category>, with <category> a
## route-level proportion of land cover in [0,1] (confirmed by inspecting all
## 16 files: observed min = 0 and max < 1 in every file, no rescaling needed —
## unlike grass_to_anthro in the 5-model pipeline, these are standing-stock
## proportions already on the same [0,1] scale the gamma ~ normal(0,1) prior
## assumes). Column names have no spaces, so check.names = FALSE isn't
## strictly required here, but it's kept for consistency with load_covariate()
## in 1c_species_iCAR_covariates.R and in case a future covariate file does
## have a header with spaces.
##
## Route-key format: same "<StateNum>-<Route>" convention as the 5-model
## pipeline, which requires strat = "bcr" (see 1c_species_iCAR_covariates.R
## for the full explanation).
##
## Per-species data prep mirrors prepare_species_data() in
## 1c_species_iCAR_covariates.R: BBS pull -> continental-US filter -> join ALL
## 16 covariates by route + calendar year -> drop any row missing ANY of the
## 16 (shared reduced dataset, same "same data for a fair comparison"
## principle used throughout this project) -> Voronoi spatial neighbours.
##
## This is a PILOT: 10 species, one 16-covariate model each (real sampler
## settings — 4 chains x 2000 warmup + 2000 sampling — same as the production
## pipeline), not a fast unit test. Lower chains/iter below if a quick smoke
## test is needed instead.
##
## Outputs (per species), all under a dedicated subdirectory so they don't mix
## with the 5-model pipeline's output/ files:
##   output/covariates16/<species>_iCAR_16covariates_<firstYear>_<lastYear>_stanfit.rds
##   output/covariates16/<species>_iCAR_16covariates_<firstYear>_<lastYear>_summ_fit.rds
##   data/stan_data/<species>_16covariates_<firstYear>_<lastYear>_stan_data.RData
## Combined across all species:
##   output/covariates16/pilot_16covariates_gamma_summary.csv   (species x covariate)
##   output/covariates16/pilot_16covariates_diagnostics.csv     (species-level fit diagnostics)
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("tests/test_16covariates_pilot.R")
source(here::here("functions", "neighbours_define_voronoi.R"))

# Ensure cmdstanr writes output to CSV (default) and set output_dir for temp
# files, matching 1c_species_iCAR_covariates.R.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ------------------------------------------------------------------
firstYear <- 2010
lastYear  <- 2025

strat <- "bcr"   # NOT "bbs_usgs" — see 1c_species_iCAR_covariates.R for why.

# Re-fit control: if TRUE, ignore existing fits and re-run the model.
force_refit <- TRUE

# 10-species pilot panel: all previously confirmed working in the 5-model
# pipeline (Baird's Sparrow, Bobolink, Burrowing Owl, Chestnut-collared
# Longspur, Grasshopper Sparrow) plus 5 more well-sampled BBS grassland
# species, for a broader but still tractable pilot.
species_filter <- c(
  "Baird's Sparrow",
  "Bobolink",
  "Burrowing Owl",
  "Chestnut-collared Longspur",
  "Grasshopper Sparrow",
  "Dickcissel",
  "Eastern Meadowlark",
  "Western Meadowlark",
  "Vesper Sparrow",
  "Savannah Sparrow"
)

# Output directories ----------------------------------------------------------
output_dir16 <- here::here("output", "covariates16")
if (!dir.exists(output_dir16)) dir.create(output_dir16, recursive = TRUE)
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

## ============================================================================
## Discover the 16 LndCov_<code>_<category>.csv files in data/ and build the
## covariate lookup tables. Doing this by pattern (rather than hard-coding 16
## filenames) means a differently-numbered or re-run export still works as
## long as the LndCov_<code>_<category>.csv naming convention holds.
## ============================================================================
lndcov_files <- list.files(here::here("data"), pattern = "^LndCov_[0-9]+_.*\\.csv$")
if (length(lndcov_files) == 0) {
  stop("No data/LndCov_<code>_<category>.csv files found — check the data/ folder.")
}

cov_meta <- tibble(file = lndcov_files) %>%
  mutate(
    code     = str_match(file, "^LndCov_([0-9]+)_")[, 2],
    category = str_match(file, "^LndCov_[0-9]+_(.*)\\.csv$")[, 2]
  ) %>%
  arrange(as.integer(code))

K <- nrow(cov_meta)
cat("Discovered", K, "land-cover covariates in data/:\n")
print(cov_meta)

if (K != 16) {
  warning("Expected 16 LndCov_*.csv files but found ", K,
          " — proceeding anyway (K is used everywhere below, not hard-coded).")
}

# Generic covariate loader, mirrors load_covariate() in
# 1c_species_iCAR_covariates.R, but keeps the value column named after the
# land-cover category itself (rather than always "value") so a left_join per
# covariate adds one distinctly-named column with no extra renaming step.
load_lndcov <- function(file, value_col) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              !!value_col := .data[[value_col]]) %>%
    filter(!is.na(.data[[value_col]])) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

cat("\nLoading", K, "land-cover covariate lookup tables...\n")
cov_lookups <- map2(cov_meta$file, cov_meta$category, load_lndcov)
names(cov_lookups) <- cov_meta$category

for (cat_name in cov_meta$category) {
  lk <- cov_lookups[[cat_name]]
  cat("  ", cat_name, ": ", nrow(lk), " route-year rows (",
      length(unique(lk$route_key)), " routes)\n", sep = "")
}

# All 16 covariates are already route-level land-cover PROPORTIONS in [0,1]
# (confirmed min = 0, max < 1 in every file when this script was written) —
# no max-rescaling needed here, unlike grass_to_anthro in the 5-model
# pipeline, whose raw range was far narrower than [0,1].

spp_df16 <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                     stringsAsFactors = FALSE)

## ============================================================================
## Function: data-prep pipeline for one species (BBS data, join all K
## covariates, drop rows missing ANY of them, spatial neighbours). Mirrors
## prepare_species_data() in 1c_species_iCAR_covariates.R, generalized from 3
## covariates to K.
## ============================================================================
prepare_species_data_16 <- function(species, species_bbs, strat,
                                    firstYear, lastYear,
                                    cov_lookups, cov_meta) {

  cat("    Preparing data for:", species, "\n")

  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data <- data_pkg[["raw_data"]]
  strata_map <- load_map(strat)

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

  # Join all K land-cover covariates by route + calendar year -----------------
  n_before <- nrow(new_data)
  for (cat_name in cov_meta$category) {
    new_data <- new_data %>%
      left_join(cov_lookups[[cat_name]], by = c("route" = "route_key", "r_year" = "year"))
  }

  match_rate <- mean(complete.cases(new_data[, cov_meta$category, drop = FALSE]))
  cat("    Covariate match rate (all", nrow(cov_meta), "land-cover covariates present):",
      round(100 * match_rate, 1), "%\n")

  cat("    Sample new_data$route values:      ",
      paste(head(unique(new_data$route), 6), collapse = ", "), "\n")
  cat("    Sample covariate route_key values: ",
      paste(head(unique(cov_lookups[[cov_meta$category[1]]]$route_key), 6), collapse = ", "), "\n")

  if (match_rate < 0.5) {
    stop("Less than half of observations matched all ", nrow(cov_meta),
         " land-cover covariates for ", species, " — new_data$route doesn't ",
         "look like '<StateNum>-<Route>'. Compare the sample values printed ",
         "above (see the analogous check in prepare_species_data() in ",
         "1c_species_iCAR_covariates.R).")
  }

  # Drop rows missing ANY of the K covariates — same reduced-dataset
  # principle as the 5-model pipeline.
  new_data <- new_data %>% filter(complete.cases(new_data[, cov_meta$category, drop = FALSE]))
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

  # K-covariate matrix, ncounts x K, column order == cov_meta$category order
  # (this is also the order gamma_out[1..K] comes back in from Stan).
  X <- as.matrix(new_data[, cov_meta$category, drop = FALSE])
  storage.mode(X) <- "double"
  dimnames(X) <- NULL  # strip row/col names — not needed by Stan and one
                        # less difference from a "plain" numeric matrix

  # Defensive check: the Stan data block declares X as
  # matrix<lower=0, upper=1>[ncounts, K], so any NA/NaN/out-of-[0,1] value
  # here would otherwise only surface as an opaque chain crash inside
  # cmdstan. Fail loudly here instead, with the exact offending cell.
  if (anyNA(X)) {
    stop("X contains NA/NaN after the complete.cases() filter for ", species,
         " — this should be impossible; check the join logic above.")
  }
  if (any(X < 0 | X > 1)) {
    bad <- which(X < 0 | X > 1, arr.ind = TRUE)
    stop("X has ", nrow(bad), " value(s) outside [0,1] for ", species,
         " (e.g. row ", bad[1, 1], ", covariate '", cov_meta$category[bad[1, 2]],
         "' = ", X[bad[1, 1], bad[1, 2]], ") — Stan's <lower=0,upper=1> bound ",
         "on X would reject this.")
  }

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
    sd_alpha_prior = sd_alpha_prior,
    K          = ncol(X),
    X          = X
  )

  if (car_stan_dat$N != base_stan_data[["nroutes"]]) {
    stop("Some routes are missing from adjacency matrix")
  }

  list(new_data = new_data, base_stan_data = base_stan_data, n_dropped = n_dropped)
}

## ============================================================================
## Function: fit the 16-covariate model for one species.
## ============================================================================
fit_one_16covariate_model <- function(species, species_f, prepped,
                                      stan_model, firstYear, lastYear, cov_meta) {

  stan_data <- prepped$base_stan_data

  # Stash the exact data list that's about to be sent to Stan in the global
  # environment BEFORE calling $sample(), so it survives even if a later
  # step in this function errors out (e.g. summary() on an all-chains-failed
  # fit — see below) — without this, the prepared data (which can be slow to
  # rebuild: BBS pull + Voronoi neighbours) would be lost, leaving nothing to
  # debug with interactively. Inspect with e.g.:
  #   str(debug_stan_data_16cov)
  #   range(debug_stan_data_16cov$X)  # confirm all covariate values in [0,1]
  #   any(is.na(debug_stan_data_16cov$X))
  # and/or re-run a fast 1-chain/low-iteration debug fit directly:
  #   fit_dbg <- stan_model$sample(data = debug_stan_data_16cov, chains = 1,
  #                                iter_warmup = 50, iter_sampling = 50,
  #                                refresh = 1, show_exceptions = TRUE)
  assign("debug_stan_data_16cov", stan_data, envir = .GlobalEnv)
  assign("debug_species_16cov", species, envir = .GlobalEnv)

  cat("    Fitting 16-covariate model for", species, "...\n")

  stanfit <- stan_model$sample(
    data = stan_data,
    refresh = 400,
    chains = 4, iter_sampling = 2000, iter_warmup = 2000,
    parallel_chains = 4,
    adapt_delta = 0.8,
    max_treedepth = 10,
    show_exceptions = TRUE,
    output_dir = cmdstanr_output_dir
  )

  # $sample() itself does NOT error even when every chain crashes — it
  # returns a CmdStanMCMC object with the (empty) results, printing "Chain N
  # finished unexpectedly!" warnings along the way. It's the *next* call,
  # $summary()/$draws(), that raises "No chains finished successfully" when
  # there's nothing to summarize. So `stanfit` is stashed here, BEFORE that
  # call, specifically so a crash still leaves something to debug with:
  #   debug_stanfit_16cov$output(1)   # chain 1's full raw console log —
  #                                   # this is where the *actual* crash
  #                                   # reason (a genuine Stan/C++ exception
  #                                   # message, not just "unexpectedly!")
  #                                   # will be printed.
  #   debug_stanfit_16cov$output(2)   # ...and so on for chains 2-4
  assign("debug_stanfit_16cov", stanfit, envir = .GlobalEnv)

  summ <- stanfit$summary()
  fit_time <- round(stanfit$time()[["total"]] / 60, 1)
  cat("    Fit time:", fit_time, "minutes\n")

  out_base <- paste0(species_f, "_iCAR_16covariates_", firstYear, "_", lastYear)
  stanfit$save_object(file.path(output_dir16, paste0(out_base, "_stanfit.rds")))
  saveRDS(summ, file.path(output_dir16, paste0(out_base, "_summ_fit.rds")))

  sp_data_file <- here::here("data", "stan_data",
                             paste0(species_f, "_16covariates_",
                                    firstYear, "_", lastYear, "_stan_data.RData"))
  save(list = c("stan_data", "cov_meta"), file = sp_data_file)

  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess  <- min(summ$ess_bulk, na.rm = TRUE)

  cat("    Max Rhat:", round(max_rhat, 4), " | Min ESS:", round(min_ess, 0), "\n")

  # Per-covariate gamma_out[k] summary, labelled with the covariate names in
  # cov_meta$category (same column order used to build X above).
  gamma_rows <- summ %>%
    filter(str_detect(variable, "^gamma_out\\[")) %>%
    mutate(k = as.integer(str_match(variable, "^gamma_out\\[([0-9]+)\\]")[, 2])) %>%
    arrange(k) %>%
    mutate(species  = species,
           covariate = cov_meta$category[k],
           code      = cov_meta$code[k]) %>%
    select(species, covariate, code, mean, sd, q5, q95, rhat, ess_bulk)

  diagnostics_row <- data.frame(
    species  = species,
    max_rhat = round(max_rhat, 4),
    min_ess  = round(min_ess, 0),
    n_obs    = nrow(prepped$new_data),
    n_routes = max(prepped$new_data$routeF),
    n_dropped_missing_cov = prepped$n_dropped,
    fit_min  = fit_time
  )

  list(gamma_rows = gamma_rows, diagnostics_row = diagnostics_row)
}

## ============================================================================
## Compile the 16-covariate Stan model once
## ============================================================================
cat("\nCompiling models/slope_iCAR_route_NB_New_16covariates.stan...\n")
model_16cov <- cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_16covariates.stan"),
                             stanc_options = list("O1"))

## ============================================================================
## Main loop: for each species, prepare data + fit the 16-covariate model
## ============================================================================
gamma_list <- list()
diagnostics_list <- list()

for (i in seq_along(species_filter)) {
  sp <- species_filter[i]
  sp_f <- species_to_f(sp)

  cat("\n================================================================\n")
  cat("  [", i, "/", length(species_filter), "]", sp, "\n")
  cat("================================================================\n")

  sp_row <- spp_df16[spp_df16$Common.Name == sp, ][1, ]
  if (is.na(sp_row$Common.Name)) {
    message("  [ERROR] Species '", sp, "' not found in spp_names_codes_group_aou.csv — skipping.")
    next
  }
  sp_bbs <- sp_row$bbs_english

  out_base  <- paste0(sp_f, "_iCAR_16covariates_", firstYear, "_", lastYear)
  summ_file <- file.path(output_dir16, paste0(out_base, "_summ_fit.rds"))

  if (file.exists(summ_file) && !force_refit) {
    cat("  Skipping (already fitted):", basename(summ_file), "\n")
    next
  }

  prepped <- tryCatch(
    prepare_species_data_16(species = sp, species_bbs = sp_bbs, strat = strat,
                            firstYear = firstYear, lastYear = lastYear,
                            cov_lookups = cov_lookups, cov_meta = cov_meta),
    error = function(e) {
      message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(prepped)) next

  fit_result <- tryCatch(
    fit_one_16covariate_model(species = sp, species_f = sp_f, prepped = prepped,
                              stan_model = model_16cov, firstYear = firstYear,
                              lastYear = lastYear, cov_meta = cov_meta),
    error = function(e) {
      message("  [ERROR] 16-covariate model fit failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(fit_result)) next

  gamma_list[[sp]] <- fit_result$gamma_rows
  diagnostics_list[[sp]] <- fit_result$diagnostics_row

  cat("  Done — 16-covariate fit + diagnostics saved for", sp, "\n")
}

## ============================================================================
## Combined pilot outputs
## ============================================================================
if (length(gamma_list) > 0) {
  gamma_all <- bind_rows(gamma_list)
  gamma_csv <- file.path(output_dir16, "pilot_16covariates_gamma_summary.csv")
  write.csv(gamma_all, gamma_csv, row.names = FALSE)
  cat("\nGamma summary (species x covariate) written to:", gamma_csv, "\n")
  print(gamma_all)
}

if (length(diagnostics_list) > 0) {
  diagnostics_all <- bind_rows(diagnostics_list)
  diag_csv <- file.path(output_dir16, "pilot_16covariates_diagnostics.csv")
  write.csv(diagnostics_all, diag_csv, row.names = FALSE)
  cat("\nPer-species diagnostics written to:", diag_csv, "\n")
  print(diagnostics_all)
}

## ============================================================================
## Post-run checks
## ============================================================================
cat("\n############################################################\n")
cat("# TEST CHECKS\n")
cat("############################################################\n\n")

all_checks <- list()

for (sp in species_filter) {
  sp_f <- species_to_f(sp)
  checks <- list()

  out_base  <- paste0(sp_f, "_iCAR_16covariates_", firstYear, "_", lastYear)
  stanfit_f <- file.path(output_dir16, paste0(out_base, "_stanfit.rds"))
  summ_f    <- file.path(output_dir16, paste0(out_base, "_summ_fit.rds"))
  stan_data_f <- here::here("data", "stan_data",
                            paste0(sp_f, "_16covariates_", firstYear, "_", lastYear, "_stan_data.RData"))

  checks[["stanfit.rds exists"]]     <- file.exists(stanfit_f)
  checks[["summ_fit.rds exists"]]    <- file.exists(summ_f)
  checks[["stan_data.RData exists"]] <- file.exists(stan_data_f)

  if (exists("gamma_all")) {
    sp_gamma <- gamma_all[gamma_all$species == sp, ]
    checks[["gamma summary has all K covariates for this species"]] <-
      nrow(sp_gamma) == K
    checks[["gamma summary has non-NA mean for every covariate"]] <-
      nrow(sp_gamma) > 0 && all(!is.na(sp_gamma$mean))
  } else {
    checks[["gamma summary exists"]] <- FALSE
  }

  n_pass  <- sum(unlist(checks) == TRUE, na.rm = TRUE)
  n_total <- length(checks)

  cat("--- Results for", sp, "---\n")
  for (nm in names(checks)) {
    cat(sprintf("  [%s] %s\n", if (isTRUE(checks[[nm]])) "PASS" else "FAIL", nm))
  }
  cat(sprintf("  %d/%d checks passed for %s\n\n", n_pass, n_total, sp))

  all_checks[[sp]] <- data.frame(species = sp, n_pass = n_pass, n_total = n_total)
}

overall <- bind_rows(all_checks)
cat("############################################################\n")
cat("# OVERALL SUMMARY\n")
cat("############################################################\n")
print(overall)

total_pass  <- sum(overall$n_pass)
total_check <- sum(overall$n_total)
cat("\n", total_pass, "/", total_check, " checks passed across ",
    length(species_filter), " species\n", sep = "")

if (all(overall$n_pass == overall$n_total)) {
  cat("\nALL CHECKS PASSED for every species.\n")
} else {
  failed_species <- overall$species[overall$n_pass != overall$n_total]
  cat("\nSOME CHECKS FAILED for:", paste(failed_species, collapse = ", "), "— see above.\n")
}
