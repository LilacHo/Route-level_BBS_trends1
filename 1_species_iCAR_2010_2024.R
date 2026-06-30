## =============================================================================
## Grassland birds: iCAR route-level trend model (fitting only)
## Time period: 2010-2024
##
## Structure mirrors 1_CH_no_habitat_routes.R:
##   - Settings block
##   - fit_species(): full data-prep + model-fit pipeline for one species
##   - Main loop: fit per species, save the fit + stan_data, collect diagnostics
##
## This script FITS the models and saves, per species:
##   output/<species>_iCAR_New_<firstYear>_<lastYear>_stanfit.rds
##   output/<species>_iCAR_New_<firstYear>_<lastYear>_summ_fit.rds
##   data/stan_data/<species>_<firstYear>_<lastYear>_stan_data.RData
##
## Per-route trend CSVs are generated separately by
## 2a_generate_route_trend_csvs.R (reads the saved fits above). Trend maps by
## 2b_plot_iCAR_trend_maps.R (reads those CSVs).
##
## Diagnostics (Rhat, ESS, leave-future-out CV) are kept: printed to console,
## stored per species, and written to a combined diagnostics CSV at the end.
##
## Based on Fitting_new_iCAR_slope_model.R (sum_to_zero_vector parameterization)
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("1_species_iCAR_2010_2024.R")
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
lastYear   <- 2024
dt         <- lastYear - firstYear

strat <- "bbs_usgs"

# Re-fit control: if TRUE, ignore existing fits and re-run the model (needed
# after changing the model, priors, or data thresholds).
force_refit <- TRUE

# Output directories
output_dir <- here::here("output")
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(here::here("data"))) dir.create(here::here("data"))
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)

# Target group species list -----------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== iCAR route-level trend model ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat(land_cover, "species to fit (n =", nrow(target_spp), "):\n")
print(target_spp %>% select(Common.Name, Code, bbs_english))

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Compile the Stan models once (reused across species) --------------------
mod.file <- "models/slope_iCAR_route_NB_New.stan"
slope_model <- cmdstan_model(mod.file, stanc_options = list("O1"))

cv_mod.file <- "models/slope_iCAR_route_NB_cv.stan"
cv_model <- cmdstan_model(cv_mod.file, stanc_options = list("O1"))

# ==========================================================================
# Function: run the full data-prep + fit pipeline for one species
# Returns a list with summ, route_map, realized_strata_map, new_data, and a
# one-row diagnostics data.frame (max_rhat, min_ess, cv_lppd, cv_cor).
# ==========================================================================
fit_species <- function(species, species_bbs, species_f, land_cover, strat,
                         firstYear, lastYear, slope_model, cv_model) {

  cat("    Running full data-prep + model fit for:", species, "\n")

  # BBS data ---------------------------------------------------------------
  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data <- data_pkg[["raw_data"]]
  strata_map <- load_map(strat)

  # Use all available routes (no spatial route filtering) ------------------
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

  # Build Stan data --------------------------------------------------------
  tmp <- data.frame(route = new_data$routeF, count = new_data$count) %>%
    group_by(route) %>%
    summarise(mean_count = mean(log(count + 1)))
  sd_alpha_prior <- sd(tmp$mean_count)

  stan_data <- list(
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

  if (car_stan_dat$N != stan_data[["nroutes"]]) {
    stop("Some routes are missing from adjacency matrix")
  }

  # Save stan_data ---------------------------------------------------------
  sp_data_file <- here::here("data", "stan_data",
                             paste0(species_f, "_",
                                    firstYear, "_", lastYear, "_stan_data.RData"))
  save(list = c("stan_data", "new_data", "route_map", "realized_strata_map",
                "car_stan_dat", "sd_alpha_prior"),
       file = sp_data_file)

  # Fit iCAR model ---------------------------------------------------------
  stanfit <- slope_model$sample(
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

  out_base <- paste0(species_f, "_iCAR_New_", firstYear, "_", lastYear)
  stanfit$save_object(file.path(output_dir, paste0(out_base, "_stanfit.rds")))
  saveRDS(summ, file.path(output_dir, paste0(out_base, "_summ_fit.rds")))

  # Convergence diagnostics ------------------------------------------------
  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess  <- min(summ$ess_bulk, na.rm = TRUE)
  cat("    Max Rhat:", round(max_rhat, 4),
      " | Min ESS:", round(min_ess, 0), "\n")

  # Leave-future-out CV (final year only) ----------------------------------
  ynext <- lastYear
  full_data <- new_data

  obs_df_fit <- full_data %>%
    filter(r_year <= ynext - 1) %>%
    mutate(observer = as.integer(factor(ObsN)))

  stan_data_cv <- list(
    count      = obs_df_fit$count,
    year       = obs_df_fit$year,
    route      = obs_df_fit$routeF,
    firstyr    = obs_df_fit$firstyr,
    observer   = obs_df_fit$observer,
    nobservers = max(obs_df_fit$observer),
    nyears     = max(obs_df_fit$year),
    nroutes    = nrow(route_map),
    ncounts    = length(obs_df_fit$count),
    fixedyear  = floor(max(obs_df_fit$year) / 2),
    N_edges    = car_stan_dat$N_edges,
    node1      = car_stan_dat$node1,
    node2      = car_stan_dat$node2,
    sd_alpha_prior = sd_alpha_prior
  )

  obs_df <- obs_df_fit %>% select(observer, ObsN) %>% distinct()

  obs_df_predict <- full_data %>%
    filter(r_year == ynext) %>%
    left_join(., obs_df, by = "ObsN") %>%
    mutate(observer = ifelse(!is.na(observer), observer, 0))

  stan_data_cv[["route_pred"]]    <- obs_df_predict$routeF
  stan_data_cv[["count_pred"]]    <- obs_df_predict$count
  stan_data_cv[["firstyr_pred"]]  <- obs_df_predict$firstyr
  stan_data_cv[["observer_pred"]] <- obs_df_predict$observer
  stan_data_cv[["ncounts_pred"]]  <- length(obs_df_predict$count)

  cv_lppd <- NA
  cv_cor  <- NA

  if (stan_data_cv[["ncounts_pred"]] > 0) {
    cv_fit <- cv_model$sample(
      data = stan_data_cv,
      refresh = 400,
      chains = 3, iter_sampling = 1000, iter_warmup = 1000,
      parallel_chains = 3,
      adapt_delta = 0.8,
      max_treedepth = 10,
      show_exceptions = FALSE,
      output_dir = cmdstanr_output_dir
    )

    log_lik_full <- posterior_samples(fit = cv_fit, parm = "log_lik", dims = "i")
    log_lik_summ <- posterior_sums(log_lik_full, quantiles = NULL, dims = "i")
    names(log_lik_summ) <- paste0("log_lik_", names(log_lik_summ))

    E_pred_full <- posterior_samples(fit = cv_fit, parm = "E_pred", dims = "i")
    E_pred_summ <- posterior_sums(E_pred_full, quantiles = NULL, dims = "i")
    names(E_pred_summ) <- paste0("E_pred_", names(E_pred_summ))

    cv_results <- bind_cols(obs_df_predict, log_lik_summ, E_pred_summ)
    cv_results$E_pred_count <- exp(cv_results$E_pred_mean)

    cv_lppd <- mean(cv_results$log_lik_mean)
    cv_cor  <- cor(cv_results$count, cv_results$E_pred_count, method = "spearman")

    cat("    CV lppd:", round(cv_lppd, 4),
        " | Spearman r:", round(cv_cor, 3), "\n")

    saveRDS(cv_results, file.path(output_dir,
                                  paste0(species_f, "_iCAR_cv_results_", ynext, ".rds")))
  } else {
    cat("    No prediction data for year", ynext, "— CV skipped.\n")
  }

  diagnostics <- data.frame(
    species   = species,
    max_rhat  = round(max_rhat, 4),
    min_ess   = round(min_ess, 0),
    cv_lppd   = ifelse(is.na(cv_lppd), NA, round(cv_lppd, 4)),
    cv_cor    = ifelse(is.na(cv_cor), NA, round(cv_cor, 3)),
    fit_min   = fit_time
  )

  list(summ = summ, route_map = route_map,
       realized_strata_map = realized_strata_map,
       new_data = new_data, diagnostics = diagnostics)
}

# ==========================================================================
# Main loop: compute per-route trends for each species in the group
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

  # Skip if this species' fit already exists
  out_base       <- paste0(sp_f, "_iCAR_New_", firstYear, "_", lastYear)
  summ_file      <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
  stan_data_file <- here::here("data", "stan_data",
                               paste0(sp_f, "_", firstYear, "_", lastYear, "_stan_data.RData"))
  if (file.exists(summ_file) && file.exists(stan_data_file) && !force_refit) {
    cat("  Skipping (already fitted):", basename(summ_file), "\n")
    next
  }

  # Run the full data-prep + model-fit pipeline
  diagnostics <- NULL
  fit_result <- tryCatch(
    fit_species(species = sp,
                species_bbs = sp_bbs,
                species_f = sp_f,
                land_cover = land_cover,
                strat = strat,
                firstYear = firstYear,
                lastYear = lastYear,
                slope_model = slope_model,
                cv_model = cv_model),
    error = function(e) {
      message("  [ERROR] Failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(fit_result)) next
  diagnostics <- fit_result$diagnostics

  results_list[[sp]] <- summ_file

  # Diagnostics (kept) -----------------------------------------------------
  if (!is.null(diagnostics)) diagnostics_list[[sp]] <- diagnostics

  # Per-route trend CSVs are generated separately in
  # 2a_generate_route_trend_csvs.R from the saved fits; trend maps in
  # 2b_plot_iCAR_trend_maps.R from those CSVs.

  cat("  Done — fit + diagnostics saved for", sp, "\n")
}

# Combined diagnostics CSV (kept) -----------------------------------------
if (length(diagnostics_list) > 0) {
  diagnostics_all <- bind_rows(diagnostics_list)
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_", land_cover, "_",
                               firstYear, "_", lastYear, ".csv"))
  write.csv(diagnostics_all, diag_csv, row.names = FALSE)
  cat("\nDiagnostics written to:", diag_csv, "\n")
  print(diagnostics_all)
}

cat("\n=== Summary ===\n")
cat("Species fitted this run:", length(results_list), "\n")
cat("Fits saved in:", output_dir, "\n")
cat("Next: run 2a_generate_route_trend_csvs.R to build per-route CSVs.\n")
