## =============================================================================
## Compare NB vs ZINB iCAR route-level trend models
## Time period: 2010-2025
##
## Motivation: 1b_posterior_predictive_checks.R found that ~27 of 39
## grassland species show a strong, consistent posterior-predictive-check
## signature of excess zeros under the plain NB model fit by
## 1_species_iCAR_2010_2025.R (p_prop_zero at or near 0, paired with
## p_sd/p_max near 1 — the model compensates for under-predicted zeros by
## over-dispersing the rest of the distribution). This script checks,
## per species, whether a zero-inflated NB (ZINB) model actually fits
## better once you account for the added complexity, rather than assuming
## it does.
##
## Runs on the full grasslands target group (all 39 species) by default —
## set `species_subset` below to a vector of Common.Name values to restrict
## to a subset (e.g. just the species flagged by the 2026-07-14 PPC run).
## For each species, this script:
##   1. Prepares the BBS data once (shared by both models).
##   2. Fits the full model with NB and with ZINB.
##   3. Runs leave-future-out CV (predict final year) for each model.
##   4. Computes comparison metrics: CV lppd, Spearman r, convergence (max
##      Rhat / min ESS), proportion of zeros in the data, and the estimated
##      zero-inflation probability (theta) for ZINB.
##
## Requires (new files, does not touch any existing model or script):
##   models/slope_iCAR_route_ZINB_New.stan
##   models/slope_iCAR_route_ZINB_cv.stan
## The NB counterparts (models/slope_iCAR_route_NB_New.stan,
## models/slope_iCAR_route_NB_cv.stan) already exist and are unchanged.
##
## Output:
##   output/<species>_iCAR_New_<firstYear>_<lastYear>_stanfit.rds   (NB)
##   output/<species>_iCAR_New_<firstYear>_<lastYear>_summ_fit.rds
##   output/<species>_iCAR_ZINB_<firstYear>_<lastYear>_stanfit.rds  (ZINB)
##   output/<species>_iCAR_ZINB_<firstYear>_<lastYear>_summ_fit.rds
##   output/model_comparison/<group>_<species>_NB_vs_ZINB.csv  (per species)
##   output/model_comparison_<group>_<firstYear>_<lastYear>.csv (combined)
##
## NB deliberately reuses the main pipeline's original "New" naming
## (1_species_iCAR_2010_2025.R): if that script has already fit a species,
## this script pulls the existing _stanfit/_summ_fit/_cv_results output from
## disk instead of refitting NB from scratch, and only fits ZINB (tagged
## "_ZINB_", which has no prior output to reuse). If no prior NB fit is
## found, it's fit here and saved under the same "New" name, so a later run
## of 1_species_iCAR_2010_2025.R would pick it up too. (NB's tag could be
## switched to a distinct "_NB_" name later; for now it must stay "New" to
## match what's already on disk.)
##
## Interpretation: prefer the model with the higher (less negative) CV lppd.
## A large theta and an NB with many excess zeros points to ZINB; if theta is
## near 0 and lppd is similar, NB is the more parsimonious choice. Note: PPC
## flagged Swainson's Hawk for a *different* reason (under-, not
## over-dispersion at the low end — p_sd/p_max near 0, not p_prop_zero near
## 0), so ZINB is not expected to help there; its result should be read with
## that caveat rather than as a failure of the comparison.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("3_compare_NB_ZINB.R")
source("functions/neighbours_define_voronoi.R")
source("functions/posterior_summary_functions.R")

cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings (match 1_species_iCAR_2010_2025.R) ------------------------------
land_cover  <- "grasslands"
firstYear   <- 2010
lastYear    <- 2025
strat       <- "bcr"
force_refit <- TRUE

# Species to compare. NULL runs the full grasslands target group (39
# species). Set to a vector of Common.Name values to restrict to a subset,
# e.g. just the species flagged by 1b_posterior_predictive_checks.R on
# 2026-07-14:
#   c("Baird's Sparrow", "Bobolink", "Cassin's Sparrow",
#     "Chestnut-collared Longspur", "Clay-colored Sparrow", "Dickcissel",
#     "Eastern Kingbird", "Eastern Meadowlark", "Grasshopper Sparrow",
#     "Gray Partridge", "Horned Lark", "Lark Bunting", "Le Conte's Sparrow",
#     "Loggerhead Shrike", "Long-billed Curlew", "McCown's Longspur",
#     "Northern Bobwhite", "Ring-necked Pheasant", "Savannah Sparrow",
#     "Scissor-tailed Flycatcher", "Sedge Wren", "Sprague's Pipit",
#     "Swainson's Hawk", "Upland Sandpiper", "Vesper Sparrow",
#     "Western Kingbird", "Western Meadowlark")
species_subset <- NULL

output_dir <- here::here("output")
if (!dir.exists(output_dir)) dir.create(output_dir)
if (!dir.exists(here::here("data", "stan_data"))) dir.create(here::here("data", "stan_data"), recursive = TRUE)
cmp_dir <- here::here("output", "model_comparison")
if (!dir.exists(cmp_dir)) dir.create(cmp_dir, recursive = TRUE)

# Target group species list -----------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

if (!is.null(species_subset)) {
  target_spp <- target_spp %>% filter(Common.Name %in% species_subset)
}

cat("=== NB vs ZINB model comparison ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("Species to compare (n =", nrow(target_spp), ")\n")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Compile all four models once --------------------------------------------
models <- list(
  NB = list(
    full = cmdstan_model("models/slope_iCAR_route_NB_New.stan", stanc_options = list("O1")),
    cv   = cmdstan_model("models/slope_iCAR_route_NB_cv.stan",  stanc_options = list("O1"))
  ),
  ZINB = list(
    full = cmdstan_model("models/slope_iCAR_route_ZINB_New.stan", stanc_options = list("O1")),
    cv   = cmdstan_model("models/slope_iCAR_route_ZINB_cv.stan",  stanc_options = list("O1"))
  )
)

# ==========================================================================
# Data prep for one species (shared by both models). Returns the pieces both
# the full fit and the CV fit need, or NULL on failure / too-few-data.
# ==========================================================================
prep_species <- function(species_bbs, species) {
  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data <- data_pkg[["raw_data"]]
  strata_map <- load_map(strat)

  # Use all available routes (no spatial route filtering)
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

  if (nrow(new_data) < 50) return(NULL)

  realized_strata_map <- strata_map %>%
    filter(strata_name %in% unique(new_data$strat_name))

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
    save_plot_data  = FALSE
  )

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

  if (car_stan_dat$N != stan_data[["nroutes"]]) return(NULL)

  list(new_data = new_data, route_map = route_map,
       car_stan_dat = car_stan_dat, stan_data = stan_data,
       sd_alpha_prior = sd_alpha_prior)
}

# ==========================================================================
# Build the CV (leave-future-out) stan data list from prepped pieces.
# ==========================================================================
build_cv_data <- function(prep) {
  new_data       <- prep$new_data
  route_map      <- prep$route_map
  car_stan_dat   <- prep$car_stan_dat
  sd_alpha_prior <- prep$sd_alpha_prior
  ynext          <- lastYear

  obs_df_fit <- new_data %>%
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

  obs_df_predict <- new_data %>%
    filter(r_year == ynext) %>%
    left_join(., obs_df, by = "ObsN") %>%
    mutate(observer = ifelse(!is.na(observer), observer, 0))

  stan_data_cv[["route_pred"]]    <- obs_df_predict$routeF
  stan_data_cv[["count_pred"]]    <- obs_df_predict$count
  stan_data_cv[["firstyr_pred"]]  <- obs_df_predict$firstyr
  stan_data_cv[["observer_pred"]] <- obs_df_predict$observer
  stan_data_cv[["ncounts_pred"]]  <- length(obs_df_predict$count)

  list(stan_data_cv = stan_data_cv, obs_df_predict = obs_df_predict)
}

# ==========================================================================
# Fit one model type (full + CV) and return its metrics.
# ==========================================================================
fit_one_model <- function(model_type, prep, cv_pieces, species_f) {
  mods <- models[[model_type]]

  # File naming: NB reuses the main pipeline's original "New" naming
  # (1_species_iCAR_2010_2025.R) so an already-fitted NB model can be pulled
  # from disk instead of refit here — NB fits are the expensive part and
  # 1_species_iCAR_2010_2025.R has usually already done them for every
  # species. ZINB has no prior fit to reuse, so it gets its own "_ZINB_"
  # tag. (NB's tag could switch to "_NB_" later; for now it must stay "New"
  # to match what's already on disk and avoid re-running NB fits.)
  name_tag        <- if (model_type == "NB") "New" else model_type
  out_base        <- paste0(species_f, "_iCAR_", name_tag, "_", firstYear, "_", lastYear)
  stanfit_file    <- file.path(output_dir, paste0(out_base, "_stanfit.rds"))
  summ_file       <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
  cv_results_file <- file.path(output_dir, paste0(species_f, "_iCAR_cv_results_", lastYear, ".rds"))

  # Full fit (for convergence diagnostics) ---------------------------------
  if (model_type == "NB" && file.exists(summ_file)) {
    cat("    [", model_type, "] Reusing existing fit:", basename(summ_file), "\n")
    summ <- readRDS(summ_file)
  } else {
    stanfit <- mods$full$sample(
      data = prep$stan_data,
      refresh = 0,
      chains = 4, iter_sampling = 2000, iter_warmup = 2000,
      parallel_chains = 4,
      adapt_delta = 0.8,
      max_treedepth = 10,
      show_exceptions = FALSE,
      output_dir = cmdstanr_output_dir
    )
    summ <- stanfit$summary()

    # Save fit + summary, tagged with model_type (New / ZINB) so this never
    # overwrites a different model's output.
    stanfit$save_object(stanfit_file)
    saveRDS(summ, summ_file)
  }

  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess  <- min(summ$ess_bulk, na.rm = TRUE)
  cat("    [", model_type, "] Max Rhat:", round(max_rhat, 4),
      " | Min ESS:", round(min_ess, 0), "\n")

  theta <- NA
  if ("theta" %in% summ$variable) {
    theta <- summ$mean[summ$variable == "theta"]
  }

  # CV fit -----------------------------------------------------------------
  stan_data_cv   <- cv_pieces$stan_data_cv
  obs_df_predict <- cv_pieces$obs_df_predict
  cv_lppd <- NA
  cv_cor  <- NA

  if (model_type == "NB" && file.exists(cv_results_file)) {
    cat("    [", model_type, "] Reusing existing CV results:", basename(cv_results_file), "\n")
    cv_results <- readRDS(cv_results_file)
    cv_lppd <- mean(cv_results$log_lik_mean)
    cv_cor  <- cor(cv_results$count, cv_results$E_pred_count, method = "spearman")
    cat("    [", model_type, "] CV lppd:", round(cv_lppd, 4),
        " | Spearman r:", round(cv_cor, 3), "\n")
  } else if (stan_data_cv[["ncounts_pred"]] > 0) {
    cv_fit <- mods$cv$sample(
      data = stan_data_cv,
      refresh = 0,
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
    cat("    [", model_type, "] CV lppd:", round(cv_lppd, 4),
        " | Spearman r:", round(cv_cor, 3), "\n")
  }

  data.frame(
    model     = model_type,
    max_rhat  = round(max_rhat, 4),
    min_ess   = round(min_ess, 0),
    cv_lppd   = ifelse(is.na(cv_lppd), NA, round(cv_lppd, 4)),
    cv_cor    = ifelse(is.na(cv_cor), NA, round(cv_cor, 3)),
    theta     = ifelse(is.na(theta), NA, round(theta, 4))
  )
}

# ==========================================================================
# Main loop
# ==========================================================================
comparison_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]
  sp_bbs  <- target_spp$bbs_english[i]

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(target_spp), "]", sp, "\n")
  cat("================================================================\n")

  sp_csv <- file.path(cmp_dir, paste0(land_cover, "_", sp_f, "_NB_vs_ZINB.csv"))
  if (file.exists(sp_csv) && !force_refit) {
    cat("  Skipping (already exists):", basename(sp_csv), "\n")
    comparison_list[[sp]] <- read.csv(sp_csv, stringsAsFactors = FALSE)
    next
  }

  res <- tryCatch({
    prep <- prep_species(sp_bbs, sp)
    if (is.null(prep)) stop("insufficient data / adjacency mismatch")

    prop_zero <- mean(prep$new_data$count == 0)
    cat("  Routes:", prep$stan_data$nroutes,
        " | Obs:", prep$stan_data$ncounts,
        " | Prop zeros:", round(prop_zero, 3), "\n")

    cv_pieces <- build_cv_data(prep)

    cat("  Fitting NB...\n")
    nb   <- fit_one_model("NB",   prep, cv_pieces, sp_f)
    cat("  Fitting ZINB...\n")
    zinb <- fit_one_model("ZINB", prep, cv_pieces, sp_f)

    cmp <- bind_rows(nb, zinb) %>%
      mutate(species      = sp,
             species_code = sp_code,
             group        = land_cover,
             prop_zero    = round(prop_zero, 3),
             nroutes      = prep$stan_data$nroutes,
             ncounts      = prep$stan_data$ncounts) %>%
      select(species, species_code, group, model, nroutes, ncounts,
             prop_zero, theta, cv_lppd, cv_cor, max_rhat, min_ess)

    # Flag the preferred model by CV lppd (higher = better)
    if (all(!is.na(cmp$cv_lppd))) {
      best <- cmp$model[which.max(cmp$cv_lppd)]
      cmp$preferred_by_lppd <- ifelse(cmp$model == best, "*", "")
      d_lppd <- round(max(cmp$cv_lppd) - min(cmp$cv_lppd), 4)
      cat("  Preferred by CV lppd:", best, " (delta =", d_lppd, ")\n")
    } else {
      cmp$preferred_by_lppd <- NA
    }

    write.csv(cmp, sp_csv, row.names = FALSE)
    cmp
  }, error = function(e) {
    message("  [ERROR] ", sp, ": ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(res)) comparison_list[[sp]] <- res
}

# Combined comparison CSV --------------------------------------------------
if (length(comparison_list) > 0) {
  comparison_all <- bind_rows(comparison_list)
  cmp_csv <- file.path(output_dir,
                       paste0("model_comparison_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))
  write.csv(comparison_all, cmp_csv, row.names = FALSE)

  cat("\n=== Summary ===\n")
  cat("Species compared:", length(comparison_list), "\n")
  cat("Combined comparison CSV:", cmp_csv, "\n")

  # Quick tally of which model wins per species
  wins <- comparison_all %>%
    filter(preferred_by_lppd == "*") %>%
    count(model, name = "n_species_preferred")
  print(wins)
} else {
  cat("\nNo species successfully compared.\n")
}
