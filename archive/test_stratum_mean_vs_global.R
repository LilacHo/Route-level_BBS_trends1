## =============================================================================
## Test: does using a STRATUM-LEVEL mean (per BCR) instead of a single GLOBAL
## mean for the route-level intercept/slope change the fitted trends?
##
## ARCHIVED (2026-07-21): tested across 6 species — route-level trend
## estimates were essentially unchanged (correlation 0.995-0.999) between
## the production NB model (models/slope_iCAR_route_NB_New.stan, the model
## used by 1_species_iCAR_2010_2025.R) and the stratum-mean alternative
## (models/slope_iCAR_route_NB_stratmean.stan). BCR strata differed
## credibly in baseline abundance (intercept) but not in trend (slope) for
## any species tested, and the stratum-mean model converged less reliably
## (up to 8.7% divergent transitions vs. none for the NB model). We
## retained the single global-mean NB model for route-level trend
## estimation — see the judgement statement near the end of this header
## block. Kept here for reference / re-running if the question is
## revisited. Originally tests/test_stratum_mean_vs_global.R; moved to
## archive/ alongside 1c_compare_NB_ZINB.R as a completed model-comparison
## test, not an active test suite file.
##
## Background:
##   models/slope_iCAR_route_NB_New.stan (the model used by
##   1_species_iCAR_2010_2025.R) gives every route the same two population-
##   level means: alpha_r = ALPHA + alpha_r' and beta_r = BETA + beta_r',
##   where alpha_r'/beta_r' are the route-level iCAR spatial deviations and
##   ALPHA/BETA are single scalars shared by ALL routes, regardless of which
##   BCR stratum they fall in.
##
##   models/slope_iCAR_route_NB_stratmean.stan replaces ALPHA/BETA with
##   stratum-indexed vectors ALPHA_strat[nstrata]/BETA_strat[nstrata] (one
##   mean per BCR), partially pooled toward a grand mean (MU_ALPHA/MU_BETA)
##   with an estimated among-stratum sd (sd_ALPHA_strat/sd_BETA_strat). This
##   nests the global-mean model as the special case where that sd -> 0, so
##   this script gives TWO complementary answers to "does it change things",
##   per species:
##     (1) directly: correlation / mean-abs-difference between the two
##         models' route-level trend & relative-abundance estimates
##     (2) implicitly: the posterior for sd_ALPHA_strat / sd_BETA_strat in
##         the stratum-mean model — if it's concentrated near zero, strata
##         don't differ meaningfully from one grand mean
##
## This script runs BOTH models, on IDENTICAL data/stan_data, for EACH
## species in `species_list` below (mirroring fit_species() in
## 1_species_iCAR_2010_2025.R for data prep), then pools the per-species
## evidence into one combined table and a short, ready-to-paste paragraph
## for the supplemental materials explaining which parameterization is
## preferred and why.
##
## Outputs (under output/stratum_mean_test/):
##   <species>_global_stanfit.rds / _stratmean_stanfit.rds   (full posterior)
##   <species>_global_summ.rds    / _stratmean_summ.rds      (summary tables)
##   <species>_route_comparison.csv   (per-route trend/rel_abundance, both
##                                      models, side by side + differences)
##   <species>_strata_summary.csv     (ALPHA_strat/BETA_strat per stratum,
##                                      plus MU_ALPHA/MU_BETA/sd_*_strat)
##   species_stratum_vs_global_summary.csv   (one row per species, combined)
##   supplemental_statement.txt              (auto-written paragraph)
##
## Judgement statement (methods-section wording, 2026-07-21):
## To evaluate whether route-level intercepts and trends should vary by
## Bird Conservation Region (BCR) rather than share a single range-wide
## mean, we compared the production NB model (single global intercept and
## trend, ALPHA/BETA) against an alternative with BCR-specific means
## (ALPHA_strat/BETA_strat) partially pooled toward a grand mean, fit to
## six species spanning a range of sample sizes (156-1,821 routes; 5-29
## realized BCR strata). The two parameterizations produced nearly
## identical route-level trend estimates (correlation = 0.995-0.999; mean
## absolute difference = 0.15-0.37 percentage points/year across species).
## While BCR strata differed credibly in baseline relative abundance
## (posterior 5th percentile of the among-stratum intercept SD exceeded
## zero for all six species), none showed credible among-stratum variation
## in trend (posterior 5th percentile of the among-stratum slope SD <=
## 0.008 for all species, versus intercept values of 0.06-0.46). The
## stratum-mean model also converged less reliably than the global-mean
## model, with lower effective sample sizes and up to 8.7% divergent
## transitions in the smallest dataset tested (Chestnut-collared Longspur,
## 156 routes), compared with no divergences under the global-mean model.
## Because trend estimates were insensitive to this parameterization and
## BCR-level heterogeneity was confined to baseline abundance rather than
## trend, we retained the single global-mean NB model (ALPHA, BETA;
## models/slope_iCAR_route_NB_New.stan) for route-level trend estimation.
## =============================================================================

rm(list = ls())

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("archive/test_stratum_mean_vs_global.R")
source("functions/neighbours_define_voronoi.R")

cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings --------------------------------------------------------------------
# List of species to test. Pick a handful that span the range of situations
# you care about (rare vs. common, few vs. many realized BCR strata, etc.) —
# a handful of species is enough to say something general in a supplement
# without paying for the full production species list.
species_list <- c(
  "Baird's Sparrow",
  "Grasshopper Sparrow",
  "Bobolink",
  "Chestnut-collared Longspur",
  "Eastern Meadowlark",
  "Dickcissel"
)

strat      <- "bcr"               # matches 1_species_iCAR_2010_2025.R
firstYear  <- 2010
lastYear   <- 2025

# MCMC settings — same as production by default; lower for a quick smoke test
chains         <- 4
iter_warmup    <- 2000
iter_sampling  <- 2000

# Thresholds used only to phrase the automated supplemental statement below;
# they don't affect the fits or the per-species CSVs. Adjust to taste.
sd_strat_threshold     <- 0.05   # among-stratum sd (log scale) considered "meaningful"
trend_cor_threshold    <- 0.98   # route-trend correlation considered "essentially identical"
trend_mad_threshold    <- 0.5    # mean abs trend difference (pct pts/yr) considered "negligible"
divergence_pct_threshold <- 1    # % divergent transitions considered a convergence problem

# Re-fit control: if TRUE, ignore existing per-species output and re-run.
force_refit <- FALSE

out_dir <- here::here("output", "stratum_mean_test")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)

spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

cat("=== Stratum-mean vs global-mean test ===\n")
cat("Species (n =", length(species_list), "):", paste(species_list, collapse = ", "), "\n")
cat("Strat:", strat, " | Period:", firstYear, "-", lastYear, "\n")

# Compile the Stan models once (reused across species) ------------------------
cat("\nCompiling models...\n")
model_global    <- cmdstan_model("models/slope_iCAR_route_NB_New.stan",
                                 stanc_options = list("O1"))
model_stratmean <- cmdstan_model("models/slope_iCAR_route_NB_stratmean.stan",
                                 stanc_options = list("O1"))

# ==============================================================================
# Function: run the full data-prep + both-model fit + comparison pipeline for
# ONE species. Returns a one-row summary data.frame (for the combined table)
# plus writes that species' full RDS/CSV outputs, exactly as the single-
# species version of this script did.
# ==============================================================================
fit_one_species <- function(species, strat, firstYear, lastYear,
                            model_global, model_stratmean,
                            chains, iter_warmup, iter_sampling,
                            cmdstanr_output_dir, out_dir, spp_df) {

  sp_f <- species_to_f(species)
  cat("\n----------------------------------------------------------------\n")
  cat("Species:", species, "\n")
  cat("----------------------------------------------------------------\n")

  # Step 1: look up the BBS name bbsBayes2 expects ----------------------------
  sp_row <- spp_df[spp_df$Common.Name == species, ][1, ]
  if (is.na(sp_row$Common.Name)) {
    stop("Species '", species, "' not found in spp_names_codes_group_aou.csv")
  }
  species_bbs <- sp_row$bbs_english

  # Step 2: data prep (mirrors fit_species() in 1_species_iCAR_2010_2025.R) ---
  cat("[1/4] Pulling BBS data + building spatial neighbours...\n")

  data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                  use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear,
                            max_year = lastYear,
                            min_n_routes = 1,
                            min_max_route_years = 1)

  raw_data   <- data_pkg[["raw_data"]]
  strata_map <- load_map(strat)

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

  if (car_stan_dat$N != max(new_data$routeF)) {
    stop("Some routes are missing from adjacency matrix")
  }

  # Step 3: build shared stan_data, plus the stratum-mean extras --------------
  tmp <- data.frame(route = new_data$routeF, count = new_data$count) %>%
    group_by(route) %>%
    summarise(mean_count = mean(log(count + 1)))
  sd_alpha_prior <- sd(tmp$mean_count)

  stan_data_base <- list(
    count      = new_data$count,
    ncounts    = nrow(new_data),
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

  # Per-route stratum index (each physical route belongs to exactly one BCR
  # stratum, so this should have exactly nroutes rows, one per routeF)
  route_strat_map <- new_data %>%
    distinct(routeF, strat_name) %>%
    arrange(routeF)

  if (nrow(route_strat_map) != stan_data_base$nroutes) {
    stop("A route maps to more than one stratum (", nrow(route_strat_map),
         " rows for ", stan_data_base$nroutes, " routes) — check new_data.")
  }

  route_strat_map$strat_idx <- as.integer(factor(route_strat_map$strat_name))
  nstrata <- length(unique(route_strat_map$strat_idx))

  cat("    Strata realized (BCRs):", nstrata, "\n")

  if (nstrata < 2) {
    stop("Only ", nstrata, " realized BCR stratum — a stratum-mean model ",
         "can't be compared to the global-mean model for this species.")
  }

  stan_data_global    <- stan_data_base
  stan_data_stratmean <- c(stan_data_base,
                           list(nstrata = nstrata,
                                route_strat = route_strat_map$strat_idx))

  # Step 4: fit both models on the identical data ------------------------------
  cat("[2/4] Fitting GLOBAL-mean model...\n")
  fit_global <- model_global$sample(
    data = stan_data_global,
    refresh = 400,
    chains = chains, iter_sampling = iter_sampling, iter_warmup = iter_warmup,
    parallel_chains = chains,
    adapt_delta = 0.8,
    max_treedepth = 10,
    show_exceptions = FALSE,
    output_dir = cmdstanr_output_dir
  )
  summ_global <- fit_global$summary()

  cat("[2/4] Fitting STRATUM-mean model...\n")
  fit_stratmean <- model_stratmean$sample(
    data = stan_data_stratmean,
    refresh = 400,
    chains = chains, iter_sampling = iter_sampling, iter_warmup = iter_warmup,
    parallel_chains = chains,
    adapt_delta = 0.8,
    max_treedepth = 10,
    show_exceptions = FALSE,
    output_dir = cmdstanr_output_dir
  )
  summ_stratmean <- fit_stratmean$summary()

  fit_global$save_object(file.path(out_dir, paste0(sp_f, "_global_stanfit.rds")))
  saveRDS(summ_global, file.path(out_dir, paste0(sp_f, "_global_summ.rds")))
  fit_stratmean$save_object(file.path(out_dir, paste0(sp_f, "_stratmean_stanfit.rds")))
  saveRDS(summ_stratmean, file.path(out_dir, paste0(sp_f, "_stratmean_summ.rds")))

  # Sampler diagnostics (divergences etc.) — the stratum-mean model has more
  # parameters (nstrata extra means + 2 more hierarchical sds) and can be
  # harder to sample than the global-mean model. Divergences here are a real
  # cost to weigh against any gain from the extra parameters, and they also
  # undermine trust in the very sd_ALPHA_strat/sd_BETA_strat estimates used
  # below to argue strata differ — so they're tracked, not just Rhat/ESS.
  diag_global    <- fit_global$diagnostic_summary(quiet = TRUE)
  diag_stratmean <- fit_stratmean$diagnostic_summary(quiet = TRUE)
  n_transitions  <- chains * iter_sampling

  num_divergent_global    <- sum(diag_global$num_divergent)
  num_divergent_stratmean <- sum(diag_stratmean$num_divergent)
  pct_divergent_global    <- 100 * num_divergent_global    / n_transitions
  pct_divergent_stratmean <- 100 * num_divergent_stratmean / n_transitions

  saveRDS(list(global = diag_global, stratmean = diag_stratmean),
          file.path(out_dir, paste0(sp_f, "_diagnostics.rds")))

  max_rhat_global <- max(summ_global$rhat, na.rm = TRUE)
  min_ess_global  <- min(summ_global$ess_bulk, na.rm = TRUE)
  max_rhat_strat  <- max(summ_stratmean$rhat, na.rm = TRUE)
  min_ess_strat   <- min(summ_stratmean$ess_bulk, na.rm = TRUE)

  cat("    Global mean    — max Rhat:", round(max_rhat_global, 4),
      " | min ESS:", round(min_ess_global, 0),
      " | divergences:", num_divergent_global,
      paste0("(", round(pct_divergent_global, 2), "%)"), "\n")
  cat("    Stratum mean   — max Rhat:", round(max_rhat_strat, 4),
      " | min ESS:", round(min_ess_strat, 0),
      " | divergences:", num_divergent_stratmean,
      paste0("(", round(pct_divergent_stratmean, 2), "%)"), "\n")

  # Step 5: compare route-level trend / relative abundance --------------------
  cat("[3/4] Comparing route-level estimates...\n")

  extract_route_trends <- function(summ, model_label) {
    beta_summ <- summ %>%
      filter(str_detect(variable, "^beta\\[")) %>%
      transmute(routeF    = as.integer(str_extract(variable, "\\d+")),
                beta      = mean,
                trend     = 100 * (exp(mean) - 1),
                trend_lci = 100 * (exp(q5)   - 1),
                trend_uci = 100 * (exp(q95)  - 1))

    alpha_summ <- summ %>%
      filter(str_detect(variable, "^alpha\\[")) %>%
      transmute(routeF        = as.integer(str_extract(variable, "\\d+")),
                alpha         = mean,
                rel_abundance = exp(mean))

    beta_summ %>%
      left_join(alpha_summ, by = "routeF") %>%
      rename_with(~ paste0(., "_", model_label), -routeF)
  }

  route_global    <- extract_route_trends(summ_global, "global")
  route_stratmean <- extract_route_trends(summ_stratmean, "stratmean")

  route_map_info <- route_map %>%
    st_drop_geometry() %>%
    select(route, routeF)

  route_comparison <- route_map_info %>%
    left_join(route_global, by = "routeF") %>%
    left_join(route_stratmean, by = "routeF") %>%
    mutate(
      trend_diff         = trend_stratmean - trend_global,
      rel_abundance_diff = rel_abundance_stratmean - rel_abundance_global
    ) %>%
    arrange(routeF)

  write.csv(route_comparison,
            file.path(out_dir, paste0(sp_f, "_route_comparison.csv")),
            row.names = FALSE)

  trend_cor      <- cor(route_comparison$trend_global, route_comparison$trend_stratmean)
  trend_mad      <- mean(abs(route_comparison$trend_diff))
  trend_max_diff <- max(abs(route_comparison$trend_diff))
  relabund_cor   <- cor(route_comparison$rel_abundance_global, route_comparison$rel_abundance_stratmean)
  relabund_mad   <- mean(abs(route_comparison$rel_abundance_diff))

  cat("    Trend correlation:", round(trend_cor, 4),
      " | mean abs diff:", round(trend_mad, 4), "pct pts/yr\n")

  # Step 6: stratum-level parameter summary ------------------------------------
  strat_lookup <- route_strat_map %>% distinct(strat_name, strat_idx) %>% arrange(strat_idx)

  alpha_strat_summ <- summ_stratmean %>%
    filter(str_detect(variable, "^ALPHA_strat\\[")) %>%
    transmute(strat_idx = as.integer(str_extract(variable, "\\d+")),
              ALPHA_strat_mean = mean, ALPHA_strat_q5 = q5, ALPHA_strat_q95 = q95)

  beta_strat_summ <- summ_stratmean %>%
    filter(str_detect(variable, "^BETA_strat\\[")) %>%
    transmute(strat_idx = as.integer(str_extract(variable, "\\d+")),
              BETA_strat_mean = mean, BETA_strat_q5 = q5, BETA_strat_q95 = q95)

  strata_summary <- strat_lookup %>%
    left_join(alpha_strat_summ, by = "strat_idx") %>%
    left_join(beta_strat_summ, by = "strat_idx")

  write.csv(strata_summary,
            file.path(out_dir, paste0(sp_f, "_strata_summary.csv")),
            row.names = FALSE)

  mu_alpha       <- summ_stratmean %>% filter(variable == "MU_ALPHA")
  mu_beta        <- summ_stratmean %>% filter(variable == "MU_BETA")
  sd_alpha_strat <- summ_stratmean %>% filter(variable == "sd_ALPHA_strat")
  sd_beta_strat  <- summ_stratmean %>% filter(variable == "sd_BETA_strat")
  global_alpha   <- summ_global %>% filter(variable == "ALPHA")
  global_beta    <- summ_global %>% filter(variable == "BETA")

  cat("    sd_ALPHA_strat (among-BCR sd):", round(sd_alpha_strat$mean, 4),
      " [", round(sd_alpha_strat$q5, 4), ",", round(sd_alpha_strat$q95, 4), "]\n")
  cat("    sd_BETA_strat  (among-BCR sd):", round(sd_beta_strat$mean, 4),
      " [", round(sd_beta_strat$q5, 4), ",", round(sd_beta_strat$q95, 4), "]\n")

  # One-row combined summary for this species ----------------------------------
  data.frame(
    species              = species,
    nroutes              = stan_data_base$nroutes,
    nstrata              = nstrata,
    max_rhat_global      = round(max_rhat_global, 4),
    min_ess_global       = round(min_ess_global, 0),
    pct_divergent_global = round(pct_divergent_global, 2),
    max_rhat_stratmean   = round(max_rhat_strat, 4),
    min_ess_stratmean    = round(min_ess_strat, 0),
    pct_divergent_stratmean = round(pct_divergent_stratmean, 2),
    trend_cor            = round(trend_cor, 4),
    trend_mad            = round(trend_mad, 4),
    trend_max_diff       = round(trend_max_diff, 4),
    relabund_cor         = round(relabund_cor, 4),
    relabund_mad         = round(relabund_mad, 4),
    ALPHA_global         = round(global_alpha$mean, 4),
    BETA_global          = round(global_beta$mean, 4),
    MU_ALPHA_stratmean   = round(mu_alpha$mean, 4),
    MU_BETA_stratmean    = round(mu_beta$mean, 4),
    sd_ALPHA_strat_mean  = round(sd_alpha_strat$mean, 4),
    sd_ALPHA_strat_q5    = round(sd_alpha_strat$q5, 4),
    sd_ALPHA_strat_q95   = round(sd_alpha_strat$q95, 4),
    sd_BETA_strat_mean   = round(sd_beta_strat$mean, 4),
    sd_BETA_strat_q5     = round(sd_beta_strat$q5, 4),
    sd_BETA_strat_q95    = round(sd_beta_strat$q95, 4),
    stringsAsFactors     = FALSE
  )
}

# ==============================================================================
# Main loop: run both models for every species in species_list
# ==============================================================================
summary_list <- list()

for (i in seq_along(species_list)) {
  sp   <- species_list[i]
  sp_f <- species_to_f(sp)

  cat("\n================================================================\n")
  cat("[", i, "/", length(species_list), "]", sp, "\n")
  cat("================================================================\n")

  # Skip if this species' comparison already exists
  route_csv <- file.path(out_dir, paste0(sp_f, "_route_comparison.csv"))
  strat_csv <- file.path(out_dir, paste0(sp_f, "_strata_summary.csv"))
  if (file.exists(route_csv) && file.exists(strat_csv) && !force_refit) {
    cat("  Skipping fit (already done) — re-reading saved CSVs for the summary table.\n")
    # Rebuild just the summary row from the saved per-species RDS summaries,
    # so re-running the script without force_refit is cheap.
    summ_global    <- readRDS(file.path(out_dir, paste0(sp_f, "_global_summ.rds")))
    summ_stratmean <- readRDS(file.path(out_dir, paste0(sp_f, "_stratmean_summ.rds")))
    route_comparison <- read.csv(route_csv)
    strata_summary   <- read.csv(strat_csv)

    trend_cor    <- cor(route_comparison$trend_global, route_comparison$trend_stratmean)
    trend_mad    <- mean(abs(route_comparison$trend_diff))
    relabund_cor <- cor(route_comparison$rel_abundance_global, route_comparison$rel_abundance_stratmean)
    relabund_mad <- mean(abs(route_comparison$rel_abundance_diff))

    mu_alpha       <- summ_stratmean %>% filter(variable == "MU_ALPHA")
    mu_beta        <- summ_stratmean %>% filter(variable == "MU_BETA")
    sd_alpha_strat <- summ_stratmean %>% filter(variable == "sd_ALPHA_strat")
    sd_beta_strat  <- summ_stratmean %>% filter(variable == "sd_BETA_strat")
    global_alpha   <- summ_global %>% filter(variable == "ALPHA")
    global_beta    <- summ_global %>% filter(variable == "BETA")

    # Divergence % requires the sampler diagnostics saved alongside the fit.
    # Species fit before this diagnostics-tracking was added won't have this
    # file yet — leave NA rather than fail, and say so.
    diag_file <- file.path(out_dir, paste0(sp_f, "_diagnostics.rds"))
    if (file.exists(diag_file)) {
      diag <- readRDS(diag_file)
      n_transitions <- chains * iter_sampling
      pct_divergent_global    <- round(100 * sum(diag$global$num_divergent)    / n_transitions, 2)
      pct_divergent_stratmean <- round(100 * sum(diag$stratmean$num_divergent) / n_transitions, 2)
    } else {
      cat("    (No saved diagnostics.rds for", sp, "— divergence % left NA;",
          "set force_refit <- TRUE to backfill.)\n")
      pct_divergent_global    <- NA
      pct_divergent_stratmean <- NA
    }

    summary_list[[sp]] <- data.frame(
      species              = sp,
      nroutes              = max(route_comparison$routeF, na.rm = TRUE),
      nstrata              = nrow(strata_summary),
      max_rhat_global      = round(max(summ_global$rhat, na.rm = TRUE), 4),
      min_ess_global       = round(min(summ_global$ess_bulk, na.rm = TRUE), 0),
      pct_divergent_global = pct_divergent_global,
      max_rhat_stratmean   = round(max(summ_stratmean$rhat, na.rm = TRUE), 4),
      min_ess_stratmean    = round(min(summ_stratmean$ess_bulk, na.rm = TRUE), 0),
      pct_divergent_stratmean = pct_divergent_stratmean,
      trend_cor            = round(trend_cor, 4),
      trend_mad            = round(trend_mad, 4),
      trend_max_diff       = round(max(abs(route_comparison$trend_diff)), 4),
      relabund_cor         = round(relabund_cor, 4),
      relabund_mad         = round(relabund_mad, 4),
      ALPHA_global         = round(global_alpha$mean, 4),
      BETA_global          = round(global_beta$mean, 4),
      MU_ALPHA_stratmean   = round(mu_alpha$mean, 4),
      MU_BETA_stratmean    = round(mu_beta$mean, 4),
      sd_ALPHA_strat_mean  = round(sd_alpha_strat$mean, 4),
      sd_ALPHA_strat_q5    = round(sd_alpha_strat$q5, 4),
      sd_ALPHA_strat_q95   = round(sd_alpha_strat$q95, 4),
      sd_BETA_strat_mean   = round(sd_beta_strat$mean, 4),
      sd_BETA_strat_q5     = round(sd_beta_strat$q5, 4),
      sd_BETA_strat_q95    = round(sd_beta_strat$q95, 4),
      stringsAsFactors     = FALSE
    )
    next
  }

  result <- tryCatch(
    fit_one_species(species = sp, strat = strat,
                    firstYear = firstYear, lastYear = lastYear,
                    model_global = model_global, model_stratmean = model_stratmean,
                    chains = chains, iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                    cmdstanr_output_dir = cmdstanr_output_dir, out_dir = out_dir, spp_df = spp_df),
    error = function(e) {
      message("  [ERROR] Failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(result)) next
  summary_list[[sp]] <- result
}

# ==============================================================================
# Combined species table + automated supplemental-materials statement
# ==============================================================================
if (length(summary_list) == 0) stop("No species fit successfully — nothing to summarize.")

species_summary <- bind_rows(summary_list)
species_summary_csv <- file.path(out_dir, "species_stratum_vs_global_summary.csv")
write.csv(species_summary, species_summary_csv, row.names = FALSE)

cat("\n\n=== Combined per-species summary ===\n")
print(species_summary %>%
        select(species, nroutes, nstrata, trend_cor, trend_mad,
               sd_ALPHA_strat_mean, sd_ALPHA_strat_q5,
               sd_BETA_strat_mean, sd_BETA_strat_q5,
               pct_divergent_global, pct_divergent_stratmean))

# A stratum's mean is treated as "credibly different" from the grand mean for
# a species if the among-stratum sd's 5th percentile clears sd_strat_threshold
# (a posterior concentrated near 0 means BCR strata aren't distinguishable).
species_summary <- species_summary %>%
  mutate(
    alpha_strat_credible = sd_ALPHA_strat_q5 > sd_strat_threshold,
    beta_strat_credible  = sd_BETA_strat_q5  > sd_strat_threshold,
    trends_essentially_identical = trend_cor > trend_cor_threshold &
                                    trend_mad < trend_mad_threshold,
    stratmean_had_convergence_issue = !is.na(pct_divergent_stratmean) &
                                       pct_divergent_stratmean > divergence_pct_threshold
  )

n_species          <- nrow(species_summary)
n_either_credible  <- sum(species_summary$alpha_strat_credible | species_summary$beta_strat_credible)
n_trends_identical <- sum(species_summary$trends_essentially_identical)
n_convergence_issue <- sum(species_summary$stratmean_had_convergence_issue)

median_trend_cor <- median(species_summary$trend_cor)
median_trend_mad <- median(species_summary$trend_mad)
range_trend_cor  <- range(species_summary$trend_cor)
range_trend_mad  <- range(species_summary$trend_mad)

median_pct_div_global    <- median(species_summary$pct_divergent_global, na.rm = TRUE)
median_pct_div_stratmean <- median(species_summary$pct_divergent_stratmean, na.rm = TRUE)
max_pct_div_stratmean    <- max(species_summary$pct_divergent_stratmean, na.rm = TRUE)
worst_divergence_species <- species_summary$species[which.max(species_summary$pct_divergent_stratmean)]
n_worse_convergence <- sum(species_summary$pct_divergent_stratmean > species_summary$pct_divergent_global,
                           na.rm = TRUE)
n_divergence_comparable <- sum(!is.na(species_summary$pct_divergent_global) &
                               !is.na(species_summary$pct_divergent_stratmean))

prefer_stratmean <- n_either_credible > n_species / 2

# Phrase "N of M species" naturally instead of e.g. "only 6 of 6 species",
# which reads oddly once N == M (or N == 0).
frac_word <- function(n, total) {
  if (total == 0) return("no")
  f <- n / total
  if (f == 1)        "all"
  else if (f == 0)    "none of"
  else if (f >= 0.5)  "most"
  else                "a minority of"
}

main_finding <- paste0(
  "Across ", n_species, " test species, allowing the route-level intercept ",
  "and/or slope to vary by BCR stratum (rather than sharing one global mean) ",
  "was credibly supported (BCR-level standard deviation exceeding ", sd_strat_threshold,
  " on the log scale at its 5th percentile, for at least one of ALPHA or BETA) in ",
  frac_word(n_either_credible, n_species), " species (", n_either_credible, " of ", n_species,
  "), while route-level trend estimates from the two parameterizations agreed closely ",
  "(correlation > ", trend_cor_threshold, ", mean absolute difference < ", trend_mad_threshold,
  " percentage points/year) in ", frac_word(n_trends_identical, n_species), " species (",
  n_trends_identical, " of ", n_species, "; median correlation = ", round(median_trend_cor, 3),
  ", range ", round(range_trend_cor[1], 3), "-", round(range_trend_cor[2], 3),
  "; median mean absolute difference = ", round(median_trend_mad, 2),
  " percentage points/year, range ", round(range_trend_mad[1], 2), "-",
  round(range_trend_mad[2], 2), ")."
)

recommendation <- if (prefer_stratmean) {
  paste0(
    "Because BCR strata differ systematically in baseline abundance and/or trend ",
    "for a majority of species tested even though route-level trends themselves ",
    "shift little, we used the stratum-mean parameterization ",
    "(models/slope_iCAR_route_NB_stratmean.stan) for route-level trend estimation. ",
    "It nests the single-global-mean model as a special case and avoids forcing ",
    "routes in ecologically distinct BCRs to share one population-level mean, even ",
    "though its practical effect on route-level trends was small for these species."
  )
} else {
  paste0(
    "Because BCR strata were not, in general, distinguishable from one grand mean, ",
    "we retained the simpler single-global-mean parameterization ",
    "(models/slope_iCAR_route_NB_New.stan) for route-level trend estimation, since ",
    "the added stratum-level parameters of models/slope_iCAR_route_NB_stratmean.stan ",
    "did not meaningfully change route-level trend or relative-abundance estimates."
  )
}

# Convergence caveat: the stratum-mean model has more parameters than the
# global-mean model and can be harder to sample. That matters because
# divergences undermine trust in exactly the sd_ALPHA_strat/sd_BETA_strat
# estimates used above to argue strata differ — so this is computed and
# reported regardless of which model is ultimately recommended, rather than
# assumed.
convergence_caveat <- paste0(
  "The stratum-mean model had a higher divergent-transition rate than the ",
  "global-mean model in ", frac_word(n_worse_convergence, n_divergence_comparable),
  " species (", n_worse_convergence, " of ", n_divergence_comparable,
  " with both models' diagnostics available; median ", round(median_pct_div_stratmean, 2),
  "% divergent transitions for the stratum-mean model vs. ", round(median_pct_div_global, 2),
  "% for the global-mean model), reaching ", round(max_pct_div_stratmean, 2),
  "% divergent transitions for ", worst_divergence_species, ". ",
  if (n_convergence_issue > 0) {
    paste0(
      n_convergence_issue, " of ", n_species, " species exceeded the ",
      divergence_pct_threshold, "% divergent-transition threshold used here to flag ",
      "a convergence problem. Divergences bias posterior exploration and can inflate ",
      "or deflate sd_ALPHA_strat/sd_BETA_strat, so the among-stratum-variation ",
      "evidence above should be treated as provisional until the stratum-mean model's ",
      "sampling is improved (e.g., a non-centered parameterization for ",
      "ALPHA_strat/BETA_strat, or a higher adapt_delta) and re-checked."
    )
  } else {
    paste0(
      "None of the ", n_species, " species exceeded the ", divergence_pct_threshold,
      "% divergent-transition threshold used here to flag a convergence problem, so ",
      "the sd_ALPHA_strat/sd_BETA_strat estimates above are reasonably trustworthy ",
      "as reported."
    )
  }
)

statement <- paste(main_finding, recommendation, convergence_caveat)

statement_file <- file.path(out_dir, "supplemental_statement.txt")
writeLines(c(main_finding, "", recommendation, "", convergence_caveat), statement_file)

cat("\n=== Supplemental-materials statement (auto-generated) ===\n\n")
cat(strwrap(main_finding, width = 80), sep = "\n")
cat("\n")
cat(strwrap(recommendation, width = 80), sep = "\n")
cat("\n")
cat(strwrap(convergence_caveat, width = 80), sep = "\n")
cat("\n")

cat("\n=== Summary ===\n")
cat("Per-species RDS/CSV outputs:      ", out_dir, "\n")
cat("Combined species table:           ", species_summary_csv, "\n")
cat("Supplemental statement:           ", statement_file, "\n")
cat("\nHow to read the combined table:\n")
cat("  - sd_ALPHA_strat / sd_BETA_strat q5 close to 0 => that species' BCR\n")
cat("    strata are NOT distinguishable from one grand mean.\n")
cat("  - trend_cor close to 1 and trend_mad small => switching parameterizations\n")
cat("    barely changes that species' route-level trends.\n")
cat("  - pct_divergent_stratmean noticeably above pct_divergent_global => the\n")
cat("    stratum-mean model is harder to sample for that species; treat its\n")
cat("    sd_ALPHA_strat/sd_BETA_strat estimates cautiously until re-checked with\n")
cat("    e.g. a higher adapt_delta or a non-centered parameterization.\n")
cat("  - The auto-generated statement above pools this evidence across species\n")
cat("    into a few sentences suitable for a supplement; edit thresholds near\n")
cat("    the top of this script (sd_strat_threshold, trend_cor_threshold,\n")
cat("    trend_mad_threshold, divergence_pct_threshold) to match your paper's\n")
cat("    standard of 'meaningful'.\n")
