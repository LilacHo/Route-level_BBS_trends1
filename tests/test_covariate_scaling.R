## =============================================================================
## Test: does transforming the land-cover covariate change the fitted effect
## or the route-level trends, compared to the raw [0,1] proportion currently
## used in 1c_species_iCAR_covariates.R?
##
## Background: 1c_species_iCAR_covariates.R feeds grassland/developed cover
## straight in as raw proportions (range [0,1]), with gamma1 ~ normal(0,1).
## Three variants of a single covariate (grassland, for simplicity — the same
## logic applies to developed) are compared here, all fit on the IDENTICAL
## reduced dataset for one species:
##
##   "raw"           grassland proportion, unchanged (current approach)
##   "standardized"  z-scored grassland: (x - mean(x)) / sd(x)
##                   a linear reparameterization only — should reproduce the
##                   same fit/predictions as "raw", just rescaling gamma1 and
##                   (hopefully) improving sampler behaviour if the real
##                   spread of grassland values is much narrower than [0,1]
##   "within_between" grassland decomposed per route into a between-route
##                   component (each route's own mean grassland value,
##                   constant over time) and a within-route component (that
##                   route-year's deviation from its own mean). alpha_r/beta_r
##                   are already free per-route parameters, so the raw
##                   covariate's between-route signal is largely redundant
##                   with them; this decomposition isolates the within-route
##                   "did habitat change over time in this route predict a
##                   trend deviation" effect from the between-route
##                   cross-sectional difference, and reduces the collinearity
##                   between gamma and beta_r that a monotonic within-route
##                   habitat trend can otherwise cause.
##
## All three are fit with models/slope_iCAR_route_NB_New_2covariates_unbounded.stan
## (identical to the 2-covariate model but without the [0,1] data bound, since
## standardized/within-route values can be negative). For "raw" and
## "standardized" the second covariate slot is filled with a constant column
## of zeros, so gamma2 has no effect on the likelihood (safe no-op) and the
## fit is equivalent to the single-covariate model.
##
## This script is self-contained ("run top to bottom" for one species): data
## prep mirrors prepare_species_data() in 1c_species_iCAR_covariates.R
## (grassland only), fits all three variants on the identical data, and
## writes comparison CSVs + prints a summary.
##
## Outputs (under output/covariate_scaling_test/):
##   <species>_<variant>_stanfit.rds / _summ.rds   (per variant, full posterior)
##   <species>_gamma_comparison.csv    (gamma1/gamma2 mean + 90% CI + Rhat/ESS/
##                                       divergences, one row per variant)
##   <species>_route_comparison.csv    (per-route trend/rel_abundance for all
##                                       three variants, side by side + diffs)
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

here::i_am("tests/test_covariate_scaling.R")
source("functions/neighbours_define_voronoi.R")

# Each variant fit below gets its own subfolder under here (created/deleted
# in the fitting loop) so raw per-chain CSVs don't pile up in %TEMP%.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings ------------------------------------------------------------------
species    <- "Baird's Sparrow"   # change to test a different species
strat      <- "bcr"               # matches 1c_species_iCAR_covariates.R
firstYear  <- 2010
lastYear   <- 2025

# MCMC settings — same as production by default; lower for a quick smoke test
chains         <- 4
iter_warmup    <- 2000
iter_sampling  <- 2000

out_dir <- here::here("output", "covariate_scaling_test")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
sp_f <- species_to_f(species)

cat("=== Covariate scaling/transformation test ===\n")
cat("Species:", species, " | Strat:", strat, " | Period:", firstYear, "-", lastYear, "\n")

# ----------------------------------------------------------------------------
# Step 1: look up the BBS name bbsBayes2 expects
# ----------------------------------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)
sp_row <- spp_df[spp_df$Common.Name == species, ][1, ]
if (is.na(sp_row$Common.Name)) {
  stop("Species '", species, "' not found in spp_names_codes_group_aou.csv")
}
species_bbs <- sp_row$bbs_english

# ----------------------------------------------------------------------------
# Step 2: grassland covariate lookup (mirrors load_covariate() in 1c)
# ----------------------------------------------------------------------------
load_covariate <- function(file, value_col) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              value     = .data[[value_col]]) %>%
    filter(!is.na(value)) %>%
    distinct(route_key, year, .keep_all = TRUE)
}
grass_lookup <- load_covariate("grasslands.csv", "grasslands") %>%
  rename(grassland = value)

# ----------------------------------------------------------------------------
# Step 3: data prep (mirrors prepare_species_data() in 1c, grassland only)
# ----------------------------------------------------------------------------
cat("\n[1/4] Pulling BBS data + joining grassland covariate...\n")

data_pkg <- bbsBayes2::stratify(by = strat, species = species_bbs,
                                use_map = FALSE) %>%
  bbsBayes2::prepare_data(min_year = firstYear,
                          max_year = lastYear,
                          min_n_routes = 1,
                          min_max_route_years = 1)

raw_data   <- data_pkg[["raw_data"]]
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

new_data <- new_data %>%
  left_join(grass_lookup, by = c("route" = "route_key", "r_year" = "year"))

match_rate <- mean(!is.na(new_data$grassland))
cat("    Grassland covariate match rate:", round(100 * match_rate, 1), "%\n")
if (match_rate < 0.5) {
  stop("Less than half of observations matched a grassland value for ", species,
       " — check strat ('bcr' expected) and the continental-US filter.")
}

new_data <- new_data %>% filter(!is.na(grassland))
if (nrow(new_data) < 50) {
  stop("Too few observations after dropping missing-covariate rows (", nrow(new_data), ")")
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

# ----------------------------------------------------------------------------
# Step 4: build the three covariate variants + shared base stan_data
# ----------------------------------------------------------------------------
tmp <- data.frame(route = new_data$routeF, count = new_data$count) %>%
  group_by(route) %>%
  summarise(mean_count = mean(log(count + 1)))
sd_alpha_prior <- sd(tmp$mean_count)

base_stan_data <- list(
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

# raw ------------------------------------------------------------------------
grass_raw <- new_data$grassland

# standardized (z-score) ------------------------------------------------------
grass_mean <- mean(new_data$grassland)
grass_sd   <- sd(new_data$grassland)
grass_z    <- (new_data$grassland - grass_mean) / grass_sd

# within/between decomposition -------------------------------------------------
route_means <- new_data %>%
  group_by(routeF) %>%
  summarise(grass_between = mean(grassland))
new_data <- new_data %>% left_join(route_means, by = "routeF")
new_data$grass_within <- new_data$grassland - new_data$grass_between

cat("\n  Covariate summary (grassland, n =", nrow(new_data), "route-years):\n")
cat("    raw:            range [", round(min(grass_raw), 3), ",", round(max(grass_raw), 3),
    "], mean", round(grass_mean, 3), ", sd", round(grass_sd, 3), "\n")
cat("    between-route:  range [", round(min(new_data$grass_between), 3), ",",
    round(max(new_data$grass_between), 3), "]\n")
cat("    within-route:   range [", round(min(new_data$grass_within), 3), ",",
    round(max(new_data$grass_within), 3), "]\n")

zeros <- rep(0, nrow(new_data))

stan_data_raw <- c(base_stan_data,
                   list(covariate1 = grass_raw, covariate2 = zeros))
stan_data_std <- c(base_stan_data,
                   list(covariate1 = grass_z, covariate2 = zeros))
stan_data_wb  <- c(base_stan_data,
                   list(covariate1 = new_data$grass_between,
                        covariate2 = new_data$grass_within))

variants <- list(raw = stan_data_raw,
                 standardized = stan_data_std,
                 within_between = stan_data_wb)

# ----------------------------------------------------------------------------
# Step 5: compile once, fit all three variants on the identical data
# ----------------------------------------------------------------------------
cat("\n[2/4] Compiling model...\n")
model_unbounded <- cmdstan_model(
  "models/slope_iCAR_route_NB_New_2covariates_unbounded.stan",
  stanc_options = list("O1"))

fits <- list()
summs <- list()
diag_rows <- list()

for (variant_name in names(variants)) {
  cat("\n[3/4] Fitting variant:", variant_name, "...\n")

  # Per-variant output subfolder so raw cmdstan CSVs stay isolated and can
  # be deleted right after this variant's fit is saved to .rds.
  variant_output_dir <- file.path(cmdstanr_output_dir, paste0(sp_f, "_", variant_name))
  if (!dir.exists(variant_output_dir)) dir.create(variant_output_dir, recursive = TRUE)

  fit <- model_unbounded$sample(
    data = variants[[variant_name]],
    refresh = 400,
    chains = chains, iter_sampling = iter_sampling, iter_warmup = iter_warmup,
    parallel_chains = chains,
    adapt_delta = 0.8,
    max_treedepth = 10,
    show_exceptions = FALSE,
    output_dir = variant_output_dir
  )
  summ <- fit$summary()

  fit$save_object(file.path(out_dir, paste0(sp_f, "_", variant_name, "_stanfit.rds")))
  saveRDS(summ, file.path(out_dir, paste0(sp_f, "_", variant_name, "_summ.rds")))

  n_divergent <- sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)

  gamma1_row <- summ %>% filter(variable == "gamma1_out")
  gamma2_row <- summ %>% filter(variable == "gamma2_out")

  diag_rows[[variant_name]] <- data.frame(
    variant     = variant_name,
    gamma1_mean = round(gamma1_row$mean, 4),
    gamma1_q5   = round(gamma1_row$q5, 4),
    gamma1_q95  = round(gamma1_row$q95, 4),
    gamma2_mean = round(gamma2_row$mean, 4),
    gamma2_q5   = round(gamma2_row$q5, 4),
    gamma2_q95  = round(gamma2_row$q95, 4),
    max_rhat    = round(max(summ$rhat, na.rm = TRUE), 4),
    min_ess     = round(min(summ$ess_bulk, na.rm = TRUE), 0),
    n_divergent = n_divergent
  )

  fits[[variant_name]]  <- fit
  summs[[variant_name]] <- summ

  cat("    gamma1:", round(gamma1_row$mean, 4),
      " [", round(gamma1_row$q5, 4), ",", round(gamma1_row$q95, 4), "]",
      " | gamma2:", round(gamma2_row$mean, 4),
      " [", round(gamma2_row$q5, 4), ",", round(gamma2_row$q95, 4), "]",
      " | max Rhat:", round(max(summ$rhat, na.rm = TRUE), 4),
      " | min ESS:", round(min(summ$ess_bulk, na.rm = TRUE), 0),
      " | divergences:", n_divergent, "\n")

  # Raw per-chain CSVs are no longer needed now that the fit is saved and
  # all diagnostics for this variant have been read — delete them instead
  # of letting them accumulate in %TEMP%.
  unlink(variant_output_dir, recursive = TRUE)
}

gamma_comparison <- bind_rows(diag_rows)
write.csv(gamma_comparison,
          file.path(out_dir, paste0(sp_f, "_gamma_comparison.csv")),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# Step 6: compare route-level trend / relative abundance across variants
# ----------------------------------------------------------------------------
cat("\n[4/4] Comparing route-level estimates across variants...\n")

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

route_map_info <- route_map %>%
  st_drop_geometry() %>%
  select(route, routeF)

route_trend_list <- lapply(names(summs), function(v) extract_route_trends(summs[[v]], v))

route_comparison <- route_map_info
for (rt in route_trend_list) {
  route_comparison <- route_comparison %>% left_join(rt, by = "routeF")
}
route_comparison <- route_comparison %>%
  mutate(
    trend_diff_std_vs_raw = trend_standardized - trend_raw,
    trend_diff_wb_vs_raw  = trend_within_between - trend_raw
  ) %>%
  arrange(routeF)

write.csv(route_comparison,
          file.path(out_dir, paste0(sp_f, "_route_comparison.csv")),
          row.names = FALSE)

cor_std <- cor(route_comparison$trend_raw, route_comparison$trend_standardized)
mad_std <- mean(abs(route_comparison$trend_diff_std_vs_raw))
cor_wb  <- cor(route_comparison$trend_raw, route_comparison$trend_within_between)
mad_wb  <- mean(abs(route_comparison$trend_diff_wb_vs_raw))

cat("\n  Route-level trend agreement (n =", nrow(route_comparison), "routes):\n")
cat("    raw vs. standardized  — correlation:", round(cor_std, 4),
    " | mean abs diff:", round(mad_std, 4), "pp\n")
cat("    raw vs. within/between — correlation:", round(cor_wb, 4),
    " | mean abs diff:", round(mad_wb, 4), "pp\n")

cat("\n  gamma comparison across variants:\n")
print(gamma_comparison)

cat("\n=== Summary ===\n")
cat("Gamma comparison CSV:", file.path(out_dir, paste0(sp_f, "_gamma_comparison.csv")), "\n")
cat("Route comparison CSV:", file.path(out_dir, paste0(sp_f, "_route_comparison.csv")), "\n")
cat("\nHow to read this:\n")
cat("  - raw vs. standardized should give NEAR-IDENTICAL route-level trends\n")
cat("    (it's a linear reparameterization); differences mainly show up in\n")
cat("    gamma1's scale/CI width and in sampler diagnostics (Rhat/ESS/\n")
cat("    divergences) — a real improvement there suggests standardizing\n")
cat("    helps the sampler without changing the substantive answer.\n")
cat("  - raw vs. within/between CAN differ in gamma1 and in route-level\n")
cat("    trends. gamma1 here is now specifically the effect of a route's\n")
cat("    grassland level DEVIATING from its own mean, isolating it from\n")
cat("    the between-route cross-sectional signal that's largely already\n")
cat("    absorbed into alpha_r/beta_r. A notable shift indicates the raw\n")
cat("    covariate was conflating within- and between-route effects, and/or\n")
cat("    that gamma1 and beta_r were competing to explain the same within-\n")
cat("    route temporal trend in the raw parameterization.\n")
