## =============================================================================
## Test: which transformation of grass_to_anthro (the grassland-to-Anthro
## transition-rate covariate) works best for the single-covariate model?
##
## Background: unlike grassland/anthro (standing-stock proportions that span
## most of [0,1]), grass_to_anthro is a year-over-year TRANSITION rate — the
## fraction of a route's area converting from grassland to anthropogenic
## cover that year. ~31% of route-years are exact zeros, and the rest trail
## off in a thin right tail up to ~0.15 (raw, unscaled). This creates two
## distinct problems for gamma1 ~ normal(0,1) in
## slope_iCAR_route_NB_New_covariate.stan:
##   (1) MAGNITUDE: the covariate's true range is far narrower than the
##       [0,1] swing the prior implicitly expects
##   (2) SHAPE: it's zero-inflated and right-skewed, so a single linear
##       gamma1 is disproportionately influenced by the rare, large
##       conversion events
##
## 1c_species_iCAR_covariates.R currently fixes (1) only, by rescaling
## grass_to_anthro by its own observed max so it spans [0,1] (a pure linear
## stretch — order/ratios between routes unchanged, zeros stay zero). This
## script checks empirically whether further transformation changes the
## fitted effect or route-level trends, comparing FOUR variants of the SAME
## covariate, all fit on the identical reduced dataset for one species:
##
##   "raw_rescaled"    grass_to_anthro / max(grass_to_anthro) — current
##                     production approach (magnitude fix only)
##   "standardized"    z-scored grass_to_anthro: (x - mean(x)) / sd(x) — a
##                     linear reparameterization of the SAME magnitude fix;
##                     included mainly to confirm it's interchangeable with
##                     raw_rescaled (as it was for grassland)
##   "sqrt_rescaled"   sqrt(grass_to_anthro), then rescaled by its own max —
##                     a genuine SHAPE fix: sqrt is concave, so it stretches
##                     apart the near-zero mass (where ~31% of the data sits
##                     exactly, and most of the rest is concentrated) and
##                     compresses the rare large conversion events, reducing
##                     their leverage on gamma1. Zeros stay exactly zero.
##                     (log1p was considered and deliberately excluded: for
##                     values this far below 1, log1p(x) ~= x, so it would be
##                     a near no-op here — not a useful transform at this
##                     scale.)
##   "within_between"  grass_to_anthro (raw_rescaled) decomposed per route
##                     into a between-route mean and a within-route
##                     deviation. Unlike grassland (almost entirely a
##                     between-route signal for Baird's Sparrow — within-
##                     route range was only +/-0.04-0.06), grass_to_anthro is
##                     a flow/transition variable that can genuinely vary
##                     within a route across years, so this decomposition may
##                     actually separate two different signals here, rather
##                     than being redundant the way it was for grassland.
##
## All four are fit with models/slope_iCAR_route_NB_New_2covariates_unbounded.stan
## (the same unbounded 2-covariate model used in test_covariate_scaling.R),
## with the unused second covariate slot filled with zeros for the three
## single-covariate variants.
##
## This script is self-contained ("run top to bottom" for one species):
## defaults to Baird's Sparrow, matching the other test scripts. Data prep
## mirrors prepare_species_data() in 1c_species_iCAR_covariates.R (restricted
## to grass_to_anthro only, since this test doesn't need grassland/anthro).
##
## Outputs (under output/grass_to_anthro_scaling_test/):
##   <species>_<variant>_stanfit.rds / _summ.rds   (per variant, full posterior)
##   <species>_gamma_comparison.csv    (gamma1/gamma2 mean + 90% CI + Rhat/ESS/
##                                       divergences, one row per variant)
##   <species>_route_comparison.csv    (per-route trend/rel_abundance for all
##                                       four variants, side by side + diffs
##                                       vs. raw_rescaled)
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

here::i_am("tests/test_grass_to_anthro_scaling.R")
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

out_dir <- here::here("output", "grass_to_anthro_scaling_test")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
sp_f <- species_to_f(species)

cat("=== grass_to_anthro transformation test ===\n")
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
# Step 2: grass_to_anthro covariate lookup (RAW, unscaled — transforms are
# applied below, after the join, so every variant starts from the same
# joined rows)
# ----------------------------------------------------------------------------
# check.names = FALSE preserves the column name exactly as it appears in the
# CSV header (needed for "grasslands to Anthro", which contains spaces).
load_covariate <- function(file, value_col) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              value     = .data[[value_col]]) %>%
    filter(!is.na(value)) %>%
    distinct(route_key, year, .keep_all = TRUE)
}
g2a_lookup <- load_covariate("grasslands to Anthro.csv", "grasslands to Anthro") %>%
  rename(grass_to_anthro_raw = value)

# ----------------------------------------------------------------------------
# Step 3: data prep (mirrors prepare_species_data() in 1c, grass_to_anthro only)
# ----------------------------------------------------------------------------
cat("\n[1/4] Pulling BBS data + joining grass_to_anthro covariate...\n")

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
  left_join(g2a_lookup, by = c("route" = "route_key", "r_year" = "year"))

match_rate <- mean(!is.na(new_data$grass_to_anthro_raw))
cat("    grass_to_anthro covariate match rate:", round(100 * match_rate, 1), "%\n")
if (match_rate < 0.5) {
  stop("Less than half of observations matched a grass_to_anthro value for ", species,
       " — check strat ('bcr' expected) and the continental-US filter.")
}

new_data <- new_data %>% filter(!is.na(grass_to_anthro_raw))
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
# Step 4: build the four transformation variants + shared base stan_data
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

g2a_raw <- new_data$grass_to_anthro_raw

# raw_rescaled (current production approach: linear rescale by own max) -----
g2a_max      <- max(g2a_raw)
g2a_rescaled <- g2a_raw / g2a_max

# standardized (z-score of the RAW, unscaled values) -------------------------
g2a_mean <- mean(g2a_raw)
g2a_sd   <- sd(g2a_raw)
g2a_z    <- (g2a_raw - g2a_mean) / g2a_sd

# sqrt_rescaled (shape fix: sqrt, then rescaled by its own max) --------------
g2a_sqrt          <- sqrt(g2a_raw)
g2a_sqrt_rescaled <- g2a_sqrt / max(g2a_sqrt)

# within/between decomposition of raw_rescaled -------------------------------
new_data$grass_to_anthro_rescaled <- g2a_rescaled
route_means <- new_data %>%
  group_by(routeF) %>%
  summarise(g2a_between = mean(grass_to_anthro_rescaled))
new_data <- new_data %>% left_join(route_means, by = "routeF")
new_data$g2a_within <- new_data$grass_to_anthro_rescaled - new_data$g2a_between

cat("\n  Covariate summary (grass_to_anthro, n =", nrow(new_data), "route-years):\n")
cat("    raw:            range [", round(min(g2a_raw), 6), ",", round(max(g2a_raw), 6),
    "], mean", round(g2a_mean, 6), ", sd", round(g2a_sd, 6),
    ", frac zero:", round(mean(g2a_raw == 0), 3), "\n")
cat("    raw_rescaled:   range [", round(min(g2a_rescaled), 3), ",", round(max(g2a_rescaled), 3), "]\n")
cat("    sqrt_rescaled:  range [", round(min(g2a_sqrt_rescaled), 3), ",", round(max(g2a_sqrt_rescaled), 3), "]\n")
cat("    between-route:  range [", round(min(new_data$g2a_between), 3), ",",
    round(max(new_data$g2a_between), 3), "]\n")
cat("    within-route:   range [", round(min(new_data$g2a_within), 3), ",",
    round(max(new_data$g2a_within), 3), "]\n")

zeros <- rep(0, nrow(new_data))

stan_data_raw_rescaled <- c(base_stan_data,
                            list(covariate1 = g2a_rescaled, covariate2 = zeros))
stan_data_standardized  <- c(base_stan_data,
                            list(covariate1 = g2a_z, covariate2 = zeros))
stan_data_sqrt_rescaled <- c(base_stan_data,
                            list(covariate1 = g2a_sqrt_rescaled, covariate2 = zeros))
stan_data_within_between <- c(base_stan_data,
                            list(covariate1 = new_data$g2a_between,
                                 covariate2 = new_data$g2a_within))

variants <- list(raw_rescaled    = stan_data_raw_rescaled,
                 standardized    = stan_data_standardized,
                 sqrt_rescaled   = stan_data_sqrt_rescaled,
                 within_between  = stan_data_within_between)

# ----------------------------------------------------------------------------
# Step 5: compile once, fit all four variants on the identical data
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
  # be deleted right after this variant's fit + diagnostics are read.
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
    trend_diff_std_vs_raw  = trend_standardized   - trend_raw_rescaled,
    trend_diff_sqrt_vs_raw = trend_sqrt_rescaled   - trend_raw_rescaled,
    trend_diff_wb_vs_raw   = trend_within_between  - trend_raw_rescaled
  ) %>%
  arrange(routeF)

write.csv(route_comparison,
          file.path(out_dir, paste0(sp_f, "_route_comparison.csv")),
          row.names = FALSE)

cor_std  <- cor(route_comparison$trend_raw_rescaled, route_comparison$trend_standardized)
mad_std  <- mean(abs(route_comparison$trend_diff_std_vs_raw))
cor_sqrt <- cor(route_comparison$trend_raw_rescaled, route_comparison$trend_sqrt_rescaled)
mad_sqrt <- mean(abs(route_comparison$trend_diff_sqrt_vs_raw))
cor_wb   <- cor(route_comparison$trend_raw_rescaled, route_comparison$trend_within_between)
mad_wb   <- mean(abs(route_comparison$trend_diff_wb_vs_raw))

cat("\n  Route-level trend agreement vs. raw_rescaled (n =", nrow(route_comparison), "routes):\n")
cat("    standardized   — correlation:", round(cor_std, 4),  " | mean abs diff:", round(mad_std, 4), "pp\n")
cat("    sqrt_rescaled  — correlation:", round(cor_sqrt, 4), " | mean abs diff:", round(mad_sqrt, 4), "pp\n")
cat("    within_between — correlation:", round(cor_wb, 4),   " | mean abs diff:", round(mad_wb, 4), "pp\n")

cat("\n  gamma comparison across variants:\n")
print(gamma_comparison)

cat("\n=== Summary ===\n")
cat("Gamma comparison CSV:", file.path(out_dir, paste0(sp_f, "_gamma_comparison.csv")), "\n")
cat("Route comparison CSV:", file.path(out_dir, paste0(sp_f, "_route_comparison.csv")), "\n")
cat("\nHow to read this:\n")
cat("  - raw_rescaled vs. standardized should be near-identical (both are\n")
cat("    pure linear rescales of the same magnitude fix) — a real difference\n")
cat("    here would be surprising and worth double-checking.\n")
cat("  - raw_rescaled vs. sqrt_rescaled is the SHAPE check: if gamma1's\n")
cat("    credible interval narrows noticeably and/or its point estimate\n")
cat("    shifts under sqrt_rescaled, the raw linear treatment was likely\n")
cat("    being dominated by a handful of large, rare conversion events —\n")
cat("    sqrt is the better-specified transform for this covariate. If\n")
cat("    gamma1 and the route-level trends barely move, the skew isn't\n")
cat("    materially affecting the fit and the extra transform isn't needed.\n")
cat("  - raw_rescaled vs. within_between is the STRUCTURE check: unlike\n")
cat("    grassland (almost entirely a between-route signal for this\n")
cat("    species), grass_to_anthro can genuinely vary within a route across\n")
cat("    years. If the between/within gammas diverge from each other and\n")
cat("    from the pooled raw_rescaled estimate, the raw covariate was\n")
cat("    blending two different signals (which route has more conversion,\n")
cat("    on average, vs. did this route see MORE conversion than usual in a\n")
cat("    given year) that are worth reporting separately.\n")
