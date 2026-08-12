## =============================================================================
## Multi-species version of test_grass_to_anthro_scaling.R: does the choice of
## grass_to_anthro transformation (raw_rescaled / standardized / sqrt_rescaled /
## within_between) hold up across several species, not just one?
##
## Background: single-species runs of test_grass_to_anthro_scaling.R showed
## inconsistent results across species —
##   - standardized collapsed gamma1 toward an artificially tight, near-zero
##     interval in every species tried (a prior-scale artifact, not a real
##     finding — see the sd/max ratio discussion for this covariate)
##   - sqrt_rescaled sharpened the estimate for Bobolink but did nothing for
##     Burrowing Owl
##   - within_between was the only variant to find a credible (CI excluding
##     zero) effect at all — the BETWEEN-route term for Burrowing Owl — while
##     for Bobolink it was the WITHIN-route term that carried what signal
##     existed
##   - production settled on raw_rescaled (divide by the covariate's own
##     observed max) specifically for COMPARABILITY with the grassland/anthro
##     models' gamma1, not because it's the most sensitive transform
##
## This script automates that same 4-variant comparison across a small panel
## of species (grassland specialists + a habitat generalist), so the pattern
## seen in 2 species can be checked against a few more before treating it as
## general. It does NOT change the production covariate (still raw_rescaled
## in 1c_species_iCAR_covariates.R) — this is a side analysis only.
##
## For each species, fits all four variants of grass_to_anthro (single
## covariate; the unused second covariate slot is zeroed out for the three
## non-decomposed variants) on the identical reduced dataset, using
## models/slope_iCAR_route_NB_New_2covariates_unbounded.stan, exactly as
## test_grass_to_anthro_scaling.R does for one species — then pools the
## per-species gamma estimates into one combined table.
##
## Outputs (under output/grass_to_anthro_scaling_test/):
##   <species>_<variant>_stanfit.rds / _summ.rds   (per species x variant)
##   <species>_gamma_comparison.csv                (per species, same as
##                                                   single-species script)
##   <species>_route_comparison.csv                (per species, same as
##                                                   single-species script)
##   species_grass_to_anthro_transform_summary.csv (ONE combined table,
##                                                   every species x variant)
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

here::i_am("tests/test_grass_to_anthro_scaling_multi_species.R")
source("functions/neighbours_define_voronoi.R")

# Each species/variant fit gets its own subfolder under here (created/deleted
# inside fit_one_species_g2a()) so raw per-chain CSVs don't pile up in %TEMP%
# across the whole multi-species loop below.
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings --------------------------------------------------------------------
# Species already run individually (Bobolink, Burrowing Owl), plus a few more
# grassland species for broader coverage. Edit freely.
species_list <- c(
  "Baird's Sparrow",
  "Grasshopper Sparrow",
  "Bobolink",
  "Chestnut-collared Longspur",
  "Burrowing Owl"
)

strat      <- "bcr"               # matches 1c_species_iCAR_covariates.R
firstYear  <- 2010
lastYear   <- 2025

# MCMC settings — same as production by default; lower for a quick smoke test
chains         <- 4
iter_warmup    <- 2000
iter_sampling  <- 2000

# Re-fit control: if TRUE, ignore existing per-species output and re-run.
force_refit <- FALSE

out_dir <- here::here("output", "grass_to_anthro_scaling_test")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

species_to_f <- function(sp) gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)

spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

cat("=== grass_to_anthro transformation test (multi-species) ===\n")
cat("Species (n =", length(species_list), "):", paste(species_list, collapse = ", "), "\n")
cat("Strat:", strat, " | Period:", firstYear, "-", lastYear, "\n")

# ----------------------------------------------------------------------------
# grass_to_anthro covariate lookup — same for every species, loaded once
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

# Compile the model once, reused across every species x variant --------------
cat("\nCompiling model...\n")
model_unbounded <- cmdstan_model(
  "models/slope_iCAR_route_NB_New_2covariates_unbounded.stan",
  stanc_options = list("O1"))

# ==============================================================================
# Function: run the full data-prep + 4-variant fit + comparison pipeline for
# ONE species. Returns a data.frame with one row per variant (for the combined
# table) plus writes that species' full RDS/CSV outputs, exactly as the
# single-species version of this script did.
# ==============================================================================
fit_one_species_g2a <- function(species, strat, firstYear, lastYear,
                                g2a_lookup, model_unbounded,
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

  # Step 2: data prep (mirrors prepare_species_data() in 1c, grass_to_anthro only)
  cat("[1/4] Pulling BBS data + joining grass_to_anthro covariate...\n")

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

  n_routes <- max(new_data$routeF)
  n_obs    <- nrow(new_data)
  cat("    Routes:", n_routes, " | Obs:", n_obs, " | Edges:", car_stan_dat$N_edges, "\n")

  if (car_stan_dat$N != n_routes) {
    stop("Some routes are missing from adjacency matrix")
  }

  # Step 3: build the four transformation variants + shared base stan_data ----
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
    nroutes    = n_routes,
    sd_alpha_prior = sd_alpha_prior
  )

  g2a_raw <- new_data$grass_to_anthro_raw

  g2a_max      <- max(g2a_raw)
  g2a_rescaled <- g2a_raw / g2a_max

  g2a_mean <- mean(g2a_raw)
  g2a_sd   <- sd(g2a_raw)
  g2a_z    <- (g2a_raw - g2a_mean) / g2a_sd

  g2a_sqrt          <- sqrt(g2a_raw)
  g2a_sqrt_rescaled <- g2a_sqrt / max(g2a_sqrt)

  new_data$grass_to_anthro_rescaled <- g2a_rescaled
  route_means <- new_data %>%
    group_by(routeF) %>%
    summarise(g2a_between = mean(grass_to_anthro_rescaled))
  new_data <- new_data %>% left_join(route_means, by = "routeF")
  new_data$g2a_within <- new_data$grass_to_anthro_rescaled - new_data$g2a_between

  zeros <- rep(0, nrow(new_data))

  variants <- list(
    raw_rescaled   = c(base_stan_data, list(covariate1 = g2a_rescaled, covariate2 = zeros)),
    standardized   = c(base_stan_data, list(covariate1 = g2a_z, covariate2 = zeros)),
    sqrt_rescaled  = c(base_stan_data, list(covariate1 = g2a_sqrt_rescaled, covariate2 = zeros)),
    within_between = c(base_stan_data, list(covariate1 = new_data$g2a_between,
                                            covariate2 = new_data$g2a_within))
  )

  # Step 4: fit all four variants ----------------------------------------------
  diag_rows <- list()
  summs <- list()

  for (variant_name in names(variants)) {
    cat("  [2/4] Fitting variant:", variant_name, "...\n")

    # Per-species/variant output subfolder so raw cmdstan CSVs stay isolated
    # and can be deleted right after this variant's fit + diagnostics are
    # read (important here since this script loops over many species).
    variant_output_dir <- file.path(cmdstanr_output_dir, paste0(sp_f, "_", variant_name))
    if (!dir.exists(variant_output_dir)) dir.create(variant_output_dir, recursive = TRUE)

    fit <- model_unbounded$sample(
      data = variants[[variant_name]],
      refresh = 0,
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
      species      = species,
      n_obs        = n_obs,
      n_routes     = n_routes,
      variant      = variant_name,
      gamma1_mean  = round(gamma1_row$mean, 4),
      gamma1_q5    = round(gamma1_row$q5, 4),
      gamma1_q95   = round(gamma1_row$q95, 4),
      gamma1_excludes_zero = (gamma1_row$q5 > 0) || (gamma1_row$q95 < 0),
      gamma2_mean  = round(gamma2_row$mean, 4),
      gamma2_q5    = round(gamma2_row$q5, 4),
      gamma2_q95   = round(gamma2_row$q95, 4),
      gamma2_excludes_zero = (gamma2_row$q5 > 0) || (gamma2_row$q95 < 0),
      max_rhat     = round(max(summ$rhat, na.rm = TRUE), 4),
      min_ess      = round(min(summ$ess_bulk, na.rm = TRUE), 0),
      n_divergent  = n_divergent
    )

    summs[[variant_name]] <- summ

    cat("    gamma1:", round(gamma1_row$mean, 4),
        " [", round(gamma1_row$q5, 4), ",", round(gamma1_row$q95, 4), "]",
        " | gamma2:", round(gamma2_row$mean, 4),
        " [", round(gamma2_row$q5, 4), ",", round(gamma2_row$q95, 4), "]",
        " | Rhat:", round(max(summ$rhat, na.rm = TRUE), 4),
        " | ESS:", round(min(summ$ess_bulk, na.rm = TRUE), 0),
        " | div:", n_divergent, "\n")

    # Raw per-chain CSVs are no longer needed now that the fit is saved and
    # all diagnostics for this variant have been read — delete them instead
    # of letting them accumulate in %TEMP% across the whole species loop.
    unlink(variant_output_dir, recursive = TRUE)
  }

  gamma_comparison <- bind_rows(diag_rows)
  write.csv(gamma_comparison,
            file.path(out_dir, paste0(sp_f, "_gamma_comparison.csv")),
            row.names = FALSE)

  # Step 5: per-route comparison (same structure as single-species script) ----
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

  gamma_comparison
}

# ==============================================================================
# Main loop: run every species, pool the per-variant gamma rows into one table
# ==============================================================================
all_rows <- list()

for (species in species_list) {

  sp_f <- species_to_f(species)
  already_done <- all(vapply(
    c("raw_rescaled", "standardized", "sqrt_rescaled", "within_between"),
    function(v) file.exists(file.path(out_dir, paste0(sp_f, "_", v, "_summ.rds"))),
    logical(1)))

  if (already_done && !force_refit) {
    cat("\nSkipping (already fitted, all 4 variants):", species, "\n")
    gc_file <- file.path(out_dir, paste0(sp_f, "_gamma_comparison.csv"))
    if (file.exists(gc_file)) {
      all_rows[[species]] <- read.csv(gc_file, stringsAsFactors = FALSE)
    }
    next
  }

  result <- tryCatch(
    fit_one_species_g2a(species = species, strat = strat,
                        firstYear = firstYear, lastYear = lastYear,
                        g2a_lookup = g2a_lookup, model_unbounded = model_unbounded,
                        chains = chains, iter_warmup = iter_warmup,
                        iter_sampling = iter_sampling,
                        cmdstanr_output_dir = cmdstanr_output_dir,
                        out_dir = out_dir, spp_df = spp_df),
    error = function(e) {
      message("[ERROR] Failed for ", species, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (!is.null(result)) all_rows[[species]] <- result
}

# ==============================================================================
# Combined summary across every species x variant
# ==============================================================================
if (length(all_rows) == 0) {
  stop("No species completed successfully — nothing to summarize.")
}

summary_all <- bind_rows(all_rows)
summary_csv <- here::here("output", "grass_to_anthro_scaling_test",
                          "species_grass_to_anthro_transform_summary.csv")
write.csv(summary_all, summary_csv, row.names = FALSE)

cat("\n\n=== Combined summary across", length(all_rows), "species ===\n")
print(summary_all %>%
        select(species, variant, gamma1_mean, gamma1_q5, gamma1_q95, gamma1_excludes_zero,
               gamma2_mean, gamma2_excludes_zero, n_divergent))

cat("\nCombined CSV written to:", summary_csv, "\n")

# Quick tally: how often did each variant find a credible (non-null) effect? --
cat("\n--- How often did each variant's gamma1 exclude zero? ---\n")
tally1 <- summary_all %>%
  group_by(variant) %>%
  summarise(n_species = n(),
            n_gamma1_excludes_zero = sum(gamma1_excludes_zero),
            .groups = "drop")
print(tally1)

cat("\n--- within_between only: how often did gamma2 (within-route) exclude zero? ---\n")
tally2 <- summary_all %>%
  filter(variant == "within_between") %>%
  select(species, gamma1_excludes_zero, gamma2_excludes_zero)
print(tally2)

cat("\nHow to read this: if 'standardized' rarely/never excludes zero across\n")
cat("species while still having the tightest CIs, that's the prior-scale\n")
cat("artifact discussed earlier, not evidence of no effect. If\n")
cat("'within_between' finds credible effects (in either gamma1 or gamma2)\n")
cat("more often than the pooled variants (raw_rescaled/sqrt_rescaled), that\n")
cat("supports reporting the decomposition as a side analysis even though\n")
cat("raw_rescaled remains the production choice for cross-model comparability.\n")
