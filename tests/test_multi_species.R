## =============================================================================
## Multi-species version of test_single_species.R: run
## 1c_species_iCAR_covariates.R and 2c_generate_route_trend_csvs_covariates.R
## end-to-end for SEVERAL species at once (all FOUR covariate models — base,
## grassland_habitat, grassland_anthro, grassland_habitat_to_anthro — per
## species), so the production pipeline can be validated across more than one
## species without running the full grassland species list.
##
## Species: set via species_filter below (a vector now, not a single string).
##
## Structurally identical to test_single_species.R, just pluralized:
##   - STEP 0 (below) builds and saves the joined BBS + covariate "datafile"
##     for EACH species in species_filter, standalone, before touching 1c —
##     same join logic as prepare_species_data() in 1c ("grassland
##     habitat.csv" + "grassland anthro.csv" + "grassland habitat to
##     anthro.csv", strat = "bcr", continental-US filter), just looped per
##     species so a bad join for any one species is caught before spending
##     time on the Stan fits.
##   - 1c_species_iCAR_covariates.R and 2c_generate_route_trend_csvs_covariates.R
##     are each sourced ONCE (not once per species) — both already support
##     `species_filter` as a vector (they do
##     `target_spp %>% filter(Common.Name %in% species_filter)` right after
##     building their target species list), so a single source() call fits
##     all species in the list.
##   - Post-run checks loop over species_filter, running the same per-species
##     checks test_single_species.R does, then report pass/fail per species
##     plus an overall tally.
##
## Note: this runs the real Stan sampler (4 chains x 2000 warmup + 2000
## sampling per model, 4 models, for EVERY species in species_filter) — for
## N species that's 4N model fits. This is a correctness/wiring smoke test
## across species, not a fast unit test; keep species_filter short for a
## reasonable runtime, or temporarily lower chains/iter_warmup/iter_sampling
## inside 1c_species_iCAR_covariates.R's fit_one_covariate_model().
## =============================================================================

library(here)
library(tidyverse)
library(bbsBayes2)

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
strat          <- "bcr"
firstYear      <- 2010
lastYear       <- 2025

here::i_am("tests/test_multi_species.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

cat("\n############################################################\n")
cat("# STEP 0: Build & verify the joined datafile for each species\n")
cat("# Species (n =", length(species_filter), "):", paste(species_filter, collapse = ", "), "\n")
cat("############################################################\n\n")

spp_df_test <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                        stringsAsFactors = FALSE)

# Covariate lookups: "<StateNum>-<Route>" key, matches 1_baseline_NB_fit.R.
# Loaded ONCE, shared across every species' join below (mirrors 1c, which
# also builds these once outside its per-species loop).
# check.names = FALSE preserves value_col exactly as it appears in the CSV
# header (needed since all three covariate names contain spaces).
load_covariate_test <- function(file, value_col) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              value     = .data[[value_col]]) %>%
    filter(!is.na(value)) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

grass_habitat_lookup_test <- load_covariate_test("grassland habitat.csv", "grassland habitat") %>%
  rename(grassland_habitat = value)
grass_anthro_lookup_test  <- load_covariate_test("grassland anthro.csv", "grassland anthro") %>%
  rename(grassland_anthro = value)
ghta_lookup_test          <- load_covariate_test("grassland habitat to anthro.csv", "grassland habitat to anthro") %>%
  rename(grassland_habitat_to_anthro = value)

# grassland habitat to anthro.csv is a year-over-year transition rate (mostly
# zero/near-zero), so it only spans [0, ~0.15] raw — rescale by its own
# observed max (matches 1c_species_iCAR_covariates.R) so it spans [0,1] like
# grassland_habitat/grassland_anthro do, keeping gamma1 ~ normal(0,1)
# sensible for it.
ghta_scale_test <- max(ghta_lookup_test$grassland_habitat_to_anthro, na.rm = TRUE)
ghta_lookup_test <- ghta_lookup_test %>%
  mutate(grassland_habitat_to_anthro = grassland_habitat_to_anthro / ghta_scale_test)
cat("Grassland-habitat-to-Anthro covariate rescaled by its observed max (",
    round(ghta_scale_test, 6), ") to span [0,1]\n\n", sep = "")

step0_summary <- list()

for (species_filter_sp in species_filter) {

  cat("--- STEP 0:", species_filter_sp, "---\n")

  sp_row <- spp_df_test[spp_df_test$Common.Name == species_filter_sp, ][1, ]
  if (is.na(sp_row$Common.Name)) {
    stop("Species '", species_filter_sp, "' not found in spp_names_codes_group_aou.csv")
  }
  sp_bbs_test <- sp_row$bbs_english

  data_pkg_test <- bbsBayes2::stratify(by = strat, species = sp_bbs_test,
                                       use_map = FALSE) %>%
    bbsBayes2::prepare_data(min_year = firstYear, max_year = lastYear,
                            min_n_routes = 1, min_max_route_years = 1)

  raw_data_test <- data_pkg_test[["raw_data"]] %>%
    filter(country_num == 840) %>%
    filter(state_num != 3)

  new_data_test <- data.frame(
    strat_name = raw_data_test$strata_name,
    strat      = raw_data_test$strata,
    route      = raw_data_test$route,
    latitude   = raw_data_test$latitude,
    longitude  = raw_data_test$longitude,
    count      = raw_data_test$count,
    year       = raw_data_test$year_num,
    firstyr    = raw_data_test$first_year,
    ObsN       = raw_data_test$observer,
    r_year     = raw_data_test$year
  )

  cat("  BBS observations before covariate join:", nrow(new_data_test),
      "(", length(unique(new_data_test$route)), "routes )\n")

  # Inner-join all three covariates: any route-year without grassland_habitat,
  # grassland_anthro, AND grassland_habitat_to_anthro values is dropped, since
  # it never joins in.
  datafile <- new_data_test %>%
    inner_join(grass_habitat_lookup_test, by = c("route" = "route_key", "r_year" = "year")) %>%
    inner_join(grass_anthro_lookup_test,  by = c("route" = "route_key", "r_year" = "year")) %>%
    inner_join(ghta_lookup_test,          by = c("route" = "route_key", "r_year" = "year"))

  match_rate_test <- nrow(datafile) / nrow(new_data_test)
  cat("  After inner-joining grassland_habitat + grassland_anthro + grassland_habitat_to_anthro covariates:\n")
  cat("    Observations:", nrow(datafile), "of", nrow(new_data_test),
      " (", round(100 * match_rate_test, 1), "% kept )\n")
  cat("    Routes:", length(unique(datafile$route)), "of",
      length(unique(new_data_test$route)), "\n")

  step0_ok <- TRUE
  if (nrow(datafile) == 0) {
    step0_ok <- FALSE
    message("  STEP 0 FAILED for ", species_filter_sp, ": the covariate join drops everything.")
  } else if (match_rate_test < 0.5) {
    step0_ok <- FALSE
    warning("  STEP 0: only ", round(100 * match_rate_test, 1), "% of observations for ",
            species_filter_sp, " matched all three covariates — lower than expected for strat = 'bcr'.")
  }

  datafile_path <- here::here("data", "stan_data",
                              paste0(species_to_f(species_filter_sp), "_test_datafile.rds"))
  saveRDS(datafile, datafile_path)
  cat("  Datafile saved to:", datafile_path, "\n\n")

  step0_summary[[species_filter_sp]] <- data.frame(
    species = species_filter_sp,
    n_obs_before_join = nrow(new_data_test),
    n_obs_after_join  = nrow(datafile),
    match_rate        = round(100 * match_rate_test, 1),
    step0_ok          = step0_ok
  )
}

step0_summary_df <- bind_rows(step0_summary)
cat("--- STEP 0 summary, all species ---\n")
print(step0_summary_df)

if (any(!step0_summary_df$step0_ok)) {
  stop("STEP 0 failed for at least one species (see above) — fix the join before running 1c/2c.")
}
cat("\nSTEP 0 complete for all species — proceeding to the full pipeline.\n")

cat("\n############################################################\n")
cat("# TEST: multi-species run —", length(species_filter), "species\n")
cat("############################################################\n\n")

cat(">>> Sourcing 1c_species_iCAR_covariates.R (all species in species_filter) ...\n\n")
source(here::here("1c_species_iCAR_covariates.R"))

cat("\n>>> Sourcing 2c_generate_route_trend_csvs_covariates.R (all species in species_filter) ...\n\n")
source(here::here("2c_generate_route_trend_csvs_covariates.R"))

## ===========================================================================
## Post-run checks — same checks as test_single_species.R, looped per species
## ===========================================================================
cat("\n############################################################\n")
cat("# TEST CHECKS\n")
cat("############################################################\n\n")

model_tags_check <- c("base", "grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro")

route_csv <- file.path(here::here("output", "species_routes_covariates"),
                       paste0("all_route_trends_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))
model_csv <- file.path(here::here("output", "species_routes_covariates"),
                       paste0("all_model_level_summary_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))

route_out <- if (file.exists(route_csv)) read.csv(route_csv, stringsAsFactors = FALSE) else NULL
model_out <- if (file.exists(model_csv)) read.csv(model_csv, stringsAsFactors = FALSE) else NULL

all_checks <- list()

for (sp in species_filter) {
  sp_f <- species_to_f(sp)
  checks <- list()

  # 1) Per-model stanfit + summary outputs from 1c ---------------------------
  for (tag in model_tags_check) {
    out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    stanfit_f <- file.path(output_dir, paste0(out_base, "_stanfit.rds"))
    summ_f    <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
    stan_data_f <- here::here("data", "stan_data",
                              paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_stan_data.RData"))

    checks[[paste0("1c: ", tag, " stanfit.rds exists")]]     <- file.exists(stanfit_f)
    checks[[paste0("1c: ", tag, " summ_fit.rds exists")]]    <- file.exists(summ_f)
    checks[[paste0("1c: ", tag, " stan_data.RData exists")]] <- file.exists(stan_data_f)
  }

  # 2) This species' rows in the combined 2c outputs -------------------------
  if (!is.null(route_out)) {
    sp_rows <- route_out[route_out$species == sp, ]
    checks[["2c: per-route CSV has rows for this species"]] <- nrow(sp_rows) > 0
    checks[["2c: per-route CSV has all 4 models for this species"]] <-
      all(model_tags_check %in% sp_rows$model)
    checks[["2c: base and covariate models have the SAME route count (shared reduced dataset)"]] <-
      length(unique(table(sp_rows$model))) == 1
  } else {
    checks[["2c: per-route combined CSV exists"]] <- FALSE
  }

  if (!is.null(model_out)) {
    sp_model_rows <- model_out[model_out$species == sp, ]
    checks[["2c: model-level CSV has rows for this species"]] <- nrow(sp_model_rows) > 0
    checks[["2c: model-level CSV has non-NA gamma1 for grassland_habitat/grassland_anthro/grassland_habitat_to_anthro"]] <-
      all(!is.na(sp_model_rows$gamma1[sp_model_rows$model %in%
                                      c("grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro")]))
    checks[["2c: model-level CSV has NA gamma1 for base"]] <-
      all(is.na(sp_model_rows$gamma1[sp_model_rows$model == "base"]))
  } else {
    checks[["2c: model-level combined CSV exists"]] <- FALSE
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

## ===========================================================================
## Overall summary across every species
## ===========================================================================
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
