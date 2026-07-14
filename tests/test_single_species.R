## =============================================================================
## Smoke test: run 1c_species_iCAR_covariates.R and
## 2c_generate_route_trend_csvs_covariates.R end-to-end for ONE species only,
## so the covariate pipeline can be validated quickly instead of waiting on
## the whole grassland species list.
##
## Species: Baird's Sparrow (change species_filter below to test another).
##
## STEP 0 (below) builds and saves the joined BBS + covariate "datafile" for
## this species DIRECTLY in this script, standalone, before touching 1c at
## all. This mirrors 1_baseline_NB_fit.R (a working reference script that
## joins these same grasslands.csv/developed.csv covariate files
## successfully) and is the same logic now embedded in 1c's
## prepare_species_data(). Building it here too means we can inspect and
## confirm the join is correct in isolation — fast, no Stan compile/sample —
## before spending time on the full 4-model Stan fit.
##
## The earlier version of this test hit a real bug: 1c used
## strat = "bbs_usgs", which makes bbsBayes2 rename BBS routes differently
## than the "<StateNum>-<Route>" format the covariate CSVs are keyed by,
## giving only a ~38% join match rate. 1_baseline_NB_fit.R uses
## strat = "bcr" plus a continental-US filter (country_num == 840,
## state_num != 3) and joins cleanly — both fixes are applied below and in
## 1c_species_iCAR_covariates.R.
##
## How the rest of this script works: 1c_species_iCAR_covariates.R and
## 2c_generate_route_trend_csvs_covariates.R both check
## `if (exists("species_filter"))` right after building their target species
## list, and if it's set, restrict target_spp to just those names. Setting
## species_filter here, before sourcing them, is what limits both scripts to
## one species — no other changes to their logic.
##
## Assumes the R working directory is the project root (Route-level_BBS_trends1),
## same assumption the sourced scripts already make via here::i_am().
##
## Note: this still runs the real Stan sampler (4 chains x 2000 warmup + 2000
## sampling per model, 4 models: base, grassland, developed,
## grassland_developed) for the one species — it is a correctness/
## wiring smoke test, not a fast unit test. If you want a quicker syntax-only
## check, temporarily lower chains/iter_warmup/iter_sampling inside
## 1c_species_iCAR_covariates.R's fit_one_covariate_model() before running.
## =============================================================================

library(here)
library(tidyverse)
library(bbsBayes2)

species_filter <- "Baird's Sparrow"
strat          <- "bcr"
firstYear      <- 2010
lastYear       <- 2025

here::i_am("tests/test_single_species.R")

cat("\n############################################################\n")
cat("# STEP 0: Build & verify the joined datafile for", species_filter, "\n")
cat("############################################################\n\n")

# Species lookup (need the bbs_english name bbsBayes2 expects) --------------
spp_df_test <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                        stringsAsFactors = FALSE)
sp_row <- spp_df_test[spp_df_test$Common.Name == species_filter, ][1, ]
sp_bbs_test <- sp_row$bbs_english

# BBS data, continental US only (matches 1_baseline_NB_fit.R) ---------------
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

cat("BBS observations before covariate join:", nrow(new_data_test),
    "(", length(unique(new_data_test$route)), "routes )\n")

# Covariate lookups: "<StateNum>-<Route>" key, matches 1_baseline_NB_fit.R --
load_covariate_test <- function(file, value_col) {
  read.csv(here::here("data", file), stringsAsFactors = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              value     = .data[[value_col]]) %>%
    filter(!is.na(value)) %>%
    distinct(route_key, year, .keep_all = TRUE)
}

grass_lookup_test <- load_covariate_test("grasslands.csv", "grasslands") %>%
  rename(grassland = value)
dev_lookup_test   <- load_covariate_test("developed.csv", "developed") %>%
  rename(developed = value)

cat("Sample new_data_test$route values:      ",
    paste(head(unique(new_data_test$route), 6), collapse = ", "), "\n")
cat("Sample covariate route_key values:      ",
    paste(head(unique(grass_lookup_test$route_key), 6), collapse = ", "), "\n")

# Inner-join both covariates: this IS the "drop all unnecessary
# observations" step — any route-year without both grassland and developed
# values is dropped, since it never joins in.
datafile <- new_data_test %>%
  inner_join(grass_lookup_test, by = c("route" = "route_key", "r_year" = "year")) %>%
  inner_join(dev_lookup_test,   by = c("route" = "route_key", "r_year" = "year"))

match_rate_test <- nrow(datafile) / nrow(new_data_test)
cat("\nAfter inner-joining grassland + developed covariates:\n")
cat("  Observations:", nrow(datafile), "of", nrow(new_data_test),
    " (", round(100 * match_rate_test, 1), "% kept )\n")
cat("  Routes:", length(unique(datafile$route)), "of",
    length(unique(new_data_test$route)), "\n")

if (nrow(datafile) == 0) {
  stop("STEP 0 FAILED: the covariate join still drops everything. Compare ",
       "the sample route values printed above against grasslands.csv / ",
       "developed.csv's StateNum-Route combinations directly before going ",
       "any further.")
}
if (match_rate_test < 0.5) {
  warning("STEP 0: only ", round(100 * match_rate_test, 1), "% of ",
          "observations matched both covariates — proceeding, but this is ",
          "lower than expected for strat = 'bcr'. Inspect the datafile ",
          "before trusting the fit.")
}

# Save the verified datafile so it can be inspected directly ----------------
datafile_path <- here::here("data", "stan_data",
                            paste0(gsub("'", "", gsub(" ", "_", species_filter)),
                                   "_test_datafile.rds"))
saveRDS(datafile, datafile_path)
cat("\nDatafile saved to:", datafile_path, "\n")
cat("STEP 0 complete — join looks correct, proceeding to the full pipeline.\n")

cat("\n############################################################\n")
cat("# TEST: single-species run — ", species_filter, "\n")
cat("############################################################\n\n")

cat(">>> Sourcing 1c_species_iCAR_covariates.R ...\n\n")
source(here::here("1c_species_iCAR_covariates.R"))

cat("\n>>> Sourcing 2c_generate_route_trend_csvs_covariates.R ...\n\n")
source(here::here("2c_generate_route_trend_csvs_covariates.R"))

## ===========================================================================
## Post-run checks
## ===========================================================================
cat("\n############################################################\n")
cat("# TEST CHECKS\n")
cat("############################################################\n\n")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}
sp_f <- species_to_f(species_filter)

checks <- list()

# 1) Per-model stanfit + summary outputs from 1c ----------------------------
for (tag in c("base", "grassland", "developed", "grassland_developed")) {
  out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
  stanfit_f <- file.path(output_dir, paste0(out_base, "_stanfit.rds"))
  summ_f    <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
  stan_data_f <- here::here("data", "stan_data",
                            paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_stan_data.RData"))

  checks[[paste0("1c: ", tag, " stanfit.rds exists")]]   <- file.exists(stanfit_f)
  checks[[paste0("1c: ", tag, " summ_fit.rds exists")]]  <- file.exists(summ_f)
  checks[[paste0("1c: ", tag, " stan_data.RData exists")]] <- file.exists(stan_data_f)
}

# 2) Combined 2c outputs exist and contain this species ---------------------
route_csv <- file.path(here::here("output", "species_routes_covariates"),
                       paste0("all_route_trends_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))
model_csv <- file.path(here::here("output", "species_routes_covariates"),
                       paste0("all_model_level_summary_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))

checks[["2c: per-route combined CSV exists"]] <- file.exists(route_csv)
checks[["2c: model-level combined CSV exists"]] <- file.exists(model_csv)

if (file.exists(route_csv)) {
  route_out <- read.csv(route_csv, stringsAsFactors = FALSE)
  sp_rows <- route_out[route_out$species == species_filter, ]
  checks[["2c: per-route CSV has rows for this species"]] <- nrow(sp_rows) > 0
  checks[["2c: per-route CSV has all 4 models for this species"]] <-
    all(c("base", "grassland", "developed", "grassland_developed") %in% sp_rows$model)
  checks[["2c: base and covariate models have the SAME route count (shared reduced dataset)"]] <-
    length(unique(table(sp_rows$model))) == 1
}

if (file.exists(model_csv)) {
  model_out <- read.csv(model_csv, stringsAsFactors = FALSE)
  sp_model_rows <- model_out[model_out$species == species_filter, ]
  checks[["2c: model-level CSV has rows for this species"]] <- nrow(sp_model_rows) > 0
  checks[["2c: model-level CSV has non-NA gamma1 for grassland/developed/both"]] <-
    all(!is.na(sp_model_rows$gamma1[sp_model_rows$model %in%
                                    c("grassland", "developed", "grassland_developed")]))
  checks[["2c: model-level CSV has NA gamma1/gamma2 for base"]] <-
    all(is.na(sp_model_rows$gamma1[sp_model_rows$model == "base"])) &&
    all(is.na(sp_model_rows$gamma2[sp_model_rows$model == "base"]))
  checks[["2c: model-level CSV has non-NA gamma2 only for grassland_developed"]] <-
    !is.na(sp_model_rows$gamma2[sp_model_rows$model == "grassland_developed"]) &&
    all(is.na(sp_model_rows$gamma2[sp_model_rows$model %in% c("grassland", "developed", "base")]))
}

# Report ---------------------------------------------------------------------
cat("\n--- Results ---\n")
for (nm in names(checks)) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(checks[[nm]])) "PASS" else "FAIL", nm))
}

n_pass <- sum(unlist(checks) == TRUE, na.rm = TRUE)
n_total <- length(checks)
cat("\n", n_pass, "/", n_total, " checks passed\n", sep = "")

if (n_pass == n_total) {
  cat("\nALL CHECKS PASSED for", species_filter, "\n")
} else {
  cat("\nSOME CHECKS FAILED for", species_filter, "— see above.\n")
}
