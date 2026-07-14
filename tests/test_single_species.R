## =============================================================================
## Smoke test: run 1c_species_iCAR_covariates.R and
## 2c_generate_route_trend_csvs_covariates.R end-to-end for ONE species only,
## so the covariate pipeline can be validated quickly instead of waiting on
## the whole grassland species list.
##
## Species: Baird's Sparrow (change species_filter below to test another).
##
## How this works: 1c_species_iCAR_covariates.R and
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
## sampling per model, 3 models) for the one species — it is a correctness/
## wiring smoke test, not a fast unit test. If you want a quicker syntax-only
## check, temporarily lower chains/iter_warmup/iter_sampling inside
## 1c_species_iCAR_covariates.R's fit_one_covariate_model() before running.
## =============================================================================

library(here)

species_filter <- "Baird's Sparrow"

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
for (tag in c("grassland", "developed", "grassland_developed")) {
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
  checks[["2c: per-route CSV has all 3 covariate models for this species"]] <-
    all(c("grassland", "developed", "grassland_developed") %in% sp_rows$model)
}

if (file.exists(model_csv)) {
  model_out <- read.csv(model_csv, stringsAsFactors = FALSE)
  sp_model_rows <- model_out[model_out$species == species_filter, ]
  checks[["2c: model-level CSV has rows for this species"]] <- nrow(sp_model_rows) > 0
  checks[["2c: model-level CSV has non-NA gamma1 for grassland/developed/both"]] <-
    all(!is.na(sp_model_rows$gamma1[sp_model_rows$model %in%
                                    c("grassland", "developed", "grassland_developed")]))
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
