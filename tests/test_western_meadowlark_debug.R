## =============================================================================
## tests/test_western_meadowlark_debug.R
##
## Standalone diagnostic re-run of ONE species/model -- Western Meadowlark,
## "anthro" -- which is the one case out of the 12 still-failing
## species/tags where the offending parameter is gamma1 ITSELF (not a
## sparse route's alpha/beta or a hyperparameter funnel; see
## helper/diagnose_nonconvergence.R's output):
##   Western Meadowlark | anthro | gamma1 | rhat 1.0124 | ess_bulk 453.7
##   (3 offending parameters total -- gamma1 plus a couple of closely
##   related/correlated quantities)
##
## This already went through 1d_refit_nonconverged_species.R's standard
## refit (iter_warmup/iter_sampling = 4000, adapt_delta = 0.9) and STILL
## came back marginally over threshold -- so this script pushes further
## (8000/8000 iterations, adapt_delta = 0.95) specifically for this one
## isolated, marginal case, since it's cheap to try harder on a single
## species and gamma1 is the actual reported covariate effect (unlike the
## route-alpha/hyperparameter-only cases, which don't need this).
##
## "base" isn't included by default -- Western Meadowlark's base fit
## already converges fine; only anthro is flagged. Add "base" to
## model_tags_to_test below if you want to double-check it anyway.
##
## Data prep and fitting both come from functions/covariate_model_fitting.R
## (shared with 1c/1d/the Green Heron debug test -- no duplicated copies
## here). Prints gamma1's own posterior summary row specifically (not just
## the whole-model max_rhat/min_ess) so it's easy to see directly whether
## the extra push actually resolved gamma1's convergence, separate from
## whatever else in the fit might still be borderline.
##
## This is READ-ONLY with respect to the production pipeline: it does not
## touch 1c/1d, model_tags, force_refit, or any file the real run reads/
## writes. Everything this script saves goes to the same *_debug subfolders
## the Green Heron debug script uses (output/rds_debug/, etc. -- safe to
## share since output filenames are species-specific) so it can never
## collide with or overwrite real run output.
##
## Usage: source this file, or Rscript tests/test_western_meadowlark_debug.R.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("tests/test_western_meadowlark_debug.R")

# species_to_f()/load_covariate()/prepare_species_data()/
# fit_one_covariate_model() -- and the neighbours_define_voronoi()/
# posterior_summary_functions.R sources they depend on -- all come from
# here, shared with 1c/1d/the Green Heron debug test.
source(here::here("functions", "covariate_model_fitting.R"))

# Settings -- match 1c_species_iCAR_covariates.R, except what this script
# exists to change ----------------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025
strat      <- "bcr"
species_name       <- "Western Meadowlark"
model_tags_to_test <- c("anthro")   # base already converges; add "base" here too if desired

# Pushed further than 1d_refit_nonconverged_species.R's standard refit
# (which used 4000/4000, adapt_delta = 0.9 and still left this marginal) --
# reasonable to try harder given it's a single, isolated, borderline case.
debug_iter_warmup    <- 8000
debug_iter_sampling  <- 8000
debug_adapt_delta    <- 0.95
debug_max_treedepth  <- 10
debug_show_exceptions <- TRUE

# Debug-only output locations -- SAME shared debug dirs the Green Heron
# debug script uses (filenames are species-specific, so this is safe to
# share); separate from production so nothing here can ever overwrite/
# collide with a real 1c/1d run ----------------------------------------------
cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output_debug")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)
rds_dir <- here::here("output", "rds_debug")
if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
route_info_dir <- here::here("data", "route_info_debug")
if (!dir.exists(route_info_dir)) dir.create(route_info_dir, recursive = TRUE)
stan_data_dir <- here::here("data", "stan_data_debug")
if (!dir.exists(stan_data_dir)) dir.create(stan_data_dir, recursive = TRUE)

# Species lookup -- same source table as 1c ---------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_row <- spp_df %>%
  filter(Common.Name == species_name, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE)

if (nrow(target_row) != 1) {
  stop("Expected exactly one matching row for '", species_name,
       "' in spp_names_codes_group_aou.csv -- found ", nrow(target_row),
       ". Check spelling / in_bbs filter.")
}

sp       <- target_row$Common.Name[1]
sp_f     <- species_to_f(sp)
sp_code  <- target_row$Code[1]
sp_bbs   <- target_row$bbs_english[1]
sp_group <- target_row$Group[1]

cat("=== Debug run:", sp, "(", sp_group, ") | tags:",
    paste(model_tags_to_test, collapse = ", "), "===\n")
cat("    iter_warmup =", debug_iter_warmup, "| iter_sampling =", debug_iter_sampling,
    "| adapt_delta =", debug_adapt_delta, "\n")

# Covariate spec -- same as 1c_species_iCAR_covariates.R --------------------
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# "base" still needs the anthro lookup available (fair-comparison
# principle), even if only testing "anthro" here -- load it unconditionally.
covariate_lookups <- list(
  anthro = load_covariate(covariate_specs$anthro$file,
                          covariate_specs$anthro$value_col, "anthro")
)
cat("  anthro covariate:", nrow(covariate_lookups$anthro), "route-year rows (",
    length(unique(covariate_lookups$anthro$route_key)), "routes)\n")

# Compile only the Stan model(s) actually needed for the tag(s) under test --
model_base   <- if ("base" %in% model_tags_to_test) {
  cmdstan_model(here::here("models", "slope_iCAR_route_NB_New.stan"), stanc_options = list("O1"))
} else NULL

model_single <- if (any(model_tags_to_test != "base")) {
  cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_covariate.stan"), stanc_options = list("O1"))
} else NULL

# --- Run it ----------------------------------------------------------------
prepped <- prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                                firstYear = firstYear, lastYear = lastYear,
                                covariate_lookups = covariate_lookups)

debug_results <- list()
for (tag in model_tags_to_test) {
  cat("\n================================================================\n")
  cat("  Testing tag:", tag, "(iter_warmup =", debug_iter_warmup,
      ", iter_sampling =", debug_iter_sampling, ", adapt_delta =", debug_adapt_delta, ")\n")
  cat("================================================================\n")
  debug_results[[tag]] <- tryCatch(
    fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                            prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                            model_base = model_base, model_single = model_single,
                            rds_dir = rds_dir, route_info_dir = route_info_dir,
                            cmdstanr_output_dir = cmdstanr_output_dir,
                            stan_data_dir = stan_data_dir,
                            chains = 4, iter_warmup = debug_iter_warmup,
                            iter_sampling = debug_iter_sampling,
                            adapt_delta = debug_adapt_delta,
                            max_treedepth = debug_max_treedepth,
                            show_exceptions = debug_show_exceptions,
                            overwrite_route_info = TRUE),
    error = function(e) {
      cat("\n    [FAILED] tag:", tag, "\n")
      cat("    Error message:", conditionMessage(e), "\n")
      NULL
    }
  )

  # gamma1's own posterior row specifically -- the parameter actually in
  # question here, separate from the whole-model max_rhat/min_ess summary.
  if (!is.null(debug_results[[tag]])) {
    out_base <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    if (file.exists(summ_file)) {
      summ <- readRDS(summ_file)
      gamma_row <- summ[summ$variable == "gamma1", ]
      cat("\n    --- gamma1 posterior summary (this run) ---\n")
      old_na_print <- getOption("na.print")
      options(na.print = "NA")
      print(as.data.frame(gamma_row))
      options(na.print = old_na_print)
      cat("    gamma1 converged (Rhat < 1.01 & ESS > 400)?",
          !is.na(gamma_row$rhat) && !is.na(gamma_row$ess_bulk) &&
            gamma_row$rhat < 1.01 && gamma_row$ess_bulk > 400, "\n")
    }
  }
}

cat("\n\n=== DEBUG RUN DONE ===\n")
cat("Results saved (if successful) to:", rds_dir, "\n")
cat("Compare gamma1's rhat/ess above against the pre-push values (rhat 1.0124, ess_bulk 453.7\n")
cat("after 1d's standard 4000/4000 refit) to see whether pushing to", debug_iter_warmup, "/",
    debug_iter_sampling, "iterations and adapt_delta =", debug_adapt_delta, "actually helped.\n")
