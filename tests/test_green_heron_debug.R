## =============================================================================
## tests/test_green_heron_debug.R
##
## Standalone diagnostic re-run of ONE species/model -- Green Heron, "base" --
## to chase down the error seen during the full 1c_species_iCAR_covariates.R
## run:
##   [ERROR] Model 'base' failed for Green Heron: arguments imply differing
##   number of rows: 4, 3
##
## Working theory (see chat): chains = 4 is requested, but Green Heron's fit
## prints "Route: 1465 | Obs: 16550 | Edges: 4163" before failing -- a large,
## well-powered CAR model, so this doesn't look like a too-little-data
## problem. More likely one of the four chains crashed/terminated during
## warmup or sampling while the other three finished fine.
##
## Data prep and fitting both come from functions/covariate_model_fitting.R
## (shared with 1c_species_iCAR_covariates.R and
## 1d_refit_nonconverged_species.R -- no duplicated copies here). That
## shared fit_one_covariate_model() now includes a chain-completion gate: it
## checks stanfit$return_codes() right after sampling, BEFORE $summary()/
## $time() are touched, and prints which chain(s) failed if any did -- so
## this script's only real job is to call it with show_exceptions = TRUE
## (which 1c's production run leaves FALSE, suppressing each chain's own
## crash/exception text) so the real per-chain error surfaces in console
## output instead of being swallowed.
##
## This is READ-ONLY with respect to the production pipeline: it does not
## touch 1c_species_iCAR_covariates.R, model_tags, force_refit, or any file
## the real run reads/writes. Everything this script saves goes to separate
## *_debug subfolders (output/rds_debug/, a debug cmdstan temp dir) so it can
## never collide with or overwrite real run output.
##
## Usage: source this file, or Rscript tests/test_green_heron_debug.R.
## Edit species_name / model_tags_to_test below to point at a different
## species or add "anthro" if you want to check that tag too.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("tests/test_green_heron_debug.R")

# species_to_f()/load_covariate()/prepare_species_data()/
# fit_one_covariate_model() -- and the neighbours_define_voronoi()/
# posterior_summary_functions.R sources they depend on -- all come from
# here, shared with 1c_species_iCAR_covariates.R and
# 1d_refit_nonconverged_species.R.
source(here::here("functions", "covariate_model_fitting.R"))

# Settings -- match 1c_species_iCAR_covariates.R, except what this script
# exists to change ----------------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025
strat      <- "bcr"
species_name       <- "Green Heron"
model_tags_to_test <- c("base")   # add "anthro" here too if you want to test that tag

# Debug-only output locations -- separate from production dirs so nothing
# here can ever overwrite/collide with a real 1c run ------------------------
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

# Covariate spec -- same as 1c_species_iCAR_covariates.R --------------------
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# 1c fits BOTH tags (including "base") on the same anthro-non-NA reduced
# dataset, so the anthro lookup is needed regardless of which tags are being
# tested here -- load it unconditionally.
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
  cat("  Testing tag:", tag, "(DEBUG: show_exceptions = TRUE)\n")
  cat("================================================================\n")
  debug_results[[tag]] <- tryCatch(
    fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                            prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                            model_base = model_base, model_single = model_single,
                            rds_dir = rds_dir, route_info_dir = route_info_dir,
                            cmdstanr_output_dir = cmdstanr_output_dir,
                            stan_data_dir = stan_data_dir,
                            chains = 4, iter_warmup = 2000, iter_sampling = 2000,
                            adapt_delta = 0.8, max_treedepth = 10,
                            show_exceptions = TRUE,
                            overwrite_route_info = TRUE),
    error = function(e) {
      cat("\n    [STILL FAILED] tag:", tag, "\n")
      cat("    Error message:", conditionMessage(e), "\n")
      NULL
    }
  )
}

cat("\n\n=== DEBUG RUN DONE ===\n")
cat("Results saved (if successful) to:", rds_dir, "\n")
cat("If a chain failed, its return code and any surfaced exception text should\n")
cat("be visible above (printed by fit_one_covariate_model() itself) -- that's\n")
cat("the real root cause to address, if it turns out to recur beyond this one\n")
cat("species (e.g. via 1d_refit_nonconverged_species.R's higher adapt_delta).\n")
