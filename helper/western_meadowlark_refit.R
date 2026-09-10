## =============================================================================
## helper/western_meadowlark_refit.R
##
## Targeted extra-push refit for Western Meadowlark | anthro specifically --
## previously tests/test_western_meadowlark_debug.R (a standalone, manually
## run debug script). Moved into helper/ and revised into a callable
## function, refit_western_meadowlark(), so
## 1d_refit_nonconverged_species.R can invoke it directly as an automatic
## follow-up step instead of it being a separate script someone has to
## remember to run by hand.
##
## Why this species gets its own targeted step: Western Meadowlark | anthro
## is a known marginal case where the offending parameter is gamma1 ITSELF
## (not a sparse route's alpha/beta or a hyperparameter funnel -- see
## helper/diagnose_nonconvergence.R's output). It already went through
## 1d_refit_nonconverged_species.R's STANDARD refit (iter_warmup/
## iter_sampling = 4000, adapt_delta = 0.9) and still came back marginally
## over threshold (rhat 1.0124, ess_bulk 453.7). Since it's cheap to try
## harder on a single, isolated, marginal species, this pushes further:
## iter_warmup/iter_sampling = 8000, adapt_delta = 0.95 by default.
##
## Intended call site: 1d_refit_nonconverged_species.R, AFTER its standard
## refit loop and post-refit convergence re-check (Step 3), and ONLY if
## Western Meadowlark | anthro is still in the "still failing" set at that
## point -- so this never runs, or overwrites anything, on a run where the
## standard refit already fixed it.
##
## Requires functions/covariate_model_fitting.R (for species_to_f(),
## prepare_species_data(), fit_one_covariate_model()) and
## helper/model_convergence.R (for model_convergence_check()) to already be
## sourced in the calling script -- this file does not re-source them
## itself, since 1d_refit_nonconverged_species.R already has both in scope
## by the time this runs.
##
## Two ways to use this (same pattern as helper/gamma_lookup.R and
## helper/model_convergence.R):
##   1. Source with `western_meadowlark_refit_skip_autorun <- TRUE` set
##      first (as 1d does), then call refit_western_meadowlark(...)
##      directly with 1d's own rds_dir/route_info_dir/stan_data_dir/
##      covariate_lookups/model_single/model_base -- no settings to
##      duplicate, and it can never drift out of sync with what 1d was
##      actually just run with.
##   2. Standalone: Rscript helper/western_meadowlark_refit.R (or source()
##      from the console with nothing pre-set) sources its own dependencies
##      and runs the extra-push refit against this project's normal
##      production output/rds, data/route_info, data/stan_data directories
##      using the defaults above -- equivalent to the old debug script.
## =============================================================================

#' Extra-push refit for Western Meadowlark, targeting whichever of its
#' model tags are still failing whole-model convergence at the time this is
#' called. Skips any tag that's already converged (e.g. if this is being
#' re-run after an earlier attempt already fixed it).
#'
#' @param rds_dir,route_info_dir,stan_data_dir,cmdstanr_output_dir output
#'   locations -- pass the SAME production directories the calling script
#'   (1d_refit_nonconverged_species.R) is already using, so this overwrites
#'   the same files 2c reads, exactly like 1d's own standard refit does.
#' @param firstYear,lastYear,strat must match how the existing fit was
#'   produced
#' @param covariate_lookups named list of covariate lookup tables (as built
#'   by load_covariate()/1d's own setup) -- reused rather than rebuilt here
#' @param model_single,model_base compiled cmdstanr models -- reused rather
#'   than recompiled here (model_base only needed if "base" is in model_tags)
#' @param model_tags which tag(s) of THIS species to check/refit (default
#'   just "anthro" -- Western Meadowlark's base fit already converges fine)
#' @param iter_warmup,iter_sampling,adapt_delta,max_treedepth,show_exceptions
#'   passed straight through to fit_one_covariate_model() -- defaults here
#'   are pushed further than 1d's standard refit (4000/4000, adapt_delta =
#'   0.9), since this is specifically for the case where that wasn't enough
#' @return a named list of per-tag diagnostics (invisibly), or NULL entries
#'   for any tag that was skipped or failed
refit_western_meadowlark <- function(rds_dir, route_info_dir, stan_data_dir,
                                     cmdstanr_output_dir, firstYear, lastYear,
                                     strat, covariate_lookups,
                                     model_single, model_base = NULL,
                                     model_tags = c("anthro"),
                                     iter_warmup = 8000, iter_sampling = 8000,
                                     adapt_delta = 0.95, max_treedepth = 10,
                                     show_exceptions = TRUE) {

  species_name <- "Western Meadowlark"

  spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                     stringsAsFactors = FALSE)
  target_row <- spp_df %>%
    filter(Common.Name == species_name, in_bbs == TRUE) %>%
    distinct(Common.Name, Code, .keep_all = TRUE)

  if (nrow(target_row) != 1) {
    message("  [SKIP] Expected exactly one matching row for '", species_name,
            "' in spp_names_codes_group_aou.csv -- found ", nrow(target_row),
            ". Skipping Western Meadowlark extra-push refit.")
    return(invisible(NULL))
  }

  sp       <- target_row$Common.Name[1]
  sp_f     <- species_to_f(sp)
  sp_bbs   <- target_row$bbs_english[1]
  sp_group <- target_row$Group[1]

  # Nothing to do if every requested tag already converges -- check BEFORE
  # paying for prepare_species_data() (BBS pull + Voronoi neighbours).
  tags_to_refit <- Filter(function(tag) {
    current <- model_convergence_check(sp_f, tag, firstYear, lastYear, rds_dir)
    !(!is.null(current) && isTRUE(current$model_converged))
  }, model_tags)

  if (length(tags_to_refit) == 0) {
    cat("\n[Western Meadowlark extra-push refit] All requested tag(s) (",
        paste(model_tags, collapse = ", "), ") already converged -- nothing to do.\n", sep = "")
    return(invisible(NULL))
  }

  cat("\n================================================================\n")
  cat("  [Western Meadowlark extra-push refit] tags:", paste(tags_to_refit, collapse = ", "), "\n")
  cat("  iter_warmup =", iter_warmup, "| iter_sampling =", iter_sampling,
      "| adapt_delta =", adapt_delta, "\n")
  cat("================================================================\n")

  prepped <- tryCatch(
    prepare_species_data(species = sp, species_bbs = sp_bbs, strat = strat,
                         firstYear = firstYear, lastYear = lastYear,
                         covariate_lookups = covariate_lookups),
    error = function(e) {
      message("  [ERROR] Data prep failed for ", sp, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(prepped)) return(invisible(NULL))

  results <- list()
  for (tag in tags_to_refit) {
    current <- model_convergence_check(sp_f, tag, firstYear, lastYear, rds_dir)

    diagnostics <- tryCatch(
      fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                              prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                              model_base = model_base, model_single = model_single,
                              rds_dir = rds_dir, route_info_dir = route_info_dir,
                              cmdstanr_output_dir = cmdstanr_output_dir,
                              stan_data_dir = stan_data_dir,
                              chains = 4, iter_warmup = iter_warmup,
                              iter_sampling = iter_sampling,
                              adapt_delta = adapt_delta,
                              max_treedepth = max_treedepth,
                              show_exceptions = show_exceptions,
                              overwrite_route_info = TRUE),
      error = function(e) {
        message("  [ERROR] Extra-push refit '", tag, "' failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(diagnostics)) next

    diagnostics$group <- sp_group
    if (!is.null(current)) {
      diagnostics$before_max_rhat     <- current$model_max_rhat
      diagnostics$before_min_ess_bulk <- current$model_min_ess_bulk
    }
    results[[tag]] <- diagnostics

    # gamma1's own posterior row specifically -- the parameter actually in
    # question here, separate from the whole-model max_rhat/min_ess summary.
    out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    if (file.exists(summ_file)) {
      summ <- readRDS(summ_file)
      gamma_row <- summ[summ$variable == "gamma1", ]
      cat("\n    --- gamma1 posterior summary (extra-push refit) ---\n")
      old_na_print <- getOption("na.print")
      options(na.print = "NA")
      print(as.data.frame(gamma_row))
      options(na.print = old_na_print)
      if (nrow(gamma_row) == 1 && !is.na(gamma_row$rhat) && !is.na(gamma_row$ess_bulk)) {
        cat("    gamma1 converged (Rhat < 1.01 & ESS > 400)?",
            gamma_row$rhat < 1.01 && gamma_row$ess_bulk > 400, "\n")
      }
    }

    cat("  Done — extra-push refit", tag, "for", sp, "\n")
  }

  invisible(results)
}

## ==========================================================================
## Auto-run: only when sourced standalone (e.g. Rscript
## helper/western_meadowlark_refit.R), NOT when 1d_refit_nonconverged_species.R
## sources this file (1d sets western_meadowlark_refit_skip_autorun <- TRUE
## first, matching the pattern used for helper/gamma_lookup.R and
## helper/model_convergence.R). Builds its own dependencies/settings from
## this project's normal production defaults, equivalent to the old
## tests/test_western_meadowlark_debug.R script.
## ==========================================================================
if (!exists("western_meadowlark_refit_skip_autorun") || !isTRUE(western_meadowlark_refit_skip_autorun)) {

  library(bbsBayes2)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(sf)
  library(spdep)
  library(concaveman)
  library(here)

  here::i_am("helper/western_meadowlark_refit.R")

  source(here::here("functions", "covariate_model_fitting.R"))
  model_convergence_skip_autorun_wm <- exists("model_convergence_skip_autorun")
  model_convergence_skip_autorun <- TRUE
  source(here::here("helper", "model_convergence.R"))
  if (!model_convergence_skip_autorun_wm) rm(model_convergence_skip_autorun)

  if (!exists("firstYear")) firstYear <- 2010
  if (!exists("lastYear"))  lastYear  <- 2025
  if (!exists("strat"))     strat     <- "bcr"

  if (!exists("rds_dir"))              rds_dir              <- here::here("output", "rds")
  if (!exists("route_info_dir"))       route_info_dir       <- here::here("data", "route_info")
  if (!exists("stan_data_dir"))        stan_data_dir        <- here::here("data", "stan_data")
  if (!exists("cmdstanr_output_dir")) {
    cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output_refit")
    if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)
  }

  if (!exists("covariate_lookups")) {
    covariate_lookups <- list(
      anthro = load_covariate("Anthro.csv", "Anthro", "anthro")
    )
  }

  wm_model_tags <- if (exists("model_tags")) intersect(model_tags, c("base", "anthro")) else c("anthro")
  if (length(wm_model_tags) == 0) wm_model_tags <- c("anthro")

  model_base_wm   <- if ("base" %in% wm_model_tags) {
    cmdstan_model(here::here("models", "slope_iCAR_route_NB_New.stan"), stanc_options = list("O1"))
  } else NULL
  model_single_wm <- cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_covariate.stan"),
                                   stanc_options = list("O1"))

  refit_western_meadowlark(rds_dir = rds_dir, route_info_dir = route_info_dir,
                           stan_data_dir = stan_data_dir,
                           cmdstanr_output_dir = cmdstanr_output_dir,
                           firstYear = firstYear, lastYear = lastYear, strat = strat,
                           covariate_lookups = covariate_lookups,
                           model_single = model_single_wm, model_base = model_base_wm,
                           model_tags = wm_model_tags)
}
