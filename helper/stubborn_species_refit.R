## =============================================================================
## helper/stubborn_species_refit.R
##
## Official, generalized extra-push refit for the named set of "stubborn"
## species/tag combinations that have been confirmed (via
## tests/test_stubborn_species.R -- all 11 species/tag combos passed, per
## output/files/test_stubborn_species log 20260911.txt) to need more than
## 1d_refit_nonconverged_species.R's standard refit (iter_warmup/
## iter_sampling = 4000, adapt_delta = 0.9) to reach whole-model convergence:
##
##   Sharp-shinned Hawk  -- base, anthro
##   Hooded Merganser    -- base, anthro
##   Broad-winged Hawk   -- base, anthro
##   Belted Kingfisher   -- base, anthro
##   Western Meadowlark  -- anthro only (its base fit already converges fine)
##   Gray Catbird        -- anthro only
##   Red Crossbill       -- anthro only
##   Wood Duck           -- anthro only
##
## Uses the SAME extra-push settings validated in
## tests/test_stubborn_species.R: iter_warmup/iter_sampling = 8000,
## adapt_delta = 0.95 (vs. 1d's standard 4000/4000, adapt_delta = 0.9). Since
## these are isolated, already-identified species (not a ~600-species blanket
## run), it's cheap to try harder on just them.
##
## Now that all species/tags tests/test_stubborn_species.R covered are
## folded in here as defaults, that test script is no longer needed and can
## be removed.
##
## This supersedes helper/western_meadowlark_refit.R, which was originally
## written for Western Meadowlark | anthro alone and then generalized to a
## refit_stubborn_species() function under that (by-then-misleading)
## filename. That file has since been removed; this file is its replacement.
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
##   1. Source with `stubborn_species_refit_skip_autorun <- TRUE` set first
##      (as 1d does), then call refit_stubborn_species(...) directly with
##      1d's own rds_dir/route_info_dir/stan_data_dir/covariate_lookups/
##      model_single/model_base -- no settings to duplicate, and it can never
##      drift out of sync with what 1d was actually just run with.
##   2. Standalone: Rscript helper/stubborn_species_refit.R (or source() from
##      the console with nothing pre-set) sources its own dependencies and
##      runs the extra-push refit for the eight default stubborn species/tags
##      against this project's normal production output/rds, data/route_info,
##      data/stan_data directories.
## =============================================================================

#' Extra-push refit for a named set of species/tag combinations, targeting
#' whichever requested tags are still failing whole-model convergence at the
#' time this is called. Skips any species/tag that's already converged (e.g.
#' if this is being re-run after an earlier attempt already fixed it), and
#' skips a whole species cleanly if it isn't found in the species list.
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
#'   than recompiled here (model_base only needed if "base" is requested for
#'   any species)
#' @param species_tags named list, names are Common.Name values, each element
#'   a character vector of which tag(s) of THAT species to check/refit.
#'   Defaults to this project's eight confirmed stubborn species/tags (see
#'   header above).
#' @param iter_warmup,iter_sampling,adapt_delta,max_treedepth,show_exceptions
#'   passed straight through to fit_one_covariate_model() -- defaults here
#'   are pushed further than 1d's standard refit (4000/4000, adapt_delta =
#'   0.9), since this is specifically for cases where that wasn't enough
#' @return a named list of per-species-per-tag diagnostics (invisibly)
refit_stubborn_species <- function(rds_dir, route_info_dir, stan_data_dir,
                                   cmdstanr_output_dir, firstYear, lastYear,
                                   strat, covariate_lookups,
                                   model_single, model_base = NULL,
                                   species_tags = list(
                                     "Sharp-shinned Hawk" = c("base", "anthro"),
                                     "Hooded Merganser"   = c("base", "anthro"),
                                     "Broad-winged Hawk"  = c("base", "anthro"),
                                     "Belted Kingfisher"  = c("base", "anthro"),
                                     "Western Meadowlark" = c("anthro"),
                                     "Gray Catbird"       = c("anthro"),
                                     "Red Crossbill"      = c("anthro"),
                                     "Wood Duck"          = c("anthro")
                                   ),
                                   iter_warmup = 8000, iter_sampling = 8000,
                                   adapt_delta = 0.95, max_treedepth = 10,
                                   show_exceptions = TRUE) {

  spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                     stringsAsFactors = FALSE)

  results <- list()

  for (species_name in names(species_tags)) {
    requested_tags <- species_tags[[species_name]]

    target_row <- spp_df %>%
      filter(Common.Name == species_name, in_bbs == TRUE) %>%
      distinct(Common.Name, Code, .keep_all = TRUE)

    if (nrow(target_row) != 1) {
      message("  [SKIP] Expected exactly one matching row for '", species_name,
              "' in spp_names_codes_group_aou.csv -- found ", nrow(target_row),
              ". Skipping this species' extra-push refit.")
      next
    }

    sp       <- target_row$Common.Name[1]
    sp_f     <- species_to_f(sp)
    sp_bbs   <- target_row$bbs_english[1]
    sp_group <- target_row$Group[1]

    # Nothing to do for this species if every requested tag already
    # converges -- check BEFORE paying for prepare_species_data() (BBS pull +
    # Voronoi neighbours).
    tags_to_refit <- Filter(function(tag) {
      current <- model_convergence_check(sp_f, tag, firstYear, lastYear, rds_dir)
      !(!is.null(current) && isTRUE(current$model_converged))
    }, requested_tags)

    if (length(tags_to_refit) == 0) {
      cat("\n[Extra-push refit]", sp, "-- all requested tag(s) (",
          paste(requested_tags, collapse = ", "), ") already converged -- nothing to do.\n")
      next
    }

    cat("\n================================================================\n")
    cat("  [Extra-push refit]", sp, "-- tags:", paste(tags_to_refit, collapse = ", "), "\n")
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
    if (is.null(prepped)) next

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
      results[[paste(sp, tag, sep = " | ")]] <- diagnostics

      # gamma1's own posterior row specifically, when present (anthro only --
      # base has no gamma1) -- the parameter often in question for these
      # marginal fits, separate from the whole-model max_rhat/min_ess summary.
      out_base  <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
      summ_file <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
      if (file.exists(summ_file)) {
        summ <- readRDS(summ_file)
        gamma_row <- summ[summ$variable == "gamma1", ]
        if (nrow(gamma_row) == 1) {
          cat("\n    --- gamma1 posterior summary (extra-push refit) ---\n")
          old_na_print <- getOption("na.print")
          options(na.print = "NA")
          print(as.data.frame(gamma_row))
          options(na.print = old_na_print)
          if (!is.na(gamma_row$rhat) && !is.na(gamma_row$ess_bulk)) {
            cat("    gamma1 converged (Rhat < 1.01 & ESS > 400)?",
                gamma_row$rhat < 1.01 && gamma_row$ess_bulk > 400, "\n")
          }
        }
      }

      cat("  Done — extra-push refit", tag, "for", sp, "\n")
    }
  }

  invisible(results)
}

## ==========================================================================
## Auto-run: only when sourced standalone (e.g. Rscript
## helper/stubborn_species_refit.R), NOT when 1d_refit_nonconverged_species.R
## sources this file (1d sets stubborn_species_refit_skip_autorun <- TRUE
## first, matching the pattern used for helper/gamma_lookup.R and
## helper/model_convergence.R). Builds its own dependencies/settings from
## this project's normal production defaults, and runs the extra-push refit
## for the eight default stubborn species/tags.
## ==========================================================================
if (!exists("stubborn_species_refit_skip_autorun") || !isTRUE(stubborn_species_refit_skip_autorun)) {

  library(bbsBayes2)
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(sf)
  library(spdep)
  library(concaveman)
  library(here)

  here::i_am("helper/stubborn_species_refit.R")

  source(here::here("functions", "covariate_model_fitting.R"))
  model_convergence_skip_autorun_stub <- exists("model_convergence_skip_autorun")
  model_convergence_skip_autorun <- TRUE
  source(here::here("helper", "model_convergence.R"))
  if (!model_convergence_skip_autorun_stub) rm(model_convergence_skip_autorun)

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

  if (!exists("species_tags")) {
    species_tags <- list(
      "Sharp-shinned Hawk" = c("base", "anthro"),
      "Hooded Merganser"   = c("base", "anthro"),
      "Broad-winged Hawk"  = c("base", "anthro"),
      "Belted Kingfisher"  = c("base", "anthro"),
      "Western Meadowlark" = c("anthro"),
      "Gray Catbird"       = c("anthro"),
      "Red Crossbill"      = c("anthro"),
      "Wood Duck"          = c("anthro")
    )
  }

  model_base_stub   <- cmdstan_model(here::here("models", "slope_iCAR_route_NB_New.stan"), stanc_options = list("O1"))
  model_single_stub <- cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_covariate.stan"), stanc_options = list("O1"))

  refit_stubborn_species(rds_dir = rds_dir, route_info_dir = route_info_dir,
                         stan_data_dir = stan_data_dir,
                         cmdstanr_output_dir = cmdstanr_output_dir,
                         firstYear = firstYear, lastYear = lastYear, strat = strat,
                         covariate_lookups = covariate_lookups,
                         model_single = model_single_stub, model_base = model_base_stub,
                         species_tags = species_tags)
}
