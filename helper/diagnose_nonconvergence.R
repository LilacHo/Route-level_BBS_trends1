## =============================================================================
## helper/diagnose_nonconvergence.R
##
## For every species/tag combo that CURRENTLY fails the whole-model
## convergence criterion (Rhat < 1.01 and bulk ESS > 400, from
## helper/model_convergence.R), identify WHICH SPECIFIC PARAMETER(S) are
## driving the failure -- not just the aggregate max_rhat/min_ess_bulk
## number model_convergence_table() reports.
##
## Why this matters: model_convergence_table()'s max_rhat/min_ess_bulk is
## the single WORST value across potentially 1000+ parameters per fit (every
## route's alpha[r] and beta[r], plus gamma1 and a handful of
## hyperparameters). A fit can get flagged as "non-converged" because of ONE
## sparsely-sampled route's alpha -- which may never fully converge no
## matter how long you run it, and doesn't call the species-level trend
## (beta) or covariate effect (gamma1) into question -- or it can be flagged
## because gamma1 or a route's beta itself (the actual quantities of
## scientific interest) are poorly mixed, which is a more serious problem.
## This script tells you which situation you're actually in, for each
## currently-failing species/tag, before spending more compute on another
## blanket refit.
##
## For each failing species/tag, reads its FULL summ_fit.rds (every
## parameter, not just gamma1 or the whole-model max/min) and reports:
##   - every parameter that individually fails Rhat < 1.01 or bulk ESS > 400
##     (not just the single worst one)
##   - each offending parameter's TYPE, via classify_parameter(): route-level
##     alpha, route-level beta, gamma1, or "other" (hyperparameters like
##     BETA/sdbeta/sdobs/observer effects/etc.) -- so you can see at a
##     glance whether beta/gamma1 (the parameters this project reports on)
##     are actually implicated, or whether it's confined to a handful of
##     individual route intercepts/slopes
##   - how MANY parameters are offending (one lagging route vs. many),
##     which distinguishes an isolated, likely-benign case from a fit that's
##     more broadly unconverged
##
## Writes output/files/nonconvergence_diagnosis_<firstYear>_<lastYear>.csv
## (one row per offending parameter, so you can filter/sort freely) and
## prints a per-species summary plus an overall cross-tab of which parameter
## TYPE is most often the worst offender across all currently-failing fits.
##
## Read-only: does not fit or refit anything, does not modify any existing
## output. Safe to re-run any time; it just reflects whatever is currently
## on disk.
##
## Two ways to use this (same pattern as helper/gamma_lookup.R and
## helper/model_convergence.R):
##   1. Standalone: Rscript helper/diagnose_nonconvergence.R (or source()
##      with nothing pre-set) falls back to this project's current
##      full-species, base+anthro defaults.
##   2. Insert into another script that already has target_spp/model_tags/
##      firstYear/lastYear/rds_dir in scope (e.g. after
##      1d_refit_nonconverged_species.R's final convergence re-check) --
##      set those variables and source this file.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/diagnose_nonconvergence.R")

# Pull in model_convergence_check()/model_convergence_table() to identify
# the current failing list, without triggering its own auto-run block.
model_convergence_skip_autorun <- TRUE
source(here::here("helper", "model_convergence.R"))

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Classify a Stan parameter name (as it appears in
#' posterior::summarise_draws()'s "variable" column) into a broad category,
#' so offending parameters can be grouped by what they actually represent.
#' Regex-based and deliberately falls back to "other" for anything
#' unrecognized, rather than guessing -- so this stays correct even if the
#' Stan model's hyperparameter names change later.
#'
#' @param varname character vector of variable names, e.g. "alpha[142]",
#'   "beta[57]", "gamma1", "sdbeta", "BETA", "sdobs[3]"
#' @return character vector of categories: "route_alpha", "route_beta",
#'   "hyperparameter_sd", "gamma1", or "other"
classify_parameter <- function(varname) {
  # Matches whatever the actual per-route parameter is called, including
  # this model's non-centered parameterization (alpha_raw[r], beta_raw[r],
  # beta_raw_space[r], etc.) -- anything starting with "alpha"/"beta" AND
  # indexed by route (has a "[") is route-level, regardless of the exact
  # "_raw"/"_space" suffix naming. sd* hyperparameters (sdalpha,
  # sdbeta_space, ...) get their own bucket since a slow-mixing variance
  # hyperparameter is a different (and common, often benign) issue from a
  # slow-mixing route-level or gamma1 parameter -- classic "funnel"
  # geometry in hierarchical models when a species has little real
  # spatial/temporal structure for that variance component to explain.
  dplyr::case_when(
    grepl("^alpha", varname) & grepl("\\[", varname) ~ "route_alpha",
    grepl("^beta",  varname) & grepl("\\[", varname) ~ "route_beta",
    varname == "gamma1"                              ~ "gamma1",
    grepl("^sd",    varname)                          ~ "hyperparameter_sd",
    TRUE                                              ~ "other"
  )
}

#' Read one species/model's *_summ_fit.rds and report every parameter that
#' individually fails Rhat < rhat_threshold or bulk ESS > ess_threshold,
#' classified by type. Returns NULL (with a message) if the file is missing.
#'
#' @param species_f file-safe species name
#' @param tag model tag, e.g. "base" or "anthro"
#' @param firstYear,lastYear must match the summ_fit.rds filename
#' @param rds_dir directory containing the *_summ_fit.rds files
#' @param rhat_threshold,ess_threshold same criterion as
#'   helper/model_convergence.R (defaults: 1.01 / 400)
#' @return a data.frame, one row per offending parameter (species_f, model,
#'   variable, category, rhat, ess_bulk, fails_rhat, fails_ess), or NULL
diagnose_convergence_failure <- function(species_f, tag, firstYear = 2010, lastYear = 2025,
                                         rds_dir = here::here("output", "rds"),
                                         rhat_threshold = 1.01, ess_threshold = 400) {
  summ_file <- file.path(rds_dir,
                         paste0(species_f, "_iCAR_", tag, "_", firstYear, "_", lastYear, "_summ_fit.rds"))

  if (!file.exists(summ_file)) {
    message("  [MISSING] ", basename(summ_file))
    return(NULL)
  }

  summ <- readRDS(summ_file)

  offenders <- summ %>%
    filter((!is.na(rhat) & rhat >= rhat_threshold) |
           (!is.na(ess_bulk) & ess_bulk <= ess_threshold)) %>%
    transmute(species_f  = species_f,
             model      = tag,
             variable   = variable,
             category   = classify_parameter(variable),
             rhat       = rhat,
             ess_bulk   = ess_bulk,
             fails_rhat = !is.na(rhat) & rhat >= rhat_threshold,
             fails_ess  = !is.na(ess_bulk) & ess_bulk <= ess_threshold)

  if (nrow(offenders) == 0) {
    # Shouldn't normally happen if this species/tag was pulled from a
    # failing list, but guard against a threshold mismatch between here and
    # model_convergence.R rather than erroring.
    message("  [NO OFFENDERS FOUND] ", basename(summ_file),
            " -- thresholds here may not match how it was flagged as failing.")
    return(NULL)
  }

  offenders %>% arrange(desc(rhat))
}

#' Build the combined non-convergence diagnosis across every CURRENTLY
#' failing species/tag (as determined fresh from model_convergence_table()),
#' write it to
#' output/files/nonconvergence_diagnosis_<firstYear>_<lastYear>.csv, and
#' print a per-species summary plus an overall cross-tab of worst-offender
#' category.
#'
#' @param target_spp,model_tags,firstYear,lastYear,rds_dir same as
#'   model_convergence_table()
#' @param out_dir where to write the combined CSV (default output/files)
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the combined offenders table (invisibly), or NULL if nothing to
#'   diagnose (every fit already converged)
diagnose_nonconvergence_table <- function(target_spp, model_tags, firstYear, lastYear,
                                          rds_dir, out_dir = here::here("output", "files"),
                                          write_csv = TRUE) {

  cat("=== Diagnosing currently-failing species/tags: which parameter(s) are driving it? ===\n")

  convergence_now <- model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                                             firstYear = firstYear, lastYear = lastYear,
                                             rds_dir = rds_dir, bird_group = "diagnosis_scratch",
                                             write_csv = FALSE)

  if (is.null(convergence_now)) {
    message("No existing fits found -- nothing to diagnose.")
    return(invisible(NULL))
  }

  failing <- convergence_now %>% filter(!model_converged)

  if (nrow(failing) == 0) {
    cat("Nothing currently fails the convergence criterion -- nothing to diagnose.\n")
    return(invisible(NULL))
  }

  cat("\n", nrow(failing), "species/tag combo(s) currently fail convergence; reading full",
      "posterior summaries to find the offending parameter(s)...\n", sep = "")

  rows <- list()
  for (i in seq_len(nrow(failing))) {
    sp      <- failing$species[i]
    sp_f    <- species_to_f(sp)
    sp_code <- failing$species_code[i]
    tag     <- failing$model[i]

    r <- diagnose_convergence_failure(sp_f, tag, firstYear, lastYear, rds_dir)
    if (!is.null(r)) {
      r$species      <- sp
      r$species_code <- sp_code
      rows[[paste(sp, tag, sep = " | ")]] <- r
    }
  }

  if (length(rows) == 0) {
    message("No offending parameters found for any currently-failing species/tag -- ",
            "check that rhat_threshold/ess_threshold match model_convergence.R's criterion.")
    return(invisible(NULL))
  }

  offenders_all <- bind_rows(rows) %>%
    select(species, species_code, model, variable, category, rhat, ess_bulk,
           fails_rhat, fails_ess, everything(), -species_f)

  # Per-species/tag summary: worst rhat param, worst (lowest) ess param, and
  # how many parameters total are implicated -- one lagging route reads very
  # differently from dozens of them ------------------------------------------
  per_species_summary <- offenders_all %>%
    group_by(species, species_code, model) %>%
    summarise(
      n_offending_params = n(),
      worst_rhat_variable = variable[which.max(rhat)],
      worst_rhat_category = category[which.max(rhat)],
      worst_rhat_value    = max(rhat, na.rm = TRUE),
      worst_ess_variable  = variable[which.min(ess_bulk)],
      worst_ess_category  = category[which.min(ess_bulk)],
      worst_ess_value     = min(ess_bulk, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(worst_rhat_value))

  cat("\n=== Per-species/tag: worst offending parameter ===\n")
  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(per_species_summary))
  options(na.print = old_na_print)

  cat("\n=== Worst-offender category, across all currently-failing species/tags ===\n")
  cat("(route_alpha/route_beta = one specific route's intercept/slope -- often just data-sparsity\n")
  cat(" in that one route, not a problem with the species-level estimate. gamma1 = the covariate\n")
  cat(" effect itself. 'other' = a hyperparameter, e.g. BETA/sdbeta/sdobs/observer effects.)\n")
  print(per_species_summary %>% count(worst_rhat_category, name = "n_species_tags_worst_rhat_here"))

  cat("\n=== How many parameters are implicated per species/tag (1 = isolated, many = broader issue) ===\n")
  print(as.data.frame(per_species_summary %>%
                        select(species, species_code, model, n_offending_params) %>%
                        arrange(desc(n_offending_params))))

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("nonconvergence_diagnosis_", firstYear, "_", lastYear, ".csv"))
    write.csv(offenders_all, out_csv, row.names = FALSE)
    cat("\nFull per-parameter diagnosis written to:", out_csv, "\n")

    summary_csv <- file.path(out_dir, paste0("nonconvergence_diagnosis_summary_", firstYear, "_", lastYear, ".csv"))
    write.csv(per_species_summary, summary_csv, row.names = FALSE)
    cat("Per-species/tag summary written to:", summary_csv, "\n")
  }

  invisible(offenders_all)
}

## ==========================================================================
## Auto-run: same `if (!exists(...))` soft-coded default pattern as
## helper/gamma_lookup.R and helper/model_convergence.R. Set
## `diagnose_nonconvergence_skip_autorun <- TRUE` before sourcing to load
## just the functions above without running anything.
## ==========================================================================
if (!exists("diagnose_nonconvergence_skip_autorun") || !isTRUE(diagnose_nonconvergence_skip_autorun)) {

  if (!exists("firstYear"))  firstYear  <- 2010
  if (!exists("lastYear"))   lastYear   <- 2025
  if (!exists("model_tags")) model_tags <- c("base", "anthro")
  if (!exists("rds_dir"))    rds_dir    <- here::here("output", "rds")

  if (!exists("target_spp")) {
    spp_df_dc <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                          stringsAsFactors = FALSE)
    target_spp <- spp_df_dc %>%
      filter(in_bbs == TRUE) %>%
      distinct(Common.Name, Code, .keep_all = TRUE) %>%
      arrange(Common.Name)
  }

  diagnose_nonconvergence_table(target_spp = target_spp, model_tags = model_tags,
                                firstYear = firstYear, lastYear = lastYear,
                                rds_dir = rds_dir)
}
