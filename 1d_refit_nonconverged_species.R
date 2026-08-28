## =============================================================================
## 1d_refit_nonconverged_species.R
##
## Targeted re-fit of ONLY the species/tag combos that fail the whole-model
## convergence criterion introduced in helper/model_convergence.R:
##   "species are accepted in models converged based on Rhat < 1.01 and bulk
##   effective sample sizes > 400" -- checked across EVERY parameter in the
##   fit (alpha, beta, gamma1 if present), for EVERY model tag, including
##   "base".
##
## Rather than re-running all ~600 species x 2 models
## (1c_species_iCAR_covariates.R) with longer iterations across the board,
## this script:
##   1. Reads every existing *_summ_fit.rds DIRECTLY off disk (via
##      helper/model_convergence.R's model_convergence_table() -- NOT the
##      diagnostics CSV; see that file's header comment for why that CSV is
##      an incomplete source across a resumable, multi-session run) to find
##      exactly which species/tag combos currently fail the criterion.
##   2. Re-fits ONLY those, with increased iter_warmup/iter_sampling and a
##      higher adapt_delta (more numerical stability during warmup), on the
##      SAME reduced (anthro-non-NA) dataset 1c itself uses.
##   3. Before refitting, COPIES the existing (under-converged) *_summ_fit.rds/
##      *_stanfit.rds/*_route_info.rds/*_stan_data.RData into parallel
##      output/rds_under-converged/, data/route_info_under-converged/, and
##      data/stan_data_under-converged/ backup folders (see
##      archive_existing_fit() below), THEN refits and overwrites the
##      originals at the SAME filenames 1c_species_iCAR_covariates.R uses --
##      so the old fit is preserved for comparison, but 2c still picks up
##      the new one from the normal production location without any change.
##
## Data prep and fitting both come from functions/covariate_model_fitting.R
## (shared with 1c_species_iCAR_covariates.R and
## tests/test_green_heron_debug.R -- no duplicated copies here). This
## script's only fitting-specific behavior is passing DIFFERENT settings
## into the shared fit_one_covariate_model(): more iterations, higher
## adapt_delta, show_exceptions = TRUE, and overwrite_route_info = TRUE.
##
## Resumable: right before refitting each species/tag, re-checks its CURRENT
## convergence status. If an earlier, interrupted run of THIS script already
## fixed it, it's skipped -- same "safe to re-run, picks up where it left
## off" principle used throughout this project (1c/2c/3c).
##
## NOTE on Green Heron specifically: the chain-crash issue discussed
## separately (a chain terminating with a nonzero return code, producing the
## "differing number of rows: 4, 3" error) is a DIFFERENT failure mode from
## slow mixing/convergence, and increasing iterations alone is not expected
## to fix it. The shared fit_one_covariate_model() now checks
## stanfit$return_codes() right after sampling and prints which chain failed
## if any did, so that should be visible in this script's console output too
## if it's the same issue. See tests/test_green_heron_debug.R for a
## single-species deep-dive on that specific failure mode.
##
## IMPORTANT: after this script finishes, re-run
## 2c_generate_route_trend_csvs_covariates.R so the improved fits propagate
## into all_route_trends_*.csv / all_model_level_summary_*.csv -- this
## script only touches output/rds/, data/route_info/, and data/stan_data/,
## not 2c's combined output.
## =============================================================================

library(bbsBayes2)
library(tidyverse)
library(cmdstanr)
library(posterior)
library(sf)
library(spdep)
library(concaveman)
library(here)

here::i_am("1d_refit_nonconverged_species.R")

# species_to_f()/load_covariate()/prepare_species_data()/
# fit_one_covariate_model() -- and the neighbours_define_voronoi()/
# posterior_summary_functions.R sources they depend on -- all come from
# here, shared with 1c_species_iCAR_covariates.R and
# tests/test_green_heron_debug.R.
source(here::here("functions", "covariate_model_fitting.R"))

# Pull in model_convergence_check()/model_convergence_table() to identify
# targets, without triggering its own auto-run block (we call it ourselves
# below with this script's own settings) ------------------------------------
model_convergence_skip_autorun <- TRUE
source(here::here("helper", "model_convergence.R"))

cmdstanr_output_dir <- file.path(tempdir(), "cmdstan_output_refit")
if (!dir.exists(cmdstanr_output_dir)) dir.create(cmdstanr_output_dir, recursive = TRUE)

# Settings -- match 1c_species_iCAR_covariates.R except the refit-specific
# knobs called out below -----------------------------------------------------
firstYear  <- 2010
lastYear   <- 2025
strat      <- "bcr"
model_tags <- c("base", "anthro")   # same tags 1c fits; base included since
                                     # model_convergence_table() checks it too

# Refit-specific settings -- bumped relative to 1c's production defaults
# (chains=4, iter_sampling=2000, iter_warmup=2000, adapt_delta=0.8,
# max_treedepth=10, show_exceptions=FALSE) to give these specific
# non-converged fits a better shot. Adjust here if a first pass isn't enough.
refit_iter_warmup    <- 4000   # was 2000
refit_iter_sampling  <- 4000   # was 2000
refit_adapt_delta    <- 0.9    # was 0.8 -- more conservative step-size
                                # adaptation, targets numerical stability
refit_max_treedepth  <- 10     # unchanged; raise if treedepth warnings show up
refit_show_exceptions <- TRUE  # unlike 1c's production FALSE -- these are
                                # known-problematic fits, so surface any
                                # per-chain crash/exception text instead of
                                # suppressing it

# Output directories -- SAME as 1c_species_iCAR_covariates.R, since this
# script overwrites those fits in place rather than creating a parallel set
# of files ------------------------------------------------------------------
output_dir     <- here::here("output")
rds_dir        <- here::here("output", "rds")
route_info_dir <- here::here("data", "route_info")
stan_data_dir  <- here::here("data", "stan_data")
if (!dir.exists(stan_data_dir)) dir.create(stan_data_dir, recursive = TRUE)
if (!dir.exists(here::here("data", "maps"))) dir.create(here::here("data", "maps"), recursive = TRUE)

# Backup directories -- mirror the production structure with an
# "_under-converged" suffix. Before each refit, the CURRENT (about-to-be-
# overwritten) fit files are copied here first, so the under-converged fit
# isn't lost if you ever want to compare before/after, or if the refit
# itself doesn't actually improve things -----------------------------------
rds_backup_dir        <- here::here("output", "rds_under-converged")
route_info_backup_dir <- here::here("data", "route_info_under-converged")
stan_data_backup_dir  <- here::here("data", "stan_data_under-converged")
if (!dir.exists(rds_backup_dir))        dir.create(rds_backup_dir, recursive = TRUE)
if (!dir.exists(route_info_backup_dir)) dir.create(route_info_backup_dir, recursive = TRUE)
if (!dir.exists(stan_data_backup_dir))  dir.create(stan_data_backup_dir, recursive = TRUE)

#' Archive the CURRENT (about-to-be-overwritten) fit files for one
#' species/tag into the *_under-converged backup directories, before
#' fit_one_covariate_model() overwrites them with a refit.
#'
#' Deliberately COPIES rather than moves (the original stays in place in
#' rds_dir/route_info_dir/stan_data_dir until fit_one_covariate_model()
#' itself overwrites it upon a SUCCESSFUL refit): if the refit throws an
#' error partway through (e.g. the same kind of chain crash discussed for
#' Green Heron), the original under-converged fit is never touched, so
#' nothing is ever lost to a failed refit attempt. The backup copy is
#' redundant in that case, which is a fine trade for never risking deleting
#' the only existing fit for a species/tag.
#'
#' Overwrites any PREVIOUS backup for the same species/tag -- this keeps
#' only the most recent pre-refit snapshot per species/tag, not a full
#' history across multiple refit attempts.
archive_existing_fit <- function(species_f, model_tag, firstYear, lastYear) {
  out_base <- paste0(species_f, "_iCAR_", model_tag, "_", firstYear, "_", lastYear)

  copy_if_exists <- function(from, to_dir) {
    if (file.exists(from)) {
      to <- file.path(to_dir, basename(from))
      file.copy(from, to, overwrite = TRUE)
      cat("    Backed up pre-refit file:", basename(from), "->", to_dir, "\n")
    }
  }

  copy_if_exists(file.path(rds_dir, paste0(out_base, "_stanfit.rds")), rds_backup_dir)
  copy_if_exists(file.path(rds_dir, paste0(out_base, "_summ_fit.rds")), rds_backup_dir)
  copy_if_exists(file.path(route_info_dir,
                           paste0(species_f, "_", model_tag, "_", firstYear, "_", lastYear, "_route_info.rds")),
                 route_info_backup_dir)
  copy_if_exists(file.path(stan_data_dir,
                           paste0(species_f, "_", model_tag, "_", firstYear, "_", lastYear, "_stan_data.RData")),
                 stan_data_backup_dir)
}

# Covariate spec -- same as 1c_species_iCAR_covariates.R --------------------
covariate_specs <- list(
  anthro = list(file = "Anthro.csv", value_col = "Anthro", rescale = FALSE)
)

# "base" still needs the anthro lookup available, since 1c fits base on the
# SAME anthro-non-NA reduced dataset as anthro (fair-comparison principle
# used throughout this project) -- load it unconditionally, exactly like 1c
# does.
cat("Loading covariates...\n")
covariate_lookups <- list(
  anthro = load_covariate(covariate_specs$anthro$file, covariate_specs$anthro$value_col, "anthro")
)
cat("  anthro covariate:", nrow(covariate_lookups$anthro), "route-year rows (",
    length(unique(covariate_lookups$anthro$route_key)), "routes)\n")

# Species list -- every in_bbs species, same as 1c --------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)
target_spp <- spp_df %>%
  filter(in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

# ==========================================================================
# Step 1: find who currently fails whole-model convergence ------------------
# ==========================================================================
cat("\n=== Checking current convergence status across all species/tags ===\n")
convergence_now <- model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                                           firstYear = firstYear, lastYear = lastYear,
                                           rds_dir = rds_dir, bird_group = "all_species_anthro",
                                           write_csv = FALSE)

if (is.null(convergence_now)) {
  stop("No existing fits found at all -- run 1c_species_iCAR_covariates.R first.")
}

failing <- convergence_now %>% filter(!model_converged)

cat("\n", nrow(failing), "species/tag combo(s) currently fail the convergence criterion ",
    "(Rhat < 1.01 and bulk ESS > 400) and will be refit with increased iterations.\n", sep = "")

if (nrow(failing) == 0) {
  cat("Nothing to refit -- every existing fit already meets the criterion. Exiting.\n")
  quit(save = "no", status = 0)
}

# print.data.frame()'s internal print.default(m, ...) call errors with
# "invalid 'na.print' specification" if getOption("na.print") has been set
# to something invalid elsewhere in the session -- reset it to a known-good
# value just for this print call, then restore whatever it was.
old_na_print <- getOption("na.print")
options(na.print = "NA")
print(as.data.frame(failing %>% select(species, species_code, model, model_max_rhat, model_min_ess_bulk)))
options(na.print = old_na_print)

# Compile only the Stan model(s) actually needed for the failing tags -------
tags_needed  <- unique(failing$model)
model_base   <- if ("base" %in% tags_needed) {
  cmdstan_model(here::here("models", "slope_iCAR_route_NB_New.stan"), stanc_options = list("O1"))
} else NULL
model_single <- if (any(tags_needed != "base")) {
  cmdstan_model(here::here("models", "slope_iCAR_route_NB_New_covariate.stan"), stanc_options = list("O1"))
} else NULL

# ==========================================================================
# Step 2: refit each failing species/tag -------------------------------------
# ==========================================================================
results_list <- list()
refit_diagnostics_list <- list()

# Group failing rows by species so prepare_species_data() (expensive: BBS
# pull + Voronoi neighbours) runs once per species even if BOTH its tags
# failed, not twice ----------------------------------------------------------
failing_species <- failing %>% distinct(species, species_code) %>% arrange(species)

for (i in seq_len(nrow(failing_species))) {
  sp      <- failing_species$species[i]
  sp_code <- failing_species$species_code[i]
  sp_f    <- species_to_f(sp)
  sp_row  <- target_spp %>% filter(Common.Name == sp, Code == sp_code)
  sp_bbs   <- sp_row$bbs_english[1]
  sp_group <- sp_row$Group[1]

  tags_for_this_species <- failing %>%
    filter(species == sp, species_code == sp_code) %>%
    pull(model)

  cat("\n================================================================\n")
  cat("  [", i, "/", nrow(failing_species), "]", sp, " (", sp_group, ") -- tags:",
      paste(tags_for_this_species, collapse = ", "), "\n")
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

  for (tag in tags_for_this_species) {

    # Resumability: re-check THIS species/tag's current status right before
    # refitting, in case an earlier (interrupted) run of this same script
    # already fixed it -----------------------------------------------------
    current <- model_convergence_check(sp_f, tag, firstYear, lastYear, rds_dir)
    if (!is.null(current) && isTRUE(current$model_converged)) {
      cat("  Skipping", tag, "-- already converged (fixed by a previous run of this script).\n")
      next
    }

    # Preserve the current (under-converged) fit before it gets overwritten
    archive_existing_fit(sp_f, tag, firstYear, lastYear)

    diagnostics <- tryCatch(
      fit_one_covariate_model(species = sp, species_f = sp_f, model_tag = tag,
                              prepped = prepped, firstYear = firstYear, lastYear = lastYear,
                              model_base = model_base, model_single = model_single,
                              rds_dir = rds_dir, route_info_dir = route_info_dir,
                              cmdstanr_output_dir = cmdstanr_output_dir,
                              stan_data_dir = stan_data_dir,
                              chains = 4, iter_warmup = refit_iter_warmup,
                              iter_sampling = refit_iter_sampling,
                              adapt_delta = refit_adapt_delta,
                              max_treedepth = refit_max_treedepth,
                              show_exceptions = refit_show_exceptions,
                              overwrite_route_info = TRUE),
      error = function(e) {
        message("  [ERROR] Refit '", tag, "' failed for ", sp, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(diagnostics)) next

    diagnostics$group <- sp_group
    if (!is.null(current)) {
      diagnostics$before_max_rhat     <- current$model_max_rhat
      diagnostics$before_min_ess_bulk <- current$model_min_ess_bulk
    }

    results_list[[paste(sp, tag, sep = " | ")]] <- TRUE
    refit_diagnostics_list[[paste(sp, tag, sep = " | ")]] <- diagnostics

    cat("  Done — refit", tag, "for", sp, "\n")
  }

  if (length(refit_diagnostics_list) > 0 && i %% 10 == 0) {
    diag_csv <- file.path(output_dir,
                          paste0("diagnostics_refit_nonconverged_", firstYear, "_", lastYear, ".csv"))
    write.csv(bind_rows(refit_diagnostics_list), diag_csv, row.names = FALSE)
    cat("  [checkpoint] refit diagnostics written to:", diag_csv,
        "(", i, "/", nrow(failing_species), "species processed )\n")
  }
}

# Final refit diagnostics CSV -------------------------------------------------
if (length(refit_diagnostics_list) > 0) {
  diag_csv <- file.path(output_dir,
                        paste0("diagnostics_refit_nonconverged_", firstYear, "_", lastYear, ".csv"))
  write.csv(bind_rows(refit_diagnostics_list), diag_csv, row.names = FALSE)
  cat("\nRefit diagnostics written to:", diag_csv, "\n")
}

# ==========================================================================
# Step 3: re-check convergence after refitting, so you can see immediately
# whether it worked -----------------------------------------------------------
# ==========================================================================
cat("\n=== Re-checking convergence after refit ===\n")
convergence_after <- model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                                             firstYear = firstYear, lastYear = lastYear,
                                             rds_dir = rds_dir,
                                             bird_group = "all_species_anthro_postrefit")

still_failing <- convergence_after %>% filter(!model_converged)
cat("\n", nrow(still_failing), "species/tag combo(s) still fail the convergence criterion after refit.\n", sep = "")
if (nrow(still_failing) > 0) {
  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(still_failing %>% select(species, species_code, model, model_max_rhat, model_min_ess_bulk)))
  options(na.print = old_na_print)
  cat("\nThese may need an even longer refit, a different adapt_delta, or individual investigation",
      "-- e.g. a chain that's outright crashing (see tests/test_green_heron_debug.R) is a different",
      "failure mode than slow mixing, and more iterations alone won't necessarily fix it.\n")
}

cat("\n\nIMPORTANT: re-run 2c_generate_route_trend_csvs_covariates.R next so these improved fits",
    "propagate into all_route_trends_*.csv / all_model_level_summary_*.csv.\n")

cat("\n=== DONE ===\n")
cat("Species/tag combos refit this run:", length(results_list), "\n")
