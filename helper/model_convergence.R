## =============================================================================
## helper/model_convergence.R
##
## Whole-model convergence check across EVERY species x EVERY model tag,
## INCLUDING "base" -- unlike helper/gamma_lookup.R's gamma_lookup_table(),
## which only covers non-"base" tags because it's structurally keyed off
## gamma1's own posterior row (base has no gamma1, so it's silently skipped
## there).
##
## Criterion (per user spec): "species are accepted in models converged
## based on Rhat < 1.01 and bulk effective sample sizes > 400" -- applied to
## the WORST rhat and LOWEST bulk ESS found across EVERY parameter in a fit
## (alpha, beta, gamma1 if present, everything), not just one parameter. A
## species/tag fails if ANY parameter is still poorly converged.
##
## Reads each *_summ_fit.rds file DIRECTLY off disk -- deliberately NOT the
## diagnostics_covariates_all_species_anthro_*.csv that
## 1c_species_iCAR_covariates.R writes. That CSV only logs species actually
## fit during the CURRENT run session (diagnostics_list starts empty every
## time 1c runs, and species skipped because they were already fit in an
## earlier session never get appended -- see 1c's main loop). Since 1c is
## explicitly designed to be interrupted and resumed across many sessions
## for a ~600-species run, that CSV is very likely incomplete as a record of
## "every species that's ever been fit." Reading summ_fit.rds files directly
## is complete regardless of how many separate sessions produced the
## current set of saved fits.
##
## Two ways to use this (same pattern as helper/gamma_lookup.R):
##   1. Insert into 1c_species_iCAR_covariates.R: add
##        source(here::here("helper", "model_convergence.R"))
##      after the main loop, once target_spp/model_tags/firstYear/lastYear/
##      rds_dir are in scope -- auto-runs using 1c's own live settings.
##   2. Standalone: Rscript helper/model_convergence.R (or source() with
##      nothing pre-set) falls back to this project's current full-species,
##      base+anthro defaults. Set `model_tags <- c(...)` / `target_spp <-
##      ...` before sourcing to check a different set.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/model_convergence.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Read one species/model's *_summ_fit.rds and compute WHOLE-MODEL
#' convergence diagnostics -- max rhat and min bulk ESS across EVERY
#' parameter in the fit -- independent of whether the model has a gamma1
#' parameter at all. Works for ANY model tag, including "base". Returns NULL
#' (with a message) if the file is missing, instead of erroring, so a batch
#' loop can skip cleanly.
#'
#' @param species_f file-safe species name, e.g. "Green_Heron" (matches the
#'   species_to_f() naming used by 1c/2c)
#' @param tag model tag, e.g. "base" or "anthro"
#' @param firstYear,lastYear must match the summ_fit.rds filename
#' @param rds_dir directory containing the *_summ_fit.rds files
#' @return a one-row data.frame (species_f/model/file/model_max_rhat/
#'   model_min_ess_bulk/model_converged), or NULL if the file is missing
model_convergence_check <- function(species_f, tag, firstYear = 2010, lastYear = 2025,
                                    rds_dir = here::here("output", "rds")) {
  summ_file <- file.path(rds_dir,
                         paste0(species_f, "_iCAR_", tag, "_", firstYear, "_", lastYear, "_summ_fit.rds"))

  if (!file.exists(summ_file)) {
    message("  [MISSING] ", basename(summ_file))
    return(NULL)
  }

  summ <- readRDS(summ_file)

  model_max_rhat     <- suppressWarnings(max(summ$rhat, na.rm = TRUE))
  model_min_ess_bulk <- suppressWarnings(min(summ$ess_bulk, na.rm = TRUE))
  model_converged <- is.finite(model_max_rhat) && is.finite(model_min_ess_bulk) &&
    model_max_rhat < 1.01 && model_min_ess_bulk > 400

  data.frame(
    species_f          = species_f,
    model              = tag,
    file               = basename(summ_file),
    model_max_rhat     = model_max_rhat,
    model_min_ess_bulk = model_min_ess_bulk,
    model_converged    = model_converged,
    stringsAsFactors   = FALSE
  )
}

#' Build the combined whole-model convergence table across EVERY species x
#' EVERY model tag (base included), write it to
#' output/files/model_convergence_<bird_group>_<firstYear>_<lastYear>.csv,
#' and print a pass/fail summary. Fully soft-coded: takes target_spp/
#' model_tags/firstYear/lastYear/rds_dir/bird_group as arguments instead of
#' hardcoding or re-deriving them, matching the same convention
#' helper/gamma_lookup.R uses.
#'
#' @param target_spp data.frame with at least Common.Name, Code columns
#' @param model_tags character vector of ALL model tags to check -- pass
#'   1c's own model_tags as-is, base included, no filtering needed here
#'   (unlike gamma_lookup_table(), which drops "base" itself)
#' @param firstYear,lastYear,rds_dir must match how the *_summ_fit.rds files
#'   were named/saved
#' @param bird_group used only to label the output CSV filename
#' @param out_dir where to write the combined CSV (default output/files)
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the combined convergence table (invisibly), or NULL if nothing found
model_convergence_table <- function(target_spp, model_tags, firstYear, lastYear,
                                    rds_dir, bird_group,
                                    out_dir = here::here("output", "files"),
                                    write_csv = TRUE) {

  cat("=== Whole-model convergence check (Rhat < 1.01 and bulk ESS > 400, ALL tags) ===\n")
  cat("Group:", bird_group, " | Period:", firstYear, "-", lastYear, "\n")
  cat("Species (n =", nrow(target_spp), ") x models:", paste(model_tags, collapse = ", "), "\n\n")

  rows <- list()
  for (i in seq_len(nrow(target_spp))) {
    sp      <- target_spp$Common.Name[i]
    sp_f    <- species_to_f(sp)
    sp_code <- target_spp$Code[i]

    for (tag in model_tags) {
      r <- model_convergence_check(sp_f, tag, firstYear, lastYear, rds_dir)
      if (!is.null(r)) {
        r$species      <- sp
        r$species_code <- sp_code
        rows[[paste(sp, tag, sep = " | ")]] <- r
      }
    }
  }

  if (length(rows) == 0) {
    message("No convergence rows found -- check that ", rds_dir,
            "/*_summ_fit.rds exist for these species/model_tags (run ",
            "1c_species_iCAR_covariates.R first) and that model_tags match your setup.")
    return(invisible(NULL))
  }

  convergence_table <- bind_rows(rows) %>%
    select(species, species_code, model, model_max_rhat, model_min_ess_bulk,
           model_converged, everything(), -species_f, -file)

  cat("\n=== Whole-model convergence summary, by model ===\n")
  print(convergence_table %>%
          group_by(model) %>%
          summarise(n_species   = n(),
                    n_converged = sum(model_converged, na.rm = TRUE),
                    n_failed    = sum(!model_converged, na.rm = TRUE),
                    .groups = "drop"))

  cat("\nSpecies FAILING whole-model convergence, by model:\n")
  # print.data.frame()'s internal print.default(m, ...) call errors with
  # "invalid 'na.print' specification" if getOption("na.print") has been set
  # to something invalid elsewhere in the session (confirmed in practice --
  # not a bug in this table itself: the very similar all-numeric summary
  # table printed just above succeeds because tibble's own print method
  # doesn't touch na.print at all, but a plain data.frame's print does).
  # Reset it to a known-good value just for this print call, then restore
  # whatever it was.
  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(convergence_table %>%
          filter(!model_converged) %>%
          arrange(model, desc(model_max_rhat)) %>%
          select(species, species_code, model, model_max_rhat, model_min_ess_bulk)))
  options(na.print = old_na_print)

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("model_convergence_", bird_group, "_",
                                         firstYear, "_", lastYear, ".csv"))
    write.csv(convergence_table, out_csv, row.names = FALSE)
    cat("\nFull table written to:", out_csv, "\n")
  }

  invisible(convergence_table)
}

## ==========================================================================
## Auto-run: call model_convergence_table() using whatever target_spp/
## model_tags/firstYear/lastYear/rds_dir/bird_group ALREADY EXIST in the
## calling environment (e.g. 1c_species_iCAR_covariates.R's, if this file is
## source()'d from the end of a 1c run) -- filling in this project's current
## full-species, base+anthro defaults via `if (!exists(...))` for anything
## not already set, so this also works standalone (Rscript
## helper/model_convergence.R) exactly as before.
## Set `model_convergence_skip_autorun <- TRUE` before sourcing to load just
## the two functions above without running anything.
## ==========================================================================
if (!exists("model_convergence_skip_autorun") || !isTRUE(model_convergence_skip_autorun)) {

  if (!exists("bird_group")) bird_group <- "all_species_anthro"
  if (!exists("firstYear"))  firstYear  <- 2010
  if (!exists("lastYear"))   lastYear   <- 2025
  if (!exists("model_tags")) model_tags <- c("base", "anthro")
  if (!exists("rds_dir"))    rds_dir    <- here::here("output", "rds")

  if (!exists("target_spp")) {
    spp_df_mc <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                          stringsAsFactors = FALSE)
    target_spp <- spp_df_mc %>%
      filter(in_bbs == TRUE) %>%
      distinct(Common.Name, Code, .keep_all = TRUE) %>%
      arrange(Common.Name)
  }

  model_convergence_table(target_spp = target_spp, model_tags = model_tags,
                          firstYear = firstYear, lastYear = lastYear,
                          rds_dir = rds_dir, bird_group = bird_group)
}
