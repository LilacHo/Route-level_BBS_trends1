## =============================================================================
## helper/gamma_lookup.R
##
## Pull the full posterior summary row for `gamma1` (mean, median, sd, mad,
## q5, q95, rhat, ess_bulk, ess_tail -- whatever posterior::summarise_draws()
## produced) out of the EXISTING *_summ_fit.rds files written by
## 1c_species_iCAR_covariates.R, for every species x covariate-model
## combination -- without re-running 1c (no re-fitting, no stanfit.rds
## loading; *_summ_fit.rds is a small tibble, not the full draws object).
##
## This is a lightweight look-and-compare tool, separate from
## 2c_generate_route_trend_csvs_covariates.R's main pipeline output: it's
## meant for quickly inspecting gamma1's full posterior summary (not just
## mean/median but also sd, rhat, ess_bulk/tail) across every species and
## model in one table, e.g. to decide mean-vs-median reporting or to spot
## poorly-converged gamma1 estimates, without touching 2c's CSVs.
##
## Nothing here is hardcoded to grasslands (or any particular bird_group/
## model_tags) -- both gamma_lookup() and gamma_lookup_table() below take
## every setting as an argument.
##
## Two ways to use this:
##
##   1. Insert into 1c_species_iCAR_covariates.R (or any script that already
##      has target_spp/model_tags/firstYear/lastYear/rds_dir/bird_group in
##      scope, which 1c does by the time its main loop finishes). Just add
##      one line at the end of 1c:
##        source(here::here("helper", "gamma_lookup.R"))
##      Sourcing this file automatically calls gamma_lookup_table() using
##      1c's OWN live target_spp/model_tags/bird_group/etc. -- no settings to
##      edit, and it can never drift out of sync with whatever 1c was
##      actually just run with (see the "Auto-run" block at the bottom).
##
##   2. Standalone: Rscript helper/gamma_lookup.R (or source() from the
##      console with nothing pre-set) falls back to reading
##      spp_names_codes_group_aou.csv itself with default settings
##      (bird_group = "grasslands", the grasslands 4-tag model_tags), using
##      the same `if (!exists(...))` override pattern already used elsewhere
##      in this project (species_filter in 1c/2c, model_tag in 4c) -- set
##      `bird_group <- "..."` / `model_tags <- c(...)` before sourcing to
##      run it for a different group.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/gamma_lookup.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Read one species/model's *_summ_fit.rds and pull out its gamma1 row
#' verbatim (every column posterior::summarise_draws() produced), tagged
#' with species/model/file for identification. Returns NULL (with a message)
#' if the file is missing or has no gamma1 row, instead of erroring, so a
#' batch loop can skip cleanly.
#'
#' @param species_f file-safe species name, e.g. "Bairds_Sparrow" (matches
#'   the species_to_f() naming used by 1c/2c)
#' @param tag covariate model tag, e.g. "grassland_habitat" or "wetlands"
#' @param firstYear,lastYear must match the summ_fit.rds filename
#' @param rds_dir directory containing the *_summ_fit.rds files
#' @return a one-row data.frame (gamma1's full summary row + species_f/tag/
#'   file columns), or NULL if unavailable
gamma_lookup <- function(species_f, tag, firstYear = 2010, lastYear = 2025,
                          rds_dir = here::here("output", "rds")) {
  summ_file <- file.path(rds_dir,
                         paste0(species_f, "_iCAR_", tag, "_", firstYear, "_", lastYear, "_summ_fit.rds"))

  if (!file.exists(summ_file)) {
    message("  [MISSING] ", basename(summ_file))
    return(NULL)
  }

  summ <- readRDS(summ_file)

  gamma_row <- summ[summ$variable == "gamma1", ]
  if (nrow(gamma_row) == 0) {
    message("  [NO gamma1] ", basename(summ_file))
    return(NULL)
  }

  gamma_row$species_f <- species_f
  gamma_row$model     <- tag
  gamma_row$file       <- basename(summ_file)

  # species_f/model/file first, then whatever summarise_draws() columns exist
  gamma_row %>% select(species_f, model, file, everything(), -variable)
}

#' Build the combined gamma1 table across every species x model_tag
#' combination, write it to
#' output/files/gamma_lookup_<bird_group>_<firstYear>_<lastYear>.csv, and
#' print a quick convergence/credibility summary. Fully soft-coded: takes
#' target_spp/model_tags/firstYear/lastYear/rds_dir/bird_group as arguments
#' instead of hardcoding or re-deriving them, so it can be called directly
#' with whatever a caller (e.g. 1c_species_iCAR_covariates.R) already built,
#' for any land cover group / covariate-tag set.
#'
#' @param target_spp data.frame with at least Common.Name, Code columns (as
#'   built by 1c_species_iCAR_covariates.R from
#'   data/spp_names_codes_group_aou.csv)
#' @param model_tags character vector of model tags to look up. "base" is
#'   automatically dropped if present (it has no gamma1) -- pass 1c's own
#'   model_tags as-is, no need to filter it yourself.
#' @param firstYear,lastYear,rds_dir must match how the *_summ_fit.rds files
#'   were named/saved
#' @param bird_group used only to label the output CSV filename
#' @param out_dir where to write the combined CSV (default output/files)
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the combined gamma_table (invisibly), or NULL if nothing found
gamma_lookup_table <- function(target_spp, model_tags, firstYear, lastYear,
                               rds_dir, bird_group,
                               out_dir = here::here("output", "files"),
                               write_csv = TRUE) {

  model_tags_gamma <- setdiff(model_tags, "base")  # base has no gamma1

  cat("=== gamma1 lookup ===\n")
  cat("Group:", bird_group, " | Period:", firstYear, "-", lastYear, "\n")
  cat("Species (n =", nrow(target_spp), ") x models:", paste(model_tags_gamma, collapse = ", "), "\n\n")

  rows <- list()
  for (i in seq_len(nrow(target_spp))) {
    sp      <- target_spp$Common.Name[i]
    sp_f    <- species_to_f(sp)
    sp_code <- target_spp$Code[i]

    for (tag in model_tags_gamma) {
      r <- gamma_lookup(sp_f, tag, firstYear, lastYear, rds_dir)
      if (!is.null(r)) {
        r$species      <- sp
        r$species_code <- sp_code
        rows[[paste(sp, tag, sep = " | ")]] <- r
      }
    }
  }

  if (length(rows) == 0) {
    message("No gamma1 rows found -- check that ", rds_dir,
            "/*_summ_fit.rds exist for these species/model_tags (run ",
            "1c_species_iCAR_covariates.R first) and that bird_group/",
            "model_tags match your setup.")
    return(invisible(NULL))
  }

  gamma_table <- bind_rows(rows) %>%
    select(species, species_code, model, everything(), -species_f, -file)

  # Reorder columns so the diagnostics people usually eyeball first (mean,
  # median, q5, q95, rhat, ess_bulk, ess_tail) come right after model, if
  # present -- falls back gracefully if posterior::summarise_draws() used a
  # different default measure set.
  preferred_order <- c("species", "species_code", "model",
                       "mean", "median", "sd", "mad", "q5", "q95",
                       "rhat", "ess_bulk", "ess_tail")
  gamma_table <- gamma_table %>%
    select(any_of(preferred_order), everything())

  cat("\n=== gamma1 table ===\n")
  print(gamma_table, n = Inf, width = Inf)

  cat("\n=== Quick summary ===\n")
  if (all(c("mean", "median") %in% names(gamma_table))) {
    cat("\nMean vs. median gamma1, per model (species-averaged, and largest |mean-median| gap):\n")
    print(gamma_table %>%
            mutate(mean_median_gap = abs(mean - median)) %>%
            group_by(model) %>%
            summarise(n              = n(),
                      mean_of_means   = mean(mean, na.rm = TRUE),
                      mean_of_medians = mean(median, na.rm = TRUE),
                      max_gap         = max(mean_median_gap, na.rm = TRUE),
                      .groups = "drop"))
  }

  if (all(c("q5", "q95") %in% names(gamma_table))) {
    cat("\nSpecies with gamma1's 90% CI excluding zero (credible effect), by model:\n")
    print(gamma_table %>%
            mutate(excludes_zero = !is.na(q5) & !is.na(q95) & (q5 > 0 | q95 < 0)) %>%
            group_by(model) %>%
            summarise(n_species       = n(),
                      n_excludes_zero = sum(excludes_zero, na.rm = TRUE),
                      .groups = "drop"))
  }

  if ("rhat" %in% names(gamma_table)) {
    cat("\nSpecies with gamma1 rhat > 1.01 (possible convergence concern), by model:\n")
    print(gamma_table %>%
            filter(rhat > 1.01) %>%
            select(species, species_code, model, rhat, any_of(c("ess_bulk", "ess_tail"))))
  }

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("gamma_lookup_", bird_group, "_",
                                         firstYear, "_", lastYear, ".csv"))
    write.csv(gamma_table, out_csv, row.names = FALSE)
    cat("\nFull table written to:", out_csv, "\n")
  }

  invisible(gamma_table)
}

## ==========================================================================
## Auto-run: call gamma_lookup_table() using whatever target_spp/model_tags/
## firstYear/lastYear/rds_dir/bird_group ALREADY EXIST in the calling
## environment (e.g. 1c_species_iCAR_covariates.R's, if this file is
## source()'d from the end of a 1c run) -- filling in grasslands defaults
## via `if (!exists(...))` for any not already set, so this also still works
## standalone (Rscript helper/gamma_lookup.R) exactly as before.
## Set `gamma_lookup_skip_autorun <- TRUE` before sourcing to load just the
## two functions above without running anything.
## ==========================================================================
if (!exists("gamma_lookup_skip_autorun") || !isTRUE(gamma_lookup_skip_autorun)) {

  if (!exists("bird_group")) bird_group <- "grasslands"
  if (!exists("firstYear"))  firstYear  <- 2010
  if (!exists("lastYear"))   lastYear   <- 2025
  if (!exists("model_tags")) model_tags <- c("base", "grassland_habitat",
                                             "grassland_anthro",
                                             "grassland_habitat_to_anthro")
  if (!exists("rds_dir"))    rds_dir    <- here::here("output", "rds")

  if (!exists("target_spp")) {
    spp_df_gl <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                          stringsAsFactors = FALSE)
    target_spp <- spp_df_gl %>%
      filter(Group == bird_group, in_bbs == TRUE) %>%
      distinct(Common.Name, Code, .keep_all = TRUE) %>%
      arrange(Common.Name)
  }

  gamma_lookup_table(target_spp = target_spp, model_tags = model_tags,
                     firstYear = firstYear, lastYear = lastYear,
                     rds_dir = rds_dir, bird_group = bird_group)
}
