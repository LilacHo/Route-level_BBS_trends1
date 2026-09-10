## =============================================================================
## helper/compare_credible_intervals.R
##
## Compares 90% vs 95% credible intervals -- for gamma1 (the covariate
## effect) AND for every route's own alpha[r]/beta[r] -- WITHOUT re-fitting
## anything.
##
## Why this needs the FULL posterior draws, not *_summ_fit.rds: summ_fit.rds
## only has whatever quantiles posterior::summarise_draws() computed by
## DEFAULT at fit time (q5/q95, a 90% CI -- see
## functions/covariate_model_fitting.R's `stanfit$summary()` call). There is
## no q2.5/q97.5 in there. Getting a 95% CI for anything already fit means
## re-deriving it from that species/tag's FULL posterior draws, saved
## separately as *_stanfit.rds (via stanfit$save_object()).
##
## This makes every read in this file much heavier than every other helper/
## script in this project (all of which only read the small summ_fit.rds):
## each *_stanfit.rds holds the complete draws for every parameter, not just
## a summary table. Expect this to be slow and memory-heavy across ~600
## species x 2 tags -- pass a filtered `target_spp` (e.g. just the species
## you actually care about) rather than running the full species list if you
## don't need everything at once.
##
## Two separate comparisons, since they answer different questions:
##
##   1. gamma1 (species/tag level, non-"base" tags only -- base has no
##      gamma1): does the "credible effect" call (90% CI excludes zero, the
##      criterion already used in helper/gamma_lookup.R and
##      2c_generate_route_trend_csvs_covariates.R's gamma1_excludes_zero)
##      change if you use a 95% CI instead? -> compare_gamma1_ci_table()
##
##   2. Per-route alpha[r]/beta[r] (every tag, including "base"): same 90%-
##      vs-95% question, but per route. The "excludes zero" / conclusion-
##      change flag is only computed for BETA (the trend) -- that's the
##      parameter where "credible effect" is a meaningful question, same as
##      gamma1. Alpha is a log-abundance intercept, not meaningfully
##      centered at zero, so no excludes_zero/flip flag is computed for it;
##      its raw 90%/95% CIs are still reported for completeness.
##      -> compare_route_ci_table()
##
## Both write a "did the credibility conclusion flip between 90% and 95%"
## flag, print a summary of how many species/routes flip, and are
## checkpointed/resumable (periodic CSV writes) given how slow this can be
## over a large species list. Both functions in this file read each
## species/tag's *_stanfit.rds AT MOST ONCE per call to
## compare_credible_intervals_table() (the combined driver below), so
## running gamma1 + route comparisons together doesn't pay the I/O cost
## twice for the same file.
##
## Two ways to use this (same soft-coded pattern as the rest of helper/):
##   1. Source with `compare_credible_intervals_skip_autorun <- TRUE` set
##      first, then call compare_credible_intervals_table(...) (or the two
##      individual table functions) directly with your own target_spp/
##      model_tags/etc.
##   2. Standalone: Rscript helper/compare_credible_intervals.R (or
##      source() with nothing pre-set) falls back to this project's current
##      full-species, base+anthro defaults -- likely to take a very long
##      time across ~600 species; consider setting `target_spp <- ...`
##      (a filtered subset) before sourcing instead.
## =============================================================================

library(dplyr)
library(posterior)
library(here)

here::i_am("helper/compare_credible_intervals.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Compute BOTH a 90% (q5/q95) and 95% (q2.5/q97.5) credible interval, plus
#' the mean, from a posterior::draws object for ONE OR MORE parameters (a
#' single parameter like "gamma1", or every indexed element of a route-level
#' parameter like "alpha"/"beta") in one vectorized call -- same quantile
#' machinery posterior::summarise_draws() uses internally by default (for
#' q5/q95), so these numbers are directly comparable to what's already in
#' *_summ_fit.rds.
#'
#' @return a data.frame with one row per variable: variable, mean, q2.5, q5,
#'   q95, q97.5
summarise_both_widths <- function(draws_obj) {
  posterior::summarise_draws(draws_obj, "mean",
                             ~posterior::quantile2(.x, probs = c(0.025, 0.05, 0.95, 0.975)))
}

#' TRUE if the interval [lo, hi] excludes zero (i.e. doesn't span it) --
#' same "credible effect" criterion already used for gamma1's 90% CI
#' elsewhere in this project (helper/gamma_lookup.R,
#' 2c_generate_route_trend_csvs_covariates.R's gamma1_excludes_zero).
ci_excludes_zero <- function(lo, hi) {
  !is.na(lo) & !is.na(hi) & (lo > 0 | hi < 0)
}

#' Extract the integer route index from a route-level Stan variable name,
#' e.g. "alpha[348]" -> 348L, "beta[27]" -> 27L.
extract_routeF <- function(varname) {
  m <- regmatches(varname, regexpr("\\[([0-9]+)\\]", varname))
  suppressWarnings(as.integer(gsub("\\[|\\]", "", m)))
}

#' Read ONE species/tag's FULL posterior (*_stanfit.rds) ONCE, and compute
#' both the gamma1 comparison (if this tag has gamma1) and the per-route
#' alpha/beta comparison from it. Returns NULL (with a message) if the
#' *_stanfit.rds file is missing, so a batch loop can skip cleanly.
#'
#' @param species_f file-safe species name
#' @param tag model tag, e.g. "base" or "anthro"
#' @param firstYear,lastYear must match the stanfit.rds filename
#' @param rds_dir directory containing the *_stanfit.rds files
#' @param route_info_dir directory containing the *_route_info.rds files
#'   (route, routeF, latitude, longitude) -- used to attach real route IDs
#'   to the per-route table; if missing, routeF alone is still returned
#' @return list(gamma1 = one-row data.frame or NULL, routes = data.frame),
#'   or NULL if the stanfit.rds itself is missing
compare_ci_for_combo <- function(species_f, tag, firstYear = 2010, lastYear = 2025,
                                 rds_dir = here::here("output", "rds"),
                                 route_info_dir = here::here("data", "route_info")) {

  stanfit_file <- file.path(rds_dir,
                            paste0(species_f, "_iCAR_", tag, "_", firstYear, "_", lastYear, "_stanfit.rds"))
  if (!file.exists(stanfit_file)) {
    message("  [MISSING stanfit] ", basename(stanfit_file))
    return(NULL)
  }

  stanfit <- readRDS(stanfit_file)

  # --- gamma1 (non-"base" tags only) -----------------------------------------
  gamma1_row <- NULL
  if (tag != "base") {
    gamma1_draws <- tryCatch(stanfit$draws(variables = "gamma1"), error = function(e) NULL)
    if (!is.null(gamma1_draws) && length(gamma1_draws) > 0) {
      g <- summarise_both_widths(gamma1_draws)
      gamma1_row <- data.frame(
        species_f              = species_f,
        model                  = tag,
        mean                   = g$mean[1],
        q2.5                   = g$q2.5[1],
        q5                     = g$q5[1],
        q95                    = g$q95[1],
        q97.5                  = g$q97.5[1],
        excludes_zero_90       = ci_excludes_zero(g$q5[1],   g$q95[1]),
        excludes_zero_95       = ci_excludes_zero(g$q2.5[1], g$q97.5[1]),
        stringsAsFactors       = FALSE
      )
      gamma1_row$ci_conclusion_changes <- gamma1_row$excludes_zero_90 != gamma1_row$excludes_zero_95
    } else {
      message("  [NO gamma1] ", basename(stanfit_file))
    }
  }

  # --- Per-route alpha[r]/beta[r] (every tag) ---------------------------------
  alpha_draws <- tryCatch(stanfit$draws(variables = "alpha"), error = function(e) NULL)
  beta_draws  <- tryCatch(stanfit$draws(variables = "beta"),  error = function(e) NULL)

  route_rows <- NULL
  if (!is.null(alpha_draws) && !is.null(beta_draws)) {
    alpha_summ <- summarise_both_widths(alpha_draws) %>%
      transmute(routeF      = extract_routeF(variable),
               alpha_mean  = mean,
               alpha_q2.5  = q2.5,
               alpha_q5    = q5,
               alpha_q95   = q95,
               alpha_q97.5 = q97.5)

    beta_summ <- summarise_both_widths(beta_draws) %>%
      transmute(routeF                    = extract_routeF(variable),
               beta_mean                 = mean,
               beta_q2.5                 = q2.5,
               beta_q5                   = q5,
               beta_q95                  = q95,
               beta_q97.5                = q97.5,
               beta_excludes_zero_90     = ci_excludes_zero(beta_q5,   beta_q95),
               beta_excludes_zero_95     = ci_excludes_zero(beta_q2.5, beta_q97.5),
               beta_ci_conclusion_changes = beta_excludes_zero_90 != beta_excludes_zero_95)

    route_info_file <- file.path(route_info_dir,
                                 paste0(species_f, "_", tag, "_", firstYear, "_", lastYear, "_route_info.rds"))
    route_info <- if (file.exists(route_info_file)) readRDS(route_info_file) else NULL

    route_rows <- alpha_summ %>%
      left_join(beta_summ, by = "routeF")
    if (!is.null(route_info)) {
      route_rows <- route_rows %>% left_join(route_info, by = "routeF")
    }
    route_rows <- route_rows %>%
      mutate(species_f = species_f, model = tag) %>%
      arrange(routeF)
  } else {
    message("  [NO alpha/beta draws] ", basename(stanfit_file))
  }

  list(gamma1 = gamma1_row, routes = route_rows)
}

#' Build BOTH combined tables (gamma1-level and route-level) across every
#' species x model_tag, reading each *_stanfit.rds only once. This is the
#' function you actually want to call for "compare 90% vs 95%" -- the two
#' individual table functions below are thin wrappers around this for when
#' you only want one or the other.
#'
#' @param target_spp data.frame with at least Common.Name, Code columns
#' @param model_tags character vector of model tags to check, e.g.
#'   c("base", "anthro") -- "base" is automatically skipped for the gamma1
#'   table (no gamma1 there) but still included for the route table
#' @param firstYear,lastYear,rds_dir,route_info_dir must match how the
#'   *_stanfit.rds/*_route_info.rds files were saved
#' @param bird_group used only to label the output CSV filenames
#' @param out_dir where to write the combined CSVs (default output/files)
#' @param write_csv if FALSE, skip writing CSVs (just return the tables)
#' @param checkpoint_every write/refresh the CSVs every N species (default
#'   10 -- lower than other helper/ scripts' 25, since this is much slower
#'   per species)
#' @return list(gamma1 = combined gamma1 table, routes = combined route
#'   table), invisibly
compare_credible_intervals_table <- function(target_spp, model_tags, firstYear, lastYear,
                                             rds_dir, route_info_dir, bird_group,
                                             out_dir = here::here("output", "files"),
                                             write_csv = TRUE, checkpoint_every = 10) {

  cat("=== Comparing 90% vs 95% credible intervals (gamma1 + per-route alpha/beta) ===\n")
  cat("Group:", bird_group, " | Period:", firstYear, "-", lastYear, "\n")
  cat("Species (n =", nrow(target_spp), ") x models:", paste(model_tags, collapse = ", "), "\n")
  cat("NOTE: this reads each species/tag's FULL *_stanfit.rds (not the small summ_fit.rds) --\n")
  cat("      expect this to be slow. Species/tags without a *_stanfit.rds (not yet fit by\n")
  cat("      1c_species_iCAR_covariates.R) are skipped up front, below. This script does NOT\n")
  cat("      skip species based on ITS OWN previous output though (unlike 1c/1d's skip-if-\n")
  cat("      already-fit) -- re-running re-derives every already-fit combo again from scratch,\n")
  cat("      so pass a smaller target_spp if you only need a subset.\n\n")

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  gamma1_csv <- file.path(out_dir, paste0("gamma1_ci_comparison_", bird_group, "_", firstYear, "_", lastYear, ".csv"))
  route_csv  <- file.path(out_dir, paste0("route_ci_comparison_",  bird_group, "_", firstYear, "_", lastYear, ".csv"))

  # Build every species x tag combo up front and check which ones actually
  # have a *_stanfit.rds -- i.e. were actually fit by
  # 1c_species_iCAR_covariates.R -- BEFORE doing any expensive work. With a
  # ~600-species, resumable, multi-session run, it's completely normal for
  # most of the list to not be fit yet; this skips those cleanly and
  # quietly (one summary line) instead of looping over every species,
  # calling compare_ci_for_combo() for each, and printing a "[MISSING
  # stanfit]" message per miss.
  combos <- expand.grid(row = seq_len(nrow(target_spp)), tag = model_tags,
                        stringsAsFactors = FALSE) %>%
    mutate(species      = target_spp$Common.Name[row],
           species_code = target_spp$Code[row],
           species_f    = species_to_f(species)) %>%
    mutate(stanfit_file = file.path(rds_dir, paste0(species_f, "_iCAR_", tag, "_",
                                                     firstYear, "_", lastYear, "_stanfit.rds")),
           has_stanfit  = file.exists(stanfit_file)) %>%
    select(-row)

  n_total <- nrow(combos)
  n_found <- sum(combos$has_stanfit)
  cat(n_found, "of", n_total, "species x tag combo(s) have a *_stanfit.rds (fit by 1c) -- skipping the other ",
      n_total - n_found, " (not yet fit by 1c_species_iCAR_covariates.R).\n\n", sep = "")

  combos_found <- combos %>% filter(has_stanfit)

  if (nrow(combos_found) == 0) {
    message("No *_stanfit.rds files found at all in ", rds_dir, " for these species/model_tags -- ",
            "run 1c_species_iCAR_covariates.R first.")
    return(invisible(list(gamma1 = NULL, routes = NULL)))
  }

  gamma1_rows <- list()
  route_rows  <- list()

  for (i in seq_len(nrow(combos_found))) {
    sp      <- combos_found$species[i]
    sp_f    <- combos_found$species_f[i]
    sp_code <- combos_found$species_code[i]
    tag     <- combos_found$tag[i]

    if (i %% 10 == 0 || i == 1) cat("[", i, "/", nrow(combos_found), "]", sp, "|", tag, "\n")

    r <- compare_ci_for_combo(sp_f, tag, firstYear, lastYear, rds_dir, route_info_dir)
    if (is.null(r)) next   # defensive -- shouldn't happen given the has_stanfit filter above

    if (!is.null(r$gamma1)) {
      r$gamma1$species      <- sp
      r$gamma1$species_code <- sp_code
      gamma1_rows[[paste(sp, tag, sep = " | ")]] <- r$gamma1
    }
    if (!is.null(r$routes)) {
      r$routes$species      <- sp
      r$routes$species_code <- sp_code
      route_rows[[paste(sp, tag, sep = " | ")]] <- r$routes
    }

    if (i %% checkpoint_every == 0 && write_csv) {
      if (length(gamma1_rows) > 0) write.csv(bind_rows(gamma1_rows), gamma1_csv, row.names = FALSE)
      if (length(route_rows)  > 0) write.csv(bind_rows(route_rows),  route_csv,  row.names = FALSE)
      cat("  [checkpoint] written at", i, "/", nrow(combos_found), "combo(s)\n")
    }
  }

  gamma1_table <- if (length(gamma1_rows) > 0) {
    bind_rows(gamma1_rows) %>%
      select(species, species_code, model, mean, q2.5, q5, q95, q97.5,
             excludes_zero_90, excludes_zero_95, ci_conclusion_changes, everything(), -species_f)
  } else NULL

  route_table <- if (length(route_rows) > 0) {
    bind_rows(route_rows) %>%
      select(species, species_code, model, routeF, any_of("route"),
             any_of("latitude"), any_of("longitude"),
             alpha_mean, alpha_q2.5, alpha_q5, alpha_q95, alpha_q97.5,
             beta_mean, beta_q2.5, beta_q5, beta_q95, beta_q97.5,
             beta_excludes_zero_90, beta_excludes_zero_95, beta_ci_conclusion_changes,
             everything(), -species_f)
  } else NULL

  # --- gamma1 summary ---------------------------------------------------------
  if (!is.null(gamma1_table)) {
    cat("\n=== gamma1: 90% vs 95% CI comparison ===\n")
    old_na_print <- getOption("na.print")
    options(na.print = "NA")
    print(as.data.frame(gamma1_table))
    options(na.print = old_na_print)

    n_flip <- sum(gamma1_table$ci_conclusion_changes, na.rm = TRUE)
    cat("\n", n_flip, "of", nrow(gamma1_table),
        "species/tag combo(s) FLIP their 'credible effect' conclusion when switching from a 90% to a 95% CI.\n", sep = "")
    if (n_flip > 0) {
      cat("(These were credible at 90% but NOT at 95% -- a 95% CI is wider, so this is the only",
          "direction a flip can go.)\n")
      old_na_print <- getOption("na.print")
      options(na.print = "NA")
      print(as.data.frame(gamma1_table %>% filter(ci_conclusion_changes) %>%
                            select(species, species_code, model, q5, q95, q2.5, q97.5)))
      options(na.print = old_na_print)
    }

    if (write_csv) {
      write.csv(gamma1_table, gamma1_csv, row.names = FALSE)
      cat("\ngamma1 comparison written to:", gamma1_csv, "\n")
    }
  } else {
    message("No gamma1 rows produced -- check that *_stanfit.rds files exist for these ",
            "species/model_tags (run 1c_species_iCAR_covariates.R first) and that at least ",
            "one non-'base' tag is in model_tags.")
  }

  # --- route summary -----------------------------------------------------------
  if (!is.null(route_table)) {
    n_route_flip <- sum(route_table$beta_ci_conclusion_changes, na.rm = TRUE)
    cat("\n=== Per-route beta: 90% vs 95% CI comparison ===\n")
    cat(n_route_flip, "of", nrow(route_table),
        "route/species/tag row(s) FLIP beta's 'credible trend' conclusion when switching from a 90% to a 95% CI.\n")
    cat("\nBy model:\n")
    print(route_table %>%
            group_by(model) %>%
            summarise(n_routes         = n(),
                      n_credible_90    = sum(beta_excludes_zero_90, na.rm = TRUE),
                      n_credible_95    = sum(beta_excludes_zero_95, na.rm = TRUE),
                      n_flip           = sum(beta_ci_conclusion_changes, na.rm = TRUE),
                      .groups = "drop"))

    if (write_csv) {
      write.csv(route_table, route_csv, row.names = FALSE)
      cat("\nPer-route comparison written to:", route_csv, "\n")
    }
  } else {
    message("No route rows produced -- check that *_stanfit.rds files exist for these species/model_tags.")
  }

  invisible(list(gamma1 = gamma1_table, routes = route_table))
}

#' Thin wrapper around compare_credible_intervals_table() for when you only
#' want the gamma1 comparison (still reads the full stanfit.rds -- there's
#' no cheaper way to get gamma1 alone once route-level data isn't needed,
#' but skipping the route extraction saves some compute per file).
compare_gamma1_ci_table <- function(target_spp, model_tags, firstYear, lastYear,
                                    rds_dir, bird_group,
                                    out_dir = here::here("output", "files"),
                                    write_csv = TRUE) {
  res <- compare_credible_intervals_table(target_spp = target_spp, model_tags = model_tags,
                                          firstYear = firstYear, lastYear = lastYear,
                                          rds_dir = rds_dir, route_info_dir = here::here("data", "route_info"),
                                          bird_group = bird_group, out_dir = out_dir, write_csv = write_csv)
  invisible(res$gamma1)
}

#' Thin wrapper around compare_credible_intervals_table() for when you only
#' want the per-route alpha/beta comparison.
compare_route_ci_table <- function(target_spp, model_tags, firstYear, lastYear,
                                   rds_dir, route_info_dir, bird_group,
                                   out_dir = here::here("output", "files"),
                                   write_csv = TRUE) {
  res <- compare_credible_intervals_table(target_spp = target_spp, model_tags = model_tags,
                                          firstYear = firstYear, lastYear = lastYear,
                                          rds_dir = rds_dir, route_info_dir = route_info_dir,
                                          bird_group = bird_group, out_dir = out_dir, write_csv = write_csv)
  invisible(res$routes)
}

## ==========================================================================
## Auto-run with this project's current full-species, base+anthro defaults --
## set `compare_credible_intervals_skip_autorun <- TRUE` before sourcing to
## load just the functions above without running anything. Given how slow
## this can be across ~600 species, strongly consider setting a filtered
## `target_spp` before sourcing standalone rather than using the full list.
## ==========================================================================
if (!exists("compare_credible_intervals_skip_autorun") || !isTRUE(compare_credible_intervals_skip_autorun)) {

  if (!exists("bird_group"))      bird_group      <- "all_species_anthro"
  if (!exists("firstYear"))       firstYear       <- 2010
  if (!exists("lastYear"))        lastYear        <- 2025
  if (!exists("model_tags"))      model_tags      <- c("base", "anthro")
  if (!exists("rds_dir"))         rds_dir         <- here::here("output", "rds")
  if (!exists("route_info_dir"))  route_info_dir  <- here::here("data", "route_info")

  if (!exists("target_spp")) {
    spp_df_ci <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                          stringsAsFactors = FALSE)
    target_spp <- spp_df_ci %>%
      filter(in_bbs == TRUE) %>%
      distinct(Common.Name, Code, .keep_all = TRUE) %>%
      arrange(Common.Name)
  }

  compare_credible_intervals_table(target_spp = target_spp, model_tags = model_tags,
                                   firstYear = firstYear, lastYear = lastYear,
                                   rds_dir = rds_dir, route_info_dir = route_info_dir,
                                   bird_group = bird_group)
}
