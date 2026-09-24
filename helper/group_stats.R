## =============================================================================
## helper/group_stats.R
##
## Descriptive statistics for gamma1 AGGREGATED BY BIRD GROUP, built on top
## of the per-species table helper/gamma_lookup.R produces (gamma_lookup.R
## is deliberately kept to per-species lookup only -- this file owns
## everything about combining species WITHIN a group).
##
## This is the reproducible replacement for the ad hoc, hand-generated
## "gamma1_ci_comparison_all_species_*.csv" files that
## 4d_visualization_publishable.R previously had to fall back on for
## per-species CI (those aren't produced by any script in this project and
## lack sd/median).
##
## For each model x group combination, gamma_group_summary() reports:
##   - n, mean_of_means, mean_of_medians, mean_median_gap (species-averaged
##     |mean - median|; a large gap flags skewed individual posteriors,
##     where a species' posterior MEAN and MEDIAN meaningfully disagree --
##     see the "MEAN vs. MEDIAN" note on gamma_group_summary() below)
##   - pct_credible_90 / pct_credible_95: % of species whose own 90%/95% CI
##     excludes zero (a CONSISTENCY/PREVALENCE measure, independent of
##     effect size -- different from the pooled means below)
##   - fe_mean/fe_ci_lo/fe_ci_hi: inverse-variance-weighted (fixed-effect)
##     pooled mean, assuming one true gamma1 per group
##   - re_mean/re_ci_lo/re_ci_hi/tau2/i2: DerSimonian-Laird random-effects
##     pooled mean, additionally allowing real between-species heterogeneity
##     (i2 = % of variance that's real heterogeneity rather than sampling
##     noise). This CI is never narrower than the fixed-effect one and is
##     the more honest summary whenever species genuinely differ (typical
##     for gamma1 -- check i2 before trusting fe_* over re_*).
## dl_random_effects() (the DerSimonian & Laird 1986 method-of-moments
## estimator) is implemented here in base R rather than via the metafor
## package, since nothing else in this project's toolchain currently
## depends on metafor.
##
## Two ways to use this:
##
##   1. Right after helper/gamma_lookup.R, in the same session: assign its
##      result first (gamma_lookup.R's own auto-run doesn't, since it's
##      meant to stay decoupled from this file) --
##        gamma_table_result <- gamma_lookup_table(target_spp = ..., ...)
##        source(here::here("helper", "group_stats.R"))
##      This file's auto-run picks up that live gamma_table_result directly
##      -- no disk round-trip needed.
##
##   2. Standalone: Rscript helper/group_stats.R (or source() with nothing
##      pre-set) reads the CSV gamma_lookup_table() already wrote to disk --
##      output/files/gamma_lookup_<run_label>_<firstYear>_<lastYear>.csv --
##      using the same `if (!exists(...))` override pattern used elsewhere
##      in this project. Defaults: run_label = "all_species", firstYear =
##      2010, lastYear = 2025 (matching gamma_lookup.R's own defaults) --
##      set `run_label <- "..."` / `firstYear <-` / `lastYear <-` before
##      sourcing to match a different gamma_lookup_*.csv.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/group_stats.R")

#' DerSimonian-Laird random-effects pooling for one set of estimates.
#'
#' Base-R implementation (no metafor dependency) of the standard
#' DerSimonian & Laird (1986) method-of-moments estimator. Used by
#' gamma_group_summary() below to pool species-level gamma1 estimates WITHIN
#' a bird group two ways: a fixed-effect mean (every species assumed to
#' share one true gamma1, so only each species' own sampling/estimation
#' variance matters) and a random-effects mean (species allowed to have
#' genuinely different true gamma1 values, with the extra between-species
#' variance tau^2 estimated from how much the species disagree beyond what
#' their individual SEs alone would predict). The random-effects CI is
#' never narrower than the fixed-effect one, and is the honest choice
#' whenever species-level heterogeneity is plausible (which it visibly is
#' for gamma1 -- see 4d_visualization_publishable.R's Fig 1 violins).
#'
#' @param yi numeric vector of point estimates (here, each species'
#'   posterior mean gamma1)
#' @param vi numeric vector of sampling VARIANCES (i.e. SE^2, or posterior
#'   sd^2), same length/order as yi
#' @return a one-row tibble: k, fe_mean, fe_se, re_mean, re_se, tau2, i2
#'   (i2 = % of total variance attributable to real between-species
#'   heterogeneity rather than each species' own sampling error; 0-100)
dl_random_effects <- function(yi, vi) {
  keep <- is.finite(yi) & is.finite(vi) & vi > 0
  yi <- yi[keep]
  vi <- vi[keep]
  k  <- length(yi)

  if (k == 0) {
    return(tibble(k = 0L, fe_mean = NA_real_, fe_se = NA_real_,
                  re_mean = NA_real_, re_se = NA_real_,
                  tau2 = NA_real_, i2 = NA_real_))
  }
  if (k == 1) {
    return(tibble(k = 1L, fe_mean = yi, fe_se = sqrt(vi),
                  re_mean = yi, re_se = sqrt(vi), tau2 = 0, i2 = 0))
  }

  w_fe    <- 1 / vi
  fe_mean <- sum(w_fe * yi) / sum(w_fe)
  fe_se   <- sqrt(1 / sum(w_fe))

  # Cochran's Q + method-of-moments tau^2 (DerSimonian-Laird), clamped at 0
  # -- a negative method-of-moments estimate just means the observed spread
  # is no more than sampling error alone would predict.
  Q    <- sum(w_fe * (yi - fe_mean)^2)
  df   <- k - 1
  C    <- sum(w_fe) - sum(w_fe^2) / sum(w_fe)
  tau2 <- max(0, (Q - df) / C)
  i2   <- max(0, 100 * (Q - df) / Q)

  w_re    <- 1 / (vi + tau2)
  re_mean <- sum(w_re * yi) / sum(w_re)
  re_se   <- sqrt(1 / sum(w_re))

  tibble(k = k, fe_mean = fe_mean, fe_se = fe_se,
         re_mean = re_mean, re_se = re_se, tau2 = tau2, i2 = i2)
}

#' Descriptive statistics for gamma1 BY BIRD GROUP, built on top of
#' helper/gamma_lookup.R's per-species output.
#'
#' On MEAN vs. MEDIAN: fe_mean/re_mean use each species' posterior MEAN (not
#' median) because inverse-variance weighting's optimality relies on a
#' roughly Gaussian per-species posterior, which pairs with mean + SE/sd,
#' not median. mean_median_gap is reported alongside so a large gap (skewed
#' posteriors, typically from sparse/weak per-species data) can be flagged
#' rather than silently trusted -- consider treating that species' pooled
#' contribution cautiously if its gap is large relative to its own sd.
#'
#' @param gamma_table output of gamma_lookup_table() (or any data.frame with
#'   at least model, group, mean, median, sd, q5, q95 columns) -- typically
#'   either the live object from gamma_lookup_table(), or the CSV it wrote
#'   read back with read.csv()
#' @param firstYear,lastYear used only to label the output CSV filename
#' @param out_dir where to write the summary CSV (default output/files)
#' @param run_label used only to label the output CSV filename (e.g.
#'   "all_species" or one bird_group)
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the group summary tibble (invisibly), or NULL if gamma_table has
#'   no usable rows/columns
gamma_group_summary <- function(gamma_table, firstYear, lastYear,
                                out_dir = here::here("output", "files"),
                                run_label = "all_species",
                                write_csv = TRUE) {

  required_cols <- c("model", "group", "mean", "median", "sd", "q5", "q95")
  missing_cols  <- setdiff(required_cols, names(gamma_table))
  if (length(missing_cols) > 0) {
    message("gamma_group_summary() skipped -- gamma_table is missing column(s): ",
            paste(missing_cols, collapse = ", "),
            ". Make sure target_spp passed to gamma_lookup_table() has a Group column.")
    return(invisible(NULL))
  }

  d <- gamma_table %>%
    filter(!is.na(group), !is.na(mean), !is.na(sd), sd > 0)

  if (nrow(d) == 0) {
    message("gamma_group_summary() skipped -- no rows with non-missing group/mean/sd.")
    return(invisible(NULL))
  }

  cat("\n=== gamma1 descriptive statistics BY BIRD GROUP ===\n")

  group_summary <- d %>%
    group_by(model, group) %>%
    summarise(
      n               = n(),
      mean_of_means   = mean(mean, na.rm = TRUE),
      mean_of_medians = mean(median, na.rm = TRUE),
      mean_median_gap = mean(abs(mean - median), na.rm = TRUE),
      pct_credible_90 = 100 * mean(!is.na(q5) & !is.na(q95) & (q5 > 0 | q95 < 0), na.rm = TRUE),
      pct_credible_95 = if (all(c("q2.5", "q97.5") %in% names(d))) {
        100 * mean(!is.na(q2.5) & !is.na(q97.5) & (q2.5 > 0 | q97.5 < 0), na.rm = TRUE)
      } else NA_real_,
      pooled          = list(dl_random_effects(mean, sd^2)),
      .groups = "drop"
    ) %>%
    tidyr::unnest(pooled) %>%
    select(-k) %>%
    mutate(
      fe_ci_lo = fe_mean - qnorm(0.95) * fe_se,   # 90% CI, matching the
      fe_ci_hi = fe_mean + qnorm(0.95) * fe_se,   # per-species q5/q95 convention
      re_ci_lo = re_mean - qnorm(0.95) * re_se,   # used throughout this project
      re_ci_hi = re_mean + qnorm(0.95) * re_se
    ) %>%
    select(model, group, n, mean_of_means, mean_of_medians, mean_median_gap,
           pct_credible_90, pct_credible_95,
           fe_mean, fe_ci_lo, fe_ci_hi,
           re_mean, re_ci_lo, re_ci_hi, tau2, i2) %>%
    arrange(model, desc(mean_of_means))

  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(group_summary))
  options(na.print = old_na_print)

  cat("\nNote: fe_* assumes one true gamma1 per group (species differ only by",
      "sampling noise); re_* (DerSimonian-Laird random-effects) additionally",
      "allows real between-species heterogeneity -- prefer re_* whenever i2 is",
      "not small (i2 = % of variance that's real heterogeneity, not noise).\n")

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("gamma_group_summary_", run_label, "_",
                                         firstYear, "_", lastYear, ".csv"))
    write.csv(group_summary, out_csv, row.names = FALSE)
    cat("\nGroup summary written to:", out_csv, "\n")
  }

  invisible(group_summary)
}

## ==========================================================================
## Auto-run: call gamma_group_summary() using either a live gamma_table_result
## already in scope (e.g. just built by sourcing helper/gamma_lookup.R and
## assigning its result) or, failing that, the CSV helper/gamma_lookup.R
## already wrote to disk. Set `group_stats_skip_autorun <- TRUE` before
## sourcing to load just the two functions above without running anything.
## ==========================================================================
if (!exists("group_stats_skip_autorun") || !isTRUE(group_stats_skip_autorun)) {

  if (!exists("firstYear")) firstYear <- 2010
  if (!exists("lastYear"))  lastYear  <- 2025
  if (!exists("run_label")) run_label <- "all_species"

  if (!exists("gamma_table_result")) {
    gamma_csv <- here::here("output", "files",
                            paste0("gamma_lookup_", run_label, "_", firstYear, "_", lastYear, ".csv"))
    if (!file.exists(gamma_csv)) {
      stop("No per-species gamma table found at ", gamma_csv, " -- run helper/gamma_lookup.R ",
           "first with matching run_label/firstYear/lastYear (or bird_group, which becomes ",
           "run_label), or assign its result to gamma_table_result before sourcing this file.")
    }
    cat("Reading per-species gamma table from:", gamma_csv, "\n")
    gamma_table_result <- read.csv(gamma_csv, stringsAsFactors = FALSE)
  }

  gamma_group_summary(gamma_table_result, firstYear = firstYear, lastYear = lastYear,
                      run_label = run_label)
}
