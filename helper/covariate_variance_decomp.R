## =============================================================================
## helper/covariate_variance_decomp.R
##
## For a route-year land-cover covariate (the same CSV shape used throughout
## this project: CountryNum/StateNum/Route/.../year/<value column>), decompose
## its total variance into a BETWEEN-route piece (differences in each route's
## average value -- already absorbed by the route-specific intercept alpha[r]
## in the iCAR models) and a WITHIN-route piece (how much a route's own value
## moves from year to year -- the only part gamma1 actually has to work with,
## since alpha[r] soaks up everything time-invariant per route).
##
## Motivation: a covariate with almost all of its variance between routes
## (ICC close to 1) is close to time-invariant at any given route over the
## study window, so a route-intercept model has very little power to
## estimate its effect on trend, no matter how ecologically meaningful the
## covariate is cross-sectionally. This is the diagnostic used to compare
## "aridlands"/"Anthro"-type level covariates (found to have ICC ~0.995-0.9996,
## i.e. almost no within-route movement over 2010-2025) against the
## "*_to_anthro" transition covariates (found to have ICC ~0.36-0.38, i.e.
## most of their variance IS within-route/temporal) -- see
## output/files/aridlands_2010_2025.txt for the worked example.
##
## Nothing here is hardcoded to aridlands (or any particular covariate list)
## -- both functions below take every setting as an argument, same
## soft-coding pattern as helper/gamma_lookup.R.
##
## Two ways to use this:
##
##   1. Standalone / interactive: source this file, then call
##      covariate_variance_decomp() for one covariate or
##      covariate_variance_decomp_table() for several at once, e.g.:
##        source(here::here("helper", "covariate_variance_decomp.R"))
##        covariate_variance_decomp_table(list(
##          anthro   = list(file = "Anthro.csv", value_col = "Anthro"),
##          habitat1 = list(file = "aridland habitat1.csv", value_col = "aridland habitat1")
##        ))
##
##   2. Insert into a 1c-style script that already has a `covariates` (or
##      `covariate_lookups`) named list built from group_config -- source
##      this file (it only defines functions, nothing runs automatically),
##      then call covariate_variance_decomp_table() directly with that same
##      list before fitting, as a quick check of which covariates actually
##      have enough within-route movement to be worth fitting.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/covariate_variance_decomp.R")

#' Within/between-route variance decomposition for one route-year covariate.
#'
#' @param file covariate CSV filename, read from data_dir (same shape as
#'   every other covariate CSV in this project: CountryNum/StateNum/Route/
#'   .../year/<value_col>)
#' @param value_col exact column name to decompose (case/spacing must match
#'   the CSV header exactly -- check.names = FALSE is used, so a column like
#'   "aridland habitat to anthro1" must be passed with its literal spaces)
#' @param firstYear,lastYear restrict to this year range (should match
#'   whatever the fitting script uses)
#' @param data_dir directory containing the covariate CSV (default data/)
#' @return invisibly, a list with the per-route decomposition (`by_route`)
#'   and the summary row (`summary`)
covariate_variance_decomp <- function(file, value_col, firstYear = 2010, lastYear = 2025,
                                      data_dir = here::here("data")) {

  dat <- read.csv(file.path(data_dir, file), stringsAsFactors = FALSE, check.names = FALSE) %>%
    transmute(route_key = paste(StateNum, Route, sep = "-"),
              year      = year,
              val       = .data[[value_col]]) %>%
    filter(!is.na(val), year >= firstYear, year <= lastYear)

  by_route <- dat %>%
    group_by(route_key) %>%
    filter(n() > 1) %>%          # need >= 2 years to have a within-route variance
    summarise(within_var  = var(val),
              route_mean  = mean(val),
              route_range = max(val) - min(val),
              n_years     = n(),
              .groups = "drop")

  between_var <- var(by_route$route_mean)
  within_var  <- mean(by_route$within_var, na.rm = TRUE)
  icc         <- between_var / (between_var + within_var)

  summ <- data.frame(
    value_col                 = value_col,
    file                      = file,
    n_routes_multiyear        = nrow(by_route),
    mean_level                = mean(dat$val),
    between_route_var         = between_var,
    avg_within_route_var      = within_var,
    icc_between_share         = icc,
    median_within_route_range = median(by_route$route_range)
  )

  invisible(list(by_route = by_route, summary = summ))
}

#' Run covariate_variance_decomp() across several covariates and print/return
#' a combined, ICC-sorted summary table -- so covariates with little
#' within-route temporal signal (high ICC, close to 1) are easy to spot
#' against ones with real within-route movement (lower ICC).
#'
#' @param covariate_specs named list where each element is
#'   list(file = "...", value_col = "..."), same shape as a group_config
#'   entry's `covariates` list in the 1c scripts (the `rescale`/other keys,
#'   if present, are ignored here)
#' @param firstYear,lastYear,data_dir passed through to
#'   covariate_variance_decomp()
#' @param out_dir where to write the combined CSV (default output/files)
#' @param label used only to name the output CSV (e.g. "aridlands")
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the combined summary table (invisibly)
covariate_variance_decomp_table <- function(covariate_specs, firstYear = 2010, lastYear = 2025,
                                            data_dir = here::here("data"),
                                            out_dir = here::here("output", "files"),
                                            label = "covariates",
                                            write_csv = TRUE) {

  cat("=== within/between-route variance decomposition ===\n")
  cat("Period:", firstYear, "-", lastYear, " | covariates:",
      paste(names(covariate_specs), collapse = ", "), "\n\n")

  rows <- list()
  for (cn in names(covariate_specs)) {
    spec <- covariate_specs[[cn]]
    res  <- covariate_variance_decomp(spec$file, spec$value_col, firstYear, lastYear, data_dir)
    res$summary$covariate <- cn
    rows[[cn]] <- res$summary
  }

  out <- bind_rows(rows) %>%
    select(covariate, value_col, file, everything()) %>%
    arrange(desc(icc_between_share))

  cat("Sorted worst (most between-route/static) to best (most within-route/temporal):\n")
  print(out, row.names = FALSE)

  cat("\nRule of thumb: ICC close to 1 means the covariate barely moves within a route\n",
      "over the study window, so a route-intercept model has little power to estimate\n",
      "its effect on trend, whatever its cross-sectional/ecological relevance. Lower ICC\n",
      "means more real within-route temporal signal for gamma1 to use.\n", sep = "")

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("covariate_variance_decomp_", label, "_",
                                         firstYear, "_", lastYear, ".csv"))
    write.csv(out, out_csv, row.names = FALSE)
    cat("\nTable written to:", out_csv, "\n")
  }

  invisible(out)
}
