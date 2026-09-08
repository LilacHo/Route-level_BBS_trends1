## =============================================================================
## helper/translate_offending_routes.R
##
## Translates the ROUTE-LEVEL offending parameters found by
## helper/diagnose_nonconvergence.R (alpha_raw[r]/beta_raw_space[r], where r
## is the internal Stan route index "routeF" -- NOT a real BBS route ID)
## into the actual route IDs and coordinates, by joining against each
## species/tag's own data/route_info/*_route_info.rds (route, routeF,
## latitude, longitude).
##
## This is the concrete follow-through on the "flag/exclude at the route
## level, not the whole species" approach discussed in chat: rather than
## dropping an entire species (e.g. Sharp-shinned Hawk, Hooded Merganser)
## because a chunk of its routes' alpha/beta failed to converge, this tells
## you EXACTLY which routes those are and what fraction of that species'
## total routes they represent -- so you can flag/exclude just those rows
## from all_route_trends_*.csv while keeping the rest of that species'
## (converged) routes as normal.
##
## Reads output/files/nonconvergence_diagnosis_<firstYear>_<lastYear>.csv
## (written by helper/diagnose_nonconvergence.R -- run that first if this
## doesn't exist yet), filters to category %in% c("route_alpha",
## "route_beta") (i.e. skips gamma1 and hyperparameter_sd/other rows, which
## have no single route to point at), extracts the routeF index from each
## variable name (e.g. "alpha_raw[348]" -> 348), and joins against the
## matching route_info.rds for that EXACT species/tag -- routeF is only
## meaningful within one species/tag's own fit (it's assigned via
## as.integer(factor(new_data$route)) on that species' own reduced dataset),
## so the same routeF number means a different physical route for a
## different species. Never mix routeF across species/tags.
##
## Writes:
##   output/files/nonconvergence_flagged_routes_<firstYear>_<lastYear>.csv
##     -- one row per offending route-level parameter, with the real route
##     ID and coordinates attached.
##   output/files/nonconvergence_flagged_routes_summary_<firstYear>_<lastYear>.csv
##     -- one row per species/tag: how many distinct routes are flagged, out
##     of how many total, and what percentage that is (a handful of edge
##     routes reads very differently from a large chunk of a species'
##     coverage).
##
## Read-only: does not modify any fit output, all_route_trends, or anything
## else. Safe to re-run any time.
##
## Usage: run helper/diagnose_nonconvergence.R first (or make sure its
## output CSV already exists), then source this file or
## Rscript helper/translate_offending_routes.R.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/translate_offending_routes.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Extract the integer routeF index from a route-level Stan variable name,
#' e.g. "alpha_raw[348]" -> 348L, "beta_raw_space[1393]" -> 1393L. Returns
#' NA_integer_ for anything that doesn't match a "[<digits>]" pattern.
extract_routeF <- function(varname) {
  m <- regmatches(varname, regexpr("\\[([0-9]+)\\]", varname))
  suppressWarnings(as.integer(gsub("\\[|\\]", "", m)))
}

#' For one species/tag, read its route_info.rds (routeF -> real route ID +
#' latitude/longitude). Returns NULL (with a message) if the file is
#' missing, instead of erroring, so a batch loop can skip cleanly.
read_route_info <- function(species_f, model_tag, firstYear, lastYear,
                            route_info_dir = here::here("data", "route_info")) {
  route_info_file <- file.path(route_info_dir,
                               paste0(species_f, "_", model_tag, "_",
                                      firstYear, "_", lastYear, "_route_info.rds"))
  if (!file.exists(route_info_file)) {
    message("  [MISSING route_info] ", basename(route_info_file))
    return(NULL)
  }
  readRDS(route_info_file)
}

#' Main: read helper/diagnose_nonconvergence.R's output CSV, filter to
#' route-level offenders (route_alpha/route_beta), and join each
#' species/tag's routeF indices against its OWN route_info.rds to attach
#' the real route ID and coordinates.
#'
#' @param firstYear,lastYear must match nonconvergence_diagnosis_*.csv and
#'   the route_info.rds filenames
#' @param diagnosis_csv path to helper/diagnose_nonconvergence.R's output
#'   (default output/files/nonconvergence_diagnosis_<firstYear>_<lastYear>.csv)
#' @param route_info_dir directory containing the *_route_info.rds files
#' @param out_dir where to write the combined CSV (default output/files)
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the combined flagged-routes table (invisibly), or NULL if there's
#'   nothing to translate (no route_alpha/route_beta offenders at all)
translate_offending_routes <- function(firstYear = 2010, lastYear = 2025,
                                       diagnosis_csv = here::here("output", "files",
                                                                  paste0("nonconvergence_diagnosis_",
                                                                         firstYear, "_", lastYear, ".csv")),
                                       route_info_dir = here::here("data", "route_info"),
                                       out_dir = here::here("output", "files"),
                                       write_csv = TRUE) {

  if (!file.exists(diagnosis_csv)) {
    stop("Can't find ", diagnosis_csv, " -- run helper/diagnose_nonconvergence.R first.")
  }

  diagnosis <- read.csv(diagnosis_csv, stringsAsFactors = FALSE)

  route_level <- diagnosis %>%
    filter(category %in% c("route_alpha", "route_beta")) %>%
    mutate(routeF = extract_routeF(variable))

  if (nrow(route_level) == 0) {
    cat("No route_alpha/route_beta offending parameters found in", diagnosis_csv,
        "-- nothing to translate (gamma1/hyperparameter_sd/other rows, if any, don't map to a route).\n")
    return(invisible(NULL))
  }

  cat("=== Translating", nrow(route_level), "route-level offending parameter row(s) to real route IDs ===\n")

  # Join species/tag by species/tag, one at a time -- routeF is only
  # meaningful WITHIN a single species/tag's own fit, never across them.
  combos <- route_level %>% distinct(species, model)

  translated_list <- list()
  route_info_cache <- list()

  for (i in seq_len(nrow(combos))) {
    sp   <- combos$species[i]
    tag  <- combos$model[i]
    sp_f <- species_to_f(sp)

    cache_key <- paste(sp_f, tag, sep = "|")
    if (is.null(route_info_cache[[cache_key]])) {
      route_info_cache[[cache_key]] <- read_route_info(sp_f, tag, firstYear, lastYear, route_info_dir)
    }
    ri <- route_info_cache[[cache_key]]
    if (is.null(ri)) next

    this_combo <- route_level %>% filter(species == sp, model == tag)
    joined <- this_combo %>% left_join(ri, by = "routeF")

    n_total_routes <- nrow(ri)
    n_flagged <- length(unique(joined$routeF))
    cat("  ", sp, "(", tag, "): ", n_flagged, "of", n_total_routes, "routes flagged (",
        round(100 * n_flagged / n_total_routes, 1), "%)\n", sep = "")

    joined$n_total_routes_this_fit <- n_total_routes
    translated_list[[cache_key]] <- joined
  }

  if (length(translated_list) == 0) {
    message("No route_info.rds files could be read for any flagged species/tag -- ",
            "check that data/route_info/ still has the matching files.")
    return(invisible(NULL))
  }

  translated_all <- bind_rows(translated_list) %>%
    select(species, species_code, model, variable, category, rhat, ess_bulk,
           routeF, route, latitude, longitude, n_total_routes_this_fit, everything())

  route_summary <- translated_all %>%
    group_by(species, species_code, model) %>%
    summarise(n_flagged_routes = n_distinct(routeF),
             n_total_routes   = dplyr::first(n_total_routes_this_fit),
             pct_flagged      = round(100 * n_flagged_routes / n_total_routes, 1),
             .groups = "drop") %>%
    arrange(desc(pct_flagged))

  cat("\n=== Per-species/tag: distinct routes flagged (out of that fit's total routes) ===\n")
  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(route_summary))
  options(na.print = old_na_print)

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("nonconvergence_flagged_routes_", firstYear, "_", lastYear, ".csv"))
    write.csv(translated_all, out_csv, row.names = FALSE)
    cat("\nFlagged routes (with real route IDs + coordinates) written to:", out_csv, "\n")

    summary_csv <- file.path(out_dir, paste0("nonconvergence_flagged_routes_summary_", firstYear, "_", lastYear, ".csv"))
    write.csv(route_summary, summary_csv, row.names = FALSE)
    cat("Per-species/tag flagged-route summary written to:", summary_csv, "\n")
  }

  invisible(translated_all)
}

## ==========================================================================
## Auto-run with defaults, same soft-coded pattern as the rest of helper/ --
## set `translate_offending_routes_skip_autorun <- TRUE` before sourcing to
## load just the functions above without running anything.
## ==========================================================================
if (!exists("translate_offending_routes_skip_autorun") || !isTRUE(translate_offending_routes_skip_autorun)) {
  if (!exists("firstYear")) firstYear <- 2010
  if (!exists("lastYear"))  lastYear  <- 2025
  translate_offending_routes(firstYear = firstYear, lastYear = lastYear)
}
