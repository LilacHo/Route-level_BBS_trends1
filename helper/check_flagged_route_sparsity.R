## =============================================================================
## helper/check_flagged_route_sparsity.R
##
## Tests TWO hypotheses for why certain routes' alpha_raw[r]/beta_raw_space[r]
## fail whole-model convergence even for large, well-sampled species
## (Sharp-shinned Hawk, Hooded Merganser, Broad-winged Hawk, Belted
## Kingfisher) that already went through 1d_refit_nonconverged_species.R's
## standard refit and still failed:
##
##   1. THIN-DATA hypothesis (per-route count level): flagged routes simply
##      have less information (lower mean count) than this species' other
##      routes. TESTED AND REJECTED in chat -- flagged routes did not have
##      lower mean counts; if anything two species/tags showed the OPPOSITE
##      (flagged routes had significantly HIGHER mean counts).
##
##   2. COUNT-VOLATILITY hypothesis: it isn't the LEVEL of a route's counts
##      that matters, but how ERRATIC/spiky they are year to year (high
##      coefficient of variation) -- a volatile route is harder for the
##      model's own likelihood to summarize with a single smooth trend.
##
##   3. NEIGHBOR-MISMATCH hypothesis: a route whose own average count level
##      is very different from its immediate spatial neighbors' (on the
##      model's log scale) creates tension between what the route's own
##      data wants and what the CAR/ICAR spatial-smoothing prior wants to
##      pull it toward -- that tension is exactly the kind of thing that
##      produces slow mixing without divergences or treedepth warnings
##      (which is what we saw in the log for these species).
##
## Both (2) and (3) are computed per flagged species/tag from the SAME
## already-saved data/stan_data/<species>_<tag>_..._stan_data.RData used for
## the (rejected) thin-data check -- (2) needs only new_data$count per
## route (already loaded there); (3) additionally needs the CAR adjacency
## edge list (stan_data$node1/node2) to look up each route's actual spatial
## neighbors.
##
## For every species/tag combo in
## output/files/nonconvergence_flagged_routes_<firstYear>_<lastYear>.csv
## (written by helper/translate_offending_routes.R), this script:
##   - computes per-route mean count, count coefficient of variation (CV =
##     sd/mean across years), and (via the adjacency edges) each route's
##     mismatch from its neighbors' average log-count level
##   - compares FLAGGED (non-converged) routes against the SAME fit's OTHER
##     routes on both metrics via two-sample Wilcoxon rank-sum tests
##
## Writes:
##   output/files/flagged_route_sparsity_check_<firstYear>_<lastYear>.csv
##     -- one row per species/tag, with both comparisons side by side
##   output/files/flagged_route_sparsity_detail_<firstYear>_<lastYear>.csv
##     -- one row per route, with its own stats + neighbor-mismatch value
##
## Read-only: does not modify any fit output. Safe to re-run any time.
##
## Usage: run helper/diagnose_nonconvergence.R and
## helper/translate_offending_routes.R first (needs their output CSV), then
## source this file or Rscript helper/check_flagged_route_sparsity.R.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/check_flagged_route_sparsity.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Build a routeF -> list-of-neighbor-routeF lookup from the CAR/ICAR
#' adjacency edge list saved in stan_data (node1, node2 -- one row per
#' undirected edge, 1-indexed to match routeF). Symmetric: if node1==r
#' contributes node2 as a neighbor of r, and vice versa.
#'
#' @param node1,node2 integer vectors of equal length, from stan_data
#' @return a named list, names are routeF (as character), values are
#'   integer vectors of that route's neighboring routeF indices
build_neighbor_lookup <- function(node1, node2) {
  neighbors <- list()
  add_edge <- function(a, b) {
    key <- as.character(a)
    neighbors[[key]] <<- c(neighbors[[key]], b)
  }
  for (i in seq_along(node1)) {
    add_edge(node1[i], node2[i])
    add_edge(node2[i], node1[i])
  }
  neighbors
}

#' For one species/tag, load its saved stan_data.RData, compute per-route
#' mean count, count coefficient of variation, and neighbor-mismatch (own
#' vs. neighbors' average log-count level), then compare the FLAGGED routes
#' against the rest on both metrics via two-sample Wilcoxon rank-sum tests.
#'
#' @param species_f file-safe species name
#' @param model_tag "base" or "anthro"
#' @param flagged_routeFs integer vector of routeF indices flagged as
#'   offending for this exact species/tag
#' @param firstYear,lastYear must match the stan_data.RData filename
#' @param stan_data_dir directory containing the *_stan_data.RData files
#' @return a list with $summary (one-row data.frame) and $per_route (full
#'   per-route detail data.frame), or NULL if the file is missing
check_route_sparsity_for_combo <- function(species_f, model_tag, flagged_routeFs,
                                           firstYear, lastYear,
                                           stan_data_dir = here::here("data", "stan_data")) {

  sp_data_file <- file.path(stan_data_dir,
                            paste0(species_f, "_", model_tag, "_",
                                   firstYear, "_", lastYear, "_stan_data.RData"))
  if (!file.exists(sp_data_file)) {
    message("  [MISSING stan_data] ", basename(sp_data_file))
    return(NULL)
  }

  # Loads new_data (route, routeF, count, year, ...) and stan_data (has the
  # CAR adjacency edges node1/node2) into a throwaway environment.
  e <- new.env()
  load(sp_data_file, envir = e)
  new_data  <- e$new_data
  stan_data <- e$stan_data

  # --- Per-route level + volatility stats -----------------------------------
  per_route <- new_data %>%
    group_by(routeF) %>%
    summarise(n_years        = n(),
             mean_count     = mean(count, na.rm = TRUE),
             sd_count       = sd(count, na.rm = TRUE),
             pct_zero_years = mean(count == 0, na.rm = TRUE) * 100,
             .groups = "drop") %>%
    mutate(
      # CV undefined (0/0) when mean_count is 0 -- leave NA rather than Inf.
      count_cv = ifelse(mean_count > 0, sd_count / mean_count, NA_real_),
      log_mean_count = log1p(mean_count),
      flagged = routeF %in% flagged_routeFs
    )

  # --- Neighbor-mismatch: each route's own log(mean count) vs. the average
  # of its immediate CAR-graph neighbors' log(mean count) -----------------
  neighbor_lookup <- build_neighbor_lookup(stan_data$node1, stan_data$node2)
  log_mean_by_route <- setNames(per_route$log_mean_count, per_route$routeF)

  neighbor_mismatch <- vapply(per_route$routeF, function(r) {
    nbrs <- neighbor_lookup[[as.character(r)]]
    if (is.null(nbrs) || length(nbrs) == 0) return(NA_real_)
    nbr_vals <- log_mean_by_route[as.character(nbrs)]
    nbr_vals <- nbr_vals[!is.na(nbr_vals)]
    if (length(nbr_vals) == 0) return(NA_real_)
    abs(log_mean_by_route[as.character(r)] - mean(nbr_vals))
  }, numeric(1))

  per_route$n_neighbors <- vapply(per_route$routeF, function(r) {
    length(neighbor_lookup[[as.character(r)]])
  }, integer(1))
  per_route$neighbor_mismatch <- neighbor_mismatch

  n_flagged <- sum(per_route$flagged)
  n_other   <- sum(!per_route$flagged)

  if (n_flagged == 0 || n_other == 0) {
    message("  [SKIP] ", species_f, " (", model_tag, ") -- need both flagged AND non-flagged ",
            "routes to compare (flagged=", n_flagged, ", other=", n_other, ")")
    return(NULL)
  }

  wilcox_safe <- function(x, y) {
    tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  }

  cv_flagged <- per_route$count_cv[per_route$flagged]
  cv_other   <- per_route$count_cv[!per_route$flagged]
  mismatch_flagged <- per_route$neighbor_mismatch[per_route$flagged]
  mismatch_other   <- per_route$neighbor_mismatch[!per_route$flagged]

  summary_row <- data.frame(
    species_f          = species_f,
    model               = model_tag,
    n_flagged_routes    = n_flagged,
    n_other_routes      = n_other,
    # count-volatility comparison
    cv_flagged_mean     = round(mean(cv_flagged, na.rm = TRUE), 3),
    cv_other_mean       = round(mean(cv_other, na.rm = TRUE), 3),
    cv_wilcox_p         = round(wilcox_safe(cv_flagged, cv_other), 4),
    flagged_more_volatile = mean(cv_flagged, na.rm = TRUE) > mean(cv_other, na.rm = TRUE),
    # neighbor-mismatch comparison
    mismatch_flagged_mean = round(mean(mismatch_flagged, na.rm = TRUE), 3),
    mismatch_other_mean   = round(mean(mismatch_other, na.rm = TRUE), 3),
    mismatch_wilcox_p     = round(wilcox_safe(mismatch_flagged, mismatch_other), 4),
    flagged_more_mismatched = mean(mismatch_flagged, na.rm = TRUE) > mean(mismatch_other, na.rm = TRUE),
    n_flagged_with_no_neighbors = sum(per_route$flagged & per_route$n_neighbors == 0)
  )

  list(summary = summary_row, per_route = per_route %>% mutate(species_f = species_f, model = model_tag))
}

#' Main: read helper/translate_offending_routes.R's flagged-routes CSV,
#' group by species/tag, and run check_route_sparsity_for_combo() for each.
#'
#' @param firstYear,lastYear must match the flagged-routes CSV and
#'   stan_data.RData filenames
#' @param flagged_routes_csv path to translate_offending_routes.R's output
#' @param stan_data_dir directory containing the *_stan_data.RData files
#' @param out_dir where to write the combined CSVs (default output/files)
#' @param write_csv if FALSE, skip writing CSVs (just return the tables)
#' @return list(summary = combined summary table, per_route = combined
#'   per-route detail table), invisibly, or NULL if nothing to check
check_flagged_route_sparsity_table <- function(firstYear = 2010, lastYear = 2025,
                                               flagged_routes_csv = here::here("output", "files",
                                                                               paste0("nonconvergence_flagged_routes_",
                                                                                      firstYear, "_", lastYear, ".csv")),
                                               stan_data_dir = here::here("data", "stan_data"),
                                               out_dir = here::here("output", "files"),
                                               write_csv = TRUE) {

  if (!file.exists(flagged_routes_csv)) {
    stop("Can't find ", flagged_routes_csv, " -- run helper/translate_offending_routes.R first.")
  }

  flagged <- read.csv(flagged_routes_csv, stringsAsFactors = FALSE)

  combos <- flagged %>% distinct(species, species_code, model)
  cat("=== Checking count-volatility and neighbor-mismatch: flagged (non-converged) routes",
      "vs. others, for", nrow(combos), "species/tag combo(s) ===\n")

  summary_list <- list()
  per_route_list <- list()

  for (i in seq_len(nrow(combos))) {
    sp      <- combos$species[i]
    sp_code <- combos$species_code[i]
    tag     <- combos$model[i]
    sp_f    <- species_to_f(sp)

    flagged_routeFs <- flagged %>%
      filter(species == sp, model == tag) %>%
      pull(routeF) %>%
      unique()

    r <- check_route_sparsity_for_combo(sp_f, tag, flagged_routeFs, firstYear, lastYear, stan_data_dir)
    if (is.null(r)) next

    r$summary$species      <- sp
    r$summary$species_code <- sp_code
    summary_list[[paste(sp, tag, sep = " | ")]] <- r$summary
    per_route_list[[paste(sp, tag, sep = " | ")]] <- r$per_route %>% mutate(species = sp, species_code = sp_code)
  }

  if (length(summary_list) == 0) {
    message("Nothing could be checked -- no matching stan_data.RData files found, or no ",
            "species/tag had both flagged and non-flagged routes.")
    return(invisible(NULL))
  }

  summary_all <- bind_rows(summary_list) %>%
    select(species, species_code, model, n_flagged_routes, n_other_routes,
           cv_flagged_mean, cv_other_mean, cv_wilcox_p, flagged_more_volatile,
           mismatch_flagged_mean, mismatch_other_mean, mismatch_wilcox_p, flagged_more_mismatched,
           n_flagged_with_no_neighbors, everything(), -species_f)

  per_route_all <- bind_rows(per_route_list) %>%
    select(species, species_code, model, routeF, flagged, n_years, mean_count, count_cv,
           neighbor_mismatch, n_neighbors, everything(), -species_f)

  cat("\n=== Count volatility (CV = sd/mean across years): flagged vs. other routes ===\n")
  cat("(flagged_more_volatile = TRUE means non-converged routes have MORE erratic year-to-year\n")
  cat(" counts than this species' other routes -- supports the volatility hypothesis.)\n")
  old_na_print <- getOption("na.print")
  options(na.print = "NA")
  print(as.data.frame(summary_all %>%
                        select(species, model, n_flagged_routes, n_other_routes,
                               cv_flagged_mean, cv_other_mean, cv_wilcox_p, flagged_more_volatile)))

  cat("\n=== Neighbor mismatch (own vs. spatial-neighbors' log-count level): flagged vs. other ===\n")
  cat("(flagged_more_mismatched = TRUE means non-converged routes look MORE different from their\n")
  cat(" immediate spatial neighbors than this species' other routes do -- supports the\n")
  cat(" neighbor-mismatch/CAR-tension hypothesis.)\n")
  print(as.data.frame(summary_all %>%
                        select(species, model, mismatch_flagged_mean, mismatch_other_mean,
                               mismatch_wilcox_p, flagged_more_mismatched, n_flagged_with_no_neighbors)))
  options(na.print = old_na_print)

  n_vol_supporting <- sum(summary_all$flagged_more_volatile, na.rm = TRUE)
  n_vol_sig        <- sum(summary_all$cv_wilcox_p < 0.05, na.rm = TRUE)
  n_mis_supporting <- sum(summary_all$flagged_more_mismatched, na.rm = TRUE)
  n_mis_sig        <- sum(summary_all$mismatch_wilcox_p < 0.05, na.rm = TRUE)
  cat("\nVolatility hypothesis: ", n_vol_supporting, " of ", nrow(summary_all),
      " combo(s) support it (", n_vol_sig, " significant, p<0.05).\n", sep = "")
  cat("Neighbor-mismatch hypothesis: ", n_mis_supporting, " of ", nrow(summary_all),
      " combo(s) support it (", n_mis_sig, " significant, p<0.05).\n", sep = "")

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    summary_csv <- file.path(out_dir, paste0("flagged_route_sparsity_check_", firstYear, "_", lastYear, ".csv"))
    write.csv(summary_all, summary_csv, row.names = FALSE)
    cat("\nSummary written to:", summary_csv, "\n")

    detail_csv <- file.path(out_dir, paste0("flagged_route_sparsity_detail_", firstYear, "_", lastYear, ".csv"))
    write.csv(per_route_all, detail_csv, row.names = FALSE)
    cat("Per-route detail written to:", detail_csv, "\n")
  }

  invisible(list(summary = summary_all, per_route = per_route_all))
}

## ==========================================================================
## Auto-run with defaults, same soft-coded pattern as the rest of helper/ --
## set `check_flagged_route_sparsity_skip_autorun <- TRUE` before sourcing
## to load just the functions above without running anything.
## ==========================================================================
if (!exists("check_flagged_route_sparsity_skip_autorun") || !isTRUE(check_flagged_route_sparsity_skip_autorun)) {
  if (!exists("firstYear")) firstYear <- 2010
  if (!exists("lastYear"))  lastYear  <- 2025
  check_flagged_route_sparsity_table(firstYear = firstYear, lastYear = lastYear)
}
