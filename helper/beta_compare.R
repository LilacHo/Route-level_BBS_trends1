## =============================================================================
## helper/beta_compare.R
##
## Compare route-level trend (beta[r]) between two already-fitted models for
## the same species -- typically "base" (no covariate) vs. a single-covariate
## model like "anthro" -- to see how much of the raw population trend gets
## reassigned to the covariate (gamma1) once it's added, and how much is left
## over in beta.
##
## Reads directly from the EXISTING *_summ_fit.rds files written by a 1c
## script (no re-fitting, no stanfit.rds/full-draws loading -- summ_fit.rds is
## the small posterior::summarise_draws() tibble, same file gamma_lookup.R
## reads from). This only needs the summ_fit.rds files to be present locally
## -- it does not need cmdstanr or a Stan installation.
##
## Nothing here is hardcoded to a species or a pair of tags -- every setting
## is an argument, same soft-coding pattern as gamma_lookup.R and
## covariate_variance_decomp.R.
##
## Interpretation reminder (see chat discussion): beta with a covariate in
## the model is NOT "the climate-only trend" -- it's "the trend left over
## after the model reassigned whatever it could explain to the covariate."
## If beta shrinks noticeably from tag_a to tag_b (e.g. base -> anthro),
## that's evidence the covariate is explaining part of the raw decline. If it
## barely moves, the covariate isn't doing much of that species' story, even
## if its own gamma1 came back credible.
##
## Usage:
##   source(here::here("helper", "beta_compare.R"))
##
##   # single species, one call:
##   res <- beta_compare("Bells_Vireo", tag_a = "base", tag_b = "anthro")
##   res$summary      # one-row summary (mean/median beta both models, etc.)
##   res$by_route      # full route-by-route comparison
##
##   # several species at once, printed + written to CSV:
##   target_spp <- data.frame(Common.Name = c("Bell's Vireo", "Sage Thrasher",
##                                            "Ladder-backed Woodpecker", "Lark Sparrow"))
##   beta_compare_table(target_spp, tag_a = "base", tag_b = "anthro",
##                      label = "aridlands_credible4")
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/beta_compare.R")

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

#' Pull every beta[route_idx] row out of a summ_fit.rds's posterior summary.
#' @param summ the data.frame read from a *_summ_fit.rds file
#' @param prefix the Stan parameter name for the route-level trend (default
#'   "beta" -- matches slope_iCAR_route_NB_New.stan / _covariate.stan)
extract_beta <- function(summ, prefix = "beta") {
  beta_rows <- summ %>% filter(grepl(paste0("^", prefix, "\\["), variable))
  beta_rows$route_idx <- as.integer(gsub(paste0(prefix, "\\[(\\d+)\\]"), "\\1", beta_rows$variable))
  beta_rows
}

#' Compare route-level beta between two fitted models for ONE species.
#'
#' @param species_f file-safe species name (matches species_to_f() naming
#'   used by 1c/2c/gamma_lookup, e.g. "Bells_Vireo")
#' @param tag_a,tag_b the two model_tags to compare (e.g. "base", "anthro")
#' @param firstYear,lastYear,rds_dir must match how the *_summ_fit.rds files
#'   were named/saved
#' @param beta_prefix Stan parameter name for route-level trend (default "beta")
#' @return a list(species_f, tag_a, tag_b, n_routes, by_route, summary), or
#'   NULL (with a message) if either file/parameter set is missing -- so a
#'   batch loop (beta_compare_table()) can skip cleanly
beta_compare <- function(species_f, tag_a, tag_b, firstYear = 2010, lastYear = 2025,
                         rds_dir = here::here("output", "rds"),
                         beta_prefix = "beta") {

  file_a <- file.path(rds_dir, paste0(species_f, "_iCAR_", tag_a, "_", firstYear, "_", lastYear, "_summ_fit.rds"))
  file_b <- file.path(rds_dir, paste0(species_f, "_iCAR_", tag_b, "_", firstYear, "_", lastYear, "_summ_fit.rds"))

  if (!file.exists(file_a)) { message("  [MISSING] ", basename(file_a)); return(NULL) }
  if (!file.exists(file_b)) { message("  [MISSING] ", basename(file_b)); return(NULL) }

  summ_a <- readRDS(file_a)
  summ_b <- readRDS(file_b)

  beta_a <- extract_beta(summ_a, beta_prefix) %>%
    select(route_idx, mean, median, sd, rhat, ess_bulk) %>%
    rename(mean_a = mean, median_a = median, sd_a = sd, rhat_a = rhat, ess_bulk_a = ess_bulk)
  beta_b <- extract_beta(summ_b, beta_prefix) %>%
    select(route_idx, mean, median, sd, rhat, ess_bulk) %>%
    rename(mean_b = mean, median_b = median, sd_b = sd, rhat_b = rhat, ess_bulk_b = ess_bulk)

  if (nrow(beta_a) == 0) { message("  [NO ", beta_prefix, " rows] ", basename(file_a)); return(NULL) }
  if (nrow(beta_b) == 0) { message("  [NO ", beta_prefix, " rows] ", basename(file_b)); return(NULL) }

  # Route indexing (routeF) is only comparable between two fits of the SAME
  # species if they were built from the same prepare_species_data() call
  # (same route set/order) -- true for base vs. a covariate model in the
  # current 1c scripts, since both draw from the same `prepped` object. If
  # nrow differs between the two fits, routes were NOT identical (e.g. one
  # tag dropped more rows for missing covariate data) -- flagged below rather
  # than silently mismatched.
  if (nrow(beta_a) != nrow(beta_b)) {
    message("  [WARNING] ", species_f, ": ", tag_a, " has ", nrow(beta_a),
            " routes, ", tag_b, " has ", nrow(beta_b),
            " -- route sets may not align 1:1 (joins below use route_idx, ",
            "mismatches will just drop out of the inner_join).")
  }

  joined <- beta_a %>%
    inner_join(beta_b, by = "route_idx") %>%
    mutate(diff        = mean_b - mean_a,
           pct_trend_a = (exp(mean_a) - 1) * 100,
           pct_trend_b = (exp(mean_b) - 1) * 100)

  summ_row <- data.frame(
    species_f        = species_f,
    tag_a             = tag_a,
    tag_b             = tag_b,
    n_routes          = nrow(joined),
    mean_beta_a       = mean(joined$mean_a),
    mean_beta_b       = mean(joined$mean_b),
    median_beta_a     = median(joined$mean_a),
    median_beta_b     = median(joined$mean_b),
    mean_pct_trend_a  = mean(joined$pct_trend_a),   # avg annual % change, model a
    mean_pct_trend_b  = mean(joined$pct_trend_b),   # avg annual % change, model b
    mean_abs_diff     = mean(abs(joined$diff)),
    route_beta_corr   = cor(joined$mean_a, joined$mean_b),
    max_rhat_a        = max(joined$rhat_a, na.rm = TRUE),
    max_rhat_b        = max(joined$rhat_b, na.rm = TRUE)
  )

  list(species_f = species_f, tag_a = tag_a, tag_b = tag_b,
       n_routes = nrow(joined), by_route = joined, summary = summ_row)
}

#' Run beta_compare() across several species and print/return a combined
#' summary table, comparing tag_a vs. tag_b's route-level trend estimates.
#'
#' @param target_spp data.frame with at least a Common.Name column
#' @param tag_a,tag_b,firstYear,lastYear,rds_dir passed through to
#'   beta_compare()
#' @param out_dir where to write the combined CSV (default output/files)
#' @param label used only to name the output CSV
#' @param write_csv if FALSE, skip writing the CSV (just return the table)
#' @return the combined summary table (invisibly)
beta_compare_table <- function(target_spp, tag_a = "base", tag_b,
                               firstYear = 2010, lastYear = 2025,
                               rds_dir = here::here("output", "rds"),
                               out_dir = here::here("output", "files"),
                               label = "comparison", write_csv = TRUE) {

  cat("=== beta comparison:", tag_a, "vs", tag_b, "===\n")
  cat("Species (n =", nrow(target_spp), ") | Period:", firstYear, "-", lastYear, "\n\n")

  rows <- list()
  for (i in seq_len(nrow(target_spp))) {
    sp   <- target_spp$Common.Name[i]
    sp_f <- species_to_f(sp)
    res  <- beta_compare(sp_f, tag_a, tag_b, firstYear, lastYear, rds_dir)
    if (!is.null(res)) {
      res$summary$species <- sp
      rows[[sp]] <- res$summary
    }
  }

  if (length(rows) == 0) {
    message("No comparisons could be built -- check that *_summ_fit.rds files ",
            "exist for both '", tag_a, "' and '", tag_b, "' for these species.")
    return(invisible(NULL))
  }

  out <- bind_rows(rows) %>% select(species, everything(), -species_f)

  cat("\n=== Summary table ===\n")
  print(out, row.names = FALSE)

  cat("\nmean_pct_trend_a/b = average annual %% population change implied by\n",
      "each route's beta, averaged across routes -- i.e. (exp(beta)-1)*100,\n",
      "then averaged. This is a quick, interpretable summary, NOT a full\n",
      "composite/weighted trend index (see 4c if this project already computes one).\n",
      "route_beta_corr close to 1 means the two models' route-level trends track\n",
      "each other closely (the covariate barely changed the trend story); lower\n",
      "values mean the covariate meaningfully reshuffled which routes look like\n",
      "they're declining vs. increasing once it's accounted for.\n", sep = "")

  if (write_csv) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    out_csv <- file.path(out_dir, paste0("beta_compare_", tag_a, "_vs_", tag_b, "_",
                                         label, "_", firstYear, "_", lastYear, ".csv"))
    write.csv(out, out_csv, row.names = FALSE)
    cat("\nTable written to:", out_csv, "\n")
  }

  invisible(out)
}
