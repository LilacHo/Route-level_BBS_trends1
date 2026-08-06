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
## Usage:
##   - source("helper/gamma_lookup.R") to get the gamma_lookup() function for
##     use on a single species/model, e.g.:
##       gamma_lookup("Bairds_Sparrow", "grassland_habitat")
##   - Rscript helper/gamma_lookup.R (or source() from the console) also runs
##     the batch section below, which reads every species x non-base model
##     combination for `land_cover` and writes one combined table to
##     output/files/gamma_lookup_<land_cover>_<firstYear>_<lastYear>.csv,
##     plus prints a quick convergence/credibility summary to the console.
## =============================================================================

library(dplyr)
library(here)

here::i_am("helper/gamma_lookup.R")

# Settings (match 1c_species_iCAR_covariates.R / 2c_generate_route_trend_csvs_covariates.R)
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2025
model_tags <- c("grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro")
# ("base" is intentionally excluded -- it has no gamma1.)

rds_dir <- here::here("output", "rds")

#' Read one species/model's *_summ_fit.rds and pull out its gamma1 row
#' verbatim (every column posterior::summarise_draws() produced), tagged
#' with species/model/file for identification. Returns NULL (with a message)
#' if the file is missing or has no gamma1 row, instead of erroring, so a
#' batch loop can skip cleanly.
#'
#' @param species_f file-safe species name, e.g. "Bairds_Sparrow" (matches
#'   the species_to_f() naming used by 1c/2c)
#' @param tag covariate model tag, e.g. "grassland_habitat"
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

## ==========================================================================
## Batch mode: every species in `land_cover` x every tag in `model_tags`.
## Runs automatically whenever this file is executed top-to-bottom
## (Rscript helper/gamma_lookup.R, or source("helper/gamma_lookup.R")).
## ==========================================================================

species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== gamma1 lookup ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("Species (n =", nrow(target_spp), ") x models:", paste(model_tags, collapse = ", "), "\n\n")

rows <- list()
for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]

  for (tag in model_tags) {
    r <- gamma_lookup(sp_f, tag, firstYear, lastYear, rds_dir)
    if (!is.null(r)) {
      r$species      <- sp
      r$species_code <- sp_code
      rows[[paste(sp, tag, sep = " | ")]] <- r
    }
  }
}

if (length(rows) == 0) {
  stop("No gamma1 rows found -- check that output/rds/*_summ_fit.rds exist ",
       "(run 1c_species_iCAR_covariates.R first) and that land_cover/model_tags ",
       "above match your setup.")
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

out_dir <- here::here("output", "files")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_csv <- file.path(out_dir, paste0("gamma_lookup_", land_cover, "_",
                                     firstYear, "_", lastYear, ".csv"))
write.csv(gamma_table, out_csv, row.names = FALSE)

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

cat("\nFull table written to:", out_csv, "\n")
