## =============================================================================
## Combine per-route trend output across ALL fitted models into one giant CSV,
## for cross-model comparison.
## Time period: 2010-2025
##
## Combines FOUR "model" sources per species, all fit by
## 1c_species_iCAR_covariates.R on the SAME covariate-reduced dataset (so
## they're directly comparable to each other): "base" (no covariates),
## "grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro".
## All three non-base models are single-covariate fits (there is no joint
## two-covariate model in this design — "grassland_anthro" here is a single
## pre-composited covariate, Developed*+CultivatedCrops, not grassland and
## anthro fit together).
##
## NOTE: "base" here is 1c's own no-covariate refit on the reduced dataset —
## NOT the full-dataset base fit from 1_species_iCAR_2010_2025.R /
## 2_generate_route_trend_csvs.R. That full-dataset fit uses more routes and
## is intentionally not comparable to the covariate models, so it is not
## included in this combined file.
##
## For all four models, reads the summary fit + lightweight route lookup
## saved by 1c_species_iCAR_covariates.R:
##   output/rds/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/route_info/<species>_<tag>_<firstYear>_<lastYear>_route_info.rds
##     (route, routeF, latitude, longitude only — NOT the full stan_data.RData,
##     which also exists but isn't needed here)
## and derives the same per-route quantities as 2_generate_route_trend_csvs.R
## (trend, trend_lci, trend_uci, rel_abundance from beta/alpha).
##
## Writes TWO combined CSVs across every species and every model:
##
##   1) Per-route file — one row per route/species/model:
##      output/species_routes_covariates/all_route_trends_<land_cover>_<firstYear>_<lastYear>.csv
##      columns: species, species_code, model, route, latitude, longitude,
##               alpha, beta, trend, trend_lci, trend_uci, rel_abundance
##
##   2) Model-level file — one row per species/model (NOT per route), for
##      quantities that are constant across a species' routes:
##      output/species_routes_covariates/all_model_level_summary_<land_cover>_<firstYear>_<lastYear>.csv
##      columns: species, species_code, model, n_routes, gamma1, gamma1_lci,
##               gamma1_uci, gamma1_excludes_zero
##
## gamma1 is the covariate effect on log(lambda) — a single value per
## species+model, not per-route, which is why it lives in the separate
## model-level file rather than being repeated down every route row. It is
## NA (all four gamma1* columns) for the "base" model (no covariate).
##
## gamma1_lci/gamma1_uci are the 5th/95th posterior percentiles (q5/q95) of
## gamma1 — i.e. a 90% credible interval, on the same log(lambda) scale as
## gamma1 itself (no exp() transform, unlike trend_lci/trend_uci). This
## mirrors how beta's 90% CI is pulled from summ's q5/q95 columns below.
## gamma1_excludes_zero is TRUE when that 90% CI does not span zero (the
## "credible effect" criterion described in the manuscript Methods).
##
## ALSO writes ONE per-route CSV per species-per-model (same route_trends
## data as above, just split apart instead of combined), so it can be fed to
## 3c_add_SDM_covariates.R the same way 2_generate_route_trend_csvs.R's
## per-species CSVs feed 3_add_SDM.R — that script requires exactly one
## species_code per input file, which the single combined CSV above doesn't
## satisfy:
##   output/species_routes_covariates/per_species/<species>_<model_tag>_route_trends.csv
##   columns: species, species_code, model, route, latitude, longitude,
##            alpha, beta, trend, trend_lci, trend_uci, rel_abundance
## (same columns as the combined per-route file, including "model", which
## 3c_add_SDM_covariates.R/4c_statistical_analysis_and_visualization_covariates.R
## just carry through unchanged since they have no model-specific logic.)
##
## This is a post-processing / combination step only; it does not fit or
## re-derive anything beyond what 2_generate_route_trend_csvs.R and
## 1c_species_iCAR_covariates.R already produced.
## =============================================================================

library(tidyverse)
library(here)

here::i_am("2c_generate_route_trend_csvs_covariates.R")

# Settings (match 1c_species_iCAR_covariates.R) ----------------------------
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2025
model_tags <- c("base", "grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro")

output_dir       <- here::here("output")
rds_dir          <- here::here("output", "rds")         # matches 1c_species_iCAR_covariates.R
route_info_dir   <- here::here("data", "route_info")    # matches 1c_species_iCAR_covariates.R
combined_out_dir <- here::here("output", "species_routes_covariates")
if (!dir.exists(combined_out_dir)) dir.create(combined_out_dir, recursive = TRUE)
per_species_out_dir <- here::here("output", "species_routes_covariates", "per_species")
if (!dir.exists(per_species_out_dir)) dir.create(per_species_out_dir, recursive = TRUE)

# Target group species list -------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Combine route-trend CSVs across all models ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("Species in group (n =", nrow(target_spp), ")\n")
cat("Models:", paste(model_tags, collapse = ", "), "\n")

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Column order every row is coerced to before combining --------------------
route_cols <- c("species", "species_code", "model", "route",
                "latitude", "longitude", "alpha", "beta",
                "trend", "trend_lci", "trend_uci", "rel_abundance")
model_cols <- c("species", "species_code", "model", "n_routes",
                "gamma1", "gamma1_lci", "gamma1_uci", "gamma1_excludes_zero")

all_route_rows <- list()
all_model_rows <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]

  cat("\n[", i, "/", nrow(target_spp), "]", sp, "\n")

  # All four models are fit by 1c_species_iCAR_covariates.R on the identical
  # reduced dataset, so all four are read the same way here.
  for (tag in model_tags) {

    out_base        <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file       <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    route_info_file <- file.path(route_info_dir,
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_route_info.rds"))

    if (!file.exists(summ_file) || !file.exists(route_info_file)) {
      cat("  [", tag, "] No fitted output found — run 1c_species_iCAR_covariates.R first. Skipping.\n")
      next
    }

    summ       <- readRDS(summ_file)
    route_info <- readRDS(route_info_file)   # route, routeF, latitude, longitude

    # Extract beta (slope) per route -> raw value + annual % trend + 90% CI --
    beta_summ <- summ %>%
      filter(str_detect(variable, "^beta\\[")) %>%
      transmute(routeF    = as.integer(str_extract(variable, "\\d+")),
                beta      = mean,
                trend     = 100 * (exp(mean) - 1),
                trend_lci = 100 * (exp(q5)   - 1),
                trend_uci = 100 * (exp(q95)  - 1))

    # Extract alpha (intercept) per route -> raw value + relative abundance --
    alpha_summ <- summ %>%
      filter(str_detect(variable, "^alpha\\[")) %>%
      transmute(routeF        = as.integer(str_extract(variable, "\\d+")),
                alpha         = mean,
                rel_abundance = exp(mean))

    # Species+model-level covariate effect, not per-route. "base" has no
    # gamma1; each of the three covariate models has exactly one. Pull the
    # 90% credible interval (q5/q95) alongside the mean, same way beta_summ
    # does above, so downstream analyses can assess whether gamma1 is
    # credibly different from zero instead of only seeing the point estimate.
    if ("gamma1" %in% summ$variable) {
      gamma1_row        <- summ[summ$variable == "gamma1", ]
      gamma1_val        <- gamma1_row$mean[1]
      gamma1_lci_val    <- gamma1_row$q5[1]
      gamma1_uci_val    <- gamma1_row$q95[1]
      gamma1_excludes_zero_val <- !is.na(gamma1_lci_val) && !is.na(gamma1_uci_val) &&
        (gamma1_lci_val > 0 || gamma1_uci_val < 0)
    } else {
      gamma1_val               <- NA_real_
      gamma1_lci_val            <- NA_real_
      gamma1_uci_val            <- NA_real_
      gamma1_excludes_zero_val <- NA
    }

    # route_info (route, routeF, latitude, longitude) already loaded above
    # from route_info.rds (already WGS84, no reprojection needed since 1c
    # builds it directly from raw_data$latitude/longitude).

    route_trends <- beta_summ %>%
      left_join(alpha_summ, by = "routeF") %>%
      left_join(route_info, by = "routeF") %>%
      mutate(species      = sp,
             species_code = sp_code,
             model        = tag) %>%
      arrange(routeF) %>%
      select(all_of(route_cols))

    all_route_rows[[paste(sp, tag, sep = " | ")]] <- route_trends

    all_model_rows[[paste(sp, tag, sep = " | ")]] <- data.frame(
      species              = sp,
      species_code         = sp_code,
      model                = tag,
      n_routes             = nrow(route_trends),
      gamma1               = gamma1_val,
      gamma1_lci           = gamma1_lci_val,
      gamma1_uci           = gamma1_uci_val,
      gamma1_excludes_zero = gamma1_excludes_zero_val
    )

    # Per-species-per-model CSV — same rows as route_trends above, just
    # written individually (one species_code per file) instead of only into
    # the combined all_route_rows list, so 3c_add_SDM_covariates.R can read
    # them the same way 3_add_SDM.R reads 2_generate_route_trend_csvs.R's
    # per-species files.
    per_species_csv <- file.path(per_species_out_dir,
                                 paste0(sp_f, "_", tag, "_route_trends.csv"))
    write.csv(route_trends, per_species_csv, row.names = FALSE)

    cat("  [", tag, "] ", nrow(route_trends), "routes -> ", basename(per_species_csv), "\n")
  }
}

if (length(all_route_rows) == 0) {
  stop("Nothing to combine — run 1c_species_iCAR_covariates.R first.")
}

all_route_trends <- bind_rows(all_route_rows)
all_model_summary <- bind_rows(all_model_rows)

route_csv <- file.path(combined_out_dir,
                       paste0("all_route_trends_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))
model_csv <- file.path(combined_out_dir,
                       paste0("all_model_level_summary_", land_cover, "_",
                              firstYear, "_", lastYear, ".csv"))

write.csv(all_route_trends, route_csv, row.names = FALSE)
write.csv(all_model_summary, model_csv, row.names = FALSE)

cat("\n=== Summary ===\n")
cat("Per-route rows:", nrow(all_route_trends), "\n")
cat("Species x model combinations:", length(all_route_rows), "\n")
print(table(all_route_trends$model))
cat("Per-route CSV written to:", route_csv, "\n")
cat("Model-level CSV written to:", model_csv, "\n")
cat("Per-species-per-model CSVs written to:", per_species_out_dir,
    "(", length(all_route_rows), "files )\n")

# Quick check on the new gamma1 CI columns: how many species have a 90%
# credible interval that excludes zero, per covariate model.
cat("\nSpecies with gamma1's 90% CI excluding zero (credible effect), by model:\n")
print(all_model_summary %>%
        filter(model != "base") %>%
        group_by(model) %>%
        summarise(n_species          = n(),
                  n_excludes_zero    = sum(gamma1_excludes_zero, na.rm = TRUE),
                  .groups = "drop"))

cat("Next: run 3c_add_SDM_covariates.R to add climate-suitability columns.\n")
