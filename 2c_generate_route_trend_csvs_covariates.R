## =============================================================================
## 2c_generate_route_trend_csvs_covariates.R
##
## Post-processing only, no fitting. Converts the raw alpha/beta posteriors
## from 1c into interpretable per-route trend estimates, for all four model
## tags (base, grassland_habitat, grassland_anthro,
## grassland_habitat_to_anthro), and combines them into single CSVs across
## every species.
##
## Input: output/rds/*_summ_fit.rds and data/route_info/*_route_info.rds
## (both from 1c), data/spp_names_codes_group_aou.csv.
##
## Output:
##   - output/species_routes_covariates/all_route_trends_grasslands_*.csv
##     (one row per route x species x model: alpha, beta, trend, trend_lci,
##     trend_uci, rel_abundance)
##   - output/species_routes_covariates/all_model_level_summary_grasslands_*.csv
##     (one row per species x model: gamma1 + its 90% CI)
##   - output/species_routes_covariates/per_species/<species>_<tag>_route_trends.csv
##     (same per-route data, split one file per species x model tag -- this
##     is what 3c reads)
##
## trend = 100 * (exp(beta) - 1) (annual % change; trend_lci/trend_uci are
## the same transform on the posterior's 5th/95th percentiles).
## rel_abundance = exp(alpha) (expected count at the fixed mid-series year).
## gamma1_excludes_zero is TRUE when its 90% CI doesn't span zero.
## =============================================================================

library(tidyverse)
library(here)

here::i_am("2c_generate_route_trend_csvs_covariates.R")

# Settings (match 1c_species_iCAR_covariates_grasslands.R) -----------------
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2025
model_tags <- c("base", "grassland_habitat", "grassland_anthro", "grassland_habitat_to_anthro")

output_dir       <- here::here("output")
rds_dir          <- here::here("output", "rds")
route_info_dir   <- here::here("data", "route_info")
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

  for (tag in model_tags) {

    out_base        <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file       <- file.path(rds_dir, paste0(out_base, "_summ_fit.rds"))
    route_info_file <- file.path(route_info_dir,
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_route_info.rds"))

    if (!file.exists(summ_file) || !file.exists(route_info_file)) {
      cat("  [", tag, "] No fitted output found — run 1c_species_iCAR_covariates_grasslands.R first. Skipping.\n")
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

    # gamma1 (covariate effect): one value per species+model, NA for "base".
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

    # Per-species-per-model CSV, so 3c can read one species/model per file.
    per_species_csv <- file.path(per_species_out_dir,
                                 paste0(sp_f, "_", tag, "_route_trends.csv"))
    write.csv(route_trends, per_species_csv, row.names = FALSE)

    cat("  [", tag, "] ", nrow(route_trends), "routes -> ", basename(per_species_csv), "\n")
  }
}

if (length(all_route_rows) == 0) {
  stop("Nothing to combine — run 1c_species_iCAR_covariates_grasslands.R first.")
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

cat("\nSpecies with gamma1's 90% CI excluding zero (credible effect), by model:\n")
print(all_model_summary %>%
        filter(model != "base") %>%
        group_by(model) %>%
        summarise(n_species          = n(),
                  n_excludes_zero    = sum(gamma1_excludes_zero, na.rm = TRUE),
                  .groups = "drop"))

cat("Next: run 3c_add_SDM_covariates.R to add climate-suitability columns.\n")
