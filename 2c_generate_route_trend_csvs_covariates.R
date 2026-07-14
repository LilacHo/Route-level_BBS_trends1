## =============================================================================
## Combine per-route trend output across ALL fitted models into one giant CSV,
## for cross-model comparison.
## Time period: 2010-2025
##
## Combines four "model" sources per species:
##   base                 — no covariates, from 1_species_iCAR_2010_2025.R,
##                          read from the CSV already written by
##                          2_generate_route_trend_csvs.R
##                          (output/species_routes/<species>_route_trends.csv)
##   grassland            — from 1c_species_iCAR_covariates.R
##   developed            — from 1c_species_iCAR_covariates.R
##   grassland_developed  — from 1c_species_iCAR_covariates.R
##
## For the three covariate models, reads the summary fit + stan_data saved by
## 1c_species_iCAR_covariates.R:
##   output/<species>_iCAR_<tag>_<firstYear>_<lastYear>_summ_fit.rds
##   data/stan_data/<species>_<tag>_<firstYear>_<lastYear>_stan_data.RData
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
##      columns: species, species_code, model, n_routes, gamma1, gamma2
##
## gamma1/gamma2 are the covariate effect(s) on log(lambda) — a single value
## per species+model, not per-route, which is why they live in the separate
## model-level file rather than being repeated down every route row. They are
## NA for the "base" model (no covariates) and gamma2 is NA for the two
## single-covariate models (grassland, developed).
##
## This is a post-processing / combination step only; it does not fit or
## re-derive anything beyond what 2_generate_route_trend_csvs.R and
## 1c_species_iCAR_covariates.R already produced.
## =============================================================================

library(tidyverse)
library(here)

here::i_am("2c_generate_route_trend_csvs_covariates.R")

# Settings (match 1_species_iCAR_2010_2025.R / 1c_species_iCAR_covariates.R) -
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2025
covariate_model_tags <- c("grassland", "developed", "grassland_developed")

output_dir      <- here::here("output")
base_route_dir  <- here::here("output", "species_routes")             # from step 2
combined_out_dir <- here::here("output", "species_routes_covariates")
if (!dir.exists(combined_out_dir)) dir.create(combined_out_dir, recursive = TRUE)

# Target group species list -------------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

# Optional external override — set `species_filter <- c("Some Species", ...)`
# in the calling environment before sourcing this script to restrict the run
# to just those species (used by tests/test_single_species.R).
if (exists("species_filter")) {
  target_spp <- target_spp %>% filter(Common.Name %in% species_filter)
  cat("species_filter applied — running for:", paste(species_filter, collapse = ", "), "\n")
}

cat("=== Combine route-trend CSVs across all models ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("Species in group (n =", nrow(target_spp), ")\n")
cat("Models: base,", paste(covariate_model_tags, collapse = ", "), "\n")

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# Column order every row is coerced to before combining --------------------
route_cols <- c("species", "species_code", "model", "route",
                "latitude", "longitude", "alpha", "beta",
                "trend", "trend_lci", "trend_uci", "rel_abundance")
model_cols <- c("species", "species_code", "model", "n_routes",
                "gamma1", "gamma2")

all_route_rows <- list()
all_model_rows <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]

  cat("\n[", i, "/", nrow(target_spp), "]", sp, "\n")

  # --- base model: read the CSV already written by 2_generate_route_trend_csvs.R
  base_csv <- file.path(base_route_dir, paste0(sp_f, "_route_trends.csv"))
  if (file.exists(base_csv)) {
    base_trends <- read.csv(base_csv, stringsAsFactors = FALSE) %>%
      mutate(model = "base") %>%
      select(all_of(route_cols))
    all_route_rows[[paste(sp, "base", sep = " | ")]] <- base_trends

    all_model_rows[[paste(sp, "base", sep = " | ")]] <- data.frame(
      species      = sp,
      species_code = sp_code,
      model        = "base",
      n_routes     = nrow(base_trends),
      gamma1       = NA_real_,
      gamma2       = NA_real_
    )
    cat("  [base] ", nrow(base_trends), "routes\n")
  } else {
    cat("  [base] No CSV found (", basename(base_csv),
        ") — run 2_generate_route_trend_csvs.R first. Skipping base for this species.\n")
  }

  # --- covariate models: read summ_fit + stan_data written by 1c
  for (tag in covariate_model_tags) {

    out_base       <- paste0(sp_f, "_iCAR_", tag, "_", firstYear, "_", lastYear)
    summ_file      <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
    stan_data_file <- here::here("data", "stan_data",
                                 paste0(sp_f, "_", tag, "_", firstYear, "_", lastYear, "_stan_data.RData"))

    if (!file.exists(summ_file) || !file.exists(stan_data_file)) {
      cat("  [", tag, "] No fitted output found — run 1c_species_iCAR_covariates.R first. Skipping.\n")
      next
    }

    summ <- readRDS(summ_file)
    load(stan_data_file)   # loads stan_data, new_data, model_tag

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

    # Species+model-level covariate effect(s), not per-route -----------------
    gamma1_val <- summ$mean[summ$variable == "gamma1"]
    gamma2_val <- if (tag == "grassland_developed") {
      summ$mean[summ$variable == "gamma2"]
    } else {
      NA_real_
    }

    # Route lat/lon straight from new_data (already WGS84, no reprojection
    # needed since 1c builds it directly from raw_data$latitude/longitude)
    route_info <- new_data %>%
      distinct(route, routeF, latitude, longitude)

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
      species      = sp,
      species_code = sp_code,
      model        = tag,
      n_routes     = nrow(route_trends),
      gamma1       = gamma1_val,
      gamma2       = gamma2_val
    )
    cat("  [", tag, "] ", nrow(route_trends), "routes\n")
  }
}

if (length(all_route_rows) == 0) {
  stop("Nothing to combine — run 2_generate_route_trend_csvs.R and/or ",
       "1c_species_iCAR_covariates.R first.")
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
