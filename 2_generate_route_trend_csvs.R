## =============================================================================
## Generate per-route trend CSVs from fitted iCAR models
## Time period: 2010-2025
##
## Reads the per-species model output saved by 1_species_iCAR_2010_2025.R:
##   output/<species>_iCAR_NB_<firstYear>_<lastYear>_summ_fit.rds   (summary)
##   data/stan_data/<species>_<firstYear>_<lastYear>_stan_data.RData
##     (loads route_map (sf), new_data, realized_strata_map, ...)
##
## and writes one per-route CSV per species (mirrors output/species_routes
## reference layout):
##   output/species_routes/<species>_route_trends.csv  with columns:
##     species, species_code, route, latitude, longitude,
##     alpha, beta, trend, trend_lci, trend_uci, rel_abundance
##
## where, per route s (alpha[s] and beta[s] are the raw Stan posterior means,
## on the model's log scale; trend and rel_abundance are their back-transformed,
## interpretable versions):
##   alpha         = mean(alpha[s])                the fitted log-scale intercept
##                                                  (mean log-count) for route s
##   beta          = mean(beta[s])                 the fitted log-scale slope
##                                                  (annual rate of change) for route s
##   trend     = 100 * (exp(beta[s])  - 1)   annual % change, i.e. beta converted
##                                            from a log-scale slope to a
##                                            percent-per-year change
##   trend_lci = 100 * (exp(q5(beta[s]))    - 1)   90% CI lower bound on trend
##   trend_uci = 100 * (exp(q95(beta[s]))   - 1)   90% CI upper bound on trend
##   rel_abundance = exp(alpha[s])            relative abundance, i.e. alpha
##                                            converted from a log-scale mean
##                                            back to the count scale (expected
##                                            count for an average observer at
##                                            the route's fixed year)
##
## This is a post-processing step only; it does not fit or load any models'
## posterior draws beyond the saved summary table.
## =============================================================================

library(tidyverse)
library(sf)
library(here)

here::i_am("2_generate_route_trend_csvs.R")

# Settings (match 1_species_iCAR_2010_2025.R) -----------------------------
land_cover <- "grasslands"
firstYear  <- 2010
lastYear   <- 2025

# Overwrite existing per-species CSVs?
overwrite <- TRUE

# Directories
output_dir <- here::here("output")
sp_out_dir <- here::here("output", "species_routes")
if (!dir.exists(sp_out_dir)) dir.create(sp_out_dir, recursive = TRUE)

# Target group species list -----------------------------------------------
spp_df <- read.csv(here::here("data", "spp_names_codes_group_aou.csv"),
                   stringsAsFactors = FALSE)

target_spp <- spp_df %>%
  filter(Group == land_cover, in_bbs == TRUE) %>%
  distinct(Common.Name, Code, .keep_all = TRUE) %>%
  arrange(Common.Name)

cat("=== Generate per-route trend CSVs ===\n")
cat("Group:", land_cover, " | Period:", firstYear, "-", lastYear, "\n")
cat("Species in group (n =", nrow(target_spp), ")\n")

# Helper: convert species name to file-safe format ------------------------
species_to_f <- function(sp) {
  gsub("'", "", gsub(" ", "_", sp, fixed = TRUE), fixed = TRUE)
}

# ==========================================================================
# Main loop: one CSV per species
# ==========================================================================
results_list <- list()

for (i in seq_len(nrow(target_spp))) {
  sp      <- target_spp$Common.Name[i]
  sp_f    <- species_to_f(sp)
  sp_code <- target_spp$Code[i]

  cat("\n[", i, "/", nrow(target_spp), "]", sp, "\n")

  sp_csv <- file.path(sp_out_dir, paste0(sp_f, "_route_trends.csv"))
  if (file.exists(sp_csv) && !overwrite) {
    cat("  Skipping (already exists):", basename(sp_csv), "\n")
    next
  }

  # Paths to pre-fitted output written by script 1
  out_base       <- paste0(sp_f, "_iCAR_NB_", firstYear, "_", lastYear)
  summ_file      <- file.path(output_dir, paste0(out_base, "_summ_fit.rds"))
  stan_data_file <- here::here("data", "stan_data",
                               paste0(sp_f, "_", firstYear, "_", lastYear, "_stan_data.RData"))

  if (!file.exists(summ_file) || !file.exists(stan_data_file)) {
    cat("  No fitted output found — run 1_species_iCAR_2010_2025.R first. Skipping.\n")
    next
  }

  summ <- readRDS(summ_file)
  load(stan_data_file)   # loads route_map (sf), new_data, realized_strata_map, ...

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

  # Lat/lon from route_map (sf object) -------------------------------------
  route_info <- route_map %>%
    st_transform(4326) %>%
    mutate(longitude = st_coordinates(.)[, 1],
           latitude  = st_coordinates(.)[, 2]) %>%
    st_drop_geometry() %>%
    select(route, routeF, latitude, longitude)

  # Assemble per-route table -----------------------------------------------
  route_trends <- beta_summ %>%
    left_join(alpha_summ, by = "routeF") %>%
    left_join(route_info, by = "routeF") %>%
    mutate(species      = sp,
           species_code = sp_code) %>%
    arrange(routeF) %>%
    select(species, species_code, route,
           latitude, longitude, alpha, beta, trend, trend_lci, trend_uci, rel_abundance)

  write.csv(route_trends, sp_csv, row.names = FALSE)
  results_list[[sp]] <- route_trends

  cat("  Median trend:", round(median(route_trends$trend), 2), "% per year",
      " | Done:", nrow(route_trends), "routes ->", basename(sp_csv), "\n")
}

cat("\n=== Summary ===\n")
cat("CSVs generated this run:", length(results_list), "\n")
cat("Per-species trend CSVs in:", sp_out_dir, "\n")
